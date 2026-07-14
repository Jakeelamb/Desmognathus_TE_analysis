#!/usr/bin/env python3
"""Join collaborator SVL measurements to the finalized cell/genome panel.

This is a descriptive species-level analysis. It does not treat the species
as independent experimental replicates or make a causal claim.
"""

from pathlib import Path
from zipfile import ZipFile
from xml.etree import ElementTree

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
SVL_PATH = ROOT / "input_data/morphological_data/desmog_SVL.xlsx"
TRAIT_PATH = ROOT / "results/data/research_review/cell_nucleus_genome_path/primary_traits_all_finalized.csv"
OUT_DIR = ROOT / "results/data/research_review/svl_cell_nucleus_genome"


def rank_average(values: np.ndarray) -> np.ndarray:
    return pd.Series(values).rank(method="average").to_numpy()


def corr(x: np.ndarray, y: np.ndarray) -> float:
    return float(np.corrcoef(x, y)[0, 1])


def permutation_p(x: np.ndarray, y: np.ndarray, statistic: str, seed: int = 20260714) -> float:
    observed = abs(corr(x, y) if statistic == "pearson" else corr(rank_average(x), rank_average(y)))
    rng = np.random.default_rng(seed)
    permuted = np.empty(20_000)
    for i in range(permuted.size):
        shuffled = rng.permutation(y)
        permuted[i] = abs(corr(x, shuffled) if statistic == "pearson" else corr(rank_average(x), rank_average(shuffled)))
    return float((np.count_nonzero(permuted >= observed) + 1) / (permuted.size + 1))


def bootstrap_ci(x: np.ndarray, y: np.ndarray, statistic: str, seed: int = 20260714) -> tuple[float, float]:
    rng = np.random.default_rng(seed)
    estimates = np.empty(20_000)
    for i in range(estimates.size):
        indices = rng.integers(0, len(x), len(x))
        xb, yb = x[indices], y[indices]
        estimates[i] = corr(xb, yb) if statistic == "pearson" else corr(rank_average(xb), rank_average(yb))
    return tuple(np.quantile(estimates, [0.025, 0.975]))


def read_first_xlsx_sheet(path: Path) -> pd.DataFrame:
    """Read the simple, value-only collaborator workbook without an Excel engine."""
    ns = {"m": "http://schemas.openxmlformats.org/spreadsheetml/2006/main",
          "r": "http://schemas.openxmlformats.org/officeDocument/2006/relationships"}
    with ZipFile(path) as archive:
        shared = []
        if "xl/sharedStrings.xml" in archive.namelist():
            root = ElementTree.fromstring(archive.read("xl/sharedStrings.xml"))
            for item in root.findall("m:si", ns):
                shared.append("".join(node.text or "" for node in item.iter() if node.tag.endswith("}t")))
        workbook = ElementTree.fromstring(archive.read("xl/workbook.xml"))
        rels = ElementTree.fromstring(archive.read("xl/_rels/workbook.xml.rels"))
        rel_map = {rel.attrib["Id"]: rel.attrib["Target"] for rel in rels}
        sheet = workbook.find("m:sheets/m:sheet", ns)
        target = rel_map[sheet.attrib["{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id"]]
        sheet_path = "xl/" + target.lstrip("/")
        root = ElementTree.fromstring(archive.read(sheet_path))
        rows = []
        for row in root.findall(".//m:sheetData/m:row", ns):
            values = {}
            for cell in row.findall("m:c", ns):
                ref = cell.attrib["r"]
                col = 0
                for char in ref:
                    if char.isalpha():
                        col = col * 26 + ord(char.upper()) - ord("A") + 1
                    else:
                        break
                value = cell.find("m:v", ns)
                text = "" if value is None else value.text
                if cell.attrib.get("t") == "s" and text:
                    text = shared[int(text)]
                values[col] = text
            rows.append([values.get(i, "") for i in range(1, max(values, default=0) + 1)])
    header, body = rows[0], rows[1:]
    return pd.DataFrame(body, columns=header)


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    svl = read_first_xlsx_sheet(SVL_PATH)
    svl["maxSVL"] = pd.to_numeric(svl["maxSVL"], errors="coerce")
    svl = svl.loc[svl["Genus"].eq("Desmognathus"), ["Species", "maxSVL"]]
    svl = svl.rename(columns={"Species": "tree_tip", "maxSVL": "svl_mm"})
    traits = pd.read_csv(TRAIT_PATH)
    traits["tree_tip"] = traits["tree_tip"].str.replace(r"^D\\. ", "", regex=True)
    data = traits.merge(svl, on="tree_tip", how="left", validate="one_to_one")
    if data["svl_mm"].isna().any():
        missing = ", ".join(data.loc[data["svl_mm"].isna(), "tree_tip"])
        raise ValueError(f"Missing collaborator SVL for: {missing}")

    keep = ["species", "tree_tip", "svl_mm", "genome_size_pg", "nucleus_area_um2", "cell_area_um2",
            "n_genome_nuclei", "n_genome_images", "n_size_cells", "n_size_images"]
    data[keep].to_csv(OUT_DIR / "svl_cell_nucleus_genome_species_panel.csv", index=False)

    relationships = [
        ("genome_size_pg", "Genome size", "pg", "genome"),
        ("cell_area_um2", "Cell area", "µm²", "cell"),
        ("nucleus_area_um2", "Nucleus area", "µm²", "nucleus"),
    ]
    rows = []
    x = data["svl_mm"].to_numpy(float)
    for column, label, units, short in relationships:
        y = data[column].to_numpy(float)
        pearson = corr(x, y)
        spearman = corr(rank_average(x), rank_average(y))
        p_pearson = permutation_p(x, y, "pearson", seed=20260714)
        p_spearman = permutation_p(x, y, "spearman", seed=20260715)
        pearson_low, pearson_high = bootstrap_ci(x, y, "pearson", seed=20260716)
        spearman_low, spearman_high = bootstrap_ci(x, y, "spearman", seed=20260717)
        rows.append({
            "trait": short,
            "trait_label": label,
            "trait_units": units,
            "n_species": len(data),
            "pearson_r": pearson,
            "pearson_permutation_p": p_pearson,
            "pearson_bootstrap_ci_low": pearson_low,
            "pearson_bootstrap_ci_high": pearson_high,
            "spearman_rho": spearman,
            "spearman_permutation_p": p_spearman,
            "spearman_bootstrap_ci_low": spearman_low,
            "spearman_bootstrap_ci_high": spearman_high,
        })
    summary = pd.DataFrame(rows)
    summary.to_csv(OUT_DIR / "svl_correlation_summary.csv", index=False)

    plt.style.use("seaborn-v0_8-whitegrid")
    fig, axes = plt.subplots(1, 3, figsize=(15, 5), constrained_layout=True)
    colors = {"genome": "#496a8f", "cell": "#4b8b73", "nucleus": "#a65d5d"}
    for ax, (column, label, units, short), result in zip(axes, relationships, rows):
        y = data[column].to_numpy(float)
        color = colors[short]
        ax.scatter(x, y, s=42, color=color, edgecolor="white", linewidth=0.7, zorder=3)
        slope, intercept = np.polyfit(x, y, 1)
        xx = np.linspace(x.min() - 3, x.max() + 3, 100)
        ax.plot(xx, slope * xx + intercept, color=color, linewidth=1.6, alpha=0.8)
        for xi, yi, name in zip(x, y, data["tree_tip"]):
            ax.annotate(name, (xi, yi), xytext=(3, 3), textcoords="offset points", fontsize=6.5, alpha=0.82)
        ax.set_xlabel("Maximum SVL (mm)")
        ax.set_ylabel(f"{label} ({units})")
        ax.set_title(label)
        ax.text(0.04, 0.96,
                f"Pearson r = {result['pearson_r']:.2f}\n"
                f"Spearman ρ = {result['spearman_rho']:.2f}",
                transform=ax.transAxes, va="top", fontsize=9,
                bbox={"facecolor": "white", "edgecolor": "0.8", "alpha": 0.9})
    fig.suptitle("Species-level SVL associations with genome and cell traits", fontsize=14)
    fig.text(0.5, -0.02,
             "Points are the species in the current primary cell/genome panel; lines are descriptive OLS fits."
             " SVL is collaborator-provided maximum SVL.", ha="center", fontsize=9)
    fig.savefig(OUT_DIR / "svl_cell_nucleus_genome_correlations.png", dpi=300, bbox_inches="tight")
    fig.savefig(OUT_DIR / "svl_cell_nucleus_genome_correlations.pdf", bbox_inches="tight")

    print(summary.to_string(index=False, float_format=lambda value: f"{value:.4f}"))
    print(f"Wrote {len(data)} species and {len(summary)} correlations to {OUT_DIR}")


if __name__ == "__main__":
    main()
