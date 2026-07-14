#!/usr/bin/env python3
"""Build TE34 descriptive products from frozen upstream outputs only.

This script deliberately separates the vetted active genomic-resource panel
from the linked-cell and integrated-path panels.  It reads existing dnaPipeTE
and corrected RepeatMasker tables; it never launches RepeatModeler,
RepeatMasker, dnaPipeTE, or read-mapping work.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
PANEL = ROOT / "path_analysis/data/derived/panels/study_te_resource_panel34_v1.csv"
LOOKUP = ROOT / "input_data/lookup_table.txt"
MERGED = ROOT / "results/data/dnaPipeTE_merged_classifications.csv"
ORDER = ROOT / "results/data/dnaPipeTE_order_breakdown.csv"
SUPERFAMILY = ROOT / "results/data/dnaPipeTE_superfamily_breakdown.csv"
HITS = ROOT / "results/data/corrected/repeatmasker_detailed_classification_hit_level_v1.csv"
HITS_MANIFEST = Path(str(HITS) + ".manifest.json")
OUT = ROOT / "results/data/corrected/te34"
FIG = ROOT / "results/figures/corrected/te34"


def canonical(value: object) -> str:
    return str(value).strip().replace("Desmognathus ", "").replace("D. ", "").replace("D.", "").lower()


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def diversity(values: np.ndarray) -> dict[str, float]:
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values) & (values > 0)]
    p = values / values.sum()
    h = float(-(p * np.log(p)).sum())
    dominance = float((p**2).sum())
    richness = len(p)
    return {
        "observed_richness": richness,
        "shannon_entropy": h,
        "gini_simpson": 1 - dominance,
        "hill_q1": float(np.exp(h)),
        "hill_q2": 1 / dominance,
        "pielou_evenness": h / np.log(richness) if richness > 1 else 0.0,
    }


def run_pca(composition: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    table = composition.copy().astype(float)
    zeros = table.eq(0)
    if zeros.any().any():
        table = table.mask(zeros, float(table.where(table.gt(0)).min().min() / 2))
        table = table.div(table.sum(axis=1), axis=0)
    clr = np.log(table).sub(np.log(table).mean(axis=1), axis=0)
    centered = clr.to_numpy() - clr.to_numpy().mean(axis=0)
    u, singular, vt = np.linalg.svd(centered, full_matrices=False)
    n_components = min(len(clr) - 1, clr.shape[1] - 1)
    scores = pd.DataFrame(u[:, :n_components] * singular[:n_components], index=clr.index)
    scores.columns = [f"PC{i + 1}" for i in range(n_components)]
    loadings = pd.DataFrame(vt[:n_components].T, index=clr.columns, columns=scores.columns)
    for column in scores:
        anchor = loadings[column].abs().idxmax()
        if loadings.loc[anchor, column] < 0:
            loadings[column] *= -1
            scores[column] *= -1
    return scores, loadings


def read_panel() -> tuple[list[str], pd.DataFrame]:
    panel = pd.read_csv(PANEL)
    panel["species"] = panel["species"].map(canonical)
    if len(panel) != 34 or panel["species"].duplicated().any():
        raise ValueError("TE resource panel must have exactly 34 unique species")
    if not panel[["has_tree_tip", "has_te"]].all().all():
        raise ValueError("TE34 panel has incomplete tree or TE evidence")
    return panel["species"].sort_values().tolist(), panel.set_index("species")


def mass_ledger(species: set[str], panel: pd.DataFrame) -> pd.DataFrame:
    headers = {
        level: set(pd.read_csv(path, nrows=0, index_col=0).columns)
        for level, path in [("Order", ORDER), ("Superfamily", SUPERFAMILY)]
    }
    totals: dict[str, dict[str, object]] = {
        item: {"n_component_rows": 0, "total_aligned_bases": 0, "sources": set()}
        for item in species
    }
    retained = {(item, level): 0 for item in species for level in headers}
    usecols = ["Species", "Source", "aligned_bases", "Order", "Superfamily"]
    for chunk in pd.read_csv(MERGED, usecols=usecols, chunksize=400_000, low_memory=False):
        chunk["species"] = chunk["Species"].map(canonical)
        chunk = chunk.loc[chunk["species"].isin(species)].copy()
        # dnaPipeTE retry directories can append labels such as ``R2`` to the
        # same SRX.  The active lookup accession is the identity boundary;
        # never sum a retry alongside its original run.
        expected_source = chunk["species"].map(panel["te_sra_accession"])
        chunk = chunk.loc[chunk["Source"].astype(str).eq(expected_source)].copy()
        if chunk.empty:
            continue
        chunk["aligned_bases"] = pd.to_numeric(chunk["aligned_bases"], errors="raise")
        if chunk["aligned_bases"].lt(0).any():
            raise ValueError("dnaPipeTE aligned bases cannot be negative")
        for item, group in chunk.groupby("species"):
            totals[item]["n_component_rows"] += len(group)
            totals[item]["total_aligned_bases"] += int(group["aligned_bases"].sum())
            totals[item]["sources"].update(group["Source"].dropna().astype(str))
        for level, labels in headers.items():
            grouped = (
                chunk.loc[chunk[level].astype("string").isin(labels)]
                .groupby("species")["aligned_bases"].sum()
            )
            for item, value in grouped.items():
                retained[(item, level)] += int(value)
    rows = []
    for item in sorted(species):
        source = totals[item]["sources"]
        if len(source) != 1:
            raise ValueError(f"{item} maps to {len(source)} dnaPipeTE SRX resources")
        total = int(totals[item]["total_aligned_bases"])
        if total <= 0:
            raise ValueError(f"{item} has no dnaPipeTE aligned bases")
        row = {
            "species": item,
            "te_sra_accession": panel.loc[item, "te_sra_accession"],
            "te_assembly_accession": panel.loc[item, "te_assembly_accession"],
            "observed_dnapipete_sra_accession": next(iter(source)),
            "n_component_rows": int(totals[item]["n_component_rows"]),
            "total_aligned_bases": total,
        }
        if row["observed_dnapipete_sra_accession"] != row["te_sra_accession"]:
            raise ValueError(f"dnaPipeTE accession mismatch for {item}")
        for level in headers:
            prefix = level.lower()
            kept = retained[(item, level)]
            row[f"{prefix}_retained_aligned_bases"] = kept
            row[f"{prefix}_unresolved_aligned_bases"] = total - kept
            row[f"{prefix}_retained_fraction"] = kept / total
            row[f"{prefix}_unresolved_fraction"] = 1 - kept / total
        rows.append(row)
    return pd.DataFrame(rows)


def compositions_and_pca(species: list[str], mass: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows: list[dict[str, object]] = []
    score_rows: list[pd.DataFrame] = []
    mass = mass.set_index("species")
    for level, path in [("order", ORDER), ("superfamily", SUPERFAMILY)]:
        table = pd.read_csv(path, index_col=0)
        table.index = table.index.map(canonical)
        table = table.reindex(species).apply(pd.to_numeric, errors="raise")
        if table.isna().any().any() or table.sum(axis=1).le(0).any():
            raise ValueError(f"{level} composition is missing a TE34 resource")
        classified = table.div(table.sum(axis=1), axis=0)
        mass_aware = classified.mul(mass[f"{level}_retained_fraction"], axis=0)
        mass_aware["Unresolved"] = mass[f"{level}_unresolved_fraction"]
        for mode, composition in [("classified_conditional", classified), ("mass_aware_unresolved_bin", mass_aware)]:
            for item, values in composition.iterrows():
                rows.append({"species": item, "te_level": level, "composition_mode": mode, **diversity(values.to_numpy())})
        scores, _ = run_pca(classified)
        scores = scores.reset_index().rename(columns={"index": "species"})
        scores.insert(1, "te_level", level)
        scores.insert(2, "composition_mode", "classified_conditional")
        score_rows.append(scores)
    return pd.DataFrame(rows), pd.concat(score_rows, ignore_index=True)


def repeatmasker_products(species: set[str]) -> tuple[pd.DataFrame, pd.DataFrame]:
    summaries = []
    columns = ["Desmognathus_Species", "percent_divergence", "query_start", "query_end", "Order"]
    for chunk in pd.read_csv(HITS, usecols=columns, chunksize=400_000, low_memory=False):
        chunk["species"] = chunk["Desmognathus_Species"].map(canonical)
        chunk = chunk.loc[chunk["species"].isin(species)].copy()
        chunk["percent_divergence"] = pd.to_numeric(chunk["percent_divergence"], errors="coerce")
        chunk["hit_bp"] = (pd.to_numeric(chunk["query_end"], errors="coerce") - pd.to_numeric(chunk["query_start"], errors="coerce")).abs() + 1
        chunk = chunk.loc[chunk["percent_divergence"].ge(0) & chunk["hit_bp"].gt(0)].copy()
        chunk["order"] = chunk["Order"].fillna("Unclassified").astype(str).str.strip().replace("", "Unclassified")
        chunk["divergence_bin_start_pct"] = np.floor(chunk["percent_divergence"].clip(upper=50)).astype(int)
        summaries.append(chunk.groupby(["species", "order", "divergence_bin_start_pct"], as_index=False).agg(hit_count=("hit_bp", "size"), hit_bp=("hit_bp", "sum")))
    landscape = pd.concat(summaries, ignore_index=True).groupby(["species", "order", "divergence_bin_start_pct"], as_index=False)[["hit_count", "hit_bp"]].sum()
    landscape["species_total_hit_bp"] = landscape.groupby("species")["hit_bp"].transform("sum")
    landscape["fraction_species_hit_bp"] = landscape["hit_bp"] / landscape["species_total_hit_bp"]
    landscape["percent_species_hit_bp"] = 100 * landscape["fraction_species_hit_bp"]
    inventory = landscape.groupby("species", as_index=False).agg(hit_count=("hit_count", "sum"), aligned_hit_bp=("hit_bp", "sum"), n_orders=("order", "nunique"))
    if set(inventory["species"]) != species:
        raise ValueError("RepeatMasker TE34 species are incomplete")
    source = json.loads(HITS_MANIFEST.read_text())
    if int(inventory["hit_count"].sum()) != int(source["n_hits"]) or int(inventory["aligned_hit_bp"].sum()) != int(source["inclusive_aligned_bp"]):
        raise ValueError("TE34 RepeatMasker accounting does not conserve frozen all-resource table")
    if not np.allclose(landscape.groupby("species")["fraction_species_hit_bp"].sum(), 1.0, atol=1e-12):
        raise ValueError("TE34 divergence landscape fractions do not close within species")
    return inventory.sort_values("species"), landscape.sort_values(["species", "order", "divergence_bin_start_pct"])


def write_figures(diversity_table: pd.DataFrame, scores: pd.DataFrame, landscape: pd.DataFrame) -> list[Path]:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    FIG.mkdir(parents=True, exist_ok=True)
    paths = []
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    for ax, level in zip(axes, ["order", "superfamily"]):
        data = diversity_table.loc[(diversity_table.te_level == level) & (diversity_table.composition_mode == "classified_conditional")]
        ax.hist(data["pielou_evenness"], bins=10, color="#0f766e", edgecolor="white")
        ax.set(title=f"TE34 {level} Pielou evenness", xlabel="evenness", ylabel="species")
    path = FIG / "te34_diversity_distributions_v1.png"; fig.tight_layout(); fig.savefig(path, dpi=300); plt.close(fig); paths.append(path)
    fig, ax = plt.subplots(figsize=(8, 6))
    data = scores.loc[scores.te_level.eq("order")]
    ax.scatter(data.PC1, data.PC2, color="#2563eb")
    for row in data.itertuples(): ax.annotate(row.species, (row.PC1, row.PC2), fontsize=7, xytext=(2, 2), textcoords="offset points")
    ax.set(title="TE34 order-composition CLR PCA", xlabel="PC1", ylabel="PC2"); fig.tight_layout()
    path = FIG / "te34_order_clr_pca_v1.png"; fig.savefig(path, dpi=300); plt.close(fig); paths.append(path)
    fig, ax = plt.subplots(figsize=(10, 6))
    low = landscape.loc[landscape.divergence_bin_start_pct.le(10)].groupby("species").hit_bp.sum().sort_values()
    ax.barh(low.index, low.values / 1_000_000, color="#9333ea")
    ax.set(title="TE34 young-repeat landscape summary (0-10% divergence)", xlabel="inclusive aligned RepeatMasker hit bp (Mb)")
    path = FIG / "te34_young_repeat_landscape_v1.png"; fig.tight_layout(); fig.savefig(path, dpi=300); plt.close(fig); paths.append(path)
    return paths


def main() -> None:
    outputs = [
        OUT / "dnapipete_mass_accounting_te34_v1.csv", OUT / "te_diversity_mass_sensitivity_te34_v1.csv",
        OUT / "te_pca_scores_te34_v1.csv", OUT / "repeatmasker_hit_inventory_te34_v1.csv",
        OUT / "repeatmasker_divergence_landscape_te34_v1.csv", OUT / "te34_descriptive_bundle_v1.manifest.json",
    ]
    if any(path.exists() for path in outputs):
        raise FileExistsError("TE34 outputs already exist; refusing to overwrite a descriptive audit")
    species, panel = read_panel()
    mass = mass_ledger(set(species), panel)
    diversity_table, scores = compositions_and_pca(species, mass)
    inventory, landscape = repeatmasker_products(set(species))
    OUT.mkdir(parents=True, exist_ok=True)
    for path, frame in zip(outputs[:5], [mass, diversity_table, scores, inventory, landscape]): frame.to_csv(path, index=False)
    figures = write_figures(diversity_table, scores, landscape)
    manifest = {
        "analysis_id": "te34_descriptive_bundle_v1", "analysis_scope": "active_vetted_te_resource_panel34_descriptive",
        "n_te_species": 34, "species": species, "panel": str(PANEL.relative_to(ROOT)), "panel_sha256": digest(PANEL),
        "inputs": {str(path.relative_to(ROOT)): digest(path) for path in [MERGED, ORDER, SUPERFAMILY, HITS, HITS_MANIFEST]},
        "repeatmasker": {"n_hits": int(inventory.hit_count.sum()), "inclusive_aligned_bp": int(inventory.aligned_hit_bp.sum())},
        "expensive_upstream_tools_executed": False, "analysis18_outputs_preserved": True,
        "eligible_for_integrated_phylogenetic_path_analysis": False,
        "path_analysis_scope": "study_integrated_path_panel18_v1 only; no cell or IOD traits imputed for TE-only species",
        "dnapipete_retry_policy": "Only the exact active lookup SRX is included; retry-labelled source directories (for example SRXR2) are not summed as a second biological resource.",
        "outputs": {str(path.relative_to(ROOT)): digest(path) for path in outputs[:5]},
        "figures": [{"path": str(path.relative_to(ROOT)), "sha256": digest(path)} for path in figures],
    }
    outputs[-1].write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps({"n_te_species": 34, **manifest["repeatmasker"]}, indent=2))


if __name__ == "__main__":
    main()
