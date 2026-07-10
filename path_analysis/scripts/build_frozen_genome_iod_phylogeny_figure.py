#!/usr/bin/env python3
"""Plot the measured genome-IOD, nucleus, and cell panels beside the time tree."""

from __future__ import annotations

import copy
import hashlib
import json
import sys
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import numpy as np
import pandas as pd
from scipy.stats import spearmanr


PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT / "path_analysis" / "scripts"))
import build_frozen_genome_iod_notebook as genome_analysis  # noqa: E402


GENOME_PATH = genome_analysis.FROZEN_PATH
GENOME_SUMMARY_PATH = genome_analysis.SPECIES_SUMMARY_PATH
CELL_PATH = (
    PROJECT_ROOT
    / "path_analysis"
    / "data"
    / "external"
    / "derived"
    / "largest_cell_mask_review"
    / "frozen_largest_cell_mask_top50.csv.gz"
)
TREE_PATH = PROJECT_ROOT / "input_data" / "phylogeny" / "desmo900dated_test.tre"
OUTPUT_DIR = genome_analysis.REPORT_DIR
FIGURE_DIR = genome_analysis.FIGURE_DIR
PNG_PATH = FIGURE_DIR / "07_measured_phylogeny_genome_nucleus_cell.png"
PDF_PATH = FIGURE_DIR / "07_measured_phylogeny_genome_nucleus_cell.pdf"
SUMMARY_PATH = OUTPUT_DIR / "phylogeny_genome_nucleus_cell_summary.csv"
CORRELATION_PATH = OUTPUT_DIR / "phylogeny_genome_nucleus_cell_correlations.csv"
MANIFEST_PATH = OUTPUT_DIR / "phylogeny_genome_nucleus_cell_manifest.json"
EXECUTED_NOTEBOOK_PATH = genome_analysis.EXECUTED_NOTEBOOK_PATH
METRIC_ORDER = (
    "relative_iod_index",
    "nucleus_area_um2",
    "cell_area_um2",
)
METRIC_STYLE = {
    "relative_iod_index": {
        "title": "Relative nuclear-IOD",
        "xlabel": "Relative IOD index\n(species median = 1.0)",
        "color": "#78A6C0",
        "format": ".2f",
    },
    "nucleus_area_um2": {
        "title": "Nucleus area",
        "xlabel": "Nucleus area (µm²)\nmedian of vetted top-50 cells",
        "color": "#9A8CB4",
        "format": ".1f",
    },
    "cell_area_um2": {
        "title": "Cell area",
        "xlabel": "Cell area (µm²)\nmedian of vetted top-50 cells",
        "color": "#D39768",
        "format": ".1f",
    },
}


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def canonical_species(value: Any) -> str:
    text = str(value).strip()
    if text.startswith("D. "):
        return text
    if text.startswith("Desmognathus "):
        return "D. " + text.split(" ", 1)[1]
    return "D. " + text


def bootstrap_equal_image_iod(
    frame: pd.DataFrame,
    *,
    n_bootstrap: int,
    seed: int,
) -> np.ndarray:
    """Draw the mean of image medians by resampling images then nuclei."""
    images = {
        str(filename): group["nuc_iod"].dropna().to_numpy(dtype=float)
        for filename, group in frame.groupby("filename", sort=True)
    }
    images = {name: values for name, values in images.items() if len(values)}
    if not images:
        raise ValueError("No usable image-specific IOD values.")
    names = np.array(sorted(images), dtype=object)
    rng = np.random.default_rng(seed)
    draws = np.empty(n_bootstrap, dtype=float)
    for replicate in range(n_bootstrap):
        sampled_names = rng.choice(names, size=len(names), replace=True)
        medians = []
        for name in sampled_names:
            values = images[str(name)]
            medians.append(
                float(np.median(rng.choice(values, size=len(values), replace=True)))
            )
        draws[replicate] = float(np.mean(medians))
    return draws


def bootstrap_paired_size_medians(
    frame: pd.DataFrame,
    *,
    n_bootstrap: int,
    seed: int,
) -> pd.DataFrame:
    """Bootstrap paired cell and corresponding-nucleus medians."""
    values = frame[["cell_area_um2", "nuc_area_um2"]].to_numpy(dtype=float)
    if len(values) == 0 or not np.isfinite(values).all():
        raise ValueError("Size panel contains no usable paired measurements.")
    rng = np.random.default_rng(seed)
    indices = rng.integers(0, len(values), size=(n_bootstrap, len(values)))
    sampled = values[indices]
    medians = np.median(sampled, axis=1)
    return pd.DataFrame(
        {
            "cell_area_um2": medians[:, 0],
            "nuc_area_um2": medians[:, 1],
        }
    )


def build_figure_data(
    *,
    n_bootstrap: int = 2_000,
    seed: int = 20260710,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    genome = pd.read_csv(GENOME_PATH, low_memory=False)
    genome_summary = pd.read_csv(GENOME_SUMMARY_PATH, low_memory=False)
    cells = pd.read_csv(CELL_PATH, low_memory=False)

    genome_species = set(genome["species"])
    cell_species = set(cells["species"])
    figure_species = sorted(cell_species)
    if not genome_species.issubset(cell_species):
        raise ValueError(
            "Every frozen genome-panel species must have a frozen size panel."
        )
    if len(figure_species) != 21:
        raise ValueError(f"Expected 21 frozen cell-panel species, found {len(figure_species)}.")
    if not cells.groupby("species").size().eq(50).all():
        raise ValueError("Size panel must contain exactly 50 vetted cells per species.")
    if not genome["review_status"].eq("reviewed_keep").all():
        raise ValueError("Genome panel contains a row that is not reviewed keep.")

    genome_points = genome_summary.set_index("species")
    anchor_values = genome_points["relative_iod_anchor"].dropna().unique()
    if len(anchor_values) != 1:
        raise ValueError("Genome summary does not contain one fixed relative-IOD anchor.")
    anchor = float(anchor_values[0])

    summary_rows: list[dict[str, Any]] = []
    draw_frames: list[pd.DataFrame] = []
    for species_index, species in enumerate(figure_species):
        genome_group = genome.loc[genome["species"].eq(species)]
        cell_group = cells.loc[cells["species"].eq(species)]
        size_draws = bootstrap_paired_size_medians(
            cell_group,
            n_bootstrap=n_bootstrap,
            seed=seed + species_index * 2017 + 1,
        )
        metric_draws = {
            "nucleus_area_um2": size_draws["nuc_area_um2"].to_numpy(),
            "cell_area_um2": size_draws["cell_area_um2"].to_numpy(),
        }
        if species in genome_species:
            metric_draws["relative_iod_index"] = bootstrap_equal_image_iod(
                genome_group,
                n_bootstrap=n_bootstrap,
                seed=seed + species_index * 2017,
            ) / anchor
            genome_point = float(genome_points.loc[species, "relative_iod_index"])
            genome_status = "frozen_primary_common_support"
        else:
            genome_point = np.nan
            genome_status = "limited_overlap_no_frozen_primary_iod"
        points = {
            "relative_iod_index": genome_point,
            "nucleus_area_um2": float(cell_group["nuc_area_um2"].median()),
            "cell_area_um2": float(cell_group["cell_area_um2"].median()),
        }
        row: dict[str, Any] = {
            "species": species,
            "genome_panel_status": genome_status,
            "n_genome_nuclei": int(len(genome_group)),
            "n_genome_images": int(genome_group["filename"].nunique()),
            "n_size_cells": int(len(cell_group)),
            "n_size_images": int(cell_group["filename"].nunique()),
        }
        for metric in METRIC_ORDER:
            row[metric] = points[metric]
            if metric not in metric_draws:
                row[f"{metric}_ci_low"] = np.nan
                row[f"{metric}_ci_high"] = np.nan
                continue
            low, high = np.quantile(metric_draws[metric], [0.025, 0.975])
            row[f"{metric}_ci_low"] = float(low)
            row[f"{metric}_ci_high"] = float(high)
            draw_frames.append(
                pd.DataFrame(
                    {
                        "species": species,
                        "metric": metric,
                        "bootstrap_replicate": np.arange(1, n_bootstrap + 1),
                        "value": metric_draws[metric],
                    }
                )
            )
        summary_rows.append(row)

    summary = pd.DataFrame(summary_rows)
    draws = pd.concat(draw_frames, ignore_index=True)
    comparison_specs = [
        (
            "relative_iod_vs_nucleus_area",
            "relative_iod_index",
            "nucleus_area_um2",
        ),
        ("relative_iod_vs_cell_area", "relative_iod_index", "cell_area_um2"),
        ("nucleus_area_vs_cell_area", "nucleus_area_um2", "cell_area_um2"),
    ]
    correlation_rows = []
    for label, left, right in comparison_specs:
        complete = summary[[left, right]].dropna()
        rho, p_value = spearmanr(complete[left], complete[right])
        correlation_rows.append(
            {
                "comparison": label,
                "left_metric": left,
                "right_metric": right,
                "spearman_rho": float(rho),
                "p_value_descriptive_only": float(p_value),
                "n_species": int(len(complete)),
                "phylogenetically_corrected": False,
            }
        )
    return summary, draws, pd.DataFrame(correlation_rows)


def load_measured_tree(species: list[str]) -> Any:
    try:
        from Bio import Phylo
        from Bio.Phylo.BaseTree import Tree
    except ModuleNotFoundError as error:
        raise ModuleNotFoundError(
            "Biopython is required for tree pruning. Run with: uv run --with biopython"
        ) from error

    full_tree = Phylo.read(TREE_PATH, "newick")
    wanted_tips = {value.replace("D. ", "") for value in species}
    available = {terminal.name for terminal in full_tree.get_terminals()}
    missing = sorted(wanted_tips - available)
    if missing:
        raise ValueError(f"Measured species are absent from the time tree: {missing}")
    mrca = full_tree.common_ancestor(
        [terminal for terminal in full_tree.get_terminals() if terminal.name in wanted_tips]
    )
    tree = Tree(root=copy.deepcopy(mrca), rooted=True)
    tree.root.branch_length = 0.0
    for terminal in list(tree.get_terminals()):
        if terminal.name not in wanted_tips:
            tree.prune(terminal)
    observed = {terminal.name for terminal in tree.get_terminals()}
    if observed != wanted_tips:
        raise ValueError("Pruned time tree does not match the measured species set.")
    tip_depths = [tree.distance(tree.root, terminal) for terminal in tree.get_terminals()]
    if max(tip_depths) - min(tip_depths) > 1e-4:
        raise ValueError("Pruned time tree is not ultrametric within tolerance.")
    return tree


def tree_coordinates(tree: Any) -> tuple[dict[Any, float], dict[Any, float], list[Any]]:
    terminals = list(tree.get_terminals())
    depths = tree.depths()
    max_depth = max(float(depths[terminal]) for terminal in terminals)
    ages = {clade: max_depth - float(depth) for clade, depth in depths.items()}
    y_coordinates: dict[Any, float] = {
        terminal: float(len(terminals) - 1 - index)
        for index, terminal in enumerate(terminals)
    }

    def assign_internal(clade: Any) -> float:
        if clade in y_coordinates:
            return y_coordinates[clade]
        child_values = [assign_internal(child) for child in clade.clades]
        y_coordinates[clade] = float((min(child_values) + max(child_values)) / 2.0)
        return y_coordinates[clade]

    assign_internal(tree.root)
    return ages, y_coordinates, terminals


def draw_tree_axis(
    ax: Any,
    tree: Any,
) -> tuple[list[str], dict[str, float], float]:
    ages, y_coordinates, terminals = tree_coordinates(tree)
    for parent in tree.find_clades(order="preorder"):
        if not parent.clades:
            continue
        child_y = [y_coordinates[child] for child in parent.clades]
        ax.plot(
            [ages[parent], ages[parent]],
            [min(child_y), max(child_y)],
            color="#30343B",
            linewidth=1.2,
        )
        for child in parent.clades:
            ax.plot(
                [ages[parent], ages[child]],
                [y_coordinates[child], y_coordinates[child]],
                color="#30343B",
                linewidth=1.2,
            )
    species_order = []
    y_by_species = {}
    max_age = max(ages.values())
    for terminal in terminals:
        species = canonical_species(terminal.name)
        species_order.append(species)
        y_by_species[species] = y_coordinates[terminal]
        ax.text(
            -0.018 * max_age,
            y_coordinates[terminal],
            species,
            va="center",
            ha="left",
            fontsize=8.6,
            color="#25282D",
            fontstyle="italic",
        )
    ax.set_xlim(max_age * 1.025, -max_age * 0.22)
    ax.set_xticks(np.linspace(0.0, max_age, 5))
    ax.set_ylim(-0.8, len(terminals) - 0.2)
    ax.set_yticks([])
    ax.set_xlabel("Time before present (MYA)")
    ax.set_title("Time-calibrated phylogeny", fontsize=12, fontweight="bold")
    ax.grid(axis="x", alpha=0.16, color="#98A1A8")
    for spine in ["top", "right", "left"]:
        ax.spines[spine].set_visible(False)
    return species_order, y_by_species, max_age


def draw_trait_axis(
    ax: Any,
    *,
    metric: str,
    summary: pd.DataFrame,
    draws: pd.DataFrame,
    species_order: list[str],
    y_by_species: dict[str, float],
) -> None:
    style = METRIC_STYLE[metric]
    available_species = []
    distributions = []
    for species in species_order:
        values = draws.loc[
            draws["species"].eq(species) & draws["metric"].eq(metric), "value"
        ].to_numpy()
        if len(values):
            available_species.append(species)
            distributions.append(values)
    positions = [y_by_species[species] for species in available_species]
    violins = ax.violinplot(
        distributions,
        positions=positions,
        vert=False,
        widths=0.72,
        showmeans=False,
        showmedians=False,
        showextrema=False,
        bw_method=0.25,
    )
    for body in violins["bodies"]:
        body.set_facecolor(style["color"])
        body.set_edgecolor("#52616B")
        body.set_linewidth(0.8)
        body.set_alpha(0.44)

    points = summary.set_index("species")[metric]
    all_values = np.concatenate(distributions)
    lower = float(np.quantile(all_values, 0.002))
    upper = float(np.quantile(all_values, 0.998))
    span = max(upper - lower, abs(upper) * 0.05, 1e-6)
    label_offset = span * 0.018
    for species in species_order:
        value = float(points.loc[species])
        y = y_by_species[species]
        if not np.isfinite(value):
            ax.text(
                0.02,
                y,
                "limited overlap; no frozen primary IOD",
                transform=ax.get_yaxis_transform(),
                va="center",
                ha="left",
                fontsize=6.6,
                color="#7A7F84",
                fontstyle="italic",
            )
            continue
        ax.scatter(value, y, s=24, color="#24282D", zorder=4)
        ax.text(
            value + label_offset,
            y,
            format(value, style["format"]),
            va="center",
            ha="left",
            fontsize=7,
            color="#454A50",
        )
    reference = float(points.median())
    ax.axvline(reference, color="#70777E", linestyle="--", linewidth=0.9, alpha=0.65)
    ax.set_xlim(lower - span * 0.06, upper + span * 0.23)
    ax.set_ylim(-0.8, len(species_order) - 0.2)
    ax.set_yticks([])
    ax.set_title(style["title"], fontsize=12, fontweight="bold")
    ax.set_xlabel(style["xlabel"])
    ax.grid(axis="x", alpha=0.16, color="#98A1A8")
    for spine in ["top", "right", "left"]:
        ax.spines[spine].set_visible(False)


def render_figure(
    summary: pd.DataFrame,
    draws: pd.DataFrame,
    correlations: pd.DataFrame,
) -> tuple[list[str], float]:
    tree = load_measured_tree(summary["species"].tolist())
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(
        1,
        4,
        figsize=(21, 12.5),
        gridspec_kw={"width_ratios": [1.55, 1.12, 1.12, 1.12], "wspace": 0.035},
    )
    species_order, y_by_species, max_age = draw_tree_axis(axes[0], tree)
    summary_order = set(summary["species"])
    if set(species_order) != summary_order:
        raise ValueError("Tree order does not contain the exact measured species set.")
    for ax, metric in zip(axes[1:], METRIC_ORDER):
        draw_trait_axis(
            ax,
            metric=metric,
            summary=summary,
            draws=draws,
            species_order=species_order,
            y_by_species=y_by_species,
        )

    correlation_lookup = correlations.set_index("comparison")["spearman_rho"]
    fig.suptitle(
        "Time-calibrated Desmognathus phylogeny with audited bootstrap trait distributions",
        fontsize=16,
        fontweight="bold",
        y=0.978,
    )
    fig.text(
        0.5,
        0.946,
        (
            "Size panel: all 21 species, with 50 manually vetted largest cells and corresponding nuclei each. "
            "Genome panel: 20 common-support species (relative IOD, not pg); D. ochrophaeus is shown without "
            "an IOD violin because it has limited image-quality overlap. No phylogenetic fills."
        ),
        ha="center",
        va="center",
        fontsize=9.5,
        color="#4A4F55",
    )
    fig.text(
        0.5,
        0.925,
        (
            "Descriptive species-level Spearman ρ: "
            f"IOD–nucleus {correlation_lookup['relative_iod_vs_nucleus_area']:.2f}; "
            f"IOD–cell {correlation_lookup['relative_iod_vs_cell_area']:.2f}; "
            f"nucleus–cell {correlation_lookup['nucleus_area_vs_cell_area']:.2f}. "
            "These correlations are not phylogenetically corrected."
        ),
        ha="center",
        va="center",
        fontsize=9.2,
        color="#4A4F55",
    )
    legend_handles = [
        Patch(
            facecolor="#8EABC0",
            edgecolor="#52616B",
            alpha=0.44,
            label="Bootstrap estimator distribution",
        ),
        Line2D(
            [],
            [],
            marker="o",
            linestyle="",
            color="#24282D",
            label="Current point estimate",
        ),
        Line2D(
            [],
            [],
            linestyle="--",
            color="#70777E",
            label="Across-species median",
        ),
    ]
    fig.legend(
        handles=legend_handles,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.907),
        ncol=3,
        frameon=False,
        fontsize=8.7,
    )
    fig.subplots_adjust(top=0.875, bottom=0.07, left=0.035, right=0.99)
    fig.savefig(PNG_PATH, dpi=220, bbox_inches="tight", facecolor="white")
    fig.savefig(PDF_PATH, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return species_order, max_age


def build(*, n_bootstrap: int = 2_000, seed: int = 20260710) -> dict[str, Any]:
    summary, draws, correlations = build_figure_data(
        n_bootstrap=n_bootstrap,
        seed=seed,
    )
    species_order, tree_height = render_figure(summary, draws, correlations)
    summary["tree_display_order"] = summary["species"].map(
        {species: index + 1 for index, species in enumerate(species_order)}
    )
    summary = summary.sort_values("tree_display_order").reset_index(drop=True)
    summary.to_csv(SUMMARY_PATH, index=False, float_format="%.8f")
    correlations.to_csv(CORRELATION_PATH, index=False, float_format="%.8f")
    manifest = {
        "analysis": "Measured-only time-tree alignment of audited microscopy traits",
        "n_measured_species": int(len(summary)),
        "n_primary_genome_species": int(summary["relative_iod_index"].notna().sum()),
        "missing_primary_genome_species": sorted(
            summary.loc[summary["relative_iod_index"].isna(), "species"].tolist()
        ),
        "species": species_order,
        "tree_source": str(TREE_PATH.resolve()),
        "tree_source_sha256": sha256_file(TREE_PATH),
        "pruned_tree_height_mya": float(tree_height),
        "genome_source": str(GENOME_PATH.resolve()),
        "genome_source_sha256": sha256_file(GENOME_PATH),
        "size_source": str(CELL_PATH.resolve()),
        "size_source_sha256": sha256_file(CELL_PATH),
        "bootstrap_replicates": int(n_bootstrap),
        "random_seed": int(seed),
        "genome_estimator": (
            "mean of image-specific median nuclear IOD, divided by the fixed "
            "median-species IOD anchor"
        ),
        "size_estimators": (
            "pooled medians of the 50 manually vetted literal-largest cells "
            "and their corresponding nuclei"
        ),
        "bootstrap_scope": (
            "genome: images then nuclei within images; size: paired resampling "
            "of the frozen selected cell-nucleus rows"
        ),
        "phylogenetic_fills_used": False,
        "absolute_genome_size_claimed": False,
        "correlations_are_phylogenetically_corrected": False,
        "summary_csv": str(SUMMARY_PATH.resolve()),
        "summary_sha256": sha256_file(SUMMARY_PATH),
        "correlation_csv": str(CORRELATION_PATH.resolve()),
        "correlation_sha256": sha256_file(CORRELATION_PATH),
        "png": str(PNG_PATH.resolve()),
        "png_sha256": sha256_file(PNG_PATH),
        "pdf": str(PDF_PATH.resolve()),
        "pdf_sha256": sha256_file(PDF_PATH),
    }
    MANIFEST_PATH.write_text(
        json.dumps(manifest, indent=2) + "\n",
        encoding="utf-8",
    )
    return manifest


def main() -> None:
    print(json.dumps(build(), indent=2))


if __name__ == "__main__":
    main()
