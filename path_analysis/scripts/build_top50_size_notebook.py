#!/usr/bin/env python3
"""Build the top-50 size-analysis notebook."""

from __future__ import annotations

from pathlib import Path

import nbformat as nbf


PROJECT_ROOT = Path(__file__).resolve().parents[2]
NOTEBOOK_PATH = PROJECT_ROOT / "path_analysis" / "notebooks" / "top50_size_analysis_figures.ipynb"


def md(text: str) -> nbf.NotebookNode:
    return nbf.v4.new_markdown_cell(text.strip() + "\n")


def code(text: str) -> nbf.NotebookNode:
    return nbf.v4.new_code_cell(text.strip() + "\n")


def build_notebook() -> nbf.NotebookNode:
    nb = nbf.v4.new_notebook()
    nb.metadata = {
        "kernelspec": {"display_name": "Python 3 (Dusky)", "language": "python", "name": "python3"},
        "language_info": {"name": "python", "pygments_lexer": "ipython3"},
    }
    nb.cells = [
        md(
            """
            # Curated Top-50 Cell, Nucleus, Genome, and TE Analysis

            This notebook is the transparent analysis surface for the current
            curated top-50 linked cell+nucleus pairs per species. The selected
            count is an audit invariant, not a biological result, so this
            notebook focuses on the measurements that can vary: cell area,
            nucleus area, nucleus darkness/IOD quality, genome-size estimates,
            within-species spread, QC sensitivity, TE predictors, and
            phylogenetic path-model evidence.

            Interpretation boundary: cell and nucleus areas are current
            morphology estimates from curated linked pairs. Genome size is the
            current IOD-calibrated estimate with image-QC primary rows and
            sensitivity sidecars. Any `genome -> nucleus -> cell` causal model
            remains sensitivity-labeled until genome size is supported by an
            independent non-IOD assay.
            """
        ),
        md("## Setup And Source Manifest"),
        code(
            r'''
from pathlib import Path
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from IPython.display import display, HTML, IFrame, FileLink

pd.set_option("display.max_columns", 140)
pd.set_option("display.max_rows", 220)
pd.set_option("display.width", 180)


def find_project_root(start=None) -> Path:
    current = (start or Path.cwd()).resolve()
    for candidate in [current, *current.parents]:
        if (candidate / "paths.yaml").exists():
            return candidate
    raise FileNotFoundError("Could not find Desmognathus_TE project root containing paths.yaml")


PROJECT_ROOT = find_project_root()
PATH_ANALYSIS = PROJECT_ROOT / "path_analysis"
RESULTS = PATH_ANALYSIS / "results"
DERIVED = PATH_ANALYSIS / "data" / "derived"
EXTERNAL = PATH_ANALYSIS / "data" / "external" / "derived"
CELLPROFILER_ROOT = Path.home() / "Projects" / "cellprofiler_test"
RUN_ROOT = CELLPROFILER_ROOT / "output" / "runs" / "mixed_cellpose_yolo_full_dataset_v1_bgclean"
TOP50_FINAL_DIR = RUN_ROOT / "top50_linked_pair_review" / "final_curated_top50_latest"
VERIFIED_TOP50_DIR = RUN_ROOT / "verified_species_dataset_top50_latest"
FIG_DIR = RESULTS / "top50_notebook_figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)
for stale_path in FIG_DIR.glob("*"):
    if stale_path.is_file() and stale_path.suffix.lower() in {".png", ".csv"}:
        stale_path.unlink()
REVIEW_APP_URL = "http://127.0.0.1:8774/"

paths = {
    "final_curated_top50_linked_pairs": (TOP50_FINAL_DIR / "final_curated_top50_linked_pairs.csv", True),
    "verified_top50_selected_pairs_with_iod_qc": (VERIFIED_TOP50_DIR / "selected_high_quality_linked_pairs_with_iod_qc.csv.gz", True),
    "verified_top50_species_estimates": (VERIFIED_TOP50_DIR / "species_estimates_verified.csv", True),
    "imported_genome_species_results": (EXTERNAL / "cellprofiler_final_species_results.csv", True),
    "imported_morphology_summary": (EXTERNAL / "cellprofiler_species_morphology_summary.csv", True),
    "genome_sensitivity": (EXTERNAL / "cellprofiler_genome_sensitivity.csv", True),
    "image_iod_quality_summary": (EXTERNAL / "cellprofiler_image_iod_quality_summary.csv", True),
    "path_input_master": (DERIVED / "path_input_master.csv", True),
    "analysis_panel_summary": (DERIVED / "analysis_panel_summary.csv", True),
    "top50_size_spread_summary": (RESULTS / "top50_size_estimate_spread_summary.csv", True),
    "top50_spearman_correlations": (RESULTS / "top50_size_spearman_correlations.csv", False),
    "top50_model_ranking_summary": (RESULTS / "top50_model_ranking_summary.csv", True),
    "genome_morphology_model_ranking": (RESULTS / "genome_morphology_model_ranking.csv", True),
    "te_genome_morphology_model_ranking": (RESULTS / "te_genome_morphology_model_ranking.csv", True),
}

source_table = pd.DataFrame(
    [
        {
            "name": name,
            "path": str(path),
            "required": required,
            "exists": path.exists(),
            "size_mb": path.stat().st_size / 1_000_000 if path.exists() else np.nan,
        }
        for name, (path, required) in paths.items()
    ]
)
display(source_table)
missing_required = source_table[source_table["required"] & ~source_table["exists"]]
assert missing_required.empty, "Required source files are missing: " + ", ".join(missing_required["name"].tolist())
print("Project root:", PROJECT_ROOT)
print("Figure output directory:", FIG_DIR)
'''
        ),
        md("## Load Exact Inputs"),
        code(
            r'''
top50_pairs = pd.read_csv(paths["final_curated_top50_linked_pairs"][0], keep_default_na=False, low_memory=False)
verified_pairs = pd.read_csv(paths["verified_top50_selected_pairs_with_iod_qc"][0], low_memory=False)
verified_species = pd.read_csv(paths["verified_top50_species_estimates"][0], low_memory=False)
genome = pd.read_csv(paths["imported_genome_species_results"][0], low_memory=False)
morphology = pd.read_csv(paths["imported_morphology_summary"][0], low_memory=False)
genome_sensitivity = pd.read_csv(paths["genome_sensitivity"][0], low_memory=False)
image_qc = pd.read_csv(paths["image_iod_quality_summary"][0], low_memory=False)
path_input = pd.read_csv(paths["path_input_master"][0], low_memory=False)
panel_summary = pd.read_csv(paths["analysis_panel_summary"][0], low_memory=False)
size_spread = pd.read_csv(paths["top50_size_spread_summary"][0], low_memory=False)
model_summary = pd.read_csv(paths["top50_model_ranking_summary"][0], low_memory=False)

numeric_candidates = [
    "cell_area_um2",
    "nuc_area_um2",
    "nucleus_area_um2",
    "nuc_iod",
    "nuc_mean_od",
    "cell_mean_od",
    "nc_area_ratio",
    "cytoplasm_area_um2",
    "quality_rank_score",
    "pair_quality_probability",
    "keep_probability",
    "quality_score",
    "iod_final_quality_score",
    "image_iod_quality_score",
    "raw_focus_lap_var",
    "raw_focus_grad_mean",
    "raw_focus_intensity_iqr",
    "nucleus_darkness_signal",
    "clarity_signal",
    "large_dark_signal",
    "large_clear_signal",
]
for df in [top50_pairs, verified_pairs]:
    for col in numeric_candidates:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")


def species_short_value(value) -> str:
    text = str(value)
    return text.replace("D. ", "").replace("Desmognathus ", "")


def species_key(series: pd.Series) -> pd.Series:
    return (
        series.astype(str)
        .str.replace("D. ", "", regex=False)
        .str.replace("Desmognathus ", "", regex=False)
        .str.strip()
        .str.lower()
    )


def savefig(fig, filename: str) -> Path:
    path = FIG_DIR / filename
    fig.tight_layout()
    fig.savefig(path, dpi=190, bbox_inches="tight")
    return path


inventory = pd.DataFrame(
    {
        "table": [
            "top50_pairs",
            "verified_pairs",
            "verified_species",
            "genome",
            "morphology",
            "genome_sensitivity",
            "image_qc",
            "path_input",
            "panel_summary",
            "size_spread",
            "model_summary",
        ],
        "rows": [
            len(top50_pairs),
            len(verified_pairs),
            len(verified_species),
            len(genome),
            len(morphology),
            len(genome_sensitivity),
            len(image_qc),
            len(path_input),
            len(panel_summary),
            len(size_spread),
            len(model_summary),
        ],
        "species": [
            top50_pairs["species"].nunique(),
            verified_pairs["species"].nunique(),
            verified_species["species"].nunique(),
            genome["species"].nunique(),
            morphology["species"].nunique(),
            genome_sensitivity["species"].nunique(),
            image_qc["species"].nunique(),
            path_input["species"].nunique(),
            np.nan,
            size_spread["species"].nunique(),
            np.nan,
        ],
    }
)
display(inventory)
'''
        ),
        md("## Selection Audit"),
        code(
            r'''
counts = top50_pairs.groupby("species", sort=True).size().rename("n_top50_pairs").reset_index()
decisions = top50_pairs.get("grid_effective_decision", pd.Series("", index=top50_pairs.index)).astype(str).str.lower()
reject_leakage = int(decisions.isin({"discard", "nucleus_only", "reject"}).sum())
audit = pd.DataFrame(
    {
        "check": [
            "species_with_selected_pairs",
            "selected_pair_rows",
            "minimum_pairs_per_species",
            "maximum_pairs_per_species",
            "explicit_reject_rows_in_final_table",
            "verified_pair_rows",
            "verified_species",
        ],
        "value": [
            top50_pairs["species"].nunique(),
            len(top50_pairs),
            int(counts["n_top50_pairs"].min()),
            int(counts["n_top50_pairs"].max()),
            reject_leakage,
            len(verified_pairs),
            verified_pairs["species"].nunique(),
        ],
    }
)
display(audit)
display(counts)
assert counts["n_top50_pairs"].eq(50).all(), "Every species should have exactly 50 selected linked pairs"
assert reject_leakage == 0, "Rejected rows leaked into final selected-pair table"
'''
        ),
        md(
            """
            ## Source Rows Feeding The Estimates

            These are the row-level selected linked pairs. The important
            measurement columns are shown directly so any figure can be traced
            back to the underlying cells and nuclei.
            """
        ),
        code(
            r'''
pair_cols = [
    "species",
    "grid_species_rank",
    "grid_selection_tier",
    "grid_effective_decision",
    "filename",
    "tile_name",
    "review_key",
    "cell_area_um2",
    "nuc_area_um2",
    "nuc_iod",
    "nuc_mean_od",
    "cell_mean_od",
    "nc_area_ratio",
    "quality_rank_score",
    "pair_quality_probability",
    "keep_probability",
    "nucleus_darkness_signal",
    "clarity_signal",
    "large_clear_signal",
]
pair_cols = [c for c in pair_cols if c in top50_pairs.columns]
display(top50_pairs[pair_cols].head(30))

qc_cols = [
    "species",
    "filename",
    "review_key",
    "cell_area_um2",
    "nuc_area_um2",
    "nuc_iod",
    "nuc_mean_od",
    "iod_final_quality_score",
    "image_iod_quality_score",
    "image_iod_qc_status",
    "image_iod_qc_pass",
    "selection_tier",
    "quality_score",
]
qc_cols = [c for c in qc_cols if c in verified_pairs.columns]
display(verified_pairs[qc_cols].head(30))
'''
        ),
        md("## Pair-Level Descriptive Statistics"),
        code(
            r'''
metric_labels = {
    "cell_area_um2": "Cell area (um2)",
    "nuc_area_um2": "Nucleus area (um2)",
    "nc_area_ratio": "Nucleus:cell area ratio",
    "cytoplasm_area_um2": "Cytoplasm area (um2)",
    "nuc_iod": "Nucleus IOD",
    "nuc_mean_od": "Nucleus mean OD / darkness",
    "cell_mean_od": "Cell mean OD",
    "iod_final_quality_score": "Final IOD quality score",
    "image_iod_quality_score": "Image IOD quality score",
    "quality_score": "Pair quality score",
    "quality_rank_score": "Grid quality rank score",
    "keep_probability": "Keep probability",
    "pair_quality_probability": "Pair quality probability",
    "raw_focus_lap_var": "Focus Laplacian variance",
    "raw_focus_grad_mean": "Focus gradient mean",
    "raw_focus_intensity_iqr": "Raw intensity IQR",
    "nucleus_darkness_signal": "Nucleus darkness signal",
    "clarity_signal": "Clarity signal",
    "large_clear_signal": "Large clear signal",
}
pair_metrics = [
    col
    for col in metric_labels
    if col in verified_pairs.columns and pd.to_numeric(verified_pairs[col], errors="coerce").notna().sum() > 0
]


def summarize_metric_by_species(df: pd.DataFrame, metric: str) -> pd.DataFrame:
    rows = []
    for species, group in df.groupby("species", sort=True):
        values = pd.to_numeric(group[metric], errors="coerce").dropna()
        if values.empty:
            continue
        mean = values.mean()
        sd = values.std(ddof=1)
        rows.append(
            {
                "species": species,
                "metric": metric,
                "metric_label": metric_labels.get(metric, metric),
                "n": int(values.size),
                "mean": mean,
                "median": values.median(),
                "sd": sd,
                "cv_pct": (sd / mean * 100.0) if mean else np.nan,
                "q1": values.quantile(0.25),
                "q3": values.quantile(0.75),
                "iqr": values.quantile(0.75) - values.quantile(0.25),
                "min": values.min(),
                "max": values.max(),
            }
        )
    return pd.DataFrame(rows)


pair_descriptive_stats = pd.concat(
    [summarize_metric_by_species(verified_pairs, metric) for metric in pair_metrics],
    ignore_index=True,
)
pair_descriptive_stats.to_csv(FIG_DIR / "top50_pair_level_descriptive_stats.csv", index=False)

key_metrics = [
    "cell_area_um2",
    "nuc_area_um2",
    "nc_area_ratio",
    "nuc_iod",
    "nuc_mean_od",
    "iod_final_quality_score",
    "quality_score",
]
key_metrics = [m for m in key_metrics if m in pair_metrics]
display(pair_descriptive_stats[pair_descriptive_stats["metric"].isin(key_metrics)])
print("Saved:", FIG_DIR / "top50_pair_level_descriptive_stats.csv")
'''
        ),
        md("## Figure: Raw Top-50 Cell, Nucleus, And Ratio Distributions"),
        code(
            r'''
def plot_species_distribution(df, metric, ax, title, xlabel, color="#345995", order_by=None):
    sub = df[["species", metric]].copy()
    sub[metric] = pd.to_numeric(sub[metric], errors="coerce")
    sub = sub.dropna()
    if sub.empty:
        ax.set_axis_off()
        ax.set_title(title + " (no data)")
        return
    order_metric = order_by or metric
    order_frame = df[["species", order_metric]].copy()
    order_frame[order_metric] = pd.to_numeric(order_frame[order_metric], errors="coerce")
    order = order_frame.groupby("species")[order_metric].median().dropna().sort_values().index.tolist()
    order = [sp for sp in order if sp in set(sub["species"])]
    labels = [species_short_value(sp) for sp in order]
    data = [sub.loc[sub["species"] == sp, metric].to_numpy(float) for sp in order]
    positions = np.arange(1, len(order) + 1)
    bp = ax.boxplot(data, vert=False, positions=positions, widths=0.58, patch_artist=True, showfliers=False)
    for patch in bp["boxes"]:
        patch.set_facecolor(color)
        patch.set_alpha(0.22)
        patch.set_edgecolor("#333333")
    for line in bp["medians"]:
        line.set_color("#111111")
        line.set_linewidth(1.4)
    rng = np.random.default_rng(41)
    for y, values in zip(positions, data):
        jitter = rng.normal(0, 0.055, size=len(values))
        ax.scatter(values, y + jitter, s=8, alpha=0.34, color=color, linewidths=0)
    ax.set_yticks(positions)
    ax.set_yticklabels(labels, fontsize=8)
    ax.set_xlabel(xlabel)
    ax.set_title(title)
    ax.grid(axis="x", alpha=0.22)


fig, axes = plt.subplots(2, 2, figsize=(15, 13), sharey=False)
plot_species_distribution(verified_pairs, "cell_area_um2", axes[0, 0], "Cell Area By Species", "cell area (um2)", "#496A81")
plot_species_distribution(verified_pairs, "nuc_area_um2", axes[0, 1], "Nucleus Area By Species", "nucleus area (um2)", "#7C6A9A")
plot_species_distribution(verified_pairs, "nc_area_ratio", axes[1, 0], "Nucleus:Cell Area Ratio By Species", "N:C area ratio", "#7A7F3F")
plot_species_distribution(verified_pairs, "cytoplasm_area_um2", axes[1, 1], "Cytoplasm Area By Species", "cytoplasm area (um2)", "#9B5C49")
savefig(fig, "top50_raw_cell_nucleus_ratio_distributions.png")
plt.show()
'''
        ),
        md("## Figure: Nuclear Darkness, Image Quality, And Cell Size"),
        code(
            r'''
quality_metrics = [
    ("nuc_mean_od", "Nucleus Mean OD / Darkness", "nucleus mean OD", "#6A4C93"),
    ("nuc_iod", "Nucleus Integrated Optical Density", "nucleus IOD", "#4D6A6D"),
    ("iod_final_quality_score", "Final IOD Quality Score", "quality score", "#8A6F3D"),
    ("quality_score", "Pair Quality Score", "quality score", "#4F6D44"),
]
quality_metrics = [m for m in quality_metrics if m[0] in verified_pairs.columns]
fig, axes = plt.subplots(2, 2, figsize=(15, 13), sharey=False)
for ax, (metric, title, xlabel, color) in zip(axes.flat, quality_metrics):
    plot_species_distribution(verified_pairs, metric, ax, title, xlabel, color, order_by="nuc_mean_od")
for ax in axes.flat[len(quality_metrics):]:
    ax.set_axis_off()
savefig(fig, "top50_nuclear_darkness_and_quality_distributions.png")
plt.show()

species_medians = (
    verified_pairs.groupby("species", sort=True)
    .agg(
        median_cell_area_um2=("cell_area_um2", "median"),
        median_nucleus_area_um2=("nuc_area_um2", "median"),
        median_nc_ratio=("nc_area_ratio", "median"),
        median_nucleus_iod=("nuc_iod", "median"),
        median_nucleus_mean_od=("nuc_mean_od", "median"),
        median_iod_quality=("iod_final_quality_score", "median"),
        median_pair_quality=("quality_score", "median"),
    )
    .reset_index()
)
display(species_medians.sort_values("median_cell_area_um2", ascending=False))

fig, axes = plt.subplots(1, 2, figsize=(14, 6))
scatter_pairs = [
    ("median_cell_area_um2", "median_nucleus_mean_od", "Median cell area vs nuclear darkness", "median cell area (um2)", "median nucleus mean OD"),
    ("median_cell_area_um2", "median_iod_quality", "Median cell area vs IOD quality", "median cell area (um2)", "median IOD quality score"),
]
for ax, (x, y, title, xlabel, ylabel) in zip(axes, scatter_pairs):
    sub = species_medians[[x, y, "species"]].dropna()
    ax.scatter(sub[x], sub[y], s=55, color="#354F52", alpha=0.82)
    for _, row in sub.iterrows():
        ax.annotate(species_short_value(row["species"]), (row[x], row[y]), xytext=(4, 2), textcoords="offset points", fontsize=7)
    rho = sub[[x, y]].corr(method="spearman").iloc[0, 1] if len(sub) >= 3 else np.nan
    ax.set_title(f"{title}\nSpearman rho={rho:.2f}" if np.isfinite(rho) else title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.grid(alpha=0.22)
savefig(fig, "top50_largest_darkest_species_medians.png")
plt.show()
'''
        ),
        md("## Species-Level Estimates, Uncertainty, And Extremes"),
        code(
            r'''
estimate_stats = size_spread.copy()
estimate_stats["genome_ci_rel_pct"] = estimate_stats["genome_ci_width_pg"] / estimate_stats["genome_pg"] * 100.0
estimate_stats["genome_iqr_rel_pct"] = estimate_stats["genome_iqr_pg"] / estimate_stats["genome_pg"] * 100.0
estimate_stats["cell_area_ci_rel_pct"] = estimate_stats["cell_area_ci_width_um2"] / estimate_stats["cell_area_um2"] * 100.0
estimate_stats["cell_area_iqr_rel_pct"] = estimate_stats["cell_area_iqr_um2"] / estimate_stats["cell_area_um2"] * 100.0
estimate_stats["nucleus_area_ci_rel_pct"] = estimate_stats["nucleus_area_ci_width_um2"] / estimate_stats["nucleus_area_um2"] * 100.0
estimate_stats["nucleus_area_iqr_rel_pct"] = estimate_stats["nucleus_area_iqr_um2"] / estimate_stats["nucleus_area_um2"] * 100.0

display(
    estimate_stats[
        [
            "species",
            "n_top50_pairs",
            "genome_pg",
            "genome_ci_low_pg",
            "genome_ci_high_pg",
            "genome_ci_rel_pct",
            "genome_support_tier",
            "genome_support_warnings",
            "cell_area_um2",
            "cell_area_ci_low_um2",
            "cell_area_ci_high_um2",
            "cell_area_ci_rel_pct",
            "nucleus_area_um2",
            "nucleus_area_ci_low_um2",
            "nucleus_area_ci_high_um2",
            "nucleus_area_ci_rel_pct",
            "nc_ratio",
        ]
    ].sort_values("genome_pg")
)


def extremes(df, metric, label, n=8):
    cols = ["species", metric]
    low = df[cols].dropna().nsmallest(n, metric).assign(which="lowest", metric_label=label)
    high = df[cols].dropna().nlargest(n, metric).assign(which="highest", metric_label=label)
    return pd.concat([low, high], ignore_index=True)


extreme_table = pd.concat(
    [
        extremes(estimate_stats, "genome_pg", "Genome size"),
        extremes(estimate_stats, "cell_area_um2", "Cell area"),
        extremes(estimate_stats, "nucleus_area_um2", "Nucleus area"),
        extremes(estimate_stats, "nc_ratio", "N:C ratio"),
        extremes(estimate_stats, "genome_ci_rel_pct", "Genome relative CI width"),
        extremes(estimate_stats, "cell_area_iqr_rel_pct", "Cell area relative IQR"),
        extremes(estimate_stats, "nucleus_area_iqr_rel_pct", "Nucleus area relative IQR"),
    ],
    ignore_index=True,
)
display(extreme_table)
estimate_stats.to_csv(FIG_DIR / "top50_species_estimate_uncertainty_stats.csv", index=False)
print("Saved:", FIG_DIR / "top50_species_estimate_uncertainty_stats.csv")
'''
        ),
        md("## Figure: Species Estimates With Bootstrap Confidence Intervals"),
        code(
            r'''
def plot_interval(ax, df, value, low, high, title, xlabel, color):
    sub = df[["species", value, low, high]].dropna().sort_values(value).reset_index(drop=True)
    y = np.arange(len(sub))
    x = sub[value].to_numpy(float)
    xlow = sub[low].to_numpy(float)
    xhigh = sub[high].to_numpy(float)
    xerr = np.vstack([x - xlow, xhigh - x])
    ax.errorbar(x, y, xerr=xerr, fmt="o", color=color, ecolor="#B8B2A7", elinewidth=2.2, capsize=3, markersize=5.5)
    ax.set_yticks(y)
    ax.set_yticklabels([species_short_value(sp) for sp in sub["species"]], fontsize=8)
    ax.set_xlabel(xlabel)
    ax.set_title(title)
    ax.grid(axis="x", alpha=0.24)


fig, axes = plt.subplots(1, 3, figsize=(18, 8), sharey=False)
plot_interval(axes[0], estimate_stats, "genome_pg", "genome_ci_low_pg", "genome_ci_high_pg", "Genome Size", "pg", "#3A5A40")
plot_interval(axes[1], estimate_stats, "cell_area_um2", "cell_area_ci_low_um2", "cell_area_ci_high_um2", "Cell Area", "um2", "#315C7A")
plot_interval(axes[2], estimate_stats, "nucleus_area_um2", "nucleus_area_ci_low_um2", "nucleus_area_ci_high_um2", "Nucleus Area", "um2", "#7A4E72")
savefig(fig, "top50_species_estimates_with_ci.png")
plt.show()
'''
        ),
        md("## Figure: Relative Uncertainty And Within-Species Spread"),
        code(
            r'''
heat_cols = [
    ("genome_ci_rel_pct", "Genome CI %"),
    ("genome_iqr_rel_pct", "Genome IQR %"),
    ("cell_area_ci_rel_pct", "Cell CI %"),
    ("cell_area_iqr_rel_pct", "Cell IQR %"),
    ("nucleus_area_ci_rel_pct", "Nucleus CI %"),
    ("nucleus_area_iqr_rel_pct", "Nucleus IQR %"),
]
heat_df = estimate_stats.set_index("species")[[c for c, _ in heat_cols]].sort_values("genome_ci_rel_pct")
matrix = heat_df.to_numpy(float)
fig, ax = plt.subplots(figsize=(11, 8))
im = ax.imshow(matrix, aspect="auto", cmap="YlOrRd")
ax.set_yticks(np.arange(len(heat_df.index)))
ax.set_yticklabels([species_short_value(sp) for sp in heat_df.index], fontsize=8)
ax.set_xticks(np.arange(len(heat_cols)))
ax.set_xticklabels([label for _, label in heat_cols], rotation=35, ha="right")
ax.set_title("Relative uncertainty and within-species spread")
for i in range(matrix.shape[0]):
    for j in range(matrix.shape[1]):
        value = matrix[i, j]
        if np.isfinite(value):
            ax.text(j, i, f"{value:.1f}", ha="center", va="center", fontsize=6.5, color="#1d1d1d")
cbar = fig.colorbar(im, ax=ax)
cbar.set_label("percent of species estimate")
savefig(fig, "top50_relative_uncertainty_heatmap.png")
plt.show()
'''
        ),
        md("## Image QC And Genome Sensitivity"),
        code(
            r'''
image_qc_display = image_qc.sort_values(["species", "image_iod_qc_status", "filename"])
display(image_qc_display)

status_counts = (
    image_qc.groupby(["species", "image_iod_qc_status"], dropna=False)
    .size()
    .rename("n_images")
    .reset_index()
)
status_wide = status_counts.pivot(index="species", columns="image_iod_qc_status", values="n_images").fillna(0)
status_wide = status_wide.loc[sorted(status_wide.index, key=species_short_value)]
display(status_wide)

fig, ax = plt.subplots(figsize=(10, max(5, 0.32 * len(status_wide))))
left = np.zeros(len(status_wide))
colors = {
    "pass": "#4F772D",
    "limited": "#C49A3A",
    "fail": "#9E2A2B",
    "nan": "#AAAAAA",
}
for col in status_wide.columns:
    values = status_wide[col].to_numpy(float)
    label = "missing" if pd.isna(col) else str(col)
    ax.barh(
        [species_short_value(sp) for sp in status_wide.index],
        values,
        left=left,
        label=label,
        color=colors.get(label, "#6C757D"),
        alpha=0.88,
    )
    left += values
ax.set_xlabel("images")
ax.set_title("Image-level IOD QC status by species")
ax.legend(loc="lower right")
ax.grid(axis="x", alpha=0.22)
savefig(fig, "top50_image_iod_qc_status_by_species.png")
plt.show()

sensitivity_view = genome_sensitivity[
    [
        "species",
        "subset_name",
        "n_selected_pairs",
        "n_selected_images",
        "n_selected_specimens",
        "effective_n",
        "support_label",
        "estimated_genome_pg_estimate",
        "estimated_genome_pg_ci_low",
        "estimated_genome_pg_ci_high",
        "pct_shift_vs_all_selected",
    ]
].sort_values(["species", "subset_name"])
display(sensitivity_view)

shift = genome_sensitivity[genome_sensitivity["subset_name"].ne("all_selected")].copy()
shift["species_short"] = shift["species"].map(species_short_value)
shift = shift.sort_values("pct_shift_vs_all_selected")
fig, ax = plt.subplots(figsize=(11, max(5.5, 0.2 * len(shift))))
colors = shift["subset_name"].map({"image_qc_pass": "#315C7A", "high_iod_qc": "#7A4E72"}).fillna("#666666")
ax.barh(shift["species_short"] + " - " + shift["subset_name"], shift["pct_shift_vs_all_selected"], color=colors)
ax.axvline(0, color="#222222", linewidth=0.9)
ax.set_xlabel("percent shift vs all selected")
ax.set_title("Genome estimate sensitivity to image/nuclear-IOD QC filtering")
ax.grid(axis="x", alpha=0.24)
savefig(fig, "top50_genome_sensitivity_qc_shift.png")
plt.show()
'''
        ),
        md("## Cell, Nucleus, And Genome Relationships"),
        code(
            r'''
def scatter_labeled(ax, df, x, y, title, xlabel, ylabel, color="#315C7A"):
    sub = df[["species", x, y]].copy()
    sub[x] = pd.to_numeric(sub[x], errors="coerce")
    sub[y] = pd.to_numeric(sub[y], errors="coerce")
    sub = sub.dropna()
    ax.scatter(sub[x], sub[y], s=55, color=color, alpha=0.82)
    for _, row in sub.iterrows():
        ax.annotate(species_short_value(row["species"]), (row[x], row[y]), xytext=(4, 2), textcoords="offset points", fontsize=7)
    rho = sub[[x, y]].corr(method="spearman").iloc[0, 1] if len(sub) >= 3 else np.nan
    if len(sub) >= 3 and sub[x].nunique() > 1:
        coeff = np.polyfit(sub[x], sub[y], 1)
        xs = np.linspace(sub[x].min(), sub[x].max(), 100)
        ax.plot(xs, coeff[0] * xs + coeff[1], color="#222222", alpha=0.42, linewidth=1.1)
    ax.set_title(f"{title}\nSpearman rho={rho:.2f}, n={len(sub)}" if np.isfinite(rho) else f"{title}\nn={len(sub)}")
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.grid(alpha=0.22)


relationships = [
    ("genome_pg", "nucleus_area_um2", "Genome vs nucleus area", "genome size (pg)", "nucleus area (um2)"),
    ("genome_pg", "cell_area_um2", "Genome vs cell area", "genome size (pg)", "cell area (um2)"),
    ("nucleus_area_um2", "cell_area_um2", "Nucleus area vs cell area", "nucleus area (um2)", "cell area (um2)"),
    ("genome_pg", "nc_ratio", "Genome vs N:C ratio", "genome size (pg)", "N:C area ratio"),
    ("cell_area_um2", "nc_ratio", "Cell area vs N:C ratio", "cell area (um2)", "N:C area ratio"),
    ("nucleus_area_um2", "nc_ratio", "Nucleus area vs N:C ratio", "nucleus area (um2)", "N:C area ratio"),
]
fig, axes = plt.subplots(2, 3, figsize=(18, 10))
for ax, args in zip(axes.flat, relationships):
    scatter_labeled(ax, estimate_stats, *args)
savefig(fig, "top50_size_relationships_detailed.png")
plt.show()

size_corr_cols = ["genome_pg", "cell_area_um2", "nucleus_area_um2", "nc_ratio"]
size_corr = estimate_stats[size_corr_cols].corr(method="spearman")
display(size_corr)
'''
        ),
        md("## TE, Organismal, Genome, Cell, And Nucleus Correlations"),
        code(
            r'''
path_view = path_input.copy()
path_view["species_key"] = species_key(path_view["species"])
estimate_view = estimate_stats.copy()
estimate_view["species_key"] = species_key(estimate_view["species"])
merged = path_view.merge(
    estimate_view[
        [
            "species_key",
            "species",
            "genome_pg",
            "cell_area_um2",
            "nucleus_area_um2",
            "nc_ratio",
            "genome_ci_rel_pct",
            "cell_area_iqr_rel_pct",
            "nucleus_area_iqr_rel_pct",
        ]
    ],
    on="species_key",
    how="inner",
    suffixes=("_path", "_top50"),
)
merged["species"] = merged["species_top50"]

candidate_vars = [
    ("genome_pg", "Genome pg"),
    ("cell_area_um2", "Cell area"),
    ("nucleus_area_um2", "Nucleus area"),
    ("nc_ratio", "N:C ratio"),
    ("order_ltr_te", "LTR fraction"),
    ("order_line_te", "LINE fraction"),
    ("order_tir_te", "TIR fraction"),
    ("ltr_line_logratio", "LTR:LINE log ratio"),
    ("retro_dna_logratio", "Retro:DNA log ratio"),
    ("order_pielou_te", "TE order evenness"),
    ("weighted_te_divergence_p90", "TE divergence p90"),
    ("weighted_te_deletions_p90", "TE deletion p90"),
    ("ectopic_log10_mean_ratio", "Ectopic log10 mean ratio"),
    ("body_size_proxy_mm", "Body size proxy"),
    ("aquaticity_index", "Aquaticity"),
    ("ltr_history_age_central_high_conf_mya", "LTR age high-conf"),
]
available_vars = [(col, label) for col, label in candidate_vars if col in merged.columns and pd.to_numeric(merged[col], errors="coerce").notna().sum() >= 5]
corr_input = merged[["species"] + [col for col, _ in available_vars]].copy()
for col, _ in available_vars:
    corr_input[col] = pd.to_numeric(corr_input[col], errors="coerce")
display(corr_input.sort_values("species"))

corr = corr_input[[col for col, _ in available_vars]].corr(method="spearman", min_periods=5)
labels = [label for _, label in available_vars]
fig, ax = plt.subplots(figsize=(12, 10))
im = ax.imshow(corr.to_numpy(float), vmin=-1, vmax=1, cmap="RdBu_r")
ax.set_xticks(np.arange(len(labels)))
ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=8)
ax.set_yticks(np.arange(len(labels)))
ax.set_yticklabels(labels, fontsize=8)
ax.set_title("Spearman correlations across morphology, genome, TE, and organismal variables")
for i in range(len(labels)):
    for j in range(len(labels)):
        value = corr.iloc[i, j]
        if np.isfinite(value):
            ax.text(j, i, f"{value:.2f}", ha="center", va="center", fontsize=6.5)
cbar = fig.colorbar(im, ax=ax, shrink=0.82)
cbar.set_label("Spearman rho")
savefig(fig, "top50_te_genome_morphology_correlation_heatmap.png")
plt.show()
'''
        ),
        md("## Figure: TE And Organismal Predictors Against Genome Size"),
        code(
            r'''
te_scatter_specs = [
    ("ltr_line_logratio", "genome_pg", "LTR:LINE balance vs genome", "LTR:LINE log ratio", "genome size (pg)", "#5C677D"),
    ("order_pielou_te", "genome_pg", "TE evenness vs genome", "TE order evenness", "genome size (pg)", "#6A994E"),
    ("weighted_te_deletions_p90", "genome_pg", "Deletion proxy vs genome", "weighted TE deletion p90", "genome size (pg)", "#9C6644"),
    ("ectopic_log10_mean_ratio", "genome_pg", "Ectopic ratio vs genome", "ectopic log10 mean ratio", "genome size (pg)", "#8E5572"),
    ("body_size_proxy_mm", "genome_pg", "Body size vs genome", "body size proxy (mm)", "genome size (pg)", "#315C7A"),
    ("ltr_history_age_central_high_conf_mya", "genome_pg", "LTR history age vs genome", "LTR high-conf age (mya)", "genome size (pg)", "#7A4E72"),
]
te_scatter_specs = [spec for spec in te_scatter_specs if spec[0] in merged.columns and pd.to_numeric(merged[spec[0]], errors="coerce").notna().sum() >= 5]
fig, axes = plt.subplots(2, 3, figsize=(18, 10))
for ax, spec in zip(axes.flat, te_scatter_specs):
    scatter_labeled(ax, merged, *spec)
for ax in axes.flat[len(te_scatter_specs):]:
    ax.set_axis_off()
savefig(fig, "top50_te_predictors_vs_genome.png")
plt.show()
'''
        ),
        md("## Panel Coverage And Model Rankings"),
        code(
            r'''
display(panel_summary.sort_values("panel_name"))
display(model_summary.sort_values(["n_species", "run"], ascending=[False, True]))

ranking_frames = []
for path in sorted(RESULTS.glob("*_model_ranking.csv")):
    if path.name == "top50_model_ranking_summary.csv":
        continue
    df = pd.read_csv(path)
    if df.empty:
        continue
    run = path.name.replace("_model_ranking.csv", "")
    top = df.sort_values("delta_CICc", na_position="last").head(1).copy()
    top["run"] = run
    top["n_models"] = len(df)
    ranking_frames.append(top)
ranking_overview = pd.concat(ranking_frames, ignore_index=True) if ranking_frames else pd.DataFrame()
ranking_overview = ranking_overview[
    [c for c in ["run", "n_models", "model", "w", "delta_CICc", "CICc", "p"] if c in ranking_overview.columns]
].sort_values("run")
display(ranking_overview)

plot_overview = ranking_overview.copy()
if "w" in plot_overview.columns:
    plot_overview["w"] = pd.to_numeric(plot_overview["w"], errors="coerce")
    plot_overview = plot_overview.dropna(subset=["w"]).sort_values("w")
    fig, ax = plt.subplots(figsize=(11, max(5, 0.28 * len(plot_overview))))
    ax.barh(plot_overview["run"], plot_overview["w"], color="#4E6E5D")
    ax.set_xlabel("top model Akaike weight")
    ax.set_title("Top model support by analysis run")
    ax.grid(axis="x", alpha=0.22)
    savefig(fig, "top50_top_model_weights_by_run.png")
    plt.show()
'''
        ),
        md("## Figure: Within-Run Model Competition"),
        code(
            r'''
selected_ranking_runs = [
    "genome_morphology",
    "te_genome_morphology",
    "te_genome_primary_mediumplus",
    "te_genome_ectopic_primary_mediumplus",
    "te_genome_ltr_history_primary_mediumplus",
    "te_genome_organismal_primary_mediumplus",
]
available_ranking_runs = [run for run in selected_ranking_runs if (RESULTS / f"{run}_model_ranking.csv").exists()]
n = len(available_ranking_runs)
fig, axes = plt.subplots(n, 1, figsize=(13, max(4, 2.1 * n)))
if n == 1:
    axes = [axes]
for ax, run in zip(axes, available_ranking_runs):
    df = pd.read_csv(RESULTS / f"{run}_model_ranking.csv")
    df["w"] = pd.to_numeric(df.get("w", np.nan), errors="coerce")
    df["delta_CICc"] = pd.to_numeric(df.get("delta_CICc", np.nan), errors="coerce")
    sub = df.sort_values("delta_CICc", na_position="last").head(6).copy()
    sub = sub.iloc[::-1]
    if sub["w"].notna().any():
        ax.barh(sub["model"], sub["w"], color="#4C6A92")
        ax.set_xlabel("Akaike weight")
    else:
        ax.barh(sub["model"], -sub["delta_CICc"], color="#4C6A92")
        ax.set_xlabel("-delta CICc")
    ax.set_title(run)
    ax.grid(axis="x", alpha=0.22)
savefig(fig, "top50_selected_model_competition.png")
plt.show()
'''
        ),
        md("## Figure: Best And Averaged Model Edges"),
        code(
            r'''
edge_label_map = {
    "gs": "genome",
    "ns": "nucleus",
    "cs": "cell",
    "ltr_balance": "LTR balance",
    "te_evenness": "TE evenness",
    "retro_dna_balance": "retro:DNA",
    "deletion_load": "deletion proxy",
    "ectopic": "ectopic",
    "body_size": "body size",
    "aquaticity": "aquaticity",
}


def edge_name(parent, child):
    return f"{edge_label_map.get(str(parent), str(parent))} -> {edge_label_map.get(str(child), str(child))}"


selected_edge_runs = [
    "genome_morphology",
    "te_genome_morphology",
    "te_genome_primary_mediumplus",
    "te_genome_ectopic_primary_mediumplus",
    "te_genome_ltr_history_primary_mediumplus",
    "te_genome_organismal_primary_mediumplus",
]
edge_frames = []
for run in selected_edge_runs:
    for kind in ["best", "average"]:
        path = RESULTS / f"{run}_{kind}_model_edges.csv"
        if not path.exists():
            continue
        df = pd.read_csv(path)
        if df.empty:
            continue
        df["run"] = run
        df["edge_set"] = kind
        edge_frames.append(df)
edges = pd.concat(edge_frames, ignore_index=True) if edge_frames else pd.DataFrame()
if not edges.empty:
    edges["coefficient"] = pd.to_numeric(edges["coefficient"], errors="coerce")
    edges["std_error"] = pd.to_numeric(edges["std_error"], errors="coerce")
    edges["edge"] = [edge_name(p, c) for p, c in zip(edges["parent"], edges["child"])]
    display(edges.sort_values(["run", "edge_set", "edge"]))

    plot_edges = edges[edges["edge_set"].eq("best")].dropna(subset=["coefficient"]).copy()
    plot_edges["label"] = plot_edges["run"] + " | " + plot_edges["edge"]
    plot_edges = plot_edges.sort_values("coefficient")
    fig, ax = plt.subplots(figsize=(12, max(5, 0.25 * len(plot_edges))))
    y = np.arange(len(plot_edges))
    xerr = plot_edges["std_error"].fillna(0).to_numpy(float)
    ax.errorbar(plot_edges["coefficient"], y, xerr=xerr, fmt="o", color="#344E41", ecolor="#B8B2A7", capsize=3)
    ax.axvline(0, color="#222222", linewidth=0.9)
    ax.set_yticks(y)
    ax.set_yticklabels(plot_edges["label"], fontsize=7)
    ax.set_xlabel("standardized coefficient")
    ax.set_title("Best-model path coefficients")
    ax.grid(axis="x", alpha=0.22)
    savefig(fig, "top50_best_model_edge_coefficients.png")
    plt.show()
else:
    print("No edge tables found.")
'''
        ),
        md("## Phylopath PDFs"),
        code(
            r'''
pdf_runs = [
    "genome_morphology",
    "te_genome_morphology",
    "te_genome_primary_mediumplus",
    "te_genome_ectopic_primary_mediumplus",
    "te_genome_ltr_history_primary_mediumplus",
    "te_genome_organismal_primary_mediumplus",
]
pdf_rows = []
for run in pdf_runs:
    for suffix in ["best_model", "average_model", "model_set", "model_summary"]:
        path = RESULTS / f"{run}_{suffix}.pdf"
        pdf_rows.append({"run": run, "artifact": suffix, "path": str(path), "exists": path.exists()})
pdf_table = pd.DataFrame(pdf_rows)
display(pdf_table)

for run in pdf_runs[:4]:
    path = RESULTS / f"{run}_best_model.pdf"
    if path.exists():
        display(HTML(f"<h4>{run}: best model PDF</h4>"))
        display(FileLink(str(path)))
        display(IFrame(src=str(path.relative_to(PROJECT_ROOT)), width="100%", height=520))
'''
        ),
        md("## Review Grid App"),
        code(
            r'''
display(HTML(f'<p>Local review app: <a href="{REVIEW_APP_URL}" target="_blank">{REVIEW_APP_URL}</a></p>'))
display(IFrame(src=REVIEW_APP_URL, width="100%", height=900))
'''
        ),
        md("## Reproducibility Commands"),
        code(
            r'''
commands = [
    "scripts/run_in_dusky.sh python path_analysis/scripts/build_top50_size_notebook.py",
    "scripts/run_in_dusky.sh jupyter nbconvert --to notebook --execute path_analysis/notebooks/top50_size_analysis_figures.ipynb --output top50_size_analysis_figures.executed.ipynb --output-dir path_analysis/notebooks --ExecutePreprocessor.timeout=300",
    "scripts/run_in_dusky.sh jupyter nbconvert --to html path_analysis/notebooks/top50_size_analysis_figures.executed.ipynb --output top50_size_analysis_figures.html --output-dir path_analysis/notebooks",
    "scripts/run_tests.sh",
]
for command in commands:
    print(command)
'''
        ),
    ]
    return nb


def main() -> None:
    NOTEBOOK_PATH.parent.mkdir(parents=True, exist_ok=True)
    nbf.write(build_notebook(), NOTEBOOK_PATH)
    print(f"Wrote {NOTEBOOK_PATH}")


if __name__ == "__main__":
    main()
