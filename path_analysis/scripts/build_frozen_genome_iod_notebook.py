#!/usr/bin/env python3
"""Build the frozen, image-quality-matched nuclear-IOD analysis notebook."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy.stats import ks_2samp, spearmanr, trim_mean


PROJECT_ROOT = Path(__file__).resolve().parents[2]
SOURCE_DIR = (
    PROJECT_ROOT
    / "path_analysis"
    / "data"
    / "external"
    / "derived"
    / "image_quality_matched_genome_iod"
)
FROZEN_PATH = SOURCE_DIR / "image_quality_matched_nuclei_frozen_reviewed.csv.gz"
MATCH_MANIFEST_PATH = SOURCE_DIR / "manifest.json"
DECISION_PATHS = [
    PROJECT_ROOT / "image_quality_matched_all_species_decisions.csv",
    PROJECT_ROOT / "image_quality_matched_replacement_review_decisions.csv",
]
REPORT_DIR = (
    PROJECT_ROOT
    / "notebooks"
    / "presentation"
    / "frozen_genome_iod_analysis"
)
NOTEBOOK_PATH = REPORT_DIR / "frozen_genome_iod_analysis.ipynb"
EXECUTED_NOTEBOOK_PATH = REPORT_DIR / "frozen_genome_iod_analysis.executed.ipynb"
FIGURE_DIR = REPORT_DIR / "figures"
HTML_PATH = REPORT_DIR / "index.html"
SPECIES_SUMMARY_PATH = REPORT_DIR / "species_relative_genome_iod_summary.csv"
IMAGE_SUMMARY_PATH = REPORT_DIR / "image_relative_genome_iod_summary.csv"
QUALITY_BALANCE_PATH = REPORT_DIR / "frozen_quality_balance.csv"
QUALITY_RESIDUAL_PATH = REPORT_DIR / "iod_quality_residual_diagnostics.csv"
ANALYSIS_MANIFEST_PATH = REPORT_DIR / "analysis_manifest.json"

MATCH_FEATURES = [
    "match_log_edge_sharpness",
    "match_log_relative_ring_noise",
]


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def hierarchical_bootstrap_equal_image_estimate(
    frame: pd.DataFrame,
    *,
    n_bootstrap: int = 2_000,
    seed: int = 20260710,
) -> tuple[float, float]:
    """Bootstrap the mean of image medians, resampling images then nuclei."""
    if frame.empty:
        return np.nan, np.nan
    image_values = {
        str(filename): group["nuc_iod"].dropna().to_numpy(dtype=float)
        for filename, group in frame.groupby("filename", sort=True)
    }
    image_values = {key: values for key, values in image_values.items() if len(values)}
    if not image_values:
        return np.nan, np.nan
    names = np.array(sorted(image_values), dtype=object)
    rng = np.random.default_rng(seed)
    estimates = np.empty(n_bootstrap, dtype=float)
    for replicate in range(n_bootstrap):
        sampled_names = rng.choice(names, size=len(names), replace=True)
        medians = []
        for name in sampled_names:
            values = image_values[str(name)]
            sampled_values = rng.choice(values, size=len(values), replace=True)
            medians.append(float(np.median(sampled_values)))
        estimates[replicate] = float(np.mean(medians))
    low, high = np.quantile(estimates, [0.025, 0.975])
    return float(low), float(high)


def summarize_species(
    frame: pd.DataFrame,
    *,
    n_bootstrap: int = 2_000,
    seed: int = 20260710,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Create equal-image species estimates and transparent sensitivity columns."""
    required = {
        "species",
        "filename",
        "specimen_group",
        "review_key",
        "nuc_iod",
        "nuc_area_um2",
        "nuc_mean_od",
    }
    missing = sorted(required - set(frame.columns))
    if missing:
        raise ValueError(f"Frozen panel is missing required columns: {missing}")

    image_aggregations: dict[str, tuple[str, str]] = {
        "n_nuclei": ("review_key", "size"),
        "image_median_iod": ("nuc_iod", "median"),
        "image_mean_iod": ("nuc_iod", "mean"),
        "image_median_nucleus_area_um2": ("nuc_area_um2", "median"),
        "image_median_mean_od": ("nuc_mean_od", "median"),
    }
    for feature in MATCH_FEATURES:
        if feature in frame.columns:
            image_aggregations[f"image_median_{feature}"] = (feature, "median")
    image_summary = (
        frame.groupby(["species", "filename", "specimen_group"], sort=True)
        .agg(**image_aggregations)
        .reset_index()
    )

    rows: list[dict[str, Any]] = []
    for species_index, (species, group) in enumerate(frame.groupby("species", sort=True)):
        images = image_summary.loc[image_summary["species"].eq(species)]
        image_medians = images["image_median_iod"].to_numpy(dtype=float)
        estimate = float(np.mean(image_medians))
        ci_low, ci_high = hierarchical_bootstrap_equal_image_estimate(
            group,
            n_bootstrap=n_bootstrap,
            seed=seed + species_index * 1009,
        )
        iod = group["nuc_iod"].dropna().to_numpy(dtype=float)
        row: dict[str, Any] = {
            "species": species,
            "n_nuclei": int(len(group)),
            "n_images": int(group["filename"].nunique()),
            "n_specimens": int(group["specimen_group"].nunique()),
            "iod_equal_image_estimate": estimate,
            "iod_equal_image_ci_low": ci_low,
            "iod_equal_image_ci_high": ci_high,
            "iod_pooled_nucleus_median": float(np.median(iod)),
            "iod_pooled_nucleus_mean": float(np.mean(iod)),
            "iod_trimmed_mean_10pct": float(trim_mean(iod, 0.10)),
            "iod_nucleus_sd": float(np.std(iod, ddof=1)),
            "iod_nucleus_iqr": float(np.quantile(iod, 0.75) - np.quantile(iod, 0.25)),
            "iod_nucleus_min": float(np.min(iod)),
            "iod_nucleus_max": float(np.max(iod)),
            "image_median_iod_min": float(np.min(image_medians)),
            "image_median_iod_max": float(np.max(image_medians)),
            "image_median_iod_range": float(np.max(image_medians) - np.min(image_medians)),
            "median_nucleus_area_um2": float(group["nuc_area_um2"].median()),
            "median_nucleus_mean_od": float(group["nuc_mean_od"].median()),
            "estimate_support": (
                "single_observed_image"
                if group["filename"].nunique() == 1
                else "multiple_observed_images"
            ),
        }
        row["pooled_median_difference_pct"] = float(
            100.0 * (row["iod_pooled_nucleus_median"] / estimate - 1.0)
        )
        row["trimmed_mean_difference_pct"] = float(
            100.0 * (row["iod_trimmed_mean_10pct"] / estimate - 1.0)
        )
        rows.append(row)

    species_summary = pd.DataFrame(rows)
    anchor = float(species_summary["iod_equal_image_estimate"].median())
    species_summary["relative_iod_anchor"] = anchor
    species_summary["relative_iod_index"] = (
        species_summary["iod_equal_image_estimate"] / anchor
    )
    species_summary["relative_iod_ci_low"] = (
        species_summary["iod_equal_image_ci_low"] / anchor
    )
    species_summary["relative_iod_ci_high"] = (
        species_summary["iod_equal_image_ci_high"] / anchor
    )
    species_summary["relative_iod_rank"] = (
        species_summary["iod_equal_image_estimate"]
        .rank(method="min", ascending=False)
        .astype(int)
    )
    species_summary = species_summary.sort_values(
        ["relative_iod_rank", "species"], kind="mergesort"
    ).reset_index(drop=True)
    return species_summary, image_summary


def frozen_quality_diagnostics(
    frame: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Measure remaining quality balance and within-species IOD associations."""
    template_species = str(frame["quality_template_species"].iloc[0])
    z_features = [f"z_{feature}" for feature in MATCH_FEATURES]
    missing = sorted(set(z_features) - set(frame.columns))
    if missing:
        raise ValueError(f"Frozen panel is missing standardized quality columns: {missing}")
    template = frame.loc[frame["species"].eq(template_species)]
    balance_rows: list[dict[str, Any]] = []
    for species, group in frame.groupby("species", sort=True):
        row: dict[str, Any] = {
            "species": species,
            "template_species": template_species,
            "n_nuclei": int(len(group)),
        }
        absolute_differences = []
        ks_distances = []
        for feature, z_feature in zip(MATCH_FEATURES, z_features):
            difference = float(group[z_feature].mean() - template[z_feature].mean())
            ks_distance = float(
                ks_2samp(group[z_feature], template[z_feature]).statistic
            )
            row[f"standardized_mean_difference_{feature}"] = difference
            row[f"ks_distance_{feature}"] = ks_distance
            absolute_differences.append(abs(difference))
            ks_distances.append(ks_distance)
        row["max_abs_standardized_mean_difference"] = max(absolute_differences)
        row["max_ks_distance"] = max(ks_distances)
        row["passes_prespecified_balance_thresholds"] = bool(
            row["max_abs_standardized_mean_difference"] <= 0.10
            and row["max_ks_distance"] <= 0.25
        )
        balance_rows.append(row)

    log_iod = np.log(frame["nuc_iod"].astype(float))
    centered_log_iod = log_iod - log_iod.groupby(frame["species"]).transform("median")
    residual_rows = []
    for feature in MATCH_FEATURES:
        values = frame[feature].astype(float)
        centered_values = values - values.groupby(frame["species"]).transform("median")
        rho, p_value = spearmanr(centered_values, centered_log_iod)
        residual_rows.append(
            {
                "quality_feature": feature,
                "within_species_centered_spearman_rho": float(rho),
                "p_value_descriptive_only": float(p_value),
                "n_nuclei": int(len(frame)),
            }
        )
    return pd.DataFrame(balance_rows), pd.DataFrame(residual_rows)


def load_latest_decisions() -> pd.DataFrame:
    parts = []
    for path in DECISION_PATHS:
        frame = pd.read_csv(path, low_memory=False)
        parts.append(frame[["species", "review_key", "decision"]].copy())
    decisions = pd.concat(parts, ignore_index=True)
    decisions["decision"] = decisions["decision"].astype(str).str.strip().str.lower()
    return decisions.drop_duplicates(["species", "review_key"], keep="last")


def build_analysis_outputs(*, n_bootstrap: int = 2_000) -> dict[str, Any]:
    REPORT_DIR.mkdir(parents=True, exist_ok=True)
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    frozen = pd.read_csv(FROZEN_PATH, low_memory=False)
    match_manifest = json.loads(MATCH_MANIFEST_PATH.read_text())
    decisions = load_latest_decisions()

    frozen_keys = set(
        map(
            tuple,
            frozen[["species", "review_key"]].itertuples(index=False, name=None),
        )
    )
    rejected_keys = set(
        map(
            tuple,
            decisions.loc[
                decisions["decision"].eq("problem"), ["species", "review_key"]
            ].itertuples(index=False, name=None),
        )
    )
    rejected_leaks = frozen_keys & rejected_keys
    counts = frozen.groupby("species").size()
    if match_manifest.get("freeze_status") != "frozen_reviewed_only_no_further_replacement":
        raise ValueError("Source panel is not marked frozen.")
    if not frozen["review_status"].eq("reviewed_keep").all():
        raise ValueError("Frozen panel contains a nucleus that was not reviewed keep.")
    if frozen.duplicated(["species", "review_key"]).any():
        raise ValueError("Frozen panel contains duplicate review keys.")
    if rejected_leaks:
        raise ValueError(f"Rejected nuclei leaked into frozen panel: {len(rejected_leaks)}")
    if counts.min() < 30:
        raise ValueError("At least one primary species has fewer than 30 reviewed nuclei.")

    species_summary, image_summary = summarize_species(
        frozen, n_bootstrap=n_bootstrap
    )
    quality_balance, quality_residual = frozen_quality_diagnostics(frozen)
    species_summary.to_csv(SPECIES_SUMMARY_PATH, index=False, float_format="%.8f")
    image_summary.to_csv(IMAGE_SUMMARY_PATH, index=False, float_format="%.8f")
    quality_balance.to_csv(QUALITY_BALANCE_PATH, index=False, float_format="%.8f")
    quality_residual.to_csv(QUALITY_RESIDUAL_PATH, index=False, float_format="%.8f")

    top = species_summary.iloc[0]
    bottom = species_summary.iloc[-1]
    max_sensitivity = float(
        species_summary[
            ["pooled_median_difference_pct", "trimmed_mean_difference_pct"]
        ]
        .abs()
        .to_numpy()
        .max()
    )
    manifest: dict[str, Any] = {
        "analysis": "Frozen image-quality-matched relative nuclear-IOD analysis",
        "frozen_source_csv": str(FROZEN_PATH.resolve()),
        "frozen_source_sha256": sha256_file(FROZEN_PATH),
        "source_match_manifest": str(MATCH_MANIFEST_PATH.resolve()),
        "source_match_manifest_sha256": sha256_file(MATCH_MANIFEST_PATH),
        "decision_files": [str(path.resolve()) for path in DECISION_PATHS],
        "decision_file_sha256": {
            str(path.resolve()): sha256_file(path) for path in DECISION_PATHS
        },
        "n_frozen_nuclei": int(len(frozen)),
        "n_primary_species": int(frozen["species"].nunique()),
        "n_images": int(frozen["filename"].nunique()),
        "n_specimens": int(frozen["specimen_group"].nunique()),
        "minimum_nuclei_per_species": int(counts.min()),
        "maximum_nuclei_per_species": int(counts.max()),
        "n_rejected_nuclei_leaked": int(len(rejected_leaks)),
        "primary_estimator": "mean across image-specific median nuclear IOD values",
        "bootstrap": {
            "method": "hierarchical resampling of images then nuclei within images",
            "replicates": int(n_bootstrap),
            "seed": 20260710,
            "interval": "2.5th and 97.5th percentiles",
            "scope": "conditional on the observed images/specimens",
        },
        "relative_index_anchor": "median of the 20 species equal-image IOD estimates",
        "relative_index_anchor_iod": float(
            species_summary["relative_iod_anchor"].iloc[0]
        ),
        "highest_relative_iod_species": str(top["species"]),
        "highest_relative_iod_index": float(top["relative_iod_index"]),
        "lowest_relative_iod_species": str(bottom["species"]),
        "lowest_relative_iod_index": float(bottom["relative_iod_index"]),
        "highest_to_lowest_fold_difference": float(
            top["iod_equal_image_estimate"] / bottom["iod_equal_image_estimate"]
        ),
        "maximum_aggregation_sensitivity_pct": max_sensitivity,
        "single_image_species": sorted(
            species_summary.loc[
                species_summary["n_images"].eq(1), "species"
            ].tolist()
        ),
        "worst_abs_standardized_quality_mean_difference": float(
            quality_balance["max_abs_standardized_mean_difference"].max()
        ),
        "worst_quality_ks_distance": float(
            quality_balance["max_ks_distance"].max()
        ),
        "all_species_pass_quality_balance_thresholds": bool(
            quality_balance["passes_prespecified_balance_thresholds"].all()
        ),
        "absolute_genome_size_claimed": False,
        "interpretation": (
            "Relative nuclear-IOD is an image-derived DNA-content proxy. "
            "It is not an absolute genome-size estimate, C-value, or picogram measurement."
        ),
        "limited_overlap_species_excluded_from_primary": match_manifest.get(
            "limited_overlap_species", []
        ),
        "outputs": {
            "species_summary_csv": str(SPECIES_SUMMARY_PATH.resolve()),
            "image_summary_csv": str(IMAGE_SUMMARY_PATH.resolve()),
            "quality_balance_csv": str(QUALITY_BALANCE_PATH.resolve()),
            "quality_residual_csv": str(QUALITY_RESIDUAL_PATH.resolve()),
            "source_notebook": str(NOTEBOOK_PATH.resolve()),
            "executed_notebook": str(EXECUTED_NOTEBOOK_PATH.resolve()),
            "html": str(HTML_PATH.resolve()),
        },
    }
    ANALYSIS_MANIFEST_PATH.write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    return manifest


def _nbformat() -> Any:
    import nbformat

    return nbformat


def md(text: str) -> Any:
    nbf = _nbformat()
    return nbf.v4.new_markdown_cell(text.strip() + "\n")


def code(text: str) -> Any:
    nbf = _nbformat()
    return nbf.v4.new_code_cell(text.strip() + "\n")


def build_notebook() -> Any:
    nbf = _nbformat()
    nb = nbf.v4.new_notebook()
    nb.metadata = {
        "kernelspec": {
            "display_name": "Python 3 (Dusky)",
            "language": "python",
            "name": "python3",
        },
        "language_info": {"name": "python", "pygments_lexer": "ipython3"},
    }
    nb.cells = [
        md(
            """
            # Frozen Quality-Matched Nuclear-IOD Genome-Size Analysis

            This notebook analyzes every nucleus in the final reviewed panel:
            721 manually accepted nuclei from 20 species. The primary estimator
            is the mean of the image-specific median IOD values, giving every
            observed image/specimen equal weight.

            **Interpretation boundary:** this is a relative nuclear-IOD
            genome-size proxy. It is not an absolute genome-size estimate,
            C-value, or picogram measurement because no independent
            DNA-content standard was imaged and calibrated with these samples.
            """
        ),
        md("## Exact frozen inputs and analysis setup"),
        code(
            r'''
from pathlib import Path
import sys
import json
import hashlib
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import spearmanr
from IPython.display import display, HTML, Image


def find_project_root(start=None) -> Path:
    current = (start or Path.cwd()).resolve()
    for candidate in [current, *current.parents]:
        if (candidate / "paths.yaml").exists():
            return candidate
    raise FileNotFoundError("Could not locate project root")


PROJECT_ROOT = find_project_root()
sys.path.insert(0, str(PROJECT_ROOT / "path_analysis" / "scripts"))
import build_frozen_genome_iod_notebook as analysis

FROZEN_PATH = analysis.FROZEN_PATH
MATCH_MANIFEST_PATH = analysis.MATCH_MANIFEST_PATH
REPORT_DIR = analysis.REPORT_DIR
FIGURE_DIR = analysis.FIGURE_DIR
FIGURE_DIR.mkdir(parents=True, exist_ok=True)

frozen = pd.read_csv(FROZEN_PATH, low_memory=False)
match_manifest = json.loads(MATCH_MANIFEST_PATH.read_text())
species_summary, image_summary = analysis.summarize_species(
    frozen, n_bootstrap=2000, seed=20260710
)
quality_balance, quality_residual = analysis.frozen_quality_diagnostics(frozen)
anchor = float(species_summary["relative_iod_anchor"].iloc[0])

plt.rcParams.update({
    "figure.dpi": 115,
    "savefig.dpi": 190,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "font.size": 10,
})
pd.set_option("display.max_columns", 80)
pd.set_option("display.max_rows", 800)
pd.set_option("display.width", 180)


def short_species(value):
    return str(value).replace("D. ", "")


def save_figure(fig, filename):
    path = FIGURE_DIR / filename
    fig.savefig(path, bbox_inches="tight", facecolor="white")
    return path


print("Frozen source:", FROZEN_PATH)
print("Frozen SHA-256:", analysis.sha256_file(FROZEN_PATH))
print("Rows:", len(frozen), "| Species:", frozen["species"].nunique())
print("Relative-index anchor IOD:", round(anchor, 3))
'''
        ),
        md(
            """
            ## Freeze and source audit

            These checks are hard failures: all rows must be reviewed keeps,
            no rejected key can appear, every species must retain at least 30
            nuclei, and the source manifest must mark the panel frozen.
            """
        ),
        code(
            r'''
decisions = analysis.load_latest_decisions()
problem_keys = set(map(
    tuple,
    decisions.loc[decisions["decision"].eq("problem"), ["species", "review_key"]]
    .itertuples(index=False, name=None),
))
frozen_keys = set(map(
    tuple,
    frozen[["species", "review_key"]].itertuples(index=False, name=None),
))
counts = frozen.groupby("species").size()
audit = pd.DataFrame({
    "check": [
        "manifest freeze status",
        "reviewed nuclei",
        "primary species",
        "images/specimens",
        "minimum nuclei per species",
        "maximum nuclei per species",
        "duplicate keys",
        "rejected-key leaks",
        "all rows explicitly reviewed keep",
    ],
    "value": [
        match_manifest["freeze_status"],
        len(frozen),
        frozen["species"].nunique(),
        f'{frozen["filename"].nunique()} / {frozen["specimen_group"].nunique()}',
        int(counts.min()),
        int(counts.max()),
        int(frozen.duplicated(["species", "review_key"]).sum()),
        len(frozen_keys & problem_keys),
        bool(frozen["review_status"].eq("reviewed_keep").all()),
    ],
})
display(audit)
assert match_manifest["freeze_status"] == "frozen_reviewed_only_no_further_replacement"
assert frozen.duplicated(["species", "review_key"]).sum() == 0
assert not (frozen_keys & problem_keys)
assert frozen["review_status"].eq("reviewed_keep").all()
assert counts.min() >= 30
'''
        ),
        md("## Complete species-level results"),
        code(
            r'''
summary_columns = [
    "relative_iod_rank",
    "species",
    "n_nuclei",
    "n_images",
    "n_specimens",
    "iod_equal_image_estimate",
    "iod_equal_image_ci_low",
    "iod_equal_image_ci_high",
    "relative_iod_index",
    "relative_iod_ci_low",
    "relative_iod_ci_high",
    "iod_pooled_nucleus_median",
    "iod_trimmed_mean_10pct",
    "pooled_median_difference_pct",
    "trimmed_mean_difference_pct",
    "median_nucleus_area_um2",
    "median_nucleus_mean_od",
    "image_median_iod_min",
    "image_median_iod_max",
    "estimate_support",
]
display(
    species_summary[summary_columns]
    .style.format({
        "iod_equal_image_estimate": "{:.2f}",
        "iod_equal_image_ci_low": "{:.2f}",
        "iod_equal_image_ci_high": "{:.2f}",
        "relative_iod_index": "{:.3f}",
        "relative_iod_ci_low": "{:.3f}",
        "relative_iod_ci_high": "{:.3f}",
        "iod_pooled_nucleus_median": "{:.2f}",
        "iod_trimmed_mean_10pct": "{:.2f}",
        "pooled_median_difference_pct": "{:+.1f}%",
        "trimmed_mean_difference_pct": "{:+.1f}%",
        "median_nucleus_area_um2": "{:.2f}",
        "median_nucleus_mean_od": "{:.4f}",
        "image_median_iod_min": "{:.2f}",
        "image_median_iod_max": "{:.2f}",
    })
)
'''
        ),
        md(
            """
            ## All 721 raw nuclei

            Every accepted nucleus is included below. IOD, its two algebraic
            components, image identity, and the two matching-quality features
            are visible. The frozen source download retains all 189 provenance
            and QC columns.
            """
        ),
        code(
            r'''
raw_columns = [
    "species",
    "filename",
    "specimen_group",
    "nucleus_label",
    "review_key",
    "nuc_iod",
    "nuc_area_um2",
    "nuc_mean_od",
    "match_log_edge_sharpness",
    "match_log_relative_ring_noise",
    "quality_match_distance",
    "review_status",
]
raw_display = frozen[raw_columns].sort_values(
    ["species", "filename", "nuc_iod"], kind="mergesort"
)
relative_source_link = (
    "../../../path_analysis/data/external/derived/image_quality_matched_genome_iod/"
    "image_quality_matched_nuclei_frozen_reviewed.csv.gz"
)
display(HTML(
    f'<p><a href="{relative_source_link}">Download the complete frozen 189-column CSV</a></p>'
    '<div style="max-height:620px;overflow:auto;border:1px solid #ccd3da;padding:4px">'
    + raw_display.to_html(index=False, float_format=lambda value: f"{value:.5f}")
    + "</div>"
))
print("Displayed raw nuclei:", len(raw_display))
'''
        ),
        md("## Raw nucleus-level IOD distributions"),
        code(
            r'''
order = species_summary.sort_values("iod_equal_image_estimate")["species"].tolist()
rng = np.random.default_rng(20260710)
fig, ax = plt.subplots(figsize=(13, 12))
for y, species in enumerate(order):
    values = frozen.loc[frozen["species"].eq(species), "nuc_iod"].to_numpy()
    jitter = rng.normal(0, 0.065, size=len(values))
    ax.scatter(values, y + jitter, s=17, alpha=0.52, color="#426B8A", edgecolor="none")
    estimate = species_summary.loc[
        species_summary["species"].eq(species), "iod_equal_image_estimate"
    ].iloc[0]
    ax.scatter(estimate, y, marker="D", s=48, color="#B33A3A", zorder=5)
ax.set_yticks(range(len(order)), [short_species(value) for value in order])
ax.set_xlabel("Nuclear IOD (raw image units)")
ax.set_ylabel("Species")
ax.set_title("Every reviewed nucleus; red diamond = equal-image species estimate")
ax.grid(axis="x", alpha=0.18)
save_figure(fig, "01_all_raw_nuclear_iod.png")
plt.show()
'''
        ),
        md(
            """
            ## Primary relative IOD estimates and conditional uncertainty

            The species index is normalized so the median species estimate is
            1.0. Confidence intervals use hierarchical resampling of images and
            then nuclei within images. They are conditional on the images and
            specimens observed here; the three single-image species cannot
            express between-image uncertainty. The index intervals divide the
            raw-IOD intervals by the observed anchor and therefore treat that
            anchor as fixed.
            """
        ),
        code(
            r'''
plot_data = species_summary.sort_values("relative_iod_index")
y = np.arange(len(plot_data))
x = plot_data["relative_iod_index"].to_numpy()
low = plot_data["relative_iod_ci_low"].to_numpy()
high = plot_data["relative_iod_ci_high"].to_numpy()
colors = np.where(plot_data["n_images"].eq(1), "#D17B29", "#2E6F95")
fig, ax = plt.subplots(figsize=(11, 11))
ax.hlines(y, low, high, color=colors, linewidth=2)
ax.scatter(x, y, color=colors, s=55, zorder=3)
ax.axvline(1.0, color="#555555", linestyle="--", linewidth=1)
ax.set_yticks(y, [short_species(value) for value in plot_data["species"]])
ax.set_xlabel("Relative nuclear-IOD index (median species = 1.0)")
ax.set_title("Equal-image species estimates with 95% hierarchical bootstrap intervals")
ax.grid(axis="x", alpha=0.18)
ax.text(
    0.01,
    0.01,
    "Orange = one observed image/specimen",
    transform=ax.transAxes,
    color="#9A541C",
)
save_figure(fig, "02_relative_iod_estimates.png")
plt.show()
'''
        ),
        md(
            """
            ## Measured-only time-calibrated phylogeny

            This panel includes all 21 species in the frozen literal-largest
            cell panel. Every species has observed cell and corresponding
            nucleus distributions. Relative nuclear IOD is shown for the 20
            common-support genome species; *D. ochrophaeus* is retained on the
            tree with its observed size data and an explicitly empty genome
            slot because its image-quality overlap was limited. Every violin
            is a bootstrap distribution of the plotted estimator, genome is
            never expressed in picograms, and no species is filled by
            phylogenetic imputation.
            """
        ),
        code(
            r'''
phylogeny_figure = FIGURE_DIR / "07_measured_phylogeny_genome_nucleus_cell.png"
phylogeny_summary_path = REPORT_DIR / "phylogeny_genome_nucleus_cell_summary.csv"
phylogeny_correlation_path = REPORT_DIR / "phylogeny_genome_nucleus_cell_correlations.csv"
assert phylogeny_figure.exists(), "Run the measured phylogeny figure builder first"
assert phylogeny_summary_path.exists()
assert phylogeny_correlation_path.exists()
display(Image(filename=str(phylogeny_figure)))
display(pd.read_csv(phylogeny_summary_path))
display(pd.read_csv(phylogeny_correlation_path))
'''
        ),
        md(
            """
            ## Pairwise genome-proxy, nucleus, and cell relationships

            These three panels show every unique pair of traits. Points are
            species estimates and bars are their bootstrap intervals. Positive
            relationships are biologically consistent with the classical
            nucleotypic expectation, but they are descriptive rather than
            causal: IOD contains nuclear area algebraically, the size panel is
            an upper-tail top-50 estimand, and shared ancestry is not corrected
            in these correlations.
            """
        ),
        code(
            r'''
pairwise_figure = FIGURE_DIR / "08_pairwise_genome_nucleus_cell_relationships.png"
assert pairwise_figure.exists(), "Run the measured phylogeny figure builder first"
display(Image(filename=str(pairwise_figure)))
display(pd.read_csv(phylogeny_correlation_path))
'''
        ),
        md("## Image-to-image variation within species"),
        code(
            r'''
species_order = species_summary.sort_values("iod_equal_image_estimate")["species"].tolist()
fig, ax = plt.subplots(figsize=(12, 11))
for y, species in enumerate(species_order):
    image_values = image_summary.loc[
        image_summary["species"].eq(species), "image_median_iod"
    ].to_numpy()
    estimate = species_summary.loc[
        species_summary["species"].eq(species), "iod_equal_image_estimate"
    ].iloc[0]
    if len(image_values) > 1:
        ax.hlines(y, image_values.min(), image_values.max(), color="#9AA6B2", linewidth=2)
    ax.scatter(image_values, np.full(len(image_values), y), s=52, color="#365F78", alpha=0.85)
    ax.scatter(estimate, y, s=58, marker="D", color="#B33A3A", zorder=4)
ax.set_yticks(range(len(species_order)), [short_species(value) for value in species_order])
ax.set_xlabel("Image-specific median nuclear IOD")
ax.set_title("Observed image/specimen medians; red diamond = equal-image estimate")
ax.grid(axis="x", alpha=0.18)
save_figure(fig, "03_image_level_iod.png")
plt.show()
display(image_summary.sort_values(["species", "filename"]))
'''
        ),
        md(
            """
            ## What makes up IOD?

            Nuclear IOD is algebraically the nuclear mask area in pixels
            multiplied by mean optical density. Its association with nuclear
            area therefore cannot independently validate a biological
            genome-size mechanism. Both components are shown explicitly.
            """
        ),
        code(
            r'''
component_data = species_summary.copy()
rho_area, p_area = spearmanr(
    component_data["median_nucleus_area_um2"],
    component_data["iod_equal_image_estimate"],
)
rho_od, p_od = spearmanr(
    component_data["median_nucleus_mean_od"],
    component_data["iod_equal_image_estimate"],
)
fig, axes = plt.subplots(1, 2, figsize=(15, 7))
specs = [
    (
        "median_nucleus_area_um2",
        "Median nucleus area (µm²)",
        rho_area,
        p_area,
        "#477998",
    ),
    (
        "median_nucleus_mean_od",
        "Median nuclear mean optical density",
        rho_od,
        p_od,
        "#7A5195",
    ),
]
for ax, (column, xlabel, rho, p_value, color) in zip(axes, specs):
    ax.scatter(
        component_data[column],
        component_data["relative_iod_index"],
        s=48,
        color=color,
        alpha=0.85,
    )
    for row in component_data.itertuples():
        ax.annotate(
            short_species(row.species),
            (getattr(row, column), row.relative_iod_index),
            xytext=(3, 3),
            textcoords="offset points",
            fontsize=7,
        )
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Relative nuclear-IOD index")
    ax.set_title(f"Spearman ρ = {rho:.2f}; descriptive p = {p_value:.3g}")
    ax.grid(alpha=0.16)
fig.suptitle("Species IOD estimates against their two measurement components", y=1.02)
save_figure(fig, "04_iod_components.png")
plt.show()
print(f"IOD vs median nucleus area: rho={rho_area:.3f}, p={p_area:.4g}")
print(f"IOD vs median mean OD: rho={rho_od:.3f}, p={p_od:.4g}")
'''
        ),
        md("## Aggregation-method sensitivity"),
        code(
            r'''
sensitivity = species_summary[[
    "species",
    "iod_equal_image_estimate",
    "iod_pooled_nucleus_median",
    "iod_trimmed_mean_10pct",
]].copy()
method_columns = [
    "iod_equal_image_estimate",
    "iod_pooled_nucleus_median",
    "iod_trimmed_mean_10pct",
]
method_labels = ["Equal-image", "Pooled median", "10% trimmed mean"]
for column in method_columns:
    sensitivity[column + "_index"] = (
        sensitivity[column] / sensitivity[column].median()
    )
order = species_summary.sort_values("relative_iod_index")["species"].tolist()
fig, ax = plt.subplots(figsize=(12, 11))
offsets = [-0.18, 0.0, 0.18]
markers = ["D", "o", "s"]
colors = ["#B33A3A", "#365F78", "#6A994E"]
for y, species in enumerate(order):
    row = sensitivity.loc[sensitivity["species"].eq(species)].iloc[0]
    values = [row[column + "_index"] for column in method_columns]
    ax.plot(values, [y] * 3, color="#B8C0C8", linewidth=1, zorder=1)
    for offset, value, marker, color in zip(offsets, values, markers, colors):
        ax.scatter(value, y + offset, marker=marker, s=42, color=color, zorder=2)
ax.axvline(1.0, color="#555555", linestyle="--", linewidth=1)
ax.set_yticks(range(len(order)), [short_species(value) for value in order])
ax.set_xlabel("Within-method relative index (method median = 1.0)")
ax.set_title("Species estimates are shown under three transparent aggregation choices")
handles = [
    plt.Line2D([], [], marker=marker, linestyle="", color=color, label=label)
    for marker, color, label in zip(markers, colors, method_labels)
]
ax.legend(handles=handles, frameon=False, loc="lower right")
ax.grid(axis="x", alpha=0.18)
save_figure(fig, "05_aggregation_sensitivity.png")
plt.show()
display(
    species_summary[[
        "species",
        "pooled_median_difference_pct",
        "trimmed_mean_difference_pct",
    ]]
    .assign(max_abs_shift_pct=lambda x: x[[
        "pooled_median_difference_pct", "trimmed_mean_difference_pct"
    ]].abs().max(axis=1))
    .sort_values("max_abs_shift_pct", ascending=False)
)
'''
        ),
        md("## Technical-quality balance and residual IOD associations"),
        code(
            r'''
display(quality_balance)
display(quality_residual)

balance_plot = quality_balance.sort_values("species")
mean_columns = [
    "standardized_mean_difference_match_log_edge_sharpness",
    "standardized_mean_difference_match_log_relative_ring_noise",
]
matrix = balance_plot[mean_columns].to_numpy()
fig, axes = plt.subplots(1, 2, figsize=(15, 10), gridspec_kw={"width_ratios": [1.0, 1.25]})
image = axes[0].imshow(matrix, cmap="coolwarm", vmin=-0.10, vmax=0.10, aspect="auto")
axes[0].set_yticks(range(len(balance_plot)), [short_species(v) for v in balance_plot["species"]])
axes[0].set_xticks([0, 1], ["Edge sharpness", "Ring noise"], rotation=25, ha="right")
axes[0].set_title("Standardized mean difference vs reviewed template")
fig.colorbar(image, ax=axes[0], shrink=0.7, label="Standardized mean difference")

centered_log_iod = np.log(frozen["nuc_iod"]) - np.log(frozen["nuc_iod"]).groupby(
    frozen["species"]
).transform("median")
feature_colors = ["#477998", "#9C6644"]
for feature, color in zip(analysis.MATCH_FEATURES, feature_colors):
    centered_quality = frozen[feature] - frozen[feature].groupby(
        frozen["species"]
    ).transform("median")
    rho = quality_residual.loc[
        quality_residual["quality_feature"].eq(feature),
        "within_species_centered_spearman_rho",
    ].iloc[0]
    axes[1].scatter(
        centered_quality,
        centered_log_iod,
        s=13,
        alpha=0.30,
        color=color,
        label=f"{feature.replace('match_log_', '')}: ρ={rho:.2f}",
    )
axes[1].axhline(0, color="#777777", linewidth=0.8)
axes[1].axvline(0, color="#777777", linewidth=0.8)
axes[1].set_xlabel("Within-species centered quality feature")
axes[1].set_ylabel("Within-species centered log IOD")
axes[1].set_title("Residual technical-quality association")
axes[1].legend(frameon=False)
axes[1].grid(alpha=0.14)
save_figure(fig, "06_quality_diagnostics.png")
plt.show()
print(
    "Worst absolute standardized mean difference:",
    round(quality_balance["max_abs_standardized_mean_difference"].max(), 3),
)
print(
    "Worst KS distance:",
    round(quality_balance["max_ks_distance"].max(), 3),
)
'''
        ),
        md("## Direct reading of the current results"),
        code(
            r'''
top = species_summary.iloc[0]
bottom = species_summary.iloc[-1]
fold = top["iod_equal_image_estimate"] / bottom["iod_equal_image_estimate"]
single_image = species_summary.loc[
    species_summary["n_images"].eq(1), "species"
].tolist()
intervals_crossing_anchor = int((
    species_summary["relative_iod_ci_low"].le(1.0)
    & species_summary["relative_iod_ci_high"].ge(1.0)
).sum())
method_shift = species_summary.assign(
    max_shift=lambda x: x[[
        "pooled_median_difference_pct", "trimmed_mean_difference_pct"
    ]].abs().max(axis=1)
).sort_values("max_shift", ascending=False).iloc[0]
rho_area, p_area = spearmanr(
    species_summary["median_nucleus_area_um2"],
    species_summary["iod_equal_image_estimate"],
)
rho_od, p_od = spearmanr(
    species_summary["median_nucleus_mean_od"],
    species_summary["iod_equal_image_estimate"],
)

display(HTML(f"""
<div style="border-left:5px solid #365F78;background:#f4f7f9;padding:14px 18px">
<p><b>Largest relative IOD:</b> {top['species']} ({top['relative_iod_index']:.3f}× the species-median anchor).</p>
<p><b>Smallest relative IOD:</b> {bottom['species']} ({bottom['relative_iod_index']:.3f}×).</p>
<p><b>Observed range:</b> {fold:.2f}-fold from highest to lowest.</p>
<p><b>Rank resolution:</b> {intervals_crossing_anchor} of {len(species_summary)}
intervals cross the 1.0 anchor, so most exact middle ranks should not be treated as
sharply separated.</p>
<p><b>Aggregation sensitivity:</b> the largest shift from equal-image weighting is
{method_shift['max_shift']:.1f}% for {method_shift['species']}.</p>
<p><b>IOD components:</b> species IOD correlates with median nucleus area
(ρ={rho_area:.2f}) and median mean optical density (ρ={rho_od:.2f}). This is expected
because both are algebraic components of IOD.</p>
<p><b>Single-image limitation:</b> {', '.join(single_image)} have only one observed
image/specimen, so their interval cannot measure between-image variation.</p>
</div>
"""))
'''
        ),
        md(
            """
            ## Conclusions and next validation step

            1. The frozen panel supports a stable **relative ranking** of
               nuclear IOD across the 20 image-quality-compatible species.
            2. Equal-image weighting is the primary estimate because it avoids
               letting images with more accepted nuclei dominate.
            3. The raw values, image-specific medians, bootstrap intervals,
               alternative aggregations, and technical-quality diagnostics
               remain visible rather than collapsing the evidence to one
               number per species.
            4. Most middle-ranked species have overlapping intervals. The point
               ranking is descriptive, not evidence that every adjacent pair
               differs biologically.
            5. This analysis cannot convert IOD to absolute genome size. That
               requires an independently measured DNA-content standard
               processed under the same staining and imaging protocol.
            6. *D. ochrophaeus* remains excluded from the primary comparison
               because its technical image-quality distribution had limited
               overlap; it belongs only in a labeled sensitivity analysis.
            """
        ),
        md("## Reproducibility and exported files"),
        code(
            r'''
exports = pd.DataFrame([
    {"artifact": "Frozen 721-nucleus source", "path": str(analysis.FROZEN_PATH)},
    {"artifact": "Species result table", "path": str(analysis.SPECIES_SUMMARY_PATH)},
    {"artifact": "Image result table", "path": str(analysis.IMAGE_SUMMARY_PATH)},
    {"artifact": "Frozen quality balance", "path": str(analysis.QUALITY_BALANCE_PATH)},
    {"artifact": "Residual quality diagnostics", "path": str(analysis.QUALITY_RESIDUAL_PATH)},
    {"artifact": "Measured phylogeny figure", "path": str(analysis.FIGURE_DIR / "07_measured_phylogeny_genome_nucleus_cell.png")},
    {"artifact": "Pairwise trait figure", "path": str(analysis.FIGURE_DIR / "08_pairwise_genome_nucleus_cell_relationships.png")},
    {"artifact": "Measured phylogeny source table", "path": str(analysis.REPORT_DIR / "phylogeny_genome_nucleus_cell_summary.csv")},
    {"artifact": "Measured phylogeny correlation table", "path": str(analysis.REPORT_DIR / "phylogeny_genome_nucleus_cell_correlations.csv")},
    {"artifact": "Analysis manifest", "path": str(analysis.ANALYSIS_MANIFEST_PATH)},
    {"artifact": "Executed notebook", "path": str(analysis.EXECUTED_NOTEBOOK_PATH)},
    {"artifact": "Rendered HTML", "path": str(analysis.HTML_PATH)},
])
display(exports)
print("Build source notebook:")
print("  uv run python path_analysis/scripts/build_frozen_genome_iod_notebook.py")
print("Build measured phylogeny figure:")
print(
    "  uv run --with biopython python "
    "path_analysis/scripts/build_frozen_genome_iod_phylogeny_figure.py"
)
print("Execute notebook:")
print(
    "  uv run jupyter nbconvert --to notebook --execute "
    "notebooks/presentation/frozen_genome_iod_analysis/frozen_genome_iod_analysis.ipynb "
    "--output frozen_genome_iod_analysis.executed.ipynb "
    "--output-dir notebooks/presentation/frozen_genome_iod_analysis "
    "--ExecutePreprocessor.timeout=300"
)
print("Render HTML:")
print(
    "  uv run jupyter nbconvert --to html "
    "notebooks/presentation/frozen_genome_iod_analysis/"
    "frozen_genome_iod_analysis.executed.ipynb "
    "--output index.html "
    "--output-dir notebooks/presentation/frozen_genome_iod_analysis"
)
'''
        ),
    ]
    return nb


def write_notebook() -> None:
    nbf = _nbformat()
    NOTEBOOK_PATH.parent.mkdir(parents=True, exist_ok=True)
    nbf.write(build_notebook(), NOTEBOOK_PATH)


def finalize_manifest() -> dict[str, Any]:
    manifest = json.loads(ANALYSIS_MANIFEST_PATH.read_text())
    for label, path in [
        ("source_notebook", NOTEBOOK_PATH),
        ("executed_notebook", EXECUTED_NOTEBOOK_PATH),
        ("html", HTML_PATH),
        ("species_summary_csv", SPECIES_SUMMARY_PATH),
        ("image_summary_csv", IMAGE_SUMMARY_PATH),
        ("quality_balance_csv", QUALITY_BALANCE_PATH),
        ("quality_residual_csv", QUALITY_RESIDUAL_PATH),
    ]:
        if path.exists():
            manifest["outputs"][f"{label}_sha256"] = sha256_file(path)
    ANALYSIS_MANIFEST_PATH.write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bootstrap-replicates", type=int, default=2_000)
    parser.add_argument("--finalize", action="store_true")
    args = parser.parse_args()
    if args.finalize:
        print(json.dumps(finalize_manifest(), indent=2))
        return
    manifest = build_analysis_outputs(n_bootstrap=args.bootstrap_replicates)
    write_notebook()
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
