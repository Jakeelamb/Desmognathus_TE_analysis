#!/usr/bin/env python3
"""Build apples-to-apples genome-IOD sensitivity panels.

This script is intentionally additive. It does not replace the current top-50
cell/nucleus results. It builds a separate nucleus-IOD genome panel from the
larger strict linked-pair candidate table using QC/focus/linkage criteria that
do not rank by nucleus darkness, cell size, nucleus size, or nucleus IOD.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[2]
CELLPROFILER_ROOT = PROJECT_ROOT.parent / "cellprofiler_test"
RUN_ROOT = CELLPROFILER_ROOT / "output" / "runs" / "mixed_cellpose_yolo_full_dataset_v1_bgclean"
DEFAULT_CANDIDATES_CSV = RUN_ROOT / "pair_review_postrepair_strict" / "triage" / "all_scored_masks.csv"
DEFAULT_CURRENT_TOP50_CSV = (
    RUN_ROOT
    / "verified_species_dataset_top50_latest"
    / "selected_high_quality_linked_pairs_with_iod_qc.csv.gz"
)
DEFAULT_OUTPUT_DIR = PROJECT_ROOT / "path_analysis" / "data" / "external" / "derived" / "balanced_genome_iod_sensitivity"
DEFAULT_VIEWER_DECISIONS_CSV = (
    DEFAULT_OUTPUT_DIR
    / "viewer_balanced_qc_only"
    / "recovered_localstorage_decisions.csv"
)
DEFAULT_REFERENCE_SPECIES = "D. fuscus"
DEFAULT_REFERENCE_GENOME_PG = 16.36

BOOL_COLUMNS = [
    "has_cell_match",
    "physical_pair_ok",
    "one_to_one_cell",
    "keep_mask_pair",
    "keep_strict_core",
    "cell_edge_touch",
]

NUMERIC_COLUMNS = [
    "cell_area_um2",
    "nuc_area_um2",
    "nuc_iod",
    "nuc_mean_od",
    "nucleus_mask_overlap_label_count",
    "nucleus_mask_best_overlap_fraction",
    "distance_over_cell_radius",
    "raw_focus_lap_var",
    "raw_focus_grad_mean",
    "raw_focus_intensity_iqr",
    "cell_shape_smoothness",
    "cell_shape_ellipse_iou",
    "cell_shape_solidity",
    "nucleus_shape_smoothness",
    "nucleus_shape_ellipse_iou",
    "nucleus_shape_solidity",
    "shape_is_suspect",
]

SELECTED_PAIR_COLUMNS = [
    "panel",
    "species",
    "filename",
    "specimen_group",
    "review_key",
    "tile_name",
    "mask_label_id",
    "nucleus_label",
    "nuc_iod",
    "nuc_mean_od",
    "nuc_area_um2",
    "cell_area_um2",
    "nc_area_ratio",
    "balanced_genome_qc_score",
    "clarity_signal",
    "shape_signal",
    "pair_distance_signal",
    "overlap_signal",
    "quality_floor_pass",
    "quality_floor_relaxed",
    "selection_rank",
    "within_image_rank",
    "nucleus_mask_best_overlap_fraction",
    "nucleus_mask_overlap_label_count",
    "distance_over_cell_radius",
    "raw_focus_lap_var",
    "raw_focus_grad_mean",
    "raw_focus_intensity_iqr",
    "cell_shape_smoothness",
    "cell_shape_ellipse_iou",
    "cell_shape_solidity",
    "nucleus_shape_smoothness",
    "nucleus_shape_ellipse_iou",
    "nucleus_shape_solidity",
    "nucleus_source_image_path",
    "nucleus_mask_path",
    "cell_source_image_path",
    "cell_mask_path",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidates-csv", type=Path, default=DEFAULT_CANDIDATES_CSV)
    parser.add_argument("--current-top50-csv", type=Path, default=DEFAULT_CURRENT_TOP50_CSV)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--target-per-species", type=int, default=50)
    parser.add_argument("--min-clarity-percentile", type=float, default=0.25)
    parser.add_argument("--reference-species", default=DEFAULT_REFERENCE_SPECIES)
    parser.add_argument("--reference-genome-pg", type=float, default=DEFAULT_REFERENCE_GENOME_PG)
    parser.add_argument("--max-overlap-label-count", type=int, default=1)
    parser.add_argument("--min-overlap-fraction", type=float, default=0.98)
    parser.add_argument(
        "--viewer-decisions-csv",
        type=Path,
        default=DEFAULT_VIEWER_DECISIONS_CSV,
        help="Optional exported/recovered viewer decisions. Problem rows are excluded from the curated panel.",
    )
    return parser.parse_args()


def require_exists(path: Path) -> None:
    if not path.exists():
        raise FileNotFoundError(path)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def normalize_species(value: Any) -> str:
    text = "" if value is None else str(value).strip()
    if not text or text.lower() == "nan":
        return ""
    if text.startswith("Desmognathus "):
        return "D. " + text.split(" ", 1)[1]
    if text.startswith("D. "):
        return text
    return "D. " + text


def parse_bool_series(series: pd.Series) -> pd.Series:
    if pd.api.types.is_bool_dtype(series):
        return series.fillna(False)
    numeric = pd.to_numeric(series, errors="coerce")
    text = series.fillna("").astype(str).str.strip().str.lower()
    return numeric.fillna(0).ne(0) | text.isin({"true", "t", "yes", "y"})


def as_clean_group(series: pd.Series, fallback: pd.Series) -> pd.Series:
    values = series.fillna("").astype(str).str.strip()
    fallback_values = fallback.fillna("").astype(str).str.strip()
    return values.where(values.str.len().gt(0), fallback_values)


def weighted_quantile(values: pd.Series, weights: pd.Series, quantile: float) -> float:
    vals = pd.to_numeric(values, errors="coerce").to_numpy(dtype=float)
    wts = pd.to_numeric(weights, errors="coerce").to_numpy(dtype=float)
    mask = np.isfinite(vals) & np.isfinite(wts) & (wts > 0)
    if not mask.any():
        return float("nan")
    vals = vals[mask]
    wts = wts[mask]
    order = np.argsort(vals)
    vals = vals[order]
    wts = wts[order]
    cumulative = np.cumsum(wts)
    cutoff = quantile * cumulative[-1]
    return float(vals[np.searchsorted(cumulative, cutoff, side="left")])


def effective_sample_size(weights: pd.Series) -> float:
    wts = pd.to_numeric(weights, errors="coerce").fillna(0.0).to_numpy(dtype=float)
    denom = float(np.sum(wts * wts))
    if denom <= 0:
        return 0.0
    total = float(np.sum(wts))
    return (total * total) / denom


def percentile_by_species(values: pd.Series, species: pd.Series) -> pd.Series:
    numeric = pd.to_numeric(values, errors="coerce")
    ranked = numeric.groupby(species).rank(method="average", pct=True)
    med = ranked.groupby(species).transform("median")
    return ranked.fillna(med).fillna(0.5).clip(0.0, 1.0)


def scaled_quality(values: pd.Series, species: pd.Series, *, higher_is_better: bool = True) -> pd.Series:
    ranked = percentile_by_species(values, species)
    if higher_is_better:
        return ranked
    return (1.0 - ranked).clip(0.0, 1.0)


def add_qc_scores(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    species = out["species"]
    out["clarity_signal"] = (
        0.45 * scaled_quality(out["raw_focus_intensity_iqr"], species)
        + 0.35 * scaled_quality(out["raw_focus_grad_mean"], species)
        + 0.20 * scaled_quality(out["raw_focus_lap_var"], species)
    )
    shape_parts = [
        scaled_quality(out[col], species)
        for col in [
            "cell_shape_smoothness",
            "cell_shape_ellipse_iou",
            "cell_shape_solidity",
            "nucleus_shape_smoothness",
            "nucleus_shape_ellipse_iou",
            "nucleus_shape_solidity",
        ]
        if col in out.columns
    ]
    out["shape_signal"] = sum(shape_parts) / len(shape_parts) if shape_parts else 0.5
    out["pair_distance_signal"] = scaled_quality(out["distance_over_cell_radius"], species, higher_is_better=False)
    out["overlap_signal"] = scaled_quality(out["nucleus_mask_best_overlap_fraction"], species)
    out["balanced_genome_qc_score"] = (
        0.45 * out["clarity_signal"]
        + 0.35 * out["shape_signal"]
        + 0.15 * out["pair_distance_signal"]
        + 0.05 * out["overlap_signal"]
    )
    return out


def prepare_candidate_frame(path: Path, *, max_overlap_label_count: int, min_overlap_fraction: float) -> pd.DataFrame:
    require_exists(path)
    df = pd.read_csv(path, low_memory=False)
    required = {
        "species",
        "filename",
        "nuc_iod",
        "nuc_area_um2",
        "cell_area_um2",
        "review_key",
        "nucleus_mask_overlap_label_count",
        "nucleus_mask_best_overlap_fraction",
        *BOOL_COLUMNS,
    }
    missing = sorted(required - set(df.columns))
    if missing:
        raise ValueError(f"Missing required candidate columns: {missing}")

    out = df.copy()
    out["species"] = out["species"].map(normalize_species)
    for col in BOOL_COLUMNS:
        out[col] = parse_bool_series(out[col])
    for col in NUMERIC_COLUMNS:
        if col in out.columns:
            out[col] = pd.to_numeric(out[col], errors="coerce")
        else:
            out[col] = np.nan
    if "decision" in out.columns:
        decision = out["decision"].fillna("").astype(str).str.strip().str.lower()
    else:
        decision = pd.Series("", index=out.index)
    out["specimen_group"] = as_clean_group(out.get("specimen_id", pd.Series("", index=out.index)), out["filename"])
    out["image_group"] = out["filename"].fillna("").astype(str).str.strip()

    strict = (
        out["species"].str.len().gt(0)
        & out["filename"].fillna("").astype(str).str.len().gt(0)
        & out["has_cell_match"]
        & out["physical_pair_ok"]
        & out["one_to_one_cell"]
        & out["keep_mask_pair"]
        & out["keep_strict_core"]
        & ~out["cell_edge_touch"]
        & out["cell_area_um2"].gt(0)
        & out["nuc_area_um2"].gt(0)
        & out["nuc_iod"].gt(0)
        & out["nucleus_mask_overlap_label_count"].le(max_overlap_label_count)
        & out["nucleus_mask_best_overlap_fraction"].ge(min_overlap_fraction)
        & ~decision.isin({"discard", "nucleus_only", "maybe"})
    )
    if "shape_is_suspect" in out.columns:
        strict &= ~out["shape_is_suspect"].fillna(0).astype(float).astype(bool)
    out = out.loc[strict].copy()
    if out.empty:
        raise ValueError("No candidate rows passed strict balanced genome-IOD eligibility.")
    return add_qc_scores(out)


def apply_quality_floor(df: pd.DataFrame, *, target_per_species: int, min_clarity_percentile: float) -> pd.DataFrame:
    out = df.copy()
    out["quality_floor_pass"] = out["clarity_signal"].ge(min_clarity_percentile)
    pass_counts = out.loc[out["quality_floor_pass"]].groupby("species").size()
    total_counts = out.groupby("species").size()
    shortfall_species = set(total_counts.index[pass_counts.reindex(total_counts.index, fill_value=0).lt(target_per_species)])
    out["quality_floor_relaxed"] = out["species"].isin(shortfall_species) & ~out["quality_floor_pass"]
    return out.loc[out["quality_floor_pass"] | out["species"].isin(shortfall_species)].copy()


def select_balanced_by_image(df: pd.DataFrame, *, target_per_species: int, panel: str) -> pd.DataFrame:
    selected_groups: list[pd.DataFrame] = []
    df = df.copy()
    if "image_group" not in df.columns:
        df["image_group"] = df["filename"].fillna("").astype(str).str.strip()
    if "quality_floor_pass" not in df.columns:
        df["quality_floor_pass"] = True
    sort_cols = [
        "species",
        "image_group",
        "balanced_genome_qc_score",
        "clarity_signal",
        "shape_signal",
        "pair_distance_signal",
        "overlap_signal",
        "review_key",
    ]
    ascending = [True, True, False, False, False, False, False, True]
    ranked = df.sort_values(sort_cols, ascending=ascending, kind="mergesort").copy()
    ranked["within_image_rank"] = ranked.groupby(["species", "image_group"], sort=False).cumcount() + 1
    for species, group in ranked.groupby("species", sort=True):
        ordered = group.sort_values(
            [
                "within_image_rank",
                "quality_floor_pass",
                "balanced_genome_qc_score",
                "clarity_signal",
                "shape_signal",
                "review_key",
            ],
            ascending=[True, False, False, False, False, True],
            kind="mergesort",
        ).head(target_per_species)
        ordered = ordered.copy()
        ordered["panel"] = panel
        ordered["selection_rank"] = np.arange(1, len(ordered) + 1)
        selected_groups.append(ordered)
    selected = pd.concat(selected_groups, ignore_index=True) if selected_groups else pd.DataFrame()
    return selected


def load_problem_review_keys(path: Path | None) -> pd.DataFrame:
    if path is None or not path.exists():
        return pd.DataFrame(columns=["species", "review_key", "decision"])
    decisions = pd.read_csv(path, low_memory=False)
    required = {"species", "review_key", "decision"}
    missing = sorted(required - set(decisions.columns))
    if missing:
        raise ValueError(f"Missing required viewer decision columns: {missing}")
    out = decisions.copy()
    out["species"] = out["species"].map(normalize_species)
    out["review_key"] = out["review_key"].fillna("").astype(str).str.strip()
    out["decision"] = out["decision"].fillna("").astype(str).str.strip().str.lower()
    out = out.loc[out["decision"].eq("problem") & out["species"].str.len().gt(0) & out["review_key"].str.len().gt(0)]
    return out[["species", "review_key", "decision"]].drop_duplicates().reset_index(drop=True)


def exclude_problem_review_keys(df: pd.DataFrame, problem_rows: pd.DataFrame) -> pd.DataFrame:
    if problem_rows.empty:
        return df.copy()
    problem_pairs = set(zip(problem_rows["species"], problem_rows["review_key"]))
    mask = [
        (species, review_key) not in problem_pairs
        for species, review_key in zip(
            df["species"].fillna("").astype(str),
            df["review_key"].fillna("").astype(str).str.strip(),
        )
    ]
    return df.loc[mask].copy()


def prepare_current_top50(path: Path) -> pd.DataFrame:
    require_exists(path)
    out = pd.read_csv(path, low_memory=False)
    out["species"] = out["species"].map(normalize_species)
    out["specimen_group"] = as_clean_group(out.get("specimen_id", pd.Series("", index=out.index)), out["filename"])
    out["image_group"] = out["filename"].fillna("").astype(str).str.strip()
    for col in NUMERIC_COLUMNS:
        if col in out.columns:
            out[col] = pd.to_numeric(out[col], errors="coerce")
    for col in ["balanced_genome_qc_score", "clarity_signal", "shape_signal", "pair_distance_signal", "overlap_signal"]:
        if col not in out.columns:
            out[col] = np.nan
    out["quality_floor_pass"] = pd.NA
    out["quality_floor_relaxed"] = pd.NA
    out["within_image_rank"] = out.groupby(["species", "image_group"], sort=False).cumcount() + 1
    out["selection_rank"] = out.groupby("species", sort=False).cumcount() + 1
    out["panel"] = "current_top50_all_selected"
    return out


def current_top50_primary_subset(current_top50: pd.DataFrame) -> pd.DataFrame:
    """Match the current primary genome rule: image-QC-pass rows when present, else all rows."""
    frames: list[pd.DataFrame] = []
    for _species, group in current_top50.groupby("species", sort=True):
        if "image_iod_qc_pass" in group.columns:
            qc_pass = parse_bool_series(group["image_iod_qc_pass"])
            chosen = group.loc[qc_pass].copy()
            if chosen.empty:
                chosen = group.copy()
        else:
            chosen = group.copy()
        chosen["panel"] = "current_top50_primary"
        chosen["selection_rank"] = np.arange(1, len(chosen) + 1)
        frames.append(chosen)
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


def assign_analysis_weights(df: pd.DataFrame) -> pd.Series:
    weights = pd.Series(0.0, index=df.index, dtype=float)
    if df.empty:
        return weights
    specimens = list(dict.fromkeys(df["specimen_group"].fillna("").astype(str)))
    n_specimens = max(len(specimens), 1)
    for specimen in specimens:
        spec_mask = df["specimen_group"].fillna("").astype(str).eq(specimen)
        spec_df = df.loc[spec_mask]
        images = list(dict.fromkeys(spec_df["image_group"].fillna("").astype(str)))
        n_images = max(len(images), 1)
        for image in images:
            img_mask = spec_mask & df["image_group"].fillna("").astype(str).eq(image)
            n_rows = max(int(img_mask.sum()), 1)
            weights.loc[img_mask] = (1.0 / n_specimens) * (1.0 / n_images) * (1.0 / n_rows)
    total = float(weights.sum())
    if total > 0:
        weights /= total
    return weights


def support_label(n_pairs: int, n_images: int, n_specimens: int, eff_n: float) -> str:
    if n_pairs >= 50 and n_images >= 2 and n_specimens >= 2 and eff_n >= 25:
        return "medium"
    if n_pairs >= 25 and n_images >= 2 and n_specimens >= 2:
        return "limited"
    return "low"


def support_warnings(n_pairs: int, n_images: int, n_specimens: int, panel: str) -> str:
    warnings: list[str] = []
    if n_pairs < 25:
        warnings.append("low_pair_count")
    if n_images < 2:
        warnings.append("low_image_count")
    if n_specimens < 2:
        warnings.append("low_specimen_count")
    if panel.startswith("current_top50"):
        warnings.append("size_darkness_selected_top50")
    return "; ".join(warnings)


def optional_numeric_median(df: pd.DataFrame, column: str) -> float:
    if column not in df.columns:
        return float("nan")
    return float(pd.to_numeric(df[column], errors="coerce").median())


def summarize_panel(panel_df: pd.DataFrame, *, reference_species: str, reference_genome_pg: float) -> pd.DataFrame:
    if panel_df.empty:
        return pd.DataFrame()
    reference = panel_df.loc[panel_df["species"].eq(reference_species)].copy()
    if reference.empty:
        raise ValueError(f"Panel {panel_df['panel'].iloc[0]!r} has no rows for reference species {reference_species!r}.")
    ref_weights = assign_analysis_weights(reference)
    reference_iod = weighted_quantile(reference["nuc_iod"], ref_weights, 0.5)
    if not np.isfinite(reference_iod) or reference_iod <= 0:
        raise ValueError(f"Panel {panel_df['panel'].iloc[0]!r} has invalid reference IOD.")
    scale = float(reference_genome_pg / reference_iod)
    rows: list[dict[str, Any]] = []
    for species, group in panel_df.groupby("species", sort=True):
        weights = assign_analysis_weights(group)
        estimated_pg = pd.to_numeric(group["nuc_iod"], errors="coerce") * scale
        n_pairs = int(len(group))
        n_images = int(group["image_group"].nunique())
        n_specimens = int(group["specimen_group"].nunique())
        eff_n = effective_sample_size(weights)
        panel = str(group["panel"].iloc[0])
        rows.append(
            {
                "panel": panel,
                "species": species,
                "n_selected_pairs": n_pairs,
                "n_selected_images": n_images,
                "n_selected_specimens": n_specimens,
                "effective_n": eff_n,
                "weighted_median_nuc_iod": weighted_quantile(group["nuc_iod"], weights, 0.5),
                "nuc_iod_q1": weighted_quantile(group["nuc_iod"], weights, 0.25),
                "nuc_iod_q3": weighted_quantile(group["nuc_iod"], weights, 0.75),
                "estimated_genome_pg": weighted_quantile(estimated_pg, weights, 0.5),
                "estimated_genome_pg_q1": weighted_quantile(estimated_pg, weights, 0.25),
                "estimated_genome_pg_q3": weighted_quantile(estimated_pg, weights, 0.75),
                "reference_species": reference_species,
                "reference_genome_pg": reference_genome_pg,
                "reference_weighted_median_nuc_iod": reference_iod,
                "reference_scale_pg_per_iod": scale,
                "iod_ratio_to_reference": weighted_quantile(group["nuc_iod"], weights, 0.5) / reference_iod,
                "median_clarity_signal": optional_numeric_median(group, "clarity_signal"),
                "median_shape_signal": optional_numeric_median(group, "shape_signal"),
                "support_label": support_label(n_pairs, n_images, n_specimens, eff_n),
                "support_warnings": support_warnings(n_pairs, n_images, n_specimens, panel),
            }
        )
    return pd.DataFrame(rows)


def build_comparison(
    species_estimates: pd.DataFrame,
    *,
    baseline_panel: str = "current_top50_primary",
    comparison_panel: str = "balanced_qc_only",
) -> pd.DataFrame:
    current = species_estimates.loc[species_estimates["panel"].eq(baseline_panel)].copy()
    balanced = species_estimates.loc[species_estimates["panel"].eq(comparison_panel)].copy()
    baseline_suffix = "_" + baseline_panel
    balanced_suffix = "_" + comparison_panel
    merged = current.merge(balanced, on="species", suffixes=(baseline_suffix, balanced_suffix))
    merged[f"{comparison_panel}_minus_current_pg"] = (
        merged[f"estimated_genome_pg_{comparison_panel}"] - merged[f"estimated_genome_pg_{baseline_panel}"]
    )
    merged[f"{comparison_panel}_vs_current_pct"] = (
        merged[f"{comparison_panel}_minus_current_pg"] / merged[f"estimated_genome_pg_{baseline_panel}"] * 100.0
    )
    merged[f"{comparison_panel}_iod_ratio_minus_current"] = (
        merged[f"iod_ratio_to_reference_{comparison_panel}"] - merged[f"iod_ratio_to_reference_{baseline_panel}"]
    )
    keep = [
        "species",
        f"estimated_genome_pg_{baseline_panel}",
        f"estimated_genome_pg_{comparison_panel}",
        f"{comparison_panel}_minus_current_pg",
        f"{comparison_panel}_vs_current_pct",
        f"iod_ratio_to_reference_{baseline_panel}",
        f"iod_ratio_to_reference_{comparison_panel}",
        f"{comparison_panel}_iod_ratio_minus_current",
        f"n_selected_pairs_{baseline_panel}",
        f"n_selected_pairs_{comparison_panel}",
        f"n_selected_images_{baseline_panel}",
        f"n_selected_images_{comparison_panel}",
        f"n_selected_specimens_{baseline_panel}",
        f"n_selected_specimens_{comparison_panel}",
        f"support_label_{baseline_panel}",
        f"support_label_{comparison_panel}",
        f"support_warnings_{baseline_panel}",
        f"support_warnings_{comparison_panel}",
    ]
    out = merged[keep].sort_values(f"{comparison_panel}_vs_current_pct").reset_index(drop=True)
    if comparison_panel == "balanced_qc_only":
        return out.rename(
            columns={
                "balanced_qc_only_minus_current_pg": "balanced_minus_current_pg",
                "balanced_qc_only_vs_current_pct": "balanced_vs_current_pct",
                "balanced_qc_only_iod_ratio_minus_current": "balanced_iod_ratio_minus_current",
            }
        )
    return out


def build_curated_comparison(species_estimates: pd.DataFrame) -> pd.DataFrame:
    original = species_estimates.loc[species_estimates["panel"].eq("balanced_qc_only")].copy()
    curated = species_estimates.loc[species_estimates["panel"].eq("balanced_qc_curated")].copy()
    merged = original.merge(curated, on="species", suffixes=("_balanced_qc_only", "_balanced_qc_curated"))
    merged["curated_minus_original_pg"] = (
        merged["estimated_genome_pg_balanced_qc_curated"] - merged["estimated_genome_pg_balanced_qc_only"]
    )
    merged["curated_vs_original_pct"] = (
        merged["curated_minus_original_pg"] / merged["estimated_genome_pg_balanced_qc_only"] * 100.0
    )
    merged["curated_iod_ratio_minus_original"] = (
        merged["iod_ratio_to_reference_balanced_qc_curated"] - merged["iod_ratio_to_reference_balanced_qc_only"]
    )
    keep = [
        "species",
        "estimated_genome_pg_balanced_qc_only",
        "estimated_genome_pg_balanced_qc_curated",
        "curated_minus_original_pg",
        "curated_vs_original_pct",
        "iod_ratio_to_reference_balanced_qc_only",
        "iod_ratio_to_reference_balanced_qc_curated",
        "curated_iod_ratio_minus_original",
        "n_selected_pairs_balanced_qc_only",
        "n_selected_pairs_balanced_qc_curated",
        "n_selected_images_balanced_qc_only",
        "n_selected_images_balanced_qc_curated",
        "n_selected_specimens_balanced_qc_only",
        "n_selected_specimens_balanced_qc_curated",
        "support_label_balanced_qc_only",
        "support_label_balanced_qc_curated",
        "support_warnings_balanced_qc_only",
        "support_warnings_balanced_qc_curated",
    ]
    return merged[keep].sort_values("curated_vs_original_pct").reset_index(drop=True)


def write_outputs(
    *,
    selected_pairs: pd.DataFrame,
    species_estimates: pd.DataFrame,
    comparison: pd.DataFrame,
    curated_comparison: pd.DataFrame,
    problem_rows: pd.DataFrame,
    output_dir: Path,
    args: argparse.Namespace,
    n_strict_candidates: int,
) -> dict[str, Any]:
    output_dir.mkdir(parents=True, exist_ok=True)
    selected_path = output_dir / "balanced_genome_iod_selected_pairs.csv.gz"
    estimates_path = output_dir / "balanced_genome_iod_species_estimates.csv"
    comparison_path = output_dir / "balanced_genome_iod_comparison.csv"
    curated_comparison_path = output_dir / "balanced_genome_iod_curated_comparison.csv"
    decisions_applied_path = output_dir / "viewer_problem_decisions_applied.csv"
    summary_path = output_dir / "summary.json"

    present_pair_cols = [col for col in SELECTED_PAIR_COLUMNS if col in selected_pairs.columns]
    selected_pairs[present_pair_cols].to_csv(selected_path, index=False)
    species_estimates.to_csv(estimates_path, index=False, float_format="%.6f")
    comparison.to_csv(comparison_path, index=False, float_format="%.6f")
    if not curated_comparison.empty:
        curated_comparison.to_csv(curated_comparison_path, index=False, float_format="%.6f")
    if not problem_rows.empty:
        problem_rows.to_csv(decisions_applied_path, index=False)

    summary = {
        "purpose": (
            "Additive apples-to-apples genome-IOD sensitivity panel. Balanced panels are selected by "
            "focus/shape/linkage QC, not nucleus darkness, cell size, nucleus size, or IOD."
        ),
        "comparison_baseline_panel": "current_top50_primary",
        "candidates_csv": str(args.candidates_csv.resolve()),
        "candidates_sha256": sha256_file(args.candidates_csv),
        "current_top50_csv": str(args.current_top50_csv.resolve()),
        "current_top50_sha256": sha256_file(args.current_top50_csv),
        "output_dir": str(output_dir.resolve()),
        "selected_pairs_csv": str(selected_path.resolve()),
        "species_estimates_csv": str(estimates_path.resolve()),
        "comparison_csv": str(comparison_path.resolve()),
        "curated_comparison_csv": str(curated_comparison_path.resolve()) if not curated_comparison.empty else None,
        "viewer_decisions_csv": str(args.viewer_decisions_csv.resolve()) if args.viewer_decisions_csv.exists() else None,
        "viewer_problem_decisions_applied_csv": str(decisions_applied_path.resolve()) if not problem_rows.empty else None,
        "n_viewer_problem_decisions_applied": int(len(problem_rows)),
        "viewer_problem_decisions_by_species": problem_rows.groupby("species").size().to_dict()
        if not problem_rows.empty
        else {},
        "n_strict_candidates": int(n_strict_candidates),
        "target_per_species": int(args.target_per_species),
        "min_clarity_percentile": float(args.min_clarity_percentile),
        "reference_species": args.reference_species,
        "reference_genome_pg": float(args.reference_genome_pg),
        "panels": species_estimates.groupby("panel").size().to_dict(),
        "largest_abs_balanced_vs_current_pct": comparison.assign(
            abs_pct=comparison["balanced_vs_current_pct"].abs()
        )
        .sort_values("abs_pct", ascending=False)
        .head(8)[["species", "balanced_vs_current_pct"]]
        .to_dict(orient="records"),
    }
    summary_path.write_text(json.dumps(summary, indent=2) + "\n")
    return summary


def build(args: argparse.Namespace) -> dict[str, Any]:
    candidates = prepare_candidate_frame(
        args.candidates_csv,
        max_overlap_label_count=args.max_overlap_label_count,
        min_overlap_fraction=args.min_overlap_fraction,
    )
    problem_rows = load_problem_review_keys(args.viewer_decisions_csv)
    qc_candidates = apply_quality_floor(
        candidates,
        target_per_species=args.target_per_species,
        min_clarity_percentile=args.min_clarity_percentile,
    )
    balanced_qc = select_balanced_by_image(qc_candidates, target_per_species=args.target_per_species, panel="balanced_qc_only")
    curated_qc = pd.DataFrame()
    if not problem_rows.empty:
        curated_candidates = exclude_problem_review_keys(qc_candidates, problem_rows)
        curated_qc = select_balanced_by_image(
            curated_candidates,
            target_per_species=args.target_per_species,
            panel="balanced_qc_curated",
        )
    balanced_all = select_balanced_by_image(candidates, target_per_species=args.target_per_species, panel="balanced_all_strict")
    current_all = prepare_current_top50(args.current_top50_csv)
    current_primary = current_top50_primary_subset(current_all)
    selected_frames = [current_primary, current_all, balanced_qc, balanced_all]
    if not curated_qc.empty:
        selected_frames.append(curated_qc)
    selected_pairs = pd.concat(selected_frames, ignore_index=True, sort=False)
    species_estimates = pd.concat(
        [
            summarize_panel(panel_df, reference_species=args.reference_species, reference_genome_pg=args.reference_genome_pg)
            for _panel, panel_df in selected_pairs.groupby("panel", sort=True)
        ],
        ignore_index=True,
    ).sort_values(["panel", "species"]).reset_index(drop=True)
    comparison = build_comparison(species_estimates, baseline_panel="current_top50_primary")
    curated_comparison = (
        build_curated_comparison(species_estimates)
        if "balanced_qc_curated" in set(species_estimates["panel"])
        else pd.DataFrame()
    )
    return write_outputs(
        selected_pairs=selected_pairs,
        species_estimates=species_estimates,
        comparison=comparison,
        curated_comparison=curated_comparison,
        problem_rows=problem_rows,
        output_dir=args.output_dir,
        args=args,
        n_strict_candidates=len(candidates),
    )


def main() -> None:
    args = parse_args()
    summary = build(args)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
