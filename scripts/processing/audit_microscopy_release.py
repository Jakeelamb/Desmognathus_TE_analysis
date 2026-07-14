#!/usr/bin/env python3
"""Build the final-panel microscopy and image-IOD release audit.

This audit preserves every historical result.  It derives a final-18 evidence
surface that keeps four questions separate:

1. Were the cell and nucleus segmentation models validated independently?
2. What erythrocyte-size estimand is represented by the frozen top-50 rows?
3. How sensitive are species ranks to the selection rule and image weighting?
4. Can nuclear IOD currently be interpreted as an absolute genome-size assay?

The corrected release deliberately reports nuclear IOD relative to the focal
fuscus images.  It never converts IOD to picograms because the repository does
not contain the calibration and wet-lab evidence needed for that claim.
"""

from __future__ import annotations

import hashlib
import json
import math
import subprocess
from itertools import combinations
from pathlib import Path
from typing import Iterable

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import spearmanr


PROJECT_ROOT = Path(__file__).resolve().parents[2]
CELL_REPO = PROJECT_ROOT.parent / "cellprofiler_test"
CELL_RUN = CELL_REPO / "output/runs/mixed_cellpose_yolo_full_dataset_v1_bgclean"

SELECTED = (
    CELL_RUN
    / "top50_linked_pair_review/final_curated_top50_latest/final_curated_top50_linked_pairs.csv"
)
CANDIDATES = (
    CELL_RUN
    / "top50_linked_pair_review/pair_quality_model/quality_scored_candidates.csv"
)
IOD_ROWS = (
    CELL_RUN
    / "verified_species_dataset_top50_latest/selected_high_quality_linked_pairs_with_iod_qc.csv.gz"
)
CELL_TRAINING = CELL_REPO / "output/tile_training_round_v1/training_manifest.csv"
NUCLEUS_TRAINING = CELL_REPO / "output/nucleus_label_manual_round1/manifest.csv"
CELL_MODEL = Path.home() / ".cellpose/models/cpsam"
CELL_RUN_PLAN = CELL_REPO / "output/runs/full_dataset_v1/logs/run_three_pipelines_plan.sh"
CELL_RUN_MANIFEST = CELL_REPO / "output/runs/full_dataset_v1/cell_size_segmentation/run_manifest.json"
CELL_RUN_LOG = CELL_REPO / "output/runs/full_dataset_v1/logs/cell_size.log"
NUCLEUS_MODEL = (
    CELL_REPO
    / "runs/segment/output/yolo_nucleus_training_final_area500_shape/weights/best.pt"
)
CELL_BENCHMARK = (
    CELL_REPO
    / "output/benchmarks/archived_production_cell_masks_tile_holdout_v1"
)
NUCLEUS_BENCHMARK = (
    CELL_REPO
    / "output/benchmarks/current_models_tile_holdout_v1_nucleus"
)

OUTPUT_DIR = PROJECT_ROOT / "results/data/corrected/microscopy"
FIGURE_DIR = PROJECT_ROOT / "results/figures/corrected/microscopy"
AUDIT_DIR = PROJECT_ROOT / "plans/publication-readiness-deep-audit"

SUPPORT_OUTPUT = OUTPUT_DIR / "microscopy_panel_support_analysis18_v1.csv"
ESTIMATOR_OUTPUT = OUTPUT_DIR / "microscopy_estimand_sensitivity_analysis18_v1.csv"
RANK_OUTPUT = OUTPUT_DIR / "microscopy_estimator_rank_stability_analysis18_v1.csv"
BOOTSTRAP_OUTPUT = OUTPUT_DIR / "microscopy_selected50_hierarchical_bootstrap_analysis18_v1.csv"
IOD_SPECIES_OUTPUT = OUTPUT_DIR / "microscopy_relative_iod_sensitivity_analysis18_v1.csv"
IOD_IMAGE_OUTPUT = OUTPUT_DIR / "microscopy_image_iod_qc_analysis18_v1.csv"
VALIDATION_OUTPUT = OUTPUT_DIR / "microscopy_segmentation_validation_analysis18_v1.csv"
SOURCE_OUTPUT = OUTPUT_DIR / "microscopy_source_registry_analysis18_v1.csv"
MANIFEST_OUTPUT = OUTPUT_DIR / "microscopy_release_audit_analysis18_v1.manifest.json"
REPORT_OUTPUT = AUDIT_DIR / "microscopy_release_audit_analysis18_v1.md"

ESTIMATOR_FIGURE = FIGURE_DIR / "microscopy_estimand_sensitivity_analysis18_v1.png"
SELECTION_FIGURE = FIGURE_DIR / "microscopy_selection_depth_analysis18_v1.png"
SUPPORT_FIGURE = FIGURE_DIR / "microscopy_support_and_review_analysis18_v1.png"
IOD_FIGURE = FIGURE_DIR / "microscopy_relative_iod_qc_sensitivity_analysis18_v1.png"
IOD_DECOMPOSITION_FIGURE = FIGURE_DIR / "microscopy_iod_decomposition_analysis18_v1.png"
IMAGE_QC_FIGURE = FIGURE_DIR / "microscopy_image_iod_qc_analysis18_v1.png"
VALIDATION_FIGURE = FIGURE_DIR / "microscopy_segmentation_validation_analysis18_v1.png"
RANK_FIGURE = FIGURE_DIR / "microscopy_estimator_rank_stability_analysis18_v1.png"

FINAL_SPECIES = [
    "amphileucus",
    "anicetus",
    "apalachicolae",
    "auriculatus",
    "bairdi",
    "campi",
    "fuscus",
    "gvnigeusgwotli",
    "intermedius",
    "kanawha",
    "marmoratus",
    "mavrokoilius",
    "monticola",
    "ocoee",
    "perlapsus",
    "tilleyi",
    "valtos",
    "welteri",
]

METRICS = {
    "cell_area_um2": "Cell area (um2)",
    "nuc_area_um2": "Nucleus area (um2)",
    "nc_area_ratio": "Nucleus/cell area ratio",
}

ESTIMATOR_ORDER = [
    "composite_selected50",
    "image_balanced_selected50",
    "literal_largest50_eligible",
    "literal_largest100_eligible",
    "eligible_candidate_pool",
    "manual_keep_pool",
]

ESTIMATOR_LABELS = {
    "composite_selected50": "Frozen composite top 50",
    "image_balanced_selected50": "Equal-image frozen top 50",
    "literal_largest50_eligible": "Literal largest 50",
    "literal_largest100_eligible": "Literal largest 100",
    "eligible_candidate_pool": "All eligible candidates",
    "manual_keep_pool": "Manual keeps",
}

PALETTE = {
    "composite_selected50": "#1B6CA8",
    "image_balanced_selected50": "#5AA6C8",
    "literal_largest50_eligible": "#C65D21",
    "literal_largest100_eligible": "#E39B45",
    "eligible_candidate_pool": "#607D6D",
    "manual_keep_pool": "#7A5195",
}


def canonical_species(value: object) -> str:
    text = str(value).strip()
    if text.lower().startswith("desmognathus "):
        text = text.split(" ", 1)[1]
    elif text.startswith("D. "):
        text = text[3:]
    elif text.startswith("D."):
        text = text[2:]
    return text.strip().lower()


def require_exists(paths: Iterable[Path]) -> None:
    missing = [str(path) for path in paths if not path.exists()]
    if missing:
        raise FileNotFoundError("Required microscopy audit inputs are missing: " + ", ".join(missing))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def relative_path(path: Path) -> str:
    resolved = path.resolve()
    try:
        return str(resolved.relative_to(PROJECT_ROOT.resolve()))
    except ValueError:
        return str(resolved)


def safe_spearman(left: pd.Series, right: pd.Series) -> tuple[float, float, int]:
    frame = pd.DataFrame({"left": left, "right": right}).replace([np.inf, -np.inf], np.nan).dropna()
    if len(frame) < 3 or frame["left"].nunique() < 2 or frame["right"].nunique() < 2:
        return float("nan"), float("nan"), int(len(frame))
    result = spearmanr(frame["left"], frame["right"])
    return float(result.statistic), float(result.pvalue), int(len(frame))


def load_panel_frame(path: Path) -> pd.DataFrame:
    frame = pd.read_csv(path, low_memory=False)
    frame = frame.copy()
    frame["species"] = frame["species"].map(canonical_species)
    frame = frame[frame["species"].isin(FINAL_SPECIES)].copy()
    return frame


def eligible_candidates(frame: pd.DataFrame) -> pd.DataFrame:
    out = frame.copy()
    if "current_decision" in out:
        decisions = out["current_decision"].fillna("").astype(str).str.strip().str.lower()
        out = out[~decisions.isin({"discard", "nucleus_only"})].copy()
    return out


def image_balanced_point(group: pd.DataFrame, metric: str) -> float:
    image_points = group.groupby("filename", sort=False)[metric].median()
    return float(image_points.median()) if len(image_points) else float("nan")


def estimator_groups(
    species: str,
    selected: pd.DataFrame,
    candidates: pd.DataFrame,
) -> dict[str, pd.DataFrame]:
    selected_group = selected[selected["species"].eq(species)].copy()
    candidate_group = candidates[candidates["species"].eq(species)].copy()
    manual = candidate_group[
        candidate_group.get("current_decision", pd.Series(index=candidate_group.index, dtype=object))
        .fillna("")
        .astype(str)
        .str.lower()
        .eq("keep")
    ].copy()
    return {
        "composite_selected50": selected_group,
        "image_balanced_selected50": selected_group,
        "literal_largest50_eligible": candidate_group.nlargest(min(50, len(candidate_group)), "cell_area_um2"),
        "literal_largest100_eligible": candidate_group.nlargest(min(100, len(candidate_group)), "cell_area_um2"),
        "eligible_candidate_pool": candidate_group,
        "manual_keep_pool": manual,
    }


def build_estimators(selected: pd.DataFrame, candidates: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for species in FINAL_SPECIES:
        groups = estimator_groups(species, selected, candidates)
        for estimator, group in groups.items():
            for metric, metric_label in METRICS.items():
                point = (
                    image_balanced_point(group, metric)
                    if estimator == "image_balanced_selected50"
                    else float(group[metric].median()) if len(group) else float("nan")
                )
                rows.append(
                    {
                        "species": species,
                        "estimator": estimator,
                        "estimator_label": ESTIMATOR_LABELS[estimator],
                        "metric": metric,
                        "metric_label": metric_label,
                        "estimate": point,
                        "n_pairs": int(len(group)),
                        "n_images": int(group["filename"].nunique()) if len(group) else 0,
                        "n_specimens": int(group["specimen_id"].nunique()) if len(group) else 0,
                        "estimand": (
                            "equal-image median of within-image medians among frozen composite top-50 rows"
                            if estimator == "image_balanced_selected50"
                            else "median among " + ESTIMATOR_LABELS[estimator].lower()
                        ),
                    }
                )
    return pd.DataFrame(rows)


def build_support(selected: pd.DataFrame, candidates: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for species in FINAL_SPECIES:
        chosen = selected[selected["species"].eq(species)].copy()
        pool = candidates[candidates["species"].eq(species)].copy()
        if len(chosen) != 50:
            raise ValueError(f"{species} does not have exactly 50 frozen rows: {len(chosen)}")
        literal = pool.nlargest(min(50, len(pool)), "cell_area_um2")
        chosen_keys = set(chosen["review_key"].astype(str))
        literal_keys = set(literal["review_key"].astype(str))
        percentiles = [float((pool["cell_area_um2"] <= value).mean()) for value in chosen["cell_area_um2"]]
        tiers = chosen["grid_selection_tier"].fillna("unknown").astype(str)
        decisions = chosen["grid_effective_decision"].fillna("unlabeled").astype(str).str.lower()
        rows.append(
            {
                "species": species,
                "n_frozen_pairs": int(len(chosen)),
                "n_eligible_candidates": int(len(pool)),
                "n_images": int(chosen["filename"].nunique()),
                "n_specimens": int(chosen["specimen_id"].nunique()),
                "n_manual_keep": int(decisions.eq("keep").sum()),
                "n_manual_maybe": int(decisions.eq("maybe").sum()),
                "n_model_ranked_unlabeled": int(decisions.eq("unlabeled").sum()),
                "manual_keep_fraction": float(decisions.eq("keep").mean()),
                "n_tier_model_ranked": int(tiers.eq("model_ranked").sum()),
                "n_tier_reviewed_keep": int(tiers.eq("reviewed_keep").sum()),
                "literal_largest50_overlap_n": int(len(chosen_keys.intersection(literal_keys))),
                "literal_largest50_overlap_fraction": float(len(chosen_keys.intersection(literal_keys)) / 50),
                "selected_mean_cell_area_percentile": float(np.mean(percentiles)),
                "selected_min_cell_area_percentile": float(np.min(percentiles)),
                "selected_max_cell_area_percentile": float(np.max(percentiles)),
                "one_to_one_cell_fraction": float(pd.to_numeric(chosen["one_to_one_cell"], errors="coerce").mean()),
                "physical_pair_ok_fraction": float(pd.to_numeric(chosen["physical_pair_ok"], errors="coerce").mean()),
                "publication_morphology_status": "sensitivity_only",
                "publication_genome_size_status": "not_approved",
            }
        )
    return pd.DataFrame(rows)


def build_rank_stability(estimators: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for metric in METRICS:
        frame = estimators[estimators["metric"].eq(metric)].pivot(
            index="species", columns="estimator", values="estimate"
        )
        for left, right in combinations(ESTIMATOR_ORDER, 2):
            rho, pvalue, n_species = safe_spearman(frame.get(left), frame.get(right))
            paired = frame[[left, right]].dropna() if left in frame and right in frame else pd.DataFrame()
            relative = (
                100.0 * (paired[right] - paired[left]) / paired[left]
                if len(paired)
                else pd.Series(dtype=float)
            )
            rows.append(
                {
                    "metric": metric,
                    "estimator_left": left,
                    "estimator_right": right,
                    "spearman_rho": rho,
                    "spearman_p": pvalue,
                    "n_species": n_species,
                    "median_pct_shift_right_vs_left": float(relative.median()) if len(relative) else np.nan,
                    "maximum_absolute_pct_shift": float(relative.abs().max()) if len(relative) else np.nan,
                }
            )
    return pd.DataFrame(rows)


def hierarchical_median_draws(
    group: pd.DataFrame,
    metric: str,
    rng: np.random.Generator,
    repetitions: int,
) -> np.ndarray:
    specimen_groups = {
        str(specimen): values[metric].dropna().to_numpy(dtype=float)
        for specimen, values in group.groupby("specimen_id", sort=False)
        if values[metric].notna().any()
    }
    specimen_ids = list(specimen_groups)
    if not specimen_ids:
        return np.asarray([], dtype=float)
    draws = np.empty(repetitions, dtype=float)
    for offset in range(repetitions):
        sampled_specimens = rng.choice(specimen_ids, size=len(specimen_ids), replace=True)
        sampled_values: list[np.ndarray] = []
        for specimen in sampled_specimens:
            values = specimen_groups[str(specimen)]
            sampled_values.append(rng.choice(values, size=len(values), replace=True))
        draws[offset] = float(np.median(np.concatenate(sampled_values)))
    return draws


def build_bootstrap(selected: pd.DataFrame, repetitions: int = 2000, seed: int = 8301) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    rows: list[dict[str, object]] = []
    for species in FINAL_SPECIES:
        group = selected[selected["species"].eq(species)]
        for metric, label in METRICS.items():
            draws = hierarchical_median_draws(group, metric, rng, repetitions)
            rows.append(
                {
                    "species": species,
                    "metric": metric,
                    "metric_label": label,
                    "estimate": float(group[metric].median()),
                    "bootstrap_ci_low": float(np.quantile(draws, 0.025)),
                    "bootstrap_ci_high": float(np.quantile(draws, 0.975)),
                    "bootstrap_sd": float(draws.std(ddof=1)),
                    "bootstrap_repetitions": int(repetitions),
                    "n_pairs": int(len(group)),
                    "n_images": int(group["filename"].nunique()),
                    "n_specimens": int(group["specimen_id"].nunique()),
                    "uncertainty_scope": (
                        "conditional_object_only_single_specimen"
                        if group["specimen_id"].nunique() == 1
                        else "hierarchical_specimen_and_object_conditional_on_frozen_selection"
                    ),
                }
            )
    return pd.DataFrame(rows)


def build_iod(iod: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    image = (
        iod.groupby(["species", "filename"], as_index=False)
        .agg(
            specimen_id=("specimen_id", "first"),
            n_selected_nuclei=("species", "size"),
            image_iod_qc_status=("image_iod_qc_status", "first"),
            image_iod_qc_score=("image_iod_quality_score", "first"),
            median_nucleus_iod=("nuc_iod", "median"),
            median_nucleus_area_um2=("nuc_area_um2", "median"),
            median_nucleus_mean_od=("nuc_mean_od", "median"),
            image_background_intensity=("nuc_i_bg", "median"),
        )
        .sort_values(["species", "filename"])
        .reset_index(drop=True)
    )

    rows: list[dict[str, object]] = []
    for species in FINAL_SPECIES:
        group = iod[iod["species"].eq(species)].copy()
        passed = group[group["image_iod_qc_pass"].fillna(False)].copy()
        high = passed[passed["iod_final_quality_score"].ge(0.60)].copy()
        image_medians = group.groupby("filename")["nuc_iod"].median()
        all_iod = float(group["nuc_iod"].median())
        pass_iod = float(passed["nuc_iod"].median()) if len(passed) else np.nan
        high_iod = float(high["nuc_iod"].median()) if len(high) else np.nan
        identity = group["nuc_iod"] / (group["nuc_area_px"] * group["nuc_mean_od"])
        rows.append(
            {
                "species": species,
                "n_selected_nuclei": int(len(group)),
                "n_images": int(group["filename"].nunique()),
                "n_specimens": int(group["specimen_id"].nunique()),
                "n_pass_images": int(passed["filename"].nunique()),
                "n_high_qc_nuclei": int(len(high)),
                "median_nucleus_iod_all": all_iod,
                "median_nucleus_iod_image_qc_pass": pass_iod,
                "median_nucleus_iod_high_qc": high_iod,
                "image_qc_pass_shift_pct": 100.0 * (pass_iod / all_iod - 1.0),
                "between_image_median_iod_cv_pct": (
                    100.0 * float(image_medians.std(ddof=1) / image_medians.mean())
                    if len(image_medians) > 1
                    else np.nan
                ),
                "median_nucleus_area_um2": float(group["nuc_area_um2"].median()),
                "median_nucleus_mean_od": float(group["nuc_mean_od"].median()),
                "iod_area_times_mean_od_max_abs_error": float(np.nanmax(np.abs(identity - 1.0))),
                "iod_is_algebraically_area_times_mean_od": bool(np.nanmax(np.abs(identity - 1.0)) < 1e-10),
                "absolute_genome_size_status": "not_approved",
                "approved_release_variable": "relative_nuclear_iod_index_sensitivity_only",
            }
        )
    species = pd.DataFrame(rows)
    fuscus = species.set_index("species").loc["fuscus"]
    for suffix in ["all", "image_qc_pass", "high_qc"]:
        column = f"median_nucleus_iod_{suffix}"
        species[f"relative_iod_to_fuscus_{suffix}"] = species[column] / float(fuscus[column])
    return species, image


def manifest_validation(
    task: str,
    manifest_path: Path,
    benchmark_dir: Path,
    *,
    production_model_uses_local_training: bool,
) -> dict[str, object]:
    manifest = pd.read_csv(manifest_path)
    manifest["species"] = manifest["species"].map(canonical_species)
    train = manifest[manifest["split"].eq("train")]
    test = manifest[manifest["split"].eq("test")]
    train_species = sorted(set(train["species"]))
    test_species = sorted(set(test["species"]))
    train_images = sorted(set(train["filename"].astype(str)))
    test_images = sorted(set(test["filename"].astype(str)))

    summary_path = benchmark_dir / "summary.csv"
    if summary_path.exists():
        summaries = pd.read_csv(summary_path)
        current = summaries[(summaries["task"].eq(task)) & (summaries["model"].eq("current"))]
        if current.empty:
            current = summaries[summaries["task"].eq(task)]
        metrics = current.iloc[0].to_dict() if len(current) else {}
    else:
        metrics = {}

    per_image_path = benchmark_dir / "per_image_metrics.csv"
    if per_image_path.exists():
        per_image = pd.read_csv(per_image_path)
        task_rows = per_image[per_image["task"].eq(task)]
        current_rows = task_rows[task_rows["model"].eq("current")]
        per_image = current_rows if len(current_rows) else task_rows
        total_gt_objects = float(pd.to_numeric(per_image.get("n_gt_objects"), errors="coerce").sum())
        total_pred_objects = float(pd.to_numeric(per_image.get("n_pred_objects"), errors="coerce").sum())
        total_gt_pixels = float(pd.to_numeric(per_image.get("gt_foreground_px"), errors="coerce").sum())
        total_pred_pixels = float(pd.to_numeric(per_image.get("pred_foreground_px"), errors="coerce").sum())
        object_count_bias_pct = 100.0 * (total_pred_objects / total_gt_objects - 1.0) if total_gt_objects else np.nan
        foreground_area_bias_pct = 100.0 * (total_pred_pixels / total_gt_pixels - 1.0) if total_gt_pixels else np.nan
        median_absolute_count_error = float(pd.to_numeric(per_image.get("count_error"), errors="coerce").abs().median())
    else:
        object_count_bias_pct = foreground_area_bias_pct = median_absolute_count_error = np.nan

    return {
        "task": task,
        "training_manifest": relative_path(manifest_path),
        "benchmark_directory": relative_path(benchmark_dir),
        "n_training_tiles_in_annotation_manifest": int(len(train)),
        "n_local_training_tiles_used_by_production_model": int(len(train)) if production_model_uses_local_training else 0,
        "n_test_tiles": int(len(test)),
        "n_training_objects_recorded_in_manifest": int(pd.to_numeric(train["n_masks"], errors="coerce").sum()),
        "n_test_objects_recorded_in_manifest": int(pd.to_numeric(test["n_masks"], errors="coerce").sum()),
        "training_species": ";".join(train_species),
        "test_species": ";".join(test_species),
        "n_training_species": int(len(train_species)),
        "n_test_species": int(len(test_species)),
        "source_image_overlap_train_test": ";".join(sorted(set(train_images).intersection(test_images))),
        "species_overlap_train_test": ";".join(sorted(set(train_species).intersection(test_species))),
        "production_model_uses_local_training": bool(production_model_uses_local_training),
        "production_model_training_source": (
            "local_manual_tiles"
            if production_model_uses_local_training
            else "external_Cellpose_cpsam_pretraining"
        ),
        "n_final18_species_in_training": int(len(set(train_species).intersection(FINAL_SPECIES))),
        "n_final18_species_in_test": int(len(set(test_species).intersection(FINAL_SPECIES))),
        "heldout_unit": "tile_only",
        "status": metrics.get("status", "missing"),
        "n_images_scored": metrics.get("n_images_scored", np.nan),
        "n_gt_objects": metrics.get("n_gt_objects", np.nan),
        "test_manifest_minus_scored_gt_objects": (
            int(pd.to_numeric(test["n_masks"], errors="coerce").sum())
            - float(metrics.get("n_gt_objects", np.nan))
        ),
        "n_pred_objects": metrics.get("n_pred_objects", np.nan),
        "foreground_iou": metrics.get("foreground_iou", np.nan),
        "foreground_dice": metrics.get("foreground_dice", np.nan),
        "precision_50": metrics.get("precision_50", np.nan),
        "recall_50": metrics.get("recall_50", np.nan),
        "f1_50": metrics.get("f1_50", np.nan),
        "precision_75": metrics.get("precision_75", np.nan),
        "recall_75": metrics.get("recall_75", np.nan),
        "f1_75": metrics.get("f1_75", np.nan),
        "object_count_bias_pct": object_count_bias_pct,
        "foreground_area_bias_pct": foreground_area_bias_pct,
        "median_absolute_count_error_per_tile": median_absolute_count_error,
        "publication_validation_status": "not_approved",
        "reason": (
            "Only two tile-level test images; source images and species overlap local training; "
            "no final-18 analysis species occurs in training or test ground truth."
            if production_model_uses_local_training
            else "Only two tile-level test images from two excluded species; no final-18 analysis "
            "species occurs in ground truth. Production cpsam did not use the local train split."
        ),
    }


def query_package_version(python_path: Path, package: str) -> str:
    if not python_path.exists():
        return "interpreter_missing"
    code = f"from importlib.metadata import version; print(version({package!r}))"
    try:
        return subprocess.check_output([str(python_path), "-c", code], text=True, timeout=30).strip()
    except (subprocess.SubprocessError, OSError):
        return "unavailable"


def build_sources() -> pd.DataFrame:
    paths = [
        ("frozen_top50_pairs", SELECTED),
        ("quality_scored_candidate_pool", CANDIDATES),
        ("selected_pairs_with_iod_qc", IOD_ROWS),
        ("cell_training_manifest", CELL_TRAINING),
        ("nucleus_training_manifest", NUCLEUS_TRAINING),
        ("cell_model", CELL_MODEL),
        ("production_cell_run_plan", CELL_RUN_PLAN),
        ("production_cell_run_manifest", CELL_RUN_MANIFEST),
        ("production_cell_run_log", CELL_RUN_LOG),
        ("nucleus_model", NUCLEUS_MODEL),
        ("cell_benchmark_summary", CELL_BENCHMARK / "summary.csv"),
        ("cell_benchmark_manifest", CELL_BENCHMARK / "benchmark_manifest.csv"),
        ("cell_benchmark_per_image", CELL_BENCHMARK / "per_image_metrics.csv"),
        ("nucleus_benchmark_summary", NUCLEUS_BENCHMARK / "summary.csv"),
        ("nucleus_benchmark_manifest", NUCLEUS_BENCHMARK / "benchmark_manifest.csv"),
        ("nucleus_benchmark_per_image", NUCLEUS_BENCHMARK / "per_image_metrics.csv"),
    ]
    rows = []
    for role, path in paths:
        rows.append(
            {
                "source_role": role,
                "path": relative_path(path),
                "sha256": sha256(path) if path.exists() else "missing",
                "size_bytes": int(path.stat().st_size) if path.exists() else np.nan,
            }
        )
    environment = {
        "cellpose_validation_runtime": query_package_version(CELL_REPO / ".venv/bin/python", "cellpose"),
        "yolo_validation_runtime": query_package_version(CELL_REPO / ".venv-yolo/bin/python", "ultralytics"),
    }
    for role, version in environment.items():
        rows.append({"source_role": role, "path": "installed environment", "sha256": "not_applicable", "size_bytes": np.nan, "version": version})
    return pd.DataFrame(rows)


def species_order(frame: pd.DataFrame, metric: str, estimator: str) -> list[str]:
    subset = frame[(frame["metric"].eq(metric)) & (frame["estimator"].eq(estimator))]
    return subset.sort_values("estimate")["species"].tolist()


def plot_estimators(estimators: pd.DataFrame) -> None:
    order = species_order(estimators, "cell_area_um2", "composite_selected50")
    fig, axes = plt.subplots(1, 2, figsize=(15, 9), constrained_layout=True)
    for axis, metric, title in zip(
        axes,
        ["cell_area_um2", "nuc_area_um2"],
        ["Cell-area estimand sensitivity", "Nucleus-area estimand sensitivity"],
    ):
        subset = estimators[
            estimators["metric"].eq(metric)
            & estimators["estimator"].isin(
                ["composite_selected50", "literal_largest50_eligible", "literal_largest100_eligible", "eligible_candidate_pool"]
            )
        ].copy()
        for offset, estimator in enumerate(
            ["eligible_candidate_pool", "literal_largest100_eligible", "literal_largest50_eligible", "composite_selected50"]
        ):
            points = subset[subset["estimator"].eq(estimator)].set_index("species").reindex(order)
            axis.scatter(
                points["estimate"],
                np.arange(len(order)) + (offset - 1.5) * 0.13,
                s=42,
                color=PALETTE[estimator],
                label=ESTIMATOR_LABELS[estimator],
                zorder=3,
            )
        axis.set_yticks(np.arange(len(order)), labels=[name.replace("gvnigeusgwotli", "gvnigeusgwotli") for name in order])
        axis.grid(axis="x", color="#D9DED9", linewidth=0.7)
        axis.set_title(title, loc="left", weight="bold")
        axis.set_xlabel(METRICS[metric])
        axis.set_ylabel("")
    axes[1].legend(frameon=False, loc="lower right")
    fig.suptitle(
        "Frozen composite top-50 values are neither pool medians nor literal largest-50 values",
        fontsize=15,
        weight="bold",
    )
    fig.savefig(ESTIMATOR_FIGURE, dpi=220)
    plt.close(fig)


def plot_selection(support: pd.DataFrame) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.5), constrained_layout=True)
    sns.regplot(
        data=support,
        x="n_eligible_candidates",
        y="selected_mean_cell_area_percentile",
        ax=axes[0],
        scatter_kws={"s": 55, "color": "#1B6CA8"},
        line_kws={"color": "#A44A3F"},
        ci=None,
    )
    rho, _, _ = safe_spearman(support["n_eligible_candidates"], support["selected_mean_cell_area_percentile"])
    for row in support.itertuples():
        if row.n_eligible_candidates in {support.n_eligible_candidates.min(), support.n_eligible_candidates.max()}:
            axes[0].annotate(row.species, (row.n_eligible_candidates, row.selected_mean_cell_area_percentile), xytext=(5, 4), textcoords="offset points", fontsize=8)
    axes[0].set_title(f"Selection depth changes with pool size (rho={rho:.3f})", loc="left", weight="bold")
    axes[0].set_xlabel("Eligible linked-pair pool size")
    axes[0].set_ylabel("Mean cell-area percentile of frozen rows")
    order = support.sort_values("literal_largest50_overlap_n")["species"]
    points = support.set_index("species").loc[order]
    axes[1].barh(order, points["literal_largest50_overlap_n"], color="#C65D21")
    axes[1].axvline(50, color="#333333", linewidth=1, linestyle="--")
    axes[1].set_xlim(0, 52)
    axes[1].set_xlabel("Frozen rows also in literal largest 50 (of 50)")
    axes[1].set_title("Only a minority of frozen rows are literal largest-50 cells", loc="left", weight="bold")
    fig.savefig(SELECTION_FIGURE, dpi=220)
    plt.close(fig)


def plot_support(support: pd.DataFrame) -> None:
    order = support.sort_values(["n_specimens", "manual_keep_fraction"])["species"]
    frame = support.set_index("species").loc[order]
    fig, axes = plt.subplots(1, 2, figsize=(14, 7.5), constrained_layout=True)
    axes[0].barh(order, frame["n_manual_keep"], color="#247B6B", label="Manual keep")
    axes[0].barh(order, frame["n_manual_maybe"], left=frame["n_manual_keep"], color="#D7A84B", label="Manual maybe")
    axes[0].barh(
        order,
        frame["n_model_ranked_unlabeled"],
        left=frame["n_manual_keep"] + frame["n_manual_maybe"],
        color="#A7AFB2",
        label="Model-ranked / unlabeled",
    )
    axes[0].set_xlim(0, 50)
    axes[0].set_xlabel("Frozen linked pairs")
    axes[0].set_title("Manual decision coverage", loc="left", weight="bold")
    axes[0].legend(frameon=False, loc="lower right")
    axes[1].barh(order, frame["n_specimens"], color="#466B7A", label="Specimens")
    axes[1].scatter(frame["n_images"], np.arange(len(order)), color="#D46A42", s=35, label="Images", zorder=3)
    axes[1].set_xlabel("Independent images / specimen IDs")
    axes[1].set_title("Biological and image support", loc="left", weight="bold")
    axes[1].legend(frameon=False, loc="lower right")
    fig.suptitle("Frozen top-50 support is uneven across species", fontsize=15, weight="bold")
    fig.savefig(SUPPORT_FIGURE, dpi=220)
    plt.close(fig)


def plot_iod(species: pd.DataFrame) -> None:
    order = species.sort_values("relative_iod_to_fuscus_all")["species"]
    frame = species.set_index("species").loc[order]
    fig, axis = plt.subplots(figsize=(10, 8), constrained_layout=True)
    y = np.arange(len(order))
    axis.scatter(frame["relative_iod_to_fuscus_all"], y - 0.15, color="#607D6D", s=48, label="All frozen nuclei")
    axis.scatter(frame["relative_iod_to_fuscus_image_qc_pass"], y, color="#1B6CA8", s=48, label="Image-QC pass")
    axis.scatter(frame["relative_iod_to_fuscus_high_qc"], y + 0.15, color="#C65D21", s=48, label="High-QC nuclei")
    axis.axvline(1.0, color="#333333", linestyle="--", linewidth=1)
    axis.set_yticks(y, labels=order)
    axis.set_xlabel("Relative nuclear IOD index (fuscus = 1 within subset)")
    axis.set_ylabel("")
    axis.grid(axis="x", color="#D9DED9", linewidth=0.7)
    axis.set_title("IOD sensitivity is relative—not an absolute C-value calibration", loc="left", weight="bold")
    axis.legend(frameon=False, loc="lower right")
    fig.savefig(IOD_FIGURE, dpi=220)
    plt.close(fig)


def plot_iod_decomposition(species: pd.DataFrame) -> None:
    rho_area, _, _ = safe_spearman(species["median_nucleus_iod_all"], species["median_nucleus_area_um2"])
    rho_od, _, _ = safe_spearman(species["median_nucleus_iod_all"], species["median_nucleus_mean_od"])
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.3), constrained_layout=True)
    scatter = axes[0].scatter(
        species["median_nucleus_area_um2"],
        species["median_nucleus_iod_all"],
        c=species["median_nucleus_mean_od"],
        cmap="viridis",
        s=75,
        edgecolor="white",
        linewidth=0.6,
    )
    for row in species.itertuples():
        if row.species in {"fuscus", "perlapsus", "valtos", "welteri"}:
            axes[0].annotate(row.species, (row.median_nucleus_area_um2, row.median_nucleus_iod_all), xytext=(4, 4), textcoords="offset points", fontsize=8)
    axes[0].set_xlabel("Median nucleus area (um2)")
    axes[0].set_ylabel("Median nuclear IOD")
    axes[0].set_title(f"IOD versus area (rho={rho_area:.3f})", loc="left", weight="bold")
    fig.colorbar(scatter, ax=axes[0], label="Median mean OD")
    axes[1].scatter(species["median_nucleus_mean_od"], species["median_nucleus_iod_all"], color="#7A5195", s=70)
    for row in species.itertuples():
        if row.species in {"fuscus", "perlapsus", "valtos", "welteri"}:
            axes[1].annotate(row.species, (row.median_nucleus_mean_od, row.median_nucleus_iod_all), xytext=(4, 4), textcoords="offset points", fontsize=8)
    axes[1].set_xlabel("Median nucleus mean OD")
    axes[1].set_ylabel("Median nuclear IOD")
    axes[1].set_title(f"IOD versus mean OD (rho={rho_od:.3f})", loc="left", weight="bold")
    fig.suptitle("Nuclear IOD is exactly pixel area x mean OD at object level", fontsize=15, weight="bold")
    fig.savefig(IOD_DECOMPOSITION_FIGURE, dpi=220)
    plt.close(fig)


def plot_image_qc(image: pd.DataFrame) -> None:
    order = FINAL_SPECIES
    max_images = int(image.groupby("species").size().max())
    matrix = np.full((len(order), max_images), np.nan)
    annotations = np.full((len(order), max_images), "", dtype=object)
    codes = {"fail": 0.0, "limited": 1.0, "pass": 2.0}
    for i, species in enumerate(order):
        rows = image[image["species"].eq(species)].sort_values("filename")
        for j, row in enumerate(rows.itertuples()):
            matrix[i, j] = codes.get(str(row.image_iod_qc_status), np.nan)
            annotations[i, j] = str(row.filename).replace("Process_", "P").replace("_raw_green.ome.tiff", "")
    fig, axis = plt.subplots(figsize=(9, 7.2), constrained_layout=True)
    cmap = sns.color_palette(["#B85450", "#D8AA48", "#3E8E74"], as_cmap=True)
    sns.heatmap(
        matrix,
        mask=np.isnan(matrix),
        cmap=cmap,
        vmin=0,
        vmax=2,
        annot=annotations,
        fmt="",
        cbar=False,
        linewidths=1,
        linecolor="white",
        yticklabels=order,
        xticklabels=[f"Image {i + 1}" for i in range(max_images)],
        ax=axis,
    )
    axis.set_title("IOD image-QC status (red fail, amber limited, green pass)", loc="left", weight="bold")
    axis.set_xlabel("")
    axis.set_ylabel("")
    fig.savefig(IMAGE_QC_FIGURE, dpi=220)
    plt.close(fig)


def plot_validation(validation: pd.DataFrame) -> None:
    metrics = ["foreground_iou", "f1_50", "f1_75"]
    labels = ["Foreground IoU", "Instance F1 @ 0.50", "Instance F1 @ 0.75"]
    fig, axes = plt.subplots(1, 2, figsize=(12, 5.3), constrained_layout=True)
    for axis, row in zip(axes, validation.itertuples()):
        values = [float(getattr(row, metric)) if pd.notna(getattr(row, metric)) else 0.0 for metric in metrics]
        bars = axis.bar(labels, values, color=["#466B7A", "#1B6CA8", "#C65D21"])
        axis.bar_label(bars, fmt="%.3f", padding=3)
        axis.set_ylim(0, 1.03)
        axis.tick_params(axis="x", rotation=20)
        axis.set_title(f"{row.task.capitalize()} model: two tile-level test images", loc="left", weight="bold")
        axis.text(
            0.02,
            0.04,
            (
                "0/18 analysis species in ground truth\nlocal training/test source images overlap"
                if row.production_model_uses_local_training
                else "0/18 analysis species in ground truth\nexternal pretrained production model"
            ),
            transform=axis.transAxes,
            color="#8E3B36",
            fontsize=9,
            va="bottom",
        )
    fig.suptitle("Current segmentation benchmark is diagnostic, not publication-grade validation", fontsize=15, weight="bold")
    fig.savefig(VALIDATION_FIGURE, dpi=220)
    plt.close(fig)


def plot_rank_stability(rank: pd.DataFrame) -> None:
    subset = rank[rank["metric"].isin(["cell_area_um2", "nuc_area_um2"])].copy()
    subset["comparison"] = subset["estimator_right"].map(ESTIMATOR_LABELS) + " vs " + subset["estimator_left"].map(ESTIMATOR_LABELS)
    matrix = subset.pivot(index="comparison", columns="metric", values="spearman_rho")
    matrix = matrix.rename(columns={"cell_area_um2": "Cell area", "nuc_area_um2": "Nucleus area"})
    fig, axis = plt.subplots(figsize=(8.5, 9), constrained_layout=True)
    sns.heatmap(matrix, annot=True, fmt=".3f", cmap="crest", vmin=0, vmax=1, linewidths=0.5, ax=axis)
    axis.set_title("Species-rank stability across microscopy estimators", loc="left", weight="bold")
    axis.set_xlabel("")
    axis.set_ylabel("")
    fig.savefig(RANK_FIGURE, dpi=220)
    plt.close(fig)


def report_text(
    support: pd.DataFrame,
    estimators: pd.DataFrame,
    rank: pd.DataFrame,
    iod_species: pd.DataFrame,
    validation: pd.DataFrame,
) -> str:
    total_manual = int(support["n_manual_keep"].sum())
    total_maybe = int(support["n_manual_maybe"].sum())
    total_unlabeled = int(support["n_model_ranked_unlabeled"].sum())
    one_specimen = support.loc[support["n_specimens"].eq(1), "species"].tolist()
    overlap_range = (int(support["literal_largest50_overlap_n"].min()), int(support["literal_largest50_overlap_n"].max()))
    pool_rho, _, _ = safe_spearman(support["n_eligible_candidates"], support["selected_mean_cell_area_percentile"])

    cell_rank = rank[
        rank["metric"].eq("cell_area_um2")
        & rank["estimator_left"].eq("composite_selected50")
        & rank["estimator_right"].eq("literal_largest100_eligible")
    ].iloc[0]
    nuc_rank = rank[
        rank["metric"].eq("nuc_area_um2")
        & rank["estimator_left"].eq("composite_selected50")
        & rank["estimator_right"].eq("literal_largest100_eligible")
    ].iloc[0]
    iod_area_rho, iod_area_p, _ = safe_spearman(iod_species["median_nucleus_iod_all"], iod_species["median_nucleus_area_um2"])
    iod_od_rho, iod_od_p, _ = safe_spearman(iod_species["median_nucleus_iod_all"], iod_species["median_nucleus_mean_od"])

    validation_lines = []
    for row in validation.itertuples():
        validation_lines.append(
            f"- **{row.task.capitalize()}:** foreground IoU {row.foreground_iou:.3f}, "
            f"instance F1@0.50 {row.f1_50:.3f}, F1@0.75 {row.f1_75:.3f}; "
            f"object-count bias {row.object_count_bias_pct:+.1f}%, foreground-area bias "
            f"{row.foreground_area_bias_pct:+.1f}%; {int(row.n_images_scored)} test tiles, "
            "0/18 focal species represented."
        )

    return f"""# Microscopy and image-IOD release audit — final 18 species, v1

## Release verdict

| Layer | Verdict | Approved interpretation |
| --- | --- | --- |
| Object trace and cell–nucleus linkage | **Approved with frozen-source provenance** | Exact source image, tile, cell mask, nucleus mask, object IDs, and one-to-one geometry are retained. |
| Segmentation accuracy | **Not approved as publication-grade validation** | The current two-tile check is diagnostic only. |
| Cell and nucleus morphology | **Sensitivity-only** | Median area of the 50 highest composite-ranked eligible linked erythrocytes per species; not typical cell size and not literal largest 50. |
| Nuclear IOD | **Sensitivity-only** | Relative nuclear IOD index, with fuscus set to 1 within each QC subset. |
| Absolute genome size from images | **Not approved** | No picogram/C-value claim enters the corrected release. |

This audit does not delete or overwrite the historical 21-species top-50 analysis.  It builds a fail-closed final-18 surface for the eventual path-analysis sensitivity suite.

## 1. Segmentation validation

The current manually labeled ground truth contains 12 large tiles (10 local-train, 2 test) from only *D. ochrophaeus* and *D. folkertsi*.  Neither species belongs to the focal 18-species analysis.  For the locally trained nucleus model, the two test tiles come from the same source images and species used for training; that split tests new tiles within two images, not new images or specimens.  The production cell model was externally pretrained `cpsam` and did not use the ten local training tiles, but its local evaluation is still limited to the same two test images and two excluded species.

The production run plan and March 7 run manifest show that the frozen cell masks were generated with the default Cellpose model (`cpsam`): no custom `--cellpose-model` argument was supplied.  The Desmognathus fine-tuned checkpoint was created on March 8 and is not the source of the frozen cell masks.  This audit scores the exact archived production masks, including their hashes and March 7 modification times.  The cached `cpsam` binary predates the run, but the original run did not store its model hash, so model-binary identity remains reconstructive rather than contemporaneously frozen.

{chr(10).join(validation_lines)}

The nucleus result is materially lower than the YOLO training dashboard's patch-validation mAP and is the more relevant full-tile diagnostic.  It still cannot estimate error in the focal panel.  The publication gate remains: annotate an image/specimen-held-out, taxonomically and technically stratified set; report mask overlap, boundary/area bias, split/merge rates, pairing error, and downstream species-summary bias.

The nucleus test manifest records 138 masks across the two rows, but the label arrays actually scored contain 243 instances.  Thus its per-row `n_masks` field is not a reliable nucleus-object count and must not be used as the validation denominator.

Figure: [`{relative_path(VALIDATION_FIGURE)}`](../../{relative_path(VALIDATION_FIGURE)})

## 2. What the frozen top 50 actually estimate

The frozen final-18 table has 900 linked pairs.  Of these, {total_manual} are explicit manual keeps, {total_maybe} are manual maybes, and {total_unlabeled} are model-ranked rows without an individual keep decision.  All frozen rows satisfy the recorded one-to-one/physical linkage gates, but a grid-selected model rank is not equivalent to an independent manual mask validation.

The ranking is composite: pair-quality score, model keep probability, the earlier top-50 score, cell area, and nucleus darkness all contribute.  Consequently:

- only {overlap_range[0]}–{overlap_range[1]} of each species' 50 frozen rows are also in its literal largest 50;
- the frozen rows occupy different cell-area quantiles across species;
- candidate-pool size correlates strongly with selection depth (Spearman rho = {pool_rho:.3f}); and
- three species have one contributing image/specimen: {', '.join(one_specimen)}.

The honest label is **median area of the 50 highest composite-ranked eligible linked erythrocytes**, or more briefly **quality-screened upper-tail morphology**.  It is not an unbiased estimate of species-average or species-median erythrocyte morphology.  A recent salamander hematology study sampled 50 erythrocytes randomly per individual, while the closest genome/cell/nucleus comparative study used 50 nuclei per individual and species medians across 1–4 individuals ([Liu et al. 2023](https://doi.org/10.7717/peerj.15446); [Mueller et al. 2008](https://doi.org/10.1016/j.zool.2007.07.010)).

Species ranks are fairly but not perfectly stable against a literal-largest-100 sensitivity: cell-area rho = {cell_rank.spearman_rho:.3f}; nucleus-area rho = {nuc_rank.spearman_rho:.3f}.  The full estimator table, not one preferred summary, must enter the downstream sensitivity analysis.

Figures:

- [`{relative_path(ESTIMATOR_FIGURE)}`](../../{relative_path(ESTIMATOR_FIGURE)})
- [`{relative_path(SELECTION_FIGURE)}`](../../{relative_path(SELECTION_FIGURE)})
- [`{relative_path(SUPPORT_FIGURE)}`](../../{relative_path(SUPPORT_FIGURE)})
- [`{relative_path(RANK_FIGURE)}`](../../{relative_path(RANK_FIGURE)})

## 3. Uncertainty and replication

The hierarchical bootstrap resamples specimens and then cells within specimen, but it remains conditional on the frozen selection and segmentation.  For a one-specimen species it is explicitly labeled `conditional_object_only_single_specimen`; it is not a population-level species interval.  None of these intervals include stain batch, acquisition, segmentation-model, reviewer, or selection-rule uncertainty.

## 4. Nuclear IOD is not currently an absolute genome-size assay

The pixel equation is implemented correctly: object IOD equals `sum(log10(I_bg / I_pixel))`, and in the exported object table `IOD = area_px × mean_OD` to numerical precision (maximum relative identity error {iod_species.iod_area_times_mean_od_max_abs_error.max():.2e}).  That is necessary but not sufficient for genome-size densitometry.

The repository does not provide a slide-linked record of fixation, acid hydrolysis, Schiff/Feulgen batch, staining time, co-processed DNA standard, microscope/camera settings, exposure/linearity/saturation tests, or the provenance and 1C/2C interpretation of the historical `fuscus = 16.36 pg` constant.  The exact methods benchmark requires controlled Feulgen staining, optical-density conversion, and a same-batch/same-slide standard ([Hardie et al. 2002](https://doi.org/10.1177/002215540205000601)).  The closest *Desmognathus* precedent used *Xenopus laevis* erythrocyte nuclei as an internal standard ([Sessions and Kezer 1986](https://doi.org/10.1007/BF00494802)).

The current [Animal Genome Size Database *D. fuscus* record](https://genomesize.com/result_species.php?id=553) also points to heterogeneous published Feulgen estimates rather than documenting the exact 16.36-pg constant.  Recovering a plausible species-level value would still not establish comparability among independently stained and imaged slides.

There is also built-in measurement dependence: species median IOD correlates with nucleus area (rho = {iod_area_rho:.3f}, p = {iod_area_p:.4g}) and mean OD (rho = {iod_od_rho:.3f}, p = {iod_od_p:.4g}), because IOD contains nuclear area algebraically.  A path `image-IOD genome size -> nucleus area` therefore cannot be treated as an independent measurement test.

The corrected release withdraws the picogram conversion and retains only relative-IOD sensitivity columns.  Absolute genome size and causal genome-to-nucleus claims require an independently calibrated assay or recovered same-batch Feulgen-standard evidence.

Figures:

- [`{relative_path(IOD_FIGURE)}`](../../{relative_path(IOD_FIGURE)})
- [`{relative_path(IOD_DECOMPOSITION_FIGURE)}`](../../{relative_path(IOD_DECOMPOSITION_FIGURE)})
- [`{relative_path(IMAGE_QC_FIGURE)}`](../../{relative_path(IMAGE_QC_FIGURE)})

## 5. Files for the path-analysis audit

- [`{relative_path(SUPPORT_OUTPUT)}`](../../{relative_path(SUPPORT_OUTPUT)}): review and specimen/image support.
- [`{relative_path(ESTIMATOR_OUTPUT)}`](../../{relative_path(ESTIMATOR_OUTPUT)}): morphology estimator sensitivity.
- [`{relative_path(RANK_OUTPUT)}`](../../{relative_path(RANK_OUTPUT)}): species-rank stability.
- [`{relative_path(BOOTSTRAP_OUTPUT)}`](../../{relative_path(BOOTSTRAP_OUTPUT)}): conditional hierarchical intervals.
- [`{relative_path(IOD_SPECIES_OUTPUT)}`](../../{relative_path(IOD_SPECIES_OUTPUT)}): relative-IOD sensitivity only.
- [`{relative_path(IOD_IMAGE_OUTPUT)}`](../../{relative_path(IOD_IMAGE_OUTPUT)}): image-level IOD QC.
- [`{relative_path(VALIDATION_OUTPUT)}`](../../{relative_path(VALIDATION_OUTPUT)}): segmentation validation and leakage audit.
- [`{relative_path(SOURCE_OUTPUT)}`](../../{relative_path(SOURCE_OUTPUT)}): immutable source hashes and runtime versions.

## Bottom line

The microscopy layer is now auditable but not wholly approved.  Linkage provenance is strong.  Morphology is usable only as a declared upper-tail sensitivity estimand with estimator and specimen-support sensitivity.  Segmentation generalization to the focal species remains unvalidated.  Image IOD is not an absolute genome-size measurement under the evidence currently stored in the repository.
"""


def main() -> None:
    require_exists(
        [
            SELECTED,
            CANDIDATES,
            IOD_ROWS,
            CELL_TRAINING,
            NUCLEUS_TRAINING,
            CELL_MODEL,
            CELL_RUN_PLAN,
            CELL_RUN_MANIFEST,
            CELL_RUN_LOG,
            NUCLEUS_MODEL,
            CELL_BENCHMARK / "summary.csv",
            NUCLEUS_BENCHMARK / "summary.csv",
        ]
    )
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    AUDIT_DIR.mkdir(parents=True, exist_ok=True)

    sns.set_theme(style="whitegrid", context="notebook")
    plt.rcParams.update({"font.family": "DejaVu Sans", "axes.spines.top": False, "axes.spines.right": False})

    selected = load_panel_frame(SELECTED)
    candidates = eligible_candidates(load_panel_frame(CANDIDATES))
    iod = load_panel_frame(IOD_ROWS)

    if set(selected["species"]) != set(FINAL_SPECIES):
        raise ValueError("Frozen microscopy source does not exactly cover the final 18 species")
    if len(selected) != 900:
        raise ValueError(f"Expected 900 final-panel frozen pairs, observed {len(selected)}")

    support = build_support(selected, candidates)
    estimators = build_estimators(selected, candidates)
    rank = build_rank_stability(estimators)
    bootstrap = build_bootstrap(selected)
    iod_species, iod_image = build_iod(iod)
    validation = pd.DataFrame(
        [
            manifest_validation(
                "cell",
                CELL_TRAINING,
                CELL_BENCHMARK,
                production_model_uses_local_training=False,
            ),
            manifest_validation(
                "nucleus",
                NUCLEUS_TRAINING,
                NUCLEUS_BENCHMARK,
                production_model_uses_local_training=True,
            ),
        ]
    )
    sources = build_sources()

    support.to_csv(SUPPORT_OUTPUT, index=False, float_format="%.9g")
    estimators.to_csv(ESTIMATOR_OUTPUT, index=False, float_format="%.9g")
    rank.to_csv(RANK_OUTPUT, index=False, float_format="%.9g")
    bootstrap.to_csv(BOOTSTRAP_OUTPUT, index=False, float_format="%.9g")
    iod_species.to_csv(IOD_SPECIES_OUTPUT, index=False, float_format="%.9g")
    iod_image.to_csv(IOD_IMAGE_OUTPUT, index=False, float_format="%.9g")
    validation.to_csv(VALIDATION_OUTPUT, index=False, float_format="%.9g")
    sources.to_csv(SOURCE_OUTPUT, index=False, float_format="%.9g")

    plot_estimators(estimators)
    plot_selection(support)
    plot_support(support)
    plot_iod(iod_species)
    plot_iod_decomposition(iod_species)
    plot_image_qc(iod_image)
    plot_validation(validation)
    plot_rank_stability(rank)

    REPORT_OUTPUT.write_text(report_text(support, estimators, rank, iod_species, validation))

    manifest = {
        "analysis_id": "microscopy_release_audit_analysis18_v1",
        "panel_species": FINAL_SPECIES,
        "n_species": len(FINAL_SPECIES),
        "n_frozen_pairs": int(len(selected)),
        "n_candidate_pairs_after_explicit_reject_exclusion": int(len(candidates)),
        "outputs": [
            relative_path(path)
            for path in [
                SUPPORT_OUTPUT,
                ESTIMATOR_OUTPUT,
                RANK_OUTPUT,
                BOOTSTRAP_OUTPUT,
                IOD_SPECIES_OUTPUT,
                IOD_IMAGE_OUTPUT,
                VALIDATION_OUTPUT,
                SOURCE_OUTPUT,
                REPORT_OUTPUT,
                ESTIMATOR_FIGURE,
                SELECTION_FIGURE,
                SUPPORT_FIGURE,
                IOD_FIGURE,
                IOD_DECOMPOSITION_FIGURE,
                IMAGE_QC_FIGURE,
                VALIDATION_FIGURE,
                RANK_FIGURE,
            ]
        ],
        "release_gates": {
            "object_trace_and_pair_linkage": "approved_with_frozen_provenance",
            "segmentation_accuracy": "not_approved_no_independent_focal_validation",
            "morphology": "sensitivity_only_quality_screened_upper_tail",
            "relative_nuclear_iod": "sensitivity_only",
            "absolute_genome_size_from_images": "not_approved",
        },
        "historical_outputs_removed": False,
    }
    MANIFEST_OUTPUT.write_text(json.dumps(manifest, indent=2) + "\n")

    print(f"Wrote {REPORT_OUTPUT}")
    print(f"Wrote {len(manifest['outputs'])} microscopy audit artifacts")


if __name__ == "__main__":
    main()
