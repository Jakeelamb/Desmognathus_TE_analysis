#!/usr/bin/env python3
"""Match nuclei across species on technical image quality before IOD analysis.

Matching uses contrast-normalized nuclear edge sharpness and relative local
background noise. It never uses IOD, nuclear darkness/contrast magnitude,
nucleus area, or cell area.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import sys
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist
from scipy.optimize import linear_sum_assignment
from scipy.stats import ks_2samp, spearmanr, wasserstein_distance


PROJECT_ROOT = Path(__file__).resolve().parents[2]
CELL_REPO = PROJECT_ROOT.parent / "cellprofiler_test"
sys.path.insert(0, str(PROJECT_ROOT / "path_analysis" / "scripts"))
sys.path.insert(0, str(CELL_REPO / "scripts"))


RUN_ROOT = CELL_REPO / "output" / "runs" / "mixed_cellpose_yolo_full_dataset_v1_bgclean"
DEFAULT_CANDIDATES = RUN_ROOT / "pair_review_postrepair_strict" / "triage" / "all_scored_masks.csv"
DEFAULT_OUTPUT_DIR = PROJECT_ROOT / "path_analysis" / "data" / "external" / "derived" / "image_quality_matched_genome_iod"
MATCH_FEATURES = [
    "match_log_edge_sharpness",
    "match_log_relative_ring_noise",
]
FORBIDDEN_MATCH_FEATURES = [
    "nuc_iod",
    "nuc_mean_od",
    "iod_inner_intensity_median",
    "iod_local_contrast",
    "iod_local_od_contrast",
    "nuc_area_um2",
    "cell_area_um2",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidates-csv", type=Path, default=DEFAULT_CANDIDATES)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--target-per-species", type=int, default=40)
    parser.add_argument("--force-pixel-qc", action="store_true")
    parser.add_argument(
        "--decisions-csv",
        type=Path,
        action="append",
        default=[],
        help="Review decision export. Repeat for replacement-review rounds.",
    )
    return parser.parse_args()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def add_matching_features(scored: pd.DataFrame) -> pd.DataFrame:
    out = scored.copy()
    edge = pd.to_numeric(out["iod_edge_sharpness"], errors="coerce").clip(lower=0)
    ring_noise = pd.to_numeric(out["iod_ring_noise_iqr"], errors="coerce").clip(lower=0)
    ring_level = pd.to_numeric(out["iod_ring_intensity_median"], errors="coerce").clip(lower=1e-6)
    out["match_log_edge_sharpness"] = np.log1p(edge)
    out["match_log_relative_ring_noise"] = np.log1p(ring_noise / ring_level)
    image_summary = (
        out.groupby(["species", "filename"], sort=False)
        .agg(
            match_image_median_log_edge_sharpness=("match_log_edge_sharpness", "median"),
            match_image_median_log_relative_ring_noise=("match_log_relative_ring_noise", "median"),
        )
        .reset_index()
    )
    out = out.merge(image_summary, on=["species", "filename"], how="left", validate="many_to_one")
    finite = np.ones(len(out), dtype=bool)
    for column in MATCH_FEATURES:
        finite &= np.isfinite(pd.to_numeric(out[column], errors="coerce").to_numpy(dtype=float))
    return out.loc[finite].copy()


def robust_standardize(frame: pd.DataFrame, feature_columns: list[str]) -> tuple[pd.DataFrame, dict[str, dict[str, float]]]:
    out = frame.copy()
    scaling: dict[str, dict[str, float]] = {}
    for feature in feature_columns:
        values = pd.to_numeric(out[feature], errors="coerce")
        low, high = values.quantile([0.01, 0.99])
        clipped = values.clip(lower=low, upper=high)
        median = float(clipped.median())
        q1, q3 = clipped.quantile([0.25, 0.75])
        scale = float(q3 - q1)
        if not np.isfinite(scale) or scale <= 0:
            scale = float(clipped.std(ddof=0))
        if not np.isfinite(scale) or scale <= 0:
            raise ValueError(f"Matching feature {feature!r} has no usable variation.")
        z_column = f"z_{feature}"
        out[z_column] = (clipped - median) / scale
        scaling[feature] = {"clip_low": float(low), "clip_high": float(high), "median": median, "scale_iqr": scale}
    return out, scaling


def quality_medoid_species(frame: pd.DataFrame, z_columns: list[str]) -> tuple[str, pd.DataFrame]:
    species = sorted(frame["species"].unique())
    distances = pd.DataFrame(0.0, index=species, columns=species)
    for left_idx, left in enumerate(species):
        for right in species[left_idx + 1 :]:
            distance = float(
                np.mean(
                    [
                        wasserstein_distance(frame.loc[frame["species"].eq(left), col], frame.loc[frame["species"].eq(right), col])
                        for col in z_columns
                    ]
                )
            )
            distances.loc[left, right] = distance
            distances.loc[right, left] = distance
    medoid = str(distances.sum(axis=1).sort_values(kind="mergesort").index[0])
    long = distances.rename_axis("species").reset_index().melt(id_vars="species", var_name="other_species", value_name="wasserstein_distance")
    return medoid, long


def image_balanced_lowest(group: pd.DataFrame, score_column: str, target: int) -> pd.DataFrame:
    queues = {
        str(image): image_rows.sort_values([score_column, "review_key"], kind="mergesort").index.tolist()
        for image, image_rows in group.groupby("filename", sort=True)
    }
    chosen: list[int] = []
    while len(chosen) < target and any(queues.values()):
        for image in sorted(queues):
            if queues[image] and len(chosen) < target:
                chosen.append(queues[image].pop(0))
    return group.loc[chosen].copy()


def build_common_template(frame: pd.DataFrame, *, medoid_species: str, z_columns: list[str], target: int) -> pd.DataFrame:
    medoid = frame.loc[frame["species"].eq(medoid_species)].copy()
    medoid_values = medoid[z_columns].to_numpy(dtype=float)
    nearest_by_species: list[np.ndarray] = []
    for species, group in frame.groupby("species", sort=True):
        if species == medoid_species:
            nearest_by_species.append(np.zeros(len(medoid), dtype=float))
            continue
        candidate_values = group[z_columns].to_numpy(dtype=float)
        nearest_by_species.append(np.sqrt(cdist(medoid_values, candidate_values, metric="sqeuclidean").min(axis=1)))
    medoid["common_support_worst_distance"] = np.max(np.vstack(nearest_by_species), axis=0)
    return image_balanced_lowest(medoid, "common_support_worst_distance", min(target, len(medoid)))


def optimal_quality_match(
    template: pd.DataFrame,
    candidates: pd.DataFrame,
    *,
    z_columns: list[str],
) -> pd.DataFrame:
    target_values = template[z_columns].to_numpy(dtype=float)
    candidate_values = candidates[z_columns].to_numpy(dtype=float)
    distances = np.sqrt(cdist(target_values, candidate_values, metric="sqeuclidean"))
    target_indices, candidate_indices = linear_sum_assignment(distances)
    assignments = {
        int(target_idx): (int(candidate_idx), float(distances[target_idx, candidate_idx]))
        for target_idx, candidate_idx in zip(target_indices, candidate_indices)
    }
    ordered_rows: list[pd.DataFrame] = []
    for target_idx in range(len(template)):
        candidate_idx, distance = assignments[target_idx]
        row = candidates.iloc[[candidate_idx]].copy()
        row["quality_template_id"] = target_idx + 1
        row["quality_match_distance"] = distance
        ordered_rows.append(row)
    return pd.concat(ordered_rows, ignore_index=True)


def match_panel(frame: pd.DataFrame, *, target_per_species: int, feature_columns: list[str] = MATCH_FEATURES) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    standardized, scaling = robust_standardize(frame, feature_columns)
    z_columns = [f"z_{column}" for column in feature_columns]
    counts = standardized.groupby("species").size()
    target = min(target_per_species, int(counts.min()))
    medoid_species, pairwise = quality_medoid_species(standardized, z_columns)
    template = build_common_template(standardized, medoid_species=medoid_species, z_columns=z_columns, target=target)
    matched_groups: list[pd.DataFrame] = []
    for species, candidates in standardized.groupby("species", sort=True):
        if species == medoid_species:
            matched = template.copy()
            matched["quality_template_id"] = range(1, len(template) + 1)
            matched["quality_match_distance"] = 0.0
        else:
            matched = optimal_quality_match(template, candidates, z_columns=z_columns)
        matched["quality_template_species"] = medoid_species
        matched_groups.append(matched)
    panel = pd.concat(matched_groups, ignore_index=True).sort_values(["species", "quality_template_id"]).reset_index(drop=True)
    panel["selection_rank"] = panel["quality_template_id"].astype(int)
    balance_rows: list[dict[str, Any]] = []
    template_values = template[z_columns]
    for species, group in panel.groupby("species", sort=True):
        row: dict[str, Any] = {
            "species": species,
            "n_matched": len(group),
            "n_images": int(group["filename"].nunique()),
            "n_specimens": int(group["specimen_group"].nunique()),
            "median_match_distance": float(group["quality_match_distance"].median()),
            "max_match_distance": float(group["quality_match_distance"].max()),
        }
        mean_diffs = []
        ks_values = []
        for feature, z_column in zip(feature_columns, z_columns):
            mean_diff = float(group[z_column].mean() - template_values[z_column].mean())
            ks = float(ks_2samp(group[z_column], template_values[z_column]).statistic)
            row[f"mean_difference_{feature}"] = mean_diff
            row[f"ks_{feature}"] = ks
            mean_diffs.append(abs(mean_diff))
            ks_values.append(ks)
        row["max_abs_standardized_mean_difference"] = max(mean_diffs)
        row["max_ks_distance"] = max(ks_values)
        row["quality_match_status"] = (
            "common_support"
            if row["max_abs_standardized_mean_difference"] <= 0.10 and row["max_ks_distance"] <= 0.25
            else "limited_overlap"
        )
        balance_rows.append(row)
    metadata = {"target_per_species": target, "template_species": medoid_species, "scaling": scaling}
    return panel, pd.DataFrame(balance_rows), {"metadata": metadata, "pairwise": pairwise}


def target_leakage_audit(frame: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for feature in MATCH_FEATURES:
        species_rhos = []
        for _species, group in frame.groupby("species", sort=True):
            rho, _p = spearmanr(group[feature], group["nuc_iod"], nan_policy="omit")
            if np.isfinite(rho):
                species_rhos.append(float(rho))
        rows.append(
            {
                "match_feature": feature,
                "global_spearman_vs_nuc_iod": float(frame[feature].corr(frame["nuc_iod"], method="spearman")),
                "median_within_species_spearman_vs_nuc_iod": float(np.median(species_rhos)),
                "max_abs_within_species_spearman_vs_nuc_iod": float(np.max(np.abs(species_rhos))),
            }
        )
    return pd.DataFrame(rows)


def load_review_registry(paths: list[Path]) -> pd.DataFrame:
    if not paths:
        return pd.DataFrame(columns=["species", "review_key", "decision"])
    parts: list[pd.DataFrame] = []
    for path in paths:
        if not path.exists():
            raise FileNotFoundError(path)
        frame = pd.read_csv(path, low_memory=False)
        required = {"species", "review_key", "decision"}
        missing = sorted(required - set(frame.columns))
        if missing:
            raise ValueError(f"Decision file {path} is missing columns: {missing}")
        frame = frame[["species", "review_key", "decision"]].copy()
        frame["decision"] = frame["decision"].fillna("").astype(str).str.strip().str.lower()
        invalid = sorted(set(frame["decision"]) - {"keep", "problem", "unsure"})
        if invalid:
            raise ValueError(f"Decision file {path} contains unsupported decisions: {invalid}")
        parts.append(frame)
    return pd.concat(parts, ignore_index=True).drop_duplicates(["species", "review_key"], keep="last")


def reviewed_panel_freeze_ready(
    frozen: pd.DataFrame,
    balance: pd.DataFrame,
    *,
    minimum_per_species: int = 30,
) -> bool:
    """Return whether reviewed keeps provide adequate common-species support."""
    expected_species = set(
        balance.loc[balance["quality_match_status"].eq("common_support"), "species"]
    )
    frozen_species = set(frozen["species"])
    if not expected_species or frozen_species != expected_species:
        return False
    counts = frozen.groupby("species", sort=False).size()
    return bool(counts.ge(minimum_per_species).all())


def build(args: argparse.Namespace) -> dict[str, Any]:
    args.output_dir.mkdir(parents=True, exist_ok=True)
    source_hash = sha256_file(args.candidates_csv)
    pixel_cache = args.output_dir / "strict_nuclei_with_pixel_quality.csv.gz"
    cache_manifest = args.output_dir / "pixel_quality_cache_manifest.json"
    cache_valid = False
    if pixel_cache.exists() and cache_manifest.exists() and not args.force_pixel_qc:
        cached = json.loads(cache_manifest.read_text())
        cache_valid = cached.get("candidates_sha256") == source_hash
    if cache_valid:
        scored = pd.read_csv(pixel_cache, low_memory=False)
    else:
        # These helpers pull in the raw-image stack (including scikit-image).
        # Import only when rebuilding pixel QC so pure matching tests and
        # hash-valid cached builds can load without raw-image dependencies.
        import audit_iod_measurement_quality as pixel_qc  # noqa: E402
        import build_balanced_genome_iod_sensitivity as source_qc  # noqa: E402

        strict = source_qc.prepare_candidate_frame(args.candidates_csv, max_overlap_label_count=1, min_overlap_fraction=0.98)
        # The legacy pixel-QC helper summarizes this historical column even
        # though it is not used by our matching features. The full strict pool
        # predates that derived rank, so use a neutral constant rather than
        # deriving it from intensity or any outcome variable.
        if "focus_rank" not in strict.columns:
            strict["focus_rank"] = 0.5
        scored = add_matching_features(pixel_qc.add_iod_quality_scores(strict))
        scored.to_csv(pixel_cache, index=False)
        cache_manifest.write_text(
            json.dumps({"candidates_csv": str(args.candidates_csv.resolve()), "candidates_sha256": source_hash, "n_rows": len(scored)}, indent=2) + "\n"
        )
    if not set(MATCH_FEATURES).issubset(scored.columns):
        scored = add_matching_features(scored)
    decisions = load_review_registry(args.decisions_csv)
    problem_keys = set(
        zip(
            decisions.loc[decisions["decision"].eq("problem"), "species"],
            decisions.loc[decisions["decision"].eq("problem"), "review_key"],
        )
    )
    if problem_keys:
        keep_mask = [(species, key) not in problem_keys for species, key in zip(scored["species"], scored["review_key"])]
        match_pool = scored.loc[keep_mask].copy()
    else:
        match_pool = scored
    panel, balance, details = match_panel(match_pool, target_per_species=args.target_per_species)
    leakage = target_leakage_audit(scored)
    panel = panel.merge(balance[["species", "quality_match_status"]], on="species", how="left", validate="many_to_one")
    panel["panel"] = "image_quality_matched_all_species"
    reviewed_keep_keys = set(
        zip(
            decisions.loc[decisions["decision"].eq("keep"), "species"],
            decisions.loc[decisions["decision"].eq("keep"), "review_key"],
        )
    )
    panel["review_status"] = [
        "reviewed_keep" if (species, key) in reviewed_keep_keys else "unreviewed_replacement"
        for species, key in zip(panel["species"], panel["review_key"])
    ]
    panel_path = args.output_dir / "image_quality_matched_nuclei_all_species_sensitivity.csv.gz"
    primary_path = args.output_dir / "image_quality_matched_nuclei_common_support.csv.gz"
    replacement_path = args.output_dir / "image_quality_matched_replacement_review.csv.gz"
    frozen_path = args.output_dir / "image_quality_matched_nuclei_frozen_reviewed.csv.gz"
    frozen_support_path = args.output_dir / "image_quality_matched_nuclei_frozen_support.csv"
    balance_path = args.output_dir / "image_quality_match_balance.csv"
    leakage_path = args.output_dir / "image_quality_target_leakage_audit.csv"
    pairwise_path = args.output_dir / "image_quality_species_pairwise_distance.csv"
    manifest_path = args.output_dir / "manifest.json"
    panel.to_csv(panel_path, index=False)
    primary = panel.loc[panel["quality_match_status"].eq("common_support")].copy()
    primary["panel"] = "image_quality_matched_common_support"
    primary.to_csv(primary_path, index=False)
    replacements = primary.loc[primary["review_status"].eq("unreviewed_replacement")].copy()
    replacements = replacements.sort_values(["species", "quality_match_distance", "review_key"], kind="mergesort")
    replacements["selection_rank"] = replacements.groupby("species").cumcount() + 1
    replacements["panel"] = "image_quality_matched_replacement_review"
    replacements.to_csv(replacement_path, index=False)
    frozen = primary.loc[primary["review_status"].eq("reviewed_keep")].copy()
    frozen_support = (
        frozen.groupby("species", sort=True)
        .agg(
            n_frozen_nuclei=("review_key", "size"),
            n_images=("filename", "nunique"),
            n_specimens=("specimen_group", "nunique"),
        )
        .reset_index()
    )
    freeze_ready = reviewed_panel_freeze_ready(frozen, balance, minimum_per_species=30)
    if freeze_ready:
        frozen["panel"] = "image_quality_matched_frozen_reviewed"
        frozen.to_csv(frozen_path, index=False)
        frozen_support.to_csv(frozen_support_path, index=False)
    balance.to_csv(balance_path, index=False, float_format="%.8f")
    leakage.to_csv(leakage_path, index=False, float_format="%.8f")
    details["pairwise"].to_csv(pairwise_path, index=False, float_format="%.8f")
    manifest = {
        "purpose": "Nuclei matched across species on technical image quality before any IOD comparison.",
        "source_candidates_csv": str(args.candidates_csv.resolve()),
        "source_candidates_sha256": source_hash,
        "strict_pixel_quality_rows": int(len(scored)),
        "match_features": MATCH_FEATURES,
        "forbidden_match_features": FORBIDDEN_MATCH_FEATURES,
        "target_requested_per_species": args.target_per_species,
        "target_matched_per_species": details["metadata"]["target_per_species"],
        "quality_template_species": details["metadata"]["template_species"],
        "matched_panel_csv": str(panel_path.resolve()),
        "matched_panel_sha256": sha256_file(panel_path),
        "common_support_panel_csv": str(primary_path.resolve()),
        "common_support_panel_sha256": sha256_file(primary_path),
        "common_support_species": sorted(primary["species"].unique().tolist()),
        "limited_overlap_species": sorted(balance.loc[balance["quality_match_status"].eq("limited_overlap"), "species"].tolist()),
        "decision_files": [str(path.resolve()) for path in args.decisions_csv],
        "decision_file_sha256": {str(path.resolve()): sha256_file(path) for path in args.decisions_csv},
        "n_problem_decisions": int(len(problem_keys)),
        "n_reviewed_keeps_in_common_support_panel": int(primary["review_status"].eq("reviewed_keep").sum()),
        "n_unreviewed_replacements": int(len(replacements)),
        "replacement_review_csv": str(replacement_path.resolve()),
        "freeze_status": "frozen_reviewed_only_no_further_replacement" if freeze_ready else "not_ready",
        "frozen_reviewed_csv": str(frozen_path.resolve()) if freeze_ready else None,
        "frozen_reviewed_sha256": sha256_file(frozen_path) if freeze_ready else None,
        "frozen_support_csv": str(frozen_support_path.resolve()) if freeze_ready else None,
        "n_frozen_reviewed": int(len(frozen)) if freeze_ready else 0,
        "minimum_frozen_nuclei_per_species": int(frozen_support["n_frozen_nuclei"].min()) if freeze_ready else None,
        "maximum_frozen_nuclei_per_species": int(frozen_support["n_frozen_nuclei"].max()) if freeze_ready else None,
        "balance_csv": str(balance_path.resolve()),
        "target_leakage_audit_csv": str(leakage_path.resolve()),
        "pairwise_quality_distance_csv": str(pairwise_path.resolve()),
        "feature_scaling": details["metadata"]["scaling"],
    }
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    return manifest


def main() -> None:
    args = parse_args()
    print(json.dumps(build(args), indent=2))


if __name__ == "__main__":
    main()
