#!/usr/bin/env python3
"""Extend the exact finalized morphology and IOD review pages.

The historical frozen tables are never modified. Existing reviewed rows enter
the extended pages as per-record default keeps; new rows remain unreviewed.
Morphology candidates use the literal-largest hard gates. IOD candidates are
matched to the original frozen technical-quality template without using IOD,
nuclear darkness, nucleus area, or cell area as matching features.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment
from scipy.spatial.distance import cdist
from scipy.stats import ks_2samp


PROJECT_ROOT = Path(__file__).resolve().parents[2]
CELL_REPO = PROJECT_ROOT.parent / "cellprofiler_test"
sys.path.insert(0, str(PROJECT_ROOT / "path_analysis/scripts"))
sys.path.insert(0, str(CELL_REPO / "scripts"))

import build_balanced_genome_iod_sensitivity as source_qc  # noqa: E402
import build_balanced_genome_iod_viewer as viewer  # noqa: E402
import build_image_quality_matched_genome_iod_panel as quality_match  # noqa: E402
import build_largest_cell_mask_review as largest_review  # noqa: E402


DERIVED = PROJECT_ROOT / "path_analysis/data/external/derived"
MORPH_ROOT = DERIVED / "largest_cell_mask_review"
IOD_ROOT = DERIVED / "image_quality_matched_genome_iod"
DEFAULT_NEW_CANDIDATES = (
    CELL_REPO
    / "output/runs/life_history_target_vsi_v1/top50_yolo_exact_v1/scored/all_scored_masks.csv"
)
DEFAULT_MORPH_FROZEN = MORPH_ROOT / "frozen_largest_cell_mask_top50.csv.gz"
DEFAULT_MORPH_EXTENSION = MORPH_ROOT / "largest_cell_mask_review_extension_20260714.csv.gz"
DEFAULT_MORPH_VIEWER = MORPH_ROOT / "viewer_frozen_analysis_top50"
DEFAULT_IOD_FROZEN = IOD_ROOT / "image_quality_matched_nuclei_frozen_reviewed.csv.gz"
DEFAULT_IOD_MATCHED = IOD_ROOT / "image_quality_matched_nuclei_all_species_sensitivity.csv.gz"
DEFAULT_IOD_MANIFEST = IOD_ROOT / "manifest.json"
DEFAULT_IOD_EXTENSION = IOD_ROOT / "image_quality_matched_review_extension_20260714.csv.gz"
DEFAULT_IOD_NEW_SCORED = IOD_ROOT / "new_species_strict_nuclei_with_pixel_quality_20260714.csv.gz"
DEFAULT_IOD_VIEWER = IOD_ROOT / "viewer"
MORPH_PANEL = "literal_largest_cell_mask_review_extension_20260714"
IOD_PANEL = "image_quality_matched_review_extension_20260714"
NEW_SPECIES = ("D. aeneus", "D. orestes", "D. wrighti")
LEGACY_IOD_REVIEW_SPECIES = ("D. ochrophaeus",)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_json(path: Path, payload: dict[str, Any]) -> None:
    path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")


def build_new_morphology_queue(
    candidates: pd.DataFrame,
    *,
    target_per_species: int,
    review_depth: int,
) -> tuple[pd.DataFrame, dict[str, dict[str, int | bool]]]:
    queues: list[pd.DataFrame] = []
    support: dict[str, dict[str, int | bool]] = {}
    for species in NEW_SPECIES:
        group = candidates.loc[candidates["species"].eq(species)].copy()
        ranked = group.sort_values(
            ["cell_area_um2", "review_key"],
            ascending=[False, True],
            kind="mergesort",
        )
        queued = ranked.head(review_depth).copy()
        queued["selection_rank"] = np.arange(1, len(queued) + 1)
        queued["selection_metric"] = "cell_area_um2_descending_after_fixed_hard_eligibility"
        queues.append(queued)
        support[species] = {
            "hard_eligible": int(len(ranked)),
            "rendered_for_review": int(len(queued)),
            "target_per_species": int(target_per_species),
            "shortfall_to_target": int(max(0, target_per_species - len(ranked))),
            "target_available": bool(len(ranked) >= target_per_species),
        }
    return pd.concat(queues, ignore_index=True), support


def apply_frozen_scaling(
    frame: pd.DataFrame,
    scaling: dict[str, dict[str, float]],
) -> pd.DataFrame:
    out = frame.copy()
    for feature in quality_match.MATCH_FEATURES:
        params = scaling[feature]
        values = pd.to_numeric(out[feature], errors="coerce")
        clipped = values.clip(lower=params["clip_low"], upper=params["clip_high"])
        out[f"z_{feature}"] = (clipped - params["median"]) / params["scale_iqr"]
    return out


def localize_tile_coordinates_for_pixel_qc(frame: pd.DataFrame) -> pd.DataFrame:
    """Translate slide-global centroids onto their tile-backed image arrays."""
    out = frame.copy()
    path_text = out["nucleus_source_image_path"].fillna("").astype(str)
    is_tile = path_text.str.contains(r"tile_y\d+_x\d+", regex=True)
    for axis in ("x", "y"):
        centroid_column = f"nuc_centroid_{axis}"
        origin_column = f"tile_{axis}0"
        centroid = pd.to_numeric(out[centroid_column], errors="coerce")
        origin = pd.to_numeric(out[origin_column], errors="coerce").fillna(0.0)
        out[f"source_global_{centroid_column}"] = centroid
        out[f"source_global_{origin_column}"] = origin
        recorded_globally = is_tile & centroid.ge(origin)
        out.loc[recorded_globally, centroid_column] = (
            centroid.loc[recorded_globally] - origin.loc[recorded_globally]
        )
        out.loc[is_tile, origin_column] = 0.0
    return out


def match_species_to_template(
    template: pd.DataFrame,
    candidates: pd.DataFrame,
) -> tuple[pd.DataFrame, dict[str, float | int | str]]:
    z_columns = [f"z_{feature}" for feature in quality_match.MATCH_FEATURES]
    if template.empty or candidates.empty:
        return candidates.iloc[0:0].copy(), {
            "n_matched": 0,
            "max_abs_standardized_mean_difference": float("nan"),
            "max_ks_distance": float("nan"),
            "quality_match_status": "insufficient_candidates",
        }
    template_values = template[z_columns].to_numpy(dtype=float)
    candidate_values = candidates[z_columns].to_numpy(dtype=float)
    distances = np.sqrt(cdist(template_values, candidate_values, metric="sqeuclidean"))
    template_indices, candidate_indices = linear_sum_assignment(distances)
    rows: list[pd.DataFrame] = []
    for template_idx, candidate_idx in sorted(zip(template_indices, candidate_indices)):
        row = candidates.iloc[[int(candidate_idx)]].copy()
        row["quality_template_id"] = int(template.iloc[int(template_idx)]["quality_template_id"])
        row["quality_match_distance"] = float(distances[int(template_idx), int(candidate_idx)])
        rows.append(row)
    matched = pd.concat(rows, ignore_index=True).sort_values("quality_template_id").reset_index(drop=True)
    matched_template = template.loc[
        template["quality_template_id"].isin(matched["quality_template_id"])
    ].sort_values("quality_template_id")
    mean_diffs = [
        abs(float(matched[column].mean() - matched_template[column].mean()))
        for column in z_columns
    ]
    ks_values = [
        float(ks_2samp(matched[column], matched_template[column]).statistic)
        for column in z_columns
    ]
    max_mean = max(mean_diffs)
    max_ks = max(ks_values)
    status = "common_support" if max_mean <= 0.10 and max_ks <= 0.25 else "limited_overlap"
    matched["quality_match_status"] = status
    matched["selection_rank"] = np.arange(1, len(matched) + 1)
    return matched, {
        "n_matched": int(len(matched)),
        "median_match_distance": float(matched["quality_match_distance"].median()),
        "max_match_distance": float(matched["quality_match_distance"].max()),
        "max_abs_standardized_mean_difference": max_mean,
        "max_ks_distance": max_ks,
        "quality_match_status": status,
    }


def compose_iod_review_extension(
    frozen: pd.DataFrame,
    original_matched: pd.DataFrame,
    new_matched: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Combine frozen, recovered legacy, and newly matched IOD review cohorts."""
    frozen = frozen.copy()
    recovered = original_matched.loc[
        original_matched["species"].isin(LEGACY_IOD_REVIEW_SPECIES)
    ].copy()
    new_matched = new_matched.copy()

    frozen["panel"] = IOD_PANEL
    frozen["viewer_default_decision"] = "keep"
    frozen["review_cohort"] = "previously_frozen"
    recovered["panel"] = IOD_PANEL
    recovered["viewer_default_decision"] = ""
    recovered["review_cohort"] = "legacy_limited_overlap_unreviewed"
    recovered["review_status"] = "unreviewed"
    new_matched["panel"] = IOD_PANEL
    new_matched["viewer_default_decision"] = ""
    new_matched["review_cohort"] = "new_species_unreviewed"
    new_matched["review_status"] = "unreviewed"

    combined = pd.concat([frozen, recovered, new_matched], ignore_index=True, sort=False)
    duplicate_keys = combined.loc[combined["review_key"].duplicated(keep=False), "review_key"]
    if not duplicate_keys.empty:
        examples = sorted(duplicate_keys.astype(str).unique())[:5]
        raise ValueError(f"Duplicate IOD review keys across cohorts: {examples}")
    return combined, recovered


def extract_existing_new_iod_matches(existing_extension: pd.DataFrame) -> pd.DataFrame:
    """Recover already computed July matches without rerunning pixel-level QC."""
    matched = existing_extension.loc[existing_extension["species"].isin(NEW_SPECIES)].copy()
    missing = sorted(set(NEW_SPECIES) - set(matched["species"].unique()))
    if missing:
        raise ValueError(f"Existing IOD extension is missing new species: {missing}")
    duplicate_keys = matched.loc[matched["review_key"].duplicated(keep=False), "review_key"]
    if not duplicate_keys.empty:
        examples = sorted(duplicate_keys.astype(str).unique())[:5]
        raise ValueError(f"Duplicate review keys in existing new-species IOD matches: {examples}")
    return matched


def render_rows(
    rows: pd.DataFrame,
    *,
    output_dir: Path,
    start_index: int,
    force: bool,
    fixed_indices: dict[str, int] | None = None,
) -> list[dict[str, Any]]:
    records: list[dict[str, Any]] = []
    for offset, (_idx, row) in enumerate(rows.iterrows()):
        index = fixed_indices.get(str(row["review_key"])) if fixed_indices else None
        if index is None:
            index = start_index + offset
        records.append(
            viewer.process_row(
                row,
                int(index),
                output_dir,
                thumb_size=384,
                crop_margin=64,
                force=force,
            )
        )
    return records


def build_morphology_extension(args: argparse.Namespace) -> dict[str, Any]:
    frozen = pd.read_csv(args.morph_frozen, low_memory=False).sort_values(
        ["species", "selection_rank", "review_key"], kind="mergesort"
    ).reset_index(drop=True)
    eligible = largest_review.prepare_candidates(args.new_candidates)
    new_queue, support = build_new_morphology_queue(
        eligible,
        target_per_species=args.morph_target,
        review_depth=args.morph_review_depth,
    )
    frozen["panel"] = MORPH_PANEL
    frozen["viewer_default_decision"] = "keep"
    frozen["review_cohort"] = "previously_frozen"
    new_queue["panel"] = MORPH_PANEL
    new_queue["viewer_default_decision"] = ""
    new_queue["review_cohort"] = "new_species_unreviewed"
    combined = pd.concat([frozen, new_queue], ignore_index=True, sort=False)
    combined.to_csv(args.morph_extension, index=False)

    args.morph_viewer.mkdir(parents=True, exist_ok=True)
    old_records = render_rows(
        frozen,
        output_dir=args.morph_viewer,
        start_index=1,
        force=False,
    )
    new_records = render_rows(
        new_queue,
        output_dir=args.morph_viewer,
        start_index=len(frozen) + 1,
        force=args.force_new_assets,
    )
    records = old_records + new_records
    viewer.write_index(
        args.morph_viewer,
        records,
        panel=MORPH_PANEL,
        title="Frozen Cell and Nucleus Masks + New Species Review",
        storage_key="largest-cell-mask-review-extension-20260714-v1",
        default_decision=None,
        blind_target_metrics=False,
        common_support_only=False,
        auto_fill_target_per_species=args.morph_target,
    )
    summary = {
        "panel": MORPH_PANEL,
        "purpose": "Preserve the finalized 21-species morphology masks and review literal-largest candidates for new species.",
        "frozen_source_csv": str(args.morph_frozen.resolve()),
        "frozen_source_sha256": sha256_file(args.morph_frozen),
        "new_candidates_csv": str(args.new_candidates.resolve()),
        "new_candidates_sha256": sha256_file(args.new_candidates),
        "extension_csv": str(args.morph_extension.resolve()),
        "extension_sha256": sha256_file(args.morph_extension),
        "index_html": str((args.morph_viewer / "index.html").resolve()),
        "existing_frozen_records": int(len(frozen)),
        "new_review_records": int(len(new_queue)),
        "auto_fill_target_per_species": int(args.morph_target),
        "species_counts": combined.groupby("species").size().to_dict(),
        "new_species_support": support,
        "freeze_status": "new_species_unreviewed",
    }
    write_json(args.morph_viewer / "summary.json", summary)
    return summary


def build_iod_extension(args: argparse.Namespace) -> dict[str, Any]:
    original_manifest = json.loads(args.iod_manifest.read_text())
    frozen = pd.read_csv(args.iod_frozen, low_memory=False)
    original_matched = pd.read_csv(args.iod_matched, low_memory=False)
    template_species = str(original_manifest["quality_template_species"])
    if args.reuse_existing_iod_extension:
        existing_extension = pd.read_csv(args.iod_extension, low_memory=False)
        new_matched = extract_existing_new_iod_matches(existing_extension)
        previous_summary_path = args.iod_viewer / "summary.json"
        previous_summary = (
            json.loads(previous_summary_path.read_text()) if previous_summary_path.exists() else {}
        )
        support = previous_summary.get("new_species_support", {})
    else:
        # Importing this helper loads the raw-image stack; avoid it for viewer-only rebuilds.
        import audit_iod_measurement_quality as pixel_qc

        strict = source_qc.prepare_candidate_frame(
            args.new_candidates,
            max_overlap_label_count=1,
            min_overlap_fraction=0.98,
        )
        strict = strict.loc[strict["species"].isin(NEW_SPECIES)].copy()
        if "focus_rank" not in strict:
            strict["focus_rank"] = 0.5
        pixel_qc_input = localize_tile_coordinates_for_pixel_qc(strict)
        scored = quality_match.add_matching_features(pixel_qc.add_iod_quality_scores(pixel_qc_input))
        scored = apply_frozen_scaling(scored, original_manifest["feature_scaling"])
        scored.to_csv(args.iod_new_scored, index=False)

        template = original_matched.loc[original_matched["species"].eq(template_species)].copy()
        template = template.sort_values("quality_template_id", kind="mergesort")
        matched_groups: list[pd.DataFrame] = []
        support = {}
        for species in NEW_SPECIES:
            group = scored.loc[scored["species"].eq(species)].copy()
            matched, metrics = match_species_to_template(template, group)
            matched["quality_template_species"] = template_species
            matched_groups.append(matched)
            support[species] = {
                "strict_candidates": int(len(group)),
                "target_per_species": int(args.iod_target),
                "shortfall_to_target": int(max(0, args.iod_target - len(matched))),
                **metrics,
            }
        new_matched = pd.concat(matched_groups, ignore_index=True, sort=False)

    combined, recovered = compose_iod_review_extension(frozen, original_matched, new_matched)
    frozen = combined.loc[combined["review_cohort"].eq("previously_frozen")].copy()
    new_matched = combined.loc[combined["review_cohort"].eq("new_species_unreviewed")].copy()
    combined.to_csv(args.iod_extension, index=False)

    args.iod_viewer.mkdir(parents=True, exist_ok=True)
    original_order = original_matched.sort_values(
        ["species", "selection_rank", "review_key"], kind="mergesort"
    ).reset_index(drop=True)
    original_indices = {
        str(review_key): int(index + 1)
        for index, review_key in enumerate(original_order["review_key"])
    }
    frozen = frozen.sort_values(["species", "selection_rank", "review_key"], kind="mergesort")
    old_records = render_rows(
        frozen,
        output_dir=args.iod_viewer,
        start_index=1,
        force=False,
        fixed_indices=original_indices,
    )
    recovered_records = render_rows(
        recovered.sort_values(["species", "selection_rank", "review_key"], kind="mergesort"),
        output_dir=args.iod_viewer,
        start_index=1,
        force=False,
        fixed_indices=original_indices,
    )
    new_records = render_rows(
        new_matched.sort_values(["species", "selection_rank", "review_key"], kind="mergesort"),
        output_dir=args.iod_viewer,
        start_index=len(original_matched) + 1,
        force=args.force_new_assets,
    )
    records = old_records + recovered_records + new_records
    viewer.write_index(
        args.iod_viewer,
        records,
        panel=IOD_PANEL,
        title="Finalized + Extended Quality-Matched IOD Mask Review",
        storage_key="image-quality-matched-nucleus-extension-20260714-v1",
        default_decision=None,
        blind_target_metrics=True,
        common_support_only=False,
        sort_by_match_distance=False,
        auto_fill_target_per_species=args.iod_target,
    )
    summary = {
        "panel": IOD_PANEL,
        "purpose": "Preserve finalized reviewed IOD nuclei, recover legacy limited-overlap species for explicit review, and review newly matched species.",
        "matching_features": quality_match.MATCH_FEATURES,
        "forbidden_matching_features": quality_match.FORBIDDEN_MATCH_FEATURES,
        "quality_template_species": template_species,
        "quality_template_source": str(args.iod_matched.resolve()),
        "frozen_source_csv": str(args.iod_frozen.resolve()),
        "frozen_source_sha256": sha256_file(args.iod_frozen),
        "new_pixel_quality_csv": str(args.iod_new_scored.resolve()),
        "new_pixel_quality_sha256": sha256_file(args.iod_new_scored),
        "extension_csv": str(args.iod_extension.resolve()),
        "extension_sha256": sha256_file(args.iod_extension),
        "index_html": str((args.iod_viewer / "index.html").resolve()),
        "existing_frozen_records": int(len(frozen)),
        "recovered_legacy_review_records": int(len(recovered)),
        "recovered_legacy_review_species": sorted(recovered["species"].unique().tolist()),
        "new_review_records": int(len(new_matched)),
        "auto_fill_target_per_species": int(args.iod_target),
        "species_counts": combined.groupby("species").size().to_dict(),
        "new_species_support": support,
        "freeze_status": "legacy_limited_overlap_and_new_species_unreviewed",
    }
    write_json(args.iod_viewer / "summary.json", summary)
    return summary


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--new-candidates", type=Path, default=DEFAULT_NEW_CANDIDATES)
    parser.add_argument("--morph-frozen", type=Path, default=DEFAULT_MORPH_FROZEN)
    parser.add_argument("--morph-extension", type=Path, default=DEFAULT_MORPH_EXTENSION)
    parser.add_argument("--morph-viewer", type=Path, default=DEFAULT_MORPH_VIEWER)
    parser.add_argument("--morph-target", type=int, default=50)
    parser.add_argument("--morph-review-depth", type=int, default=100)
    parser.add_argument("--iod-frozen", type=Path, default=DEFAULT_IOD_FROZEN)
    parser.add_argument("--iod-matched", type=Path, default=DEFAULT_IOD_MATCHED)
    parser.add_argument("--iod-manifest", type=Path, default=DEFAULT_IOD_MANIFEST)
    parser.add_argument("--iod-extension", type=Path, default=DEFAULT_IOD_EXTENSION)
    parser.add_argument("--iod-new-scored", type=Path, default=DEFAULT_IOD_NEW_SCORED)
    parser.add_argument("--iod-viewer", type=Path, default=DEFAULT_IOD_VIEWER)
    parser.add_argument("--iod-target", type=int, default=40)
    parser.add_argument("--force-new-assets", action="store_true")
    parser.add_argument(
        "--reuse-existing-iod-extension",
        action="store_true",
        help="Reuse already computed July matches and only rebuild the combined viewer.",
    )
    parser.add_argument("--skip-morphology", action="store_true")
    parser.add_argument("--skip-iod", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.skip_morphology and args.skip_iod:
        raise ValueError("Cannot skip both morphology and IOD builds")
    required = [args.new_candidates]
    if not args.skip_morphology:
        required.append(args.morph_frozen)
    if not args.skip_iod:
        required.extend([args.iod_frozen, args.iod_matched, args.iod_manifest])
    for path in required:
        if not path.exists():
            raise FileNotFoundError(path)
    morphology = None if args.skip_morphology else build_morphology_extension(args)
    iod = None if args.skip_iod else build_iod_extension(args)
    print(json.dumps({"morphology": morphology, "iod": iod}, indent=2))


if __name__ == "__main__":
    main()
