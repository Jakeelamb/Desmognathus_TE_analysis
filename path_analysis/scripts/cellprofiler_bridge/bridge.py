from __future__ import annotations

from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

from .common import (
    GENOME_GB_PER_PG,
    REFERENCE_GENOME_PG,
    canonical_species,
    existing_pct,
    file_record,
    nonempty_pct,
    sha256_for_file,
    write_json,
)
from .pipeline_runs import get_run_paths, metadata_path


PROJECT_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_ACTIVE_RAW_RUN_TAG = "full_dataset_v1"
DEFAULT_ACTIVE_MIXED_RUN_TAG = "mixed_cellpose_yolo_full_dataset_v1"
GENOME_CALIBRATION_IMAGE_TYPE = "brightfield"


def load_metadata(cellprofiler_root: Path) -> pd.DataFrame:
    path = metadata_path(cellprofiler_root)
    df = pd.read_csv(path)
    df["filename"] = df["filename"].astype(str)
    df["species"] = df["species"].map(canonical_species)
    return df


def _required_paths_missing(paths: list[Path]) -> list[str]:
    return [str(path.resolve()) for path in paths if not path.exists()]


def build_morphology_bundle(
    cellprofiler_root: Path,
    out_dir: Path,
    summary_rows: list[dict[str, object]],
    gap_rows: list[dict[str, object]],
    mixed_run_tag: str,
) -> dict[str, object]:
    linkage_root = get_run_paths(cellprofiler_root, mixed_run_tag).root / "linkage"
    species_summary_path = linkage_root / "species_summary.csv"
    size_summary_path = linkage_root / "species_size_correlation_summary.csv"
    nucleus_iod_summary_path = linkage_root / "species_nucleus_iod_summary.csv"
    linked_pairs_path = linkage_root / "linked_nucleus_pairs.csv.gz"
    summary_json_path = linkage_root / "summary.json"

    required = [species_summary_path, linked_pairs_path]
    missing = _required_paths_missing(required)
    if missing:
        gap_rows.append(
            {
                "bundle": "morphology",
                "gap_type": "missing_morphology_bundle",
                "species": pd.NA,
                "details": f"One or more mixed-linkage morphology files are missing: {missing}",
            }
        )
        return {"present": False, "summary_path": pd.NA, "image_trace_path": pd.NA}

    linked_pairs = pd.read_csv(linked_pairs_path)
    species_summary = pd.read_csv(species_summary_path)

    for frame in [species_summary, linked_pairs]:
        frame["species"] = frame["species"].map(canonical_species)

    strict_pairs = linked_pairs[linked_pairs["keep_strict_core"].fillna(False)].copy()
    strict_pairs["cytoplasm_area_um2"] = pd.to_numeric(strict_pairs["cytoplasm_area_um2"], errors="coerce")

    if size_summary_path.exists():
        size_summary = pd.read_csv(size_summary_path)
        size_summary["species"] = size_summary["species"].map(canonical_species)
    else:
        size_summary = (
            strict_pairs.groupby("species", dropna=True)
            .agg(
                n_strict_core_pairs=("species", "size"),
                cell_area_median_um2=("cell_area_um2", "median"),
                nucleus_area_median_um2=("nuc_area_um2", "median"),
                median_nc_area_ratio=("nc_area_ratio", "median"),
            )
            .reset_index()
        )

    cytoplasm_summary = (
        strict_pairs.groupby("species", dropna=True)
        .agg(species_median_cytoplasm_area_um2=("cytoplasm_area_um2", "median"))
        .reset_index()
    )

    morphology = (
        species_summary[["species", "n_specimens", "n_images", "n_analysis_ready_images"]]
        .rename(
            columns={
                "n_specimens": "n_specimens_strict",
                "n_images": "n_images_strict",
            }
        )
        .merge(
            size_summary[
                [
                    "species",
                    "n_strict_core_pairs",
                    "cell_area_median_um2",
                    "nucleus_area_median_um2",
                    "median_nc_area_ratio",
                ]
            ].rename(
                columns={
                    "n_strict_core_pairs": "n_pairs_strict",
                    "cell_area_median_um2": "species_median_cell_area_um2",
                    "nucleus_area_median_um2": "species_median_nuc_area_um2",
                    "median_nc_area_ratio": "species_median_nc_ratio",
                }
            ),
            on="species",
            how="inner",
            validate="one_to_one",
        )
        .merge(cytoplasm_summary, on="species", how="left", validate="one_to_one")
        .sort_values("species")
        .reset_index(drop=True)
    )
    morphology["low_support_for_species_median"] = (
        morphology["n_images_strict"].lt(2) | morphology["n_pairs_strict"].lt(250)
    )
    morphology["source_file"] = str(linked_pairs_path.resolve())
    morphology["source_sha256"] = sha256_for_file(linked_pairs_path)

    summary_out = out_dir / "cellprofiler_species_morphology_summary.csv"
    morphology.to_csv(summary_out, index=False)

    image_trace = (
        strict_pairs.groupby(["species", "filename"], dropna=True)
        .agg(
            specimen_id=("specimen_id", "first"),
            n_pairs_strict=("filename", "size"),
            median_cell_area_um2_strict=("cell_area_um2", "median"),
            median_nuc_area_um2_strict=("nuc_area_um2", "median"),
            median_nuc_iod_strict=("nuc_iod", "median"),
            median_nc_ratio_strict=("nc_area_ratio", "median"),
            median_cytoplasm_area_um2_strict=("cytoplasm_area_um2", "median"),
            cell_mask_path=("cell_mask_path", "first"),
            nucleus_mask_path=("nucleus_mask_path", "first"),
            roi_zip_path=("roi_zip_path", "first"),
            cell_tile_manifest_path=("cell_tile_manifest_path", "first"),
            nucleus_tile_manifest_path=("nucleus_tile_manifest_path", "first"),
            cell_source_image_path=("cell_source_image_path", "first"),
            nucleus_source_image_path=("nucleus_source_image_path", "first"),
            cell_run_manifest_path=("cell_run_manifest_path", "first"),
            nucleus_run_manifest_path=("nucleus_run_manifest_path", "first"),
        )
        .reset_index()
        .sort_values(["species", "filename"])
        .reset_index(drop=True)
    )
    image_trace["linked_pairs_source_path"] = str(linked_pairs_path.resolve())
    image_trace["linked_pairs_source_sha256"] = sha256_for_file(linked_pairs_path)

    image_trace_out = out_dir / "cellprofiler_species_morphology_image_trace.csv"
    image_trace.to_csv(image_trace_out, index=False)

    summary_rows.append(
        {
            "bundle": "morphology",
            "present": True,
            "n_species": int(morphology["species"].nunique()),
            "n_images": int(image_trace["filename"].nunique()),
            "n_rows": int(len(strict_pairs)),
            "nonempty_source_image_path_pct": round(nonempty_pct(strict_pairs["cell_source_image_path"]), 2),
            "existing_source_image_path_pct": round(existing_pct(strict_pairs["cell_source_image_path"]), 2),
            "nonempty_tile_manifest_path_pct": round(nonempty_pct(strict_pairs["cell_tile_manifest_path"]), 2),
            "existing_tile_manifest_path_pct": round(existing_pct(strict_pairs["cell_tile_manifest_path"]), 2),
            "nonempty_mask_path_pct": round(nonempty_pct(strict_pairs["nucleus_mask_path"]), 2),
            "existing_mask_path_pct": round(existing_pct(strict_pairs["nucleus_mask_path"]), 2),
            "nonempty_roi_zip_path_pct": round(nonempty_pct(strict_pairs["roi_zip_path"]), 2),
            "existing_roi_zip_path_pct": round(existing_pct(strict_pairs["roi_zip_path"]), 2),
            "notes": "Strict-core linked morphology pairs are traceable to image, tile manifest, and saved cell/nucleus mask paths.",
        }
    )

    if nonempty_pct(strict_pairs["roi_zip_path"]) == 0.0:
        gap_rows.append(
            {
                "bundle": "morphology",
                "gap_type": "missing_roi_zip_trace",
                "species": pd.NA,
                "details": "Current mixed-linkage morphology rows have saved nucleus mask paths but no ROI ZIP paths.",
            }
        )

    return {
        "present": True,
        "origin": "upstream_mixed_linkage_run",
        "summary_path": str(summary_out.resolve()),
        "summary_sha256": sha256_for_file(summary_out),
        "image_trace_path": str(image_trace_out.resolve()),
        "image_trace_sha256": sha256_for_file(image_trace_out),
        "source_files": [
            file_record("species_summary", species_summary_path),
            file_record("size_summary", size_summary_path),
            file_record("species_nucleus_iod_summary", nucleus_iod_summary_path),
            file_record("linked_nucleus_pairs", linked_pairs_path),
            file_record("summary_json", summary_json_path),
        ],
    }


def build_linked_genome_trace(
    cellprofiler_root: Path,
    out_dir: Path,
    summary_rows: list[dict[str, object]],
    gap_rows: list[dict[str, object]],
    mixed_run_tag: str,
) -> dict[str, object]:
    linkage_root = get_run_paths(cellprofiler_root, mixed_run_tag).root / "linkage"
    linked_pairs_path = linkage_root / "linked_nucleus_pairs.csv.gz"
    image_summary_path = linkage_root / "image_summary.csv"
    species_summary_path = linkage_root / "species_summary.csv"
    summary_json_path = linkage_root / "summary.json"

    required = [linked_pairs_path, image_summary_path, species_summary_path]
    missing = _required_paths_missing(required)
    if missing:
        gap_rows.append(
            {
                "bundle": "genome_linked_yolo",
                "gap_type": "missing_linked_genome_bundle",
                "species": pd.NA,
                "details": f"One or more mixed-linkage genome source files are missing: {missing}",
            }
        )
        return {"present": False, "image_trace_path": pd.NA, "source_files": []}

    linked_pairs = pd.read_csv(linked_pairs_path, keep_default_na=False)
    image_summary = pd.read_csv(image_summary_path, keep_default_na=False)
    for frame in [linked_pairs, image_summary]:
        frame["species"] = frame["species"].map(canonical_species)

    strict_mask = (
        linked_pairs.get("keep_strict_core", pd.Series(False, index=linked_pairs.index))
        .fillna(False)
        .astype(str)
        .str.lower()
        .isin(["true", "1", "yes"])
    )
    strict_pairs = linked_pairs.loc[strict_mask].copy()
    if strict_pairs.empty:
        gap_rows.append(
            {
                "bundle": "genome_linked_yolo",
                "gap_type": "empty_linked_genome_trace",
                "species": pd.NA,
                "details": "The mixed-linkage run contains no strict-core linked nucleus rows, so no linked genome trace can be built.",
            }
        )
        return {"present": False, "image_trace_path": pd.NA, "source_files": []}

    strict_pairs["nuc_area_um2"] = pd.to_numeric(strict_pairs["nuc_area_um2"], errors="coerce")
    strict_pairs["nuc_iod"] = pd.to_numeric(strict_pairs["nuc_iod"], errors="coerce")
    strict_pairs["analysis_ready_image"] = False

    image_flags = image_summary[["filename", "analysis_ready_image"]].copy()
    image_flags["analysis_ready_image"] = (
        image_flags["analysis_ready_image"].fillna(False).astype(str).str.lower().isin(["true", "1", "yes"])
    )

    image_trace = (
        strict_pairs.groupby(["species", "filename", "specimen_id", "image_type"], dropna=False)
        .agg(
            n_nuclei=("nuc_iod", "size"),
            median_nucleus_area_um2=("nuc_area_um2", "median"),
            median_nucleus_iod=("nuc_iod", "median"),
            cell_source_image_path=("cell_source_image_path", "first"),
            nucleus_source_image_path=("nucleus_source_image_path", "first"),
            cell_tile_manifest_path=("cell_tile_manifest_path", "first"),
            nucleus_tile_manifest_path=("nucleus_tile_manifest_path", "first"),
            cell_mask_path=("cell_mask_path", "first"),
            nucleus_mask_path=("nucleus_mask_path", "first"),
            roi_zip_path=("roi_zip_path", "first"),
            cell_run_manifest_path=("cell_run_manifest_path", "first"),
            nucleus_run_manifest_path=("nucleus_run_manifest_path", "first"),
            raw_imagej_results_path=("raw_imagej_results_path", "first"),
        )
        .reset_index()
        .merge(image_flags, on="filename", how="left", validate="one_to_one")
        .sort_values(["species", "filename"])
        .reset_index(drop=True)
    )
    image_trace["analysis_ready_image"] = image_trace["analysis_ready_image"].fillna(False).astype(bool)
    image_trace["cell_source_image_exists"] = image_trace["cell_source_image_path"].fillna("").astype(str).map(
        lambda value: Path(value).exists() if value else False
    )
    image_trace["nucleus_source_image_exists"] = image_trace["nucleus_source_image_path"].fillna("").astype(str).map(
        lambda value: Path(value).exists() if value else False
    )
    image_trace["source_image_exists"] = (
        image_trace["cell_source_image_exists"] & image_trace["nucleus_source_image_exists"]
    )
    image_trace["linked_pairs_source_path"] = str(linked_pairs_path.resolve())
    image_trace["linked_pairs_source_sha256"] = sha256_for_file(linked_pairs_path)

    trace_out = out_dir / "cellprofiler_linked_genome_image_trace.csv"
    image_trace.to_csv(trace_out, index=False)

    summary_rows.append(
        {
            "bundle": "genome_linked_yolo",
            "present": True,
            "n_species": int(image_trace["species"].dropna().nunique()),
            "n_images": int(image_trace["filename"].nunique()),
            "n_rows": int(len(strict_pairs)),
            "nonempty_source_image_path_pct": round(nonempty_pct(image_trace["nucleus_source_image_path"]), 2),
            "existing_source_image_path_pct": round(existing_pct(image_trace["nucleus_source_image_path"]), 2),
            "nonempty_tile_manifest_path_pct": round(nonempty_pct(image_trace["nucleus_tile_manifest_path"]), 2),
            "existing_tile_manifest_path_pct": round(existing_pct(image_trace["nucleus_tile_manifest_path"]), 2),
            "nonempty_mask_path_pct": round(nonempty_pct(strict_pairs["nucleus_mask_path"]), 2),
            "existing_mask_path_pct": round(existing_pct(strict_pairs["nucleus_mask_path"]), 2),
            "nonempty_roi_zip_path_pct": round(nonempty_pct(strict_pairs["roi_zip_path"]), 2),
            "existing_roi_zip_path_pct": round(existing_pct(strict_pairs["roi_zip_path"]), 2),
            "notes": "Linked strict-core YOLO nuclei are the authoritative publication nucleus/IOD trace.",
        }
    )

    return {
        "present": True,
        "origin": "upstream_mixed_linkage_run",
        "image_trace_path": str(trace_out.resolve()),
        "image_trace_sha256": sha256_for_file(trace_out),
        "source_files": [
            file_record("linked_nucleus_pairs", linked_pairs_path),
            file_record("image_summary", image_summary_path),
            file_record("species_summary", species_summary_path),
            file_record("summary_json", summary_json_path),
        ],
    }


def build_raw_genome_trace(
    cellprofiler_root: Path,
    out_dir: Path,
    metadata: pd.DataFrame,
    summary_rows: list[dict[str, object]],
    gap_rows: list[dict[str, object]],
    raw_run_tag: str,
) -> dict[str, object]:
    run_root = get_run_paths(cellprofiler_root, raw_run_tag).nucleus_iod
    run_specs = [
        {
            "bundle": "genome_raw_brightfield",
            "image_type": "brightfield",
            "measurement_path": run_root / "brightfield" / "measurements" / "nucleus_iod_measurements.csv",
            "image_index_path": run_root / "brightfield" / "image_index.csv",
            "run_manifest_path": run_root / "brightfield" / "run_manifest.json",
            "expected_preservation": "Dried Blood",
        },
        {
            "bundle": "genome_raw_pmount",
            "image_type": "pmount",
            "measurement_path": run_root / "pmount" / "measurements" / "nucleus_iod_measurements.csv",
            "image_index_path": run_root / "pmount" / "image_index.csv",
            "run_manifest_path": run_root / "pmount" / "run_manifest.json",
            "expected_preservation": pd.NA,
        },
    ]

    trace_frames: list[pd.DataFrame] = []
    source_files: list[dict[str, object]] = []
    for spec in run_specs:
        measurement_path = spec["measurement_path"]
        image_index_path = spec["image_index_path"]
        run_manifest_path = spec["run_manifest_path"]
        required = [measurement_path, image_index_path, run_manifest_path]
        missing = _required_paths_missing(required)
        if missing:
            gap_rows.append(
                {
                    "bundle": spec["bundle"],
                    "gap_type": "missing_raw_genome_run",
                    "species": pd.NA,
                    "details": f"One or more raw nucleus-IOD run files are missing: {missing}",
                }
            )
            continue

        measurements = pd.read_csv(measurement_path)
        image_index = pd.read_csv(image_index_path)
        measurements["filename"] = measurements["filename"].astype(str)
        measurements["species"] = measurements["species"].map(canonical_species)
        image_index["filename"] = image_index["filename"].astype(str)

        merged = measurements.merge(
            metadata[
                [
                    "filename",
                    "species",
                    "specimen_id",
                    "preservation",
                    "mounting",
                    "image_type",
                ]
            ],
            on=["filename", "species", "specimen_id", "image_type"],
            how="left",
            validate="many_to_one",
        )

        image_trace = (
            merged.groupby(
                ["species", "filename", "specimen_id", "preservation", "mounting", "image_type"],
                dropna=False,
            )
            .agg(
                n_nuclei=("iod", "size"),
                median_nucleus_area_um2=("area_um2", "median"),
                median_nucleus_iod=("iod", "median"),
            )
            .reset_index()
            .merge(
                image_index[
                    [
                        "filename",
                        "source_image_path",
                        "artifact_dir",
                        "tile_manifest_path",
                        "raw_imagej_results_path",
                        "status",
                    ]
                ],
                on="filename",
                how="left",
                validate="one_to_one",
            )
            .sort_values(["species", "filename"])
            .reset_index(drop=True)
        )
        image_trace["bundle"] = spec["bundle"]
        image_trace["measurement_source_path"] = str(measurement_path.resolve())
        image_trace["measurement_source_sha256"] = sha256_for_file(measurement_path)
        image_trace["image_index_path"] = str(image_index_path.resolve())
        image_trace["image_index_sha256"] = sha256_for_file(image_index_path)
        image_trace["run_manifest_path"] = str(run_manifest_path.resolve())
        image_trace["run_manifest_sha256"] = sha256_for_file(run_manifest_path)

        trace_frames.append(image_trace)
        source_files.extend(
            [
                file_record(f"{spec['bundle']}_measurements", measurement_path),
                file_record(f"{spec['bundle']}_image_index", image_index_path),
                file_record(f"{spec['bundle']}_run_manifest", run_manifest_path),
            ]
        )

        summary_rows.append(
            {
                "bundle": spec["bundle"],
                "present": True,
                "n_species": int(image_trace["species"].dropna().nunique()),
                "n_images": int(image_trace["filename"].nunique()),
                "n_rows": int(len(merged)),
                "nonempty_source_image_path_pct": round(nonempty_pct(image_trace["source_image_path"]), 2),
                "existing_source_image_path_pct": round(existing_pct(image_trace["source_image_path"]), 2),
                "nonempty_tile_manifest_path_pct": round(nonempty_pct(image_trace["tile_manifest_path"]), 2),
                "existing_tile_manifest_path_pct": round(existing_pct(image_trace["tile_manifest_path"]), 2),
                "nonempty_mask_path_pct": 0.0,
                "existing_mask_path_pct": 0.0,
                "nonempty_roi_zip_path_pct": 0.0,
                "existing_roi_zip_path_pct": 0.0,
                "notes": "Raw nucleus-IOD runs are traceable to image and tile-manifest level, but not to per-object mask or ROI artifacts.",
            }
        )

        source_exists_pct = existing_pct(image_trace["source_image_path"])
        raw_results_exists_pct = existing_pct(image_trace["raw_imagej_results_path"])
        gap_rows.append(
            {
                "bundle": spec["bundle"],
                "gap_type": "missing_object_mask_trace",
                "species": pd.NA,
                "details": "Current raw nucleus-IOD tables do not carry per-object mask_path or roi_zip_path fields.",
            }
        )
        if source_exists_pct < 100.0:
            gap_rows.append(
                {
                    "bundle": spec["bundle"],
                    "gap_type": "missing_source_image_files",
                    "species": pd.NA,
                    "details": "One or more source_image_path entries are populated in the raw nucleus-IOD index, but the referenced files are no longer present on disk.",
                }
            )
        if raw_results_exists_pct < 100.0:
            gap_rows.append(
                {
                    "bundle": spec["bundle"],
                    "gap_type": "missing_raw_imagej_results",
                    "species": pd.NA,
                    "details": "One or more raw_imagej_results_path entries are populated in the raw nucleus-IOD index, but the referenced files are missing on disk.",
                }
            )

        if spec["image_type"] == "brightfield":
            expected = set(
                metadata.loc[
                    metadata["image_type"].eq("brightfield") & metadata["preservation"].eq("Dried Blood"),
                    "species",
                ]
                .dropna()
                .map(canonical_species)
            )
            observed = set(image_trace["species"].dropna())
            for species in sorted(expected - observed):
                gap_rows.append(
                    {
                        "bundle": spec["bundle"],
                        "gap_type": "metadata_species_missing_from_run",
                        "species": species,
                        "details": "Brightfield metadata exists for this species, but no brightfield nucleus-IOD measurements were present in the current run.",
                    }
                )

    if trace_frames:
        combined = pd.concat(trace_frames, ignore_index=True).sort_values(["bundle", "species", "filename"])
        trace_out = out_dir / "cellprofiler_raw_genome_image_trace.csv"
        combined.to_csv(trace_out, index=False)
        return {
            "present": True,
            "origin": "upstream_raw_nucleus_iod_runs",
            "image_trace_path": str(trace_out.resolve()),
            "image_trace_sha256": sha256_for_file(trace_out),
            "source_files": source_files,
        }

    return {"present": False, "image_trace_path": pd.NA, "source_files": source_files}


def _round_or_na(value: object, digits: int = 4) -> float | pd.NA:
    if pd.isna(value):
        return pd.NA
    return round(float(value), digits)


def _bool_or_false(value: object) -> bool:
    if pd.isna(value):
        return False
    if isinstance(value, str):
        return value.strip().lower() in {"true", "1", "yes"}
    return bool(value)


def _sem_from_values(values: pd.Series) -> float | pd.NA:
    numeric = pd.to_numeric(values, errors="coerce").dropna()
    if len(numeric) < 2:
        return pd.NA
    return float(numeric.std(ddof=1) / (len(numeric) ** 0.5))


def _cv_pct(values: pd.Series) -> float | pd.NA:
    numeric = pd.to_numeric(values, errors="coerce").dropna()
    if len(numeric) < 2:
        return pd.NA
    mean = float(numeric.mean())
    if mean <= 0:
        return pd.NA
    return float(numeric.std(ddof=1) / mean * 100.0)


def _load_convergence_summary(path: Path) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame(columns=["species", "state_converged", "state_capped_not_converged"])
    df = pd.read_csv(path)
    df["species"] = df["species"].map(canonical_species)
    df["state_converged"] = df["converged"].fillna(False).astype(bool)
    df["state_capped_not_converged"] = df["capped_not_converged"].fillna(False).astype(bool)
    return df[["species", "state_converged", "state_capped_not_converged"]]


def _build_genome_state_summary(
    *,
    raw_trace_path: Path,
    cellprofiler_root: Path,
    out_dir: Path,
    gap_rows: list[dict[str, object]],
    raw_run_tag: str,
) -> tuple[pd.DataFrame, Path, list[dict[str, object]]]:
    raw_trace = pd.read_csv(raw_trace_path)
    if raw_trace.empty:
        gap_rows.append(
            {
                "bundle": "genome_species_bundle",
                "gap_type": "empty_raw_genome_trace",
                "species": pd.NA,
                "details": "The imported raw genome image trace is empty, so a raw-run species bundle cannot be reconstructed.",
            }
        )
        return pd.DataFrame(), out_dir / "cellprofiler_genome_state_summary.csv", []

    raw_trace["species"] = raw_trace["species"].map(canonical_species)
    raw_trace["n_nuclei"] = pd.to_numeric(raw_trace["n_nuclei"], errors="coerce")
    raw_trace["median_nucleus_area_um2"] = pd.to_numeric(raw_trace["median_nucleus_area_um2"], errors="coerce")
    raw_trace["median_nucleus_iod"] = pd.to_numeric(raw_trace["median_nucleus_iod"], errors="coerce")
    raw_trace["source_image_exists"] = raw_trace["source_image_path"].fillna("").astype(str).map(
        lambda value: Path(value).exists() if value else False
    )

    specimen_summary = (
        raw_trace.groupby(["species", "image_type", "specimen_id"], dropna=False)
        .agg(
            n_images=("filename", "nunique"),
            n_nuclei=("n_nuclei", "sum"),
            specimen_mean_area_um2=("median_nucleus_area_um2", "mean"),
            specimen_median_area_um2=("median_nucleus_area_um2", "median"),
            specimen_mean_iod=("median_nucleus_iod", "mean"),
            specimen_median_iod=("median_nucleus_iod", "median"),
        )
        .reset_index()
    )

    convergence_frames: list[pd.DataFrame] = []
    convergence_records: list[dict[str, object]] = []
    run_root = get_run_paths(cellprofiler_root, raw_run_tag).nucleus_iod
    for image_type in ["brightfield", "pmount"]:
        convergence_path = run_root / image_type / "convergence_summary.csv"
        convergence_records.append(file_record(f"{image_type}_convergence_summary", convergence_path))
        if convergence_path.exists():
            df = _load_convergence_summary(convergence_path)
            if not df.empty:
                df["image_type"] = image_type
                convergence_frames.append(df)

    convergence = (
        pd.concat(convergence_frames, ignore_index=True)
        if convergence_frames
        else pd.DataFrame(columns=["species", "image_type", "state_converged", "state_capped_not_converged"])
    )

    rows: list[dict[str, object]] = []
    for (species, image_type), image_group in raw_trace.groupby(["species", "image_type"], dropna=False):
        specimen_group = specimen_summary[
            specimen_summary["species"].eq(species) & specimen_summary["image_type"].eq(image_type)
        ].copy()
        if specimen_group.empty:
            continue

        specimen_units = specimen_group["specimen_mean_area_um2"].dropna()
        image_units = image_group["median_nucleus_area_um2"].dropna()
        if len(specimen_units) >= 2:
            support_unit = "specimen_mean_image_median_area_um2"
            se_area_um2 = _sem_from_values(specimen_units)
            cv_area_pct = _cv_pct(specimen_units)
        else:
            support_unit = "image_median_nucleus_area_um2"
            se_area_um2 = _sem_from_values(image_units)
            cv_area_pct = _cv_pct(image_units)

        row = {
            "species": species,
            "image_type": image_type,
            "n_images": int(image_group["filename"].nunique()),
            "n_specimens": int(specimen_group["specimen_id"].nunique()),
            "n_nuclei": int(pd.to_numeric(image_group["n_nuclei"], errors="coerce").fillna(0).sum()),
            "state_mean_area_um2": _round_or_na(specimen_group["specimen_mean_area_um2"].mean()),
            "state_median_area_um2": _round_or_na(specimen_group["specimen_median_area_um2"].median()),
            "state_mean_iod": _round_or_na(specimen_group["specimen_mean_iod"].mean()),
            "state_median_iod": _round_or_na(specimen_group["specimen_median_iod"].median()),
            "state_se_area_um2": _round_or_na(se_area_um2),
            "state_cv_area_pct": _round_or_na(cv_area_pct),
            "support_unit": support_unit,
            "source_image_existing_pct": _round_or_na(image_group["source_image_exists"].mean() * 100.0, 2),
            "tile_manifest_existing_pct": _round_or_na(existing_pct(image_group["tile_manifest_path"]), 2),
            "raw_imagej_results_existing_pct": _round_or_na(existing_pct(image_group["raw_imagej_results_path"]), 2),
        }
        rows.append(row)

    state_summary = pd.DataFrame(rows)
    if state_summary.empty:
        gap_rows.append(
            {
                "bundle": "genome_species_bundle",
                "gap_type": "missing_genome_state_summary",
                "species": pd.NA,
                "details": "No species/image-type state summaries could be derived from the raw genome trace.",
            }
        )
        return state_summary, out_dir / "cellprofiler_genome_state_summary.csv", convergence_records

    if not convergence.empty:
        state_summary = state_summary.merge(
            convergence,
            on=["species", "image_type"],
            how="left",
            validate="one_to_one",
        )
    state_summary["state_converged"] = state_summary["state_converged"].fillna(False).astype(bool)
    state_summary["state_capped_not_converged"] = state_summary["state_capped_not_converged"].fillna(False).astype(bool)

    ref_lookup = (
        state_summary.loc[state_summary["species"].eq("fuscus"), ["image_type", "state_mean_area_um2"]]
        .dropna(subset=["state_mean_area_um2"])
        .drop_duplicates(subset=["image_type"])
        .set_index("image_type")["state_mean_area_um2"]
        .to_dict()
    )

    missing_ref_states = sorted(set(state_summary["image_type"]) - set(ref_lookup))
    if missing_ref_states:
        gap_rows.append(
            {
                "bundle": "genome_species_bundle",
                "gap_type": "missing_reference_state",
                "species": "fuscus",
                "details": f"Reference species is missing one or more preparation states needed for raw-run scaling: {missing_ref_states}",
            }
        )

    genome_pg = []
    genome_se_pg = []
    genome_gb = []
    state_se_pct = []
    for _, row in state_summary.iterrows():
        ref_area = ref_lookup.get(row["image_type"])
        mean_area = row["state_mean_area_um2"]
        se_area = row["state_se_area_um2"]
        if ref_area is None or pd.isna(ref_area) or pd.isna(mean_area) or float(ref_area) <= 0:
            genome_pg.append(pd.NA)
            genome_se_pg.append(pd.NA)
            genome_gb.append(pd.NA)
            state_se_pct.append(pd.NA)
            continue
        estimate_pg = REFERENCE_GENOME_PG * float(mean_area) / float(ref_area)
        estimate_se_pg = pd.NA
        estimate_se_pct = pd.NA
        if pd.notna(se_area):
            estimate_se_pg = REFERENCE_GENOME_PG * float(se_area) / float(ref_area)
            if estimate_pg > 0:
                estimate_se_pct = float(estimate_se_pg) / float(estimate_pg) * 100.0
        genome_pg.append(_round_or_na(estimate_pg))
        genome_se_pg.append(_round_or_na(estimate_se_pg))
        genome_gb.append(_round_or_na(float(estimate_pg) * GENOME_GB_PER_PG))
        state_se_pct.append(_round_or_na(estimate_se_pct))

    state_summary["reference_species"] = "D. fuscus"
    state_summary["reference_genome_pg"] = REFERENCE_GENOME_PG
    state_summary["reference_state_mean_area_um2"] = state_summary["image_type"].map(ref_lookup).map(_round_or_na)
    state_summary["genome_pg_estimate"] = genome_pg
    state_summary["genome_se_pg_estimate"] = genome_se_pg
    state_summary["genome_gb_estimate"] = genome_gb
    state_summary["state_se_pct"] = state_se_pct
    state_summary = state_summary.sort_values(["species", "image_type"]).reset_index(drop=True)

    out_path = out_dir / "cellprofiler_genome_state_summary.csv"
    state_summary.to_csv(out_path, index=False)
    return state_summary, out_path, convergence_records


def _build_brightfield_only_genome_bundle(
    *,
    state_summary: pd.DataFrame,
    raw_trace_path: Path,
    state_summary_path: Path,
) -> pd.DataFrame:
    brightfield = state_summary[state_summary["image_type"].eq(GENOME_CALIBRATION_IMAGE_TYPE)].copy()
    if brightfield.empty:
        return brightfield

    bundle = pd.DataFrame(
        {
            "species": brightfield["species"].map(lambda value: f"D. {value}"),
            "primary_genome_pg": brightfield["genome_pg_estimate"],
            "primary_genome_se_pg": brightfield["genome_se_pg_estimate"],
            "primary_genome_gb": brightfield["genome_gb_estimate"],
            "primary_n_images": brightfield["n_images"].astype(int),
            "primary_n_specimens": brightfield["n_specimens"].astype(int),
            "primary_n_nuclei": brightfield["n_nuclei"].astype(int),
            "primary_cv_area_pct": brightfield["state_cv_area_pct"],
            "primary_state": brightfield["image_type"],
            "primary_state_converged": brightfield["state_converged"].astype(bool),
            "primary_state_capped_not_converged": brightfield["state_capped_not_converged"].astype(bool),
            "primary_support_unit": brightfield["support_unit"],
            "primary_source_image_existing_pct": brightfield["source_image_existing_pct"],
            "primary_tile_manifest_existing_pct": brightfield["tile_manifest_existing_pct"],
            "primary_raw_imagej_results_existing_pct": brightfield["raw_imagej_results_existing_pct"],
            "alternate_state": pd.NA,
            "alternate_genome_pg": pd.NA,
            "alternate_n_images": pd.NA,
            "alternate_n_specimens": pd.NA,
            "alternate_source_image_existing_pct": pd.NA,
            "cross_state_pct_diff": pd.NA,
            "primary_selection_reason": "brightfield_only_policy",
            "bundle_origin": "reconstructed_from_upstream_raw_runs",
            "reconstruction_method": "brightfield_only_specimen_mean_of_image_median_nucleus_area_scaled_to_fuscus",
            "reconstruction_source_path": str(raw_trace_path.resolve()),
            "reconstruction_source_sha256": sha256_for_file(raw_trace_path),
            "raw_image_trace_path": str(raw_trace_path.resolve()),
            "raw_image_trace_sha256": sha256_for_file(raw_trace_path),
            "genome_state_summary_path": str(state_summary_path.resolve()),
            "genome_state_summary_sha256": sha256_for_file(state_summary_path),
        }
    )

    statuses: list[str] = []
    flags: list[str] = []
    for _, row in bundle.iterrows():
        status, flag_summary = _flag_summary_for_primary_state(row)
        statuses.append(status)
        flags.append(flag_summary)
    bundle["result_status"] = statuses
    bundle["flag_summary"] = flags
    return bundle.sort_values("species").reset_index(drop=True)


def _build_linked_genome_state_summary(
    *,
    linked_trace_path: Path,
    out_dir: Path,
    gap_rows: list[dict[str, object]],
) -> tuple[pd.DataFrame, Path]:
    linked_trace = pd.read_csv(linked_trace_path)
    if linked_trace.empty:
        gap_rows.append(
            {
                "bundle": "genome_species_bundle",
                "gap_type": "empty_linked_genome_trace",
                "species": pd.NA,
                "details": "The imported linked genome image trace is empty, so a linked-run species bundle cannot be reconstructed.",
            }
        )
        return pd.DataFrame(), out_dir / "cellprofiler_genome_state_summary.csv"

    linked_trace["species"] = linked_trace["species"].map(canonical_species)
    linked_trace["n_nuclei"] = pd.to_numeric(linked_trace["n_nuclei"], errors="coerce")
    linked_trace["median_nucleus_area_um2"] = pd.to_numeric(linked_trace["median_nucleus_area_um2"], errors="coerce")
    linked_trace["median_nucleus_iod"] = pd.to_numeric(linked_trace["median_nucleus_iod"], errors="coerce")
    linked_trace["analysis_ready_image"] = (
        linked_trace["analysis_ready_image"].fillna(False).astype(str).str.lower().isin(["true", "1", "yes"])
    )
    linked_trace["source_image_exists"] = linked_trace["source_image_exists"].fillna(False).astype(bool)

    rows: list[dict[str, object]] = []
    for (species, image_type), image_group in linked_trace.groupby(["species", "image_type"], dropna=False):
        selected = image_group[image_group["analysis_ready_image"]].copy()
        support_tier = "analysis_ready_strict_core"
        if selected.empty:
            selected = image_group.copy()
            support_tier = "all_strict_core_images"

        specimen_summary = (
            selected.groupby(["species", "image_type", "specimen_id"], dropna=False)
            .agg(
                n_images=("filename", "nunique"),
                n_nuclei=("n_nuclei", "sum"),
                specimen_mean_area_um2=("median_nucleus_area_um2", "mean"),
                specimen_median_area_um2=("median_nucleus_area_um2", "median"),
                specimen_mean_iod=("median_nucleus_iod", "mean"),
                specimen_median_iod=("median_nucleus_iod", "median"),
            )
            .reset_index()
        )
        if specimen_summary.empty:
            continue

        specimen_iod_units = specimen_summary["specimen_mean_iod"].dropna()
        image_iod_units = selected["median_nucleus_iod"].dropna()
        if len(specimen_iod_units) >= 2:
            support_unit = "specimen_mean_image_median_nucleus_iod"
            se_measurement = _sem_from_values(specimen_iod_units)
            cv_measurement_pct = _cv_pct(specimen_iod_units)
        else:
            support_unit = "image_median_nucleus_iod"
            se_measurement = _sem_from_values(image_iod_units)
            cv_measurement_pct = _cv_pct(image_iod_units)

        row = {
            "species": species,
            "image_type": image_type,
            "n_images": int(selected["filename"].nunique()),
            "n_total_strict_images": int(image_group["filename"].nunique()),
            "n_analysis_ready_images": int(image_group["analysis_ready_image"].sum()),
            "n_specimens": int(specimen_summary["specimen_id"].nunique()),
            "n_nuclei": int(pd.to_numeric(selected["n_nuclei"], errors="coerce").fillna(0).sum()),
            "state_mean_area_um2": _round_or_na(specimen_summary["specimen_mean_area_um2"].mean()),
            "state_median_area_um2": _round_or_na(specimen_summary["specimen_median_area_um2"].median()),
            "state_mean_iod": _round_or_na(specimen_summary["specimen_mean_iod"].mean()),
            "state_median_iod": _round_or_na(specimen_summary["specimen_median_iod"].median()),
            "state_se_measurement": _round_or_na(se_measurement),
            "state_cv_measurement_pct": _round_or_na(cv_measurement_pct),
            "support_unit": support_unit,
            "support_tier": support_tier,
            "source_image_existing_pct": _round_or_na(selected["source_image_exists"].mean() * 100.0, 2),
            "tile_manifest_existing_pct": _round_or_na(existing_pct(selected["nucleus_tile_manifest_path"]), 2),
            "mask_existing_pct": _round_or_na(existing_pct(selected["nucleus_mask_path"]), 2),
            "roi_zip_existing_pct": _round_or_na(existing_pct(selected["roi_zip_path"]), 2),
        }
        rows.append(row)

    state_summary = pd.DataFrame(rows)
    if state_summary.empty:
        gap_rows.append(
            {
                "bundle": "genome_species_bundle",
                "gap_type": "missing_linked_genome_state_summary",
                "species": pd.NA,
                "details": "No species/image-type summaries could be derived from the linked strict-core genome trace.",
            }
        )
        return state_summary, out_dir / "cellprofiler_genome_state_summary.csv"

    ref_lookup = (
        state_summary.loc[state_summary["species"].eq("fuscus"), ["image_type", "state_mean_iod"]]
        .dropna(subset=["state_mean_iod"])
        .drop_duplicates(subset=["image_type"])
        .set_index("image_type")["state_mean_iod"]
        .to_dict()
    )

    missing_ref_states = sorted(set(state_summary["image_type"]) - set(ref_lookup))
    if missing_ref_states:
        gap_rows.append(
            {
                "bundle": "genome_species_bundle",
                "gap_type": "missing_linked_reference_state",
                "species": "fuscus",
                "details": f"Reference species is missing linked IOD summaries for one or more preparation states: {missing_ref_states}",
            }
        )

    genome_pg = []
    genome_se_pg = []
    genome_gb = []
    state_se_pct = []
    reference_signal = []
    for _, row in state_summary.iterrows():
        ref_iod = ref_lookup.get(row["image_type"])
        mean_iod = row["state_mean_iod"]
        se_measurement = row["state_se_measurement"]
        if ref_iod is None or pd.isna(ref_iod) or pd.isna(mean_iod) or float(ref_iod) <= 0:
            reference_signal.append(pd.NA)
            genome_pg.append(pd.NA)
            genome_se_pg.append(pd.NA)
            genome_gb.append(pd.NA)
            state_se_pct.append(pd.NA)
            continue
        estimate_pg = REFERENCE_GENOME_PG * float(mean_iod) / float(ref_iod)
        estimate_se_pg = pd.NA
        estimate_se_pct = pd.NA
        if pd.notna(se_measurement):
            estimate_se_pg = REFERENCE_GENOME_PG * float(se_measurement) / float(ref_iod)
            if estimate_pg > 0:
                estimate_se_pct = float(estimate_se_pg) / float(estimate_pg) * 100.0
        reference_signal.append(_round_or_na(ref_iod))
        genome_pg.append(_round_or_na(estimate_pg))
        genome_se_pg.append(_round_or_na(estimate_se_pg))
        genome_gb.append(_round_or_na(float(estimate_pg) * GENOME_GB_PER_PG))
        state_se_pct.append(_round_or_na(estimate_se_pct))

    state_summary["reference_species"] = "D. fuscus"
    state_summary["reference_genome_pg"] = REFERENCE_GENOME_PG
    state_summary["reference_state_mean_iod"] = reference_signal
    state_summary["genome_pg_estimate"] = genome_pg
    state_summary["genome_se_pg_estimate"] = genome_se_pg
    state_summary["genome_gb_estimate"] = genome_gb
    state_summary["state_se_pct"] = state_se_pct
    state_summary = state_summary.sort_values(["species", "image_type"]).reset_index(drop=True)

    out_path = out_dir / "cellprofiler_genome_state_summary.csv"
    state_summary.to_csv(out_path, index=False)
    return state_summary, out_path


def _build_brightfield_only_linked_genome_bundle(
    *,
    state_summary: pd.DataFrame,
    linked_trace_path: Path,
    state_summary_path: Path,
) -> pd.DataFrame:
    brightfield = state_summary[state_summary["image_type"].eq(GENOME_CALIBRATION_IMAGE_TYPE)].copy()
    if brightfield.empty:
        return brightfield

    bundle = pd.DataFrame(
        {
            "species": brightfield["species"].map(lambda value: f"D. {value}"),
            "primary_genome_pg": brightfield["genome_pg_estimate"],
            "primary_genome_se_pg": brightfield["genome_se_pg_estimate"],
            "primary_genome_gb": brightfield["genome_gb_estimate"],
            "primary_n_images": brightfield["n_images"].astype(int),
            "primary_n_specimens": brightfield["n_specimens"].astype(int),
            "primary_n_nuclei": brightfield["n_nuclei"].astype(int),
            "primary_n_analysis_ready_images": brightfield["n_analysis_ready_images"].astype(int),
            "primary_total_strict_images": brightfield["n_total_strict_images"].astype(int),
            "primary_cv_measurement_pct": brightfield["state_cv_measurement_pct"],
            "primary_cv_area_pct": brightfield["state_cv_measurement_pct"],
            "primary_measurement_kind": "linked_nucleus_iod",
            "primary_state": brightfield["image_type"],
            "primary_state_converged": False,
            "primary_state_capped_not_converged": False,
            "primary_support_unit": brightfield["support_unit"],
            "primary_support_tier": brightfield["support_tier"],
            "primary_source_image_existing_pct": brightfield["source_image_existing_pct"],
            "primary_tile_manifest_existing_pct": brightfield["tile_manifest_existing_pct"],
            "primary_mask_existing_pct": brightfield["mask_existing_pct"],
            "primary_roi_zip_existing_pct": brightfield["roi_zip_existing_pct"],
            "alternate_state": pd.NA,
            "alternate_genome_pg": pd.NA,
            "alternate_n_images": pd.NA,
            "alternate_n_specimens": pd.NA,
            "alternate_source_image_existing_pct": pd.NA,
            "cross_state_pct_diff": pd.NA,
            "primary_selection_reason": brightfield["support_tier"].map(
                lambda value: "analysis_ready_images_preferred"
                if value == "analysis_ready_strict_core"
                else "fallback_all_strict_core_images"
            ),
            "bundle_origin": "reconstructed_from_upstream_linked_runs",
            "reconstruction_method": "brightfield_only_analysis_ready_preferred_specimen_mean_of_image_median_linked_nucleus_iod_scaled_to_fuscus",
            "reconstruction_source_path": str(linked_trace_path.resolve()),
            "reconstruction_source_sha256": sha256_for_file(linked_trace_path),
            "raw_image_trace_path": str(linked_trace_path.resolve()),
            "raw_image_trace_sha256": sha256_for_file(linked_trace_path),
            "genome_state_summary_path": str(state_summary_path.resolve()),
            "genome_state_summary_sha256": sha256_for_file(state_summary_path),
        }
    )

    statuses: list[str] = []
    flags: list[str] = []
    for _, row in bundle.iterrows():
        status, flag_summary = _flag_summary_for_primary_state(row)
        statuses.append(status)
        flags.append(flag_summary)
    bundle["result_status"] = statuses
    bundle["flag_summary"] = flags
    return bundle.sort_values("species").reset_index(drop=True)


def _flag_summary_for_primary_state(row: pd.Series) -> tuple[str, str]:
    flags: list[str] = []
    se_pct = pd.to_numeric(pd.Series([row.get("state_se_pct")]), errors="coerce").iloc[0]
    cv_pct = pd.to_numeric(
        pd.Series(
            [
                row.get(
                    "primary_cv_measurement_pct",
                    row.get("state_cv_measurement_pct", row.get("state_cv_area_pct")),
                )
            ]
        ),
        errors="coerce",
    ).iloc[0]
    cross_state_pct_diff = pd.to_numeric(pd.Series([row.get("cross_state_pct_diff")]), errors="coerce").iloc[0]

    if row.get("primary_n_images", 0) < 3:
        flags.append("low_image_count")
    if row.get("primary_n_specimens", 0) < 2:
        flags.append("low_specimen_count")
    if row.get("primary_n_analysis_ready_images", 0) < 1:
        flags.append("no_analysis_ready_images")
    if pd.notna(se_pct) and se_pct >= 10.0:
        flags.append("high_standard_error")
    if pd.notna(cv_pct) and cv_pct >= 25.0:
        flags.append("high_measurement_cv")
    if row.get("primary_state_capped_not_converged", False):
        flags.append("state_not_converged")
    if row.get("primary_source_image_existing_pct", 0.0) < 100.0:
        flags.append("missing_source_image_files")
    if row.get("primary_mask_existing_pct", 100.0) < 100.0:
        flags.append("missing_nucleus_mask_files")
    if pd.notna(cross_state_pct_diff) and cross_state_pct_diff >= 10.0:
        flags.append("preservation_shift")

    if "state_not_converged" in flags or "preservation_shift" in flags:
        status = "sensitivity_limited"
    elif (row.get("primary_n_images", 0) < 2 and row.get("primary_n_specimens", 0) < 2) or (
        pd.notna(se_pct) and se_pct >= 20.0
    ):
        status = "sensitivity_limited"
    elif len(flags) >= 2:
        status = "caution"
    elif len(flags) == 1:
        status = "minor_caution"
    else:
        status = "stable"

    return status, "; ".join(flags)


def _selection_reason(primary: pd.Series, alternate: pd.Series | None) -> str:
    if alternate is None:
        return "only_available_state"
    if primary["n_specimens"] > alternate["n_specimens"]:
        return "more_specimens"
    if primary["n_images"] > alternate["n_images"]:
        return "more_images"
    if _bool_or_false(primary["state_converged"]) and not _bool_or_false(alternate["state_converged"]):
        return "converged_state"
    primary_source_pct = primary.get("source_image_existing_pct")
    alternate_source_pct = alternate.get("source_image_existing_pct")
    if (
        pd.notna(primary_source_pct)
        and pd.notna(alternate_source_pct)
        and float(primary_source_pct) > float(alternate_source_pct)
    ):
        return "better_source_image_retention"
    return "lower_standard_error"


def _reconstruct_genome_species_bundle_from_raw_runs(
    cellprofiler_root: Path,
    out_dir: Path,
    raw_genome_info: dict[str, object],
    summary_rows: list[dict[str, object]],
    gap_rows: list[dict[str, object]],
    raw_run_tag: str,
) -> dict[str, object]:
    raw_trace_path = Path(str(raw_genome_info.get("image_trace_path", "")))
    if not raw_trace_path.exists():
        gap_rows.append(
            {
                "bundle": "genome_species_bundle",
                "gap_type": "missing_raw_genome_trace",
                "species": pd.NA,
                "details": "The raw genome image trace is missing, so the genome species bundle cannot be rebuilt from upstream raw runs.",
            }
        )
        return {"present": False, "summary_path": pd.NA}

    state_summary, state_summary_path, convergence_records = _build_genome_state_summary(
        raw_trace_path=raw_trace_path,
        cellprofiler_root=cellprofiler_root,
        out_dir=out_dir,
        gap_rows=gap_rows,
        raw_run_tag=raw_run_tag,
    )
    if state_summary.empty:
        return {"present": False, "summary_path": pd.NA}

    bundle = _build_brightfield_only_genome_bundle(
        state_summary=state_summary,
        raw_trace_path=raw_trace_path,
        state_summary_path=state_summary_path,
    )
    if bundle.empty:
        gap_rows.append(
            {
                "bundle": "genome_species_bundle",
                "gap_type": "missing_brightfield_genome_states",
                "species": pd.NA,
                "details": "No brightfield genome states were available for the imported genome bundle reconstruction.",
            }
        )
        return {"present": False, "summary_path": pd.NA}

    out_path = out_dir / "cellprofiler_final_species_results.csv"
    bundle.to_csv(out_path, index=False)

    summary_rows.append(
        _genome_bundle_summary_row(
            bundle,
            "Reconstructed from the active brightfield raw nucleus-IOD runs using specimen-level image-median nucleus area scaled to D. fuscus.",
        )
    )
    gap_rows.append(
        {
            "bundle": "genome_species_bundle",
            "gap_type": "upstream_bundle_missing_reconstructed_from_raw_runs",
            "species": pd.NA,
            "details": "The legacy CellProfiler species bundle is missing, so the imported genome snapshot was rebuilt directly from the current brightfield raw nucleus-IOD runs.",
        }
    )

    audit_payload = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "bundle_path": str(out_path.resolve()),
        "bundle_sha256": sha256_for_file(out_path),
        "raw_image_trace_path": str(raw_trace_path.resolve()),
        "raw_image_trace_sha256": sha256_for_file(raw_trace_path),
        "genome_state_summary_path": str(state_summary_path.resolve()),
        "genome_state_summary_sha256": sha256_for_file(state_summary_path),
        "reference_genome_pg": REFERENCE_GENOME_PG,
        "reference_gb_per_pg": GENOME_GB_PER_PG,
        "selection_rule": "restrict to brightfield states only",
        "estimation_rule": "Scale specimen-level means of image-median brightfield nucleus area to D. fuscus; fall back to image-level SEM when only one specimen is available.",
        "calibration_image_type": GENOME_CALIBRATION_IMAGE_TYPE,
        "n_species": int(len(bundle)),
    }
    audit_path = out_dir / "cellprofiler_final_species_results_reconstruction.json"
    write_json(audit_path, audit_payload)

    source_files = list(raw_genome_info.get("source_files", [])) + convergence_records + [file_record("genome_state_summary", state_summary_path)]
    return {
        "present": True,
        "origin": "reconstructed_from_upstream_raw_runs",
        "summary_path": str(out_path.resolve()),
        "summary_sha256": sha256_for_file(out_path),
        "state_summary_path": str(state_summary_path.resolve()),
        "state_summary_sha256": sha256_for_file(state_summary_path),
        "reconstruction_audit_path": str(audit_path.resolve()),
        "reconstruction_audit_sha256": sha256_for_file(audit_path),
        "source_files": source_files,
        "recovered_from_downstream": False,
    }


def _reconstruct_genome_species_bundle_from_linked_runs(
    cellprofiler_root: Path,
    out_dir: Path,
    summary_rows: list[dict[str, object]],
    gap_rows: list[dict[str, object]],
    mixed_run_tag: str,
) -> dict[str, object]:
    linked_genome_info = build_linked_genome_trace(
        cellprofiler_root,
        out_dir,
        summary_rows,
        gap_rows,
        mixed_run_tag,
    )
    linked_trace_path = Path(str(linked_genome_info.get("image_trace_path", "")))
    if not linked_trace_path.exists():
        gap_rows.append(
            {
                "bundle": "genome_species_bundle",
                "gap_type": "missing_linked_genome_trace",
                "species": pd.NA,
                "details": "The linked genome image trace is missing, so the genome species bundle cannot be rebuilt from the mixed YOLO linkage run.",
            }
        )
        return {"present": False, "summary_path": pd.NA, "source_files": []}

    state_summary, state_summary_path = _build_linked_genome_state_summary(
        linked_trace_path=linked_trace_path,
        out_dir=out_dir,
        gap_rows=gap_rows,
    )
    if state_summary.empty:
        return {"present": False, "summary_path": pd.NA, "source_files": linked_genome_info.get("source_files", [])}

    bundle = _build_brightfield_only_linked_genome_bundle(
        state_summary=state_summary,
        linked_trace_path=linked_trace_path,
        state_summary_path=state_summary_path,
    )
    if bundle.empty:
        gap_rows.append(
            {
                "bundle": "genome_species_bundle",
                "gap_type": "missing_linked_brightfield_states",
                "species": pd.NA,
                "details": "No brightfield linked nucleus-IOD states were available for the imported genome bundle reconstruction.",
            }
        )
        return {"present": False, "summary_path": pd.NA, "source_files": linked_genome_info.get("source_files", [])}

    out_path = out_dir / "cellprofiler_final_species_results.csv"
    bundle.to_csv(out_path, index=False)

    summary_rows.append(
        _genome_bundle_summary_row(
            bundle,
            "Reconstructed from strict-core linked YOLO nuclei, preferring analysis-ready images, using specimen-level image-median nucleus IOD scaled to D. fuscus.",
        )
    )

    audit_payload = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "bundle_path": str(out_path.resolve()),
        "bundle_sha256": sha256_for_file(out_path),
        "linked_image_trace_path": str(linked_trace_path.resolve()),
        "linked_image_trace_sha256": sha256_for_file(linked_trace_path),
        "genome_state_summary_path": str(state_summary_path.resolve()),
        "genome_state_summary_sha256": sha256_for_file(state_summary_path),
        "reference_genome_pg": REFERENCE_GENOME_PG,
        "reference_gb_per_pg": GENOME_GB_PER_PG,
        "selection_rule": "restrict to brightfield strict-core linked nuclei and prefer analysis-ready images when available",
        "estimation_rule": "Scale specimen-level means of image-median linked nucleus IOD to D. fuscus.",
        "calibration_image_type": GENOME_CALIBRATION_IMAGE_TYPE,
        "calibration_signal": "linked_nucleus_iod",
        "mixed_linkage_run_tag": mixed_run_tag,
        "n_species": int(len(bundle)),
    }
    audit_path = out_dir / "cellprofiler_final_species_results_reconstruction.json"
    write_json(audit_path, audit_payload)

    source_files = list(linked_genome_info.get("source_files", [])) + [
        file_record("genome_state_summary", state_summary_path)
    ]
    return {
        "present": True,
        "origin": "reconstructed_from_upstream_linked_runs",
        "summary_path": str(out_path.resolve()),
        "summary_sha256": sha256_for_file(out_path),
        "state_summary_path": str(state_summary_path.resolve()),
        "state_summary_sha256": sha256_for_file(state_summary_path),
        "reconstruction_audit_path": str(audit_path.resolve()),
        "reconstruction_audit_sha256": sha256_for_file(audit_path),
        "source_files": source_files,
        "recovered_from_downstream": False,
    }


def _genome_bundle_summary_row(bundle: pd.DataFrame, notes: str) -> dict[str, object]:
    image_total = pd.to_numeric(bundle.get("primary_n_images"), errors="coerce").fillna(0).sum()
    return {
        "bundle": "genome_species_bundle",
        "present": True,
        "n_species": int(bundle["species"].nunique()),
        "n_images": int(image_total),
        "n_rows": int(len(bundle)),
        "nonempty_source_image_path_pct": pd.NA,
        "existing_source_image_path_pct": pd.NA,
        "nonempty_tile_manifest_path_pct": pd.NA,
        "existing_tile_manifest_path_pct": pd.NA,
        "nonempty_mask_path_pct": pd.NA,
        "existing_mask_path_pct": pd.NA,
        "nonempty_roi_zip_path_pct": pd.NA,
        "existing_roi_zip_path_pct": pd.NA,
        "notes": notes,
    }


def _import_upstream_genome_species_bundle(
    upstream_path: Path,
    out_dir: Path,
    raw_genome_info: dict[str, object],
    summary_rows: list[dict[str, object]],
) -> dict[str, object]:
    bundle = pd.read_csv(upstream_path)
    bundle["bundle_origin"] = "upstream_cellprofiler_qc_report_blockbalanced"
    bundle["upstream_source_path"] = str(upstream_path.resolve())
    bundle["upstream_source_sha256"] = sha256_for_file(upstream_path)
    bundle["raw_image_trace_path"] = raw_genome_info.get("image_trace_path", pd.NA)
    bundle["raw_image_trace_sha256"] = raw_genome_info.get("image_trace_sha256", pd.NA)

    out_path = out_dir / "cellprofiler_final_species_results.csv"
    bundle.to_csv(out_path, index=False)

    summary_rows.append(
        _genome_bundle_summary_row(
            bundle,
            "Imported directly from the upstream CellProfiler species bundle.",
        )
    )

    return {
        "present": True,
        "origin": "upstream_cellprofiler_qc_report_blockbalanced",
        "summary_path": str(out_path.resolve()),
        "summary_sha256": sha256_for_file(out_path),
        "source_files": [file_record("upstream_genome_species_bundle", upstream_path)],
        "recovered_from_downstream": False,
    }


def _recover_genome_species_bundle(
    out_dir: Path,
    raw_genome_info: dict[str, object],
    summary_rows: list[dict[str, object]],
    gap_rows: list[dict[str, object]],
) -> dict[str, object]:
    master_path = PROJECT_ROOT / "path_analysis" / "data" / "derived" / "master_species_table.csv"
    path_input_path = PROJECT_ROOT / "path_analysis" / "data" / "derived" / "path_input_master.csv"

    if not master_path.exists():
        gap_rows.append(
            {
                "bundle": "genome_species_bundle",
                "gap_type": "missing_recovery_source_table",
                "species": pd.NA,
                "details": "master_species_table.csv is missing, so the genome species bundle cannot be recovered locally.",
            }
        )
        return {"present": False, "summary_path": pd.NA}

    master = pd.read_csv(master_path)
    genome = (
        master[
            [
                "species",
                "genome_size_pg",
                "genome_size_se_pg",
                "genome_size_gb",
                "genome_n_images",
                "genome_n_specimens",
                "genome_n_nuclei",
                "genome_cv_area_pct",
                "genome_result_status",
                "genome_flag_summary",
            ]
        ]
        .dropna(subset=["genome_size_pg"])
        .drop_duplicates(subset=["species"])
        .sort_values("species")
        .reset_index(drop=True)
    )

    verified_against_path_input = False
    if path_input_path.exists():
        path_input = pd.read_csv(path_input_path)
        compare_cols = [
            "species",
            "genome_size_pg",
            "genome_size_se_pg",
            "genome_size_gb",
            "genome_result_status",
            "genome_flag_summary",
        ]
        overlap = (
            path_input[compare_cols]
            .dropna(subset=["genome_size_pg"])
            .drop_duplicates(subset=["species"])
            .sort_values("species")
            .reset_index(drop=True)
        )
        verified_against_path_input = genome[compare_cols].equals(overlap[compare_cols])

    bundle = pd.DataFrame(
        {
            "species": genome["species"].map(lambda value: f"D. {value}"),
            "primary_genome_pg": genome["genome_size_pg"],
            "primary_genome_se_pg": genome["genome_size_se_pg"],
            "primary_genome_gb": genome["genome_size_gb"],
            "primary_n_images": genome["genome_n_images"],
            "primary_n_specimens": genome["genome_n_specimens"],
            "primary_n_nuclei": genome["genome_n_nuclei"],
            "primary_cv_area_pct": genome["genome_cv_area_pct"],
            "result_status": genome["genome_result_status"],
            "flag_summary": genome["genome_flag_summary"],
            "bundle_origin": "recovered_from_path_analysis_tables",
            "reconstruction_source_path": str(master_path.resolve()),
            "reconstruction_source_sha256": sha256_for_file(master_path),
            "reconstruction_verified_against_path_input_master": verified_against_path_input,
            "raw_image_trace_path": raw_genome_info.get("image_trace_path", pd.NA),
            "raw_image_trace_sha256": raw_genome_info.get("image_trace_sha256", pd.NA),
        }
    )

    out_path = out_dir / "cellprofiler_final_species_results.csv"
    bundle.to_csv(out_path, index=False)

    summary_rows.append(
        _genome_bundle_summary_row(
            bundle,
            "Recovered locally from existing path-analysis tables because the upstream CellProfiler species bundle is missing.",
        )
    )
    gap_rows.append(
        {
            "bundle": "genome_species_bundle",
            "gap_type": "upstream_bundle_missing_recovered_locally",
            "species": pd.NA,
            "details": "The upstream cellprofiler_test species bundle is missing, so this local snapshot was reconstructed from current path-analysis tables.",
        }
    )

    audit_payload = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "bundle_path": str(out_path.resolve()),
        "bundle_sha256": sha256_for_file(out_path),
        "reconstruction_source_path": str(master_path.resolve()),
        "reconstruction_source_sha256": sha256_for_file(master_path),
        "verified_against_path_input_master": verified_against_path_input,
        "raw_image_trace_path": raw_genome_info.get("image_trace_path"),
        "raw_image_trace_sha256": raw_genome_info.get("image_trace_sha256"),
        "n_species": int(len(bundle)),
    }
    audit_path = out_dir / "cellprofiler_final_species_results_reconstruction.json"
    write_json(audit_path, audit_payload)

    return {
        "present": True,
        "origin": "recovered_from_path_analysis_tables",
        "summary_path": str(out_path.resolve()),
        "summary_sha256": sha256_for_file(out_path),
        "reconstruction_audit_path": str(audit_path.resolve()),
        "reconstruction_audit_sha256": sha256_for_file(audit_path),
        "source_files": [file_record("reconstruction_source_table", master_path)],
        "recovered_from_downstream": True,
    }


def build_genome_species_bundle(
    cellprofiler_root: Path,
    out_dir: Path,
    summary_rows: list[dict[str, object]],
    gap_rows: list[dict[str, object]],
    mixed_run_tag: str,
) -> dict[str, object]:
    return _reconstruct_genome_species_bundle_from_linked_runs(
        cellprofiler_root,
        out_dir,
        summary_rows,
        gap_rows,
        mixed_run_tag,
    )


def write_discovery_file(
    cellprofiler_root: Path,
    out_dir: Path,
    genome_info: dict[str, object],
    morphology_info: dict[str, object],
    raw_genome_info: dict[str, object],
    raw_run_tag: str,
    mixed_run_tag: str,
) -> Path:
    expected_genome_bundle_path = cellprofiler_root / "output" / "qc_report_blockbalanced" / "final_species_results.csv"
    active_inventory_path = cellprofiler_root / "docs" / "ACTIVE_WORKFLOW_INVENTORY.md"

    payload = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "cellprofiler_root": str(cellprofiler_root.resolve()),
        "reference_genome_pg": REFERENCE_GENOME_PG,
        "active_workflow_inventory": file_record("active_workflow_inventory", active_inventory_path),
        "active_run_tags": {
            "raw_auxiliary_run": raw_run_tag,
            "mixed_linkage_run": mixed_run_tag,
        },
        "expected_genome_species_bundle": file_record("expected_genome_species_bundle", expected_genome_bundle_path),
        "genome_species_bundle": genome_info,
        "morphology_bundle": morphology_info,
        "raw_genome_trace": raw_genome_info,
        "native_bridge_modules": [
            file_record("bridge", Path(__file__)),
            file_record("pipeline_runs", Path(__file__).with_name("pipeline_runs.py")),
            file_record("convergence", Path(__file__).with_name("convergence.py")),
        ],
    }
    out_path = out_dir / "cellprofiler_source_discovery.json"
    write_json(out_path, payload)
    return out_path


def rebuild_imported_artifacts(
    *,
    cellprofiler_root: Path,
    out_dir: Path,
    raw_run_tag: str = DEFAULT_ACTIVE_RAW_RUN_TAG,
    mixed_run_tag: str = DEFAULT_ACTIVE_MIXED_RUN_TAG,
) -> dict[str, object]:
    cellprofiler_root = cellprofiler_root.resolve()
    out_dir = out_dir.resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    summary_rows: list[dict[str, object]] = []
    gap_rows: list[dict[str, object]] = []

    raw_genome_info = {
        "present": False,
        "origin": "not_used_in_linked_yolo_publication_default",
        "details": "The publication default derives the genome bundle from linked YOLO nuclei rather than the standalone raw nucleus-IOD run.",
        "source_files": [],
    }
    genome_info = build_genome_species_bundle(
        cellprofiler_root,
        out_dir,
        summary_rows,
        gap_rows,
        mixed_run_tag,
    )
    morphology_info = build_morphology_bundle(
        cellprofiler_root,
        out_dir,
        summary_rows,
        gap_rows,
        mixed_run_tag,
    )

    discovery_path = write_discovery_file(
        cellprofiler_root,
        out_dir,
        genome_info,
        morphology_info,
        raw_genome_info,
        raw_run_tag,
        mixed_run_tag,
    )

    summary = pd.DataFrame(summary_rows).sort_values("bundle").reset_index(drop=True)
    gaps = (
        pd.DataFrame(gap_rows)
        .sort_values(["bundle", "gap_type", "species"], na_position="last")
        .reset_index(drop=True)
    )

    summary_path = out_dir / "cellprofiler_traceability_audit_summary.csv"
    gaps_path = out_dir / "cellprofiler_traceability_audit_gaps.csv"
    summary.to_csv(summary_path, index=False)
    gaps.to_csv(gaps_path, index=False)

    return {
        "discovery_path": str(discovery_path.resolve()),
        "summary_path": str(summary_path.resolve()),
        "gaps_path": str(gaps_path.resolve()),
        "morphology_present": bool(morphology_info.get("present")),
        "raw_genome_trace_present": bool(raw_genome_info.get("present")),
        "genome_species_bundle_present": bool(genome_info.get("present")),
    }
