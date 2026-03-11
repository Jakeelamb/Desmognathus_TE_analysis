#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import importlib
import json
import os
import shutil
import subprocess
import sys
import time
import traceback
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pandas as pd
import tifffile


PROJECT_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_CELLPROFILER_ROOT = Path.home() / "Projects" / "cellprofiler_test"
DEFAULT_MANIFEST = DEFAULT_CELLPROFILER_ROOT / "data" / "manifests" / "full_dataset_v1.csv"
DEFAULT_CELL_MODEL = (
    DEFAULT_CELLPROFILER_ROOT
    / "output"
    / "tile_training_round_v1"
    / "train"
    / "models"
    / "desmognathus_tile_round1"
)
DEFAULT_NUCLEUS_MODEL = (
    DEFAULT_CELLPROFILER_ROOT
    / "runs"
    / "segment"
    / "output"
    / "yolo_nucleus_training_final_area500_shape"
    / "weights"
    / "best.pt"
)
DEFAULT_YOLO_PYTHON = DEFAULT_CELLPROFILER_ROOT / ".venv-yolo" / "bin" / "python"
DEFAULT_LOG_ROOT = PROJECT_ROOT / "path_analysis" / "logs" / "publication_rebuild"
DEFAULT_OPTIMIZED_CELL_RUNNER = PROJECT_ROOT / "path_analysis" / "scripts" / "run_optimized_cell_size_from_manifest.py"

PATH_MODEL_RUNS: list[tuple[str, str]] = [
    ("te_genome", "te_genome_primary_mediumplus"),
    ("te_genome", "te_genome_primary_phylofill"),
    ("te_genome", "te_genome_primary_strict_body"),
    ("te_genome_organismal", "te_genome_organismal_primary_mediumplus"),
    ("te_genome_organismal", "te_genome_organismal_primary_phylofill"),
    ("te_genome_organismal", "te_genome_organismal_primary_strict_body"),
    ("te_genome_ltr_history", "te_genome_ltr_history_primary_mediumplus"),
    ("te_genome_ltr_history", "te_genome_ltr_history_primary_strict_body"),
    ("te_genome_ectopic", "te_genome_ectopic_primary_mediumplus"),
    ("te_genome_ectopic", "te_genome_ectopic_primary_phylofill"),
    ("te_genome_ectopic", "te_genome_ectopic_primary_strict_body"),
    ("te_genome_ectopic_organismal", "te_genome_ectopic_organismal_primary_mediumplus"),
    ("te_genome_ectopic_organismal", "te_genome_ectopic_organismal_primary_phylofill"),
    ("te_genome_ectopic_organismal", "te_genome_ectopic_organismal_primary_strict_body"),
    ("te_genome_morphology", "te_genome_morphology_primary_mediumplus"),
    ("te_genome_morphology", "te_genome_morphology_primary_phylofill"),
    ("te_genome_morphology", "te_genome_morphology_primary_strict_body"),
]


def utc_stamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")


def default_raw_run_tag(stamp: str) -> str:
    return f"full_dataset_pubrebuild_{stamp}"


def default_mixed_run_tag(stamp: str) -> str:
    return f"mixed_cellpose_yolo_pubrebuild_{stamp}"


def git_hash(path: Path) -> str:
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "--short", "HEAD"],
            cwd=path,
            text=True,
            timeout=10,
        ).strip()
    except Exception:
        return ""


def parse_args() -> argparse.Namespace:
    stamp = utc_stamp()
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cellprofiler-root", type=Path, default=DEFAULT_CELLPROFILER_ROOT)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--cell-model", type=Path, default=DEFAULT_CELL_MODEL)
    parser.add_argument("--nucleus-model", type=Path, default=DEFAULT_NUCLEUS_MODEL)
    parser.add_argument("--yolo-python", type=Path, default=DEFAULT_YOLO_PYTHON)
    parser.add_argument("--optimized-cell-runner", type=Path, default=DEFAULT_OPTIMIZED_CELL_RUNNER)
    parser.add_argument("--raw-run-tag", default=default_raw_run_tag(stamp))
    parser.add_argument("--mixed-run-tag", default=default_mixed_run_tag(stamp))
    parser.add_argument(
        "--reuse-cell-run-tag",
        default="",
        help="Reuse the existing cell-size segmentation run from this tag instead of rerunning Cellpose.",
    )
    parser.add_argument(
        "--reuse-nucleus-run-tag",
        default="",
        help="Reuse the existing nucleus-IOD run from this tag instead of rerunning ImageJ nucleus measurements.",
    )
    parser.add_argument(
        "--reuse-prepare-run-tag",
        default="",
        help="Reuse the prepare stage from this mixed run tag instead of rebuilding mixed linkage tiles.",
    )
    parser.add_argument(
        "--reuse-yolo-run-tag",
        default="",
        help="Reuse the YOLO nucleus-measurement stage from this mixed run tag instead of rerunning nucleus segmentation.",
    )
    parser.add_argument("--image-type", default="brightfield", choices=["brightfield"])
    parser.add_argument("--tile-format", default="png", choices=["png", "tiff"])
    parser.add_argument("--yolo-device", default=None)
    parser.add_argument("--yolo-batch", type=int, default=8)
    parser.add_argument("--yolo-conf", type=float, default=0.10)
    parser.add_argument("--yolo-iou", type=float, default=0.50)
    parser.add_argument("--yolo-max-det", type=int, default=512)
    parser.add_argument("--yolo-min-mask-area", type=int, default=500)
    parser.add_argument("--yolo-min-circularity", type=float, default=0.55)
    parser.add_argument("--yolo-min-solidity", type=float, default=0.92)
    parser.add_argument("--yolo-max-aspect-ratio", type=float, default=2.8)
    parser.add_argument("--nucleus-workers", type=int, default=1)
    parser.add_argument("--log-root", type=Path, default=DEFAULT_LOG_ROOT)
    parser.add_argument(
        "--wait-for-cell-run-complete",
        action="store_true",
        help="When reusing a cell run, wait until its image index covers the full manifest before continuing.",
    )
    parser.add_argument(
        "--wait-for-nucleus-run-complete",
        action="store_true",
        help="When reusing a nucleus-IOD run, wait until its completion artifacts are present before continuing.",
    )
    parser.add_argument(
        "--wait-poll-seconds",
        type=int,
        default=300,
        help="Polling interval in seconds while waiting for a reused cell run to finish.",
    )
    parser.add_argument(
        "--allow-existing-run-tags",
        action="store_true",
        help="Allow reusing existing upstream run directories instead of requiring a fresh tag.",
    )
    return parser.parse_args()


def validate_inputs(args: argparse.Namespace) -> None:
    required = [
        args.cellprofiler_root,
        args.manifest,
        args.cell_model,
        args.nucleus_model,
        args.yolo_python,
        args.optimized_cell_runner,
    ]
    missing = [str(path) for path in required if not Path(path).exists()]
    if missing:
        raise FileNotFoundError(f"Missing required inputs: {missing}")


def ensure_fresh_dir(path: Path, allow_existing: bool) -> None:
    if path.exists():
        if allow_existing:
            return
        if any(path.iterdir()):
            raise FileExistsError(f"Refusing to reuse non-empty directory: {path}")
    path.mkdir(parents=True, exist_ok=True)


def log_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")


def existing_nucleus_stage_dir(base_dir: Path, image_type: str) -> Path:
    candidate = base_dir / image_type
    return candidate if candidate.exists() else base_dir


def expected_manifest_rows(manifest: Path, image_type: str) -> int:
    rows = pd.read_csv(manifest)
    if image_type != "all" and "image_type" in rows.columns:
        rows = rows.loc[rows["image_type"].astype(str) == image_type].copy()
    return int(len(rows))


def wait_for_cell_run_completion(
    *,
    manifest: Path,
    image_type: str,
    cell_run_dir: Path,
    poll_seconds: int,
) -> None:
    expected_rows = expected_manifest_rows(manifest, image_type)
    index_path = cell_run_dir / "image_index.csv"
    print(
        f"[wait] Waiting for cell run completion in {cell_run_dir} "
        f"(expected_rows={expected_rows}, poll_seconds={poll_seconds})"
    )
    while True:
        if not index_path.exists():
            print(f"[wait] image_index.csv not present yet at {index_path}")
            time.sleep(poll_seconds)
            continue

        try:
            index = pd.read_csv(index_path)
        except Exception as exc:
            print(f"[wait] Could not read {index_path}: {exc}")
            time.sleep(poll_seconds)
            continue

        if index.empty:
            print(f"[wait] {index_path} is empty")
            time.sleep(poll_seconds)
            continue

        statuses = index.get("status", pd.Series(dtype="object")).astype(str)
        completed_rows = int(len(index))
        planned_rows = int((statuses == "planned").sum())
        pending_rows = max(expected_rows - completed_rows, 0)
        print(
            f"[wait] cell index rows={completed_rows}/{expected_rows}, "
            f"planned={planned_rows}, pending={pending_rows}"
        )
        if completed_rows >= expected_rows and planned_rows == 0:
            print("[wait] Cell run index is complete; continuing.")
            return
        time.sleep(poll_seconds)


def wait_for_nucleus_run_completion(
    *,
    nucleus_stage_dir: Path,
    poll_seconds: int,
) -> None:
    required = [
        nucleus_stage_dir / "measurements" / "nucleus_iod_measurements.csv",
        nucleus_stage_dir / "image_index.csv",
        nucleus_stage_dir / "run_manifest.json",
    ]
    print(
        f"[wait] Waiting for nucleus run completion in {nucleus_stage_dir} "
        f"(poll_seconds={poll_seconds})"
    )
    while True:
        missing = [path for path in required if not path.exists()]
        if not missing:
            print("[wait] Nucleus run completion artifacts are present; continuing.")
            return
        print("[wait] nucleus completion artifacts missing:", ", ".join(str(path) for path in missing))
        time.sleep(poll_seconds)


class RunRecorder:
    def __init__(self, log_dir: Path) -> None:
        self.log_dir = log_dir
        self.records: list[dict[str, Any]] = []

    def append(self, payload: dict[str, Any]) -> None:
        self.records.append(payload)
        log_json(self.log_dir / "stage_status.json", {"stages": self.records})


def run_logged_command(
    *,
    label: str,
    cmd: list[str],
    cwd: Path,
    log_path: Path,
    recorder: RunRecorder,
    env: dict[str, str] | None = None,
    continue_on_error: bool = False,
) -> int:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    started = time.time()
    print(f"\n[{label}] cwd={cwd}")
    print(f"[{label}] $ {' '.join(cmd)}")
    with log_path.open("w") as handle:
        handle.write(f"$ {' '.join(cmd)}\n\n")
        proc = subprocess.Popen(
            cmd,
            cwd=cwd,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            env=env,
            bufsize=1,
        )
        assert proc.stdout is not None
        for line in proc.stdout:
            sys.stdout.write(f"[{label}] {line}")
            handle.write(line)
        rc = proc.wait()
    elapsed = round(time.time() - started, 1)
    status = {
        "label": label,
        "command": cmd,
        "cwd": str(cwd),
        "log_path": str(log_path),
        "returncode": rc,
        "elapsed_seconds": elapsed,
        "continued_after_error": bool(rc != 0 and continue_on_error),
    }
    recorder.append(status)
    if rc != 0 and not continue_on_error:
        raise RuntimeError(f"{label} failed with exit code {rc}")
    return rc


def setup_cellprofiler_imports(cellprofiler_root: Path) -> None:
    additions = [
        cellprofiler_root / "src",
        cellprofiler_root / "scripts",
    ]
    for path in additions:
        text = str(path)
        if text not in sys.path:
            sys.path.insert(0, text)


def robust_prepare_tiles(
    *,
    cellprofiler_root: Path,
    cell_csv: Path,
    cell_run_manifest: Path,
    output_dir: Path,
    image_type: str,
    tile_format: str,
) -> dict[str, Any]:
    setup_cellprofiler_imports(cellprofiler_root)
    prep = importlib.import_module("prepare_mixed_linkage_tiles")

    ensure_fresh_dir(output_dir, allow_existing=False)
    tiles_dir = output_dir / "tiles"
    tiles_dir.mkdir(parents=True, exist_ok=True)

    cells = prep.load_cells(cell_csv, image_type)
    cells = prep.backfill_cell_traceability(cells, cell_run_manifest)
    tiles = prep.unique_tiles(cells)

    tile_rows: list[dict[str, Any]] = []
    error_rows: list[dict[str, Any]] = []
    grouped = list(tiles.groupby("source_image_path", sort=True))
    for idx, (source_path_str, group) in enumerate(grouped, start=1):
        source_path = Path(str(source_path_str))
        print(f"[prepare] [{idx}/{len(grouped)}] {source_path.name} ({len(group)} tiles)")
        try:
            chunk_rows = prep.extract_tiles(group.reset_index(drop=True), tiles_dir, tile_format)
            tile_rows.extend(chunk_rows)
        except Exception as exc:
            print(f"[prepare] ERROR {source_path.name}: {exc}")
            error_rows.append(
                {
                    "source_image_path": str(source_path),
                    "filename": str(group["filename"].iloc[0]) if len(group) else "",
                    "species": str(group["species"].iloc[0]) if len(group) else "",
                    "n_tiles_requested": int(len(group)),
                    "error": str(exc),
                }
            )

    if not tile_rows:
        raise RuntimeError("Tile preparation produced no successful rows")

    backfilled_csv = output_dir / "cell_measurements_backfilled.csv"
    cells.to_csv(backfilled_csv, index=False)
    tile_manifest = output_dir / "tile_manifest.csv"
    prep.write_tile_manifest(tile_manifest, tile_rows)

    error_csv = output_dir / "image_errors.csv"
    if error_rows:
        pd.DataFrame(error_rows).to_csv(error_csv, index=False)

    summary = {
        "cell_csv": str(cell_csv.resolve()),
        "cell_run_manifest": str(cell_run_manifest.resolve()),
        "output_dir": str(output_dir.resolve()),
        "image_type": image_type,
        "tile_format": tile_format,
        "n_cell_rows": int(len(cells)),
        "n_unique_tiles_requested": int(len(tiles)),
        "n_unique_images_requested": int(tiles["filename"].nunique()),
        "n_unique_tiles_prepared": int(len(tile_rows)),
        "n_images_failed": int(len(error_rows)),
        "cell_measurements_backfilled_csv": str(backfilled_csv.resolve()),
        "tile_manifest_csv": str(tile_manifest.resolve()),
        "image_errors_csv": str(error_csv.resolve()) if error_rows else "",
    }
    log_json(output_dir / "summary.json", summary)
    return summary


def robust_run_yolo_measurements(
    *,
    cellprofiler_root: Path,
    manifest_path: Path,
    model_path: Path,
    yolo_python: Path,
    output_dir: Path,
    background_cache: Path | None,
    object_kind: str,
    device_request: str | None,
    batch: int,
    conf: float,
    iou: float,
    max_det: int,
    min_mask_area: int,
    min_circularity: float,
    min_solidity: float,
    max_aspect_ratio: float,
) -> dict[str, Any]:
    ensure_fresh_dir(output_dir, allow_existing=False)
    with manifest_path.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        raise RuntimeError(f"Manifest has no rows: {manifest_path}")
    chunk_manifest_dir = output_dir / "chunk_manifests"
    chunk_manifest_dir.mkdir(parents=True, exist_ok=True)
    chunk_output_root = output_dir / "chunks"
    chunk_output_root.mkdir(parents=True, exist_ok=True)

    all_measurements: list[dict[str, Any]] = []
    tile_manifest_rows: list[dict[str, Any]] = []
    error_rows: list[dict[str, Any]] = []
    summary_rows: list[dict[str, Any]] = []

    grouped = list(pd.DataFrame(rows).groupby("filename", sort=True))
    for idx, (filename, group) in enumerate(grouped, start=1):
        image_stem = Path(str(filename)).stem
        chunk_manifest = chunk_manifest_dir / f"{image_stem}.csv"
        chunk_output_dir = chunk_output_root / image_stem
        chunk_tmp_dir = Path("/tmp") / f"{output_dir.name}_{image_stem}"
        group.to_csv(chunk_manifest, index=False)
        print(f"[yolo] [{idx}/{len(grouped)}] {filename} ({len(group)} tiles)")
        shutil.rmtree(chunk_tmp_dir, ignore_errors=True)

        cmd = [
            str(yolo_python),
            str(cellprofiler_root / "scripts" / "run_yolo_tile_measurements.py"),
            "--manifest",
            str(chunk_manifest),
            "--model",
            str(model_path),
            "--output-dir",
            str(chunk_tmp_dir),
            "--object-kind",
            object_kind,
            "--batch",
            str(batch),
            "--conf",
            str(conf),
            "--iou",
            str(iou),
            "--max-det",
            str(max_det),
            "--min-mask-area",
            str(min_mask_area),
            "--min-circularity",
            str(min_circularity),
            "--min-solidity",
            str(min_solidity),
            "--max-aspect-ratio",
            str(max_aspect_ratio),
        ]
        if background_cache is not None:
            cmd.extend(["--background-cache", str(background_cache)])
        if device_request:
            cmd.extend(["--device", device_request])

        # Do not resolve the interpreter symlink here; we need the venv root,
        # not the base interpreter path that the symlink may target.
        yolo_venv_root = yolo_python.parent.parent
        clean_env = os.environ.copy()
        for key in [
            "PYTHONHOME",
            "PYTHONPATH",
            "UV_ACTIVE",
            "UV_PROJECT",
            "UV_PYTHON",
            "CONDA_PREFIX",
            "CONDA_DEFAULT_ENV",
            "CONDA_SHLVL",
            "PYTHONEXECUTABLE",
        ]:
            clean_env.pop(key, None)
        clean_env["VIRTUAL_ENV"] = str(yolo_venv_root)
        clean_env["PATH"] = f"{yolo_venv_root / 'bin'}:{clean_env.get('PATH', '')}"
        result = subprocess.run(
            cmd,
            cwd=cellprofiler_root,
            capture_output=True,
            text=True,
            env=clean_env,
        )
        if result.returncode != 0:
            print(f"[yolo] ERROR {filename}: exit {result.returncode}")
            error_text = (result.stdout + "\n" + result.stderr).strip()
            if error_text:
                print(error_text[:1200])
            error_rows.append(
                {
                    "filename": str(filename),
                    "species": str(group["species"].iloc[0]) if len(group) else "",
                    "n_tiles_requested": int(len(group)),
                    "returncode": int(result.returncode),
                    "error": error_text[:4000],
                }
            )
            continue

        shutil.rmtree(chunk_output_dir, ignore_errors=True)
        shutil.copytree(chunk_tmp_dir, chunk_output_dir)
        chunk_tile_manifest = chunk_output_dir / "tile_manifest.csv"
        chunk_measurements = chunk_output_dir / "all_measurements.csv"
        chunk_summary = chunk_output_dir / "summary.json"
        if not chunk_tile_manifest.exists() or not chunk_measurements.exists() or not chunk_summary.exists():
            error_rows.append(
                {
                    "filename": str(filename),
                    "species": str(group["species"].iloc[0]) if len(group) else "",
                    "n_tiles_requested": int(len(group)),
                    "returncode": 0,
                    "error": "Chunk run completed but expected outputs were missing",
                }
            )
            continue

        with chunk_tile_manifest.open(newline="") as handle:
            tile_manifest_rows.extend(list(csv.DictReader(handle)))
        with chunk_measurements.open(newline="") as handle:
            all_measurements.extend(list(csv.DictReader(handle)))
        summary_rows.append(json.loads(chunk_summary.read_text()))

    error_csv = output_dir / "tile_errors.csv"
    if error_rows:
        pd.DataFrame(error_rows).to_csv(error_csv, index=False)

    if not tile_manifest_rows:
        raise RuntimeError("YOLO measurement stage produced no successful tiles")
    with (output_dir / "tile_manifest.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(tile_manifest_rows[0].keys()))
        writer.writeheader()
        writer.writerows(tile_manifest_rows)
    measurement_fieldnames = list(all_measurements[0].keys()) if all_measurements else []
    with (output_dir / "all_measurements.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=measurement_fieldnames)
        writer.writeheader()
        writer.writerows(all_measurements)

    summary = {
        "manifest": str(manifest_path.resolve()),
        "yolo_python": str(yolo_python.resolve()),
        "model": str(model_path.resolve()),
        "output_dir": str(output_dir.resolve()),
        "object_kind": object_kind,
        "background_cache": str(background_cache.resolve()) if background_cache is not None else "",
        "device": device_request or "auto",
        "n_images_requested": int(len(grouped)),
        "n_images_processed": int(len(summary_rows)),
        "n_images_failed": int(len(error_rows)),
        "n_tiles_requested": int(len(rows)),
        "n_tiles_processed": int(len(tile_manifest_rows)),
        "n_tiles_failed": int(len(error_rows)),
        "n_objects": int(len(all_measurements)),
        "n_total_predicted_instances": int(sum(int(row.get("n_total_predicted_instances", 0)) for row in summary_rows)),
        "n_total_positive_pixels": int(sum(int(row.get("n_total_positive_pixels", 0)) for row in summary_rows)),
        "min_mask_area": int(min_mask_area),
        "min_circularity": float(min_circularity),
        "min_solidity": float(min_solidity),
        "max_aspect_ratio": float(max_aspect_ratio),
        "tile_errors_csv": str(error_csv.resolve()) if error_rows else "",
    }
    log_json(output_dir / "summary.json", summary)
    return summary


def _path_readable(path_str: str, suffix: str = "") -> bool:
    value = str(path_str).strip()
    if not value:
        return False
    path = Path(value)
    if not path.exists():
        return False
    try:
        if suffix == ".zip":
            import zipfile

            with zipfile.ZipFile(path):
                return True
        with tifffile.TiffFile(path) as tif:
            _ = tif.pages[0].shape
        return True
    except Exception:
        return False


def robust_linkage(
    *,
    cellprofiler_root: Path,
    cell_csv: Path,
    nucleus_csv: Path,
    output_dir: Path,
    image_type: str,
    window_radius: int = 1,
) -> dict[str, Any]:
    setup_cellprofiler_imports(cellprofiler_root)
    linkage = importlib.import_module("build_cell_nucleus_linkage_report")

    ensure_fresh_dir(output_dir, allow_existing=False)

    cells = linkage.normalize_cells(linkage.load_measurements(cell_csv), image_type)
    nuclei = linkage.normalize_nuclei(linkage.load_measurements(nucleus_csv), image_type)

    cell_mask_ok = {
        path: _path_readable(path)
        for path in sorted({str(value).strip() for value in cells["cell_mask_path"] if str(value).strip()})
    }
    nucleus_mask_ok = {
        path: _path_readable(path)
        for path in sorted({str(value).strip() for value in nuclei["nucleus_mask_path"] if str(value).strip()})
    }
    roi_ok = {
        path: _path_readable(path, suffix=".zip")
        for path in sorted({str(value).strip() for value in nuclei["roi_zip_path"] if str(value).strip()})
    }

    cells = cells[cells["cell_mask_path"].map(lambda value: cell_mask_ok.get(str(value).strip(), False))].copy()
    nuclei = nuclei[
        nuclei["nucleus_mask_path"].map(lambda value: nucleus_mask_ok.get(str(value).strip(), False))
    ].copy()
    nuclei.loc[
        ~nuclei["roi_zip_path"].map(lambda value: roi_ok.get(str(value).strip(), False) if str(value).strip() else False),
        "roi_zip_path",
    ] = ""

    shared_files = sorted(set(cells["filename"]) & set(nuclei["filename"]))
    cells = cells[cells["filename"].isin(shared_files)].copy()
    nuclei = nuclei[nuclei["filename"].isin(shared_files)].copy()
    if cells.empty or nuclei.empty:
        raise RuntimeError("No shared cell/nucleus files remain after mask validation")

    linked, cell_linkage, link_summary = linkage.link_nuclei_to_cells(cells, nuclei, window_radius)
    linked = linkage.annotate_pair_qc(linked, cell_linkage)
    image_summary = linkage.summarize_images(linked, cells, nuclei)
    image_summary, image_qc_config = linkage.annotate_image_qc(image_summary)
    species_summary = linkage.summarize_species(image_summary, linked)
    summary = linkage.global_summary(linked, image_summary, cell_linkage, link_summary, image_qc_config)
    summary["shared_images"] = len(shared_files)
    summary["shared_species"] = int(cells["species"].nunique())
    summary["cell_csv"] = str(cell_csv.resolve())
    summary["nucleus_csv"] = str(nucleus_csv.resolve())
    summary["validated_cell_mask_paths"] = int(sum(cell_mask_ok.values()))
    summary["validated_nucleus_mask_paths"] = int(sum(nucleus_mask_ok.values()))
    summary["validated_roi_zip_paths"] = int(sum(roi_ok.values()))
    summary["rejected_cell_mask_paths"] = int(len(cell_mask_ok) - sum(cell_mask_ok.values()))
    summary["rejected_nucleus_mask_paths"] = int(len(nucleus_mask_ok) - sum(nucleus_mask_ok.values()))
    summary["rejected_roi_zip_paths"] = int(len(roi_ok) - sum(roi_ok.values()))

    linked.to_csv(output_dir / "linked_nucleus_pairs.csv.gz", index=False)
    cell_linkage.to_csv(output_dir / "cell_linkage_summary.csv.gz", index=False)
    image_summary.to_csv(output_dir / "image_summary.csv", index=False, float_format="%.6f")
    species_summary.to_csv(output_dir / "species_summary.csv", index=False, float_format="%.6f")
    (output_dir / "summary.json").write_text(json.dumps(summary, indent=2))
    linkage.write_methods_summary(output_dir, summary)

    strict = linked[linked["keep_strict_core"]].copy()
    linkage.plot_retention(summary, output_dir / "retention_by_tier.png")
    if len(strict):
        linkage.plot_scatter(
            strict,
            "cell_area_um2",
            "nuc_area_um2",
            "Cell area (um^2)",
            "Nucleus area (um^2)",
            "Strict-core cell area vs nucleus area",
            output_dir / "cell_vs_nucleus_area.png",
        )
        linkage.plot_scatter(
            strict,
            "cell_area_um2",
            "nuc_iod",
            "Cell area (um^2)",
            "Nucleus IOD",
            "Strict-core cell area vs nucleus IOD",
            output_dir / "cell_vs_nucleus_iod.png",
        )
        linkage.plot_ratio_hist(strict, output_dir / "nc_ratio_hist.png")

    html_text = linkage.build_html(output_dir, summary, image_summary, species_summary)
    (output_dir / "index.html").write_text(html_text)

    validation_payload = {
        "cell_mask_paths_checked": len(cell_mask_ok),
        "cell_mask_paths_valid": int(sum(cell_mask_ok.values())),
        "nucleus_mask_paths_checked": len(nucleus_mask_ok),
        "nucleus_mask_paths_valid": int(sum(nucleus_mask_ok.values())),
        "roi_zip_paths_checked": len(roi_ok),
        "roi_zip_paths_valid": int(sum(roi_ok.values())),
    }
    log_json(output_dir / "mask_validation_summary.json", validation_payload)
    return summary


def run_downstream_rebuild(
    *,
    args: argparse.Namespace,
    log_dir: Path,
    recorder: RunRecorder,
) -> None:
    commands: list[tuple[str, list[str], bool]] = [
        (
            "pull_cellprofiler_estimates",
            [
                "python3",
                "path_analysis/scripts/pull_cellprofiler_estimates.py",
                "--cellprofiler-root",
                str(args.cellprofiler_root),
                "--raw-run-tag",
                args.raw_run_tag,
                "--mixed-run-tag",
                args.mixed_run_tag,
            ],
            False,
        ),
        ("build_master_dataset", ["python3", "path_analysis/scripts/build_master_dataset.py"], False),
        ("prepare_te_features", ["python3", "path_analysis/scripts/prepare_te_features.py"], False),
        ("prepare_amphibio_traits", ["python3", "path_analysis/scripts/prepare_amphibio_traits.py"], False),
        ("prepare_external_morphometrics", ["python3", "path_analysis/scripts/prepare_external_morphometrics.py"], False),
        ("prepare_nc_biodiversity_traits", ["python3", "path_analysis/scripts/prepare_nc_biodiversity_traits.py"], False),
        ("prepare_curated_organismal_traits", ["python3", "path_analysis/scripts/prepare_curated_organismal_traits.py"], False),
        (
            "build_phylogenetic_trait_imputation",
            [
                "bash",
                "-lc",
                'source "$HOME/miniconda3/etc/profile.d/conda.sh" && conda activate Dusky'
                ' && export R_LIBS_USER='
                ' && export R_PROFILE_USER=/dev/null'
                ' && export R_ENVIRON_USER=/dev/null'
                ' && Rscript path_analysis/scripts/build_phylogenetic_trait_imputation.R',
            ],
            False,
        ),
        ("build_te_model_panel", ["python3", "path_analysis/scripts/build_te_model_panel.py"], False),
        ("build_path_input_master", ["python3", "path_analysis/scripts/build_path_input_master.py"], False),
        ("build_analysis_panels", ["python3", "path_analysis/scripts/build_analysis_panels.py"], False),
        ("audit_source_traceability", ["python3", "path_analysis/scripts/audit_source_traceability.py"], False),
        ("audit_publication_readiness", ["python3", "path_analysis/scripts/audit_publication_readiness.py"], False),
    ]

    for idx, (label, cmd, continue_on_error) in enumerate(commands, start=1):
        run_logged_command(
            label=f"{idx:02d}_{label}",
            cmd=cmd,
            cwd=PROJECT_ROOT,
            log_path=log_dir / f"{idx:02d}_{label}.log",
            recorder=recorder,
            continue_on_error=continue_on_error,
        )

    offset = len(commands)
    for pos, (family, panel) in enumerate(PATH_MODEL_RUNS, start=1):
        label = f"{offset + pos:02d}_path_model_{panel}"
        cmd = [
            "bash",
            "-lc",
            'source "$HOME/miniconda3/etc/profile.d/conda.sh" && conda activate Dusky'
            + " && export R_LIBS_USER="
            + " && export R_PROFILE_USER=/dev/null"
            + " && export R_ENVIRON_USER=/dev/null"
            + f" && Rscript path_analysis/scripts/path_model_scaffold.R --family {family} --panel {panel}",
        ]
        run_logged_command(
            label=label,
            cmd=cmd,
            cwd=PROJECT_ROOT,
            log_path=log_dir / f"{label}.log",
            recorder=recorder,
            continue_on_error=True,
        )

    final_cmds: list[tuple[str, list[str], bool]] = [
        ("audit_publication_readiness_post_models", ["python3", "path_analysis/scripts/audit_publication_readiness.py"], False),
        ("build_paper_freeze_manifest", ["python3", "scripts/processing/build_paper_freeze_manifest.py"], False),
    ]
    for pos, (label, cmd, continue_on_error) in enumerate(final_cmds, start=1):
        index = offset + len(PATH_MODEL_RUNS) + pos
        run_logged_command(
            label=f"{index:02d}_{label}",
            cmd=cmd,
            cwd=PROJECT_ROOT,
            log_path=log_dir / f"{index:02d}_{label}.log",
            recorder=recorder,
            continue_on_error=continue_on_error,
        )


def main() -> None:
    args = parse_args()
    validate_inputs(args)

    cellprofiler_root = args.cellprofiler_root.resolve()
    manifest = args.manifest.resolve()
    raw_run_root = cellprofiler_root / "output" / "runs" / args.raw_run_tag
    mixed_run_root = cellprofiler_root / "output" / "runs" / args.mixed_run_tag
    log_dir = args.log_root.resolve() / utc_stamp()
    log_dir.mkdir(parents=True, exist_ok=True)
    recorder = RunRecorder(log_dir)

    ensure_fresh_dir(raw_run_root, args.allow_existing_run_tags)
    ensure_fresh_dir(mixed_run_root, args.allow_existing_run_tags)

    gpu_available = False
    try:
        result = subprocess.run(
            ["nvidia-smi", "-L"],
            capture_output=True,
            text=True,
            timeout=10,
        )
        gpu_available = result.returncode == 0
    except Exception:
        gpu_available = False
    if not gpu_available:
        try:
            import torch

            gpu_available = bool(torch.cuda.is_available())
        except Exception:
            gpu_available = False

    context = {
        "started_at": datetime.now(timezone.utc).isoformat(),
        "project_root": str(PROJECT_ROOT),
        "project_git_hash": git_hash(PROJECT_ROOT),
        "cellprofiler_root": str(cellprofiler_root),
        "cellprofiler_git_hash": git_hash(cellprofiler_root),
        "manifest": str(manifest),
        "raw_run_tag": args.raw_run_tag,
        "mixed_run_tag": args.mixed_run_tag,
        "raw_run_root": str(raw_run_root),
        "mixed_run_root": str(mixed_run_root),
        "cell_model": str(args.cell_model.resolve()),
        "optimized_cell_runner": str(args.optimized_cell_runner.resolve()),
        "nucleus_model": str(args.nucleus_model.resolve()),
        "reuse_cell_run_tag": args.reuse_cell_run_tag,
        "reuse_nucleus_run_tag": args.reuse_nucleus_run_tag,
        "reuse_prepare_run_tag": args.reuse_prepare_run_tag,
        "reuse_yolo_run_tag": args.reuse_yolo_run_tag,
        "gpu_available": gpu_available,
    }
    log_json(log_dir / "context.json", context)
    print(json.dumps(context, indent=2))

    if args.reuse_cell_run_tag:
        raw_cell_dir = cellprofiler_root / "output" / "runs" / args.reuse_cell_run_tag / "cell_size_segmentation"
    else:
        raw_cell_dir = raw_run_root / "cell_size_segmentation"
    if args.reuse_nucleus_run_tag:
        raw_nucleus_dir = cellprofiler_root / "output" / "runs" / args.reuse_nucleus_run_tag / "nucleus_iod"
        raw_nucleus_stage_dir = existing_nucleus_stage_dir(raw_nucleus_dir, args.image_type)
    else:
        raw_nucleus_dir = raw_run_root / "nucleus_iod"
        raw_nucleus_stage_dir = raw_nucleus_dir / args.image_type
    if args.reuse_prepare_run_tag:
        mixed_prepare_dir = (
            cellprofiler_root / "output" / "runs" / args.reuse_prepare_run_tag / "prepare"
        )
    else:
        mixed_prepare_dir = mixed_run_root / "prepare"
    if args.reuse_yolo_run_tag:
        mixed_nucleus_dir = (
            cellprofiler_root / "output" / "runs" / args.reuse_yolo_run_tag / "nucleus_measurements"
        )
    else:
        mixed_nucleus_dir = mixed_run_root / "nucleus_measurements"
    mixed_linkage_dir = mixed_run_root / "linkage"

    if args.reuse_cell_run_tag:
        prewait_required = [raw_cell_dir / "run_manifest.json"]
        if not args.wait_for_cell_run_complete:
            prewait_required.append(raw_cell_dir / "all_measurements.csv")
        missing = [str(path) for path in prewait_required if not path.exists()]
        if missing:
            raise FileNotFoundError(f"Missing reused cell-size run inputs: {missing}")
        reuse_note = {
            "label": "01_reuse_cell_size_segmentation",
            "returncode": 0,
            "elapsed_seconds": 0.0,
            "reused_from_run_tag": args.reuse_cell_run_tag,
            "cell_run_dir": str(raw_cell_dir.resolve()),
        }
        recorder.append(reuse_note)
        log_json(log_dir / "01_reuse_cell_size_segmentation.json", reuse_note)
        if args.wait_for_cell_run_complete:
            wait_started = time.time()
            wait_for_cell_run_completion(
                manifest=manifest,
                image_type=args.image_type,
                cell_run_dir=raw_cell_dir,
                poll_seconds=max(int(args.wait_poll_seconds), 1),
            )
            recorder.append(
                {
                    "label": "01a_wait_for_reused_cell_run_completion",
                    "returncode": 0,
                    "elapsed_seconds": round(time.time() - wait_started, 1),
                    "cell_run_dir": str(raw_cell_dir.resolve()),
                }
            )
        postwait_required = [
            raw_cell_dir / "all_measurements.csv",
            raw_cell_dir / "run_manifest.json",
        ]
        missing = [str(path) for path in postwait_required if not path.exists()]
        if missing:
            raise FileNotFoundError(f"Missing completed reused cell-size run inputs: {missing}")
    else:
        cell_cmd = [
            sys.executable,
            str(args.optimized_cell_runner.resolve()),
            "--cellprofiler-root",
            str(cellprofiler_root),
            "--manifest",
            str(manifest),
            "--output-dir",
            str(raw_cell_dir),
            "--backend",
            "cellpose",
            "--cellpose-model",
            str(args.cell_model),
            "--image-type",
            args.image_type,
            "--resume",
            "--no-crops",
        ]
        if gpu_available:
            cell_cmd.append("--gpu")

        run_logged_command(
            label="01_cell_size_segmentation",
            cmd=cell_cmd,
            cwd=cellprofiler_root,
            log_path=log_dir / "01_cell_size_segmentation.log",
            recorder=recorder,
        )

    background_cache_path: Path | None = None
    if args.reuse_nucleus_run_tag:
        prewait_required: list[Path] = []
        if not args.wait_for_nucleus_run_complete:
            prewait_required = [
                raw_nucleus_stage_dir / "measurements" / "nucleus_iod_measurements.csv",
                raw_nucleus_stage_dir / "image_index.csv",
                raw_nucleus_stage_dir / "run_manifest.json",
            ]
        missing = [str(path) for path in prewait_required if not path.exists()]
        if missing:
            raise FileNotFoundError(f"Missing reused nucleus-IOD run inputs: {missing}")
        reuse_note = {
            "label": "02_reuse_nucleus_iod",
            "returncode": 0,
            "elapsed_seconds": 0.0,
            "reused_from_run_tag": args.reuse_nucleus_run_tag,
            "nucleus_run_dir": str(raw_nucleus_stage_dir.resolve()),
        }
        recorder.append(reuse_note)
        log_json(log_dir / "02_reuse_nucleus_iod.json", reuse_note)
        if args.wait_for_nucleus_run_complete:
            wait_started = time.time()
            wait_for_nucleus_run_completion(
                nucleus_stage_dir=raw_nucleus_stage_dir,
                poll_seconds=max(int(args.wait_poll_seconds), 1),
            )
            recorder.append(
                {
                    "label": "02a_wait_for_reused_nucleus_run_completion",
                    "returncode": 0,
                    "elapsed_seconds": round(time.time() - wait_started, 1),
                    "nucleus_run_dir": str(raw_nucleus_stage_dir.resolve()),
                }
            )
        postwait_required = [
            raw_nucleus_stage_dir / "measurements" / "nucleus_iod_measurements.csv",
            raw_nucleus_stage_dir / "image_index.csv",
            raw_nucleus_stage_dir / "run_manifest.json",
        ]
        missing = [str(path) for path in postwait_required if not path.exists()]
        if missing:
            raise FileNotFoundError(f"Missing completed reused nucleus-IOD run inputs: {missing}")
        background_cache_path = raw_nucleus_stage_dir / "measurements" / "nucleus_iod_measurements.csv"
    else:
        recorder.append(
            {
                "label": "02_linked_nucleus_policy",
                "returncode": 0,
                "elapsed_seconds": 0.0,
                "details": "Skipping standalone raw nucleus-IOD stage; publication default derives nucleus IOD from linked YOLO nuclei.",
            }
        )
        log_json(
            log_dir / "02_linked_nucleus_policy.json",
            {
                "policy": "linked_yolo_nucleus_iod_is_authoritative",
                "background_cache_path": "",
                "details": "No standalone background cache was supplied; YOLO tile measurements will compute i_bg directly from each tile.",
            },
        )

    if args.reuse_prepare_run_tag:
        required = [
            mixed_prepare_dir / "tile_manifest.csv",
            mixed_prepare_dir / "cell_measurements_backfilled.csv",
            mixed_prepare_dir / "summary.json",
        ]
        missing = [str(path) for path in required if not path.exists()]
        if missing:
            raise FileNotFoundError(f"Missing reused prepare inputs: {missing}")
        prepare_summary = json.loads((mixed_prepare_dir / "summary.json").read_text())
        reuse_note = {
            "label": "03_reuse_prepare_mixed_tiles",
            "returncode": 0,
            "elapsed_seconds": 0.0,
            "reused_from_run_tag": args.reuse_prepare_run_tag,
            "prepare_dir": str(mixed_prepare_dir.resolve()),
        }
        recorder.append(reuse_note)
        log_json(log_dir / "03_reuse_prepare_mixed_tiles.json", reuse_note)
    else:
        started = time.time()
        prepare_summary = robust_prepare_tiles(
            cellprofiler_root=cellprofiler_root,
            cell_csv=raw_cell_dir / "all_measurements.csv",
            cell_run_manifest=raw_cell_dir / "run_manifest.json",
            output_dir=mixed_prepare_dir,
            image_type=args.image_type,
            tile_format=args.tile_format,
        )
        recorder.append(
            {
                "label": "03_prepare_mixed_tiles",
                "returncode": 0,
                "elapsed_seconds": round(time.time() - started, 1),
                "summary_path": str((mixed_prepare_dir / "summary.json").resolve()),
            }
        )

    if args.reuse_yolo_run_tag:
        required = [
            mixed_nucleus_dir / "all_measurements.csv",
            mixed_nucleus_dir / "tile_manifest.csv",
            mixed_nucleus_dir / "summary.json",
        ]
        missing = [str(path) for path in required if not path.exists()]
        if missing:
            raise FileNotFoundError(f"Missing reused YOLO inputs: {missing}")
        yolo_summary = json.loads((mixed_nucleus_dir / "summary.json").read_text())
        reuse_note = {
            "label": "04_reuse_yolo_nucleus_measurements",
            "returncode": 0,
            "elapsed_seconds": 0.0,
            "reused_from_run_tag": args.reuse_yolo_run_tag,
            "nucleus_measurement_dir": str(mixed_nucleus_dir.resolve()),
        }
        recorder.append(reuse_note)
        log_json(log_dir / "04_reuse_yolo_nucleus_measurements.json", reuse_note)
    else:
        started = time.time()
        yolo_summary = robust_run_yolo_measurements(
            cellprofiler_root=cellprofiler_root,
            manifest_path=mixed_prepare_dir / "tile_manifest.csv",
            model_path=args.nucleus_model.resolve(),
            yolo_python=args.yolo_python.resolve(),
            output_dir=mixed_nucleus_dir,
            background_cache=background_cache_path,
            object_kind="nucleus",
            device_request=args.yolo_device,
            batch=args.yolo_batch,
            conf=args.yolo_conf,
            iou=args.yolo_iou,
            max_det=args.yolo_max_det,
            min_mask_area=args.yolo_min_mask_area,
            min_circularity=args.yolo_min_circularity,
            min_solidity=args.yolo_min_solidity,
            max_aspect_ratio=args.yolo_max_aspect_ratio,
        )
        recorder.append(
            {
                "label": "04_yolo_nucleus_measurements",
                "returncode": 0,
                "elapsed_seconds": round(time.time() - started, 1),
                "summary_path": str((mixed_nucleus_dir / "summary.json").resolve()),
            }
        )

    started = time.time()
    linkage_summary = robust_linkage(
        cellprofiler_root=cellprofiler_root,
        cell_csv=mixed_prepare_dir / "cell_measurements_backfilled.csv",
        nucleus_csv=mixed_nucleus_dir / "all_measurements.csv",
        output_dir=mixed_linkage_dir,
        image_type=args.image_type,
    )
    recorder.append(
        {
            "label": "05_cell_nucleus_linkage",
            "returncode": 0,
            "elapsed_seconds": round(time.time() - started, 1),
            "summary_path": str((mixed_linkage_dir / "summary.json").resolve()),
        }
    )

    run_downstream_rebuild(args=args, log_dir=log_dir, recorder=recorder)

    final_summary = {
        "completed_at": datetime.now(timezone.utc).isoformat(),
        "raw_run_tag": args.raw_run_tag,
        "mixed_run_tag": args.mixed_run_tag,
        "prepare_summary": prepare_summary,
        "yolo_summary": yolo_summary,
        "linkage_summary": linkage_summary,
        "log_dir": str(log_dir),
    }
    log_json(log_dir / "final_summary.json", final_summary)
    print(json.dumps(final_summary, indent=2))


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        error_payload = {
            "failed_at": datetime.now(timezone.utc).isoformat(),
            "error": str(exc),
            "traceback": traceback.format_exc(),
        }
        print(json.dumps(error_payload, indent=2), file=sys.stderr)
        raise
