#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import importlib
import json
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Any


DEFAULT_CELLPROFILER_ROOT = Path.home() / "Projects" / "cellprofiler_test"

CSV_COLUMNS = [
    "label",
    "seed_label",
    "area_px",
    "area_um2",
    "solidity",
    "circularity",
    "iod",
    "mean_od",
    "centroid_y",
    "centroid_x",
    "i_bg",
    "tile_name",
    "tile_y0",
    "tile_x0",
    "tile_h",
    "tile_w",
    "tile_score",
    "mask_path",
    "raw_mask_path",
    "mask_label_id",
    "tile_manifest_path",
    "overlay_path",
    "nucleus_area_px_seed",
    "nucleus_area_um2_seed",
    "nc_ratio_seed",
    "distance_over_cell_radius_seed",
    "cell_extent",
    "cell_edge_touch",
    "cell_threshold",
    "nucleus_threshold",
    "source_image_path",
    "run_manifest_path",
    "filename",
    "slide_id",
    "specimen_id",
    "species",
    "image_type",
]

CROP_FIELD_DEFAULTS = {
    "crop_filename": "",
    "crop_size_px": "",
    "centroid_in_crop_y": "",
    "centroid_in_crop_x": "",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cellprofiler-root", type=Path, default=DEFAULT_CELLPROFILER_ROOT)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--gpu", action="store_true")
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--debug", action="store_true")
    parser.add_argument("--no-crops", action="store_true")
    parser.add_argument("--tile-filter", choices=["green", "auto", "yellow", "all"], default="auto")
    parser.add_argument("--max-cells-image", type=int, default=500)
    parser.add_argument(
        "--max-cells-species",
        type=int,
        default=0,
        help="Stop scheduling more images for a species once this many cells have been collected; 0 disables the cap.",
    )
    parser.add_argument("--diameter", type=float, default=None)
    parser.add_argument("--backend", choices=["cellpose"], default="cellpose")
    parser.add_argument("--cellpose-model", type=Path, required=True)
    parser.add_argument("--image-type", choices=["brightfield", "pmount", "all"], default="brightfield")
    parser.add_argument("--gpu-tile-batch-images", type=int, default=1)
    parser.add_argument("--gpu-patch-batch-size", type=int, default=16)
    parser.add_argument("--cpu-tile-batch-images", type=int, default=1)
    parser.add_argument("--cpu-patch-batch-size", type=int, default=8)
    return parser.parse_args()


def setup_runtime(cellprofiler_root: Path) -> tuple[Any, Any, Any]:
    for extra in [
        cellprofiler_root / "src",
        cellprofiler_root / "scripts",
        cellprofiler_root / "Cellsize_segmentation_cellpose_pipeline" / "scripts",
    ]:
        text = str(extra.resolve())
        if text not in sys.path:
            sys.path.insert(0, text)

    pipeline_runs = importlib.import_module("cellprofiler_tools.pipeline_runs")
    convergence = importlib.import_module("cellprofiler_tools.convergence")
    segment_cells = importlib.import_module("segment_cells")
    return pipeline_runs, convergence, segment_cells


def append_rows(csv_path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        return
    write_header = not csv_path.exists() or csv_path.stat().st_size == 0
    with csv_path.open("a", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=CSV_COLUMNS, extrasaction="ignore")
        if write_header:
            writer.writeheader()
        writer.writerows(rows)


def write_index(index_path: Path, index_rows: list[dict[str, Any]]) -> None:
    with index_path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=index_rows[0].keys()
            if index_rows
            else [
                "filename",
                "image_type",
                "species",
                "slide_id",
                "specimen_id",
                "status",
                "source_image_path",
                "measurement_csv_path",
                "tile_manifest_path",
                "mask_root",
                "overlay_path",
            ],
        )
        writer.writeheader()
        writer.writerows(index_rows)


def build_run_manifest(args: argparse.Namespace, rows: list[dict[str, str]], git_hash: str) -> dict[str, Any]:
    return {
        "timestamp": time.strftime("%Y-%m-%dT%H:%M:%S"),
        "git_hash": git_hash,
        "manifest_path": str(args.manifest.resolve()),
        "output_dir": str(args.output_dir),
        "image_count": len(rows),
        "backend": args.backend,
        "gpu": args.gpu,
        "tile_filter": args.tile_filter,
        "max_cells_image": args.max_cells_image,
        "max_cells_species": args.max_cells_species,
        "diameter_override": args.diameter,
        "cellpose_model": str(args.cellpose_model.resolve()),
        "crops_enabled": not args.no_crops,
        "optimization_mode": "prefetch_gpu_finalize_pipeline",
        "gpu_tile_batch_images": args.gpu_tile_batch_images,
        "gpu_patch_batch_size": args.gpu_patch_batch_size,
        "cpu_tile_batch_images": args.cpu_tile_batch_images,
        "cpu_patch_batch_size": args.cpu_patch_batch_size,
        "traceability_note": "Optimized publication runner preserves tile manifests, masks, and measurement CSVs while allowing crop suppression for throughput.",
        "dry_run": args.dry_run,
    }


def _load_batch_tiles(segment_cells: Any, z: Any, batch_entries: list[tuple[float, int, int, int, int, str]]) -> tuple[list[Any], list[tuple[float, int, int, int, int, str]]]:
    tiles: list[Any] = []
    tile_infos: list[tuple[float, int, int, int, int, str]] = []
    for score, y0, x0, y1, x1, tile_name in batch_entries:
        tile = segment_cells.np.array(z[y0:y1, x0:x1])
        if tile.dtype == segment_cells.np.uint16:
            tile = (tile >> 8).astype(segment_cells.np.uint8)
        tiles.append(tile)
        tile_infos.append((score, y0, x0, y1, x1, tile_name))
    return tiles, tile_infos


def _finalize_batch(
    segment_cells: Any,
    *,
    tiles: list[Any],
    raw_label_images: list[Any],
    tile_infos: list[tuple[float, int, int, int, int, str]],
    artifacts: dict[str, Path],
    crop_dir: Path | None,
    crop_size: int,
    image_stem: str,
    crop_counter_start: int,
    batch_dt: float,
) -> dict[str, Any]:
    batch_measurements: list[dict[str, Any]] = []
    batch_manifest_rows: list[dict[str, Any]] = []
    processed_manifest_rows: dict[str, dict[str, Any]] = {}
    tile_summaries: list[dict[str, Any]] = []
    crop_counter = crop_counter_start

    for tile, raw_labeled, (score, y0, x0, y1, x1, tile_name) in zip(tiles, raw_label_images, tile_infos):
        tile_meas, _, raw_labeled, filtered_labeled = segment_cells.process_labeled_tile(tile, raw_labeled, y0, x0)

        raw_mask_path = artifacts["mask_dir"] / f"{tile_name[:-5]}__raw_labels.tiff"
        mask_path = artifacts["mask_dir"] / f"{tile_name[:-5]}__filtered_labels.tiff"
        segment_cells._write_label_mask(raw_labeled, raw_mask_path)
        segment_cells._write_label_mask(filtered_labeled, mask_path)

        manifest_row = {
            "tile_name": tile_name,
            "tile_y0": y0,
            "tile_x0": x0,
            "tile_h": y1 - y0,
            "tile_w": x1 - x0,
            "tile_score": score,
            "tile_status": "processed",
            "mask_path": str(mask_path.resolve()),
            "raw_mask_path": str(raw_mask_path.resolve()),
        }
        batch_manifest_rows.append(manifest_row)
        processed_manifest_rows[tile_name] = manifest_row

        for measurement in tile_meas:
            measurement["tile_score"] = score
            measurement["mask_path"] = str(mask_path.resolve())
            measurement["raw_mask_path"] = str(raw_mask_path.resolve())
            measurement["mask_label_id"] = measurement["label"]
            measurement["tile_manifest_path"] = str(artifacts["tile_manifest"].resolve())
            measurement["overlay_path"] = ""
            for key, default in CROP_FIELD_DEFAULTS.items():
                measurement.setdefault(key, default)

        if crop_dir and tile_meas:
            for measurement in tile_meas:
                crop_counter += 1
                local_y = measurement["centroid_y"] - y0
                local_x = measurement["centroid_x"] - x0
                filename = segment_cells.save_cell_crop(
                    tile,
                    local_y,
                    local_x,
                    crop_counter,
                    crop_dir,
                    image_stem,
                    crop_size,
                )
                measurement["crop_filename"] = str((crop_dir / filename).resolve())
                measurement["crop_size_px"] = crop_size
                measurement["centroid_in_crop_y"] = crop_size // 2
                measurement["centroid_in_crop_x"] = crop_size // 2

        batch_measurements.extend(tile_meas)
        tile_summaries.append(
            {
                "batch_dt": batch_dt,
                "batch_size": len(tile_infos),
                "n_cells": len(tile_meas),
                "score": score,
                "tile_name": tile_name,
                "x0": x0,
                "y0": y0,
            }
        )

    return {
        "batch_measurements": batch_measurements,
        "batch_manifest_rows": batch_manifest_rows,
        "crop_counter_end": crop_counter,
        "processed_manifest_rows": processed_manifest_rows,
        "tile_summaries": tile_summaries,
    }


def _add_budget_rows(
    budget_rows: dict[str, dict[str, Any]],
    batch_entries: list[tuple[float, int, int, int, int, str]],
) -> None:
    for score, y0, x0, y1, x1, tile_name in batch_entries:
        budget_rows[tile_name] = {
            "tile_name": tile_name,
            "tile_y0": y0,
            "tile_x0": x0,
            "tile_h": y1 - y0,
            "tile_w": x1 - x0,
            "tile_score": score,
            "tile_status": "not_processed_budget",
            "mask_path": "",
            "raw_mask_path": "",
        }


def optimized_process_tiled(
    segment_cells: Any,
    image_path: Path,
    name: str,
    img_h: int,
    img_w: int,
    output_dir: Path,
    *,
    crops: bool,
    crop_size: int,
    selected_tiles: list[tuple[int, int, int, int]] | None,
    tile_filter: str,
    max_cells: int,
) -> list[dict[str, Any]]:
    t0 = time.time()
    store = segment_cells.tifffile.imread(str(image_path), aszarr=True)
    z = segment_cells.zarr.open(store, mode="r")
    if z.ndim == 3:
        z = z[0]

    crop_dir = output_dir / "crops" / name if crops else None
    artifacts = segment_cells._artifact_paths(output_dir, name)
    partial_measurement_path = artifacts["partial_measurements"]
    existing_measurements = segment_cells._load_csv_rows(partial_measurement_path)
    existing_manifest_rows = segment_cells._load_tile_manifest_rows(artifacts["tile_manifest"])
    processed_manifest_rows = {
        row["tile_name"]: row
        for row in existing_manifest_rows
        if row.get("tile_status") == "processed" and row.get("tile_name")
    }
    processed_tile_names = set(processed_manifest_rows)
    if processed_tile_names or existing_measurements:
        print(
            f"  {name}: resuming {len(processed_tile_names)} processed tiles "
            f"and {len(existing_measurements)} saved measurements"
        )
    skipped_manifest_rows: list[dict[str, Any]] = []

    tile_coords: list[tuple[int, int, int, int]] = []
    for y0 in range(0, img_h, segment_cells.TILE_SIZE):
        for x0 in range(0, img_w, segment_cells.TILE_SIZE):
            y1 = min(y0 + segment_cells.TILE_SIZE, img_h)
            x1 = min(x0 + segment_cells.TILE_SIZE, img_w)
            if y1 - y0 >= segment_cells.EDGE_MARGIN * 2 and x1 - x0 >= segment_cells.EDGE_MARGIN * 2:
                tile_coords.append((y0, x0, y1, x1))

    total_tiles = len(tile_coords)
    if selected_tiles is not None:
        selected = {(t[0], t[1], t[2], t[3]) for t in selected_tiles}
        tile_coords = [tile for tile in tile_coords if tile in selected]
        print(f"  {name}: {len(tile_coords)}/{total_tiles} tiles selected")

    scored: list[tuple[float, int, int, int, int, str]] = []
    tiles_skipped = 0
    tiles_skipped_yellow = 0

    for y0, x0, y1, x1 in tile_coords:
        cy, cx = (y0 + y1) // 2, (x0 + x1) // 2
        half = min(128, (y1 - y0) // 2, (x1 - x0) // 2)
        probe = segment_cells.np.array(z[cy - half : cy + half, cx - half : cx + half])
        if probe.dtype == segment_cells.np.uint16:
            probe = (probe >> 8).astype(segment_cells.np.uint8)
        tmean, tstd = float(probe.mean()), float(probe.std())

        is_red = tmean > 230 or tmean < 10 or tstd < 5
        is_yellow = not is_red and (tstd < 15 or tmean > 200)
        tile_name = f"tile_y{y0:06d}_x{x0:06d}.tiff"

        if is_red:
            tiles_skipped += 1
            skipped_manifest_rows.append(
                {
                    "tile_name": tile_name,
                    "tile_y0": y0,
                    "tile_x0": x0,
                    "tile_h": y1 - y0,
                    "tile_w": x1 - x0,
                    "tile_score": "",
                    "tile_status": "skipped_empty",
                    "mask_path": "",
                    "raw_mask_path": "",
                }
            )
            continue
        if tile_filter == "green" and is_yellow:
            tiles_skipped_yellow += 1
            skipped_manifest_rows.append(
                {
                    "tile_name": tile_name,
                    "tile_y0": y0,
                    "tile_x0": x0,
                    "tile_h": y1 - y0,
                    "tile_w": x1 - x0,
                    "tile_score": "",
                    "tile_status": "skipped_borderline",
                    "mask_path": "",
                    "raw_mask_path": "",
                }
            )
            continue
        if tile_filter == "yellow" and not is_yellow:
            skipped_manifest_rows.append(
                {
                    "tile_name": tile_name,
                    "tile_y0": y0,
                    "tile_x0": x0,
                    "tile_h": y1 - y0,
                    "tile_w": x1 - x0,
                    "tile_score": "",
                    "tile_status": "skipped_filter",
                    "mask_path": "",
                    "raw_mask_path": "",
                }
            )
            continue

        score = float(segment_cells.np.exp(-((tstd - 18) ** 2) / (2 * 10**2)))
        if is_yellow:
            score *= 0.5
        scored.append((score, y0, x0, y1, x1, tile_name))

    scored.sort(reverse=True)
    skip_detail = f"{tiles_skipped} empty"
    if tiles_skipped_yellow:
        skip_detail += f", {tiles_skipped_yellow} borderline"
    print(f"  {name}: {img_w}x{img_h} ({len(scored)} tiles to process, {skip_detail}, filter={tile_filter})")

    measurements = list(existing_measurements)
    tiles_processed = len(processed_tile_names)
    crop_counter = len(measurements) if crop_dir and measurements else 0
    pending_scored = [entry for entry in scored if entry[5] not in processed_tile_names]
    budget_rows: dict[str, dict[str, Any]] = {}
    cap_reached = False

    batch_entries_list = [
        pending_scored[start : start + segment_cells._tile_batch_images()]
        for start in range(0, len(pending_scored), segment_cells._tile_batch_images())
    ]

    try:
        with ThreadPoolExecutor(max_workers=1) as prefetch_pool, ThreadPoolExecutor(max_workers=1) as finalize_pool:
            prefetch_future = None
            if batch_entries_list:
                prefetch_future = prefetch_pool.submit(_load_batch_tiles, segment_cells, z, batch_entries_list[0])
            pending_finalize = None

            for batch_index, batch_entries in enumerate(batch_entries_list):
                assert prefetch_future is not None
                tiles, tile_infos = prefetch_future.result()
                if batch_index + 1 < len(batch_entries_list):
                    prefetch_future = prefetch_pool.submit(
                        _load_batch_tiles, segment_cells, z, batch_entries_list[batch_index + 1]
                    )
                else:
                    prefetch_future = None

                t_batch = time.time()
                raw_label_images = segment_cells.segment_cells_batch(tiles)
                batch_dt = time.time() - t_batch

                if pending_finalize is not None:
                    result = pending_finalize.result()
                    segment_cells._append_rows(partial_measurement_path, result["batch_measurements"])
                    segment_cells._append_tile_manifest_rows(result["batch_manifest_rows"], artifacts["tile_manifest"])
                    start_tile_index = tiles_processed
                    measurements.extend(result["batch_measurements"])
                    tiles_processed += len(result["tile_summaries"])
                    crop_counter = result["crop_counter_end"]
                    processed_manifest_rows.update(result["processed_manifest_rows"])
                    for offset, summary in enumerate(result["tile_summaries"], start=1):
                        tile_index = start_tile_index + offset
                        if tile_index <= 3 or tile_index % 20 == 0:
                            print(
                                f"    tile {tile_index}/{len(scored)}: ({summary['y0']},{summary['x0']}) "
                                f"score={summary['score']:.2f} {summary['n_cells']} cells "
                                f"[batch {summary['batch_dt']:.1f}s/{summary['batch_size']} tiles] "
                                f"total={len(measurements)}",
                                flush=True,
                            )
                    if max_cells > 0 and len(measurements) >= max_cells:
                        print(f"    Reached {max_cells} cell cap after {tiles_processed}/{len(scored)} tiles")
                        _add_budget_rows(budget_rows, batch_entries)
                        for later_entries in batch_entries_list[batch_index + 1 :]:
                            _add_budget_rows(budget_rows, later_entries)
                        cap_reached = True
                        break

                pending_finalize = finalize_pool.submit(
                    _finalize_batch,
                    segment_cells,
                    tiles=tiles,
                    raw_label_images=raw_label_images,
                    tile_infos=tile_infos,
                    artifacts=artifacts,
                    crop_dir=crop_dir,
                    crop_size=crop_size,
                    image_stem=name,
                    crop_counter_start=crop_counter,
                    batch_dt=batch_dt,
                )

            if not cap_reached and pending_finalize is not None:
                result = pending_finalize.result()
                segment_cells._append_rows(partial_measurement_path, result["batch_measurements"])
                segment_cells._append_tile_manifest_rows(result["batch_manifest_rows"], artifacts["tile_manifest"])
                start_tile_index = tiles_processed
                measurements.extend(result["batch_measurements"])
                tiles_processed += len(result["tile_summaries"])
                crop_counter = result["crop_counter_end"]
                processed_manifest_rows.update(result["processed_manifest_rows"])
                for offset, summary in enumerate(result["tile_summaries"], start=1):
                    tile_index = start_tile_index + offset
                    if tile_index <= 3 or tile_index % 20 == 0:
                        print(
                            f"    tile {tile_index}/{len(scored)}: ({summary['y0']},{summary['x0']}) "
                            f"score={summary['score']:.2f} {summary['n_cells']} cells "
                            f"[batch {summary['batch_dt']:.1f}s/{summary['batch_size']} tiles] "
                            f"total={len(measurements)}",
                            flush=True,
                        )
                if max_cells > 0 and len(measurements) >= max_cells:
                    print(f"    Reached {max_cells} cell cap after {tiles_processed}/{len(scored)} tiles")
    finally:
        if hasattr(store, "close"):
            store.close()

    elapsed = time.time() - t0
    print(f"  {name}: {tiles_processed}/{len(scored)} tiles processed ({skip_detail}) -> {len(measurements)} isolated [{elapsed:.0f}s]")
    if crop_dir:
        print(f"  {name}: saved {len(measurements)} crops to {crop_dir}")

    final_manifest_rows = list(skipped_manifest_rows)
    for _, _, _, _, _, tile_name in scored:
        if tile_name in processed_manifest_rows:
            final_manifest_rows.append(processed_manifest_rows[tile_name])
        elif tile_name in budget_rows:
            final_manifest_rows.append(budget_rows[tile_name])
    segment_cells._write_tile_manifest(final_manifest_rows, artifacts["tile_manifest"])
    segment_cells._save_csv(measurements, output_dir, name)
    if partial_measurement_path.exists():
        partial_measurement_path.unlink()
    return measurements


def main() -> None:
    args = parse_args()
    cellprofiler_root = args.cellprofiler_root.resolve()
    pipeline_runs, convergence, segment_cells = setup_runtime(cellprofiler_root)

    segment_cells.GPU_TILE_BATCH_IMAGES = args.gpu_tile_batch_images
    segment_cells.GPU_PATCH_BATCH_SIZE = args.gpu_patch_batch_size
    segment_cells.CPU_TILE_BATCH_IMAGES = args.cpu_tile_batch_images
    segment_cells.CPU_PATCH_BATCH_SIZE = args.cpu_patch_batch_size

    rows = pipeline_runs.load_manifest_rows(args.manifest)
    if args.image_type != "all":
        rows = [row for row in rows if row["image_type"] == args.image_type]
        print(f"Cell-size: filtered to {len(rows)} {args.image_type} images")
    if args.max_cells_species > 0:
        rows = convergence.interleave_by_species(rows)
        print(f"Cell-size: interleaved images by species with cap={args.max_cells_species}")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    combined_csv = args.output_dir / "all_measurements.csv"
    run_manifest_path = args.output_dir / "run_manifest.json"
    index_path = args.output_dir / "image_index.csv"
    index_rows: list[dict[str, Any]] = []
    species_counts: dict[str, int] = {}
    run_manifest_path.write_text(
        json.dumps(build_run_manifest(args, rows, pipeline_runs.current_git_hash()), indent=2)
    )

    if args.resume and args.max_cells_species > 0 and combined_csv.exists():
        with combined_csv.open(newline="") as handle:
            for row in csv.DictReader(handle):
                species = str(row.get("species", "")).strip()
                if species:
                    species_counts[species] = species_counts.get(species, 0) + 1

    if not args.dry_run:
        segment_cells.get_model(
            gpu=args.gpu,
            diameter=args.diameter,
            pretrained_model=str(args.cellpose_model.resolve()),
        )

        original_process_tiled = segment_cells._process_tiled

        def patched_process_tiled(
            image_path: Path,
            name: str,
            img_h: int,
            img_w: int,
            output_dir: Path,
            crops: bool = True,
            crop_size: int = segment_cells.DEFAULT_CROP_SIZE,
            selected_tiles: list[tuple[int, int, int, int]] | None = None,
            tile_filter: str = "auto",
            max_cells: int = 0,
        ) -> list[dict[str, Any]]:
            return optimized_process_tiled(
                segment_cells,
                image_path=image_path,
                name=name,
                img_h=img_h,
                img_w=img_w,
                output_dir=output_dir,
                crops=crops,
                crop_size=crop_size,
                selected_tiles=selected_tiles,
                tile_filter=tile_filter,
                max_cells=max_cells,
            )

        segment_cells._process_tiled = patched_process_tiled
    else:
        original_process_tiled = None

    try:
        for row in rows:
            image_path = pipeline_runs.resolve_image_path(row)
            image_output_dir = args.output_dir / row["image_type"]
            stem = Path(row["filename"]).stem
            measurement_csv = image_output_dir / f"{stem}_measurements.csv"
            overlay_path = image_output_dir / "debug" / f"{stem}_overlay.png"
            tile_manifest_path = image_output_dir / "tile_manifests" / f"{stem}_tile_manifest.csv"
            mask_root = image_output_dir / "masks" / stem
            status = "planned"
            species = str(row.get("species", "")).strip()

            if args.max_cells_species > 0 and species and species_counts.get(species, 0) >= args.max_cells_species:
                status = "skipped_species_cap"
                index_rows.append(
                    {
                        "filename": row["filename"],
                        "image_type": row["image_type"],
                        "species": row["species"],
                        "slide_id": row["slide_id"],
                        "specimen_id": row["specimen_id"],
                        "status": status,
                        "source_image_path": str(image_path),
                        "measurement_csv_path": str(measurement_csv),
                        "tile_manifest_path": str(tile_manifest_path),
                        "mask_root": str(mask_root),
                        "overlay_path": str(overlay_path),
                    }
                )
                write_index(index_path, index_rows)
                continue

            if args.resume and measurement_csv.exists():
                status = "skipped_existing"
                index_rows.append(
                    {
                        "filename": row["filename"],
                        "image_type": row["image_type"],
                        "species": row["species"],
                        "slide_id": row["slide_id"],
                        "specimen_id": row["specimen_id"],
                        "status": status,
                        "source_image_path": str(image_path),
                        "measurement_csv_path": str(measurement_csv),
                        "tile_manifest_path": str(tile_manifest_path),
                        "mask_root": str(mask_root),
                        "overlay_path": str(overlay_path),
                    }
                )
                write_index(index_path, index_rows)
                continue

            if not args.dry_run:
                try:
                    measurements = segment_cells.process_image(
                        image_path,
                        image_output_dir,
                        debug=args.debug,
                        crops=not args.no_crops,
                        tile_filter=args.tile_filter,
                        image_type=row["image_type"],
                        max_cells=args.max_cells_image,
                        diameter_override=args.diameter,
                        pretrained_model=str(args.cellpose_model.resolve()),
                    )
                    for measurement in measurements:
                        measurement["filename"] = row["filename"]
                        measurement["slide_id"] = row["slide_id"]
                        measurement["specimen_id"] = row["specimen_id"]
                        measurement["species"] = row["species"]
                        measurement["image_type"] = row["image_type"]
                        measurement["source_image_path"] = str(image_path.resolve())
                        measurement["run_manifest_path"] = str(run_manifest_path.resolve())
                    append_rows(combined_csv, measurements)
                    if args.max_cells_species > 0 and species:
                        species_counts[species] = species_counts.get(species, 0) + len(measurements)
                    status = "completed"
                except Exception as exc:
                    print(f"  ERROR {row['filename']}: {exc}")
                    status = "error"

            index_rows.append(
                {
                    "filename": row["filename"],
                    "image_type": row["image_type"],
                    "species": row["species"],
                    "slide_id": row["slide_id"],
                    "specimen_id": row["specimen_id"],
                    "status": status,
                    "source_image_path": str(image_path),
                    "measurement_csv_path": str(measurement_csv),
                    "tile_manifest_path": str(tile_manifest_path),
                    "mask_root": str(mask_root),
                    "overlay_path": str(overlay_path),
                }
            )
            write_index(index_path, index_rows)
    finally:
        if not args.dry_run and original_process_tiled is not None:
            segment_cells._process_tiled = original_process_tiled

    print(f"Cell-size image index: {index_path}")
    print(f"Cell-size manifest: {run_manifest_path}")


if __name__ == "__main__":
    main()
