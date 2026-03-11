from __future__ import annotations

import csv
import subprocess
import time
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

from .common import write_json


REQUIRED_MANIFEST_COLUMNS = {
    "filename",
    "image_type",
    "species",
    "slide_id",
    "specimen_id",
}


@dataclass(frozen=True)
class CanonicalRunPaths:
    run_tag: str
    root: Path
    cell_size_segmentation: Path
    nucleus_iod: Path
    qc: Path
    traceability: Path
    logs: Path
    manifests: Path

    def as_dict(self) -> dict[str, str]:
        return {key: str(value) for key, value in asdict(self).items()}


def data_root(cellprofiler_root: Path) -> Path:
    root = Path(cellprofiler_root)
    candidate = root / "data"
    if candidate.exists():
        return candidate
    raise FileNotFoundError(f"Could not locate data root under {root}")


def output_root(cellprofiler_root: Path) -> Path:
    return Path(cellprofiler_root) / "output"


def runs_root(cellprofiler_root: Path) -> Path:
    return output_root(cellprofiler_root) / "runs"


def metadata_path(cellprofiler_root: Path) -> Path:
    root = Path(cellprofiler_root)
    candidates = [
        root / "data" / "metadata" / "master_image_metadata.csv",
        root / "master_image_metadata.csv",
    ]
    for candidate in candidates:
        if candidate.exists():
            return candidate
    raise FileNotFoundError(f"Could not locate master_image_metadata.csv under {root}")


def get_run_paths(cellprofiler_root: Path, run_tag: str) -> CanonicalRunPaths:
    root = runs_root(cellprofiler_root) / run_tag
    return CanonicalRunPaths(
        run_tag=run_tag,
        root=root,
        cell_size_segmentation=root / "cell_size_segmentation",
        nucleus_iod=root / "nucleus_iod",
        qc=root / "qc",
        traceability=root / "traceability",
        logs=root / "logs",
        manifests=root / "manifests",
    )


def ensure_run_directories(paths: CanonicalRunPaths) -> list[Path]:
    directories = [
        paths.root,
        paths.cell_size_segmentation,
        paths.nucleus_iod,
        paths.qc,
        paths.traceability,
        paths.logs,
        paths.manifests,
    ]
    for directory in directories:
        directory.mkdir(parents=True, exist_ok=True)
    return directories


def load_metadata_lookup(cellprofiler_root: Path) -> dict[tuple[str, str], dict[str, str]]:
    with metadata_path(cellprofiler_root).open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    return {
        (row["image_type"], row["filename"]): row
        for row in rows
        if row.get("filename") and row.get("image_type")
    }


def load_manifest_rows(cellprofiler_root: Path, manifest_path: Path) -> list[dict[str, str]]:
    manifest_path = Path(manifest_path)
    if not manifest_path.is_absolute():
        manifest_path = Path(cellprofiler_root) / manifest_path
    with manifest_path.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        raise ValueError(f"Manifest is empty: {manifest_path}")

    lookup = load_metadata_lookup(cellprofiler_root)
    enriched = []
    for row in rows:
        key = (row.get("image_type", ""), row.get("filename", ""))
        meta = lookup.get(key, {})
        merged = dict(meta)
        merged.update({key: value for key, value in row.items() if value not in (None, "")})
        enriched.append(merged)

    missing = {required for required in REQUIRED_MANIFEST_COLUMNS if any(not row.get(required) for row in enriched)}
    if missing:
        cols = ", ".join(sorted(missing))
        raise ValueError(f"Manifest rows still missing required fields after metadata join: {cols}")
    return enriched


def resolve_image_path(cellprofiler_root: Path, row: dict[str, str]) -> Path:
    root = data_root(cellprofiler_root)
    candidate = root / row["image_type"] / row["filename"]
    if candidate.exists():
        return candidate
    return candidate


def manifest_summary(rows: list[dict[str, str]]) -> dict[str, Any]:
    def counts(key: str) -> dict[str, int]:
        values: dict[str, int] = {}
        for row in rows:
            value = row.get(key, "")
            values[value] = values.get(value, 0) + 1
        return dict(sorted(values.items()))

    return {
        "image_count": len(rows),
        "image_types": counts("image_type"),
        "preservation": counts("preservation"),
        "mounting": counts("mounting"),
        "species": counts("species"),
    }


def rows_missing_images(cellprofiler_root: Path, rows: list[dict[str, str]]) -> list[dict[str, str]]:
    missing = []
    for row in rows:
        image_path = resolve_image_path(cellprofiler_root, row)
        if not image_path.exists():
            missing.append(
                {
                    "filename": row["filename"],
                    "image_type": row["image_type"],
                    "species": row.get("species", ""),
                    "expected_path": str(image_path),
                }
            )
    return missing


def current_git_hash(cellprofiler_root: Path) -> str:
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "--short", "HEAD"],
            cwd=Path(cellprofiler_root),
            text=True,
            timeout=5,
        ).strip()
    except Exception:
        return ""


def copy_manifest_into_run(manifest_path: Path, paths: CanonicalRunPaths) -> Path:
    dest = paths.manifests / Path(manifest_path).name
    dest.write_text(Path(manifest_path).read_text())
    return dest


def build_run_config(
    *,
    cellprofiler_root: Path,
    run_tag: str,
    manifest_path: Path,
    rows: list[dict[str, str]],
    paths: CanonicalRunPaths,
    extra: dict[str, Any] | None = None,
) -> dict[str, Any]:
    root = Path(cellprofiler_root).resolve()
    config = {
        "run_tag": run_tag,
        "created_at": time.strftime("%Y-%m-%dT%H:%M:%S"),
        "git_hash": current_git_hash(root),
        "repo_root": str(root),
        "manifest_path": str(Path(manifest_path).resolve()),
        "run_paths": paths.as_dict(),
        "manifest_summary": manifest_summary(rows),
        "shared_inputs": {
            "data_root": str(data_root(root)),
            "metadata_path": str(metadata_path(root)),
        },
        "pipeline_entrypoints": {
            "cell_size": str(root / "cell_size_segmentation_pipeline" / "run_from_manifest.py"),
            "nucleus_iod": str(root / "nucleus_iod_estimate_pipeline" / "run_from_manifest.py"),
        },
        "traceability_contract": [
            "result_row",
            "object_row",
            "mask_with_coordinates",
            "tile_or_block_coordinates",
            "source_image_file",
        ],
    }
    if extra:
        config.update(extra)
    return config


def write_run_readme(paths: CanonicalRunPaths, config: dict[str, Any]) -> Path:
    summary = config["manifest_summary"]
    lines = [
        f"# Canonical Run: {config['run_tag']}",
        "",
        "This run root is the shared execution/output surface for the active production pipelines.",
        "",
        "## Manifest",
        "",
        f"- Source manifest: `{config['manifest_path']}`",
        f"- Copied manifest: `{paths.manifests / Path(config['manifest_path']).name}`",
        "",
        "## Summary",
        "",
        f"- Images: {summary['image_count']}",
        f"- Image types: {summary['image_types']}",
        f"- Preservation: {summary['preservation']}",
        "",
        "## Pipeline Roots",
        "",
        f"- Cell size: `{paths.cell_size_segmentation}`",
        f"- Nucleus IOD: `{paths.nucleus_iod}`",
        f"- QC: `{paths.qc}`",
        f"- Traceability: `{paths.traceability}`",
        f"- Logs: `{paths.logs}`",
    ]
    readme = paths.root / "README.md"
    readme.write_text("\n".join(lines) + "\n")
    return readme


def write_run_config(path: Path, config: dict[str, Any]) -> None:
    write_json(path, config)
