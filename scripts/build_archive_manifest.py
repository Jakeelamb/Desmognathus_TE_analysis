#!/usr/bin/env python3
"""Build a deterministic checksum manifest for a pre-cleanup archive.

The manifest intentionally excludes Git internals, its own output files, and
unambiguous transient tool caches. Every other regular file is hashed in
relative-path order and assigned a broad retention class.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import stat
import tempfile
import time
from collections import defaultdict
from collections.abc import Iterable
from concurrent.futures import ThreadPoolExecutor
from datetime import UTC, datetime
from pathlib import Path, PurePosixPath
from typing import NamedTuple

MANIFEST_NAME = "ARCHIVE_MANIFEST.csv"
SUMMARY_NAME = "ARCHIVE_MANIFEST_SUMMARY.json"
OUTPUT_NAMES = frozenset({MANIFEST_NAME, SUMMARY_NAME})

CACHE_DIRECTORY_NAMES = frozenset(
    {
        "__pycache__",
        ".pytest_cache",
        ".ruff_cache",
        ".mypy_cache",
        ".ipynb_checkpoints",
    }
)

HPC_INPUT_PREFIXES = (
    PurePosixPath("input_data/dnaPipeTE"),
    PurePosixPath("input_data/repeatmasker"),
    PurePosixPath("input_data/ectopic_recombination"),
)
GENOME_CACHE_PREFIX = PurePosixPath("input_data/genomes")
INTERMEDIATE_PREFIX = PurePosixPath("interim")
DERIVED_PREFIXES = (
    PurePosixPath("results"),
    PurePosixPath("Publication"),
    PurePosixPath("path_analysis/data/derived"),
)

RETENTION_DESCRIPTIONS = {
    "hpc_input_preserve": (
        "Known HPC-derived dnaPipeTE, RepeatMasker, or ectopic-recombination input; "
        "preserve as primary computational evidence."
    ),
    "remote_genome_cache": (
        "Downloaded genome assembly cache; preserve until accession and upstream "
        "checksum provenance are independently frozen."
    ),
    "intermediate_regenerable": (
        "Local intermediate; eligible for regeneration after its inputs and code are verified."
    ),
    "derived_result": (
        "Generated result, publication artifact, or derived analysis table retained "
        "for recovery and audit."
    ),
    "other_pre_cleanup_worktree": (
        "Other file present in the dirty pre-cleanup working tree; retain for "
        "recovery and reference."
    ),
}


class ManifestRow(NamedTuple):
    relative_path: str
    bytes: int
    sha256: str
    retention_class: str


def has_prefix(path: PurePosixPath, prefix: PurePosixPath) -> bool:
    """Return whether *path* is equal to or contained by *prefix*."""

    return path == prefix or prefix in path.parents


def retention_class(relative_path: str) -> str:
    path = PurePosixPath(relative_path)
    if any(has_prefix(path, prefix) for prefix in HPC_INPUT_PREFIXES):
        return "hpc_input_preserve"
    if has_prefix(path, GENOME_CACHE_PREFIX):
        return "remote_genome_cache"
    if has_prefix(path, INTERMEDIATE_PREFIX):
        return "intermediate_regenerable"
    if any(has_prefix(path, prefix) for prefix in DERIVED_PREFIXES):
        return "derived_result"
    return "other_pre_cleanup_worktree"


def iter_regular_files(root: Path) -> tuple[list[Path], dict[str, object]]:
    """Discover manifestable files without following directory symlinks."""

    files: list[Path] = []
    excluded_cache_files = 0
    excluded_cache_bytes = 0
    symlinks_not_followed: list[str] = []

    for current, directory_names, file_names in os.walk(root, topdown=True, followlinks=False):
        current_path = Path(current)

        kept_directories: list[str] = []
        for name in sorted(directory_names):
            candidate = current_path / name
            if name == ".git":
                continue
            if name in CACHE_DIRECTORY_NAMES:
                for cache_current, _, cache_files in os.walk(candidate, followlinks=False):
                    for cache_name in cache_files:
                        cache_path = Path(cache_current) / cache_name
                        try:
                            cache_stat = cache_path.stat(follow_symlinks=False)
                        except FileNotFoundError:
                            continue
                        if stat.S_ISREG(cache_stat.st_mode):
                            excluded_cache_files += 1
                            excluded_cache_bytes += cache_stat.st_size
                continue
            if candidate.is_symlink():
                symlinks_not_followed.append(candidate.relative_to(root).as_posix())
                continue
            kept_directories.append(name)
        directory_names[:] = kept_directories

        for name in sorted(file_names):
            candidate = current_path / name
            relative = candidate.relative_to(root).as_posix()
            if relative in OUTPUT_NAMES:
                continue
            if candidate.parent == root and name.startswith(".ARCHIVE_MANIFEST."):
                continue
            try:
                candidate_stat = candidate.stat(follow_symlinks=False)
            except FileNotFoundError:
                continue
            if stat.S_ISLNK(candidate_stat.st_mode):
                symlinks_not_followed.append(relative)
            elif stat.S_ISREG(candidate_stat.st_mode):
                files.append(candidate)

    files.sort(key=lambda path: path.relative_to(root).as_posix())
    exclusion_details: dict[str, object] = {
        "git_directories": "excluded wherever the directory name is .git",
        "manifest_outputs": sorted(OUTPUT_NAMES),
        "transient_cache_directories": sorted(CACHE_DIRECTORY_NAMES),
        "transient_cache_file_count": excluded_cache_files,
        "transient_cache_bytes": excluded_cache_bytes,
        "symlinks_not_followed": sorted(symlinks_not_followed),
    }
    return files, exclusion_details


def sha256_file(path: Path, chunk_size: int = 8 * 1024 * 1024) -> tuple[int, str]:
    """Hash one stable regular file, failing if it changes during the read."""

    before = path.stat(follow_symlinks=False)
    digest = hashlib.sha256()
    with path.open("rb", buffering=0) as handle:
        while chunk := handle.read(chunk_size):
            digest.update(chunk)
    after = path.stat(follow_symlinks=False)

    stable_before = (before.st_dev, before.st_ino, before.st_size, before.st_mtime_ns)
    stable_after = (after.st_dev, after.st_ino, after.st_size, after.st_mtime_ns)
    if stable_before != stable_after:
        raise RuntimeError(f"file changed while hashing: {path}")
    return before.st_size, digest.hexdigest()


def build_rows(root: Path, files: Iterable[Path], workers: int) -> list[ManifestRow]:
    paths = list(files)
    with ThreadPoolExecutor(max_workers=workers) as executor:
        hashes = executor.map(sha256_file, paths)
        rows = [
            ManifestRow(
                relative_path=path.relative_to(root).as_posix(),
                bytes=size,
                sha256=digest,
                retention_class=retention_class(path.relative_to(root).as_posix()),
            )
            for path, (size, digest) in zip(paths, hashes, strict=True)
        ]
    return rows


def new_temporary_file(root: Path) -> tuple[int, Path]:
    descriptor, raw_path = tempfile.mkstemp(dir=root, prefix=".ARCHIVE_MANIFEST.", suffix=".tmp")
    return descriptor, Path(raw_path)


def write_csv_temporary(root: Path, rows: list[ManifestRow]) -> tuple[Path, str]:
    descriptor, temporary_path = new_temporary_file(root)
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, lineterminator="\n")
            writer.writerow(ManifestRow._fields)
            writer.writerows(rows)
            handle.flush()
            os.fsync(handle.fileno())
        _, digest = sha256_file(temporary_path)
        return temporary_path, digest
    except BaseException:
        temporary_path.unlink(missing_ok=True)
        raise


def write_json_temporary(root: Path, summary: dict[str, object]) -> Path:
    descriptor, temporary_path = new_temporary_file(root)
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as handle:
            json.dump(summary, handle, indent=2, sort_keys=True)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        return temporary_path
    except BaseException:
        temporary_path.unlink(missing_ok=True)
        raise


def fsync_directory(path: Path) -> None:
    descriptor = os.open(path, os.O_RDONLY | os.O_DIRECTORY)
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def write_outputs(
    root: Path,
    rows: list[ManifestRow],
    exclusions: dict[str, object],
    started_at: float,
) -> dict[str, object]:
    by_class: dict[str, dict[str, int]] = defaultdict(lambda: {"file_count": 0, "bytes": 0})
    for row in rows:
        by_class[row.retention_class]["file_count"] += 1
        by_class[row.retention_class]["bytes"] += row.bytes

    csv_temporary, manifest_sha256 = write_csv_temporary(root, rows)
    duration_seconds = time.monotonic() - started_at
    summary: dict[str, object] = {
        "archive_root": str(root),
        "duration_seconds": round(duration_seconds, 3),
        "exclusions": exclusions,
        "file_count": len(rows),
        "generated_at_utc": datetime.now(UTC).isoformat(),
        "hash_algorithm": "sha256",
        "manifest_file": MANIFEST_NAME,
        "manifest_sha256": manifest_sha256,
        "retention_class_descriptions": RETENTION_DESCRIPTIONS,
        "retention_classes": {
            name: by_class.get(name, {"file_count": 0, "bytes": 0})
            for name in RETENTION_DESCRIPTIONS
        },
        "total_bytes": sum(row.bytes for row in rows),
    }
    try:
        json_temporary = write_json_temporary(root, summary)
    except BaseException:
        csv_temporary.unlink(missing_ok=True)
        raise

    try:
        os.replace(csv_temporary, root / MANIFEST_NAME)
        os.replace(json_temporary, root / SUMMARY_NAME)
        fsync_directory(root)
    except BaseException:
        csv_temporary.unlink(missing_ok=True)
        json_temporary.unlink(missing_ok=True)
        raise
    return summary


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("archive_root", type=Path, help="archive directory to scan")
    parser.add_argument(
        "--workers",
        type=int,
        default=min(4, os.cpu_count() or 1),
        help="parallel hashing workers (default: up to 4)",
    )
    arguments = parser.parse_args()
    if arguments.workers < 1:
        parser.error("--workers must be at least 1")
    return arguments


def main() -> int:
    arguments = parse_arguments()
    root = arguments.archive_root.expanduser().resolve()
    if not root.is_dir():
        raise SystemExit(f"archive root is not a directory: {root}")

    started_at = time.monotonic()
    files, exclusions = iter_regular_files(root)
    rows = build_rows(root, files, arguments.workers)
    summary = write_outputs(root, rows, exclusions, started_at)
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
