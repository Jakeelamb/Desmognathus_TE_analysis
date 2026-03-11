#!/usr/bin/env python3
"""Pull available CellProfiler outputs into path_analysis and audit traceability."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from cellprofiler_bridge.bridge import rebuild_imported_artifacts


PROJECT_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_CELLPROFILER_ROOT = Path.home() / "Projects" / "cellprofiler_test"
DEFAULT_OUT_DIR = PROJECT_ROOT / "path_analysis" / "data" / "external" / "derived"
DEFAULT_RAW_RUN_TAG = "full_dataset_v1"
DEFAULT_MIXED_RUN_TAG = "mixed_cellpose_yolo_full_dataset_v1"
def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--cellprofiler-root",
        type=Path,
        default=DEFAULT_CELLPROFILER_ROOT,
        help="Root of the cellprofiler_test repository",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=DEFAULT_OUT_DIR,
        help="Output directory for imported CellProfiler artifacts",
    )
    parser.add_argument(
        "--raw-run-tag",
        default=DEFAULT_RAW_RUN_TAG,
        help="Raw nucleus-IOD / cell-size run tag to import",
    )
    parser.add_argument(
        "--mixed-run-tag",
        default=DEFAULT_MIXED_RUN_TAG,
        help="Mixed linkage run tag to import",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    payload = rebuild_imported_artifacts(
        cellprofiler_root=args.cellprofiler_root,
        out_dir=args.out_dir,
        raw_run_tag=args.raw_run_tag,
        mixed_run_tag=args.mixed_run_tag,
    )
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
