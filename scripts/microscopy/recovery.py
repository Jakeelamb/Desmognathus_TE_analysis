#!/usr/bin/env python3
"""Immutable, metadata-blinded review packets for existing nuclear masks."""

from __future__ import annotations

import argparse
import csv
import hashlib
import importlib
import itertools
import json
import math
import os
import re
import shlex
import shutil
import sys
import zipfile
from collections import Counter, defaultdict
from dataclasses import asdict, dataclass
from datetime import UTC, datetime
from pathlib import Path

import numpy as np
import tifffile
from PIL import Image

sys.dont_write_bytecode = True
SCHEMA = "microscopy-recovery-v1"
PILOT = "standard47_fuscus35"
STANDARD_REL = (
    "output/runs/life_history_pmount_reference_v1/nucleus_measurements/all_measurements.csv"
)
TARGET_REL = "analyses/02_morphology/data/nuclear_iod_objects.csv"
STANDARD_SHA = "33f8a90344e362e6d61ca242703a501d27bd1a0b786f86b5035e4dfc3b29bd22"
TARGET_SHA = "fffc91fa5f53a374711b7b71baf2292733e6185fa7d86bb02dc189a9cfc6d4bf"
STANDARD_PINS_SHA = "d3bfed24e791e61836818572879177ed6b93c203f4b11b344d0881da6c305bfc"
STANDARD_PINS_PATH = Path(__file__).with_name("standard_qc_source_pins.json")
UI_PATH = Path(__file__).with_name("recovery_review.html")
PROTOCOL_TEXT = (
    "Inspect the complete isolated nuclear boundary, focus, overlap and debris. "
    "Every decision needs your name and explicit reasons. Blank cases remain unreviewed. "
    "Exports are authored events; this page's local cache is a draft. "
    "Metadata blinding does not prevent recognizing morphology.\n"
)
CANONICAL_CELL = Path("/home/jake/Projects/cellprofiler_test")
COUNTS = {"413": 31, "414": 16, "314": 12, "315": 8, "320": 5, "321": 10}
SPECIMENS = {
    "413": "32469",
    "414": "32470",
    "314": "3372",
    "315": "3373",
    "320": "3401",
    "321": "3402",
}
BACKGROUNDS = {"413": 255, "414": 255, "314": 255, "315": 255, "320": 218, "321": 216}
CASE_KEYS = {"case_id", "raw", "mask", "width", "height", "state_sha256"}
EVENT_KEYS = {"case_id", "state_sha256", "reviewer", "time", "decision", "reasons", "note"}
SHA_RE = re.compile(r"[0-9a-f]{64}\Z")


@dataclass(frozen=True)
class SourcePixel:
    path: str
    sha256: str
    axes: str
    dtype: str
    width_px: int
    height_px: int
    tile_global_x0: int
    tile_global_y0: int
    centroid_frame: str
    channel: str = "green"


@dataclass(frozen=True)
class Mask:
    path: str
    sha256: str
    label: int
    width_px: int
    height_px: int
    frame: str
    lineage: str


@dataclass(frozen=True)
class PilotCase:
    case_id: str
    object_id: str
    cohort: str
    process: str
    specimen: str
    source: SourcePixel
    baseline_mask: Mask
    pre_repair_mask: Mask | None
    stored_centroid_x: float
    stored_centroid_y: float
    stored_area_px: int
    stored_iod: float
    recorded_i_bg: float
    source_bbox: tuple[int, int, int, int]
    legacy: dict[str, str]


def clock() -> str:
    return datetime.now(UTC).isoformat()


def digest(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def encoded(value) -> bytes:
    return (json.dumps(value, sort_keys=True, indent=2, allow_nan=False) + "\n").encode()


def value_hash(value) -> str:
    return hashlib.sha256(encoded(value)).hexdigest()


def write_json(path: Path, value) -> None:
    with path.open("xb") as handle:
        handle.write(encoded(value))


def no_symlinks(path: Path) -> Path:
    # Check the supplied ancestry before normalization can erase a symlink/.. pair.
    supplied = path.absolute()
    for component in (supplied, *supplied.parents):
        if component.is_symlink():
            raise ValueError(f"Symlinked input/output or ancestor: {component}")
    result = Path(os.path.abspath(path))
    for component in (result, *result.parents):
        if component.is_symlink():
            raise ValueError(f"Symlinked input/output or ancestor: {component}")
    return result


def input_file(path: Path) -> Path:
    result = no_symlinks(path)
    if not result.is_file():
        raise ValueError(f"Input is not a regular file: {result}")
    return result


def fresh_destination(path: Path, sources: list[Path]) -> Path:
    result = no_symlinks(path)
    if result.exists():
        raise ValueError(f"Output already exists: {result}")
    if not result.parent.is_dir():
        raise ValueError(f"Output parent does not exist: {result.parent}")
    for source in sources:
        source = no_symlinks(source)
        if result == source or result in source.parents or source in result.parents:
            raise ValueError(f"Input/output ancestry overlap: {source}")
    return result


def verify_pins(pins: dict[str, str]) -> None:
    if not isinstance(pins, dict) or not pins:
        raise ValueError("Missing source pins")
    for name, expected in pins.items():
        if not isinstance(expected, str) or not SHA_RE.fullmatch(expected):
            raise ValueError("Invalid source hash")
        path = input_file(Path(name))
        if digest(path) != expected:
            raise ValueError(f"Source hash mismatch: {path}")


def load_json(path: Path):
    def unique(pairs):
        result = {}
        for key, value in pairs:
            if key in result:
                raise ValueError(f"Duplicate JSON key: {key}")
            result[key] = value
        return result

    return json.loads(
        input_file(path).read_text(),
        object_pairs_hook=unique,
        parse_constant=lambda value: (_ for _ in ()).throw(ValueError(f"Nonfinite JSON: {value}")),
    )


def number(value, name: str) -> float:
    try:
        result = float(value)
    except (TypeError, ValueError) as error:
        raise ValueError(f"Invalid {name}") from error
    if not math.isfinite(result):
        raise ValueError(f"Nonfinite {name}")
    return result


def integer(value, name: str, *, positive: bool = False) -> int:
    result = number(value, name)
    if result != math.floor(result) or result < (1 if positive else 0):
        raise ValueError(f"Invalid integer {name}")
    return int(result)


def read_array(path: Path, *, raw: bool = False, binary: bool = False) -> np.ndarray:
    input_file(path)
    if path.suffix.lower() in {".tif", ".tiff"}:
        with tifffile.TiffFile(path) as image:
            if len(image.series) != 1 or image.series[0].axes != "YX":
                raise ValueError(f"Image requires exactly one YX series: {path}")
            result = image.asarray()
    elif path.suffix.lower() == ".png":
        with Image.open(path) as image:
            if raw and image.mode != "L":
                raise ValueError(f"Raw image requires native grayscale uint8: {path}")
            result = np.array(image)
    else:
        raise ValueError(f"Unsupported lossless image: {path}")
    if result.ndim != 2 or min(result.shape) <= 0:
        raise ValueError(f"Image requires nonempty YX pixels: {path}")
    if raw and result.dtype != np.uint8:
        raise ValueError(f"Raw image requires native grayscale uint8: {path}")
    if not raw and result.dtype.kind not in "uib":
        raise ValueError(f"Label mask requires integer pixels: {path}")
    if not raw and np.any(result < 0):
        raise ValueError(f"Negative mask labels: {path}")
    if binary and (result.dtype != np.uint8 or not np.isin(result, [0, 255]).all()):
        raise ValueError(f"Binary mask requires uint8 0/255: {path}")
    return result


def save_png(path: Path, array: np.ndarray) -> None:
    Image.fromarray(np.asarray(array, dtype=np.uint8)).save(path, format="PNG")


def read_csv(path: Path) -> list[dict[str, str]]:
    with input_file(path).open(newline="") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames or len(reader.fieldnames) != len(set(reader.fieldnames)):
            raise ValueError("Missing or duplicate CSV header")
        return list(reader)


def remember_pin(path: Path, pins: dict[str, str]) -> str:
    name = str(path)
    if name not in pins:
        pins[name] = digest(path)
    return pins[name]


def mapped_path(value: str, cell_root: Path) -> Path:
    path = Path(value)
    if path.is_relative_to(CANONICAL_CELL):
        path = cell_root / path.relative_to(CANONICAL_CELL)
    if not path.is_absolute() or not path.is_relative_to(cell_root):
        raise ValueError(f"Source outside the selected CellProfiler root: {path}")
    return input_file(path)


def case_from_row(
    row: dict[str, str], standard: bool, cell_root: Path, seed: str, pins: dict[str, str]
) -> PilotCase:
    process_pattern = (
        r"Process_(413|414)_green\.ome\.tiff"
        if standard
        else r"Process_(314|315|320|321)_raw_green\.ome\.tiff"
    )
    match = re.fullmatch(process_pattern, row.get("filename", ""))
    if not match or row.get("species") != "D. fuscus":
        raise ValueError("Wrong pilot process/species")
    process = match[1]
    if row.get("specimen_id") != SPECIMENS[process] or row.get("image_type") != (
        "pmount" if standard else "brightfield"
    ):
        raise ValueError("Wrong pilot specimen/image identity")
    tile = row.get("tile_name", "")
    if not re.fullmatch(r"tile_y\d+_x\d+\.tiff", tile):
        raise ValueError("Wrong tile identity")
    label = integer(row.get("label" if standard else "nucleus_label"), "label", positive=True)
    identity = f"{row['filename']}::{tile}::{label}"
    if not standard and row.get("nucleus_object_id") != identity:
        raise ValueError("Nucleus identity differs from process/tile/label")
    source = mapped_path(
        row["source_image_path" if standard else "nucleus_source_image_path"], cell_root
    )
    mask = mapped_path(row["mask_path" if standard else "nucleus_mask_path"], cell_root)
    old = None if standard else mapped_path(row["nucleus_mask_path_original"], cell_root)
    for field in (
        ("tile_manifest_path", "run_manifest_path")
        if standard
        else ("nucleus_tile_manifest_path", "nucleus_run_manifest_path")
    ):
        metadata = mapped_path(row[field], cell_root)
        remember_pin(metadata, pins)
    for path in (source, mask, old):
        if path is not None:
            remember_pin(path, pins)
    width = integer(row["tile_width_px"], "tile width", positive=True)
    height = integer(row["tile_height_px"], "tile height", positive=True)
    x0 = integer(row["tile_x0"], "tile x origin")
    y0 = integer(row["tile_y0"], "tile y origin")
    background = number(row["i_bg" if standard else "nuc_i_bg"], "recorded background")
    if background != BACKGROUNDS[process]:
        raise ValueError("Wrong recorded background for pilot process")
    return PilotCase(
        "N" + hashlib.sha256((seed + "\0" + identity).encode()).hexdigest()[:20],
        identity,
        f"standard{process}" if standard else "fuscus_target",
        process,
        row["specimen_id"],
        SourcePixel(
            str(source),
            pins[str(source)],
            "YX",
            "uint8",
            width,
            height,
            x0,
            y0,
            "slide_global" if standard else "tile_local",
        ),
        Mask(
            str(mask),
            pins[str(mask)],
            label,
            width,
            height,
            "tile_local",
            "standard_original" if standard else "released_current",
        ),
        None
        if old is None
        else Mask(
            str(old), pins[str(old)], label, width, height, "tile_local", "historical_pre_repair"
        ),
        number(row["centroid_x" if standard else "nuc_centroid_x"], "centroid x"),
        number(row["centroid_y" if standard else "nuc_centroid_y"], "centroid y"),
        integer(row["area_px" if standard else "nuc_area_px"], "stored area", positive=True),
        number(row["iod" if standard else "nuc_iod"], "stored IOD"),
        background,
        (0, 0, 0, 0),
        {
            key: row.get(key, "")
            for key in (
                "image_type",
                "review_status",
                "review_decision",
                "review_decision_source",
                "nucleus_repair_applied",
                "has_cell_match",
                "keep_strict_core",
            )
        },
    )


def group_cases(cases: list[dict]):
    ordered = sorted(
        cases,
        key=lambda case: (case["source"]["path"], case["baseline_mask"]["path"], case["object_id"]),
    )
    return itertools.groupby(
        ordered, key=lambda case: (case["source"]["path"], case["baseline_mask"]["path"])
    )


def native_membership(
    case: dict, raw: np.ndarray, labels: np.ndarray
) -> tuple[np.ndarray, tuple[int, int, int, int]]:
    source = case["source"]
    if source["axes"] != "YX" or source["dtype"] != "uint8" or source["channel"] != "green":
        raise ValueError("Wrong axes/dtype/channel contract")
    if raw.shape != labels.shape or raw.shape != (source["height_px"], source["width_px"]):
        raise ValueError("Source/mask shape differs from declared tile dimensions")
    if (
        case["baseline_mask"]["width_px"] != source["width_px"]
        or case["baseline_mask"]["height_px"] != source["height_px"]
    ):
        raise ValueError("Mask dimensions differ from source")
    membership = labels == case["baseline_mask"]["label"]
    yy, xx = np.nonzero(membership)
    if len(xx) != case["stored_area_px"]:
        raise ValueError(f"Original area mismatch: {case['case_id']}")
    expected_x = case["stored_centroid_x"]
    expected_y = case["stored_centroid_y"]
    frame = source["centroid_frame"]
    if frame == "slide_global":
        expected_x -= source["tile_global_x0"]
        expected_y -= source["tile_global_y0"]
    elif frame != "tile_local":
        raise ValueError("Unknown centroid coordinate frame")
    if abs(float(xx.mean()) - expected_x) > 0.011 or abs(float(yy.mean()) - expected_y) > 0.011:
        raise ValueError("Centroid/origin mismatch")
    return membership, (int(xx.min()), int(yy.min()), int(xx.max()) + 1, int(yy.max()) + 1)


def validate_case_set(cases: list[dict]) -> None:
    if (
        not isinstance(cases, list)
        or len(cases) != 82
        or Counter(case["process"] for case in cases) != COUNTS
    ):
        raise ValueError("Pilot requires exactly the frozen 47 standards + 35 target identities")
    for field in ("case_id", "object_id"):
        if len({case[field] for case in cases}) != len(cases):
            raise ValueError(f"Duplicate {field}")
    for case in cases:
        process = case["process"]
        standard = process in {"413", "414"}
        if (
            case["cohort"] != (f"standard{process}" if standard else "fuscus_target")
            or case["specimen"] != SPECIMENS[process]
        ):
            raise ValueError("Wrong cohort/specimen contract")
        if case["recorded_i_bg"] != BACKGROUNDS[process]:
            raise ValueError("Wrong recorded background")
        if case["source"]["centroid_frame"] != ("slide_global" if standard else "tile_local"):
            raise ValueError("Wrong pilot coordinate frame")
        if case["baseline_mask"]["lineage"] != (
            "standard_original" if standard else "released_current"
        ):
            raise ValueError("Wrong baseline mask lineage")


def inventory(args) -> dict:
    paper_root = no_symlinks(args.paper_root)
    cell_root = no_symlinks(args.cell_root)
    standard_path = input_file(cell_root / STANDARD_REL)
    target_path = input_file(paper_root / TARGET_REL)
    pins = {str(standard_path): STANDARD_SHA, str(target_path): TARGET_SHA}
    verify_pins(pins)
    standard_rows = read_csv(standard_path)
    target_rows = [row for row in read_csv(target_path) if row.get("species") == "D. fuscus"]
    cases = [asdict(case_from_row(row, True, cell_root, args.seed, pins)) for row in standard_rows]
    cases += [asdict(case_from_row(row, False, cell_root, args.seed, pins)) for row in target_rows]
    validate_case_set(cases)
    standard_pins = STANDARD_PINS_PATH
    if digest(standard_pins) != STANDARD_PINS_SHA:
        raise ValueError("Frozen standard pin manifest changed")
    frozen = {
        str(cell_root / row["path"]): row["sha256"]
        for row in load_json(standard_pins)["sources"]
        if row["root"] == "cell"
    }
    for case in cases:
        if case["process"] in {"413", "414"}:
            for source in (case["source"], case["baseline_mask"]):
                if frozen.get(source["path"]) != source["sha256"]:
                    raise ValueError("Standard source differs from frozen source pins")
    pins[str(standard_pins)] = STANDARD_PINS_SHA
    verify_pins(pins)
    for (raw_path, mask_path), group in group_cases(cases):
        raw = read_array(Path(raw_path), raw=True)
        labels = read_array(Path(mask_path))
        for case in group:
            _, bbox = native_membership(case, raw, labels)
            case["source_bbox"] = list(bbox)
    verify_pins(pins)
    protected = [str(cell_root)] + [
        str(paper_root / name)
        for name in ("analyses", "data", "input_data", "outputs", "Publication", "results", ".venv")
    ]
    out = fresh_destination(args.out, [Path(path) for path in [*pins, *protected]])
    out.mkdir()
    args._created_output = out
    result = {
        "schema": SCHEMA,
        "kind": "inventory",
        "pilot": PILOT,
        "seed": args.seed,
        "coordinate_policy": "standards_slide_global; targets_tile_local; crops_from_full_label_bbox",
        "source_pins": pins,
        "protected_roots": protected,
        "cases": cases,
    }
    write_json(out / "inventory.json", result)
    return {"inventory": str(out / "inventory.json"), "cases": len(cases)}


def load_inventory(path: Path) -> dict:
    result = load_json(path)
    if (
        result.get("schema") != SCHEMA
        or result.get("kind") != "inventory"
        or result.get("pilot") != PILOT
    ):
        raise ValueError("Unsupported inventory")
    validate_case_set(result["cases"])
    verify_pins(result["source_pins"])
    standard_paths = [
        Path(name)
        for name, pin in result["source_pins"].items()
        if name.endswith(STANDARD_REL) and pin == STANDARD_SHA
    ]
    target_paths = [
        Path(name)
        for name, pin in result["source_pins"].items()
        if name.endswith(TARGET_REL) and pin == TARGET_SHA
    ]
    if len(standard_paths) != 1 or len(target_paths) != 1:
        raise ValueError("Inventory lacks the frozen authority tables")
    cell_root = standard_paths[0].parents[len(Path(STANDARD_REL).parts) - 1]
    paper_root = target_paths[0].parents[len(Path(TARGET_REL).parts) - 1]
    expected_protected = [str(cell_root)] + [
        str(paper_root / name)
        for name in ("analyses", "data", "input_data", "outputs", "Publication", "results", ".venv")
    ]
    if (
        result["protected_roots"] != expected_protected
        or result["coordinate_policy"]
        != "standards_slide_global; targets_tile_local; crops_from_full_label_bbox"
    ):
        raise ValueError("Inventory guard/coordinate policy differs from source roots")
    pins = result["source_pins"].copy()
    expected = [
        asdict(case_from_row(row, True, cell_root, result["seed"], pins))
        for row in read_csv(standard_paths[0])
    ]
    expected += [
        asdict(case_from_row(row, False, cell_root, result["seed"], pins))
        for row in read_csv(target_paths[0])
        if row.get("species") == "D. fuscus"
    ]
    original_by_id = {case["object_id"]: case for case in expected}
    if set(original_by_id) != {case["object_id"] for case in result["cases"]}:
        raise ValueError("Inventory differs from frozen authoritative identities")
    for case in result["cases"]:
        source = original_by_id[case["object_id"]]
        if any(case.get(key) != value for key, value in source.items() if key != "source_bbox"):
            raise ValueError("Inventory scalar/source contract differs from authority tables")
    frozen = {
        str(cell_root / row["path"]): row["sha256"]
        for row in load_json(STANDARD_PINS_PATH)["sources"]
        if row["root"] == "cell"
    }
    if digest(STANDARD_PINS_PATH) != STANDARD_PINS_SHA:
        raise ValueError("Frozen standard pin manifest changed")
    for case in result["cases"]:
        if case["process"] in {"413", "414"}:
            for source in (case["source"], case["baseline_mask"]):
                if frozen.get(source["path"]) != source["sha256"]:
                    raise ValueError("Inventory standard source differs from frozen pins")
    return result


def backend():
    return importlib.import_module("recovery_backend")


def tools_for_run():
    return backend().load_toolchain(Path(__file__).resolve().parents[2])


def measure(raw: np.ndarray, membership: np.ndarray, background: float, tools) -> dict:
    result = backend().measure_pixels(
        np.asarray(raw[membership], dtype=np.uint8), background, tools
    )
    if result["count"] != int(np.count_nonzero(membership)) or not math.isfinite(result["iod"]):
        raise ValueError("Native reducer returned invalid count/IOD")
    return result


def pack(args) -> dict:
    inv_path = input_file(args.inventory)
    inv = load_inventory(inv_path)
    ui = input_file(UI_PATH)
    inputs = {**inv["source_pins"], str(inv_path): digest(inv_path), str(ui): digest(ui)}
    out = fresh_destination(args.out, [Path(path) for path in [*inputs, *inv["protected_roots"]]])
    margin = integer(args.margin, "crop margin")
    if margin > 2048:
        raise ValueError("Crop margin exceeds supported limit")
    tools = tools_for_run()
    out.mkdir()
    args._created_output = out
    (out / "cases").mkdir()
    (out / "reviewer").mkdir()
    (out / "reviewer" / "assets").mkdir()
    cases = []
    for (raw_path, mask_path), group in group_cases(inv["cases"]):
        raw = read_array(Path(raw_path), raw=True)
        labels = read_array(Path(mask_path))
        for original in group:
            membership, bbox = native_membership(original, raw, labels)
            if list(bbox) != original["source_bbox"]:
                raise ValueError("Original source bounding box changed")
            replay = measure(raw, membership, original["recorded_i_bg"], tools)
            if abs(replay["iod"] - original["stored_iod"]) > 1e-6:
                raise ValueError(f"Recorded-background IOD replay mismatch: {original['case_id']}")
            shared = measure(raw, membership, 255, tools)
            x0, y0 = max(0, bbox[0] - margin), max(0, bbox[1] - margin)
            x1, y1 = min(raw.shape[1], bbox[2] + margin), min(raw.shape[0], bbox[3] + margin)
            crop = raw[y0:y1, x0:x1]
            mask = membership[y0:y1, x0:x1]
            case_id = original["case_id"]
            private = out / "cases" / case_id
            private.mkdir()
            paths = {
                "raw": private / "raw.png",
                "mask": private / "original_mask.png",
                "known_nuclei": private / "known_nuclei.png",
            }
            save_png(paths["raw"], crop)
            save_png(paths["mask"], mask.astype(np.uint8) * 255)
            save_png(paths["known_nuclei"], (labels[y0:y1, x0:x1] > 0).astype(np.uint8) * 255)
            raw_hash, mask_hash = digest(paths["raw"]), digest(paths["mask"])
            crop_contract = {
                "case_id": case_id,
                "source_sha256": original["source"]["sha256"],
                "source_mask_sha256": original["baseline_mask"]["sha256"],
                "box_xyxy": [x0, y0, x1, y1],
                "raw_sha256": raw_hash,
                "mask_sha256": mask_hash,
                "label": original["baseline_mask"]["label"],
            }
            state = value_hash(crop_contract)
            copied = out / "reviewer" / "assets" / case_id
            copied.mkdir()
            shutil.copyfile(paths["raw"], copied / "raw.png")
            shutil.copyfile(paths["mask"], copied / "mask.png")
            case = {
                **original,
                "crop": crop_contract,
                "state_sha256": state,
                "source_label_touches_edge": bbox[0] == 0
                or bbox[1] == 0
                or bbox[2] == raw.shape[1]
                or bbox[3] == raw.shape[0],
                "context_padding_clipped": x0 == 0
                or y0 == 0
                or x1 == raw.shape[1]
                or y1 == raw.shape[0],
                "whole_label_covered": int(mask.sum()) == original["stored_area_px"],
                "recorded_background_replay": replay,
                "shared_255_sensitivity": shared,
            }
            write_json(private / "case.json", case)
            cases.append(case)
    cases.sort(
        key=lambda case: hashlib.sha256(
            (inv["seed"] + "\0shuffle\0" + case["case_id"]).encode()
        ).hexdigest()
    )
    public = [
        {
            "case_id": case["case_id"],
            "raw": f"assets/{case['case_id']}/raw.png",
            "mask": f"assets/{case['case_id']}/mask.png",
            "width": case["crop"]["box_xyxy"][2] - case["crop"]["box_xyxy"][0],
            "height": case["crop"]["box_xyxy"][3] - case["crop"]["box_xyxy"][1],
            "state_sha256": case["state_sha256"],
        }
        for case in cases
    ]
    write_json(out / "reviewer" / "cases.json", public)
    shutil.copyfile(ui, out / "reviewer" / "index.html")
    (out / "reviewer" / "protocol.txt").write_text(PROTOCOL_TEXT)
    write_json(
        out / "provenance.private.json", {"cases": cases, "legacy_fields_are_evidence_only": True}
    )
    verify_pins(inputs)
    assets = {
        str(path.relative_to(out)): digest(path)
        for path in sorted(out.rglob("*"))
        if path.is_file()
    }
    packet = {
        "schema": SCHEMA,
        "kind": "packet",
        "pilot": PILOT,
        "inventory_sha256": digest(inv_path),
        "source_pins": inputs,
        "protected_roots": inv["protected_roots"],
        "coordinate_policy": inv["coordinate_policy"],
        "crop_margin_px": margin,
        "cases": cases,
        "assets": assets,
        "toolchain_identity": tools.identity,
        "proposals": [],
        "decisions": [],
        "scientific_status": "diagnostic_original_only; no calibration or denominator promotion",
    }
    write_json(out / "packet.json", packet)
    return {
        "packet": str(out),
        "packet_sha256": digest(out / "packet.json"),
        "cases": len(cases),
        "authored_proposals": 0,
        "authored_decisions": 0,
    }


def relative_asset(root: Path, name: str) -> Path:
    if not isinstance(name, str) or Path(name).is_absolute() or ".." in Path(name).parts:
        raise ValueError("Invalid relative asset path")
    return input_file(root / name)


def load_packet(path: Path, *, arrays: bool = True) -> dict:
    root = no_symlinks(path)
    packet = load_json(root / "packet.json")
    if (
        packet.get("schema") != SCHEMA
        or packet.get("kind") != "packet"
        or packet.get("pilot") != PILOT
    ):
        raise ValueError("Unsupported packet")
    validate_case_set(packet["cases"])
    verify_pins(packet["source_pins"])
    inventories = [
        Path(name)
        for name, pin in packet["source_pins"].items()
        if pin == packet["inventory_sha256"] and Path(name).name == "inventory.json"
    ]
    if len(inventories) != 1:
        raise ValueError("Packet lacks the exact original inventory")
    inv = load_inventory(inventories[0])
    if (
        packet["protected_roots"] != inv["protected_roots"]
        or packet["coordinate_policy"] != inv["coordinate_policy"]
    ):
        raise ValueError("Packet guard/coordinate policy differs from inventory")
    original_by_id = {case["case_id"]: case for case in inv["cases"]}
    for case in packet["cases"]:
        original = original_by_id.get(case["case_id"])
        if original is None or any(case.get(key) != value for key, value in original.items()):
            raise ValueError("Packet original scalar/source contract differs from inventory")
    if packet.get("proposals") != [] or packet.get("decisions") != []:
        raise ValueError("Original packet must have blank authored state")
    allowed_assets = {
        "provenance.private.json",
        "reviewer/cases.json",
        "reviewer/index.html",
        "reviewer/protocol.txt",
    }
    for case in packet["cases"]:
        case_id = case["case_id"]
        allowed_assets.update(
            f"cases/{case_id}/{name}"
            for name in ("case.json", "raw.png", "original_mask.png", "known_nuclei.png")
        )
        allowed_assets.update(
            f"reviewer/assets/{case_id}/{name}" for name in ("raw.png", "mask.png")
        )
    if set(packet["assets"]) != allowed_assets:
        raise ValueError("Packet asset roles differ from the exact blinded 82-case layout")
    ui = input_file(UI_PATH)
    ui_hash = digest(ui)
    if packet["source_pins"].get(str(ui)) != ui_hash:
        raise ValueError("Packet lacks the mandatory current UI source pin")
    if digest(root / "reviewer" / "index.html") != ui_hash:
        raise ValueError("Reviewer HTML differs from the pinned UI source")
    if (root / "reviewer" / "protocol.txt").read_bytes() != PROTOCOL_TEXT.encode():
        raise ValueError("Reviewer protocol differs from the canonical blinded text")
    if load_json(root / "provenance.private.json") != {
        "cases": packet["cases"],
        "legacy_fields_are_evidence_only": True,
    }:
        raise ValueError("Private provenance differs from the exact packet crosswalk")
    expected = allowed_assets | {"packet.json", "receipt.json"}
    allowed_directories = {
        str(parent)
        for name in allowed_assets
        for parent in Path(name).parents
        if parent != Path(".")
    }
    actual = set()
    for entry in root.rglob("*"):
        no_symlinks(entry)
        if entry.is_file():
            actual.add(str(entry.relative_to(root)))
        elif entry.is_dir() and str(entry.relative_to(root)) not in allowed_directories:
            raise ValueError("Unexpected packet directory outside the blinded layout")
    if actual - expected or set(packet["assets"]) - actual:
        raise ValueError("Missing/extra packet assets")
    for name, expected_hash in packet["assets"].items():
        if digest(relative_asset(root, name)) != expected_hash:
            raise ValueError(f"Packet asset hash mismatch: {name}")
    public = load_json(root / "reviewer" / "cases.json")
    if (
        not isinstance(public, list)
        or len(public) != 82
        or any(set(row) != CASE_KEYS for row in public)
    ):
        raise ValueError("Reviewer JSON has invalid/unblinded fields or missing cases")
    if [row["case_id"] for row in public] != [case["case_id"] for case in packet["cases"]]:
        raise ValueError("Reviewer case set/order differs from private crosswalk")
    for case, public_case in zip(packet["cases"], public, strict=True):
        if (
            value_hash(case["crop"]) != case["state_sha256"]
            or public_case["state_sha256"] != case["state_sha256"]
        ):
            raise ValueError("Case state hash mismatch")
        if load_json(root / "cases" / case["case_id"] / "case.json") != case:
            raise ValueError("Case record differs from packet")
        if (
            public_case["raw"] != f"assets/{case['case_id']}/raw.png"
            or public_case["mask"] != f"assets/{case['case_id']}/mask.png"
        ):
            raise ValueError("Reviewer asset identity differs from case")
        for field, filename in (("raw", "raw.png"), ("mask", "original_mask.png")):
            if digest(root / "reviewer" / public_case[field]) != digest(
                root / "cases" / case["case_id"] / filename
            ):
                raise ValueError("Reviewer pixels/mask differ from private crop")
        x0, y0, x1, y1 = case["crop"]["box_xyxy"]
        if public_case["width"] != x1 - x0 or public_case["height"] != y1 - y0:
            raise ValueError("Reviewer dimensions differ from crop")
        private = root / "cases" / case["case_id"]
        if (
            digest(private / "raw.png") != case["crop"]["raw_sha256"]
            or digest(private / "original_mask.png") != case["crop"]["mask_sha256"]
        ):
            raise ValueError("Crop raw/mask hash differs from original state")
    if arrays:
        tools = tools_for_run()
        if tools.identity != packet["toolchain_identity"]:
            raise ValueError("Numeric toolchain identity differs from packet")
        for (raw_path, mask_path), group in group_cases(packet["cases"]):
            raw = read_array(Path(raw_path), raw=True)
            labels = read_array(Path(mask_path))
            for case in group:
                membership, bbox = native_membership(case, raw, labels)
                if list(bbox) != case["source_bbox"]:
                    raise ValueError("Source bounding box differs from actual full label")
                x0, y0, x1, y1 = case["crop"]["box_xyxy"]
                margin = integer(packet["crop_margin_px"], "crop margin")
                expected_box = [
                    max(0, bbox[0] - margin),
                    max(0, bbox[1] - margin),
                    min(raw.shape[1], bbox[2] + margin),
                    min(raw.shape[0], bbox[3] + margin),
                ]
                if [x0, y0, x1, y1] != expected_box:
                    raise ValueError("Crop bounds differ from full-label/fixed-margin policy")
                private = root / "cases" / case["case_id"]
                if not np.array_equal(read_array(private / "raw.png", raw=True), raw[y0:y1, x0:x1]):
                    raise ValueError("Crop pixels changed from original source")
                if not np.array_equal(
                    read_array(private / "original_mask.png", binary=True) > 0,
                    membership[y0:y1, x0:x1],
                ):
                    raise ValueError("Original mask membership changed")
                if not np.array_equal(
                    read_array(private / "known_nuclei.png", binary=True) > 0,
                    labels[y0:y1, x0:x1] > 0,
                ):
                    raise ValueError("Known-nucleus support changed")
                for field, background in (
                    ("recorded_background_replay", case["recorded_i_bg"]),
                    ("shared_255_sensitivity", 255),
                ):
                    measured = measure(
                        raw[y0:y1, x0:x1], membership[y0:y1, x0:x1], background, tools
                    )
                    if measured != case[field]:
                        raise ValueError("Stored native measurement differs from fresh reduction")
                if abs(case["recorded_background_replay"]["iod"] - case["stored_iod"]) > 1e-6:
                    raise ValueError("Original signed IOD replay mismatch")
    verify_pins(packet["source_pins"])
    return packet


def validate(args) -> dict:
    packet = load_packet(args.packet)
    verify_pins(packet_file_pins(args.packet, packet))
    return {
        "packet_sha256": digest(args.packet / "packet.json"),
        "cases": len(packet["cases"]),
        "pixel_equality": True,
        "mask_equality": True,
        "numeric_replay": True,
    }


def named_time(author, timestamp) -> None:
    if not isinstance(author, str) or not author.strip() or len(author) > 200:
        raise ValueError("Named author/reviewer is required")
    if not isinstance(timestamp, str) or not re.fullmatch(
        r"\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}(?:\.\d+)?(?:Z|[+-]\d{2}:\d{2})", timestamp
    ):
        raise ValueError("Authored time requires timezone ISO8601")
    try:
        parsed = datetime.fromisoformat(timestamp)
    except ValueError as error:
        raise ValueError("Invalid authored time") from error
    if parsed.utcoffset() is None:
        raise ValueError("Authored time requires timezone")


def proposal_state(proposal: dict) -> str:
    return value_hash(
        {
            key: value
            for key, value in proposal.items()
            if key not in {"state_sha256", "proposal_id", "root"}
        }
    )


def packet_file_pins(root: Path, packet: dict) -> dict[str, str]:
    return {
        **packet["source_pins"],
        str(root / "packet.json"): digest(root / "packet.json"),
        **{str(relative_asset(root, name)): pin for name, pin in packet["assets"].items()},
    }


def prepare_editor(args) -> dict:
    packet = load_packet(args.packet)
    matches = [case for case in packet["cases"] if case["case_id"] == args.case]
    if len(matches) != 1:
        raise ValueError("Editor workspace identifies unknown case")
    case = matches[0]
    inputs = packet_file_pins(args.packet, packet)
    out = fresh_destination(
        args.out, [args.packet, *[Path(path) for path in [*inputs, *packet["protected_roots"]]]]
    )
    tools = tools_for_run()
    out.mkdir()
    args._created_output = out
    private = args.packet / "cases" / case["case_id"]
    shutil.copyfile(private / "raw.png", out / "raw.png")
    shutil.copyfile(private / "original_mask.png", out / "original_mask.png")
    sidecar = {
        "packet_sha256": digest(args.packet / "packet.json"),
        "case_id": case["case_id"],
        "crop_sha256": case["crop"]["raw_sha256"],
        "parent_state_sha256": case["state_sha256"],
        "author": "",
        "time": "",
        "rationale": "",
    }
    write_json(out / "authorship.template.json", sidecar)
    private_home = out / ".imagej-home"
    private_prefs = out / ".imagej-prefs"
    private_home.mkdir()
    private_prefs.mkdir()
    launch = [
        str(tools.java),
        f"-Duser.home={private_home}",
        f"-Djava.util.prefs.userRoot={private_prefs}",
        "-Dsun.awt.disableGtkFileDialogs=true",
        "-cp",
        str(tools.ij_jar),
        "ij.ImageJ",
        str(out / "raw.png"),
    ]
    write_json(
        out / "editor.json",
        {
            "schema": SCHEMA,
            "kind": "editor_workspace",
            "case_id": case["case_id"],
            "launch_argv": launch,
            "toolchain_identity": tools.identity,
            "authored_proposals": 0,
            "authored_decisions": 0,
        },
    )
    launch_command = shlex.join(
        [
            "pixi",
            "run",
            "--manifest-path",
            str(Path(__file__).resolve().parents[2] / "pixi.toml"),
            "-e",
            "recovery",
            "--clean-env",
            "--locked",
            "--no-install",
            "--executable",
            "python",
            "-c",
            'import json, subprocess; subprocess.run(json.load(open("editor.json"))["launch_argv"], check=True)',
        ]
    )
    (out / "README.txt").write_text(
        "These are disposable diagnostic copies. Preserve the raw pixels and use crop-local coordinates.\n"
        "Open raw.png in the pinned ImageJ editor. Open original_mask.png separately as a reference.\n"
        f"Run this command from the workspace to launch the exact editor.json argv:\n{launch_command}\n"
        "ImageJ home and Java preferences stay inside this copy.\n"
        "In ROI Manager, author one supported closed nucleus selection named nucleus. Optional background "
        "regions must be named background_001, background_002, etc. A single selection may be saved as "
        "nucleus.roi (at most 1 MiB). Multiple selections must be saved as a ZIP containing nucleus.roi "
        "and optional background_NNN.roi entries. Do not include the reference mask as an ROI.\n"
        "Copy authorship.template.json to authorship.json and explicitly fill your name, timezone ISO8601 time "
        "and rationale. Import the authored nucleus.roi or ZIP and sidecar with import-roi. "
        "An import is a proposal only.\n"
        "No ROI or decision has been authored by this preparation command.\n"
    )
    verify_pins(inputs)
    return {
        "workspace": str(out),
        "case_id": case["case_id"],
        "launch_argv": launch,
        "authored_proposals": 0,
        "authored_decisions": 0,
    }


def import_roi(args) -> dict:
    packet = load_packet(args.packet)
    packet_hash = digest(args.packet / "packet.json")
    sidecar = load_json(args.sidecar)
    required = {
        "packet_sha256",
        "case_id",
        "crop_sha256",
        "parent_state_sha256",
        "author",
        "time",
        "rationale",
    }
    if set(sidecar) != required:
        raise ValueError("Proposal sidecar requires exact authored binding fields")
    named_time(sidecar["author"], sidecar["time"])
    if not isinstance(sidecar["rationale"], str) or not sidecar["rationale"].strip():
        raise ValueError("Authored proposal rationale is required")
    matches = [case for case in packet["cases"] if case["case_id"] == sidecar["case_id"]]
    if len(matches) != 1 or sidecar["packet_sha256"] != packet_hash:
        raise ValueError("Proposal identifies unknown packet/case")
    case = matches[0]
    if (
        sidecar["crop_sha256"] != case["crop"]["raw_sha256"]
        or sidecar["parent_state_sha256"] != case["state_sha256"]
    ):
        raise ValueError("Stale proposal state or changed crop hash")
    roi = input_file(args.roi)
    single_roi = None
    if roi.suffix.lower() == ".roi":
        if roi.name != "nucleus.roi":
            raise ValueError("A single ROI must be named exactly nucleus.roi")
        with roi.open("rb") as handle:
            single_roi = handle.read(1024 * 1024 + 1)
        if not single_roi or len(single_roi) > 1024 * 1024:
            raise ValueError("Single nucleus.roi requires nonempty bytes of at most 1 MiB")
    input_pins = {
        **packet_file_pins(args.packet, packet),
        str(args.packet / "packet.json"): packet_hash,
        str(roi): digest(roi),
        str(input_file(args.sidecar)): digest(args.sidecar),
    }
    if single_roi is not None and hashlib.sha256(single_roi).hexdigest() != input_pins[str(roi)]:
        raise ValueError("Single ROI source changed while reading")
    out = fresh_destination(
        args.out, [args.packet, *[Path(path) for path in [*input_pins, *packet["protected_roots"]]]]
    )
    tools = tools_for_run()
    out.mkdir()
    args._created_output = out
    if single_roi is None:
        shutil.copyfile(roi, out / "roi.zip")
    else:
        entry = zipfile.ZipInfo("nucleus.roi", date_time=(1980, 1, 1, 0, 0, 0))
        entry.compress_type = zipfile.ZIP_STORED
        with zipfile.ZipFile(out / "roi.zip", "x", compression=zipfile.ZIP_STORED) as archive:
            archive.writestr(entry, single_roi)
    shutil.copyfile(args.sidecar, out / "authorship.json")
    x0, y0, x1, y1 = case["crop"]["box_xyxy"]
    rasters = backend().rasterize_rois(out / "roi.zip", x1 - x0, y1 - y0, out / "rasters", tools)
    if "nucleus" not in rasters or any(
        key != "nucleus" and not re.fullmatch(r"background_\d{3}", key) for key in rasters
    ):
        raise ValueError("ROI bridge returned invalid roles")
    private = args.packet / "cases" / case["case_id"]
    raw = read_array(private / "raw.png", raw=True)
    original = read_array(private / "original_mask.png", binary=True) > 0
    known = read_array(private / "known_nuclei.png", binary=True) > 0
    nucleus = read_array(rasters["nucleus"], binary=True) > 0
    if nucleus.shape != raw.shape or not nucleus.any():
        raise ValueError("Proposed nucleus requires nonempty exact crop grid")
    backgrounds = {}
    for key, path in sorted(rasters.items()):
        if key == "nucleus":
            continue
        region = read_array(path, binary=True) > 0
        if region.shape != raw.shape or not region.any() or np.any(region & (known | nucleus)):
            raise ValueError(
                "Background ROI is empty/outside grid or overlaps a known/proposed nucleus"
            )
        pixels = raw[region]
        summary = backend().measure_pixels(np.asarray(pixels, dtype=np.uint8), 255, tools)
        candidates = {
            "mean": summary["raw_sum"] / summary["count"],
            "median": summary["median_pixel"],
        }
        if any(not 1 <= value <= 255 for value in candidates.values()):
            raise ValueError("Candidate background must remain within [1,255]")
        backgrounds[key] = {
            "area_px": summary["count"],
            "raw_mean": candidates["mean"],
            "raw_median": candidates["median"],
            "status": "authored_candidate_reference_only",
            "nucleus_sensitivity": {
                name: measure(raw, nucleus, value, tools) for name, value in candidates.items()
            },
        }
    original_area = int(original.sum())
    proposed_area = int(nucleus.sum())
    added = int((nucleus & ~original).sum())
    removed = int((original & ~nucleus).sum())
    proposal = {
        "schema": SCHEMA,
        "kind": "proposal",
        **sidecar,
        "roi_sha256": digest(out / "roi.zip"),
        "nucleus_sha256": digest(rasters["nucleus"]),
        "background_sha256": {
            key: digest(path) for key, path in sorted(rasters.items()) if key != "nucleus"
        },
        "rasters": {key: str(path.relative_to(out)) for key, path in sorted(rasters.items())},
        "geometry": {
            "original_area_px": original_area,
            "proposed_area_px": proposed_area,
            "added_px": added,
            "removed_px": removed,
        },
        "recorded_background": measure(raw, nucleus, case["recorded_i_bg"], tools),
        "shared_255_sensitivity": measure(raw, nucleus, 255, tools),
        "candidate_backgrounds": backgrounds,
        "toolchain_identity": tools.identity,
        "decision": None,
        "input_pins": input_pins,
    }
    proposal["state_sha256"] = proposal_state(proposal)
    proposal["proposal_id"] = "P" + proposal["state_sha256"][:20]
    verify_pins(input_pins)
    write_json(out / "proposal.json", proposal)
    return {
        "proposal": str(out),
        "proposal_id": proposal["proposal_id"],
        "state_sha256": proposal["state_sha256"],
        "decision": None,
    }


def load_proposals(paths: list[Path], packet_root: Path, packet: dict) -> dict[str, dict]:
    packet_hash = digest(packet_root / "packet.json")
    cases = {case["case_id"]: case for case in packet["cases"]}
    result = {}
    for path in paths:
        root = no_symlinks(path)
        proposal = load_json(root / "proposal.json")
        if (
            proposal.get("schema") != SCHEMA
            or proposal.get("kind") != "proposal"
            or proposal.get("packet_sha256") != packet_hash
        ):
            raise ValueError("Unsupported proposal or different packet")
        case = cases.get(proposal["case_id"])
        if (
            not case
            or proposal["crop_sha256"] != case["crop"]["raw_sha256"]
            or proposal["parent_state_sha256"] != case["state_sha256"]
        ):
            raise ValueError("Unknown/stale proposal identity")
        named_time(proposal["author"], proposal["time"])
        if not proposal["rationale"].strip() or proposal.get("decision") is not None:
            raise ValueError("Proposal cannot manufacture acceptance")
        verify_pins(proposal["input_pins"])
        sidecar = load_json(root / "authorship.json")
        if any(proposal.get(key) != value for key, value in sidecar.items()) or set(sidecar) != {
            "packet_sha256",
            "case_id",
            "crop_sha256",
            "parent_state_sha256",
            "author",
            "time",
            "rationale",
        }:
            raise ValueError("Proposal authorship changed")
        if digest(root / "roi.zip") != proposal["roi_sha256"]:
            raise ValueError("Proposal ROI bytes changed")
        for key, name in proposal["rasters"].items():
            expected = (
                proposal["nucleus_sha256"]
                if key == "nucleus"
                else proposal["background_sha256"][key]
            )
            if digest(relative_asset(root, name)) != expected:
                raise ValueError("Proposal raster changed")
        if (
            proposal_state(proposal) != proposal["state_sha256"]
            or proposal["proposal_id"] != "P" + proposal["state_sha256"][:20]
        ):
            raise ValueError("Proposal state/identity changed")
        if proposal["state_sha256"] in result:
            raise ValueError("Duplicate proposal state")
        result[proposal["state_sha256"]] = {**proposal, "root": str(root)}
    return result


def validated_events(value, states: dict[str, set[str]]) -> list[dict]:
    events = value if isinstance(value, list) else [value]
    result = []
    seen = set()
    for event in events:
        if not isinstance(event, dict) or set(event) != EVENT_KEYS:
            raise ValueError("Decision requires exact authored event fields")
        named_time(event["reviewer"], event["time"])
        if event["case_id"] not in states or event["state_sha256"] not in states[event["case_id"]]:
            raise ValueError("Unknown case or stale decision state")
        if event["decision"] not in {"keep", "problem", "unsure"}:
            raise ValueError("Decision must be explicitly keep/problem/unsure")
        reasons = event["reasons"]
        if (
            not isinstance(reasons, list)
            or not reasons
            or any(not isinstance(reason, str) or not reason.strip() for reason in reasons)
            or len(set(reasons)) != len(reasons)
        ):
            raise ValueError("Every decision requires explicit unique nonempty reasons")
        if not isinstance(event["note"], str):
            raise TypeError("Decision note must be text")
        key = (event["case_id"], event["state_sha256"], event["reviewer"])
        if key in seen:
            raise ValueError("Duplicate/colliding reviewer event; no last-writer-wins")
        seen.add(key)
        result.append(event)
    return result


def known_states(packet: dict, proposals: dict) -> dict[str, set[str]]:
    states = {case["case_id"]: {case["state_sha256"]} for case in packet["cases"]}
    for state, proposal in proposals.items():
        states[proposal["case_id"]].add(state)
    return states


def load_policy(path: Path | None, states: dict[str, set[str]]) -> dict:
    policy = load_json(path) if path else {"name": "original_baseline_v1", "state_selection": {}}
    if (
        not isinstance(policy, dict)
        or set(policy) != {"name", "state_selection"}
        or not isinstance(policy["name"], str)
        or not policy["name"].strip()
        or not isinstance(policy["state_selection"], dict)
    ):
        raise ValueError("Review/report policy requires name and explicit state_selection")
    for case_id, state in policy["state_selection"].items():
        if not isinstance(state, str) or case_id not in states or state not in states[case_id]:
            raise ValueError("Review/report policy selects unknown/stale state")
    return policy


def review_queue(args) -> dict:
    packet = load_packet(args.packet)
    proposals = load_proposals(args.proposals, args.packet, packet)
    policy = load_policy(args.policy, known_states(packet, proposals))
    inputs = packet_file_pins(args.packet, packet)
    for proposal in proposals.values():
        root = Path(proposal["root"])
        inputs.update(proposal["input_pins"])
        inputs[str(root / "proposal.json")] = digest(root / "proposal.json")
        inputs[str(root / "roi.zip")] = proposal["roi_sha256"]
        inputs[str(root / "authorship.json")] = digest(root / "authorship.json")
        for role, name in proposal["rasters"].items():
            inputs[str(relative_asset(root, name))] = (
                proposal["nucleus_sha256"]
                if role == "nucleus"
                else proposal["background_sha256"][role]
            )
    if args.policy:
        inputs[str(input_file(args.policy))] = digest(args.policy)
    out = fresh_destination(
        args.out,
        [
            args.packet,
            *args.proposals,
            *[Path(path) for path in [*inputs, *packet["protected_roots"]]],
        ],
    )
    out.mkdir()
    args._created_output = out
    reviewer = out / "reviewer"
    reviewer.mkdir()
    (reviewer / "assets").mkdir()
    public = []
    selected_proposals = 0
    for case in packet["cases"]:
        case_id = case["case_id"]
        state = policy["state_selection"].get(case_id, case["state_sha256"])
        proposal = proposals.get(state)
        private = args.packet / "cases" / case_id
        raw = private / "raw.png"
        mask = (
            private / "original_mask.png"
            if proposal is None
            else relative_asset(Path(proposal["root"]), proposal["rasters"]["nucleus"])
        )
        x0, y0, x1, y1 = case["crop"]["box_xyxy"]
        membership = read_array(mask, binary=True)
        if membership.shape != (y1 - y0, x1 - x0) or not membership.any():
            raise ValueError("Selected proposal mask differs from the exact crop grid")
        destination = reviewer / "assets" / case_id
        destination.mkdir()
        shutil.copyfile(raw, destination / "raw.png")
        shutil.copyfile(mask, destination / "mask.png")
        if digest(destination / "raw.png") != case["crop"]["raw_sha256"] or digest(
            destination / "mask.png"
        ) != (case["crop"]["mask_sha256"] if proposal is None else proposal["nucleus_sha256"]):
            raise ValueError("Review queue copies differ from the selected source state")
        public.append(
            {
                "case_id": case_id,
                "raw": f"assets/{case_id}/raw.png",
                "mask": f"assets/{case_id}/mask.png",
                "width": x1 - x0,
                "height": y1 - y0,
                "state_sha256": state,
            }
        )
        selected_proposals += proposal is not None
    write_json(reviewer / "cases.json", public)
    shutil.copyfile(args.packet / "reviewer" / "index.html", reviewer / "index.html")
    shutil.copyfile(args.packet / "reviewer" / "protocol.txt", reviewer / "protocol.txt")
    verify_pins(inputs)
    assets = {
        str(path.relative_to(out)): digest(path)
        for path in sorted(reviewer.rglob("*"))
        if path.is_file()
    }
    write_json(
        out / "queue.private.json",
        {
            "schema": SCHEMA,
            "kind": "review_queue",
            "packet_sha256": digest(args.packet / "packet.json"),
            "input_pins": inputs,
            "policy": policy,
            "policy_sha256": value_hash(policy),
            "assets": assets,
            "cases": len(public),
            "supplied_proposals": len(proposals),
            "selected_proposals": selected_proposals,
            "authored_proposals": 0,
            "authored_decisions": 0,
        },
    )
    return {
        "queue": str(out),
        "serve_directory": str(reviewer),
        "cases": len(public),
        "selected_proposals": selected_proposals,
        "authored_proposals": 0,
        "authored_decisions": 0,
    }


def decision(args) -> dict:
    packet = load_packet(args.packet)
    proposals = load_proposals(args.proposals, args.packet, packet)
    events_path = input_file(args.event)
    events_hash = digest(events_path)
    events = validated_events(load_json(events_path), known_states(packet, proposals))
    if not events:
        raise ValueError("Decision export contains no authored events")
    sources = [
        args.packet,
        events_path,
        *args.proposals,
        *[Path(path) for path in [*packet["source_pins"], *packet["protected_roots"]]],
    ]
    out = fresh_destination(args.out, sources)
    out.mkdir()
    args._created_output = out
    record = {
        "schema": SCHEMA,
        "kind": "decision_events",
        "packet_sha256": digest(args.packet / "packet.json"),
        "source_event_sha256": events_hash,
        "events": events,
        "proposal_pins": {
            str(path / "proposal.json"): digest(path / "proposal.json") for path in args.proposals
        },
    }
    verify_pins({**packet_file_pins(args.packet, packet), str(events_path): events_hash})
    write_json(out / "events.json", record)
    return {"events": str(out / "events.json"), "authored_events": len(events)}


def read_event_set(
    paths: list[Path], packet_root: Path, states: dict
) -> tuple[list[dict], dict[str, str]]:
    events = []
    pins = {}
    for path in paths:
        path = input_file(path / "events.json" if path.is_dir() else path)
        value = load_json(path)
        if (
            not isinstance(value, dict)
            or value.get("schema") != SCHEMA
            or value.get("kind") != "decision_events"
            or value.get("packet_sha256") != digest(packet_root / "packet.json")
        ):
            raise ValueError("Report requires immutable recorded decision events")
        verify_pins(
            value.get("proposal_pins", {})
            or {str(packet_root / "packet.json"): digest(packet_root / "packet.json")}
        )
        events.extend(value["events"])
        pins[str(path)] = digest(path)
    return validated_events(events, states), pins


def report(args) -> dict:
    packet = load_packet(args.packet)
    proposals = load_proposals(args.proposals, args.packet, packet)
    states = known_states(packet, proposals)
    events, event_pins = read_event_set(args.events, args.packet, states)
    policy = load_policy(args.policy, states)
    inputs = {**packet_file_pins(args.packet, packet), **event_pins}
    inputs.update(
        {str(path / "proposal.json"): digest(path / "proposal.json") for path in args.proposals}
    )
    if args.policy:
        inputs[str(input_file(args.policy))] = digest(args.policy)
    out = fresh_destination(
        args.out,
        [
            args.packet,
            *args.proposals,
            *[Path(path) for path in [*inputs, *packet["protected_roots"]]],
        ],
    )
    groups = defaultdict(list)
    for event in events:
        groups[(event["case_id"], event["state_sha256"])].append(event)
    rows = []
    counts = Counter()
    by_image = defaultdict(Counter)
    approved = []
    for case in packet["cases"]:
        state = policy["state_selection"].get(case["case_id"], case["state_sha256"])
        selected_events = groups[(case["case_id"], state)]
        decisions = {event["decision"] for event in selected_events}
        status = (
            "missing"
            if not decisions
            else "conflict"
            if len(decisions) > 1
            else next(iter(decisions))
        )
        if status == "keep":
            approved.append(case["case_id"])
        counts[status] += 1
        by_image[(case["cohort"], case["process"])][status] += 1
        proposal = proposals.get(state)
        rows.append(
            {
                "case_id": case["case_id"],
                "object_id": case["object_id"],
                "cohort": case["cohort"],
                "process": case["process"],
                "original_area_px": case["stored_area_px"],
                "original_stored_iod": case["stored_iod"],
                "original_recorded_background": case["recorded_i_bg"],
                "original_recorded_background_replay": case["recorded_background_replay"],
                "original_shared_255_sensitivity": case["shared_255_sensitivity"],
                "selected_state_sha256": state,
                "review_status": status,
                "reviewer_events": selected_events,
                "proposal": None
                if proposal is None
                else {
                    key: proposal[key]
                    for key in (
                        "proposal_id",
                        "geometry",
                        "recorded_background",
                        "shared_255_sensitivity",
                        "candidate_backgrounds",
                    )
                },
            }
        )
    result = {
        "schema": SCHEMA,
        "kind": "report",
        "packet_sha256": digest(args.packet / "packet.json"),
        "input_pins": inputs,
        "policy": policy,
        "policy_sha256": value_hash(policy),
        "rows": rows,
        "authored_proposals": len(proposals),
        "authored_decisions": len(events),
        "review_counts": {
            name: counts[name] for name in ("missing", "keep", "problem", "unsure", "conflict")
        },
        "counts_by_original_image_cohort": [
            {"cohort": cohort, "process": process, "counts": dict(counts)}
            for (cohort, process), counts in sorted(by_image.items())
        ],
        "review_approved_subset": {
            "name": f"{policy['name']}:unanimous_explicit_keep",
            "case_ids": approved,
        },
        "proposal_sensitivities": [
            {
                key: proposal[key]
                for key in (
                    "proposal_id",
                    "case_id",
                    "state_sha256",
                    "geometry",
                    "recorded_background",
                    "shared_255_sensitivity",
                    "candidate_backgrounds",
                )
            }
            for proposal in proposals.values()
        ],
        "unselected_state_events": [
            event
            for event in events
            if event["state_sha256"]
            != policy["state_selection"].get(
                event["case_id"],
                next(
                    case["state_sha256"]
                    for case in packet["cases"]
                    if case["case_id"] == event["case_id"]
                ),
            )
        ],
        "scientific_status": "diagnostic sensitivity only; no pg/1C, pooled denominator, calibration or promotion",
    }
    verify_pins(inputs)
    out.mkdir()
    args._created_output = out
    write_json(out / "report.json", result)
    return {
        "report": str(out / "report.json"),
        "cases": len(rows),
        "review_counts": result["review_counts"],
        "authored_proposals": len(proposals),
        "authored_decisions": len(events),
    }


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    commands = result.add_subparsers(dest="command", required=True)
    command = commands.add_parser("inventory")
    command.add_argument("--paper-root", type=Path, required=True)
    command.add_argument("--cell-root", type=Path, required=True)
    command.add_argument("--pilot", choices=[PILOT], default=PILOT)
    command.add_argument("--seed", default="microscopy-recovery-pilot-v1")
    command.add_argument("--out", type=Path, required=True)
    command.set_defaults(run=inventory)
    command = commands.add_parser("pack")
    command.add_argument("--inventory", type=Path, required=True)
    command.add_argument("--margin", type=int, default=32)
    command.add_argument("--out", type=Path, required=True)
    command.set_defaults(run=pack)
    command = commands.add_parser("validate")
    command.add_argument("--packet", type=Path, required=True)
    command.set_defaults(run=validate)
    command = commands.add_parser("prepare-editor")
    command.add_argument("--packet", type=Path, required=True)
    command.add_argument("--case", required=True)
    command.add_argument("--out", type=Path, required=True)
    command.set_defaults(run=prepare_editor)
    command = commands.add_parser("review-queue")
    command.add_argument("--packet", type=Path, required=True)
    command.add_argument("--proposals", type=Path, nargs="*", default=[])
    command.add_argument("--policy", type=Path)
    command.add_argument("--out", type=Path, required=True)
    command.set_defaults(run=review_queue)
    command = commands.add_parser("import-roi")
    command.add_argument("--packet", type=Path, required=True)
    command.add_argument("--roi", "--proposal", type=Path, required=True)
    command.add_argument("--sidecar", type=Path, required=True)
    command.add_argument("--out", type=Path, required=True)
    command.set_defaults(run=import_roi)
    command = commands.add_parser("decision")
    command.add_argument("--packet", type=Path, required=True)
    command.add_argument("--event", type=Path, required=True)
    command.add_argument("--proposals", type=Path, nargs="*", default=[])
    command.add_argument("--out", type=Path, required=True)
    command.set_defaults(run=decision)
    command = commands.add_parser("report")
    command.add_argument("--packet", type=Path, required=True)
    command.add_argument("--events", type=Path, nargs="*", default=[])
    command.add_argument("--proposals", type=Path, nargs="*", default=[])
    command.add_argument("--policy", type=Path)
    command.add_argument("--out", type=Path, required=True)
    command.set_defaults(run=report)
    return result


def main(argv: list[str] | None = None) -> int:
    args = parser().parse_args(argv)
    receipt = {
        "schema": SCHEMA,
        "command": args.command,
        "argv": argv if argv is not None else sys.argv[1:],
        "started_utc": clock(),
    }
    try:
        summary = args.run(args)
        receipt.update({"finished_utc": clock(), "exit": 0, "summary": summary})
        if hasattr(args, "_created_output"):
            write_json(args._created_output / "receipt.json", receipt)
        print(json.dumps(summary, sort_keys=True, allow_nan=False))
        return 0
    except (ValueError, OSError, KeyError, TypeError, ImportError, RuntimeError) as error:
        receipt.update({"finished_utc": clock(), "exit": 1, "error": str(error)})
        if (
            hasattr(args, "_created_output")
            and not (args._created_output / "receipt.json").exists()
        ):
            try:
                write_json(args._created_output / "receipt.json", receipt)
            except OSError:
                pass
        print(f"recovery: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
