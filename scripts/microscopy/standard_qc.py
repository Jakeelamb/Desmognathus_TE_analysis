#!/usr/bin/env python3
"""Review existing standard masks against preserved native quality helpers."""

from __future__ import annotations

import argparse
import hashlib
import html
import importlib
import importlib.metadata
import json
import os
import platform
import shutil
import subprocess
import sys
from datetime import UTC, datetime
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd

sys.dont_write_bytecode = True
STANDARDS = {"Process_413_green.ome.tiff": "32469", "Process_414_green.ome.tiff": "32470"}
COUNTS = {"Process_413_green.ome.tiff": 31, "Process_414_green.ome.tiff": 16}
MATCH_FEATURES = ("match_log_edge_sharpness", "match_log_relative_ring_noise")
PIXEL_METRICS = ("iod_edge_sharpness", "iod_ring_noise_iqr", "iod_ring_intensity_median")
SHAPE_METRICS = tuple(
    f"{part}_shape_{metric}"
    for part in ("cell", "nucleus")
    for metric in ("smoothness", "ellipse_iou", "solidity")
) + ("raw_focus_lap_var", "raw_focus_grad_mean", "raw_focus_intensity_iqr")
PINS = Path(__file__).with_name("standard_qc_source_pins.json")
PIN_MANIFEST_SHA256 = "d3bfed24e791e61836818572879177ed6b93c203f4b11b344d0881da6c305bfc"
STANDARD_RUN = "output/runs/life_history_pmount_reference_v1"
AUDIT_RELATIVE = "outputs/standard_calibration_review_20261013/pixel_audit_v2"
TEMPLATE_RELATIVE = (
    "path_analysis/data/external/derived/image_quality_matched_genome_iod/"
    "image_quality_matched_nuclei_all_species_sensitivity.csv.gz"
)


def clock() -> str:
    return datetime.now(UTC).isoformat()


def digest(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def require_hash(path: Path, expected: str) -> None:
    actual = digest(path)
    if actual != expected:
        raise ValueError(f"Source hash mismatch: {path}: {actual} != {expected}")


def fresh_destination(path: Path, sources: list[Path]) -> Path:
    out = Path(os.path.abspath(path))
    for component in (out, *out.parents):
        if component.is_symlink():
            raise ValueError(f"Symlinked output or ancestor: {component}")
    if out.exists():
        raise ValueError(f"Output already exists: {out}")
    if not out.parent.is_dir():
        raise ValueError(f"Output parent does not exist: {out.parent}")
    for source in sources:
        source = source.resolve()
        if out == source or out in source.parents or source in out.parents:
            raise ValueError(f"Input/output ancestry overlap: {source}")
    return out


def source_paths(roots: dict[str, Path]) -> dict[Path, str]:
    require_hash(PINS, PIN_MANIFEST_SHA256)
    entries = json.loads(PINS.read_text())["sources"]
    return {roots[row["root"]] / row["path"]: row["sha256"] for row in entries}


def verify_sources(pins: dict[Path, str]) -> None:
    for path, expected in pins.items():
        require_hash(path, expected)


def identity_frame(frame: pd.DataFrame, *, inventory: bool = False) -> pd.DataFrame:
    out = frame.copy()
    required = {"filename", "specimen_id", "species", "image_type", "tile_name"}
    required.add("label" if inventory else "nucleus_label")
    if missing := sorted(required - set(out.columns)):
        raise ValueError(f"Missing identity fields: {missing}")
    for name in ("filename", "specimen_id", "species", "image_type", "tile_name"):
        out[name] = out[name].fillna("").astype(str).str.strip()
    label_name = "label" if inventory else "nucleus_label"
    label = pd.to_numeric(out[label_name], errors="coerce")
    if not (np.isfinite(label) & label.gt(0) & label.eq(np.floor(label))).all():
        raise ValueError("Nucleus labels must be finite positive integers")
    expected = out["filename"].map(STANDARDS)
    if not (
        expected.notna()
        & out["specimen_id"].eq(expected)
        & out["species"].eq("D. fuscus")
        & out["image_type"].eq("pmount")
        & out["tile_name"].str.fullmatch(r"tile_y\d+_x\d+\.tiff")
    ).all():
        raise ValueError("Wrong standard process, specimen, species, image type or tile identity")
    keys = out["filename"] + "::" + out["tile_name"] + "::" + label.astype(int).astype(str)
    if "nucleus_object_id" in out and not out["nucleus_object_id"].astype(str).eq(keys).all():
        raise ValueError("Nucleus identity differs from original process/tile/label")
    if keys.duplicated().any():
        raise ValueError("Duplicate nucleus identity")
    out["nucleus_object_id"] = keys
    out["review_key"] = keys
    out["standard_process_id"] = out["filename"].str.extract(r"Process_(\d+)_", expand=False)
    return out


def validate_census(frame: pd.DataFrame, inventory: pd.DataFrame) -> None:
    if len(inventory) != 47 or inventory["filename"].value_counts().to_dict() != COUNTS:
        raise ValueError("The original standard census must contain exactly 31 + 16 objects")
    if len(frame) != 47 or set(frame["nucleus_object_id"]) != set(inventory["nucleus_object_id"]):
        raise ValueError("Linkage does not preserve the exact 47 original identities")


def validate_template(template: pd.DataFrame, scaling: dict, native) -> pd.DataFrame:
    required = {"quality_template_id", *MATCH_FEATURES, *(f"z_{name}" for name in MATCH_FEATURES)}
    if missing := sorted(required - set(template)):
        raise ValueError(f"Missing template fields: {missing}")
    ids = pd.to_numeric(template["quality_template_id"], errors="coerce")
    if len(template) != 40 or sorted(ids.tolist()) != list(range(1, 41)):
        raise ValueError("Frozen template requires the unique original IDs 1 through 40")
    if set(scaling) != set(MATCH_FEATURES):
        raise ValueError("Frozen scaling must contain exactly the two technical match features")
    for name in MATCH_FEATURES:
        params = scaling[name]
        if set(params) != {"clip_low", "clip_high", "median", "scale_iqr"}:
            raise ValueError(f"Invalid frozen scaling fields: {name}")
        values = np.asarray(list(params.values()), dtype=float)
        if (
            not np.isfinite(values).all()
            or params["scale_iqr"] <= 0
            or params["clip_low"] >= params["clip_high"]
        ):
            raise ValueError(f"Invalid frozen scaling values: {name}")
    numeric = template[list(required - {"quality_template_id"})].apply(
        pd.to_numeric, errors="coerce"
    )
    if not np.isfinite(numeric.to_numpy()).all():
        raise ValueError("Nonfinite frozen template feature")
    result = native.apply_frozen_scaling(template, scaling)
    for name in MATCH_FEATURES:
        if not np.allclose(result[f"z_{name}"], template[f"z_{name}"], rtol=0, atol=1e-12):
            raise ValueError("Frozen template scaling differs from preserved z scores")
    return template.sort_values("quality_template_id").reset_index(drop=True)


def import_native(roots: dict[str, Path], pins: dict[Path, str]):
    sys.path[:0] = [str(roots["cell"] / "scripts"), str(roots["archive"] / "path_analysis/scripts")]
    names = {
        "linkage": "build_cell_nucleus_linkage_report",
        "review": "run_linked_pair_review",
        "extension": "extend_finalized_microscopy_reviewers",
    }
    modules = SimpleNamespace(
        **{key: importlib.import_module(value) for key, value in names.items()}
    )
    for module in tuple(sys.modules.values()):
        source = getattr(module, "__file__", None)
        if not source:
            continue
        path = Path(source).resolve()
        if any(path.is_relative_to(roots[key]) for key in ("cell", "archive")):
            if path not in pins:
                raise ValueError(f"Unpinned native imported source: {path}")
            require_hash(path, pins[path])
    return modules


def load_decisions(paths: list[Path], inventory: pd.DataFrame) -> dict[str, str]:
    valid_keys = set(inventory["review_key"])
    result: dict[str, str] = {}
    for path in paths:
        frame = pd.read_csv(path, dtype=str, keep_default_na=False)
        if not {"species", "review_key", "decision"}.issubset(frame):
            raise ValueError("Manual decisions require species, review_key and decision")
        if frame["review_key"].duplicated().any():
            raise ValueError("Duplicate manual review key within a decision export")
        for row in frame.to_dict("records"):
            if row["species"] != "D. fuscus" or row["review_key"] not in valid_keys:
                raise ValueError("Manual decision identifies an unknown standard")
            if row["decision"] not in {"keep", "problem", "unsure"}:
                raise ValueError("Manual decision must be keep, problem or unsure")
            result[row["review_key"]] = row["decision"]
    return result


def finite_evidence(frame: pd.DataFrame) -> pd.Series:
    required = [
        *SHAPE_METRICS,
        *PIXEL_METRICS,
        "nuc_iod",
        "nuc_area_um2",
        "cell_area_um2",
        "nucleus_mask_overlap_label_count",
        "nucleus_mask_best_overlap_fraction",
    ]
    if missing := sorted(set(required) - set(frame)):
        raise ValueError(f"Missing required measured evidence: {missing}")
    values = frame[required].apply(pd.to_numeric, errors="coerce")
    return pd.Series(np.isfinite(values.to_numpy()).all(axis=1), index=frame.index)


def cell_boundary_evidence(row: dict, native) -> dict:
    path = native.clean_text(row.get("cell_mask_path"))
    label = native.as_int(row.get("mask_label_id"), default=0)
    if not path or label <= 0:
        return {"cell_edge_touch": None, "cell_mask_label_present": False}
    array = native.load_label_array(path)
    if array.ndim != 2 or not np.any(array == label):
        return {"cell_edge_touch": None, "cell_mask_label_present": False}
    edge = any(
        np.any(border == label) for border in (array[0], array[-1], array[:, 0], array[:, -1])
    )
    return {"cell_edge_touch": bool(edge), "cell_mask_label_present": True}


def attach_pixel_audit(frame: pd.DataFrame, inventory: pd.DataFrame, audit: Path) -> pd.DataFrame:
    receipt = json.loads((audit / "receipt.json").read_text())
    validation = json.loads((audit / "validation.json").read_text())
    if receipt["exit_status"] != 0 or receipt["changed_inputs"] or validation["exit_status"] != 0:
        raise ValueError("Pixel audit execution or source-integrity validation failed")
    for name in ("command_receipt.json", "validation_command_receipt.json"):
        if json.loads((audit / name).read_text())["exit_status"] != 0:
            raise ValueError("Pixel audit lacks a successful command receipt")
    if not all(check["pass"] for check in validation["checks"]):
        raise ValueError("Pixel audit validation has a failed acceptance check")
    pixel = pd.read_csv(audit / "objects.csv")
    pixel = pixel.loc[pixel["role"].eq("standard")].copy()
    if (
        len(pixel) != 47
        or pixel["row_key"].duplicated().any()
        or set(pixel["row_key"]) != set(inventory["review_key"])
    ):
        raise ValueError("Pixel audit differs from the original standard census")
    require_hash(audit / "objects.csv", receipt["artifacts"]["objects.csv"]["sha256"])
    details = {
        row["row_key"]: row
        for row in json.loads((audit / "objects.json").read_text())
        if row["role"] == "standard"
    }
    originals = inventory.set_index("review_key")
    for row in pixel.to_dict("records"):
        original = originals.loc[row["row_key"]]
        detail = details[row["row_key"]]
        if (
            not row["area_exact_match"]
            or not row["iod_six_decimal_match"]
            or not row["helper_complete_label_coverage"]
        ):
            raise ValueError("Pixel audit lacks exact area, rounded IOD or full-label evidence")
        if (
            int(row["stored_area_px"]) != int(original["area_px"])
            or abs(row["stored_iod"] - float(original["iod"])) > 1e-12
        ):
            raise ValueError("Pixel measurement differs from original nucleus record")
        for field, audit_field in (("mask_path", "mask_path"), ("source_image_path", "raw_path")):
            if Path(original[field]).resolve() != Path(row[audit_field]).resolve():
                raise ValueError("Pixel audit source differs from the original nucleus")
        for field in ("filename", "specimen_id", "label", "tile_name"):
            if str(detail["metric_input_record"][field]) != str(original[field]):
                raise ValueError("Pixel-local view changed an original identity")
        local = detail["metric_input_record"]
        for axis in ("x", "y"):
            if (
                abs(
                    float(local[f"nuc_centroid_{axis}"])
                    - (float(original[f"centroid_{axis}"]) - float(original[f"tile_{axis}0"]))
                )
                > 1e-9
                or float(local[f"tile_{axis}0"]) != 0
            ):
                raise ValueError("Pixel-local coordinate view is inconsistent")
        for relative in detail["gallery_images"]:
            require_hash(audit / relative, receipt["artifacts"][relative]["sha256"])
    metrics = [
        name
        for name in pixel
        if name.startswith("iod_") and name not in {"iod_delta", "iod_six_decimal_match"}
    ]
    out = frame.merge(
        pixel[["row_key", *metrics]],
        left_on="review_key",
        right_on="row_key",
        how="left",
        validate="one_to_one",
    )
    out["source_nucleus_edge_touch"] = out["review_key"].map(
        {key: value["gallery_mask_touches_tile_edge"] for key, value in details.items()}
    )
    out["pixel_coordinate_mode"] = "standard_global_to_local_copy"
    return out


def strict_candidates(frame: pd.DataFrame, directory: Path, native) -> pd.DataFrame:
    out = frame.copy()
    out = out.loc[out["finite_required_evidence"]].copy()
    for name in native.BOOL_COLUMNS:
        if name not in out:
            raise ValueError(f"Missing required linked-pair flag: {name}")
        text = out[name].astype(str).str.lower()
        if not text.isin({"true", "false", "1", "0"}).all():
            raise ValueError(f"Invalid boolean linked-pair flag: {name}")
    path = directory / "candidate_input.csv.gz"
    out.to_csv(path, index=False)
    if out.empty:
        return out
    try:
        return native.prepare_candidate_frame(
            path, max_overlap_label_count=1, min_overlap_fraction=0.98
        )
    except ValueError as exc:
        if str(exc) == "No candidate rows passed strict balanced genome-IOD eligibility.":
            return out.iloc[0:0].copy()
        raise


def full_template_balance(template: pd.DataFrame, matched: pd.DataFrame, native) -> dict:
    if matched.empty:
        return {
            "n_reference": 40,
            "status": "insufficient_candidates",
            "max_abs_standardized_mean_difference": None,
            "max_ks_distance": None,
        }
    mean = max(
        abs(float(matched[f"z_{name}"].mean() - template[f"z_{name}"].mean()))
        for name in MATCH_FEATURES
    )
    ks = max(
        float(native.ks_2samp(matched[f"z_{name}"], template[f"z_{name}"]).statistic)
        for name in MATCH_FEATURES
    )
    return {
        "n_reference": 40,
        "status": "balanced" if mean <= 0.10 and ks <= 0.25 else "limited_overlap",
        "max_abs_standardized_mean_difference": mean,
        "max_ks_distance": ks,
    }


def project_groups(
    census: pd.DataFrame, candidates: pd.DataFrame, template: pd.DataFrame, scaling: dict, native
):
    proposed_pool = candidates.loc[~candidates["manual_decision"].eq("problem")].copy()
    eligible = (
        native.quality_match.add_matching_features(proposed_pool)
        if len(proposed_pool)
        else proposed_pool.copy()
    )
    if len(eligible):
        eligible = native.apply_frozen_scaling(eligible, scaling)
    groups = []
    selections = []
    for filename, specimen in STANDARDS.items():
        rows = census.loc[census["filename"].eq(filename)]
        pool = eligible.loc[eligible["filename"].eq(filename)]
        matched, subset = native.match_species_to_template(template, pool)
        subset["n_reference"] = len(matched)
        missing = sorted(
            set(range(1, 41)) - set(matched.get("quality_template_id", pd.Series(dtype=int)))
        )
        matched["accepted_for_standard_summary"] = matched["manual_decision"].eq("keep")
        accepted = matched.loc[matched["accepted_for_standard_summary"]]
        groups.append(
            {
                "standard_process_id": rows["standard_process_id"].iloc[0],
                "specimen_id": specimen,
                "input_count": len(rows),
                "linked_count": int(rows["has_cell_match"].sum()),
                "finite_evidence_count": int(rows["finite_required_evidence"].sum()),
                "strict_count": int(candidates["filename"].eq(filename).sum()),
                "finite_match_count": len(pool),
                "matched_count": len(matched),
                "manual_keep_count": int(rows["manual_decision"].eq("keep").sum()),
                "manual_problem_count": int(rows["manual_decision"].eq("problem").sum()),
                "manual_unsure_count": int(rows["manual_decision"].eq("unsure").sum()),
                "unreviewed_count": int(rows["manual_decision"].eq("").sum()),
                "accepted_matched_count": len(accepted),
                "accepted_standard_median_iod": float(accepted["nuc_iod"].median())
                if len(accepted)
                else None,
                "target_template_count": 40,
                "missing_template_ids": missing,
                "template_shortfall": len(missing),
                "assigned_subset_balance": {
                    key: value if not isinstance(value, float) or np.isfinite(value) else None
                    for key, value in subset.items()
                },
                "full_template_balance": full_template_balance(template, matched, native),
                "calibration_status": "unresolved; review projection does not validate DNA content",
            }
        )
        selections.append(matched)
    return pd.concat(selections, ignore_index=True), groups


def write_review(census: pd.DataFrame, out: Path, native, audit: Path) -> None:
    review = out / "review"
    review.mkdir()
    records = []
    for index, row in census.reset_index(drop=True).iterrows():
        record = native.process_row(
            row, index + 1, review, thumb_size=384, crop_margin=64, force=True
        )
        record["initial_decision"] = row["manual_decision"]
        records.append(record)
    native.write_index(
        review,
        records,
        panel="designated_standards_existing_masks",
        title="Designated standards: 413 / 32469 and 414 / 32470",
        storage_key="designated-standards-existing-masks-v1",
        default_decision=None,
        blind_target_metrics=True,
        common_support_only=False,
        auto_fill_target_per_species=None,
    )
    (review / "records.json").write_text(json.dumps(records, indent=2) + "\n")
    crop_dir = review / "nuclear_audit_crops"
    crop_dir.mkdir()
    audit_rows = {
        row["row_key"]: row
        for row in json.loads((audit / "objects.json").read_text())
        if row["role"] == "standard"
    }
    cards = []
    for row in census.to_dict("records"):
        source = audit_rows[row["review_key"]]
        panels = []
        for caption, relative in zip(
            ("Raw green channel", "Whole nucleus label", "Label outline"), source["gallery_images"]
        ):
            path = audit / relative
            shutil.copy2(path, crop_dir / path.name)
            panels.append(
                f'<figure><img src="nuclear_audit_crops/{html.escape(path.name)}" alt="{caption}"><figcaption>{caption}</figcaption></figure>'
            )
        edge = "yes" if source["gallery_mask_touches_tile_edge"] else "no"
        padding = "yes" if source["gallery_clipped"] else "no"
        cards.append(
            f"<article><h2>{html.escape(row['review_key'])}</h2><p>Specimen {html.escape(str(row['specimen_id']))}; source-label edge: {edge}; crop padding clipped: {padding}.</p><div>{''.join(panels)}</div></article>"
        )
    (review / "nuclear_audit.html").write_text(
        '<!doctype html><html lang="en"><meta charset="utf-8"><title>Full-label standard nucleus crops</title><style>body{font:16px sans-serif;margin:24px}article{border-top:1px solid #bbb;padding:16px 0}h2{font-size:15px}article div{display:flex;gap:12px;flex-wrap:wrap}figure{margin:0}img{max-width:320px;height:auto}</style><h1>Full-label standard nucleus crops</h1><p>47 unchanged candidate labels from the validated pixel audit. Source-label edges and clipped crop padding are distinct. No crop is an accepted standard measurement. Record manual decisions in the <a href="index.html">paired-mask review</a>.</p>'
        + "".join(cards)
        + "</html>"
    )
    index = review / "index.html"
    index.write_text(
        index.read_text().replace(
            "</header>",
            '<p style="padding:0 20px"><a href="nuclear_audit.html">Inspect full-label nucleus crops and source-edge flags</a>. Keep / problem / unsure are manual decisions; technical matching does not accept an object.</p></header>',
            1,
        )
    )


def execute(args: argparse.Namespace) -> dict:
    roots = {
        key: Path(getattr(args, f"{key}_root")).resolve() for key in ("paper", "cell", "archive")
    }
    pins = source_paths(roots)
    pins[PINS.resolve()] = digest(PINS)
    pins[Path(__file__).resolve()] = digest(Path(__file__).resolve())
    for name in ("pixi.toml", "pixi.lock"):
        path = roots["paper"] / name
        pins[path] = digest(path)
    for path in args.decisions_csv:
        pins[path.resolve()] = digest(path)
    out = fresh_destination(args.output_dir, [roots["cell"], roots["archive"], *pins])
    verify_sources(pins)
    inventory_path = roots["cell"] / STANDARD_RUN / "nucleus_measurements/all_measurements.csv"
    inventory = identity_frame(
        pd.read_csv(inventory_path, dtype=str, keep_default_na=False), inventory=True
    )
    validate_census(inventory, inventory)
    decisions = load_decisions(args.decisions_csv, inventory)
    out.mkdir()
    os.environ["MPLCONFIGDIR"] = str(out / "matplotlib_config")
    native = import_native(roots, pins)
    linkage = out / "linkage"
    command = [
        sys.executable,
        "-B",
        str(roots["cell"] / "scripts/build_cell_nucleus_linkage_report.py"),
        "--cell-csv",
        str(roots["cell"] / STANDARD_RUN / "cell_measurements/all_measurements.csv"),
        "--nucleus-csv",
        str(inventory_path),
        "--image-type",
        "pmount",
        "--window-radius",
        "1",
        "--output-dir",
        str(linkage),
    ]
    start = clock()
    env = {**os.environ, "PYTHONDONTWRITEBYTECODE": "1"}
    with (out / "linkage.log").open("w") as log:
        result = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, env=env, check=False)
    receipt = {
        "argv": command,
        "start_utc": start,
        "finish_utc": clock(),
        "exit_status": result.returncode,
    }
    (out / "linkage_command_receipt.json").write_text(json.dumps(receipt, indent=2) + "\n")
    if result.returncode:
        raise RuntimeError(
            f"Native linkage failed with exit {result.returncode}; see {out / 'linkage.log'}"
        )
    census = identity_frame(
        pd.read_csv(linkage / "linked_nucleus_pairs.csv.gz", dtype={"specimen_id": str})
    )
    validate_census(census, inventory)
    census["review_key"] = [
        native.review.review_key_for_row(row) for row in census.to_dict("records")
    ]
    boundaries = pd.DataFrame(
        [cell_boundary_evidence(row, native.review) for row in census.to_dict("records")]
    )
    census = pd.concat([census, boundaries], axis=1)
    shapes = pd.DataFrame(
        [native.review.compute_shape_qc_for_row(row) for row in census.to_dict("records")]
    )
    census = pd.concat([census, shapes], axis=1)
    census["manual_decision"] = census["review_key"].map(decisions).fillna("")
    audit = roots["paper"] / AUDIT_RELATIVE
    census = attach_pixel_audit(census, inventory, audit)
    census["finite_required_evidence"] = finite_evidence(census) & census["cell_mask_label_present"]
    candidates = strict_candidates(census, out, native.extension.source_qc)
    census["native_strict_eligible"] = census["review_key"].isin(candidates["review_key"])
    manifest = json.loads(
        (roots["paper"] / "analyses/02_morphology/provenance/iod_source_manifest.json").read_text()
    )
    scaling = manifest["feature_scaling"]
    panel = pd.read_csv(roots["archive"] / TEMPLATE_RELATIVE)
    template = validate_template(
        panel.loc[panel["species"].eq("D. apalachicolae")].copy(), scaling, native.extension
    )
    matched, summaries = project_groups(census, candidates, template, scaling, native.extension)
    census["proposed_technical_match"] = census["review_key"].isin(matched["review_key"])
    census["accepted_for_standard_summary"] = census["proposed_technical_match"] & census[
        "manual_decision"
    ].eq("keep")
    census["review_cohort"] = "designated_standard_existing_masks"
    census["panel"] = "designated_standards_existing_masks"
    census["selection_rank"] = np.arange(1, len(census) + 1)
    census["gate_reasons"] = census.apply(
        lambda row: ";".join(
            name
            for name, passed in {
                "no_cell_link": row["has_cell_match"],
                "nonphysical_pair": row["physical_pair_ok"],
                "nonunique_cell": row["one_to_one_cell"],
                "mask_pair_rejected": row["keep_mask_pair"],
                "native_strict_core_rejected": row["keep_strict_core"],
                "cell_edge": not row["cell_edge_touch"],
                "nonfinite_evidence": row["finite_required_evidence"],
                "suspect_shape": not row["shape_is_suspect"],
                "overlap_multiple_labels": row["nucleus_mask_overlap_label_count"] <= 1,
                "overlap_below_098": row["nucleus_mask_best_overlap_fraction"] >= 0.98,
                "manual_problem": row["manual_decision"] != "problem",
            }.items()
            if not passed
        ),
        axis=1,
    )
    census.to_csv(out / "standard_census.csv", index=False)
    matched.to_csv(out / "proposed_technical_matches.csv", index=False)
    template.to_csv(out / "frozen_template.csv", index=False)
    (out / "frozen_scaling.json").write_text(json.dumps(scaling, indent=2) + "\n")
    (out / "standard_summary.json").write_text(
        json.dumps(summaries, indent=2, allow_nan=False) + "\n"
    )
    write_review(census, out, native.extension.viewer, audit)
    verify_sources(pins)
    return {
        "standards": summaries,
        "source_sha256": {str(path): value for path, value in pins.items()},
        "source_preservation": "all pinned inputs unchanged",
        "torch_loaded": "torch" in sys.modules,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    for root in ("paper", "cell", "archive"):
        parser.add_argument(f"--{root}-root", required=True, type=Path)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--decisions-csv", action="append", type=Path, default=[])
    args = parser.parse_args()
    start = clock()
    try:
        result = execute(args)
    except (ValueError, FileNotFoundError) as exc:
        print(f"REFUSAL: {exc}", file=sys.stderr)
        return 2
    result.update(
        {
            "start_utc": start,
            "finish_utc": clock(),
            "exit_status": 0,
            "argv": sys.argv,
            "cwd": str(Path.cwd()),
            "python": {
                "executable": sys.executable,
                "prefix": sys.prefix,
                "version": platform.python_version(),
            },
            "environment": {
                name: os.environ.get(name)
                for name in (
                    "PYTHONNOUSERSITE",
                    "PYTHONHOME",
                    "PYTHONPATH",
                    "R_LIBS",
                    "R_LIBS_SITE",
                    "R_LIBS_USER",
                    "CONDA_PREFIX",
                    "PIXI_ENVIRONMENT_NAME",
                )
            },
            "packages": {
                name: {
                    "version": importlib.metadata.version(name),
                    "location": str(importlib.metadata.distribution(name).locate_file("")),
                }
                for name in (
                    "numpy",
                    "pandas",
                    "scipy",
                    "Pillow",
                    "tifffile",
                    "scikit-image",
                    "matplotlib",
                    "roifile",
                )
            },
        }
    )
    (args.output_dir / "receipt.json").write_text(
        json.dumps(result, indent=2, allow_nan=False) + "\n"
    )
    print(json.dumps({"output": str(args.output_dir), "standards": result["standards"]}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
