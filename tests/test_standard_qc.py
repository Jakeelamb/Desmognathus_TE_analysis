"""Boundary controls and genuine frozen-helper acceptance for standard review."""

from __future__ import annotations

import importlib.util
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

PATH = Path(__file__).resolve().parents[1] / "scripts/microscopy/standard_qc.py"
SPEC = importlib.util.spec_from_file_location("standard_qc", PATH)
qc = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(qc)


def inventory_rows() -> pd.DataFrame:
    rows = []
    for filename, count in qc.COUNTS.items():
        for label in range(1, count + 1):
            rows.append(
                {
                    "filename": filename,
                    "specimen_id": qc.STANDARDS[filename],
                    "species": "D. fuscus",
                    "image_type": "pmount",
                    "tile_name": "tile_y000000_x000000.tiff",
                    "label": label,
                }
            )
    return pd.DataFrame(rows)


@pytest.mark.parametrize(
    "field,value",
    [
        ("filename", "Process_314_raw_green.ome.tiff"),
        ("specimen_id", "32470"),
        ("species", "D. apalachicolae"),
        ("image_type", "brightfield"),
        ("label", float("nan")),
        ("label", 1.5),
    ],
)
def test_identity_rejects_wrong_standard_or_invalid_label(field, value):
    frame = inventory_rows()
    frame[field] = frame[field].astype(object)
    frame.loc[0, field] = value
    with pytest.raises(ValueError):
        qc.identity_frame(frame, inventory=True)


def test_duplicate_and_changed_object_identity_are_rejected():
    frame = inventory_rows()
    frame.loc[1, "label"] = 1
    with pytest.raises(ValueError, match="Duplicate nucleus"):
        qc.identity_frame(frame, inventory=True)
    frame = qc.identity_frame(inventory_rows(), inventory=True).rename(
        columns={"label": "nucleus_label"}
    )
    frame.loc[0, "nucleus_object_id"] = "wrong"
    with pytest.raises(ValueError, match="differs from original"):
        qc.identity_frame(frame)


def test_missing_census_object_is_rejected():
    inventory = qc.identity_frame(inventory_rows(), inventory=True)
    with pytest.raises(ValueError, match="exact 47"):
        qc.validate_census(inventory.iloc[:-1], inventory)


def test_changed_source_is_rejected(tmp_path):
    path = tmp_path / "source"
    path.write_text("original")
    fingerprint = qc.digest(path)
    path.write_text("changed")
    with pytest.raises(ValueError, match="hash mismatch"):
        qc.require_hash(path, fingerprint)


def test_existing_destination_and_source_overlap_are_rejected(tmp_path):
    source = tmp_path / "original"
    source.mkdir()
    (source / "data").write_text("preserved")
    with pytest.raises(ValueError, match="already exists"):
        qc.fresh_destination(source, [])
    with pytest.raises(ValueError, match="ancestry overlap"):
        qc.fresh_destination(source / "new", [source])
    link = tmp_path / "link"
    link.symlink_to(source, target_is_directory=True)
    with pytest.raises(ValueError, match="Symlinked"):
        qc.fresh_destination(link / "new", [])
    assert (source / "data").read_text() == "preserved"


def test_manual_registry_does_not_accept_unreviewed(tmp_path):
    inventory = qc.identity_frame(inventory_rows(), inventory=True)
    assert qc.load_decisions([], inventory) == {}
    row = {"species": "D. fuscus", "review_key": inventory.loc[0, "review_key"], "decision": "keep"}
    path = tmp_path / "decisions.csv"
    pd.DataFrame([row]).to_csv(path, index=False)
    assert qc.load_decisions([path], inventory) == {row["review_key"]: "keep"}
    pd.DataFrame([{**row, "decision": "unreviewed"}]).to_csv(path, index=False)
    with pytest.raises(ValueError, match="keep, problem or unsure"):
        qc.load_decisions([path], inventory)
    pd.DataFrame([{**row, "review_key": "Process_314::target"}]).to_csv(path, index=False)
    with pytest.raises(ValueError, match="unknown standard"):
        qc.load_decisions([path], inventory)


def test_nonfinite_required_metrics_cannot_pass_evidence_gate():
    frame = pd.DataFrame(
        [
            {
                name: 1.0
                for name in [
                    *qc.SHAPE_METRICS,
                    *qc.PIXEL_METRICS,
                    "nuc_iod",
                    "nuc_area_um2",
                    "cell_area_um2",
                    "nucleus_mask_overlap_label_count",
                    "nucleus_mask_best_overlap_fraction",
                ]
            }
        ]
    )
    assert qc.finite_evidence(frame).tolist() == [True]
    for name in [
        *qc.SHAPE_METRICS,
        *qc.PIXEL_METRICS,
        "nucleus_mask_overlap_label_count",
        "nucleus_mask_best_overlap_fraction",
    ]:
        changed = frame.copy()
        changed.loc[0, name] = np.nan
        assert qc.finite_evidence(changed).tolist() == [False]


@pytest.fixture(scope="module")
def native():
    if os.environ.get("STANDARD_QC_NATIVE") != "1":
        pytest.skip(
            "Native acceptance requires the image-QC Pixi environment and pinned archival inputs"
        )
    roots = {
        key: Path(os.environ[f"STANDARD_QC_{key.upper()}_ROOT"]).resolve()
        for key in ("paper", "cell", "archive")
    }
    pins = qc.source_paths(roots)
    qc.verify_sources(pins)
    modules = qc.import_native(roots, pins)
    manifest = json.loads(
        (roots["paper"] / "analyses/02_morphology/provenance/iod_source_manifest.json").read_text()
    )
    panel = pd.read_csv(roots["archive"] / qc.TEMPLATE_RELATIVE)
    template = qc.validate_template(
        panel.loc[panel["species"].eq("D. apalachicolae")].copy(),
        manifest["feature_scaling"],
        modules.extension,
    )
    return modules.extension, template, manifest["feature_scaling"]


@pytest.mark.parametrize(
    "defect", ["duplicate_id", "missing_id", "nan_feature", "wrong_z", "zero_scale", "nan_scale"]
)
def test_native_template_corruption_is_rejected(native, defect):
    module, original, original_scaling = native
    template = original.copy()
    scaling = json.loads(json.dumps(original_scaling))
    if defect == "duplicate_id":
        template.loc[1, "quality_template_id"] = 1
    elif defect == "missing_id":
        template = template.iloc[:-1]
    elif defect == "nan_feature":
        template.loc[0, qc.MATCH_FEATURES[0]] = np.nan
    elif defect == "wrong_z":
        template.loc[0, f"z_{qc.MATCH_FEATURES[0]}"] += 0.1
    elif defect == "zero_scale":
        scaling[qc.MATCH_FEATURES[0]]["scale_iqr"] = 0
    else:
        scaling[qc.MATCH_FEATURES[0]]["scale_iqr"] = np.nan
    with pytest.raises(ValueError):
        qc.validate_template(template, scaling, module)


def projection_fixture(template: pd.DataFrame):
    census = qc.identity_frame(inventory_rows(), inventory=True)
    census["has_cell_match"] = True
    census["finite_required_evidence"] = True
    census["manual_decision"] = ""
    census["nuc_iod"] = 100.0
    selected = census.loc[census["standard_process_id"].eq("414")].copy().reset_index(drop=True)
    selected["iod_edge_sharpness"] = np.expm1(template.iloc[:16][qc.MATCH_FEATURES[0]].to_numpy())
    selected["iod_ring_intensity_median"] = 1.0
    selected["iod_ring_noise_iqr"] = np.expm1(template.iloc[:16][qc.MATCH_FEATURES[1]].to_numpy())
    return census, selected


def test_genuine_native_zero_and_16_candidate_shortfalls_and_manual_acceptance(native):
    module, template, scaling = native
    census, candidates = projection_fixture(template)
    matched, groups = qc.project_groups(census, candidates, template, scaling, module)
    zero, sixteen = groups
    assert zero["matched_count"] == 0
    assert zero["template_shortfall"] == 40
    assert zero["assigned_subset_balance"]["quality_match_status"] == "insufficient_candidates"
    assert zero["full_template_balance"]["status"] == "insufficient_candidates"
    assert sixteen["matched_count"] == 16
    assert sixteen["template_shortfall"] == 24
    assert len(sixteen["missing_template_ids"]) == 24
    assert sixteen["assigned_subset_balance"]["quality_match_status"] == "common_support"
    assert sixteen["full_template_balance"]["n_reference"] == 40
    assert (
        sixteen["full_template_balance"]["max_abs_standardized_mean_difference"]
        > sixteen["assigned_subset_balance"]["max_abs_standardized_mean_difference"]
    )
    assert matched["accepted_for_standard_summary"].sum() == 0
    assert sixteen["accepted_matched_count"] == 0
    assert sixteen["accepted_standard_median_iod"] is None
    candidates.loc[0, "manual_decision"] = "keep"
    census.loc[census["review_key"].eq(candidates.loc[0, "review_key"]), "manual_decision"] = "keep"
    kept, summaries = qc.project_groups(census, candidates, template, scaling, module)
    assert kept["accepted_for_standard_summary"].sum() == 1
    assert summaries[1]["accepted_matched_count"] == 1
    assert summaries[1]["accepted_standard_median_iod"] == 100.0


def test_nonfinite_evidence_never_enters_native_strict_candidates(native, tmp_path):
    module, _template, _scaling = native
    frame = pd.DataFrame(
        [
            {
                **{name: True for name in module.source_qc.BOOL_COLUMNS},
                "cell_edge_touch": False,
                "finite_required_evidence": False,
                "manual_decision": "",
                "species": "D. fuscus",
                "review_key": "candidate",
                "filename": "Process_414_green.ome.tiff",
            }
        ]
    )
    assert qc.strict_candidates(frame, tmp_path, module.source_qc).empty


def test_genuine_native_mask_loading_measures_cell_edge_and_missing_label(native, tmp_path):
    import sys

    import tifffile

    review = sys.modules["run_linked_pair_review"]
    mask = np.zeros((12, 12), dtype=np.uint32)
    mask[5:8, 4:7] = 1
    mask[0, 2:5] = 2
    path = tmp_path / "labels.tiff"
    tifffile.imwrite(path, mask)
    row = {"cell_mask_path": str(path), "mask_label_id": 1}
    assert qc.cell_boundary_evidence(row, review) == {
        "cell_edge_touch": False,
        "cell_mask_label_present": True,
    }
    assert qc.cell_boundary_evidence({**row, "mask_label_id": 2}, review) == {
        "cell_edge_touch": True,
        "cell_mask_label_present": True,
    }
    assert qc.cell_boundary_evidence({**row, "mask_label_id": 999}, review) == {
        "cell_edge_touch": None,
        "cell_mask_label_present": False,
    }


def test_manual_problem_excluded_from_match_but_keeps_strict_count(native):
    module, template, scaling = native
    census, candidates = projection_fixture(template)
    candidates.loc[0, "manual_decision"] = "problem"
    census.loc[census["review_key"].eq(candidates.loc[0, "review_key"]), "manual_decision"] = (
        "problem"
    )
    matched, summaries = qc.project_groups(census, candidates, template, scaling, module)
    assert summaries[1]["strict_count"] == 16
    assert summaries[1]["finite_match_count"] == 15
    assert summaries[1]["manual_problem_count"] == 1
    assert candidates.loc[0, "review_key"] not in set(matched["review_key"])


def test_invalid_linked_boolean_is_rejected_before_native_eligibility(native, tmp_path):
    module, _template, _scaling = native
    frame = pd.DataFrame(
        [
            {
                **{name: True for name in module.source_qc.BOOL_COLUMNS},
                "has_cell_match": "possibly",
                "finite_required_evidence": True,
            }
        ]
    )
    with pytest.raises(ValueError, match="Invalid boolean"):
        qc.strict_candidates(frame, tmp_path, module.source_qc)
    assert not (tmp_path / "candidate_input.csv.gz").exists()
