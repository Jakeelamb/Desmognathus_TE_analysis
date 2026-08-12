from __future__ import annotations

import importlib.util
import json
import math
import shutil
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
VALIDATOR_PATH = ROOT / "analyses/02_morphology/validate_release.py"
IOD_PROVENANCE_PATH = ROOT / "analyses/02_morphology/provenance/iod_source_manifest.json"
S22_PATH = ROOT / "analyses/02_morphology/data/relative_nuclear_iod_species.csv"
S01_PATH = ROOT / "data/identity/analysis_availability.csv"
CONDITIONAL_CALIBRATION_COLUMNS = (
    "iod_ratio_to_fuscus",
    "genome_size_pg_fuscus_anchored",
    "genome_size_pg_ci_low",
    "genome_size_pg_ci_high",
    "genome_size_reference_species",
    "genome_size_reference_pg",
    "genome_size_calibration_pg_per_iod",
    "genome_size_rank",
)


def load_validator():
    spec = importlib.util.spec_from_file_location("morphology_validate_release", VALIDATOR_PATH)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def copy_morphology_data(tmp_path: Path) -> None:
    destination = tmp_path / "analyses/02_morphology/data"
    destination.parent.mkdir(parents=True)
    shutil.copytree(ROOT / "analyses/02_morphology/data", destination)
    provenance_destination = tmp_path / "analyses/02_morphology/provenance"
    provenance_destination.mkdir(parents=True)
    shutil.copy2(
        IOD_PROVENANCE_PATH,
        provenance_destination / "iod_source_manifest.json",
    )
    identity_destination = tmp_path / "data/identity"
    identity_destination.mkdir(parents=True)
    for filename in [
        "analysis_availability.csv",
        "path24_panel.csv",
        "species_taxonomy_crosswalk.csv",
    ]:
        shutil.copy2(ROOT / "data/identity" / filename, identity_destination / filename)


def test_release_replay_matches_frozen_canonical_tables() -> None:
    errors = load_validator().validate_release(ROOT)
    assert errors == ()


def test_s22_conditional_calibration_replays_from_declared_assembly_anchor() -> None:
    validator = load_validator()
    released = pd.read_csv(S22_PATH)
    replayed = validator.reconstruct_s22_calibration(ROOT)
    selected_columns = ["species", *CONDITIONAL_CALIBRATION_COLUMNS]
    pd.testing.assert_frame_equal(
        released[selected_columns],
        replayed[selected_columns],
        check_exact=False,
        rtol=1e-10,
        atol=5.1e-9,
    )

    provenance = json.loads(IOD_PROVENANCE_PATH.read_text(encoding="utf-8"))
    anchor = provenance["conditional_genome_size_anchor"]
    expected_pg = anchor["reported_assembly_span_gbp"] / anchor["gbp_per_pg"]
    assert math.isclose(
        anchor["derived_anchor_pg_1c"],
        expected_pg,
        rel_tol=0,
        abs_tol=1e-12,
    )
    assert anchor["assembly_accession"] == "GCA_050004315.1"
    assert anchor["ploidy_semantics"] == "1C_equivalent"
    assert anchor["status"] == "conditional_assembly_derived_not_measured_c_value"
    assert anchor["reference_standard_processes"] == [
        "Process_413_green.ome.tiff",
        "Process_414_green.ome.tiff",
    ]
    assert anchor["reference_standard_specimen_ids"] == ["32469", "32470"]
    assert anchor["author_confirmed_same_staining_runs"] is True
    assert anchor["anchor_uncertainty_included"] is False


def test_s22_replay_rejects_tampered_assembly_anchor_arithmetic(tmp_path: Path) -> None:
    copy_morphology_data(tmp_path)
    path = tmp_path / "analyses/02_morphology/provenance/iod_source_manifest.json"
    provenance = json.loads(path.read_text(encoding="utf-8"))
    provenance["conditional_genome_size_anchor"]["derived_anchor_pg_1c"] += 0.25
    path.write_text(json.dumps(provenance, indent=2) + "\n", encoding="utf-8")

    errors = load_validator().validate_release(tmp_path)
    assert any(error.startswith("[S22_CALIBRATION_PROVENANCE]") for error in errors)


def test_s22_replay_rejects_tampered_conditional_estimate(tmp_path: Path) -> None:
    copy_morphology_data(tmp_path)
    path = tmp_path / "analyses/02_morphology/data/relative_nuclear_iod_species.csv"
    released = pd.read_csv(path)
    released.loc[
        released["species"].eq("D. valtos"),
        "genome_size_pg_fuscus_anchored",
    ] += 1.0
    released.to_csv(path, index=False)

    errors = load_validator().validate_release(tmp_path)
    assert any(error.startswith("[S22_GENOME_SIZE_PG_FUSCUS_ANCHORED]") for error in errors)


def test_s01_conditional_scale_replays_as_an_exact_s22_projection() -> None:
    validator = load_validator()
    released = pd.read_csv(S22_PATH)
    projected = validator.project_s22_availability(ROOT, released).set_index("species")
    stored = pd.read_csv(S01_PATH).set_index("species")
    column = "step06_genome_size_pg_fuscus_anchored"
    pd.testing.assert_series_equal(
        stored[column],
        projected[column],
        check_exact=False,
        check_names=False,
        rtol=1e-10,
        atol=5.1e-9,
    )


def test_s01_projection_rejects_tampered_conditional_scale(tmp_path: Path) -> None:
    copy_morphology_data(tmp_path)
    path = tmp_path / "data/identity/analysis_availability.csv"
    availability = pd.read_csv(path)
    availability.loc[
        availability["species"].eq("fuscus"),
        "step06_genome_size_pg_fuscus_anchored",
    ] += 1.0
    availability.to_csv(path, index=False)

    errors = load_validator().validate_release(tmp_path)
    assert any(error.startswith("[S01_S22_CONDITIONAL_SCALE]") for error in errors)


def test_s19_replay_rejects_changed_biological_support(tmp_path: Path) -> None:
    copy_morphology_data(tmp_path)
    path = tmp_path / "analyses/02_morphology/data/cell_nucleus_species_estimates.csv"
    released = pd.read_csv(path)
    released.loc[0, "largest_specimen_fraction"] = 0.0
    released.to_csv(path, index=False)

    errors = load_validator().validate_release(tmp_path)
    assert any(error.startswith("[S19_LARGEST_SPECIMEN_FRACTION]") for error in errors)


def test_s23_replay_rejects_changed_image_summary(tmp_path: Path) -> None:
    copy_morphology_data(tmp_path)
    path = tmp_path / "analyses/02_morphology/data/nuclear_iod_by_image.csv"
    released = pd.read_csv(path)
    released.loc[0, "image_median_iod"] += 1
    released.to_csv(path, index=False)

    errors = load_validator().validate_release(tmp_path)
    assert any(error.startswith("[S23_IMAGE_MEDIAN_IOD]") for error in errors)


def test_s22_replay_rejects_changed_path24_anchor(tmp_path: Path) -> None:
    copy_morphology_data(tmp_path)
    path = tmp_path / "analyses/02_morphology/data/relative_nuclear_iod_species.csv"
    released = pd.read_csv(path)
    released.loc[0, "relative_iod_anchor"] *= 1.01
    released.to_csv(path, index=False)

    errors = load_validator().validate_release(tmp_path)
    assert any(error.startswith("[S22_RELATIVE_IOD_ANCHOR]") for error in errors)


def test_s22_replay_rejects_invalid_bootstrap_boundaries(tmp_path: Path) -> None:
    copy_morphology_data(tmp_path)
    path = tmp_path / "analyses/02_morphology/data/relative_nuclear_iod_species.csv"
    released = pd.read_csv(path)
    released.loc[0, "iod_equal_image_ci_low"] = released.loc[0, "iod_equal_image_estimate"] + 1
    released.to_csv(path, index=False)

    errors = load_validator().validate_release(tmp_path)
    assert any(error.startswith("[S22_BOOTSTRAP_BOUNDS]") for error in errors)


def test_pixel_scale_replay_rejects_changed_physical_area(tmp_path: Path) -> None:
    copy_morphology_data(tmp_path)
    path = tmp_path / "analyses/02_morphology/data/cell_nucleus_objects.csv"
    released = pd.read_csv(path, low_memory=False)
    released.loc[0, "cell_area_um2"] += 0.0144
    released.to_csv(path, index=False)

    errors = load_validator().validate_release(tmp_path)
    assert any(error.startswith("[PIXEL_SCALE_S18_CELL]") for error in errors)


def test_s24_replay_rejects_changed_balance_decision(tmp_path: Path) -> None:
    copy_morphology_data(tmp_path)
    path = tmp_path / "analyses/02_morphology/data/nuclear_iod_quality_balance.csv"
    released = pd.read_csv(path)
    released.loc[
        released["species"].eq("D. amphileucus"),
        "passes_prespecified_balance_thresholds",
    ] = False
    released.to_csv(path, index=False)

    errors = load_validator().validate_release(tmp_path)
    assert any(error.startswith("[S24_BALANCE_BOOLEAN]") for error in errors)


def test_release_replay_rejects_changed_animal_support(tmp_path: Path) -> None:
    copy_morphology_data(tmp_path)
    path = tmp_path / "analyses/02_morphology/data/cell_nucleus_objects.csv"
    released = pd.read_csv(path, low_memory=False)
    released.loc[0, "specimen_id"] = 999_999_999
    released.to_csv(path, index=False)

    errors = load_validator().validate_release(tmp_path)
    assert any(error.startswith("[MORPHOLOGY_ANIMALS]") for error in errors)
