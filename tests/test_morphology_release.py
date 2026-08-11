from __future__ import annotations

import importlib.util
import shutil
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
VALIDATOR_PATH = ROOT / "analyses/02_morphology/validate_release.py"


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
    identity_destination = tmp_path / "data/identity"
    identity_destination.mkdir(parents=True)
    for filename in ["path24_panel.csv", "species_taxonomy_crosswalk.csv"]:
        shutil.copy2(ROOT / "data/identity" / filename, identity_destination / filename)


def test_release_replay_matches_frozen_canonical_tables() -> None:
    errors = load_validator().validate_release(ROOT)
    assert errors == ()


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
    released.loc[0, "iod_equal_image_ci_low"] = (
        released.loc[0, "iod_equal_image_estimate"] + 1
    )
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
