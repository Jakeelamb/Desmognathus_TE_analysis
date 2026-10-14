"""Portable boundary checks and explicitly qualified native recovery controls."""

from __future__ import annotations

import importlib.util
import json
import math
import os
import sys
import zipfile
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
PATH = ROOT / "scripts/microscopy/recovery_backend.py"
SPEC = importlib.util.spec_from_file_location("recovery_backend", PATH)
backend = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = backend
SPEC.loader.exec_module(backend)


def dummy_tools(tmp_path):
    return backend.Toolchain(
        tmp_path / "missing.json",
        tmp_path / "java",
        tmp_path / "ij.jar",
        tmp_path / "classes",
        tmp_path / "reducer",
        "fixture-only",
    )


@pytest.mark.parametrize(
    "pixels",
    [
        np.array([], dtype=np.uint8),
        np.array([1], dtype=np.uint16),
        np.ones((1, 1), dtype=np.uint8),
        np.array([1.0]),
        [1],
    ],
)
def test_reducer_boundary_rejects_pixel_shape_dtype(pixels, tmp_path):
    with pytest.raises(ValueError, match="uint8 vector"):
        backend.measure_pixels(pixels, 255, dummy_tools(tmp_path))


@pytest.mark.parametrize("background", [0, 256, -1, math.nan, math.inf, True, "255", 1 + 2j])
def test_reducer_boundary_rejects_background(background, tmp_path):
    with pytest.raises(ValueError, match="Background"):
        backend.measure_pixels(np.array([1], dtype=np.uint8), background, dummy_tools(tmp_path))


@pytest.mark.parametrize("dims", [(0, 1), (1, -1), (True, 1), (1.0, 1), (10001, 10000)])
def test_raster_boundary_rejects_dimensions(dims, tmp_path):
    with pytest.raises(ValueError, match="dimensions"):
        backend.rasterize_rois(
            tmp_path / "none.zip", *dims, tmp_path / "out", dummy_tools(tmp_path)
        )


@pytest.mark.parametrize(
    "names",
    [
        ["../nucleus.roi"],
        ["nested/nucleus.roi"],
        ["other.roi"],
        ["background_001.roi"],
        ["nucleus.roi", "nucleus.roi"],
        ["nucleus.roi", "background_000.roi"],
        ["NUCLEUS.roi"],
    ],
)
def test_zip_role_preflight_refuses_unsafe_or_missing_nucleus(names, tmp_path):
    archive = tmp_path / "input.zip"
    with zipfile.ZipFile(archive, "w") as out:
        for name in names:
            out.writestr(name, b"Iout" + b"\0" * 60)
    destination = tmp_path / "out"
    before = archive.read_bytes()
    with pytest.raises(ValueError):
        backend.rasterize_rois(archive, 16, 12, destination, dummy_tools(tmp_path))
    assert not destination.exists()
    assert archive.read_bytes() == before


def test_zip_malformed_oversized_nonregular_and_existing_destination(tmp_path):
    archive = tmp_path / "input.zip"
    archive.write_bytes(b"not a zip")
    with pytest.raises(ValueError, match="Malformed"):
        backend.rasterize_rois(archive, 16, 12, tmp_path / "out", dummy_tools(tmp_path))
    with zipfile.ZipFile(archive, "w") as out:
        out.writestr("nucleus.roi", b"x" * (1048576 + 1))
    with pytest.raises(ValueError, match="oversized"):
        backend.rasterize_rois(archive, 16, 12, tmp_path / "out", dummy_tools(tmp_path))
    entry = zipfile.ZipInfo("nucleus.roi")
    entry.external_attr = 0o120777 << 16
    with zipfile.ZipFile(archive, "w") as out:
        out.writestr(entry, b"Iout" + b"\0" * 60)
    with pytest.raises(ValueError, match="Non-regular"):
        backend.rasterize_rois(archive, 16, 12, tmp_path / "out", dummy_tools(tmp_path))
    destination = tmp_path / "out"
    destination.mkdir()
    (destination / "sentinel").write_text("preserve")
    with pytest.raises(ValueError, match="Fresh destination"):
        backend.rasterize_rois(archive, 16, 12, destination, dummy_tools(tmp_path))
    assert (destination / "sentinel").read_text() == "preserve"


def test_zip_symlink_and_input_ancestry_refused(tmp_path):
    archive = tmp_path / "input.zip"
    archive.write_bytes(b"placeholder")
    link = tmp_path / "linked.zip"
    link.symlink_to(archive)
    with pytest.raises(ValueError, match="Symlink"):
        backend.rasterize_rois(link, 16, 12, tmp_path / "out", dummy_tools(tmp_path))
    with pytest.raises(ValueError, match="Fresh destination"):
        backend.rasterize_rois(archive, 16, 12, tmp_path, dummy_tools(tmp_path))


def test_destination_parent_segments_cannot_escape_project_scratch(tmp_path):
    archive = tmp_path / "input.zip"
    archive.write_bytes(b"placeholder")
    disguised = ROOT / ".pixi/runs/../../scripts/recovery-refusal-fixture"
    with pytest.raises(ValueError, match="source paths are protected"):
        backend.rasterize_rois(archive, 16, 12, disguised, dummy_tools(tmp_path))
    assert not (ROOT / "scripts/recovery-refusal-fixture").exists()


def test_receipt_directory_symlink_cannot_redirect_writes(tmp_path):
    protected = tmp_path / "protected"
    protected.mkdir()
    (protected / "sentinel").write_text("preserve")
    evidence = tmp_path / "evidence"
    evidence.symlink_to(protected, target_is_directory=True)
    with pytest.raises(ValueError, match="Symlink"):
        backend.retained_run(["unexecuted-fixture-command"], evidence, "fixture")
    assert sorted(p.name for p in protected.iterdir()) == ["sentinel"]
    assert (protected / "sentinel").read_text() == "preserve"


def test_load_is_read_only_for_missing_or_changed_cache(tmp_path, monkeypatch):
    inputs = {"schema": 1, "tools": {"java": "/fixture-only/java"}, "fixture_only": True}
    monkeypatch.setattr(backend, "tool_inputs", lambda root: inputs)
    with pytest.raises(ValueError, match="missing or changed"):
        backend.load_toolchain(tmp_path)
    assert not (tmp_path / ".pixi").exists()
    identity = backend.tool_identity(inputs)
    cache = tmp_path / ".pixi/recovery-tools" / identity
    (cache / "classes").mkdir(parents=True)
    for name in ["ij.jar", "microscopy_recovery", "classes/RecoveryRoiBridge.class"]:
        (cache / name).write_bytes(b"fixture-only non-native artifact")
    artifacts = {
        name: backend.sha256(cache / name)
        for name in ["ij.jar", "microscopy_recovery", "classes/RecoveryRoiBridge.class"]
    }
    manifest = {
        "status": "complete",
        "inputs": inputs,
        "identity": identity,
        "artifacts": artifacts,
    }
    (cache / "manifest.json").write_text(json.dumps(manifest))
    (cache / "microscopy_recovery").write_bytes(b"changed")
    before = {str(p): p.read_bytes() for p in cache.rglob("*") if p.is_file()}
    with pytest.raises(ValueError, match="artifact changed"):
        backend.load_toolchain(tmp_path)
    assert before == {str(p): p.read_bytes() for p in cache.rglob("*") if p.is_file()}
    (cache / "microscopy_recovery").write_bytes(b"fixture-only non-native artifact")
    monkeypatch.setattr(backend, "IMAGEJ_SHA256", artifacts["ij.jar"])
    with pytest.raises(ValueError, match="qualification missing or changed"):
        backend.load_toolchain(tmp_path)
    qualification = cache / "qualification.json"
    qualification.write_text(json.dumps({"status": "passed", "fixture_only": True}))
    manifest["qualification_sha256"] = backend.sha256(qualification)
    (cache / "manifest.json").write_text(json.dumps(manifest))
    assert backend.load_toolchain(tmp_path).identity == identity
    qualification.write_text(json.dumps({"status": "failed", "fixture_only": True}))
    with pytest.raises(ValueError, match="qualification missing or changed"):
        backend.load_toolchain(tmp_path)


@pytest.fixture
def native_tools():
    if os.environ.get("PIXI_ENVIRONMENT_NAME", os.environ.get("PIXI_ENV_NAME")) != "recovery":
        pytest.skip(
            "native backend NOT EXECUTED: run locked recovery environment after explicit recovery-setup"
        )
    return backend.load_toolchain(ROOT)


def test_native_signed_floor_and_compensated_reduction(native_tools):
    pixels = np.tile(np.array([0, 1, 218, 255, 100, 250], dtype=np.uint8), 10000)
    before = pixels.copy()
    actual = backend.measure_pixels(pixels, 218, native_tools)
    reference = math.fsum(math.log10(218 / max(int(pixel), 1)) for pixel in pixels)
    assert actual["count"] == pixels.size
    assert actual["raw_sum"] == int(pixels.astype(np.uint64).sum())
    assert actual["min_pixel"] == 0 and actual["max_pixel"] == 255
    assert actual["zero_pixels"] == actual["one_pixels"] == actual["saturated_pixels"] == 10000
    assert actual["iod"] == pytest.approx(reference, rel=0, abs=1e-10)
    assert actual["mean_od"] == pytest.approx(reference / pixels.size, rel=0, abs=1e-14)
    assert actual["median_pixel"] == 159
    assert np.array_equal(pixels, before)
    negative = backend.measure_pixels(np.array([255], dtype=np.uint8), 218, native_tools)
    assert negative["iod"] < 0
    assert negative["median_pixel"] == 255
    assert backend.measure_pixels(np.array([0, 1], dtype=np.uint8), 1, native_tools)["iod"] == 0
    assert (
        backend.measure_pixels(np.array([0, 255], dtype=np.uint8), 255, native_tools)[
            "median_pixel"
        ]
        == 127.5
    )
    assert (
        backend.measure_pixels(np.array([1, 1, 255], dtype=np.uint8), 255, native_tools)[
            "median_pixel"
        ]
        == 1
    )


@pytest.mark.parametrize(
    "name,count,raw_sum",
    [
        ("asymmetric_integer", 67, 414876),
        ("concave_freehand", 73, 523940),
        ("subpixel_polygon", 65, 439827),
        ("composite_hole", 118, 887123),
        ("boundary_rectangle", 48, 295152),
    ],
)
def test_native_genuine_imagej_geometry_and_source_preservation(
    name, count, raw_sum, native_tools, tmp_path
):
    from PIL import Image

    fixture = native_tools.manifest_path.parent / "native_fixture"
    archive = fixture / f"{name}.zip"
    before = archive.read_bytes()
    paths = backend.rasterize_rois(archive, 16, 12, tmp_path / "out", native_tools)
    assert set(paths) == {"nucleus"}
    with Image.open(paths["nucleus"]) as image:
        assert image.mode == "L"
        actual = np.asarray(image)
    with Image.open(fixture / f"{name}.png") as image:
        reference = np.asarray(image)
    assert actual.shape == (12, 16) and actual.dtype == np.uint8
    assert np.array_equal(actual, reference)
    assert set(np.unique(actual)) == {0, 255}
    assert np.count_nonzero(actual) == count
    pixels = np.fromfunction(lambda y, x: 13 + x * 237 + y * 1051, (12, 16), dtype=int)
    assert int(pixels[actual != 0].sum()) == raw_sum
    assert archive.read_bytes() == before
    assert (tmp_path / "out/rasterization.json").is_file()
    if name == "composite_hole":
        assert not actual[4:7, 5:9].any()


@pytest.mark.parametrize(
    "name",
    [
        "line",
        "point",
        "z",
        "stack",
        "rotated",
        "rounded",
        "spline",
        "out_of_bounds",
        "subpixel_rectangle",
        "ellipse",
    ],
)
def test_native_unsupported_geometry_rejected_with_retained_staging(name, native_tools, tmp_path):
    fixture = native_tools.manifest_path.parent / "native_fixture"
    archive = fixture / f"reject_{name}.zip"
    before = archive.read_bytes()
    out = tmp_path / "out"
    with pytest.raises(ValueError, match="ImageJ ROI rejection"):
        backend.rasterize_rois(archive, 16, 12, out, native_tools)
    assert archive.read_bytes() == before
    assert out.exists() and not (out / "rasterization.json").exists()
    receipts = list((out / "evidence").glob("*.receipt.json"))
    assert len(receipts) == 1 and json.loads(receipts[0].read_text())["exit_code"] == 1


def test_native_entry_identity_owned_by_zip_and_backgrounds(native_tools, tmp_path):
    fixture = native_tools.manifest_path.parent / "native_fixture"
    with zipfile.ZipFile(fixture / "boundary_rectangle.zip") as archive:
        encoded = archive.read("nucleus.roi")
    path = tmp_path / "roles.zip"
    with zipfile.ZipFile(path, "w") as archive:
        archive.writestr("nucleus.roi", encoded)
        archive.writestr("background_001.roi", encoded)
    paths = backend.rasterize_rois(path, 16, 12, tmp_path / "out", native_tools)
    assert set(paths) == {"nucleus", "background_001"}
    assert paths["nucleus"].read_bytes() == paths["background_001"].read_bytes()
