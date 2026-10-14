"""Fresh-output controls and genuine Bio-Formats unchanged crop handoff."""

from __future__ import annotations

import importlib.util
import json
import os
import sys
import uuid
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "source_reader", ROOT / "scripts/microscopy/source_reader.py"
)
reader = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(reader)


@pytest.mark.parametrize(
    "bounds",
    [
        (-1, 0, 1, 1),
        (0, -1, 1, 1),
        (0, 0, 0, 1),
        (0, 0, 1, -1),
        (0, 0, 1001, 1000),
        (True, 0, 1, 1),
        (0, 0, 1.5, 1),
        ("0", 0, 1, 1),
    ],
)
def test_refuse_invalid_crop_bounds(bounds):
    with pytest.raises(ValueError):
        reader.crop_bounds(*bounds)


def test_fresh_output_and_source_tree_refusal(tmp_path):
    root = tmp_path / "repo"
    source = tmp_path / "archive/source.vsi"
    source.parent.mkdir()
    source.write_bytes(b"fixture")
    output = root / ".pixi/runs/new"
    assert reader.fresh_run(root, output, source) == output
    with pytest.raises(FileExistsError):
        reader.fresh_run(root, output, source)
    for path in (tmp_path / "outside", source.parent / "output", root / ".pixi/runs"):
        with pytest.raises(ValueError):
            reader.fresh_run(root, path, source)
    alias = root / ".pixi/runs/alias"
    alias.symlink_to(source.parent, target_is_directory=True)
    with pytest.raises(ValueError, match="Symlink"):
        reader.fresh_run(root, alias / "output", source)


def test_source_under_runs_still_refuses_its_companion_tree(tmp_path):
    source = tmp_path / ".pixi/runs/archive/source.vsi"
    source.parent.mkdir(parents=True)
    source.write_bytes(b"fixture")
    with pytest.raises(ValueError, match="companion"):
        reader.fresh_run(tmp_path, source.parent / "companion/output", source)


def test_missing_optional_jar_keeps_failed_receipt(tmp_path, monkeypatch):
    source = tmp_path / "archive/source.vsi"
    source.parent.mkdir()
    source.write_bytes(b"fixture")
    java = tmp_path / "java"
    monkeypatch.setattr(reader, "java_path", lambda root: java)
    monkeypatch.setattr(reader, "tool_pins", lambda root, java: {"fixture": True})
    out = tmp_path / ".pixi/runs/missing-tools"
    with pytest.raises(ValueError, match="explicit"):
        reader.read_source(tmp_path, "metadata", source, out)
    receipt = json.loads((out / "receipt.json").read_text())
    assert receipt["status"] == "failed"
    assert receipt["start_utc"] and receipt["end_utc"]
    assert receipt["source"]["sha256"] == reader.sha256(source)
    assert receipt["companion_ets_hashes_computed"] is False
    assert "not declared read-only" in receipt["source"]["access_policy"]


def test_bad_seed_is_retained_and_not_repaired(tmp_path, monkeypatch):
    seed = tmp_path / "bad.jar"
    seed.write_bytes(b"not a JAR")
    monkeypatch.setattr(reader, "java_path", lambda root: tmp_path / "java")
    monkeypatch.setattr(reader, "tool_pins", lambda root, java: {"fixture": True})
    with pytest.raises(ValueError, match="checksum"):
        reader.setup(tmp_path, seed)
    receipt = reader.cache_path(tmp_path) / "setup.json"
    before = receipt.read_bytes()
    assert json.loads(before)["status"] == "failed"
    with pytest.raises(ValueError, match="in-place repair"):
        reader.setup(tmp_path, seed)
    assert receipt.read_bytes() == before


def test_failed_native_command_keeps_raw_streams(tmp_path):
    command = reader.run_java(
        [
            sys.executable,
            "-c",
            ("import sys; print('output'); print('error', file=sys.stderr); sys.exit(7)"),
        ],
        tmp_path,
    )
    assert command["exit_code"] == 7
    assert (tmp_path / "command.stdout.log").read_text() == "output\n"
    assert (tmp_path / "command.stderr.log").read_text() == "error\n"
    assert command["stdout_sha256"] == reader.sha256(tmp_path / "command.stdout.log")
    assert command["start_utc"] and command["end_utc"]


def test_timed_out_native_command_keeps_partial_logs(tmp_path, monkeypatch):
    monkeypatch.setattr(reader, "TIMEOUT", 0.1)
    command = reader.run_java(
        [sys.executable, "-c", "import time; print('started', flush=True); time.sleep(10)"],
        tmp_path,
    )
    assert command["exit_code"] is None
    assert "TimeoutExpired" in command["error"]
    assert (tmp_path / "command.stdout.log").read_text() == "started\n"
    assert command["start_utc"] and command["end_utc"]


@pytest.mark.parametrize("kind", ["uint16", "float", "rgb", "wrong_shape", "stack", "white"])
def test_reject_bad_crop_pixel_encoding_and_shape(tmp_path, kind):
    tifffile = pytest.importorskip("tifffile", reason="Requires image-qc or recovery")
    path = tmp_path / "bad.tiff"
    array = np.zeros((4, 5), dtype=np.uint8)
    kwargs = {"photometric": "minisblack"}
    if kind == "uint16":
        array = np.full((4, 5), 256, dtype=np.uint16)
    elif kind == "float":
        array = np.full((4, 5), np.nan, dtype=np.float32)
    elif kind == "rgb":
        array = np.zeros((4, 5, 3), dtype=np.uint8)
        kwargs = {"photometric": "rgb"}
    elif kind == "wrong_shape":
        array = np.zeros((5, 4), dtype=np.uint8)
    elif kind == "stack":
        array = np.zeros((2, 4, 5), dtype=np.uint8)
    elif kind == "white":
        kwargs = {"photometric": "miniswhite"}
    tifffile.imwrite(path, array, **kwargs)
    with pytest.raises(ValueError):
        reader.validate_crop(path, 5, 4)


def native_available():
    if os.environ.get("PIXI_ENVIRONMENT_NAME") != "recovery":
        pytest.skip("Native controls require locked recovery Java21")
    try:
        reader.java_path(ROOT)
        reader.load_jar(ROOT)
    except ValueError as error:
        pytest.skip(f"Explicit optional source-reader setup required: {error}")
    return pytest.importorskip("tifffile")


def test_genuine_public_crop_cli_exact_green_pixels_and_metadata():
    tifffile = native_available()
    evidence = ROOT / ".pixi/runs/source-reader-controls" / uuid.uuid4().hex
    evidence.mkdir(parents=True)
    y, x = np.indices((17, 19))
    rgb = np.zeros((17, 19, 3), dtype=np.uint8)
    rgb[..., 0] = 51
    rgb[..., 1] = (y * 19 + x) % 256
    rgb[..., 2] = 203
    # Deliberate 0/1/255 extremes detect inversion, clamp and scale.
    rgb[3, 4:7, 1] = [0, 1, 255]
    source = evidence / "fixture.tiff"
    tifffile.imwrite(source, rgb, photometric="rgb")
    before = reader.sha256(source)
    # Output must be outside the source directory even for synthetic controls.
    crop = evidence.with_name(evidence.name + "-crop")
    assert (
        reader.main(
            [
                "crop",
                "--source",
                str(source),
                "--out",
                str(crop),
                "--x",
                "4",
                "--y",
                "3",
                "--width",
                "8",
                "--height",
                "6",
            ]
        )
        == 0
    )
    pixels = reader.validate_crop(crop / "crop.tiff", 8, 6)
    np.testing.assert_array_equal(pixels, rgb[3:9, 4:12, 1])
    assert reader.sha256(source) == before
    receipt = json.loads((crop / "receipt.json").read_text())
    assert receipt["status"] == "complete" and receipt["command"]["exit_code"] == 0
    assert receipt["transforms"] == []
    assert "-autoscale" not in receipt["command"]["argv"]
    assert "-cache" not in receipt["command"]["argv"]
    assert "-no-upgrade" in receipt["command"]["argv"]
    metadata = evidence.with_name(evidence.name + "-metadata")
    assert reader.main(["metadata", "--source", str(source), "--out", str(metadata)]) == 0
    assert (metadata / "metadata.ome.xml").is_file()
    assert sorted(path.name for path in evidence.iterdir()) == ["fixture.tiff"]


def test_genuine_uint16_crop_fails_without_coercion():
    tifffile = native_available()
    evidence = ROOT / ".pixi/runs/source-reader-controls" / uuid.uuid4().hex
    evidence.mkdir(parents=True)
    source = evidence / "uint16.tiff"
    rgb = np.full((6, 7, 3), 1000, dtype=np.uint16)
    tifffile.imwrite(source, rgb, photometric="rgb")
    out = evidence.with_name(evidence.name + "-rejected")
    assert (
        reader.main(
            [
                "crop",
                "--source",
                str(source),
                "--out",
                str(out),
                "--x",
                "0",
                "--y",
                "0",
                "--width",
                "3",
                "--height",
                "2",
            ]
        )
        == 1
    )
    receipt = json.loads((out / "receipt.json").read_text())
    assert receipt["status"] == "failed"
    assert receipt["command"]["exit_code"] == 0
    assert "uint8" in receipt["error"]
    assert tifffile.imread(out / "crop.tiff").dtype == np.uint16


def test_genuine_pyramidal_rgb_crop_contains_only_full_resolution_green_pixels():
    tifffile = native_available()
    evidence = ROOT / ".pixi/runs/source-reader-controls" / uuid.uuid4().hex
    evidence.mkdir(parents=True)
    source = evidence / "pyramid.ome.tiff"
    y, x = np.indices((32, 40))
    rgb = np.zeros((32, 40, 3), dtype=np.uint8)
    rgb[..., 0] = 51
    rgb[..., 1] = (y * 40 + x) % 256
    rgb[..., 2] = 203
    rgb[3, 4:7, 1] = [0, 1, 255]
    with tifffile.TiffWriter(source, ome=True) as writer:
        writer.write(rgb, photometric="rgb", subifds=2, metadata={"axes": "YXS"})
        # Coarse planes deliberately differ from the original pixels.
        writer.write(np.full((16, 20, 3), 17, dtype=np.uint8), photometric="rgb", subfiletype=1)
        writer.write(np.full((8, 10, 3), 29, dtype=np.uint8), photometric="rgb", subfiletype=1)
    with tifffile.TiffFile(source) as pyramid:
        assert [level.shape for level in pyramid.series[0].levels] == [
            (32, 40, 3),
            (16, 20, 3),
            (8, 10, 3),
        ]
    before = reader.sha256(source)
    out = evidence.with_name(evidence.name + "-pyramid-crop")
    assert (
        reader.main(
            [
                "crop",
                "--source",
                str(source),
                "--out",
                str(out),
                "--x",
                "4",
                "--y",
                "3",
                "--width",
                "8",
                "--height",
                "6",
            ]
        )
        == 0
    )
    np.testing.assert_array_equal(reader.validate_crop(out / "crop.tiff", 8, 6), rgb[3:9, 4:12, 1])
    with tifffile.TiffFile(out / "crop.tiff") as crop:
        assert len(crop.pages) == len(crop.series) == len(crop.series[0].levels) == 1
        assert crop.pages[0].subifds is None
    assert reader.sha256(source) == before
    assert sorted(path.name for path in evidence.iterdir()) == ["pyramid.ome.tiff"]
    assert {path.name for path in out.iterdir()} == {
        "crop.tiff",
        "command.stdout.log",
        "command.stderr.log",
        "command.json",
        "receipt.json",
    }


@pytest.fixture
def tiff_page_fixture():
    if os.environ.get("PIXI_ENVIRONMENT_NAME") != "recovery":
        pytest.skip("TIFF-page command controls require locked recovery")
    tifffile = pytest.importorskip("tifffile")
    evidence = ROOT / ".pixi/runs/source-reader-controls" / uuid.uuid4().hex
    evidence.mkdir(parents=True)
    source = evidence / "copied.ome.btf"
    yy, xx = np.indices((17, 19))
    rgb = np.zeros((17, 19, 3), dtype=np.uint8)
    rgb[..., 0] = 51
    rgb[..., 1] = (yy * 19 + xx) % 256
    rgb[..., 2] = 203
    rgb[3, 4:7, 1] = [0, 1, 255]
    # Copied OME fields intentionally disagree with the real 19x17 single RGB page.
    description = (
        '<?xml version="1.0" encoding="UTF-8"?>'
        '<OME xmlns="http://www.openmicroscopy.org/Schemas/OME/2016-06">'
        '<Image ID="Image:0"><Pixels ID="Pixels:0" DimensionOrder="XYZCT" '
        'Type="uint16" SizeX="999" SizeY="888" SizeZ="4" SizeC="3" SizeT="2">'
        '<Channel ID="Channel:0" SamplesPerPixel="1"/><TiffData IFD="4"/>'
        "</Pixels></Image></OME>"
    )
    tifffile.imwrite(
        source, rgb, bigtiff=True, photometric="rgb", description=description, metadata=None
    )
    return tifffile, evidence, source, rgb, description


def test_tiff_page_bigtiff_ignores_copied_ome_and_never_full_hashes_or_decodes(
    tiff_page_fixture, monkeypatch
):
    tifffile, evidence, source, rgb, description = tiff_page_fixture
    original_hash = reader.sha256
    before = original_hash(source)

    def guarded_hash(path):
        assert Path(path).resolve() != source.resolve(), "Stat-only attempted a payload hash"
        return original_hash(path)

    monkeypatch.setattr(reader, "sha256", guarded_hash)
    original_asarray = tifffile.TiffPage.asarray

    def guarded_asarray(page, *args, **kwargs):
        assert Path(page.parent.filehandle.path).resolve() != source.resolve(), "Full source decode"
        return original_asarray(page, *args, **kwargs)

    monkeypatch.setattr(tifffile.TiffPage, "asarray", guarded_asarray)
    original_memmap = tifffile.memmap
    mappings = []

    def guarded_memmap(path, **kwargs):
        assert kwargs == {"page": 0, "mode": "r"}
        mappings.append(path)
        return original_memmap(path, **kwargs)

    monkeypatch.setattr(tifffile, "memmap", guarded_memmap)
    out = evidence.with_name(evidence.name + "-tiff-page-crop")
    assert (
        reader.main(
            [
                "crop",
                "--source",
                str(source),
                "--out",
                str(out),
                "--reader",
                "tiff-page",
                "--source-hash",
                "stat-only",
                "--x",
                "4",
                "--y",
                "3",
                "--width",
                "8",
                "--height",
                "6",
            ]
        )
        == 0
    )
    np.testing.assert_array_equal(reader.validate_crop(out / "crop.tiff", 8, 6), rgb[3:9, 4:12, 1])
    assert len(mappings) == 1
    receipt = json.loads((out / "receipt.json").read_text())
    assert receipt["reader"] == "tiff-page" and receipt["status"] == "complete"
    assert receipt["source"]["sha256"] is None
    assert receipt["source"]["full_payload_hash_computed"] is False
    assert receipt["source"]["hash_scope"] == "filesystem stat only; source payload not hashed"
    assert receipt["selection"] == {"page": 0, "sample": 1}
    assert receipt["command"]["header"]["shape"] == [17, 19, 3]
    assert "ignored" in receipt["coordinate_authority"]
    assert (out / "embedded-description.txt").read_text() == description
    metadata = evidence.with_name(evidence.name + "-tiff-page-metadata")
    assert (
        reader.main(
            [
                "metadata",
                "--source",
                str(source),
                "--out",
                str(metadata),
                "--reader",
                "tiff-page",
                "--source-hash",
                "stat-only",
            ]
        )
        == 0
    )
    header = json.loads((metadata / "metadata.json").read_text())
    assert header["width"] == 19 and header["height"] == 17
    assert header["dtype"] == "uint8" and header["samples_per_pixel"] == 3
    assert header["embedded_description"]["contains_ome_text"] is True
    assert header["embedded_ome_series_validity"].startswith("not evaluated")
    assert (metadata / "embedded-description.txt").read_text() == description
    assert len(mappings) == 1  # metadata never maps or reads pixels
    assert original_hash(source) == before


def test_tiff_page_bounds_rejected_before_mapping_or_decoding(tiff_page_fixture, monkeypatch):
    tifffile, evidence, source, _, _ = tiff_page_fixture

    def forbidden(*args, **kwargs):
        pytest.fail("Out-of-header-bounds crop attempted a source pixel read")

    monkeypatch.setattr(tifffile, "memmap", forbidden)
    monkeypatch.setattr(tifffile.TiffPage, "asarray", forbidden)
    out = evidence.with_name(evidence.name + "-tiff-page-bounds")
    with pytest.raises(ValueError, match="actual TIFF page header bounds"):
        reader.read_source(
            ROOT,
            "crop",
            source,
            out,
            reader="tiff-page",
            source_hash="stat-only",
            x=18,
            y=0,
            width=2,
            height=1,
        )
    receipt = json.loads((out / "receipt.json").read_text())
    assert receipt["status"] == receipt["command"]["status"] == "failed"
    assert not (out / "crop.tiff").exists()


@pytest.mark.parametrize("kind", ["uint16", "compressed", "extra_samples", "planar", "pages"])
def test_tiff_page_rejects_unsupported_encoding_before_mapping(
    tiff_page_fixture, monkeypatch, kind
):
    tifffile, evidence, _, rgb, _ = tiff_page_fixture
    source = evidence / "unsupported.btf"
    options = {"bigtiff": True, "photometric": "rgb", "metadata": None}
    if kind == "uint16":
        rgb = rgb.astype(np.uint16) + 256
    elif kind == "compressed":
        options["compression"] = "deflate"
    elif kind == "extra_samples":
        rgb = np.zeros((17, 19, 4), dtype=np.uint8)
    elif kind == "planar":
        rgb = np.moveaxis(rgb, 2, 0)
        options["planarconfig"] = "separate"
    elif kind == "pages":
        rgb = np.stack([rgb, rgb])
    tifffile.imwrite(source, rgb, **options)

    def forbidden(*args, **kwargs):
        pytest.fail("Rejected encoding attempted a source pixel read")

    monkeypatch.setattr(tifffile, "memmap", forbidden)
    monkeypatch.setattr(tifffile.TiffPage, "asarray", forbidden)
    out = evidence.with_name(evidence.name + "-tiff-page-rejection")
    with pytest.raises(ValueError, match="requires"):
        reader.read_source(
            ROOT,
            "crop",
            source,
            out,
            reader="tiff-page",
            source_hash="stat-only",
            x=0,
            y=0,
            width=2,
            height=2,
        )
    assert not (out / "crop.tiff").exists()


def test_tiff_page_stat_change_refuses_completion(tiff_page_fixture, monkeypatch):
    _, evidence, source, _, _ = tiff_page_fixture
    original_pin = reader.source_pin
    calls = []

    def changed_pin(path, hash_mode="full"):
        pin = original_pin(path, hash_mode)
        calls.append(pin)
        if len(calls) > 1:
            pin["mtime_ns"] += 1
        return pin

    monkeypatch.setattr(reader, "source_pin", changed_pin)
    out = evidence.with_name(evidence.name + "-tiff-page-stat-change")
    with pytest.raises(ValueError, match="changed during read"):
        reader.read_source(
            ROOT, "metadata", source, out, reader="tiff-page", source_hash="stat-only"
        )
    assert json.loads((out / "receipt.json").read_text())["status"] == "failed"


@pytest.fixture
def tiled_page_fixture(tiff_page_fixture):
    tifffile, evidence, _, _, description = tiff_page_fixture
    source = evidence / "tiled-copied.ome.btf"
    yy, xx = np.indices((777, 913))
    rgb = np.zeros((777, 913, 3), dtype=np.uint8)
    rgb[..., 0] = 41
    rgb[..., 1] = (yy * 13 + xx * 7) % 256
    rgb[..., 2] = 197
    rgb[511, 511:514, 1] = [0, 1, 255]
    tifffile.imwrite(
        source,
        rgb,
        bigtiff=True,
        tile=(512, 512),
        photometric="rgb",
        description=description,
        metadata=None,
    )
    return tifffile, evidence, source, rgb, description


@pytest.mark.parametrize(
    "bounds,indices", [((507, 509, 13, 11), [0, 1, 2, 3]), ((905, 770, 8, 7), [3])]
)
def test_tiled_page_exact_green_seam_and_partial_image_edge(
    tiled_page_fixture, monkeypatch, bounds, indices
):
    tifffile, evidence, source, rgb, _ = tiled_page_fixture
    original_hash = reader.sha256

    def guarded_hash(path):
        assert Path(path).resolve() != source.resolve(), "Stat-only source payload hash"
        return original_hash(path)

    monkeypatch.setattr(reader, "sha256", guarded_hash)
    original_asarray = tifffile.TiffPage.asarray

    def guarded_asarray(page, *args, **kwargs):
        assert Path(page.parent.filehandle.path).resolve() != source.resolve(), "Full source decode"
        return original_asarray(page, *args, **kwargs)

    monkeypatch.setattr(tifffile.TiffPage, "asarray", guarded_asarray)

    def forbidden(*args, **kwargs):
        pytest.fail("Tiled source attempted memory-map or all-segments decoding")

    monkeypatch.setattr(tifffile, "memmap", forbidden)
    monkeypatch.setattr(tifffile.TiffPage, "segments", forbidden)
    x, y, width, height = bounds
    out = evidence.with_name(evidence.name + "-tiled-crop")
    reader.read_source(
        ROOT,
        "crop",
        source,
        out,
        reader="tiff-page",
        source_hash="stat-only",
        x=x,
        y=y,
        width=width,
        height=height,
    )
    np.testing.assert_array_equal(
        reader.validate_crop(out / "crop.tiff", width, height),
        rgb[y : y + height, x : x + width, 1],
    )
    receipt = json.loads((out / "receipt.json").read_text())
    assert receipt["source"]["sha256"] is None
    assert receipt["status"] == "complete"
    operation = receipt["command"]
    assert operation["header"]["memmappable"] is False
    assert [item["tile_index"] for item in operation["tile_reads"]] == indices
    assert operation["decoded_tile_count"] == len(indices)
    assert operation["encoded_bytes_read"] == len(indices) * 512 * 512 * 3
    assert "canonical decode" in operation["method"]


@pytest.mark.parametrize(
    "kind",
    [
        "compressed",
        "uint16",
        "orientation",
        "missing_offset",
        "missing_count",
        "truncated",
        "decoder_none",
        "outbounds",
    ],
)
def test_tiled_page_rejects_bad_encoding_or_missing_pixels(tiled_page_fixture, monkeypatch, kind):
    tifffile, evidence, _, rgb, description = tiled_page_fixture
    source = evidence / "bad-tiled.btf"
    options = {
        "bigtiff": True,
        "tile": (512, 512),
        "photometric": "rgb",
        "description": description,
        "metadata": None,
    }
    if kind == "compressed":
        options["compression"] = "deflate"
    elif kind == "uint16":
        rgb = rgb.astype(np.uint16) + 256
    elif kind == "orientation":
        options["extratags"] = [(274, "H", 1, 6, False)]
    tifffile.imwrite(source, rgb, **options)
    if kind in {"missing_offset", "missing_count"}:
        with tifffile.TiffFile(source, mode="r+b") as image:
            tag = image.pages[0].tags[
                "TileOffsets" if kind == "missing_offset" else "TileByteCounts"
            ]
            values = list(tag.value)
            values[0] = 0
            tag.overwrite(tuple(values))
    elif kind == "truncated":
        with source.open("r+b") as handle:
            handle.truncate(source.stat().st_size - 1)
    elif kind == "decoder_none":
        monkeypatch.setattr(
            tifffile.TiffPage,
            "decode",
            property(lambda page: lambda data, index: (None, (0, 0, 0, 0, 0), (1, 512, 512, 3))),
        )

    def forbidden(*args, **kwargs):
        pytest.fail("Rejected tiled source attempted full-page decoding")

    monkeypatch.setattr(tifffile.TiffPage, "asarray", forbidden)
    monkeypatch.setattr(tifffile.TiffPage, "segments", forbidden)
    monkeypatch.setattr(tifffile, "memmap", forbidden)
    out = evidence.with_name(evidence.name + "-tiled-rejected")
    bounds = (
        {"x": 907, "y": 771, "width": 6, "height": 6}
        if kind == "truncated"
        else {"x": 0, "y": 0, "width": 2, "height": 2}
    )
    if kind == "outbounds":
        bounds.update(x=912, width=2)
    with pytest.raises(ValueError):
        reader.read_source(
            ROOT, "crop", source, out, reader="tiff-page", source_hash="stat-only", **bounds
        )
    receipt = json.loads((out / "receipt.json").read_text())
    assert receipt["status"] == "failed"
    assert not (out / "crop.tiff").exists()


@pytest.mark.parametrize("kind", ["tile_size", "tile_region"])
def test_tiled_page_read_bounds_are_enforced_before_decode(tiff_page_fixture, monkeypatch, kind):
    tifffile, evidence, _, _, description = tiff_page_fixture
    source = evidence / "large-tiles.btf"
    tile = (1024, 1024) if kind == "tile_size" else (512, 512)
    shape = (1024, 1024, 3) if kind == "tile_size" else (512, 8192, 3)
    tifffile.imwrite(
        source,
        np.zeros(shape, dtype=np.uint8),
        bigtiff=True,
        tile=tile,
        photometric="rgb",
        description=description,
        metadata=None,
    )

    def forbidden(*args, **kwargs):
        pytest.fail("Oversized tiles attempted decode")

    monkeypatch.setattr(tifffile.TiffPage, "decode", property(lambda page: forbidden))
    out = evidence.with_name(evidence.name + "-tile-bound")
    width = 1 if kind == "tile_size" else 8192
    with pytest.raises(ValueError):
        reader.read_source(
            ROOT,
            "crop",
            source,
            out,
            reader="tiff-page",
            source_hash="stat-only",
            x=0,
            y=0,
            width=width,
            height=1,
        )
    assert not (out / "crop.tiff").exists()
