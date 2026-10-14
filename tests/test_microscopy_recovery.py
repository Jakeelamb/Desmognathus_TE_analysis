"""Public-command controls for immutable nuclear recovery infrastructure."""

from __future__ import annotations

import csv
import importlib.util
import math
import sys
import zipfile
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

Image = pytest.importorskip(
    "PIL.Image", reason="Recovery image controls require the image-qc or recovery Pixi environment"
)

tifffile = pytest.importorskip(
    "tifffile", reason="Recovery image controls require the image-qc or recovery Pixi environment"
)

PATH = Path(__file__).resolve().parents[1] / "scripts/microscopy/recovery.py"
SPEC = importlib.util.spec_from_file_location("microscopy_recovery", PATH)
recovery = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = recovery
SPEC.loader.exec_module(recovery)


def dump(path, value):
    path.write_bytes(recovery.encoded(value))


def rebind_crosswalk(root, value):
    path = root / "provenance.private.json"
    dump(path, {"cases": value["cases"], "legacy_fields_are_evidence_only": True})
    value["assets"][str(path.relative_to(root))] = recovery.digest(path)


def fake_measure(pixels, background, tools):
    assert pixels.dtype == np.uint8 and pixels.ndim == 1 and pixels.size
    iod = math.fsum(math.log10(background / max(int(pixel), 1)) for pixel in pixels)
    return {
        "count": pixels.size,
        "raw_sum": int(pixels.astype(np.uint64).sum()),
        "min_pixel": int(pixels.min()),
        "max_pixel": int(pixels.max()),
        "zero_pixels": int((pixels == 0).sum()),
        "one_pixels": int((pixels == 1).sum()),
        "saturated_pixels": int((pixels == 255).sum()),
        "background": background,
        "iod": iod,
        "mean_od": iod / pixels.size,
        "median_pixel": float(np.median(pixels)),
    }


@pytest.fixture
def pilot(tmp_path, monkeypatch):
    paper = tmp_path / "paper"
    cell = tmp_path / "cell"
    paper.mkdir()
    cell.mkdir()
    standards, targets, frozen = [], [], []
    arrays = {}
    for process, count in recovery.COUNTS.items():
        standard = process in {"413", "414"}
        name = f"Process_{process}_{'green' if standard else 'raw_green'}.ome.tiff"
        directory = cell / "pixels" / process
        directory.mkdir(parents=True)
        raw = np.full((80, 80), 240, dtype=np.uint8)
        labels = np.zeros(raw.shape, dtype=np.uint32)
        for label in range(1, count + 1):
            x, y = 8 + ((label - 1) % 8) * 8, 8 + ((label - 1) // 8) * 8
            labels[y : y + 2, x : x + 3] = label
            raw[y : y + 2, x : x + 3] = (100 + label) if process != "320" else 230
        raw_path, mask_path, old_path = (
            directory / "raw.tiff",
            directory / "labels.tiff",
            directory / "old.tiff",
        )
        tifffile.imwrite(raw_path, raw, metadata={"axes": "YX"})
        tifffile.imwrite(mask_path, labels, metadata={"axes": "YX"})
        tifffile.imwrite(old_path, labels, metadata={"axes": "YX"})
        manifest = directory / "manifest.json"
        dump(manifest, {"fixture": True})
        arrays[process] = (raw_path, mask_path)
        if standard:
            frozen.extend(
                {
                    "root": "cell",
                    "path": str(path.relative_to(cell)),
                    "sha256": recovery.digest(path),
                }
                for path in (raw_path, mask_path)
            )
        for label in range(1, count + 1):
            yy, xx = np.nonzero(labels == label)
            background = recovery.BACKGROUNDS[process]
            iod = fake_measure(raw[labels == label], background, None)["iod"]
            row = {
                "filename": name,
                "specimen_id": recovery.SPECIMENS[process],
                "species": "D. fuscus",
                "image_type": "pmount" if standard else "brightfield",
                "tile_name": "tile_y002000_x001000.tiff",
                "tile_x0": 1000,
                "tile_y0": 2000,
                "tile_width_px": 80,
                "tile_height_px": 80,
            }
            prefix = "" if standard else "nuc_"
            row["label" if standard else "nucleus_label"] = label
            row[prefix + "area_px"] = len(xx)
            row[prefix + "iod"] = f"{iod:.6f}"
            row[prefix + "centroid_x"] = float(xx.mean()) + (1000 if standard else 0)
            row[prefix + "centroid_y"] = float(yy.mean()) + (2000 if standard else 0)
            row[prefix + "i_bg"] = background
            row["source_image_path" if standard else "nucleus_source_image_path"] = str(raw_path)
            row["mask_path" if standard else "nucleus_mask_path"] = str(mask_path)
            row["tile_manifest_path" if standard else "nucleus_tile_manifest_path"] = str(manifest)
            row["run_manifest_path" if standard else "nucleus_run_manifest_path"] = str(manifest)
            if not standard:
                row["nucleus_mask_path_original"] = str(old_path)
                row["nucleus_object_id"] = f"{name}::{row['tile_name']}::{label}"
                row["review_decision"] = "keep"
                row["nucleus_repair_applied"] = "True"
            (standards if standard else targets).append(row)
    standard_path = cell / recovery.STANDARD_REL
    target_path = paper / recovery.TARGET_REL
    for path, rows in ((standard_path, standards), (target_path, targets)):
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
    pins = tmp_path / "frozen.json"
    dump(pins, {"sources": frozen})
    ui = tmp_path / "ui.html"
    ui.write_text("<!doctype html><title>Fixture review</title>")
    monkeypatch.setattr(recovery, "STANDARD_SHA", recovery.digest(standard_path))
    monkeypatch.setattr(recovery, "TARGET_SHA", recovery.digest(target_path))
    monkeypatch.setattr(recovery, "STANDARD_PINS_PATH", pins)
    monkeypatch.setattr(recovery, "STANDARD_PINS_SHA", recovery.digest(pins))
    monkeypatch.setattr(recovery, "UI_PATH", ui)
    native = SimpleNamespace(
        load_toolchain=lambda root: SimpleNamespace(
            identity="fixture-native", java=Path("/fixture/java"), ij_jar=Path("/fixture/ij.jar")
        ),
        measure_pixels=fake_measure,
    )
    monkeypatch.setattr(recovery, "backend", lambda: native)
    return SimpleNamespace(
        root=tmp_path,
        paper=paper,
        cell=cell,
        standard=standard_path,
        target=target_path,
        arrays=arrays,
        native=native,
    )


def run(pilot, *arguments):
    return recovery.main([str(argument) for argument in arguments])


def packet(pilot):
    inv = pilot.root / "inventory"
    out = pilot.root / "packet"
    assert (
        run(
            pilot, "inventory", "--paper-root", pilot.paper, "--cell-root", pilot.cell, "--out", inv
        )
        == 0
    )
    assert (
        run(pilot, "pack", "--inventory", inv / "inventory.json", "--margin", 4, "--out", out) == 0
    )
    return out, recovery.load_json(out / "packet.json")


def authored(case, reviewer="alice", decision="keep"):
    return {
        "case_id": case["case_id"],
        "state_sha256": case["state_sha256"],
        "reviewer": reviewer,
        "time": "2026-10-13T12:00:00Z",
        "decision": decision,
        "reasons": ["complete_isolated_nucleus" if decision == "keep" else "uncertain_boundary"],
        "note": "fixture",
    }


def test_all_82_originals_blind_lossless_and_signed_replay(pilot):
    root, value = packet(pilot)
    assert run(pilot, "validate", "--packet", root) == 0
    assert len(value["cases"]) == 82 and value["proposals"] == value["decisions"] == []
    public = recovery.load_json(root / "reviewer" / "cases.json")
    assert all(set(case) == recovery.CASE_KEYS for case in public)
    text = (root / "reviewer" / "cases.json").read_text()
    assert not any(
        word in text for word in ("Process_", "fuscus", "iod", "specimen", str(pilot.cell))
    )
    for asset in (root / "reviewer" / "assets").rglob("*.png"):
        with Image.open(asset) as image:
            assert not image.info and image.mode == "L"
    negative = [case for case in value["cases"] if case["process"] == "320"]
    assert all(case["recorded_background_replay"]["iod"] < 0 for case in negative)
    assert all(case["shared_255_sensitivity"]["iod"] > 0 for case in negative)
    report = pilot.root / "report"
    assert run(pilot, "report", "--packet", root, "--out", report) == 0
    derived = recovery.load_json(report / "report.json")
    assert (
        len(derived["rows"]) == 82
        and derived["authored_decisions"] == derived["authored_proposals"] == 0
    )
    assert derived["review_counts"] == {
        "missing": 82,
        "keep": 0,
        "problem": 0,
        "unsure": 0,
        "conflict": 0,
    }
    assert derived["review_approved_subset"]["case_ids"] == []


@pytest.mark.parametrize(
    "change",
    [
        "missing",
        "extra",
        "wrong_species",
        "wrong_label",
        "wrong_background",
        "wrong_origin",
        "wrong_cohort",
    ],
)
def test_authority_and_coordinate_errors_refuse_inventory(pilot, monkeypatch, change):
    rows = recovery.read_csv(pilot.target)
    if change == "missing":
        rows.pop()
    elif change == "extra":
        rows.append(rows[0].copy())
    elif change == "wrong_species":
        rows[0]["species"] = "D. apalachicolae"
    elif change == "wrong_label":
        rows[0]["nucleus_label"] = "1.5"
    elif change == "wrong_background":
        rows[0]["nuc_i_bg"] = "255" if rows[0]["nuc_i_bg"] != "255" else "254"
    elif change == "wrong_origin":
        rows[0]["nuc_centroid_x"] = str(float(rows[0]["nuc_centroid_x"]) + 1000)
    else:
        rows[0]["specimen_id"] = "wrong"
    with pilot.target.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    monkeypatch.setattr(recovery, "TARGET_SHA", recovery.digest(pilot.target))
    assert (
        run(
            pilot,
            "inventory",
            "--paper-root",
            pilot.paper,
            "--cell-root",
            pilot.cell,
            "--out",
            pilot.root / "bad",
        )
        == 1
    )
    assert not (pilot.root / "bad").exists()


@pytest.mark.parametrize("change", ["channel", "dtype", "axes", "mask_area"])
def test_pixel_format_and_area_errors_refuse_inventory(pilot, change):
    raw_path, mask_path = pilot.arrays["314"]
    raw = tifffile.imread(raw_path)
    if change == "channel":
        tifffile.imwrite(raw_path, np.stack([raw] * 3, axis=-1), metadata={"axes": "YXS"})
    elif change == "dtype":
        tifffile.imwrite(raw_path, raw.astype(np.uint16), metadata={"axes": "YX"})
    elif change == "axes":
        tifffile.imwrite(raw_path, raw[None], metadata={"axes": "ZYX"})
    else:
        labels = tifffile.imread(mask_path)
        labels[8, 8] = 0
        tifffile.imwrite(mask_path, labels, metadata={"axes": "YX"})
    assert (
        run(
            pilot,
            "inventory",
            "--paper-root",
            pilot.paper,
            "--cell-root",
            pilot.cell,
            "--out",
            pilot.root / "bad",
        )
        == 1
    )


def test_writers_refuse_existing_symlink_and_input_ancestry_without_writes(pilot):
    existing = pilot.root / "existing"
    existing.mkdir()
    args = ("inventory", "--paper-root", pilot.paper, "--cell-root", pilot.cell)
    assert run(pilot, *args, "--out", existing) == 1
    assert list(existing.iterdir()) == []
    assert run(pilot, *args, "--out", pilot.cell / "new") == 1
    assert not (pilot.cell / "new").exists()
    link = pilot.root / "link"
    link.symlink_to(pilot.root, target_is_directory=True)
    assert run(pilot, *args, "--out", link / "new") == 1
    assert not (pilot.root / "new").exists()
    assert run(pilot, *args, "--out", link / ".." / "escaped") == 1
    assert not (pilot.root.parent / "escaped").exists()


@pytest.mark.parametrize("tree", ["data", ".venv"])
def test_inventory_refuses_new_artifacts_inside_preserved_trees(pilot, tree):
    protected = pilot.paper / tree
    protected.mkdir()
    assert (
        run(
            pilot,
            "inventory",
            "--paper-root",
            pilot.paper,
            "--cell-root",
            pilot.cell,
            "--out",
            protected / "new",
        )
        == 1
    )
    assert list(protected.iterdir()) == []


def test_changed_source_and_rebound_inventory_scalars_rejected(pilot):
    root, _ = packet(pilot)
    inventory_path = pilot.root / "inventory" / "inventory.json"
    inv = recovery.load_json(inventory_path)
    inv["cases"][0]["stored_iod"] += 10
    dump(inventory_path, inv)
    assert run(pilot, "pack", "--inventory", inventory_path, "--out", pilot.root / "bad") == 1
    assert run(pilot, "validate", "--packet", root) == 1


def test_source_change_outside_selected_crop_is_refused(pilot):
    root, value = packet(pilot)
    raw_path = Path(value["cases"][0]["source"]["path"])
    raw = tifffile.imread(raw_path)
    raw[-1, -1] -= 1
    tifffile.imwrite(raw_path, raw, metadata={"axes": "YX"})
    assert run(pilot, "validate", "--packet", root) == 1


@pytest.mark.parametrize(
    "change",
    ["missing_case", "extra_case", "cohort", "origin", "measurement", "extra_reviewer_asset"],
)
def test_packet_identity_contract_and_measurements_cannot_be_rebound(pilot, change, capsys):
    root, value = packet(pilot)
    if change == "missing_case":
        value["cases"].pop()
    elif change == "extra_case":
        value["cases"].append(value["cases"][0].copy())
    elif change == "extra_reviewer_asset":
        (root / "reviewer" / "source-map.json").write_text("{}")
    else:
        case = value["cases"][0]
        if change == "cohort":
            case["cohort"] = "wrong"
        elif change == "origin":
            case["source"]["tile_global_x0"] += 1
        else:
            case["recorded_background_replay"]["iod"] += 1
        private = root / "cases" / case["case_id"] / "case.json"
        dump(private, case)
        value["assets"][str(private.relative_to(root))] = recovery.digest(private)
        rebind_crosswalk(root, value)
    dump(root / "packet.json", value)
    assert run(pilot, "validate", "--packet", root) == 1
    if change == "measurement":
        assert "Stored native measurement differs from fresh reduction" in capsys.readouterr().err


@pytest.mark.parametrize(
    "change", ["extra_file", "html", "protocol", "provenance", "missing_ui_pin", "extra_directory"]
)
def test_packet_blinding_refuses_rehashed_metadata_changes(pilot, change):
    root, value = packet(pilot)
    if change == "missing_ui_pin":
        del value["source_pins"][str(recovery.UI_PATH)]
    elif change == "extra_directory":
        (root / "reviewer" / "Process_320").mkdir()
    else:
        relative = {
            "extra_file": "reviewer/source-map.json",
            "html": "reviewer/index.html",
            "protocol": "reviewer/protocol.txt",
            "provenance": "provenance.private.json",
        }[change]
        artifact = root / relative
        if change == "provenance":
            crosswalk = recovery.load_json(artifact)
            crosswalk["unvalidated_crosswalk"] = "Process_320"
            dump(artifact, crosswalk)
        else:
            artifact.write_text("Process_320 specimen 3401\n")
        value["assets"][relative] = recovery.digest(artifact)
    dump(root / "packet.json", value)
    assert run(pilot, "validate", "--packet", root) == 1


@pytest.mark.parametrize("asset", ["raw.png", "original_mask.png", "known_nuclei.png"])
def test_changed_packet_pixels_rejected_even_when_file_hashes_rebound(pilot, asset, capsys):
    root, value = packet(pilot)
    case = value["cases"][0]
    path = root / "cases" / case["case_id"] / asset
    image = np.array(Image.open(path))
    image[0, 0] = (int(image[0, 0]) + 1) % 256 if asset == "raw.png" else 255 - image[0, 0]
    recovery.save_png(path, image)
    value["assets"][str(path.relative_to(root))] = recovery.digest(path)
    if asset in {"raw.png", "original_mask.png"}:
        public = (
            root
            / "reviewer"
            / "assets"
            / case["case_id"]
            / ("raw.png" if asset == "raw.png" else "mask.png")
        )
        recovery.save_png(public, image)
        value["assets"][str(public.relative_to(root))] = recovery.digest(public)
        case["crop"]["raw_sha256" if asset == "raw.png" else "mask_sha256"] = recovery.digest(path)
        case["state_sha256"] = recovery.value_hash(case["crop"])
        private_record = root / "cases" / case["case_id"] / "case.json"
        dump(private_record, case)
        value["assets"][str(private_record.relative_to(root))] = recovery.digest(private_record)
        queue_path = root / "reviewer" / "cases.json"
        queue = recovery.load_json(queue_path)
        next(row for row in queue if row["case_id"] == case["case_id"])["state_sha256"] = case[
            "state_sha256"
        ]
        dump(queue_path, queue)
        value["assets"][str(queue_path.relative_to(root))] = recovery.digest(queue_path)
        rebind_crosswalk(root, value)
    dump(root / "packet.json", value)
    assert run(pilot, "validate", "--packet", root) == 1
    expected = {
        "raw.png": "Crop pixels changed from original source",
        "original_mask.png": "Original mask membership changed",
        "known_nuclei.png": "Known-nucleus support changed",
    }
    assert expected[asset] in capsys.readouterr().err


@pytest.mark.parametrize(
    "change", ["stale", "duplicate", "blank", "unnamed", "untimed", "bad_time", "reasonless"]
)
def test_invalid_decision_exports_refused(pilot, change):
    root, value = packet(pilot)
    event = authored(value["cases"][0])
    events = [event]
    if change == "stale":
        event["state_sha256"] = "0" * 64
    elif change == "duplicate":
        events.append(event.copy())
    elif change == "blank":
        event["decision"] = ""
    elif change == "unnamed":
        event["reviewer"] = " "
    elif change == "untimed":
        event["time"] = "2026-10-13T12:00:00"
    elif change == "bad_time":
        event["time"] = "2026-99-13T12:00:00Z"
    else:
        event["reasons"] = []
    source = pilot.root / "draft.json"
    dump(source, events)
    assert (
        run(pilot, "decision", "--packet", root, "--event", source, "--out", pilot.root / "events")
        == 1
    )
    assert not (pilot.root / "events").exists()


def test_named_reviewers_keep_conflicts_and_duplicate_collisions_visible(pilot):
    root, value = packet(pilot)
    source = pilot.root / "draft.json"
    dump(source, [authored(value["cases"][0]), authored(value["cases"][0], "bob", "problem")])
    events = pilot.root / "events"
    assert run(pilot, "decision", "--packet", root, "--event", source, "--out", events) == 0
    assert run(pilot, "decision", "--packet", root, "--event", source, "--out", events) == 1
    assert (
        run(pilot, "report", "--packet", root, "--events", events, "--out", pilot.root / "report")
        == 0
    )
    report = recovery.load_json(pilot.root / "report" / "report.json")
    assert report["review_counts"]["conflict"] == 1 and report["review_counts"]["missing"] == 81
    assert report["review_approved_subset"]["case_ids"] == []
    assert (
        run(
            pilot,
            "report",
            "--packet",
            root,
            "--events",
            events,
            events,
            "--out",
            pilot.root / "duplicate_report",
        )
        == 1
    )


def proposal_inputs(pilot, root, case):
    roi = pilot.root / "authored.zip"
    roi.write_bytes(b"fixture-not-real-roi")
    sidecar = pilot.root / "authored.json"
    dump(
        sidecar,
        {
            "packet_sha256": recovery.digest(root / "packet.json"),
            "case_id": case["case_id"],
            "crop_sha256": case["crop"]["raw_sha256"],
            "parent_state_sha256": case["state_sha256"],
            "author": "alice",
            "time": "2026-10-13T12:00:00Z",
            "rationale": "Synthetic control only",
        },
    )
    return roi, sidecar


def test_prepare_editor_makes_disposable_copies_and_blank_authorship(pilot):
    root, value = packet(pilot)
    case = value["cases"][0]
    out = pilot.root / "workspace"
    assert (
        run(pilot, "prepare-editor", "--packet", root, "--case", case["case_id"], "--out", out) == 0
    )
    sidecar = recovery.load_json(out / "authorship.template.json")
    assert sidecar["author"] == sidecar["time"] == sidecar["rationale"] == ""
    assert sidecar["crop_sha256"] == recovery.digest(out / "raw.png")
    editor = recovery.load_json(out / "editor.json")
    assert editor["launch_argv"][1:3] == [
        f"-Duser.home={out / '.imagej-home'}",
        f"-Djava.util.prefs.userRoot={out / '.imagej-prefs'}",
    ]
    assert (out / ".imagej-home").is_dir() and (out / ".imagej-prefs").is_dir()
    assert "-Dsun.awt.disableGtkFileDialogs=true" in editor["launch_argv"]
    assert "subprocess.run" in (out / "README.txt").read_text()
    assert (out / "raw.png").read_bytes() == (
        root / "cases" / case["case_id"] / "raw.png"
    ).read_bytes()
    assert not list(out.glob("*.roi")) and not list(out.glob("*.zip"))
    assert (
        run(
            pilot,
            "prepare-editor",
            "--packet",
            root,
            "--case",
            "unknown",
            "--out",
            pilot.root / "unknown",
        )
        == 1
    )
    recovery.save_png(out / "raw.png", np.zeros((2, 2), dtype=np.uint8))
    assert run(pilot, "validate", "--packet", root) == 0


def fake_rasterizer(pilot, root, case, *, overlap=False, error=False):
    original = np.array(Image.open(root / "cases" / case["case_id"] / "original_mask.png"))

    def rasterize(roi, width, height, out, tools):
        if error:
            raise ValueError("Unsupported ROI geometry fixture")
        assert original.shape == (height, width)
        out.mkdir()
        nucleus = original.copy()
        nucleus[0, 0] = 255
        background = original.copy() if overlap else np.zeros_like(original)
        if not overlap:
            background[-1, -1] = 255
        for name, pixels in (("nucleus", nucleus), ("background_001", background)):
            recovery.save_png(out / f"{name}.png", pixels)
        return {name: out / f"{name}.png" for name in ("nucleus", "background_001")}

    pilot.native.rasterize_rois = rasterize


def test_import_is_proposal_only_and_candidate_reference_is_separate(pilot):
    root, value = packet(pilot)
    case = value["cases"][0]
    roi, sidecar = proposal_inputs(pilot, root, case)
    fake_rasterizer(pilot, root, case)
    out = pilot.root / "proposal"
    assert (
        run(pilot, "import-roi", "--packet", root, "--roi", roi, "--sidecar", sidecar, "--out", out)
        == 0
    )
    proposal = recovery.load_json(out / "proposal.json")
    assert proposal["decision"] is None and proposal["geometry"]["added_px"] == 1
    assert (
        proposal["candidate_backgrounds"]["background_001"]["status"]
        == "authored_candidate_reference_only"
    )
    assert (out / "roi.zip").read_bytes() == roi.read_bytes()
    assert (
        run(pilot, "report", "--packet", root, "--proposals", out, "--out", pilot.root / "report")
        == 0
    )
    report = recovery.load_json(pilot.root / "report" / "report.json")
    assert report["review_counts"]["missing"] == 82 and report["authored_proposals"] == 1
    assert len(report["proposal_sensitivities"]) == 1
    assert report["review_approved_subset"]["case_ids"] == []
    assert recovery.load_json(root / "packet.json")["proposals"] == []


def test_single_roi_wraps_exact_bytes_for_existing_backend_deterministically(pilot):
    root, value = packet(pilot)
    case = value["cases"][0]
    _, sidecar = proposal_inputs(pilot, root, case)
    source = pilot.root / "nucleus.roi"
    original_bytes = b"Synthetic ROI fixture\x00\xff\x01"
    source.write_bytes(original_bytes)
    source_hash = recovery.digest(source)
    observed = []

    def inspect_archive(roi, width, height, out, tools):
        assert roi.name == "roi.zip"
        with zipfile.ZipFile(roi) as archive:
            assert archive.namelist() == ["nucleus.roi"]
            entry = archive.getinfo("nucleus.roi")
            assert entry.date_time == (1980, 1, 1, 0, 0, 0)
            assert entry.compress_type == zipfile.ZIP_STORED
            assert archive.read(entry) == original_bytes
        observed.append(roi.read_bytes())
        out.mkdir()
        pixels = np.array(Image.open(root / "cases" / case["case_id"] / "original_mask.png"))
        assert pixels.shape == (height, width)
        pixels[0, 0] = 255
        recovery.save_png(out / "nucleus.png", pixels)
        return {"nucleus": out / "nucleus.png"}

    pilot.native.rasterize_rois = inspect_archive
    for name in ("single-proposal-one", "single-proposal-two"):
        out = pilot.root / name
        assert (
            run(
                pilot,
                "import-roi",
                "--packet",
                root,
                "--roi",
                source,
                "--sidecar",
                sidecar,
                "--out",
                out,
            )
            == 0
        )
        proposal = recovery.load_json(out / "proposal.json")
        assert proposal["input_pins"][str(source)] == source_hash
        assert proposal["roi_sha256"] == recovery.digest(out / "roi.zip")
        assert proposal["decision"] is None
        assert proposal["background_sha256"] == proposal["candidate_backgrounds"] == {}
    assert observed[0] == observed[1]
    assert source.read_bytes() == original_bytes and recovery.digest(source) == source_hash
    assert run(pilot, "validate", "--packet", root) == 0


@pytest.mark.parametrize(
    "change",
    [
        "basename",
        "case",
        "oversize",
        "empty",
        "symlink",
        "existing_destination",
        "protected_destination",
    ],
)
def test_single_roi_rejections_preserve_source_and_destination(pilot, change):
    root, value = packet(pilot)
    case = value["cases"][0]
    _, sidecar = proposal_inputs(pilot, root, case)
    source = pilot.root / (
        "other.roi"
        if change == "basename"
        else "Nucleus.roi"
        if change == "case"
        else "nucleus.roi"
    )
    source.write_bytes(
        b"x" * (1024 * 1024 + 1)
        if change == "oversize"
        else b""
        if change == "empty"
        else b"Synthetic ROI fixture"
    )
    original_hash = recovery.digest(source)
    if change == "symlink":
        saved = pilot.root / "saved.roi"
        source.rename(saved)
        source.symlink_to(saved)
    out = pilot.root / "refused-single"
    if change == "existing_destination":
        out.mkdir()
    elif change == "protected_destination":
        out = pilot.cell / "forbidden-single"
    assert (
        run(
            pilot,
            "import-roi",
            "--packet",
            root,
            "--roi",
            source,
            "--sidecar",
            sidecar,
            "--out",
            out,
        )
        == 1
    )
    assert recovery.digest(source) == original_hash
    assert not out.exists() or list(out.iterdir()) == []


def test_proposal_state_review_and_explicit_report_policy(pilot):
    root, value = packet(pilot)
    case = value["cases"][0]
    roi, sidecar = proposal_inputs(pilot, root, case)
    fake_rasterizer(pilot, root, case)
    proposal_root = pilot.root / "proposal"
    assert (
        run(
            pilot,
            "import-roi",
            "--packet",
            root,
            "--roi",
            roi,
            "--sidecar",
            sidecar,
            "--out",
            proposal_root,
        )
        == 0
    )
    proposal = recovery.load_json(proposal_root / "proposal.json")
    event = authored(case)
    event["state_sha256"] = proposal["state_sha256"]
    event_path = pilot.root / "proposal-review.json"
    dump(event_path, [event])
    assert (
        run(
            pilot,
            "decision",
            "--packet",
            root,
            "--event",
            event_path,
            "--out",
            pilot.root / "stale",
        )
        == 1
    )
    events = pilot.root / "proposal-events"
    assert (
        run(
            pilot,
            "decision",
            "--packet",
            root,
            "--event",
            event_path,
            "--proposals",
            proposal_root,
            "--out",
            events,
        )
        == 0
    )
    policy = pilot.root / "policy.json"
    dump(
        policy,
        {
            "name": "synthetic_explicit_state",
            "state_selection": {case["case_id"]: proposal["state_sha256"]},
        },
    )
    out = pilot.root / "selected-report"
    assert (
        run(
            pilot,
            "report",
            "--packet",
            root,
            "--events",
            events,
            "--proposals",
            proposal_root,
            "--policy",
            policy,
            "--out",
            out,
        )
        == 0
    )
    report = recovery.load_json(out / "report.json")
    assert report["review_counts"]["keep"] == 1 and report["review_counts"]["missing"] == 81
    assert report["review_approved_subset"]["case_ids"] == [case["case_id"]]
    assert (
        next(row for row in report["rows"] if row["case_id"] == case["case_id"])["proposal"][
            "geometry"
        ]["added_px"]
        == 1
    )
    proposal["geometry"]["added_px"] = 0
    dump(proposal_root / "proposal.json", proposal)
    assert (
        run(
            pilot,
            "report",
            "--packet",
            root,
            "--proposals",
            proposal_root,
            "--out",
            pilot.root / "tampered",
        )
        == 1
    )


def queue_proposal(pilot):
    root, value = packet(pilot)
    case = value["cases"][0]
    roi, sidecar = proposal_inputs(pilot, root, case)
    fake_rasterizer(pilot, root, case)
    proposal_root = pilot.root / "proposal"
    assert (
        run(
            pilot,
            "import-roi",
            "--packet",
            root,
            "--roi",
            roi,
            "--sidecar",
            sidecar,
            "--out",
            proposal_root,
        )
        == 0
    )
    proposal = recovery.load_json(proposal_root / "proposal.json")
    policy_path = pilot.root / "policy.json"
    dump(
        policy_path,
        {
            "name": "synthetic_selected_proposal",
            "state_selection": {case["case_id"]: proposal["state_sha256"]},
        },
    )
    return root, value, case, proposal_root, proposal, policy_path


def test_review_queue_blinds_selected_proposal_and_exports_bound_state(pilot):
    root, value, case, proposal_root, proposal, policy = queue_proposal(pilot)
    out = pilot.root / "queue"
    assert (
        run(
            pilot,
            "review-queue",
            "--packet",
            root,
            "--proposals",
            proposal_root,
            "--policy",
            policy,
            "--out",
            out,
        )
        == 0
    )
    queue = recovery.load_json(out / "reviewer" / "cases.json")
    assert len(queue) == 82 and all(set(row) == recovery.CASE_KEYS for row in queue)
    assert [row["case_id"] for row in queue] == [row["case_id"] for row in value["cases"]]
    selected = next(row for row in queue if row["case_id"] == case["case_id"])
    assert selected["state_sha256"] == proposal["state_sha256"]
    assert recovery.digest(out / "reviewer" / selected["mask"]) == proposal["nucleus_sha256"]
    for original, public in zip(value["cases"], queue, strict=True):
        assert recovery.digest(out / "reviewer" / public["raw"]) == original["crop"]["raw_sha256"]
        if original["case_id"] != case["case_id"]:
            assert public["state_sha256"] == original["state_sha256"]
            assert (
                recovery.digest(out / "reviewer" / public["mask"])
                == original["crop"]["mask_sha256"]
            )
    text = (out / "reviewer" / "cases.json").read_text()
    assert not any(
        word in text for word in ("Process_", "fuscus", "iod", "specimen", str(pilot.cell))
    )
    assert not (out / "reviewer" / "queue.private.json").exists()
    manifest = recovery.load_json(out / "queue.private.json")
    assert (
        manifest["selected_proposals"] == 1
        and manifest["authored_proposals"] == manifest["authored_decisions"] == 0
    )
    assert manifest["policy"] == recovery.load_json(policy)
    assert str(proposal_root / "proposal.json") in manifest["input_pins"]
    draft = pilot.root / "queue-event.json"
    dump(draft, [authored(selected)])
    assert (
        run(
            pilot,
            "decision",
            "--packet",
            root,
            "--event",
            draft,
            "--out",
            pilot.root / "unknown-state",
        )
        == 1
    )
    events = pilot.root / "queue-events"
    assert (
        run(
            pilot,
            "decision",
            "--packet",
            root,
            "--event",
            draft,
            "--proposals",
            proposal_root,
            "--out",
            events,
        )
        == 0
    )
    assert (
        recovery.load_json(events / "events.json")["events"][0]["state_sha256"]
        == proposal["state_sha256"]
    )
    assert run(pilot, "review-queue", "--packet", root, "--out", out) == 1


def test_review_queue_defaults_to_original_states_despite_supplied_proposals(pilot):
    root, value, _, proposal_root, _, _ = queue_proposal(pilot)
    out = pilot.root / "original-queue"
    assert (
        run(pilot, "review-queue", "--packet", root, "--proposals", proposal_root, "--out", out)
        == 0
    )
    queue = recovery.load_json(out / "reviewer" / "cases.json")
    assert [row["state_sha256"] for row in queue] == [
        case["state_sha256"] for case in value["cases"]
    ]
    manifest = recovery.load_json(out / "queue.private.json")
    assert manifest["supplied_proposals"] == 1 and manifest["selected_proposals"] == 0


@pytest.mark.parametrize(
    "change",
    [
        "stale_policy",
        "unknown_case",
        "wrong_case_state",
        "extra_policy_field",
        "stale_proposal",
        "changed_roi",
        "changed_source",
    ],
)
def test_review_queue_refuses_stale_bindings_and_changed_inputs(pilot, change):
    root, value, case, proposal_root, proposal, policy_path = queue_proposal(pilot)
    policy = recovery.load_json(policy_path)
    if change == "stale_policy":
        policy["state_selection"][case["case_id"]] = "0" * 64
    elif change == "unknown_case":
        policy["state_selection"] = {"unknown": case["state_sha256"]}
    elif change == "wrong_case_state":
        policy["state_selection"][value["cases"][1]["case_id"]] = proposal["state_sha256"]
    elif change == "extra_policy_field":
        policy["auto_keep"] = True
    elif change == "stale_proposal":
        proposal["parent_state_sha256"] = "0" * 64
        dump(proposal_root / "proposal.json", proposal)
    elif change == "changed_roi":
        (proposal_root / "roi.zip").write_bytes(b"changed")
    else:
        raw_path = Path(case["source"]["path"])
        raw = tifffile.imread(raw_path)
        raw[-1, -1] -= 1
        tifffile.imwrite(raw_path, raw, metadata={"axes": "YX"})
    dump(policy_path, policy)
    out = pilot.root / "refused-queue"
    assert (
        run(
            pilot,
            "review-queue",
            "--packet",
            root,
            "--proposals",
            proposal_root,
            "--policy",
            policy_path,
            "--out",
            out,
        )
        == 1
    )
    assert not out.exists()


@pytest.mark.parametrize("change", ["overlap", "unsupported", "stale", "crop", "packet", "author"])
def test_proposal_refusal_controls(pilot, change):
    root, value = packet(pilot)
    case = value["cases"][0]
    roi, sidecar = proposal_inputs(pilot, root, case)
    fake_rasterizer(pilot, root, case, overlap=change == "overlap", error=change == "unsupported")
    if change in {"stale", "crop", "packet", "author"}:
        metadata = recovery.load_json(sidecar)
        metadata[
            {
                "stale": "parent_state_sha256",
                "crop": "crop_sha256",
                "packet": "packet_sha256",
                "author": "author",
            }[change]
        ] = "" if change == "author" else "0" * 64
        dump(sidecar, metadata)
    out = pilot.root / "refused"
    assert (
        run(pilot, "import-roi", "--packet", root, "--roi", roi, "--sidecar", sidecar, "--out", out)
        == 1
    )
    assert not (out / "proposal.json").exists()
    if change in {"overlap", "unsupported"}:
        assert recovery.load_json(out / "receipt.json")["exit"] == 1


def test_json_duplicate_keys_and_nonfinite_values_refused(tmp_path):
    path = tmp_path / "bad.json"
    path.write_text('{"case_id":"first","case_id":"second"}')
    with pytest.raises(ValueError, match="Duplicate JSON"):
        recovery.load_json(path)
    path.write_text('{"value":NaN}')
    with pytest.raises(ValueError, match="Nonfinite"):
        recovery.load_json(path)
