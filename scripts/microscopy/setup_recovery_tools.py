"""Explicit locked-prefix acquisition, compilation and qualification of recovery tools."""

from __future__ import annotations

import argparse
import json
import math
import re
import urllib.request
from datetime import UTC, datetime
from pathlib import Path

import numpy as np
from PIL import Image
from recovery_backend import (
    IMAGEJ_SHA256,
    IMAGEJ_URL,
    Toolchain,
    java_argv,
    json_bytes,
    load_toolchain,
    no_symlink,
    retained_run,
    sha256,
    tool_identity,
    tool_inputs,
)


def check_run(
    argv: list[str], evidence: Path, label: str, *, expected: int = 0, data: bytes | None = None
):
    result = retained_run(argv, evidence, label, data=data)
    if result.returncode != expected:
        raise ValueError(f"{label} exit {result.returncode}, expected {expected}; see {evidence}")
    return result


def qualify(tools: Toolchain, cache: Path) -> dict:
    evidence = cache / "evidence"
    home = cache / "private_home"
    home.mkdir()
    base = java_argv(tools, home)
    fixture = cache / "native_fixture"
    check_run(base + ["fixture", str(fixture)], evidence, "imagej_fixture")
    expected = {
        "asymmetric_integer": (67, 414876),
        "concave_freehand": (73, 523940),
        "subpixel_polygon": (65, 439827),
        "composite_hole": (118, 887123),
        "boundary_rectangle": (48, 295152),
    }
    rows = {}
    for line in (fixture / "fixture.tsv").read_text().splitlines():
        name, count, raw_sum, imagej_raw = line.split("\t")
        rows[name] = (int(count), int(raw_sum))
        if float(imagej_raw) != int(raw_sum):
            raise ValueError("ImageJ RawIntDen equality failed")
    if rows != expected:
        raise ValueError("Native ImageJ fixture differs from retained five-geometry reference")
    for name in expected:
        out = cache / f"roundtrip_{name}"
        out.mkdir()
        check_run(
            base + ["raster", str(fixture / f"{name}.zip"), str(out), "16", "12"],
            evidence,
            f"roundtrip_{name}",
        )
        with Image.open(fixture / f"{name}.png") as image:
            original = np.asarray(image)
        with Image.open(out / "nucleus.png") as image:
            actual = np.asarray(image)
        if (
            original.dtype != np.uint8
            or not np.array_equal(original, actual)
            or not np.isin(actual, [0, 255]).all()
        ):
            raise ValueError("Canonical full-grid PNG roundtrip failed")
        if name == "composite_hole" and actual[4:7, 5:9].any():
            raise ValueError("Composite hole filled")
    for archive in sorted(fixture.glob("reject_*.zip")):
        out = cache / archive.stem
        out.mkdir()
        check_run(
            base + ["raster", str(archive), str(out), "16", "12"],
            evidence,
            archive.stem,
            expected=1,
        )
    if (fixture / "original_before.tif").read_bytes() != (
        fixture / "original_after.tif"
    ).read_bytes():
        raise ValueError("Synthetic source TIFF changed")
    pixels = bytes([0, 1, 218, 255, 255, 100, 250])
    result = check_run([str(tools.reducer), "218"], evidence, "rust_signed_floor", data=pixels)
    value = json.loads(result.stdout)
    iod = math.fsum(math.log10(218 / max(pixel, 1)) for pixel in pixels)
    if (
        value["count"] != len(pixels)
        or value["raw_sum"] != sum(pixels)
        or value["zero_pixels"] != 1
        or value["one_pixels"] != 1
        or value["saturated_pixels"] != 2
        or value["median_pixel"] != 218
        or not math.isclose(value["iod"], iod, rel_tol=0, abs_tol=1e-12)
    ):
        raise ValueError("Rust signed/floor reducer fixture failed")
    negative = json.loads(
        check_run([str(tools.reducer), "218"], evidence, "rust_negative", data=bytes([255])).stdout
    )
    if negative["iod"] >= 0:
        raise ValueError("Negative OD was dropped")
    check_run(
        [str(tools.reducer), "0"], evidence, "rust_background_reject", expected=2, data=b"\x01"
    )
    check_run([str(tools.reducer), "255"], evidence, "rust_empty_reject", expected=2, data=b"")
    summary = {
        "status": "passed",
        "fixture_only": True,
        "geometry_reference": expected,
        "native_java": str(tools.java),
        "imagej_sha256": sha256(tools.ij_jar),
        "composite_hole_preserved": True,
        "source_tiff_bytes_unchanged": True,
        "signed_floor_reducer": value,
        "negative_od": negative,
        "negative_geometry_controls": len(list(fixture.glob("reject_*.zip"))),
    }
    (cache / "qualification.json").write_bytes(json_bytes(summary))
    return summary


def setup(repo_root: Path) -> Toolchain:
    root = repo_root.resolve()
    inputs = tool_inputs(root)
    identity = tool_identity(inputs)
    cache = root / ".pixi/recovery-tools" / identity
    no_symlink(cache)
    if cache.exists():
        # An interrupted or rejected setup remains inspectable and is never repaired in place.
        return load_toolchain(root)
    cache.mkdir(parents=True, exist_ok=False)
    evidence = cache / "evidence"
    status = {
        "status": "building",
        "identity": identity,
        "inputs": inputs,
        "start_utc": datetime.now(UTC).isoformat(),
    }
    (cache / "setup_status.json").write_bytes(json_bytes(status))
    try:
        compiler = check_run([inputs["tools"]["rustc"], "--version"], evidence, "rustc_version")
        if not compiler.stdout.decode().startswith("rustc 1.99.0 "):
            raise ValueError("Locked Rust 1.99.0 required")
        runtime = check_run([inputs["tools"]["java"], "-version"], evidence, "java_version")
        if re.search(r'openjdk version "21\.0\.10(?:-internal)?"', runtime.stderr.decode()) is None:
            raise ValueError("Locked OpenJDK 21.0.10 required")
        javac = check_run([inputs["tools"]["javac"], "-version"], evidence, "javac_version")
        if "javac 21.0.10" not in (javac.stdout + javac.stderr).decode():
            raise ValueError("Locked javac 21.0.10 required")
        acquisition = {
            "url": IMAGEJ_URL,
            "expected_sha256": IMAGEJ_SHA256,
            "start_utc": datetime.now(UTC).isoformat(),
        }
        try:
            with (
                urllib.request.urlopen(IMAGEJ_URL, timeout=60) as response,
                (cache / "ij.jar").open("xb") as out,
            ):
                acquisition.update(final_url=response.url, http_status=response.status)
                total = 0
                while block := response.read(1024 * 1024):
                    total += len(block)
                    if total > 8 * 1024 * 1024:
                        raise ValueError("ImageJ download exceeded byte bound")
                    out.write(block)
            acquisition["sha256"] = sha256(cache / "ij.jar")
            if acquisition["sha256"] != IMAGEJ_SHA256:
                raise ValueError("Downloaded ImageJ checksum mismatch")
            acquisition["status"] = "verified"
        except Exception as error:
            acquisition.update(status="failed", error=repr(error))
            raise
        finally:
            acquisition["end_utc"] = datetime.now(UTC).isoformat()
            (cache / "acquisition.json").write_bytes(json_bytes(acquisition))
        classes = cache / "classes"
        classes.mkdir()
        check_run(
            [
                inputs["tools"]["javac"],
                "--release",
                "21",
                "-classpath",
                str(cache / "ij.jar"),
                "-d",
                str(classes),
                str(root / "scripts/microscopy/src/RecoveryRoiBridge.java"),
            ],
            evidence,
            "compile_java",
        )
        reducer = cache / "microscopy_recovery"
        check_run(
            [
                inputs["tools"]["rustc"],
                "--edition=2024",
                "-C",
                "opt-level=2",
                "-C",
                f"linker={inputs['tools']['x86_64-conda-linux-gnu-gcc']}",
                "-o",
                str(reducer),
                str(root / "scripts/microscopy/src/microscopy_recovery.rs"),
            ],
            evidence,
            "compile_rust",
        )
        tools = Toolchain(
            cache / "manifest.json",
            Path(inputs["tools"]["java"]),
            cache / "ij.jar",
            classes,
            reducer,
            identity,
        )
        qualify(tools, cache)
        if tool_inputs(root) != inputs:
            raise ValueError("Sources, lock or runtime changed during setup")
        artifacts = {
            name: sha256(cache / name)
            for name in ("ij.jar", "microscopy_recovery", "classes/RecoveryRoiBridge.class")
        }
        manifest = {
            "status": "complete",
            "identity": identity,
            "inputs": inputs,
            "artifacts": artifacts,
            "qualification_sha256": sha256(cache / "qualification.json"),
            "completed_utc": datetime.now(UTC).isoformat(),
        }
        with (cache / "manifest.json").open("xb") as out:
            out.write(json_bytes(manifest))
        status["status"] = "complete"
    except Exception as error:
        status.update(status="failed", error=repr(error))
        raise
    finally:
        status["end_utc"] = datetime.now(UTC).isoformat()
        (cache / "setup_status.json").write_bytes(json_bytes(status))
    return load_toolchain(root)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", type=Path, default=Path(__file__).resolve().parents[2])
    args = parser.parse_args()
    tools = setup(args.repo_root)
    print(
        json.dumps(
            {"status": "complete", "identity": tools.identity, "manifest": str(tools.manifest_path)}
        )
    )


if __name__ == "__main__":
    main()
