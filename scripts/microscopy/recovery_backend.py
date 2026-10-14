"""Pinned native measurement and ImageJ ROI exchange; no scientific decisions."""

from __future__ import annotations

import hashlib
import json
import math
import os
import re
import shutil
import subprocess
import sys
import uuid
import zipfile
from dataclasses import dataclass
from datetime import UTC, datetime
from pathlib import Path

import numpy as np

IMAGEJ_URL = "https://repo.maven.apache.org/maven2/net/imagej/ij/1.54p/ij-1.54p.jar"
IMAGEJ_SHA256 = "2e1a09961dfb41cee66ddc821b2577a41a072566ce45a49bae69267099741e20"
SOURCE_PATHS = (
    "scripts/microscopy/recovery_backend.py",
    "scripts/microscopy/setup_recovery_tools.py",
    "scripts/microscopy/src/microscopy_recovery.rs",
    "scripts/microscopy/src/RecoveryRoiBridge.java",
)
FIELDS = {
    "count",
    "raw_sum",
    "min_pixel",
    "max_pixel",
    "zero_pixels",
    "one_pixels",
    "saturated_pixels",
    "background",
    "iod",
    "mean_od",
    "median_pixel",
}


@dataclass(frozen=True)
class Toolchain:
    manifest_path: Path
    java: Path
    ij_jar: Path
    roi_class_dir: Path
    reducer: Path
    identity: str


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def json_bytes(value: object) -> bytes:
    return (json.dumps(value, sort_keys=True, indent=2, allow_nan=False) + "\n").encode()


def no_symlink(path: Path) -> None:
    for component in (path, *path.parents):
        if component.is_symlink():
            raise ValueError(f"Symlink path refused: {component}")


def active_prefix(repo_root: Path) -> Path:
    root = repo_root.resolve()
    prefix = Path(os.environ.get("CONDA_PREFIX", "/missing-pixi-prefix")).resolve()
    if (
        os.environ.get("PIXI_ENVIRONMENT_NAME", os.environ.get("PIXI_ENV_NAME")) != "recovery"
        or Path(os.environ.get("PIXI_PROJECT_ROOT", "/missing-pixi-project")).resolve() != root
        or prefix != (root / ".pixi/envs/recovery").resolve()
        or not Path(sys.executable).resolve().is_relative_to(prefix)
    ):
        raise ValueError("Run through the project's locked recovery Pixi environment")
    return prefix


def tool_inputs(repo_root: Path) -> dict:
    """Collect source, lock, runtime and compiler pins without acquiring/building."""
    root = repo_root.resolve()
    prefix = active_prefix(root)
    tools = {}
    for name in (
        "java",
        "javac",
        "rustc",
        "x86_64-conda-linux-gnu-gcc",
        "x86_64-conda-linux-gnu-ld",
    ):
        found = shutil.which(name)
        if found is None or not Path(found).resolve().is_relative_to(prefix):
            raise ValueError(f"{name} must belong to the active recovery prefix")
        tools[name] = str(Path(found).resolve())
    java_home = Path(tools["java"]).parent.parent
    runtime = [java_home / "release", java_home / "lib/modules", java_home / "lib/server/libjvm.so"]
    # Shared launch/compiler libraries participate in the runtime identity.
    runtime += sorted((java_home / "lib").glob("*.so"))
    runtime += sorted((prefix / "lib").glob("librustc_driver*.so"))
    pins = {str(path): sha256(path) for path in [*map(Path, tools.values()), *runtime]}
    sources = {name: sha256(root / name) for name in SOURCE_PATHS}
    return {
        "schema": 1,
        "sources": sources,
        "pixi_lock_sha256": sha256(root / "pixi.lock"),
        "pixi_manifest_sha256": sha256(root / "pixi.toml"),
        "tools": tools,
        "runtime_pins": pins,
        "imagej_sha256": IMAGEJ_SHA256,
        "imagej_url": IMAGEJ_URL,
    }


def tool_identity(inputs: dict) -> str:
    return hashlib.sha256(json_bytes(inputs)).hexdigest()


def load_toolchain(repo_root: Path) -> Toolchain:
    inputs = tool_inputs(repo_root)
    identity = tool_identity(inputs)
    cache = repo_root.resolve() / ".pixi/recovery-tools" / identity
    no_symlink(cache)
    manifest_path = cache / "manifest.json"
    if not manifest_path.is_file():
        raise ValueError("Recovery toolchain missing or changed; run explicit recovery-setup")
    manifest = json.loads(manifest_path.read_text())
    if (
        manifest.get("status") != "complete"
        or manifest.get("inputs") != inputs
        or manifest.get("identity") != identity
    ):
        raise ValueError("Recovery toolchain manifest mismatch")
    expected = {"ij.jar", "microscopy_recovery", "classes/RecoveryRoiBridge.class"}
    artifacts = manifest.get("artifacts", {})
    actual_classes = {
        str(p.relative_to(cache)) for p in (cache / "classes").rglob("*") if p.is_file()
    }
    if set(artifacts) != expected or actual_classes != {"classes/RecoveryRoiBridge.class"}:
        raise ValueError("Recovery toolchain artifact set mismatch")
    for name, digest in artifacts.items():
        path = cache / name
        no_symlink(path)
        if not path.is_file() or sha256(path) != digest:
            raise ValueError(f"Recovery toolchain artifact changed: {name}")
    if artifacts["ij.jar"] != IMAGEJ_SHA256:
        raise ValueError("ImageJ checksum mismatch")
    qualification_path = cache / "qualification.json"
    no_symlink(qualification_path)
    if not qualification_path.is_file() or sha256(qualification_path) != manifest.get(
        "qualification_sha256"
    ):
        raise ValueError("Recovery native qualification missing or changed")
    qualification = json.loads(qualification_path.read_text())
    if qualification.get("status") != "passed" or qualification.get("fixture_only") is not True:
        raise ValueError("Recovery native qualification did not pass")
    return Toolchain(
        manifest_path,
        Path(inputs["tools"]["java"]),
        cache / "ij.jar",
        cache / "classes",
        cache / "microscopy_recovery",
        identity,
    )


def clean_native_env() -> dict[str, str]:
    env = os.environ.copy()
    for name in (
        "JAVA_TOOL_OPTIONS",
        "_JAVA_OPTIONS",
        "JDK_JAVA_OPTIONS",
        "CLASSPATH",
        "LD_PRELOAD",
        "LD_LIBRARY_PATH",
        "GCC_EXEC_PREFIX",
        "COMPILER_PATH",
        "LIBRARY_PATH",
        "CPATH",
        "C_INCLUDE_PATH",
        "CPLUS_INCLUDE_PATH",
    ):
        env.pop(name, None)
    return env


def retained_run(
    argv: list[str], evidence: Path, label: str, *, data: bytes | None = None, timeout: int = 120
) -> subprocess.CompletedProcess:
    """Retain clocks, actual exit and raw streams even for a rejected command."""
    no_symlink(evidence)
    evidence.mkdir(parents=True, exist_ok=True)
    stem = evidence / f"{label}-{uuid.uuid4().hex}"
    receipt = {
        "argv": argv,
        "start_utc": datetime.now(UTC).isoformat(),
        "stdin_sha256": hashlib.sha256(data).hexdigest() if data is not None else None,
    }
    try:
        result = subprocess.run(
            argv,
            input=data,
            capture_output=True,
            timeout=timeout,
            env=clean_native_env(),
            check=False,
        )
        receipt["exit_code"] = result.returncode
        stdout, stderr = result.stdout, result.stderr
    except (OSError, subprocess.TimeoutExpired) as error:
        receipt.update(exit_code=None, error=repr(error))
        stdout = getattr(error, "stdout", None) or b""
        stderr = getattr(error, "stderr", None) or b""
        result = None
    receipt["end_utc"] = datetime.now(UTC).isoformat()
    stdout_path = stem.with_suffix(".stdout.log")
    stderr_path = stem.with_suffix(".stderr.log")
    stdout_path.write_bytes(stdout)
    stderr_path.write_bytes(stderr)
    receipt.update(
        stdout_path=str(stdout_path),
        stderr_path=str(stderr_path),
        stdout_sha256=sha256(stdout_path),
        stderr_sha256=sha256(stderr_path),
    )
    stem.with_suffix(".receipt.json").write_bytes(json_bytes(receipt))
    if result is None:
        raise ValueError(f"Native command failed to execute; retained receipt: {stem}")
    return result


def measure_pixels(pixels: np.ndarray, background: float, tools: Toolchain) -> dict:
    if (
        not isinstance(pixels, np.ndarray)
        or pixels.dtype != np.uint8
        or pixels.ndim != 1
        or not 0 < pixels.size <= 100_000_000
    ):
        raise ValueError("Nonempty bounded one-dimensional uint8 vector required")
    if isinstance(background, (bool, np.bool_)) or not isinstance(
        background, (int, float, np.integer, np.floating)
    ):
        raise ValueError("Background must be finite in [1,255]")  # noqa: TRY004 -- one public boundary refusal type
    background = float(background)
    if not math.isfinite(background) or not 1 <= background <= 255:
        raise ValueError("Background must be finite in [1,255]")
    manifest = json.loads(tools.manifest_path.read_text())
    if (
        sha256(tools.reducer) != manifest["artifacts"]["microscopy_recovery"]
        or manifest["identity"] != tools.identity
    ):
        raise ValueError("Reducer changed after toolchain loading")
    result = retained_run(
        [str(tools.reducer), repr(background)],
        tools.manifest_path.parent / "runs",
        "measure",
        data=pixels.tobytes(order="C"),
    )
    if result.returncode != 0:
        raise ValueError(f"Native reducer rejected input: {result.stderr.decode(errors='replace')}")
    value = json.loads(result.stdout)
    if (
        set(value) != FIELDS
        or value["count"] != pixels.size
        or value["background"] != background
        or not all(math.isfinite(v) for v in value.values())
    ):
        raise ValueError("Invalid native measurement response")
    return value


def java_argv(tools: Toolchain, home: Path) -> list[str]:
    return [
        str(tools.java),
        "-Djava.awt.headless=true",
        f"-Duser.home={home}",
        f"-Djava.util.prefs.userRoot={home / 'prefs'}",
        "-cp",
        os.pathsep.join([str(tools.roi_class_dir), str(tools.ij_jar)]),
        "RecoveryRoiBridge",
    ]


def rasterize_rois(
    roi_zip: Path, width: int, height: int, output_dir: Path, tools: Toolchain
) -> dict[str, Path]:
    if (
        type(width) is not int
        or type(height) is not int
        or width <= 0
        or height <= 0
        or width * height > 100_000_000
    ):
        raise ValueError("Bounded positive integer crop dimensions required")
    roi_zip = Path(roi_zip).absolute()
    output_dir = Path(output_dir).absolute()
    no_symlink(roi_zip)
    no_symlink(output_dir)
    roi_zip = roi_zip.resolve()
    output_dir = output_dir.resolve()
    if not roi_zip.is_file() or roi_zip.stat().st_size > 20 * 1024 * 1024:
        raise ValueError("Bounded regular ROI archive required")
    if output_dir.exists() or roi_zip.is_relative_to(output_dir):
        raise ValueError("Fresh destination outside input ancestry required")
    repo = Path(__file__).resolve().parents[2]
    if output_dir.is_relative_to(repo) and not output_dir.is_relative_to(repo / ".pixi/runs"):
        raise ValueError("Project output must be under .pixi/runs; source paths are protected")
    before = sha256(roi_zip)
    roles = []
    try:
        with zipfile.ZipFile(roi_zip) as archive:
            entries = archive.infolist()
            total = 0
            for entry in entries:
                name = entry.filename
                if (
                    entry.is_dir()
                    or (
                        name != "nucleus.roi"
                        and re.fullmatch(r"background_[0-9]{3}\.roi", name) is None
                    )
                    or name == "background_000.roi"
                    or name in roles
                    or entry.flag_bits & 1
                    or entry.file_size < 64
                    or entry.file_size > 1048576
                ):
                    raise ValueError(
                        f"Invalid, duplicate, encrypted or oversized ROI entry: {name}"
                    )
                mode = entry.external_attr >> 16
                if mode & 0o170000 not in (0, 0o100000):
                    raise ValueError("Non-regular ZIP entry refused")
                roles.append(name)
                total += entry.file_size
            if not 1 <= len(roles) <= 64 or total > 16777216 or "nucleus.roi" not in roles:
                raise ValueError(
                    "Exactly one nucleus.roi and bounded optional backgrounds required"
                )
    except zipfile.BadZipFile as error:
        raise ValueError("Malformed ROI ZIP") from error
    output_dir.mkdir(parents=True, exist_ok=False)
    home = output_dir / "private_home"
    home.mkdir()
    result = retained_run(
        java_argv(tools, home) + ["raster", str(roi_zip), str(output_dir), str(width), str(height)],
        output_dir / "evidence",
        "raster",
    )
    if sha256(roi_zip) != before:
        raise ValueError("ROI input changed during rasterization")
    if result.returncode != 0:
        raise ValueError(
            f"ImageJ ROI rejection; retained failed staging: {output_dir}: {result.stderr.decode(errors='replace')}"
        )
    from PIL import Image

    paths = {}
    for entry in roles:
        role = entry[:-4]
        path = output_dir / f"{role}.png"
        with Image.open(path) as image:
            array = np.asarray(image)
            if (
                image.mode != "L"
                or array.dtype != np.uint8
                or array.shape != (height, width)
                or not np.isin(array, [0, 255]).all()
                or not array.any()
            ):
                raise ValueError("Noncanonical or empty ImageJ PNG mask")
        paths[role] = path
    manifest = {
        "status": "complete",
        "toolchain": tools.identity,
        "roi_zip_sha256": before,
        "width": width,
        "height": height,
        "masks": {role: sha256(path) for role, path in paths.items()},
    }
    (output_dir / "rasterization.json").write_bytes(json_bytes(manifest))
    return paths
