"""Read-only checks of this project's installed Pixi packages and R overlay."""

from __future__ import annotations

import argparse
import importlib.metadata
import json
import os
import re
import shutil
import subprocess
import sys
import tomllib
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PHYLOPATH_SHA256 = "fbbfe6ffe78fe4b6d59c98202e3c6e03332d379b4a18aed62008bc6944d5227f"
OVERLAY_IDENTITY = {
    "package": "phylopath",
    "version": "1.3.1",
    "r_version": "4.4.3",
    "source_sha256": PHYLOPATH_SHA256,
}


def runtime_prefix() -> tuple[str, Path]:
    name = os.environ.get("PIXI_ENVIRONMENT_NAME", os.environ.get("PIXI_ENV_NAME", "default"))
    if name not in {"default", "bonsai"}:
        raise ValueError(f"Unsupported Pixi environment: {name}")
    prefix = ROOT / ".pixi/envs" / name
    if Path(os.environ.get("PIXI_PROJECT_ROOT", "")).resolve() != ROOT:
        raise ValueError("Activate this project's Pixi environment before checking it.")
    if Path(os.environ.get("CONDA_PREFIX", "")).resolve() != prefix.resolve():
        raise ValueError(f"Wrong runtime prefix; expected {prefix}")
    if Path(sys.prefix).resolve() != prefix.resolve():
        raise ValueError(f"Wrong Python interpreter; expected {prefix}/bin/python")
    for tool in ("python", "R", "Rscript"):
        if not os.access(prefix / "bin" / tool, os.X_OK):
            raise ValueError(f"Missing Pixi runtime: {prefix}/bin/{tool}; run make setup")
    return name, prefix


def dependencies(manifest: dict, name: str, kind: str) -> dict:
    result = dict(manifest.get(kind, {}))
    environment = manifest["environments"][name]
    features = environment if isinstance(environment, list) else environment["features"]
    for feature in features:
        result.update(manifest["feature"][feature].get(kind, {}))
    return result


def locked_packages(name: str) -> list[dict]:
    pixi = shutil.which("pixi")
    if pixi is None:
        raise ValueError("Pixi executable is required for the read-only lock check.")
    result = subprocess.run(
        [
            pixi,
            "list",
            "--manifest-path",
            str(ROOT / "pixi.toml"),
            "--environment",
            name,
            "--json",
            "--locked",
            "--no-install",
        ],
        cwd=ROOT,
        text=True,
        capture_output=True,
        check=True,
    )
    records = json.loads(result.stdout)
    if not isinstance(records, list) or not records:
        raise ValueError("Pixi returned no locked package records.")
    return records


def normalize_package_name(name: str) -> str:
    return re.sub(r"[-_.]+", "-", name).lower()


def check_packages(name: str, prefix: Path, manifest: dict) -> None:
    installed = {}
    for path in (prefix / "conda-meta").glob("*.json"):
        record = json.loads(path.read_text())
        installed[record["name"]] = record
    records = locked_packages(name)
    locked_conda = {record["name"] for record in records if record["kind"].lower() == "conda"}
    unexpected = set(installed) - locked_conda
    if unexpected:
        raise ValueError(f"Installed Conda packages absent from Pixi lock: {sorted(unexpected)}")
    for record in records:
        package, version = record["name"], record["version"]
        if record["kind"].lower() == "conda":
            actual = installed.get(package)
            if actual is None:
                raise ValueError(f"Missing installed Conda package: {package}")
            if actual["version"] != version or actual["build"] != record["build"]:
                raise ValueError(f"Installed package differs from Pixi lock: {package}")
        elif record["kind"].lower() == "pypi":
            actual = importlib.metadata.distribution(package)
            if actual.version != version:
                raise ValueError(f"Installed Python package differs from Pixi lock: {package}")
            if not Path(actual.locate_file("")).resolve().is_relative_to(prefix.resolve()):
                raise ValueError(f"Python package resolved outside Pixi runtime: {package}")
        else:
            raise ValueError(f"Unknown Pixi package kind: {record['kind']}")
    locked_names = {normalize_package_name(record["name"]) for record in records}
    for kind in ("dependencies", "pypi-dependencies"):
        for package in dependencies(manifest, name, kind):
            if normalize_package_name(package) not in locked_names:
                raise ValueError(f"Manifest dependency absent from Pixi lock: {package}")
    print(f"Installed packages match Pixi lock: PASS ({name})")


def check_python(name: str) -> None:
    observed = ".".join(map(str, sys.version_info[:3]))
    valid = observed == "3.13.11" if name == "default" else sys.version_info[:2] == (3, 14)
    if not valid:
        expected = "3.13.11" if name == "default" else "3.14.x"
        raise ValueError(f"Python version: expected {expected}, found {observed}")
    print(f"Python interpreter: PASS ({observed})")


def check_r(prefix: Path, manifest: dict, name: str, overlay: Path, without: bool) -> None:
    if not without:
        if overlay.is_symlink() or not overlay.is_dir():
            raise ValueError("Missing or unsafe private R overlay; run make setup-r")
        receipt = overlay / ".phylopath-source.json"
        if json.loads(receipt.read_text()) != OVERLAY_IDENTITY:
            raise ValueError("Damaged R overlay identity; existing overlay was not changed")
    expected = {}
    for package, specification in dependencies(manifest, name, "dependencies").items():
        if not package.startswith("r-"):
            continue
        version = specification if isinstance(specification, str) else specification["version"]
        expected[package[2:]] = version.removeprefix("==").removeprefix("=")
    # Pass package pins as arguments; no source interpolation or shell evaluation.
    arguments = [f"{package}={version}" for package, version in expected.items()]
    environment = dict(os.environ)
    for key in (
        "R_LIBS",
        "R_LIBS_SITE",
        "R_HOME",
        "R_ENVIRON",
        "R_PROFILE",
        "R_DEFAULT_PACKAGES",
        "PYTHONPATH",
        "PYTHONHOME",
        "PYTHONUSERBASE",
    ):
        environment.pop(key, None)
    environment.update(
        R_LIBS_USER=str(overlay),
        R_LIBS_SITE=str(overlay),
        R_ENVIRON="/dev/null",
        R_PROFILE="/dev/null",
        R_ENVIRON_USER="/dev/null",
        R_PROFILE_USER="/dev/null",
        CONDA_PREFIX=str(prefix),
    )
    code = r"""
args <- commandArgs(trailingOnly = TRUE)
overlay <- args[[1]]
without <- args[[2]] == "true"
prefix <- Sys.getenv("CONDA_PREFIX")
# Verify R startup honored isolation before exercising the frozen path runner reset.
expected_site <- if (dir.exists(overlay)) normalizePath(overlay) else character()
if (!identical(normalizePath(.Library.site), expected_site)) {
  stop("R startup did not retain only the private site library")
}
.libPaths(c(if (dir.exists(overlay)) overlay, file.path(prefix, "lib/R/library")))
for (pin in args[-c(1, 2, 3)]) {
  fields <- strsplit(pin, "=", fixed = TRUE)[[1]]
  package <- fields[[1]]
  expected <- fields[[2]]
  if (package == "base") {
    if (as.character(getRversion()) != expected) stop("Wrong R version: ", getRversion())
  } else {
    if (!requireNamespace(package, quietly = TRUE)) stop("Missing R package: ", package)
    expected_version <- package_version(gsub("_", "-", expected, fixed = TRUE))
    if (packageVersion(package) != expected_version) {
      stop("Wrong R package version: ", package)
    }
    location <- normalizePath(find.package(package))
    if (!startsWith(location, paste0(normalizePath(file.path(prefix, "lib/R/library")), "/"))) {
      stop("R package outside Pixi runtime: ", package)
    }
  }
}
if (!without) {
  .libPaths(file.path(prefix, "lib/R/library"))
  if (!requireNamespace("phylopath", quietly = TRUE) ||
      as.character(packageVersion("phylopath")) != "1.3.1") stop("phylopath 1.3.1 is required")
  if (normalizePath(find.package("phylopath")) !=
      normalizePath(file.path(overlay, "phylopath"))) stop("phylopath must use private overlay")
}
if (!capabilities("cairo")) stop("R Cairo support is required")
source(args[[3]])
stopifnot(identical(resolve_gbe_font(), "Arial"))
cat("R runtime, packages, Cairo, Arial and private overlay: PASS\n")
"""
    subprocess.run(
        [
            str(prefix / "bin/Rscript"),
            "--vanilla",
            "-",
            str(overlay),
            str(without).lower(),
            str(ROOT / "scripts/publication/gbe_theme.R"),
            *arguments,
        ],
        input=code,
        env=environment,
        cwd=ROOT,
        text=True,
        check=True,
    )
    converters = {
        "coreutils": ("sha256sum", "--version"),
        "ghostscript": ("gs", "--version"),
        "imagemagick": ("magick", "-version"),
        "pandoc": ("pandoc", "--version"),
    }
    for package, (tool, option) in converters.items():
        if package in dependencies(manifest, name, "dependencies"):
            executable = prefix / "bin" / tool
            if not os.access(executable, os.X_OK):
                raise ValueError(f"Missing environment-local converter: {executable}")
            subprocess.run([str(executable), option], stdout=subprocess.DEVNULL, check=True)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--python-only", action="store_true")
    group.add_argument("--r-only", action="store_true")
    parser.add_argument("--without-phylopath", action="store_true")
    parser.add_argument("--overlay", type=Path)
    args = parser.parse_args()
    try:
        name, prefix = runtime_prefix()
        manifest = tomllib.loads((ROOT / "pixi.toml").read_text())
        check_python(name)
        check_packages(name, prefix, manifest)
        if not args.python_only:
            overlay = args.overlay or ROOT / ".pixi/r-library"
            check_r(prefix, manifest, name, overlay, args.without_phylopath)
    except (
        ValueError,
        OSError,
        KeyError,
        json.JSONDecodeError,
        importlib.metadata.PackageNotFoundError,
        subprocess.CalledProcessError,
    ) as error:
        print(f"Environment: FAIL; {error}. Existing environment was not changed.", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
