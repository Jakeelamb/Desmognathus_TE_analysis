from __future__ import annotations

import importlib.util
import json
import os
import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace

import pytest

ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture
def checker(tmp_path, monkeypatch):
    spec = importlib.util.spec_from_file_location(
        "environment_checker", ROOT / "scripts/check_environment.py"
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    monkeypatch.setattr(module, "ROOT", tmp_path)
    return module


def manifest():
    return {
        "dependencies": {"r-base": "==4.4.3"},
        "feature": {
            "analysis": {
                "dependencies": {"python": "==3.13.11"},
                "pypi-dependencies": {"numpy": "==2.4.2"},
            }
        },
        "environments": {"default": {"features": ["analysis"]}},
    }


@pytest.mark.parametrize(
    "name,version,valid",
    [
        ("default", (3, 13, 11), True),
        ("default", (3, 13, 10), False),
        ("bonsai", (3, 14, 2), True),
        ("bonsai", (3, 13, 11), False),
        ("image-qc", (3, 13, 11), True),
        ("image-qc", (3, 14, 2), False),
        ("recovery", (3, 13, 11), True),
        ("recovery", (3, 14, 2), False),
    ],
)
def test_python_environment_version(checker, monkeypatch, name, version, valid):
    monkeypatch.setattr(checker.sys, "version_info", version)
    if valid:
        checker.check_python(name)
    else:
        with pytest.raises(ValueError, match="Python version"):
            checker.check_python(name)


@pytest.mark.parametrize("problem", ["project", "prefix", "interpreter", "missing"])
def test_runtime_identity_rejects_foreign_and_missing_runtime(checker, monkeypatch, problem):
    prefix = checker.ROOT / ".pixi/envs/default"
    (prefix / "bin").mkdir(parents=True)
    for tool in ("python", "R", "Rscript"):
        executable = prefix / "bin" / tool
        executable.write_text("runtime")
        executable.chmod(0o755)
    monkeypatch.setenv("PIXI_PROJECT_ROOT", str(checker.ROOT))
    monkeypatch.setenv("PIXI_ENVIRONMENT_NAME", "default")
    monkeypatch.setenv("CONDA_PREFIX", str(prefix))
    monkeypatch.setattr(checker.sys, "prefix", str(prefix))
    if problem == "project":
        monkeypatch.setenv("PIXI_PROJECT_ROOT", "/foreign")
    elif problem == "prefix":
        monkeypatch.setenv("CONDA_PREFIX", "/foreign")
    elif problem == "interpreter":
        monkeypatch.setattr(checker.sys, "prefix", "/foreign")
    else:
        (prefix / "bin/python").unlink()
    with pytest.raises(ValueError):
        checker.runtime_prefix()


@pytest.mark.parametrize("problem", [None, "version", "build", "missing", "pypi", "foreign", "extra"])
def test_installed_records_must_match_native_lock(checker, monkeypatch, tmp_path, problem):
    prefix = tmp_path / ".pixi/envs/default"
    metadata = prefix / "conda-meta"
    metadata.mkdir(parents=True)
    locked = [
        {"name": "r-base", "version": "4.4.3", "build": "r-build", "kind": "conda"},
        {"name": "python", "version": "3.13.11", "build": "py-build", "kind": "conda"},
        {"name": "numpy", "version": "2.4.2", "build": "", "kind": "pypi"},
    ]
    for record in locked[:2]:
        actual = dict(record)
        if record["name"] == "r-base":
            if problem == "missing":
                continue
            if problem in {"version", "build"}:
                actual[problem] = "wrong"
        (metadata / f"{record['name']}.json").write_text(json.dumps(actual))
    if problem == "extra":
        (metadata / "foreign.json").write_text(
            json.dumps({"name": "foreign", "version": "1", "build": "foreign"})
        )
    monkeypatch.setattr(checker, "locked_packages", lambda name: locked)
    distribution = SimpleNamespace(
        version="wrong" if problem == "pypi" else "2.4.2",
        locate_file=lambda path: (
            Path("/foreign") if problem == "foreign" else prefix / "lib/site-packages"
        ),
    )
    monkeypatch.setattr(checker.importlib.metadata, "distribution", lambda name: distribution)
    before = {path: path.read_bytes() for path in metadata.iterdir()}
    if problem is None:
        checker.check_packages("default", prefix, manifest())
    else:
        with pytest.raises(ValueError):
            checker.check_packages("default", prefix, manifest())
    assert {path: path.read_bytes() for path in metadata.iterdir()} == before


@pytest.mark.parametrize(
    "declared,valid",
    [("ruamel.yaml", True), ("RuAmEl._--YaMl", True), ("ruamel.yamls", False)],
)
def test_python_package_identity_normalizes_separators(
    checker, monkeypatch, tmp_path, declared, valid
):
    prefix = tmp_path / ".pixi/envs/default"
    locked = [{"name": "ruamel-yaml", "version": "0.19.1", "build": "", "kind": "pypi"}]
    monkeypatch.setattr(checker, "locked_packages", lambda name: locked)
    distribution = SimpleNamespace(
        version="0.19.1", locate_file=lambda path: prefix / "lib/site-packages"
    )
    monkeypatch.setattr(checker.importlib.metadata, "distribution", lambda name: distribution)
    configuration = {
        "feature": {"bonsai": {"pypi-dependencies": {declared: "==0.19.1"}}},
        "environments": {"bonsai": {"features": ["bonsai"]}},
    }
    if valid:
        checker.check_packages("bonsai", prefix, configuration)
    else:
        with pytest.raises(ValueError, match="Manifest dependency absent from Pixi lock"):
            checker.check_packages("bonsai", prefix, configuration)
    assert not prefix.exists()


@pytest.mark.parametrize("pixi_exe", [None, "/project/pixi"])
def test_native_list_is_locked_and_does_not_install(checker, monkeypatch, pixi_exe):
    calls = []
    lookups = []
    if pixi_exe:
        monkeypatch.setenv("PIXI_EXE", pixi_exe)
    else:
        monkeypatch.delenv("PIXI_EXE", raising=False)

    def which(name):
        lookups.append(name)
        return "/fake/pixi"

    monkeypatch.setattr(checker.shutil, "which", which)

    def run(arguments, **kwargs):
        calls.append((arguments, kwargs))
        return SimpleNamespace(stdout='[{"name":"python"}]')

    monkeypatch.setattr(checker.subprocess, "run", run)
    assert checker.locked_packages("default") == [{"name": "python"}]
    assert lookups == [pixi_exe or "pixi"]
    assert calls[0][0] == [
        "/fake/pixi",
        "list",
        "--manifest-path",
        str(checker.ROOT / "pixi.toml"),
        "--environment",
        "default",
        "--json",
        "--locked",
        "--no-install",
    ]


@pytest.mark.parametrize("problem", ["missing", "identity", "symlink"])
def test_private_overlay_is_validated_before_r_runs(checker, monkeypatch, tmp_path, problem):
    overlay = tmp_path / ".pixi/r-library"
    if problem == "identity":
        overlay.mkdir(parents=True)
        (overlay / ".phylopath-source.json").write_text("{}")
    elif problem == "symlink":
        overlay.parent.mkdir()
        overlay.symlink_to(tmp_path, target_is_directory=True)

    def forbidden(*args, **kwargs):
        pytest.fail("R must not run with a missing or damaged overlay")

    monkeypatch.setattr(checker.subprocess, "run", forbidden)
    with pytest.raises(ValueError):
        checker.check_r(tmp_path / "prefix", manifest(), "default", overlay, False)


def test_r_checker_isolates_libraries_and_exercises_frozen_path_reset(
    checker, monkeypatch, tmp_path
):
    overlay = tmp_path / ".pixi/r-library"
    overlay.mkdir(parents=True)
    (overlay / ".phylopath-source.json").write_text(json.dumps(checker.OVERLAY_IDENTITY))
    monkeypatch.setenv("R_LIBS", "/foreign")
    monkeypatch.setenv("R_LIBS_SITE", "/foreign-site")
    calls = []
    monkeypatch.setattr(
        checker.subprocess, "run", lambda args, **kwargs: calls.append((args, kwargs))
    )
    checker.check_r(tmp_path / "prefix", manifest(), "default", overlay, False)
    arguments, options = calls[0]
    assert arguments[:3] == [str(tmp_path / "prefix/bin/Rscript"), "--vanilla", "-"]
    assert "r-base=4.4.3" not in arguments
    assert "base=4.4.3" in arguments
    assert "R_LIBS" not in options["env"]
    assert options["env"]["R_LIBS_SITE"] == str(overlay)
    assert options["env"]["R_PROFILE_USER"] == "/dev/null"
    assert options["env"]["R_PROFILE"] == "/dev/null"
    assert options["env"]["R_ENVIRON"] == "/dev/null"
    assert ".Library.site <-" not in options["input"]
    assert "normalizePath(.Library.site)" in options["input"]
    assert '.libPaths(file.path(prefix, "lib/R/library"))' in options["input"]
    assert "resolve_gbe_font()" in options["input"]


def test_check_env_fails_when_run_with_foreign_python():
    environment = dict(os.environ, PIXI_PROJECT_ROOT=str(ROOT), CONDA_PREFIX="/foreign")
    result = subprocess.run(
        [sys.executable, "-B", str(ROOT / "scripts/check_environment.py")],
        env=environment,
        text=True,
        capture_output=True,
        check=False,
    )
    assert result.returncode == 1
    assert "Wrong runtime prefix" in result.stderr
