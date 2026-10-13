from __future__ import annotations

import hashlib
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SHA256 = "fbbfe6ffe78fe4b6d59c98202e3c6e03332d379b4a18aed62008bc6944d5227f"
FAKE_TOOL = r"""
import json, os, subprocess, sys
from pathlib import Path
name = Path(sys.argv[0]).name
args = sys.argv[1:]
entry = dict(tool=name, args=args, cwd=os.getcwd(),
             environment={k: os.environ.get(k) for k in
               ("CONDA_PREFIX", "R_LIBS_USER", "R_LIBS", "R_LIBS_SITE", "R_ENVIRON_USER",
                "R_PROFILE_USER", "R_ENVIRON", "R_PROFILE", "PYTHONPATH", "PYTHONHOME", "VIRTUAL_ENV")})
with open(os.environ["FAKE_LOG"], "a") as log:
    log.write(json.dumps(entry) + "\n")
if name == "pixi":
    assert "--locked" in args and "--no-install" in args and "--executable" in args
    root = Path(args[args.index("--manifest-path") + 1]).parent
    env = dict(os.environ, PIXI_PROJECT_ROOT=str(root),
               CONDA_PREFIX=str(root / ".pixi/envs/default"), PIXI_ENVIRONMENT_NAME="default")
    command = args[args.index("--executable") + 1:]
    sys.exit(subprocess.call(command, env=env))
elif name == "python":
    if os.environ.get("FAKE_CHECK_FAILURE"):
        print("Damaged R overlay; existing environment was not changed", file=sys.stderr)
        sys.exit(1)
    if "--without-phylopath" not in args:
        overlay = Path(args[args.index("--overlay") + 1]) if "--overlay" in args else Path(os.environ["R_LIBS_USER"])
        if not (overlay / "phylopath/DESCRIPTION").is_file() or not (overlay / ".phylopath-source.json").is_file():
            print("Missing private R overlay", file=sys.stderr)
            sys.exit(1)
elif name == "curl":
    if os.environ.get("FAKE_DOWNLOAD_FAILURE"): sys.exit(1)
    Path(args[args.index("-o") + 1]).write_text("pinned tar fixture")
elif name == "sha256sum":
    print(os.environ.get("FAKE_SHA256", "wrong") + "  " + args[0])
elif name == "R":
    if os.environ.get("FAKE_INSTALL_FAILURE"): sys.exit(1)
    library = Path(next(a.split("=", 1)[1] for a in args if a.startswith("--library=")))
    (library / "phylopath").mkdir()
    (library / "phylopath/DESCRIPTION").write_text("Version: 1.3.1\n")
"""


class Runtime:
    def __init__(self, directory: Path):
        self.root = directory / "project with spaces"
        (self.root / "scripts").mkdir(parents=True)
        for name in ("activate_pixi.sh", "run_r.sh", "setup_r.sh"):
            shutil.copy2(ROOT / "scripts" / name, self.root / "scripts" / name)
        self.prefix = self.root / ".pixi/envs/default"
        (self.prefix / "bin").mkdir(parents=True)
        self.overlay = self.root / ".pixi/r-library"
        self.tools = directory / "tools"
        self.tools.mkdir()
        self.log = directory / "calls.jsonl"
        for name in ("R", "Rscript", "python"):
            self.tool(self.prefix / "bin" / name)
        for name in ("pixi", "curl", "sha256sum"):
            self.tool(self.tools / name)
        self.env = {
            key: value
            for key, value in os.environ.items()
            if not key.startswith(("PIXI_", "CONDA_", "R_LIBS", "PYTHON"))
        }
        self.env.update(
            PATH=f"{self.tools}:/usr/bin:/bin",
            FAKE_LOG=str(self.log),
            FAKE_SHA256=SHA256,
            R_LIBS="/foreign",
            R_LIBS_SITE="/foreign-site",
            PYTHONPATH="/foreign-python",
            PYTHONHOME="/foreign-home",
            VIRTUAL_ENV="/foreign-venv",
        )
        self.outside = directory / "outside"
        self.outside.mkdir()

    @staticmethod
    def tool(path: Path):
        path.write_text(f"#!{sys.executable}\n{FAKE_TOOL}")
        path.chmod(0o755)

    def activate(self):
        self.env.update(
            PIXI_PROJECT_ROOT=str(self.root),
            CONDA_PREFIX=str(self.prefix),
            PIXI_ENVIRONMENT_NAME="default",
        )

    def installed_overlay(self):
        (self.overlay / "phylopath").mkdir(parents=True)
        (self.overlay / "phylopath/DESCRIPTION").write_text("Version: 1.3.1\n")
        (self.overlay / ".phylopath-source.json").write_text("{}")

    def run(self, script: str, *arguments: str):
        return subprocess.run(
            ["bash", str(self.root / "scripts" / script), *arguments],
            cwd=self.outside,
            env=self.env,
            text=True,
            capture_output=True,
            check=False,
        )

    def calls(self):
        return (
            [json.loads(line) for line in self.log.read_text().splitlines()]
            if self.log.exists()
            else []
        )


def snapshot(root: Path):
    return {
        str(path.relative_to(root)): hashlib.sha256(path.read_bytes()).hexdigest()
        for path in root.rglob("*")
        if path.is_file()
    }


@pytest.fixture
def runtime(tmp_path):
    return Runtime(tmp_path)


@pytest.mark.parametrize(
    "activated,ambient_conda", [(False, None), (False, "/home/jake/miniconda3"), (True, None)]
)
def test_runner_preserves_arguments_root_and_isolation(runtime, activated, ambient_conda):
    if ambient_conda:
        runtime.env["CONDA_PREFIX"] = ambient_conda
    if activated:
        runtime.activate()
    arguments = ["--no-environ", "script with spaces.R", "--output", "'quoted' $literal;echo no"]
    result = runtime.run("run_r.sh", "Rscript", *arguments)
    assert result.returncode == 0, result.stderr
    calls = runtime.calls()
    assert [call["tool"] for call in calls] == (["Rscript"] if activated else ["pixi", "Rscript"])
    call = calls[-1]
    assert call["args"] == arguments
    assert call["cwd"] == str(runtime.root)
    assert call["environment"] == {
        "CONDA_PREFIX": str(runtime.prefix),
        "R_LIBS_USER": str(runtime.overlay),
        "R_LIBS": None,
        "R_LIBS_SITE": str(runtime.overlay),
        "R_ENVIRON_USER": "/dev/null",
        "R_PROFILE_USER": "/dev/null",
        "R_ENVIRON": "/dev/null",
        "R_PROFILE": "/dev/null",
        "PYTHONPATH": None,
        "PYTHONHOME": None,
        "VIRTUAL_ENV": None,
    }


@pytest.mark.parametrize(
    "script,args",
    [("run_r.sh", ["Rscript", "a.R"]), ("setup_r.sh", []), ("setup_r.sh", ["--check"])],
)
@pytest.mark.parametrize("problem", ["missing", "foreign-prefix", "foreign-project"])
def test_runtime_guard_never_executes_or_provisions(runtime, script, args, problem):
    runtime.activate()
    if problem == "missing":
        (runtime.prefix / "bin/Rscript").unlink()
    elif problem == "foreign-prefix":
        runtime.env["CONDA_PREFIX"] = "/foreign"
    else:
        runtime.env["PIXI_PROJECT_ROOT"] = "/foreign"
    before = snapshot(runtime.root)
    result = runtime.run(script, *args)
    assert result.returncode == 1
    assert runtime.calls() == []
    assert snapshot(runtime.root) == before


@pytest.mark.parametrize("args", [[], ["--check"]])
def test_existing_setup_is_read_only(runtime, args):
    runtime.activate()
    runtime.installed_overlay()
    before = snapshot(runtime.root)
    result = runtime.run("setup_r.sh", *args)
    assert result.returncode == 0, result.stderr
    assert snapshot(runtime.root) == before
    assert [call["tool"] for call in runtime.calls()] == ["python"]


@pytest.mark.parametrize("args", [[], ["--check"]])
def test_damaged_existing_overlay_is_never_repaired(runtime, args):
    runtime.activate()
    runtime.installed_overlay()
    runtime.env["FAKE_CHECK_FAILURE"] = "1"
    before = snapshot(runtime.root)
    result = runtime.run("setup_r.sh", *args)
    assert result.returncode == 1
    assert "existing environment was not changed" in result.stderr
    assert snapshot(runtime.root) == before
    assert [call["tool"] for call in runtime.calls()] == ["python"]


def test_missing_overlay_check_is_read_only(runtime):
    runtime.activate()
    before = snapshot(runtime.root)
    result = runtime.run("setup_r.sh", "--check")
    assert result.returncode == 1
    assert snapshot(runtime.root) == before
    assert not runtime.overlay.exists()


def test_setup_requires_activated_project(runtime):
    result = runtime.run("setup_r.sh")
    assert result.returncode == 1
    assert "inside this project Pixi environment" in result.stderr
    assert runtime.calls() == []


def test_setup_stages_verifies_and_installs_only_overlay_once(runtime):
    runtime.activate()
    before = snapshot(runtime.prefix)
    first = runtime.run("setup_r.sh")
    assert first.returncode == 0, first.stderr
    overlay_before = snapshot(runtime.overlay)
    second = runtime.run("setup_r.sh")
    assert second.returncode == 0, second.stderr
    assert snapshot(runtime.prefix) == before
    assert snapshot(runtime.overlay) == overlay_before
    calls = runtime.calls()
    installs = [call for call in calls if call["tool"] == "R"]
    assert len(installs) == 1
    library = next(arg for arg in installs[0]["args"] if arg.startswith("--library="))
    assert "phylopath-setup." in library
    assert (
        json.loads((runtime.overlay / ".phylopath-source.json").read_text())["source_sha256"]
        == SHA256
    )
    assert not list((runtime.root / ".pixi").glob("phylopath-setup.*"))
    assert not (runtime.root / ".pixi/r-overlay.setup.lock").exists()


@pytest.mark.parametrize("failure", ["DOWNLOAD", "SHA256", "INSTALL", "CHECK"])
def test_failed_setup_never_changes_runtime_or_exposes_overlay(runtime, failure):
    runtime.activate()
    runtime.env[f"FAKE_{failure}_FAILURE" if failure != "SHA256" else "FAKE_SHA256"] = "bad"
    before = snapshot(runtime.prefix)
    result = runtime.run("setup_r.sh")
    assert result.returncode == 1
    assert snapshot(runtime.prefix) == before
    assert not runtime.overlay.exists()
    assert not list((runtime.root / ".pixi").glob("phylopath-setup.*"))
    if failure in {"DOWNLOAD", "SHA256", "CHECK"}:
        assert not any(call["tool"] == "R" for call in runtime.calls())


def test_other_setup_lock_is_not_removed(runtime):
    runtime.activate()
    lock = runtime.root / ".pixi/r-overlay.setup.lock"
    lock.mkdir()
    result = runtime.run("setup_r.sh")
    assert result.returncode == 1
    assert lock.is_dir()
    assert not any(call["tool"] == "R" for call in runtime.calls())


@pytest.mark.parametrize(
    "script,args", [("run_r.sh", ["R", "a.R"]), ("setup_r.sh", ["--check", "extra"])]
)
def test_invalid_usage_executes_nothing(runtime, script, args):
    result = runtime.run(script, *args)
    assert result.returncode == 2
    assert "Usage:" in result.stderr
    assert runtime.calls() == []
