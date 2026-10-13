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
PHYLOPATH_SHA256 = "fbbfe6ffe78fe4b6d59c98202e3c6e03332d379b4a18aed62008bc6944d5227f"

FAKE_TOOL = r"""
import json
import os
import shutil
import sys
from pathlib import Path

tool = Path(sys.argv[0]).name
args = sys.argv[1:]
stdin = sys.stdin.read() if tool == "Rscript" else ""
with open(os.environ["FAKE_LOG"], "a") as log:
    log.write(json.dumps({"tool": tool, "args": args, "stdin": stdin,
                          "cwd": os.getcwd(), "prefix": os.environ.get("CONDA_PREFIX"),
                          "r_environ": os.environ.get("R_ENVIRON_USER"),
                          "r_profile": os.environ.get("R_PROFILE_USER"),
                          "r_libs": os.environ.get("R_LIBS"),
                          "r_libs_user": os.environ.get("R_LIBS_USER")}) + "\n")
registered = Path(os.environ["FAKE_REGISTERED"])
if tool == "conda":
    if os.environ.get("FAKE_CONDA_FAILURE"):
        print("Conda fixture error", file=sys.stderr)
        sys.exit(1)
    if args == ["config", "--show", "envs_dirs"]:
        print("envs_dirs:")
        print("  - " + os.environ["FAKE_ENV_DIR"])
        for prefix in json.loads(registered.read_text()):
            print("  - " + str(Path(prefix).parent))
    elif args[:2] == ["env", "create"] or args[:2] == ["create", "--yes"]:
        if os.environ.get("FAKE_CREATE_FAILURE"):
            print("Conda create fixture error", file=sys.stderr)
            sys.exit(1)
        if "--prefix" in args:
            prefix = Path(args[args.index("--prefix") + 1])
        else:
            prefix = Path(os.environ["FAKE_ENV_DIR"]) / args[args.index("--name") + 1]
        if prefix.exists():
            print("Refusing to overwrite existing prefix", file=sys.stderr)
            sys.exit(1)
        spec = Path(args[args.index("--file") + 1]).read_text()
        if args[0] == "create" and "@EXPLICIT" not in spec:
            print("Expected explicit package lock", file=sys.stderr)
            sys.exit(99)
        if args[0] == "env" and not spec.startswith("name:"):
            print("Expected YAML environment specification", file=sys.stderr)
            sys.exit(99)
        shutil.copytree(os.environ["FAKE_TEMPLATE"], prefix)
        (prefix / "phylopath-installed").unlink()
        registered.write_text(json.dumps([*json.loads(registered.read_text()), str(prefix)]))
        Path(os.environ["FAKE_REGISTRY"]).write_text(str(prefix) + "\n")
    else:
        print("Unexpected Conda mutation", file=sys.stderr)
        sys.exit(99)
elif tool == "uname":
    print(os.environ["FAKE_OS"] if args == ["-s"] else os.environ["FAKE_ARCH"])
elif tool == "curl":
    if os.environ.get("FAKE_DOWNLOAD_FAILURE"):
        sys.exit(1)
    Path(args[args.index("-o") + 1]).write_text("fake pinned source")
elif tool == "sha256sum":
    print(os.environ["FAKE_SHA256"] + "  " + args[0])
elif tool == "R":
    if os.environ.get("FAKE_INSTALL_FAILURE"):
        print("R install fixture error", file=sys.stderr)
        sys.exit(1)
    Path(os.environ["CONDA_PREFIX"], "phylopath-installed").write_text("1.3.1")
elif tool == "Rscript" and stdin:
    failure = os.environ.get("FAKE_R_FAILURE")
    if not Path(os.environ["CONDA_PREFIX"], "phylopath-installed").exists():
        failure = "phylopath 1.3.1 is required"
    if failure:
        print(failure, file=sys.stderr)
        sys.exit(1)
elif tool in {"gs", "convert"} and os.environ.get("FAKE_CONVERTER_FAILURE"):
    sys.exit(1)
"""


class RFixture:
    def __init__(self, directory: Path) -> None:
        self.directory = directory
        self.project = directory / "project with spaces"
        (self.project / "scripts/publication").mkdir(parents=True)
        for name in ("run_r.sh", "setup_r.sh", "r_environment.sh"):
            shutil.copy2(ROOT / "scripts" / name, self.project / "scripts" / name)
        self.registry = directory / "environments.txt"
        self.registry.write_text("")
        # Keep the host HOME untouched; redirect only the copied registry fixture.
        helper = self.project / "scripts/r_environment.sh"
        helper.write_text(
            helper.read_text().replace(
                '"$HOME/.conda/environments.txt"',
                f'"{self.registry}"',
            )
        )
        shutil.copy2(ROOT / "environment.yml", self.project / "environment.yml")
        shutil.copy2(
            ROOT / "environment-linux-64.explicit.txt",
            self.project / "environment-linux-64.explicit.txt",
        )
        shutil.copy2(
            ROOT / "scripts/publication/gbe_theme.R",
            self.project / "scripts/publication/gbe_theme.R",
        )
        self.tools = directory / "tools"
        self.tools.mkdir()
        self.template = directory / "template"
        (self.template / "bin").mkdir(parents=True)
        (self.template / "lib/R/library").mkdir(parents=True)
        (self.template / "phylopath-installed").write_text("1.3.1")
        for name in ("conda", "curl", "sha256sum", "uname"):
            self.write_tool(self.tools / name)
        for name in ("Rscript", "R", "gs", "convert"):
            self.write_tool(self.template / "bin" / name)
        self.registered = directory / "registered.json"
        self.registered.write_text("[]")
        self.log = directory / "calls.jsonl"
        self.outside = directory / "outside"
        self.outside.mkdir()
        self.environment = {
            **os.environ,
            "PATH": f"{self.tools}:/usr/bin:/bin",
            "DESMOGNATHUS_R_ENV": "desmognathus",
            "FAKE_LOG": str(self.log),
            "FAKE_REGISTERED": str(self.registered),
            "FAKE_REGISTRY": str(self.registry),
            "FAKE_TEMPLATE": str(self.template),
            "FAKE_ENV_DIR": str(directory / "custom env directory"),
            "FAKE_SHA256": PHYLOPATH_SHA256,
            "FAKE_OS": "Linux",
            "FAKE_ARCH": "x86_64",
            "R_LIBS": "/unrelated/user/library",
        }

    @staticmethod
    def write_tool(path: Path) -> None:
        path.write_text(f"#!{sys.executable}\n{FAKE_TOOL}")
        path.chmod(0o755)

    def existing(self, name: str = "desmognathus") -> Path:
        prefix = self.directory / "registered outside default envs" / name
        shutil.copytree(self.template, prefix)
        self.registered.write_text(json.dumps([str(prefix)]))
        return prefix

    def run(self, script: str, *args: str) -> subprocess.CompletedProcess[str]:
        return subprocess.run(
            ["bash", str(self.project / "scripts" / script), *args],
            cwd=self.outside,
            env=self.environment,
            text=True,
            capture_output=True,
            check=False,
        )

    def calls(self) -> list[dict]:
        return [json.loads(line) for line in self.log.read_text().splitlines()]


@pytest.fixture
def r_fixture(tmp_path: Path) -> RFixture:
    return RFixture(tmp_path)


def snapshot(directory: Path) -> dict[str, str]:
    return {
        str(path.relative_to(directory)): hashlib.sha256(path.read_bytes()).hexdigest()
        for path in directory.rglob("*")
        if path.is_file()
    }


@pytest.mark.parametrize("arguments", [(), ("--check",)])
def test_existing_setup_is_read_only(r_fixture: RFixture, arguments: tuple[str, ...]) -> None:
    prefix = r_fixture.existing()
    before = snapshot(prefix)
    result = r_fixture.run("setup_r.sh", *arguments)
    assert result.returncode == 0, result.stderr
    assert snapshot(prefix) == before
    assert [(call["tool"], call["args"]) for call in r_fixture.calls()] == [
        ("conda", ["config", "--show", "envs_dirs"]),
        (
            "Rscript",
            [
                "--vanilla",
                "-",
                str(r_fixture.project / "environment.yml"),
                str(r_fixture.project / "scripts/publication/gbe_theme.R"),
            ],
        ),
        ("gs", ["--version"]),
        ("convert", ["-version"]),
    ]


@pytest.mark.parametrize(
    "failure",
    [
        "Missing required R package: ape",
        "R 4.4.x is required; found 4.3.3",
        "phylopath 1.3.1 is required",
        "Required rmarkdown version: 2.31; found 2.30",
        "Arial is required for GBE figures",
        "R Cairo support is required",
    ],
)
def test_existing_dependency_failure_never_repairs(r_fixture: RFixture, failure: str) -> None:
    prefix = r_fixture.existing()
    before = snapshot(prefix)
    r_fixture.environment["FAKE_R_FAILURE"] = failure
    result = r_fixture.run("setup_r.sh")
    assert result.returncode == 1
    assert failure in result.stderr
    assert "Existing environment was not changed" in result.stderr
    assert "<new-name-or-prefix>" in result.stderr
    assert snapshot(prefix) == before
    assert [call["tool"] for call in r_fixture.calls()] == ["conda", "Rscript"]


def test_missing_check_does_not_provision(r_fixture: RFixture) -> None:
    result = r_fixture.run("setup_r.sh", "--check")
    assert result.returncode == 1
    assert "--check does not create" in result.stderr
    assert [call["args"] for call in r_fixture.calls()] == [
        ["config", "--show", "envs_dirs"],
    ]


def test_unset_selection_defaults_to_desmognathus(r_fixture: RFixture) -> None:
    prefix = r_fixture.existing()
    r_fixture.environment.pop("DESMOGNATHUS_R_ENV")
    result = r_fixture.run("run_r.sh", "Rscript", "script.R")
    assert result.returncode == 0, result.stderr
    assert r_fixture.calls()[-1]["prefix"] == str(prefix)


@pytest.mark.parametrize("stage", ["CREATE", "INSTALL"])
def test_provisioning_failure_is_reported_and_retained(r_fixture: RFixture, stage: str) -> None:
    r_fixture.environment[f"FAKE_{stage}_FAILURE"] = "1"
    result = r_fixture.run("setup_r.sh")
    assert result.returncode == 1
    assert f"{stage.lower()} fixture error" in result.stderr
    assert not any(call["tool"] == "Rscript" for call in r_fixture.calls())
    if stage == "INSTALL":
        prefix = Path(json.loads(r_fixture.registered.read_text())[0])
        assert prefix.is_dir()
        assert not (prefix / "phylopath-installed").exists()


@pytest.mark.parametrize("use_prefix", [False, True])
@pytest.mark.parametrize(
    "system,architecture", [("Linux", "x86_64"), ("Linux", "aarch64"), ("Darwin", "x86_64")]
)
def test_missing_setup_creates_once_then_only_validates(
    r_fixture: RFixture,
    use_prefix: bool,
    system: str,
    architecture: str,
) -> None:
    selection = "separate-env"
    if use_prefix:
        selection = str(r_fixture.directory / "new custom prefix")
    r_fixture.environment["DESMOGNATHUS_R_ENV"] = selection
    r_fixture.environment["FAKE_OS"] = system
    r_fixture.environment["FAKE_ARCH"] = architecture
    spec_before = (r_fixture.project / "environment.yml").read_bytes()
    lock_before = (r_fixture.project / "environment-linux-64.explicit.txt").read_bytes()
    first = r_fixture.run("setup_r.sh")
    assert first.returncode == 0, first.stderr
    prefix = Path(json.loads(r_fixture.registered.read_text())[0])
    assert (prefix / "phylopath-installed").read_text() == "1.3.1"
    before = snapshot(prefix)
    second = r_fixture.run("setup_r.sh")
    assert second.returncode == 0, second.stderr
    calls = r_fixture.calls()
    creates = [call for call in calls if call["tool"] == "conda" and "create" in call["args"]]
    assert len(creates) == 1
    option = "--prefix" if use_prefix else "--name"
    locked = (system, architecture) == ("Linux", "x86_64")
    command = ["create", "--yes"] if locked else ["env", "create"]
    specification = "environment-linux-64.explicit.txt" if locked else "environment.yml"
    assert creates[0]["args"] == [
        *command,
        option,
        selection,
        "--file",
        str(r_fixture.project / specification),
    ]
    assert len([call for call in calls if call["tool"] == "R"]) == 1
    assert snapshot(prefix) == before
    assert (r_fixture.project / "environment.yml").read_bytes() == spec_before
    assert (r_fixture.project / "environment-linux-64.explicit.txt").read_bytes() == lock_before


def test_missing_linux_lock_does_not_silently_resolve_packages(r_fixture: RFixture) -> None:
    (r_fixture.project / "environment-linux-64.explicit.txt").unlink()
    result = r_fixture.run("setup_r.sh")
    assert result.returncode == 1
    assert "Missing Linux x86_64 package lock" in result.stderr
    assert not any(
        call["tool"] == "conda" and "create" in call["args"] for call in r_fixture.calls()
    )


@pytest.mark.parametrize("arguments", [(), ("--check",)])
def test_existing_environment_does_not_need_creation_lock(
    r_fixture: RFixture,
    arguments: tuple[str, ...],
) -> None:
    r_fixture.existing()
    (r_fixture.project / "environment-linux-64.explicit.txt").unlink()
    result = r_fixture.run("setup_r.sh", *arguments)
    assert result.returncode == 0, result.stderr
    assert not any(call["tool"] == "uname" for call in r_fixture.calls())


@pytest.mark.parametrize("failure", ["checksum", "download"])
def test_new_environment_source_failure_never_installs(r_fixture: RFixture, failure: str) -> None:
    if failure == "checksum":
        r_fixture.environment["FAKE_SHA256"] = "wrong"
    else:
        r_fixture.environment["FAKE_DOWNLOAD_FAILURE"] = "1"
    result = r_fixture.run("setup_r.sh")
    assert result.returncode == 1
    expected = "checksum mismatch" if failure == "checksum" else "Could not download"
    assert expected in result.stderr
    assert not any(call["tool"] == "R" for call in r_fixture.calls())
    # A failed new environment is retained and will never be repaired implicitly.
    retry = r_fixture.run("setup_r.sh", "--check")
    assert retry.returncode == 1
    assert "Existing environment was not changed" in retry.stderr
    creates = [
        call for call in r_fixture.calls() if call["tool"] == "conda" and "create" in call["args"]
    ]
    assert len(creates) == 1


@pytest.mark.parametrize("use_prefix", [False, True])
def test_runner_preserves_arguments_root_and_isolation(
    r_fixture: RFixture, use_prefix: bool
) -> None:
    prefix = r_fixture.existing("custom-name")
    r_fixture.environment["DESMOGNATHUS_R_ENV"] = str(prefix) if use_prefix else "custom-name"
    arguments = ["--no-environ", "script with spaces.R", "--output", "value 'quoted' $literal"]
    result = r_fixture.run("run_r.sh", "Rscript", *arguments)
    assert result.returncode == 0, result.stderr
    call = r_fixture.calls()[-1]
    assert call["args"] == arguments
    assert call["cwd"] == str(r_fixture.project)
    assert call["prefix"] == str(prefix)
    assert call["r_environ"] == call["r_profile"] == "/dev/null"
    assert call["r_libs"] is None
    assert call["r_libs_user"] == str(prefix / "lib/R/library")
    if use_prefix:
        assert not any(item["tool"] == "conda" for item in r_fixture.calls())


def test_relative_prefix_has_same_meaning_for_setup_and_runner(r_fixture: RFixture) -> None:
    prefix = r_fixture.existing()
    r_fixture.environment["DESMOGNATHUS_R_ENV"] = os.path.relpath(prefix, r_fixture.outside)
    check = r_fixture.run("setup_r.sh", "--check")
    run = r_fixture.run("run_r.sh", "Rscript", "script.R")
    assert check.returncode == run.returncode == 0, check.stderr + run.stderr
    assert all(call["prefix"] == str(prefix) for call in r_fixture.calls())


def test_registered_prefix_outside_env_directories_is_resolved(r_fixture: RFixture) -> None:
    prefix = r_fixture.existing()
    r_fixture.registered.write_text("[]")
    r_fixture.registry.write_text(f"{prefix}\n/stale/desmognathus\n")
    before = r_fixture.registry.read_bytes()
    result = r_fixture.run("run_r.sh", "Rscript", "script.R")
    assert result.returncode == 0, result.stderr
    assert r_fixture.calls()[-1]["prefix"] == str(prefix)
    assert r_fixture.registry.read_bytes() == before


def test_explicit_prefix_check_does_not_require_conda(r_fixture: RFixture) -> None:
    prefix = r_fixture.existing()
    (r_fixture.tools / "conda").unlink()
    r_fixture.environment["DESMOGNATHUS_R_ENV"] = str(prefix)
    result = r_fixture.run("setup_r.sh", "--check")
    assert result.returncode == 0, result.stderr
    assert not any(call["tool"] == "conda" for call in r_fixture.calls())


@pytest.mark.parametrize("script", ["setup_r.sh", "run_r.sh"])
def test_missing_conda_for_named_environment_fails_without_provisioning(
    r_fixture: RFixture,
    script: str,
) -> None:
    (r_fixture.tools / "conda").unlink()
    arguments = ["--check"] if script == "setup_r.sh" else ["Rscript", "script.R"]
    result = r_fixture.run(script, *arguments)
    assert result.returncode == 1
    assert "Conda is required to resolve" in result.stderr
    assert not r_fixture.log.exists()


def test_missing_runner_never_provisions(r_fixture: RFixture) -> None:
    result = r_fixture.run("run_r.sh", "Rscript", "script.R")
    assert result.returncode == 1
    assert "Run: make setup-r" in result.stderr
    assert [call["args"] for call in r_fixture.calls()] == [
        ["config", "--show", "envs_dirs"],
    ]


def test_failed_converter_fails_read_only(r_fixture: RFixture) -> None:
    prefix = r_fixture.existing()
    before = snapshot(prefix)
    r_fixture.environment["FAKE_CONVERTER_FAILURE"] = "1"
    result = r_fixture.run("setup_r.sh", "--check")
    assert result.returncode == 1
    assert "Existing environment was not changed" in result.stderr
    assert snapshot(prefix) == before


@pytest.mark.parametrize("script", ["setup_r.sh", "run_r.sh"])
def test_conda_error_is_not_treated_as_missing(r_fixture: RFixture, script: str) -> None:
    r_fixture.environment["FAKE_CONDA_FAILURE"] = "1"
    result = r_fixture.run(script, *(["Rscript", "script.R"] if script == "run_r.sh" else []))
    assert result.returncode == 1
    assert "Could not read Conda environment directories" in result.stderr
    assert [call["args"] for call in r_fixture.calls()] == [
        ["config", "--show", "envs_dirs"],
    ]


@pytest.mark.parametrize("missing", ["R", "Rscript", "gs", "convert"])
def test_missing_runtime_or_converter_fails_read_only(r_fixture: RFixture, missing: str) -> None:
    prefix = r_fixture.existing()
    (prefix / "bin" / missing).unlink()
    before = snapshot(prefix)
    result = r_fixture.run("setup_r.sh", "--check")
    assert result.returncode == 1
    assert "Existing environment was not changed" in result.stderr
    assert snapshot(prefix) == before
    assert not any(call["tool"] in {"R", "curl"} for call in r_fixture.calls())


def test_ambiguous_name_requires_explicit_prefix(r_fixture: RFixture) -> None:
    first = r_fixture.existing()
    second = r_fixture.directory / "another location/desmognathus"
    shutil.copytree(r_fixture.template, second)
    r_fixture.registered.write_text(json.dumps([str(first), str(second)]))
    result = r_fixture.run("setup_r.sh")
    assert result.returncode == 1
    assert "Ambiguous R environment" in result.stderr
    assert [call["tool"] for call in r_fixture.calls()] == ["conda"]


@pytest.mark.parametrize(
    "script,arguments",
    [
        ("setup_r.sh", ["--unknown"]),
        ("setup_r.sh", ["--check", "extra"]),
        ("run_r.sh", ["R", "script.R"]),
    ],
)
def test_invalid_usage_never_invokes_conda(
    r_fixture: RFixture,
    script: str,
    arguments: list[str],
) -> None:
    result = r_fixture.run(script, *arguments)
    assert result.returncode == 2
    assert "Usage:" in result.stderr
    assert not r_fixture.log.exists()
