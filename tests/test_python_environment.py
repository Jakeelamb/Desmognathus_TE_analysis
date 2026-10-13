from __future__ import annotations

import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]


@pytest.mark.parametrize("correct_version", [True, False])
def test_interpreter_check_accepts_pin_and_rejects_other_version(
    tmp_path: Path, correct_version: bool
) -> None:
    scripts = tmp_path / "scripts"
    scripts.mkdir()
    checker = scripts / "check_python_environment.py"
    shutil.copyfile(ROOT / "scripts/check_python_environment.py", checker)
    actual_version = ".".join(map(str, sys.version_info[:3]))
    pin = actual_version if correct_version else "0.0.0"
    (tmp_path / ".python-version").write_text(pin + "\n")
    result = subprocess.run(
        [sys.executable, "-B", str(checker)], capture_output=True, text=True, check=False
    )
    assert result.returncode == (0 if correct_version else 1)
    if correct_version:
        assert f"PASS; Python {actual_version}" in result.stdout
    else:
        assert f"expected {pin}, found {actual_version}" in result.stderr
        assert "new UV_PROJECT_ENVIRONMENT" in result.stderr
    assert (tmp_path / ".python-version").read_text() == pin + "\n"
    assert not (scripts / "__pycache__").exists()


@pytest.fixture
def make_fixture(tmp_path: Path) -> tuple[Path, dict[str, str], Path]:
    shutil.copyfile(ROOT / "Makefile", tmp_path / "Makefile")
    binaries = tmp_path / "bin"
    binaries.mkdir()
    scripts = tmp_path / "scripts"
    scripts.mkdir()
    uv = binaries / "uv"
    uv.write_text(
        '#!/bin/sh\nprintf "uv %s\\n" "$*" >> "$COMMAND_LOG"\n'
        'if [ "$1" = sync ] && [ "${FAIL_SYNC:-0}" = 1 ]; then exit 1; fi\n'
        'exit "${FAIL_UV:-0}"\n'
    )
    uv.chmod(0o755)
    setup_r = scripts / "setup_r.sh"
    setup_r.write_text('#!/bin/sh\nprintf "r %s\\n" "$*" >> "$COMMAND_LOG"\n')
    setup_r.chmod(0o755)
    environment = tmp_path / "env"
    (environment / "bin").mkdir(parents=True)
    (environment / "bin/python").symlink_to(sys.executable)
    log = tmp_path / "commands.log"
    variables = {
        **os.environ,
        "PATH": str(binaries) + os.pathsep + os.environ["PATH"],
        "UV_PROJECT_ENVIRONMENT": str(environment),
        "COMMAND_LOG": str(log),
    }
    return tmp_path, variables, log


def test_check_env_checks_python_lock_before_r(
    make_fixture: tuple[Path, dict[str, str], Path],
) -> None:
    directory, variables, log = make_fixture
    result = subprocess.run(
        ["make", "check-env"], cwd=directory, env=variables, capture_output=True, check=False
    )
    assert result.returncode == 0, result.stderr
    assert log.read_text().splitlines() == [
        (
            "uv run --locked --no-sync --offline --no-python-downloads "
            "python -B scripts/check_python_environment.py"
        ),
        "uv sync --locked --check --all-groups --offline --no-python-downloads",
        "r --check",
    ]


def test_failed_python_check_blocks_lock_and_r_checks(
    make_fixture: tuple[Path, dict[str, str], Path],
) -> None:
    directory, variables, log = make_fixture
    variables["FAIL_UV"] = "1"
    result = subprocess.run(
        ["make", "check-env"], cwd=directory, env=variables, capture_output=True, check=False
    )
    assert result.returncode != 0
    commands = log.read_text().splitlines()
    assert len(commands) == 1
    assert commands[0].startswith("uv run --locked --no-sync")


def test_missing_environment_does_not_call_uv_or_r(
    make_fixture: tuple[Path, dict[str, str], Path],
) -> None:
    directory, variables, log = make_fixture
    missing = directory / "missing-env"
    variables["UV_PROJECT_ENVIRONMENT"] = str(missing)
    result = subprocess.run(
        ["make", "check-env"],
        cwd=directory,
        env=variables,
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode != 0
    assert "Python environment missing" in result.stderr
    assert not log.exists()
    assert not missing.exists()


def test_unsynchronized_packages_block_r_check(
    make_fixture: tuple[Path, dict[str, str], Path],
) -> None:
    directory, variables, log = make_fixture
    variables["FAIL_SYNC"] = "1"
    result = subprocess.run(
        ["make", "check-env"], cwd=directory, env=variables, capture_output=True, check=False
    )
    assert result.returncode != 0
    commands = log.read_text().splitlines()
    assert len(commands) == 2
    assert commands[1].startswith("uv sync --locked --check --all-groups --offline")


def test_python_make_entrypoints_are_locked_and_do_not_sync() -> None:
    targets = [
        "status", "catalog", "publication", "refresh", "identity", "validate", "test", "lint",
        "report", "figure-review", "figure-data", "figures", "methods-word", "te-bonsai",
        "te-pca", "path", "path-quick",
    ]
    result = subprocess.run(
        ["make", "-n", *targets], cwd=ROOT, capture_output=True, text=True, check=False
    )
    assert result.returncode == 0, result.stderr
    commands = result.stdout.splitlines()
    python_commands = [command for command in commands if command.startswith("uv ")]
    assert python_commands
    assert all(command.startswith("uv run --locked --no-sync ") for command in python_commands)
    assert "uv run --locked --no-sync python -m pytest -q" in commands
    assert "uv run --locked --no-sync python -m scripts.project publication" in commands
    assert "uv run --locked --no-sync python -m scripts.project catalog" in commands
    assert (
        "uv run --locked --no-sync python scripts/publication/build_methods_document.py --check"
        in commands
    )


def test_setup_python_uses_existing_lock(make_fixture: tuple[Path, dict[str, str], Path]) -> None:
    directory, variables, log = make_fixture
    result = subprocess.run(
        ["make", "setup-python"], cwd=directory, env=variables, capture_output=True, check=False
    )
    assert result.returncode == 0, result.stderr
    assert log.read_text().splitlines() == ["uv sync --locked --all-groups"]
