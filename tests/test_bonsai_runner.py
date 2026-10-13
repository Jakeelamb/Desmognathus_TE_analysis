from __future__ import annotations

import importlib.util
import json
import os
import re
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("bonsai_runner", ROOT / "scripts/run_bonsai.py")
assert SPEC is not None and SPEC.loader is not None
RUNNER = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(RUNNER)
PINNED_COMMIT = RUNNER.UPSTREAM_COMMIT

UPSTREAM_PROGRAM = """import argparse, json, os
from pathlib import Path
parser = argparse.ArgumentParser()
parser.add_argument('--config_filepath', required=True)
parser.add_argument('--step', choices=['all'], required=True)
args = parser.parse_args()
config = json.loads(Path(args.config_filepath).read_text())
directory = Path(args.config_filepath).parent
record = {'step': args.step, 'python': os.sys.executable,
          'threads': {key: os.environ.get(key) for key in
          ['OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS']}}
(directory / 'upstream_call.json').write_text(json.dumps(record))
assert os.environ['PYTHONDONTWRITEBYTECODE'] == '1'
if os.environ.get('FAKE_BONSAI_MODE') == 'fail':
    print('fixture upstream failure', flush=True)
    raise SystemExit(9)
output = Path(config['results_folder']) / 'final_bonsai'
output.mkdir(parents=True, exist_ok=True)
(output / 'tree.nwk').write_text('' if os.environ.get('FAKE_BONSAI_MODE') == 'empty'
                               else '(a:1,b:1);\\n')
(output / 'metadata_bonsai.txt').write_text('reconsTime\\n1.5\\n')
print('fixture upstream success', flush=True)
"""

HELPER_PROGRAM = """import argparse, json
from pathlib import Path
parser = argparse.ArgumentParser()
parser.add_argument('mode')
parser.add_argument('--run-dir', type=Path)
parser.add_argument('--output-dir', type=Path)
parser.add_argument('--tree', type=Path)
parser.add_argument('--metadata', type=Path)
parser.add_argument('--replicates', type=int)
parser.add_argument('--seed', type=int)
parser.add_argument('--write', action='store_true')
args = parser.parse_args()
if args.mode.startswith('prepare'):
    if args.mode == 'prepare-bootstrap':
        plan = args.run_dir / 'te34_bonsai_feature_bootstrap_plan.json'
        identity = {'replicates': args.replicates, 'seed': args.seed}
        if plan.exists():
            assert json.loads(plan.read_text()) == identity
            raise SystemExit(0)
        plan.write_text(json.dumps(identity))
        (args.run_dir / 'te34_bonsai_feature_bootstrap_draws.csv').write_text('replicate_id\\n0\\n')
        directories = [args.run_dir / 'bootstrap' / f'replicate_{i:03d}'
                       for i in range(args.replicates + 1)]
    else:
        directories = [args.run_dir]
    for directory in directories:
        (directory / 'input').mkdir(parents=True)
        (directory / 'results').mkdir()
        for name in ['features.txt', 'cellID.txt', 'geneID.txt']:
            (directory / 'input' / name).write_text('fixture input\\n')
        config = {'results_folder': str(directory / 'results')}
        (directory / 'bonsai_config.yaml').write_text(json.dumps(config))
else:
    assert args.write
    if args.tree:
        assert args.tree.read_text().strip() == '(a:1,b:1);'
        assert args.metadata.is_file()
    args.output_dir.mkdir()
    (args.output_dir / 'products.json').write_text(json.dumps({'mode': args.mode}))
"""


def git(path: Path, *args: str) -> str:
    return subprocess.check_output(["git", "-C", str(path), *args], text=True).strip()


@pytest.fixture
def campaign(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> tuple[Path, Path]:
    project = tmp_path / "project with spaces"
    project.mkdir()
    comparison = project / "analyses/01_transposable_elements/bonsai_comparison"
    comparison.mkdir(parents=True)
    for name in ["compare_bonsai_to_pca.py", "validate_strategy_localization.py"]:
        (comparison / name).write_text(HELPER_PROGRAM)
    scripts = project / "scripts"
    scripts.mkdir()
    plotter = scripts / "run_r.sh"
    plotter.write_text(
        f"#!{sys.executable}\nimport sys\nfrom pathlib import Path\n"
        "(Path(sys.argv[-1]) / 'plot.pdf').write_text('fixture plot')\n"
    )
    plotter.chmod(0o755)
    tools = tmp_path / "fixture tools"
    tools.mkdir()
    pixi = tools / "pixi"
    pixi.write_text(
        f"#!{sys.executable}\nimport json, os, subprocess, sys\n"
        "from pathlib import Path\n"
        "assert '--locked' in sys.argv and '--no-install' in sys.argv\n"
        "assert sys.argv[sys.argv.index('--environment') + 1] == 'default'\n"
        "with Path(os.environ['FIXTURE_PIXI_CALLS']).open('a') as log:\n"
        "    log.write(json.dumps(sys.argv[1:]) + '\\n')\n"
        "index = sys.argv.index('--executable')\n"
        "raise SystemExit(subprocess.call([sys.executable, *sys.argv[index + 2:]]))\n"
    )
    pixi.chmod(0o755)
    monkeypatch.setenv("PATH", f"{tools}:{os.environ['PATH']}")
    monkeypatch.setenv("FIXTURE_PIXI_CALLS", str(tmp_path / "pixi calls.jsonl"))
    source = tmp_path / "upstream with spaces"
    (source / "bonsai").mkdir(parents=True)
    (source / "requirements.txt").write_text("# mpi4py\n" + "\n".join(sorted(RUNNER.REQUIREMENTS)))
    (source / "bonsai/bonsai_main.py").write_text(UPSTREAM_PROGRAM)
    git(source, "init", "-q")
    git(source, "add", ".")
    git(
        source,
        "-c",
        "user.name=Runner Fixture",
        "-c",
        "user.email=fixture@example.invalid",
        "commit",
        "-qm",
        "fixture Bonsai implementation",
    )
    monkeypatch.setattr(RUNNER, "ROOT", project)
    monkeypatch.setattr(RUNNER, "COMPARISON", comparison)
    monkeypatch.setattr(RUNNER, "UPSTREAM_COMMIT", git(source, "rev-parse", "HEAD"))
    # Execution fixtures exercise orchestration using the hydrated audit Python;
    # the production guard is checked separately below.
    monkeypatch.setattr(RUNNER, "require_pixi_environment", lambda: None)
    return project, source


def invoke(run_dir: Path, source: Path, mode: str = "bootstrap", *args: str) -> int:
    return RUNNER.main([mode, "--run-dir", str(run_dir), "--upstream", str(source), *args])


def test_refit_executes_prepare_upstream_analyze_and_plot(campaign: tuple[Path, Path]) -> None:
    project, source = campaign
    run_dir = project / ".pixi/runs/refit with spaces"
    assert invoke(run_dir, source, "refit") == 0
    assert json.loads((run_dir / "comparison/products.json").read_text()) == {"mode": "analyze"}
    assert (run_dir / "comparison/plot.pdf").is_file()
    call = json.loads((run_dir / "upstream_call.json").read_text())
    assert call["step"] == "all" and call["python"] == sys.executable
    assert "fixture upstream success" in (run_dir / "bonsai.log").read_text()
    assert git(source, "status", "--porcelain", "--untracked-files=all") == ""
    calls = Path(os.environ["FIXTURE_PIXI_CALLS"]).read_text().splitlines()
    assert len(calls) == 2


def test_bounded_bootstrap_calibrates_and_resumes_nonempty_trees(
    campaign: tuple[Path, Path], capsys: pytest.CaptureFixture[str]
) -> None:
    project, source = campaign
    run_dir = project / ".pixi/runs/smoke"
    assert invoke(run_dir, source, "bootstrap", "--replicates", "2", "--jobs", "1") == 0
    assert "execution smoke only" in capsys.readouterr().out
    assert not (run_dir / "validation").exists()
    replicates = sorted((run_dir / "bootstrap").iterdir())
    assert [path.name for path in replicates] == ["replicate_000", "replicate_001", "replicate_002"]
    plan = json.loads((run_dir / RUNNER.BOOTSTRAP_PLAN).read_text())
    assert plan == {"replicates": 2, "seed": 20260826}
    for path in replicates:
        call = json.loads((path / "upstream_call.json").read_text())
        assert set(call["threads"].values()) == {"1"}
    retained_log = (replicates[0] / "bonsai.log").read_bytes()
    (replicates[1] / "results/final_bonsai/tree.nwk").write_text("")
    (replicates[1] / "bonsai.log").write_text("old failed log")
    assert invoke(run_dir, source, "bootstrap", "--replicates", "2", "--resume") == 0
    output = capsys.readouterr().out
    assert "SKIP replicate_000" in output and "SKIP replicate_002" in output
    assert "DONE replicate_001" in output
    assert (replicates[0] / "bonsai.log").read_bytes() == retained_log
    assert git(source, "status", "--porcelain", "--untracked-files=all") == ""


def test_default_campaign_analyzes_all_200_replicates_plus_calibration(
    campaign: tuple[Path, Path],
) -> None:
    project, source = campaign
    run_dir = project / ".pixi/runs/campaign"
    assert invoke(run_dir, source) == 0
    assert (
        len(list((run_dir / "bootstrap").glob("replicate_*/results/final_bonsai/tree.nwk"))) == 201
    )
    assert json.loads((run_dir / "validation/products.json").read_text()) == {
        "mode": "analyze-bootstrap"
    }
    assert (run_dir / "validation/plot.pdf").is_file()
    assert not (run_dir / ".bonsai-run.lock").exists()


def test_run_lock_refuses_concurrent_resume_without_replacing_lock(
    campaign: tuple[Path, Path],
) -> None:
    project, source = campaign
    run_dir = project / ".pixi/runs/owned"
    assert invoke(run_dir, source, "bootstrap", "--replicates", "1") == 0
    lock = run_dir / ".bonsai-run.lock"
    lock.write_text("other runner owns this")
    assert invoke(run_dir, source, "bootstrap", "--replicates", "1", "--resume") == 1
    assert lock.read_text() == "other runner owns this"


@pytest.mark.parametrize(
    "path",
    ["data/new", "Publication/new", "analyses/new", "outputs/new", "scripts/new", ".git/new"],
)
def test_frozen_project_paths_are_rejected_without_writes(
    campaign: tuple[Path, Path], path: str
) -> None:
    project, source = campaign
    run_dir = project / path
    assert invoke(run_dir, source, "refit") == 1
    assert not run_dir.exists()
    assert not Path(os.environ["FIXTURE_PIXI_CALLS"]).exists()


def test_existing_unrelated_directories_and_refit_resume_are_rejected(
    campaign: tuple[Path, Path],
) -> None:
    project, source = campaign
    run_dir = project / ".pixi/runs/unrelated"
    run_dir.mkdir(parents=True)
    sentinel = run_dir / "important.txt"
    sentinel.write_text("keep")
    assert invoke(run_dir, source) == 1
    assert invoke(run_dir, source, "bootstrap", "--resume") == 1
    assert invoke(run_dir, source, "refit", "--resume") == 1
    assert list(run_dir.iterdir()) == [sentinel] and sentinel.read_text() == "keep"


@pytest.mark.parametrize(
    "change", ["config", "input", "manifest", "replicates", "symlink", "hardlink"]
)
def test_resume_rejects_changed_preparation_and_redirected_writes(
    campaign: tuple[Path, Path], tmp_path: Path, change: str
) -> None:
    project, source = campaign
    run_dir = project / ".pixi/runs/resume"
    assert invoke(run_dir, source, "bootstrap", "--replicates", "1") == 0
    directory = run_dir / "bootstrap/replicate_000"
    retained = tmp_path / "protected log"
    retained.write_text("must survive")
    extra = []
    if change == "config":
        (directory / "bonsai_config.yaml").write_text("{}")
    elif change == "input":
        (directory / "input/features.txt").write_text("changed")
    elif change == "manifest":
        (run_dir / RUNNER.RUN_MANIFEST).write_text("{}")
    elif change == "replicates":
        extra = ["--replicates", "2"]
    else:
        log = directory / "bonsai.log"
        log.unlink()
        if change == "symlink":
            log.symlink_to(retained)
        else:
            os.link(retained, log)
    assert invoke(run_dir, source, "bootstrap", "--replicates", "1", "--resume", *extra) == 1
    assert retained.read_text() == "must survive"


@pytest.mark.parametrize("failure", ["fail", "empty"])
def test_failed_or_empty_upstream_result_is_not_success(
    campaign: tuple[Path, Path], monkeypatch: pytest.MonkeyPatch, failure: str
) -> None:
    project, source = campaign
    run_dir = project / ".pixi/runs/failure"
    monkeypatch.setenv("FAKE_BONSAI_MODE", failure)
    assert invoke(run_dir, source, "refit") == 1
    assert not (run_dir / "comparison").exists()
    assert (run_dir / "bonsai.log").is_file()


def test_dirty_wrong_commit_and_corrupt_source_caches_are_never_repaired(
    campaign: tuple[Path, Path], monkeypatch: pytest.MonkeyPatch
) -> None:
    project, source = campaign
    (source / "requirements.txt").write_text("numpy==0.0.0")
    assert invoke(project / ".pixi/runs/dirty", source, "refit") == 1
    assert (source / "requirements.txt").read_text() == "numpy==0.0.0"
    git(source, "add", ".")
    git(
        source,
        "-c",
        "user.name=Fixture",
        "-c",
        "user.email=fixture@example.invalid",
        "commit",
        "-qm",
        "different commit",
    )
    assert invoke(project / ".pixi/runs/wrong", source, "refit") == 1
    cache = project / ".pixi/sources" / f"bonsai-{RUNNER.UPSTREAM_COMMIT}"
    cache.mkdir(parents=True)
    (cache / "sentinel").write_text("corrupt cache stays")
    monkeypatch.setattr(RUNNER, "UPSTREAM_URL", str(source))
    assert RUNNER.main(["refit", "--run-dir", str(project / ".pixi/runs/cache")]) == 1
    assert list(cache.iterdir()) == [cache / "sentinel"]
    assert (cache / "sentinel").read_text() == "corrupt cache stays"


def test_first_source_cache_fetches_exact_commit_from_explicit_origin(
    campaign: tuple[Path, Path], monkeypatch: pytest.MonkeyPatch
) -> None:
    project, source = campaign
    monkeypatch.setattr(RUNNER, "UPSTREAM_URL", str(source))
    run_dir = project / ".pixi/runs/cached"
    assert RUNNER.main(["refit", "--run-dir", str(run_dir)]) == 0
    cache = project / ".pixi/sources" / f"bonsai-{RUNNER.UPSTREAM_COMMIT}"
    assert git(cache, "rev-parse", "HEAD") == git(source, "rev-parse", "HEAD")
    assert git(cache, "status", "--porcelain", "--untracked-files=all") == ""


def test_redirected_cache_and_upstream_overlap_are_rejected_without_writes(
    campaign: tuple[Path, Path], tmp_path: Path
) -> None:
    project, source = campaign
    cache_parent = project / ".pixi/sources"
    cache_parent.parent.mkdir()
    protected = tmp_path / "preserved source location"
    protected.mkdir()
    cache_parent.symlink_to(protected, target_is_directory=True)
    assert RUNNER.main(["refit", "--run-dir", str(project / ".pixi/runs/cache")]) == 1
    assert list(protected.iterdir()) == []
    assert invoke(source / "scratch", source, "refit") == 1
    assert not (source / "scratch").exists()


def test_git_worktree_source_pointer_is_verified_read_only(
    campaign: tuple[Path, Path], tmp_path: Path
) -> None:
    _, source = campaign
    linked = tmp_path / "linked upstream"
    git(source, "worktree", "add", "--quiet", "--detach", str(linked), RUNNER.UPSTREAM_COMMIT)
    pointer = (linked / ".git").read_bytes()
    assert pointer.startswith(b"gitdir: ")
    assert RUNNER.verify_upstream(linked) == linked
    assert (linked / ".git").read_bytes() == pointer


def test_pinned_production_commit_and_environment_guard(monkeypatch: pytest.MonkeyPatch) -> None:
    assert re.fullmatch(r"[0-9a-f]{40}", PINNED_COMMIT)
    assert PINNED_COMMIT == "b16ef6f0a164610327022c7d0b92de82b3596aad"
    monkeypatch.delenv("PIXI_ENVIRONMENT_NAME", raising=False)
    with pytest.raises(ValueError, match="locked Bonsai Pixi task"):
        RUNNER.require_pixi_environment()


@pytest.fixture
def bonsai_environment(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Path:
    project = tmp_path / "project"
    prefix = project / ".pixi/envs/bonsai"
    (prefix / "bin").mkdir(parents=True)
    for name in ["python", "git"]:
        (prefix / "bin" / name).touch()
    default_python = project / ".pixi/envs/default/bin/python"
    default_python.parent.mkdir(parents=True)
    default_python.touch()
    monkeypatch.setattr(RUNNER, "ROOT", project)
    monkeypatch.setattr(RUNNER.sys, "executable", str(prefix / "bin/python"))
    monkeypatch.setattr(
        RUNNER.shutil,
        "which",
        lambda name: str(prefix / "bin/git") if name == "git" else "/fixture/pixi",
    )
    monkeypatch.setenv("CONDA_PREFIX", str(prefix))
    monkeypatch.setenv("PIXI_PROJECT_ROOT", str(project))
    monkeypatch.setenv("PIXI_ENVIRONMENT_NAME", "bonsai")
    monkeypatch.setenv("PIXI_ENV_NAME", "bonsai")
    return prefix


def test_environment_guard_accepts_realpath_alias_of_project_prefix(
    bonsai_environment: Path, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    prefix = bonsai_environment
    actual = tmp_path / "resolved Bonsai environment"
    prefix.rename(actual)
    prefix.symlink_to(actual, target_is_directory=True)
    monkeypatch.setenv("CONDA_PREFIX", str(actual))
    RUNNER.require_pixi_environment()


@pytest.mark.parametrize("environment_key", ["PIXI_ENVIRONMENT_NAME", "PIXI_ENV_NAME"])
def test_environment_guard_rejects_foreign_prefix_with_matching_pixi_metadata(
    bonsai_environment: Path,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    environment_key: str,
) -> None:
    foreign = tmp_path / "foreign environment"
    (foreign / "bin").mkdir(parents=True)
    for name in ["python", "git"]:
        (foreign / "bin" / name).touch()
    monkeypatch.setenv("CONDA_PREFIX", str(foreign))
    monkeypatch.setattr(RUNNER.sys, "executable", str(foreign / "bin/python"))
    monkeypatch.setattr(
        RUNNER.shutil,
        "which",
        lambda name: str(foreign / "bin/git") if name == "git" else "/fixture/pixi",
    )
    for key in ["PIXI_ENVIRONMENT_NAME", "PIXI_ENV_NAME"]:
        monkeypatch.delenv(key)
    monkeypatch.setenv(environment_key, "bonsai")
    with pytest.raises(ValueError, match="locked Bonsai Pixi task"):
        RUNNER.require_pixi_environment()


@pytest.mark.parametrize("binary", ["python", "git"])
def test_environment_guard_rejects_native_binary_redirected_outside_prefix(
    bonsai_environment: Path, tmp_path: Path, binary: str
) -> None:
    foreign = tmp_path / f"foreign-{binary}"
    foreign.touch()
    native = bonsai_environment / "bin" / binary
    native.unlink()
    native.symlink_to(foreign)
    with pytest.raises(ValueError, match="locked Bonsai Pixi task|Git from the active Pixi"):
        RUNNER.require_pixi_environment()
