"""Run the frozen Bonsai preparation and analysis through locked Pixi environments."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
COMPARISON = ROOT / "analyses/01_transposable_elements/bonsai_comparison"
UPSTREAM_URL = "https://github.com/dhdegroot/Bonsai-data-representation.git"
UPSTREAM_COMMIT = "b16ef6f0a164610327022c7d0b92de82b3596aad"
REQUIREMENTS = {
    "matplotlib==3.11.0",
    "numpy==2.4.6",
    "scipy==1.17.1",
    "pandas==3.0.3",
    "scikit-learn==1.9.0",
    "ruamel.yaml==0.19.1",
    "psutil==7.2.2",
}
RUN_MANIFEST = "bonsai_pixi_run.json"
BOOTSTRAP_PLAN = "te34_bonsai_feature_bootstrap_plan.json"
BOOTSTRAP_DRAWS = "te34_bonsai_feature_bootstrap_draws.csv"
SEED = 20260826


def execute(argv: list[str], **kwargs: object) -> None:
    subprocess.run(argv, check=True, **kwargs)


def git_output(upstream: Path, *args: str) -> str:
    return subprocess.check_output(["git", "-C", str(upstream), *args], text=True).strip()


def verify_upstream(upstream: Path) -> Path:
    upstream = upstream.resolve()
    if not upstream.is_dir():
        raise ValueError(f"Upstream source directory is missing: {upstream}")
    if Path(git_output(upstream, "rev-parse", "--show-toplevel")).resolve() != upstream:
        raise ValueError(f"Upstream must be the root of its Git checkout: {upstream}")
    if git_output(upstream, "rev-parse", "HEAD") != UPSTREAM_COMMIT:
        raise ValueError(f"Upstream checkout is not pinned to {UPSTREAM_COMMIT}: {upstream}")
    if git_output(upstream, "status", "--porcelain", "--untracked-files=all"):
        raise ValueError(f"Refusing to use a dirty upstream checkout: {upstream}")
    requirements = {
        line.strip()
        for line in (upstream / "requirements.txt").read_text().splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    }
    if requirements != REQUIREMENTS:
        raise ValueError("Upstream requirements do not match the seven locked Bonsai pins")
    if not (upstream / "bonsai/bonsai_main.py").is_file():
        raise ValueError("Upstream Bonsai entry point is missing")
    return upstream


def get_upstream(explicit: Path | None) -> Path:
    if explicit is not None:
        return verify_upstream(explicit)
    cache = ROOT / ".pixi/sources" / f"bonsai-{UPSTREAM_COMMIT}"
    if not cache.resolve().is_relative_to(ROOT / ".pixi/sources"):
        raise ValueError("Source cache must not redirect writes outside .pixi/sources")
    if cache.exists() or cache.is_symlink():
        return verify_upstream(cache)
    cache.mkdir(parents=True, exist_ok=False)
    # A failed acquisition is retained. Never repair or replace an existing cache.
    execute(["git", "init", "-q", str(cache)])
    execute(["git", "-C", str(cache), "remote", "add", "origin", UPSTREAM_URL])
    execute(["git", "-C", str(cache), "fetch", "-q", "--depth", "1", "origin", UPSTREAM_COMMIT])
    execute(["git", "-C", str(cache), "checkout", "-q", "--detach", "FETCH_HEAD"])
    return verify_upstream(cache)


def require_pixi_environment() -> None:
    prefix = Path(os.environ.get("CONDA_PREFIX", "/missing-pixi-environment")).resolve()
    if (
        os.environ.get("PIXI_ENVIRONMENT_NAME", os.environ.get("PIXI_ENV_NAME")) != "bonsai"
        or Path(os.environ.get("PIXI_PROJECT_ROOT", "/missing-pixi-project")).resolve() != ROOT
        or prefix != (ROOT / ".pixi/envs/bonsai").resolve()
        or not Path(sys.executable).resolve().is_relative_to(prefix)
    ):
        raise ValueError("Launch this runner through the project's locked Bonsai Pixi task")
    git = shutil.which("git")
    if git is None or not Path(git).resolve().is_relative_to(prefix):
        raise ValueError("Bonsai requires Git from the active Pixi environment")
    if shutil.which(os.environ.get("PIXI_EXE", "pixi")) is None:
        raise ValueError("Pixi is required for locked default-environment preparation and analysis")
    if not (ROOT / ".pixi/envs/default/bin/python").is_file():
        raise ValueError("Install the locked default Pixi environment before running Bonsai")


def safe_run_directory(path: Path, upstream: Path | None = None) -> Path:
    run_dir = path.expanduser().resolve()
    if ".git" in run_dir.parts:
        raise ValueError("Run directories may not be inside .git")
    if ROOT.is_relative_to(run_dir):
        raise ValueError("Run directory may not contain the project")
    if run_dir.is_relative_to(ROOT) and not run_dir.is_relative_to(ROOT / ".pixi/runs"):
        raise ValueError("Project scratch runs must be below .pixi/runs; science paths are frozen")
    if upstream is not None and (
        run_dir.is_relative_to(upstream) or upstream.is_relative_to(run_dir)
    ):
        raise ValueError("Run directory must not overlap the upstream source checkout")
    return run_dir


def check_resume_tree(run_dir: Path) -> None:
    # A reused result/log directory must not redirect writes outside this owned run.
    for path in [run_dir, *run_dir.rglob("*")]:
        if path.is_symlink():
            raise ValueError(f"Resume directory contains a symlink: {path}")
        if path.is_file() and path.stat().st_nlink != 1:
            raise ValueError(f"Resume directory contains a hard-linked file: {path}")


def prepared_files(run_dir: Path, mode: str, replicates: int) -> dict[str, str]:
    files = []
    if mode == "bootstrap":
        files.extend([run_dir / BOOTSTRAP_PLAN, run_dir / BOOTSTRAP_DRAWS])
        directories = [
            run_dir / "bootstrap" / f"replicate_{index:03d}" for index in range(replicates + 1)
        ]
    else:
        directories = [run_dir]
    for directory in directories:
        files.extend([directory / "bonsai_config.yaml"])
        files.extend(
            directory / "input" / name for name in ("features.txt", "cellID.txt", "geneID.txt")
        )
    return {
        str(path.relative_to(run_dir)): hashlib.sha256(path.read_bytes()).hexdigest()
        for path in files
    }


def manifest_identity(run_dir: Path, mode: str, replicates: int) -> dict[str, object]:
    return {
        "schema": "desmognathus-pixi-bonsai-run-v1",
        "run_dir": str(run_dir),
        "mode": mode,
        "upstream_commit": UPSTREAM_COMMIT,
        "replicates": replicates if mode == "bootstrap" else None,
        "seed": SEED if mode == "bootstrap" else None,
    }


def local_helper(name: str, *args: str) -> None:
    execute(
        [
            os.environ.get("PIXI_EXE", "pixi"),
            "run",
            "--manifest-path",
            str(ROOT / "pixi.toml"),
            "--environment",
            "default",
            "--clean-env",
            "--locked",
            "--no-install",
            "--executable",
            "python",
            str(COMPARISON / name),
            *args,
        ],
        cwd=ROOT,
        # Pixi's host glibc probe needs system tools even in --clean-env.
        env={**os.environ, "PATH": os.environ.get("PATH", "") + os.pathsep + os.defpath},
    )


def prepare(run_dir: Path, mode: str, replicates: int, resume: bool) -> None:
    identity = manifest_identity(run_dir, mode, replicates)
    if resume:
        check_resume_tree(run_dir)
        saved = json.loads((run_dir / RUN_MANIFEST).read_text())
        if saved != {**identity, "prepared_files": prepared_files(run_dir, mode, replicates)}:
            raise ValueError("Resume requires an unchanged prepared manifest from this runner")
    if mode == "refit":
        local_helper("compare_bonsai_to_pca.py", "prepare", "--run-dir", str(run_dir))
    else:
        local_helper(
            "validate_strategy_localization.py",
            "prepare-bootstrap",
            "--run-dir",
            str(run_dir),
            "--replicates",
            str(replicates),
            "--seed",
            str(SEED),
        )
    if not resume:
        saved = {**identity, "prepared_files": prepared_files(run_dir, mode, replicates)}
        (run_dir / RUN_MANIFEST).write_text(json.dumps(saved, indent=2) + "\n")


def run_one(upstream: Path, directory: Path, bootstrap: bool) -> None:
    tree = directory / "results/final_bonsai/tree.nwk"
    if bootstrap and tree.is_file() and tree.stat().st_size:
        print(f"SKIP {directory.name}", flush=True)
        return
    env = os.environ.copy()
    env["PYTHONDONTWRITEBYTECODE"] = "1"
    env.pop("PYTHONPATH", None)
    env.pop("PYTHONHOME", None)
    if bootstrap:
        env.update(
            {
                name: "1"
                for name in (
                    "OPENBLAS_NUM_THREADS",
                    "OMP_NUM_THREADS",
                    "MKL_NUM_THREADS",
                    "NUMEXPR_NUM_THREADS",
                )
            }
        )
    with (directory / "bonsai.log").open("w") as log:
        execute(
            [
                sys.executable,
                str(upstream / "bonsai/bonsai_main.py"),
                "--config_filepath",
                str(directory / "bonsai_config.yaml"),
                "--step",
                "all",
            ],
            cwd=upstream,
            env=env,
            stdout=log,
            stderr=subprocess.STDOUT,
        )
    if not tree.is_file() or not tree.stat().st_size:
        raise ValueError(f"Upstream did not produce a nonempty final tree: {tree}")
    print(f"DONE {directory.name}", flush=True)


def run(args: argparse.Namespace) -> None:
    require_pixi_environment()
    run_dir = safe_run_directory(args.run_dir)
    if args.resume and args.mode != "bootstrap":
        raise ValueError("Only bootstrap campaigns support explicit resume")
    if run_dir.exists() and not args.resume:
        raise ValueError(f"Refusing to overwrite an existing run directory: {run_dir}")
    if args.resume and not (run_dir / RUN_MANIFEST).is_file():
        raise ValueError("Resume requires an existing run prepared by this runner")
    upstream = get_upstream(args.upstream)
    safe_run_directory(run_dir, upstream)
    if not args.resume:
        # Exclusive mkdir protects existing directories even if empty or unrelated.
        run_dir.mkdir(parents=True, exist_ok=False)
    lock = run_dir / ".bonsai-run.lock"
    with lock.open("x") as handle:
        handle.write(f"{os.getpid()}\n")
        handle.flush()
        try:
            prepare(run_dir, args.mode, args.replicates, args.resume)
            run_prepared(args, run_dir, upstream)
        finally:
            lock.unlink()


def run_prepared(args: argparse.Namespace, run_dir: Path, upstream: Path) -> None:
    if args.mode == "refit":
        run_one(upstream, run_dir, False)
        destination = run_dir / "comparison"
        local_helper(
            "compare_bonsai_to_pca.py",
            "analyze",
            "--tree",
            str(run_dir / "results/final_bonsai/tree.nwk"),
            "--metadata",
            str(run_dir / "results/final_bonsai/metadata_bonsai.txt"),
            "--output-dir",
            str(destination),
            "--write",
        )
        plot = "plot_bonsai_comparison.R"
    else:
        directories = [
            run_dir / "bootstrap" / f"replicate_{index:03d}" for index in range(args.replicates + 1)
        ]
        with ThreadPoolExecutor(max_workers=args.jobs) as pool:
            list(pool.map(lambda directory: run_one(upstream, directory, True), directories))
        if args.replicates != 200:
            print(
                f"PASS: execution smoke only; calibration plus {args.replicates} replicates "
                f"under {run_dir}; full-campaign analysis and plots were not run"
            )
            return
        destination = run_dir / "validation"
        local_helper(
            "validate_strategy_localization.py",
            "analyze-bootstrap",
            "--run-dir",
            str(run_dir),
            "--output-dir",
            str(destination),
            "--write",
        )
        plot = "plot_bonsai_validation.R"
    execute(
        [str(ROOT / "scripts/run_r.sh"), "Rscript", str(COMPARISON / plot), str(destination)],
        cwd=ROOT,
    )
    print(f"PASS: completed pinned Bonsai {args.mode} under {run_dir}")


def positive_int(value: str) -> int:
    number = int(value)
    if number < 1:
        raise argparse.ArgumentTypeError("must be at least 1")
    return number


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=("refit", "bootstrap"))
    parser.add_argument(
        "--run-dir",
        type=Path,
        required=True,
        help="new scratch directory; inside this project use .pixi/runs",
    )
    parser.add_argument(
        "--upstream", type=Path, help="verify an existing pinned checkout read-only"
    )
    parser.add_argument(
        "--replicates",
        type=positive_int,
        default=200,
        help="bootstrap replicates; values other than 200 are execution smoke only",
    )
    parser.add_argument("--jobs", type=positive_int, default=8)
    parser.add_argument(
        "--resume", action="store_true", help="resume this runner's bootstrap manifest"
    )
    args = parser.parse_args(argv)
    try:
        run(args)
    except (ValueError, OSError, subprocess.CalledProcessError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
