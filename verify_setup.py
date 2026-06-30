#!/usr/bin/env python3
"""
Read-only setup verification for the Desmognathus TE analysis repository.

This script checks the active environment, local input-data contract, path
configuration, and script presence. It must not create result directories or
write generated analysis artifacts.
"""

from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Callable

import yaml


PROJECT_ROOT = Path(__file__).resolve().parent

GREEN = "\033[92m"
RED = "\033[91m"
YELLOW = "\033[93m"
BLUE = "\033[94m"
RESET = "\033[0m"


def print_header(text: str) -> None:
    print(f"\n{BLUE}{'=' * 60}{RESET}")
    print(f"{BLUE}{text}{RESET}")
    print(f"{BLUE}{'=' * 60}{RESET}")


def print_status(item: str, status: bool, note: str = "") -> None:
    mark = f"{GREEN}OK{RESET}" if status else f"{RED}FAIL{RESET}"
    note_str = f" {YELLOW}({note}){RESET}" if note else ""
    print(f"  {mark} {item}{note_str}")


def rel(path: Path) -> str:
    try:
        return str(path.relative_to(PROJECT_ROOT))
    except ValueError:
        return str(path)


def path_is_relative_to(path: Path, root: Path) -> bool:
    try:
        path.resolve().relative_to(root.resolve())
        return True
    except ValueError:
        return False


def check_execution_environment() -> bool:
    print_header("Execution Environment")

    executable = Path(sys.executable)
    print(f"  Python executable: {executable}")
    print(f"  Python version: {sys.version.split()[0]}")

    conda_prefix = os.environ.get("CONDA_PREFIX")
    conda_env = os.environ.get("CONDA_DEFAULT_ENV")
    if not conda_prefix:
        print_status("CONDA_PREFIX", False, "not set; use scripts/run_in_dusky.sh for reproducible local runs")
        return False

    prefix_path = Path(conda_prefix)
    env_label = conda_env or prefix_path.name
    print(f"  Conda environment: {env_label}")
    print(f"  Conda prefix: {prefix_path}")

    python_matches_prefix = path_is_relative_to(executable, prefix_path)
    print_status(
        "Python resolves inside CONDA_PREFIX",
        python_matches_prefix,
        "PATH likely needs CONDA_PREFIX/bin first" if not python_matches_prefix else "",
    )
    if not python_matches_prefix:
        print("    Hint: scripts/run_in_dusky.sh prepends CONDA_PREFIX/bin before running commands.")
    return python_matches_prefix


def check_python_packages() -> bool:
    print_header("Python Packages")

    required_packages = [
        ("pandas", "pandas"),
        ("numpy", "numpy"),
        ("matplotlib", "matplotlib"),
        ("seaborn", "seaborn"),
        ("tqdm", "tqdm"),
        ("yaml", "yaml"),
        ("scipy", "scipy"),
    ]
    optional_packages = [
        ("biopython", "Bio", "Used by genome and sequence helper utilities"),
        ("dask", "dask", "Used by divergence.py"),
    ]

    all_ok = True
    for display_name, import_name in required_packages:
        try:
            __import__(import_name)
            print_status(display_name, True)
        except ImportError:
            print_status(display_name, False, "not installed")
            all_ok = False

    print(f"\n  {YELLOW}Optional/runtime-specific packages:{RESET}")
    for display_name, import_name, note in optional_packages:
        try:
            __import__(import_name)
            print_status(display_name, True, note)
        except ImportError:
            print_status(display_name, False, note)
            all_ok = False

    return all_ok


def check_r_availability() -> bool:
    print_header("R Environment")

    rscript = shutil.which("Rscript")
    print_status("Rscript available", rscript is not None)
    if not rscript:
        return False

    try:
        result = subprocess.run(
            [rscript, "--version"],
            capture_output=True,
            text=True,
            timeout=5,
            check=False,
        )
        version_output = result.stderr.strip() or result.stdout.strip()
        version = version_output.splitlines()[0] if version_output else "unknown"
        print(f"    Version: {version}")
    except Exception as exc:
        print_status("Rscript version probe", False, str(exc))
        return False

    all_ok = True
    conda_prefix = os.environ.get("CONDA_PREFIX")
    if conda_prefix:
        try:
            lib_result = subprocess.run(
                [
                    rscript,
                    "-e",
                    "cat(paste(normalizePath(.libPaths(), mustWork=FALSE), collapse='\\n'))",
                ],
                capture_output=True,
                text=True,
                timeout=10,
                check=False,
            )
            if lib_result.returncode != 0:
                print_status("R library path probe", False, lib_result.stderr.strip())
                all_ok = False
            else:
                lib_paths = [Path(line.strip()) for line in lib_result.stdout.splitlines() if line.strip()]
                outside_prefix = [
                    str(path)
                    for path in lib_paths
                    if not path_is_relative_to(path, Path(conda_prefix))
                ]
                print_status(
                    "R libraries resolve inside CONDA_PREFIX",
                    len(outside_prefix) == 0,
                    "; ".join(outside_prefix),
                )
                all_ok = all_ok and len(outside_prefix) == 0
        except Exception as exc:
            print_status("R library path probe", False, str(exc))
            all_ok = False

    try:
        pkg_result = subprocess.run(
            [
                rscript,
                "-e",
                "required <- c('ape'); missing <- required[!vapply(required, requireNamespace, logical(1), quietly=TRUE)]; if (length(missing)) stop(paste(missing, collapse=', '))",
            ],
            capture_output=True,
            text=True,
            timeout=15,
            check=False,
        )
        print_status("Required R packages", pkg_result.returncode == 0, pkg_result.stderr.strip())
        all_ok = all_ok and pkg_result.returncode == 0
    except Exception as exc:
        print_status("Required R packages", False, str(exc))
        all_ok = False

    return all_ok


def check_external_tools() -> bool:
    print_header("External Tools")

    tools = ["bowtie2", "samtools", "bedtools", "RepeatMasker", "trf", "TEsorter"]
    all_ok = True
    for tool in tools:
        found = shutil.which(tool) is not None
        print_status(tool, found)
        all_ok = all_ok and found
    return all_ok


def check_path_config() -> bool:
    print_header("Path Configuration")

    root_config = PROJECT_ROOT / "paths.yaml"
    if not root_config.exists():
        print_status(rel(root_config), False)
        return False
    print_status(rel(root_config), True)

    with root_config.open("r", encoding="utf-8") as handle:
        root_data = yaml.safe_load(handle)

    all_ok = True
    required_keys = [
        ("input_data", "genomes"),
        ("input_data", "lookup_table"),
        ("results", "root"),
        ("results", "data"),
        ("results", "reports"),
    ]
    for section, key in required_keys:
        present = section in root_data and key in root_data[section]
        print_status(f"paths.{section}.{key}", present)
        all_ok = all_ok and present

    return all_ok


def manifest_entry_present(entry: dict[str, object]) -> tuple[bool, str]:
    path_text = str(entry.get("path", "")).strip()
    kind = str(entry.get("kind", "")).strip()
    if not path_text:
        return False, "missing path"
    if "*" in path_text or kind.endswith("_glob"):
        matches = sorted(PROJECT_ROOT.glob(path_text))
        return len(matches) > 0, f"{len(matches)} match(es)"
    full_path = PROJECT_ROOT / path_text
    return full_path.exists(), "present" if full_path.exists() else "missing"


def check_data_manifest(*, validate_paths: bool = True) -> bool:
    print_header("Data Manifest")

    manifest_path = PROJECT_ROOT / "DATA_MANIFEST.yml"
    exists = manifest_path.exists()
    print_status("DATA_MANIFEST.yml", exists)
    if not exists:
        return False

    with manifest_path.open("r", encoding="utf-8") as handle:
        manifest = yaml.safe_load(handle)

    required_sections = [
        "required_local_inputs",
        "path_analysis_imports",
        "tracked_small_external_snapshots",
        "generated_outputs",
        "workflow_steps",
    ]
    all_ok = True
    for section in required_sections:
        present = section in manifest and isinstance(manifest[section], list)
        print_status(section, present)
        all_ok = all_ok and present
    if not validate_paths:
        print_status("manifest path validation", True, "skipped by --skip-data")
        return all_ok

    print("\n  Required manifest paths:")
    for section in ["required_local_inputs", "path_analysis_imports", "tracked_small_external_snapshots"]:
        for entry in manifest.get(section, []):
            ok, note = manifest_entry_present(entry)
            print_status(f"{section}: {entry.get('path')}", ok, note)
            all_ok = all_ok and ok
    return all_ok


def check_directory_structure() -> bool:
    print_header("Directory Structure")

    required_dirs = {
        "Input Data": [
            "input_data",
            "input_data/dnaPipeTE",
            "input_data/repeatmasker",
            "input_data/phylogeny",
        ],
        "Scripts": [
            "scripts/processing",
            "scripts/visualization",
            "path_analysis/scripts",
        ],
    }

    all_ok = True
    for category, dirs in required_dirs.items():
        print(f"\n  {category}:")
        for dir_path in dirs:
            full_path = PROJECT_ROOT / dir_path
            exists = full_path.exists() and full_path.is_dir()
            print_status(dir_path, exists)
            all_ok = all_ok and exists

    return all_ok


def check_input_files() -> bool:
    print_header("Critical Input Files")

    files_to_check = {
        "input_data/lookup_table.txt": "Species lookup table",
        "input_data/phylogeny/desmo900dated_test.tre": "Phylogenetic tree",
    }
    data_dirs = {
        "input_data/dnaPipeTE": "dnaPipeTE files",
        "input_data/repeatmasker": ".align files",
    }

    all_ok = True
    for file_path, description in files_to_check.items():
        full_path = PROJECT_ROOT / file_path
        exists = full_path.exists() and full_path.is_file()
        print_status(f"{description} ({file_path})", exists)
        all_ok = all_ok and exists

    print("\n  Data Directory Contents:")
    for dir_path, description in data_dirs.items():
        full_path = PROJECT_ROOT / dir_path
        if full_path.exists() and full_path.is_dir():
            file_count = sum(1 for child in full_path.iterdir() if child.is_file())
            status = file_count > 0
            note = f"{file_count} files" if status else "empty"
            print_status(f"{description} ({dir_path})", status, note)
            all_ok = all_ok and status
        else:
            print_status(f"{description} ({dir_path})", False, "directory not found")
            all_ok = False

    return all_ok


def check_scripts() -> bool:
    print_header("Processing Scripts")

    scripts = [
        ("scripts/processing/dnaPipe.py", "dnaPipeTE processing"),
        ("scripts/processing/repeatmask.py", "RepeatMasker processing"),
        ("scripts/processing/ec.py", "Ectopic recombination"),
        ("scripts/processing/divergence.py", "Divergence analysis"),
        ("scripts/processing/diversity_stats.py", "Canonical diversity metrics"),
        ("scripts/processing/pca.R", "PCA analysis"),
        ("scripts/processing/clean_tree_phylo.R", "Phylogeny cleaning"),
        ("scripts/processing/phylogenetic_pca_analysis.R", "Phylogenetic PCA"),
        ("scripts/processing/analyze_phylogenetic_signal.R", "Phylogenetic signal"),
        ("scripts/processing/analyze_phylogenetic_correlogram.R", "Phylogenetic correlogram"),
        ("path_analysis/scripts/build_master_dataset.py", "Path-analysis master dataset"),
        ("path_analysis/scripts/build_analysis_panels.py", "Path-analysis panels"),
        ("path_analysis/scripts/path_model_scaffold.R", "Path-analysis model scaffold"),
        ("path_analysis/scripts/audit_source_traceability.py", "Path-analysis source traceability"),
    ]

    all_ok = True
    for script_path, description in scripts:
        full_path = PROJECT_ROOT / script_path
        exists = full_path.exists() and full_path.is_file()
        print_status(f"{description} ({script_path})", exists)
        all_ok = all_ok and exists

    return all_ok


def check_existing_outputs() -> None:
    print_header("Existing Output Files")

    results_data = PROJECT_ROOT / "results" / "data"
    results_tables = PROJECT_ROOT / "results" / "tables"

    if not results_data.exists():
        print("  results/data does not exist yet")
        return

    csv_files = list(results_data.glob("*.csv"))
    tre_files = list(results_data.glob("*.tre"))

    print(f"\n  CSV files: {len(csv_files)}")
    print(f"  Tree files: {len(tre_files)}")
    if results_tables.exists():
        table_csv_files = list(results_tables.rglob("*.csv"))
        print(f"  Table CSV files: {len(table_csv_files)}")

    if csv_files or tre_files:
        print(f"\n  {YELLOW}Note: some analyses may have been run already{RESET}")


def run_checks(checks: list[tuple[str, Callable[[], bool]]]) -> dict[str, bool]:
    results: dict[str, bool] = {}
    for check_name, check_func in checks:
        try:
            results[check_name] = bool(check_func())
        except Exception as exc:
            print(f"\n{RED}Error running {check_name}: {exc}{RESET}")
            results[check_name] = False
    return results


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--json",
        action="store_true",
        help="Print a machine-readable summary as the final line.",
    )
    parser.add_argument(
        "--skip-data",
        action="store_true",
        help="Skip local input-data checks; useful for source-only CI.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    print(f"\n{BLUE}{'=' * 60}")
    print("Desmognathus TE Analysis - Setup Verification")
    print(f"{'=' * 60}{RESET}\n")

    checks: list[tuple[str, Callable[[], bool]]] = [
        ("Execution Environment", check_execution_environment),
        ("Path Configuration", check_path_config),
        ("Data Manifest", lambda: check_data_manifest(validate_paths=not args.skip_data)),
        ("Python Packages", check_python_packages),
        ("R Environment", check_r_availability),
        ("External Tools", check_external_tools),
        ("Directory Structure", check_directory_structure),
        ("Scripts", check_scripts),
    ]
    if not args.skip_data:
        checks.append(("Critical Input Files", check_input_files))

    results = run_checks(checks)
    check_existing_outputs()

    print_header("Summary")
    all_passed = all(results.values())
    for check_name, passed in results.items():
        print_status(check_name, passed)

    if all_passed:
        print(f"\n{GREEN}All critical checks passed.{RESET}\n")
        exit_code = 0
    else:
        print(f"\n{RED}Some checks failed.{RESET}")
        print(f"{YELLOW}Address the failed checks before running analyses.{RESET}\n")
        exit_code = 1

    if args.json:
        print(json.dumps({"ok": all_passed, "checks": results}, sort_keys=True))

    return exit_code


if __name__ == "__main__":
    sys.exit(main())
