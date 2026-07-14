#!/usr/bin/env python3
"""Fail-closed audit and integrity manifest for the research-review notebooks.

This audit is intentionally lightweight. It reads notebook JSON, optional local
HTML review surfaces, and portable figure files; it does not execute a notebook
or launch any upstream bioinformatics or image-analysis program.
"""

from __future__ import annotations

import ast
import hashlib
import json
import re
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[2]
WORKSPACE = ROOT.parent
NOTEBOOK_DIR = ROOT / "notebooks/research_review"
OUTPUT = NOTEBOOK_DIR / "research_review_manifest.json"
FIGURE_DIR = ROOT / "results/figures/research_review"

EXPECTED_NOTEBOOKS = (
    "01_phylogeny_tree_trimming.ipynb",
    "02_data_tables_and_provenance.ipynb",
    "03_repeat_analysis_te34.ipynb",
    "04_ltr_deletion_footprint.ipynb",
    "05_cell_modeling_and_measurement.ipynb",
    "06_genome_size_estimation.ipynb",
    "07_cell_nucleus_genome_path_analysis.ipynb",
    "08_integrated_phylogenetic_path_analysis.ipynb",
)

VIEWER_HTML_PATHS = (
    "cellprofiler_test/output/runs/mixed_cellpose_yolo_full_dataset_v1_bgclean/linkage/index.html",
    "cellprofiler_test/output/tile_annotation_bundle_v1/index.html",
    "cellprofiler_test/output/tile_bootstrap_review_v1/index.html",
    "cellprofiler_test/output/runs/mixed_cellpose_yolo_full_dataset_v1_bgclean/linked_species_stats_reviewed/index.html",
    "cellprofiler_test/output/runs/mixed_cellpose_yolo_full_dataset_v1_bgclean/phylogenetic_analysis/index.html",
    "cellprofiler_test/output/runs/mixed_cellpose_yolo_full_dataset_v1_bgclean/integrated_multinode_phylo/index.html",
)

FIGURE_SUFFIXES = {".pdf", ".png", ".svg"}
EXPENSIVE_TOOL_PATTERN = re.compile(
    r"(?i)(?<![A-Za-z0-9_])(?:RepeatMasker|RepeatModeler|dnaPipeTE|cellpose)(?![A-Za-z0-9_])"
)
SHELL_ESCAPE_PATTERN = re.compile(
    r"(?im)^\s*(?:!|%(?:sx|system)\b|%%(?:bash|sh|script)\b)"
)
DIRECT_SHELL_TOOL_PATTERN = re.compile(
    r"(?im)^\s*!\s*(?:[\w./-]+/)?(?:RepeatMasker|RepeatModeler|dnaPipeTE|cellpose)\b"
)


def sha256(path: Path) -> str:
    """Return a streaming SHA-256 digest for *path*."""

    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _call_name(node: ast.AST) -> str:
    """Return a dotted name for a simple Python call target."""

    if isinstance(node, ast.Name):
        return node.id
    if isinstance(node, ast.Attribute):
        parent = _call_name(node.value)
        return f"{parent}.{node.attr}" if parent else node.attr
    if isinstance(node, ast.Call):
        parent = _call_name(node.func)
        return f"{parent}()" if parent else ""
    return ""


def _literal_strings(node: ast.AST) -> list[str]:
    """Collect literal strings beneath an AST node without evaluating it."""

    return [
        child.value
        for child in ast.walk(node)
        if isinstance(child, ast.Constant) and isinstance(child.value, str)
    ]


def detect_forbidden_execution(source: str) -> list[dict[str, str]]:
    """Find process-launch or Cellpose execution surfaces in Python cell code.

    Plain strings such as corrected RepeatMasker filenames are deliberately
    allowed.  They document or load frozen products and are not executions.
    Markdown is excluded by the caller.
    """

    findings: list[dict[str, str]] = []

    def add(reason: str, evidence: str) -> None:
        item = {"reason": reason, "evidence": evidence[:240]}
        if item not in findings:
            findings.append(item)

    if SHELL_ESCAPE_PATTERN.search(source):
        match = SHELL_ESCAPE_PATTERN.search(source)
        assert match is not None
        add("shell execution syntax is not allowed", match.group(0).strip())
    direct_tool = DIRECT_SHELL_TOOL_PATTERN.search(source)
    if direct_tool:
        add("expensive upstream executable invocation", direct_tool.group(0).strip())

    try:
        tree = ast.parse(source)
    except SyntaxError:
        # IPython shell syntax is handled above.  Other non-Python syntax is
        # reported by the notebook audit so it cannot silently evade checks.
        if not SHELL_ESCAPE_PATTERN.search(source):
            add("code cell is not parseable Python", "SyntaxError")
        return findings

    subprocess_aliases: set[str] = set()
    subprocess_call_aliases: set[str] = set()
    os_aliases: set[str] = set()
    os_process_aliases: set[str] = set()
    cellpose_aliases: set[str] = set()
    cellpose_symbols: set[str] = set()

    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            for alias in node.names:
                local_name = alias.asname or alias.name.split(".")[0]
                if alias.name == "subprocess" or alias.name.startswith("subprocess."):
                    subprocess_aliases.add(local_name)
                    add("subprocess import is not allowed", f"import {alias.name}")
                if alias.name == "os":
                    os_aliases.add(local_name)
                if alias.name == "cellpose" or alias.name.startswith("cellpose."):
                    cellpose_aliases.add(local_name)
                    add("Cellpose import can execute segmentation or training", f"import {alias.name}")
        elif isinstance(node, ast.ImportFrom):
            module = node.module or ""
            if module == "subprocess":
                for alias in node.names:
                    subprocess_call_aliases.add(alias.asname or alias.name)
                add("subprocess import is not allowed", f"from {module} import ...")
            if module == "os":
                process_names = {
                    "system",
                    "popen",
                    "execl",
                    "execle",
                    "execlp",
                    "execlpe",
                    "execv",
                    "execve",
                    "execvp",
                    "execvpe",
                    "spawnl",
                    "spawnle",
                    "spawnlp",
                    "spawnlpe",
                    "spawnv",
                    "spawnve",
                    "spawnvp",
                    "spawnvpe",
                }
                for alias in node.names:
                    if alias.name in process_names:
                        os_process_aliases.add(alias.asname or alias.name)
            if module == "cellpose" or module.startswith("cellpose."):
                for alias in node.names:
                    cellpose_symbols.add(alias.asname or alias.name)
                add("Cellpose import can execute segmentation or training", f"from {module} import ...")

    process_methods = {"run", "Popen", "call", "check_call", "check_output", "getoutput"}
    os_process_methods = {
        "system",
        "popen",
        "execl",
        "execle",
        "execlp",
        "execlpe",
        "execv",
        "execve",
        "execvp",
        "execvpe",
        "spawnl",
        "spawnle",
        "spawnlp",
        "spawnlpe",
        "spawnv",
        "spawnve",
        "spawnvp",
        "spawnvpe",
    }

    for node in ast.walk(tree):
        if not isinstance(node, ast.Call):
            continue
        name = _call_name(node.func)
        parts = name.split(".")
        is_process_call = (
            (len(parts) >= 2 and parts[0] in subprocess_aliases and parts[-1] in process_methods)
            or name in subprocess_call_aliases
            or (len(parts) >= 2 and parts[0] in os_aliases and parts[-1] in os_process_methods)
            or name in os_process_aliases
            or name in {"get_ipython().system", "get_ipython().getoutput"}
        )
        if is_process_call:
            literals = " ".join(_literal_strings(node))
            tool = EXPENSIVE_TOOL_PATTERN.search(literals)
            evidence = name if tool is None else f"{name}: {tool.group(0)}"
            add("external process launch is not allowed", evidence)

        root_name = parts[0] if parts else ""
        if root_name in cellpose_aliases or root_name in cellpose_symbols:
            add("Cellpose API execution is not allowed", name)
        if parts and parts[-1] in {"Cellpose", "CellposeModel"}:
            add("Cellpose model construction is not allowed", name)

    return findings


def audit_notebook(path: Path) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    """Audit one notebook and return its manifest record and failures."""

    try:
        notebook = json.loads(path.read_text())
    except (json.JSONDecodeError, UnicodeDecodeError) as exc:
        return {}, [{"gate": "valid_notebook_json", "detail": str(exc)}]

    cells = notebook.get("cells")
    if not isinstance(cells, list):
        return {}, [{"gate": "valid_notebook_schema", "detail": "cells is not a list"}]

    code_cells = [cell for cell in cells if cell.get("cell_type") == "code"]
    markdown_cells = [cell for cell in cells if cell.get("cell_type") == "markdown"]
    failures: list[dict[str, Any]] = []
    unexecuted: list[int] = []
    errors: list[dict[str, Any]] = []
    forbidden: list[dict[str, Any]] = []
    output_count = 0

    if not code_cells:
        failures.append({"gate": "has_code_cells", "detail": "no code cells"})

    for code_index, cell in enumerate(code_cells):
        if cell.get("execution_count") is None:
            unexecuted.append(code_index)
        outputs = cell.get("outputs", [])
        output_count += len(outputs)
        for output in outputs:
            if output.get("output_type") == "error":
                errors.append(
                    {
                        "code_cell_index": code_index,
                        "ename": output.get("ename"),
                        "evalue": output.get("evalue"),
                    }
                )
        source = cell.get("source", "")
        if isinstance(source, list):
            source = "".join(source)
        for finding in detect_forbidden_execution(str(source)):
            forbidden.append({"code_cell_index": code_index, **finding})

    if unexecuted:
        failures.append({"gate": "all_code_cells_executed", "detail": unexecuted})
    if errors:
        failures.append({"gate": "no_error_outputs", "detail": errors})
    if forbidden:
        failures.append(
            {"gate": "no_expensive_or_external_execution_in_code", "detail": forbidden}
        )

    kernelspec = notebook.get("metadata", {}).get("kernelspec", {})
    try:
        display_path = str(path.relative_to(ROOT))
    except ValueError:
        display_path = str(path)
    record = {
        "path": display_path,
        "sha256": sha256(path),
        "size_bytes": path.stat().st_size,
        "n_cells": len(cells),
        "n_code_cells": len(code_cells),
        "n_markdown_cells": len(markdown_cells),
        "n_outputs": output_count,
        "all_code_cells_executed": not unexecuted and bool(code_cells),
        "n_error_outputs": len(errors),
        "forbidden_execution_findings": forbidden,
        "kernel_name": kernelspec.get("name"),
        "status": "pass" if not failures else "fail",
    }
    return record, failures


def _file_record(path: Path, display_path: str) -> dict[str, Any]:
    return {
        "path": display_path,
        "sha256": sha256(path),
        "size_bytes": path.stat().st_size,
    }


def build_manifest() -> dict[str, Any]:
    """Audit the full bundle and return a release manifest or raise."""

    failures: list[dict[str, Any]] = []
    expected_paths = [NOTEBOOK_DIR / name for name in EXPECTED_NOTEBOOKS]
    missing_notebooks = [name for name, path in zip(EXPECTED_NOTEBOOKS, expected_paths) if not path.exists()]
    if missing_notebooks:
        failures.append({"gate": "expected_notebook_set_complete", "detail": missing_notebooks})

    visible_notebooks = {
        path.name for path in NOTEBOOK_DIR.glob("*.ipynb") if not path.name.startswith(".")
    }
    unexpected_notebooks = sorted(visible_notebooks - set(EXPECTED_NOTEBOOKS))
    if unexpected_notebooks:
        failures.append({"gate": "no_unexpected_notebooks", "detail": unexpected_notebooks})

    notebook_records: dict[str, dict[str, Any]] = {}
    for name, path in zip(EXPECTED_NOTEBOOKS, expected_paths):
        if not path.exists():
            continue
        record, notebook_failures = audit_notebook(path)
        notebook_records[name] = record
        for failure in notebook_failures:
            failures.append({"notebook": name, **failure})

    viewer_records: list[dict[str, Any]] = []
    missing_viewers: list[str] = []
    for relative in VIEWER_HTML_PATHS:
        path = WORKSPACE / relative
        if not path.exists():
            missing_viewers.append(relative)
            continue
        viewer_records.append(_file_record(path, f"../{relative}"))
    figure_records: list[dict[str, Any]] = []
    if FIGURE_DIR.exists():
        for path in sorted(FIGURE_DIR.rglob("*")):
            if path.is_file() and path.suffix.lower() in FIGURE_SUFFIXES:
                figure_records.append(_file_record(path, str(path.relative_to(ROOT))))
    if not figure_records:
        failures.append({"gate": "review_figures_present", "detail": "no PNG, PDF, or SVG files"})

    if failures:
        summary = json.dumps(failures, indent=2)
        raise RuntimeError(f"Research-review notebook audit failed:\n{summary}")

    all_executed = all(
        record["all_code_cells_executed"] for record in notebook_records.values()
    )
    no_errors = all(record["n_error_outputs"] == 0 for record in notebook_records.values())
    no_forbidden = all(
        not record["forbidden_execution_findings"] for record in notebook_records.values()
    )

    return {
        "schema_version": 1,
        "analysis_id": "Desmognathus_research_review_notebooks_v1",
        "audited_at_utc": datetime.now(timezone.utc).isoformat(),
        "status": "pass",
        "audit_scope": (
            "Static integrity audit of executed notebook code/output, optional local HTML "
            "review surfaces, and portable review figures; no notebook or upstream tool was run."
        ),
        "gates": {
            "expected_notebook_set_complete": len(notebook_records) == len(EXPECTED_NOTEBOOKS),
            "no_unexpected_notebooks": True,
            "all_code_cells_executed": all_executed,
            "no_error_outputs": no_errors,
            "no_expensive_or_external_execution_detected_in_code": no_forbidden,
            "optional_viewer_availability_recorded": True,
            "notebook_sha256_recorded": True,
            "review_figures_present": bool(figure_records),
            "figure_sha256_recorded": bool(figure_records),
        },
        # Required compatibility field.  This statement is scoped to notebook
        # code by audit_scope and the more precise gate immediately above.
        "expensive_upstream_tools_executed": False,
        "expensive_upstream_execution_detection_scope": "notebook code cells",
        "forbidden_execution_evidence": [],
        "notebooks": notebook_records,
        "optional_viewer_html": viewer_records,
        "optional_viewer_html_missing": missing_viewers,
        "figure_outputs": figure_records,
    }


def main() -> None:
    OUTPUT.unlink(missing_ok=True)
    manifest = build_manifest()
    temporary = OUTPUT.with_suffix(OUTPUT.suffix + ".tmp")
    temporary.write_text(json.dumps(manifest, indent=2) + "\n")
    temporary.replace(OUTPUT)
    print(
        "Research-review notebook audit passed: "
        f"{len(manifest['notebooks'])} notebooks, "
        f"{len(manifest['optional_viewer_html'])} optional viewers, "
        f"{len(manifest['figure_outputs'])} figures"
    )
    print(f"Manifest: {OUTPUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
