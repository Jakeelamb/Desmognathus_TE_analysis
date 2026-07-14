#!/usr/bin/env python3
"""Publish the full fitted path analysis as an executed notebook and index.html."""

from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import nbformat as nbf


PROJECT_ROOT = Path(__file__).resolve().parents[2]
SCRIPT_DIR = PROJECT_ROOT / "path_analysis" / "scripts"
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import build_cell_nucleus_genome_path_presentation as analysis  # noqa: E402


OUTPUT_DIR = PROJECT_ROOT / "notebooks/research_review/cell_nucleus_genome_path_analysis"
NOTEBOOK_PATH = analysis.NOTEBOOK_PATH
INDEX_HTML_PATH = OUTPUT_DIR / "index.html"
MANIFEST_PATH = OUTPUT_DIR / "presentation_manifest.json"


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def execute_notebook(notebook: nbf.NotebookNode) -> nbf.NotebookNode:
    from nbclient import NotebookClient

    client = NotebookClient(
        notebook,
        timeout=900,
        kernel_name="python3",
        resources={"metadata": {"path": str(PROJECT_ROOT)}},
        allow_errors=False,
    )
    return client.execute()


def export_html(notebook: nbf.NotebookNode) -> None:
    from nbconvert import HTMLExporter
    from traitlets.config import Config

    config = Config()
    config.HTMLExporter.exclude_input_prompt = True
    config.HTMLExporter.exclude_output_prompt = True
    exporter = HTMLExporter(config=config, template_name="lab")
    html, _ = exporter.from_notebook_node(notebook)
    INDEX_HTML_PATH.write_text(html, encoding="utf-8")


def validate_executed_notebook(notebook: nbf.NotebookNode) -> None:
    code_cells = [cell for cell in notebook.cells if cell.cell_type == "code"]
    errors = [
        output
        for cell in code_cells
        for output in cell.get("outputs", [])
        if output.get("output_type") == "error"
    ]
    if errors:
        raise RuntimeError(f"Executed notebook contains {len(errors)} error output(s)")
    if not code_cells or any(cell.get("execution_count") is None for cell in code_cells):
        raise RuntimeError("Every path-analysis code cell must be executed")
    cell_ids = [cell.get("id") for cell in notebook.cells]
    if len(cell_ids) != len(set(cell_ids)) or any(not value for value in cell_ids):
        raise RuntimeError("Notebook cell IDs must be present and unique")


def write_manifest(notebook: nbf.NotebookNode) -> None:
    source_manifest = json.loads(analysis.ANALYSIS_MANIFEST_PATH.read_text())
    figure_manifest = json.loads(analysis.FIGURE_MANIFEST_PATH.read_text())
    manifest = {
        "analysis_id": "cell_nucleus_genome_phylogenetic_path_presentation_v1",
        "status": "pass",
        "canonical_notebook": str(NOTEBOOK_PATH.relative_to(PROJECT_ROOT)),
        "canonical_html": str(INDEX_HTML_PATH.relative_to(PROJECT_ROOT)),
        "notebook_sha256": sha256_file(NOTEBOOK_PATH),
        "html_sha256": sha256_file(INDEX_HTML_PATH),
        "analysis_manifest": str(
            analysis.ANALYSIS_MANIFEST_PATH.relative_to(PROJECT_ROOT)
        ),
        "analysis_manifest_sha256": sha256_file(analysis.ANALYSIS_MANIFEST_PATH),
        "figure_manifest": str(analysis.FIGURE_MANIFEST_PATH.relative_to(PROJECT_ROOT)),
        "figure_manifest_sha256": sha256_file(analysis.FIGURE_MANIFEST_PATH),
        "reference_method_audit": str(
            analysis.REFERENCE_AUDIT_PATH.relative_to(PROJECT_ROOT)
        ),
        "gates": {
            "primary_species": source_manifest["primary_species"],
            "labeled_dags": source_manifest["labeled_dags"],
            "markov_equivalence_classes": source_manifest[
                "markov_equivalence_classes"
            ],
            "measurement_bootstrap_replicates": source_manifest[
                "measurement_bootstrap_replicates"
            ],
            "published_bootstrap_trees": source_manifest[
                "published_bootstrap_trees_analyzed"
            ],
            "simulation_replicates_per_class_and_regime": source_manifest[
                "simulation_replicates_per_class_and_regime"
            ],
            "fit_failures": source_manifest["failure_count"],
            "causal_direction_identified": source_manifest[
                "causal_direction_identified"
            ],
            "all_code_cells_executed": all(
                cell.get("execution_count") is not None
                for cell in notebook.cells
                if cell.cell_type == "code"
            ),
            "no_error_outputs": True,
            "figure_release_gates_passed": figure_manifest[
                "release_gates_passed"
            ],
        },
    }
    MANIFEST_PATH.write_text(json.dumps(manifest, indent=2) + "\n")


def build() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    analysis.build()
    conclusion = json.loads(analysis.CONCLUSION_PATH.read_text())
    source_notebook = analysis.build_notebook(conclusion)
    executed = execute_notebook(source_notebook)
    validate_executed_notebook(executed)
    nbf.write(executed, NOTEBOOK_PATH)
    export_html(executed)
    write_manifest(executed)
    print(f"Notebook: {NOTEBOOK_PATH}")
    print(f"HTML: {INDEX_HTML_PATH}")


if __name__ == "__main__":
    build()
