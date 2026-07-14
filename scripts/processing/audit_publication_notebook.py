#!/usr/bin/env python3
"""Fail-closed integrity check and manifest for the executed audit notebook."""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
NOTEBOOK = ROOT / "notebooks/Desmognathus_publication_audit_analysis18_v1.ipynb"
OUTPUT = ROOT / "notebooks/Desmognathus_publication_audit_analysis18_v1.manifest.json"

REQUIRED_ARTIFACTS = [
    "results/data/corrected/dnapipete/dnapipete_mass_accounting_analysis18_v2.csv",
    "results/data/corrected/ectopic/ectopic_species_robustness_analysis18_v1.csv",
    "results/data/corrected/diversity_pca/te_diversity_mass_sensitivity_analysis18_v1.csv",
    "results/data/corrected/phylogeny/phylogeny_tree_metrics_analysis18_v1.csv",
    "results/data/corrected/microscopy/microscopy_panel_support_analysis18_v1.csv",
    "results/data/corrected/microscopy/microscopy_segmentation_validation_analysis18_v1.csv",
    "results/data/corrected/path_analysis/corrected_path_anchor_model_comparison_analysis18_v1.csv",
    "results/data/corrected/path_analysis/corrected_path_simulation_calibration_analysis18_v1.csv",
    "results/data/corrected/path_analysis/publication_release_gate_matrix_analysis18_v1.csv",
    "plans/publication-readiness-deep-audit/corrected_path_analysis_audit_analysis18_v1.md",
    "plans/publication-readiness-deep-audit/upstream_methods_primary_literature_benchmark.md",
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> None:
    if not NOTEBOOK.exists():
        raise FileNotFoundError(NOTEBOOK)
    missing = [relative for relative in REQUIRED_ARTIFACTS if not (ROOT / relative).exists()]
    if missing:
        raise FileNotFoundError("Notebook release inputs missing: " + ", ".join(missing))

    notebook = json.loads(NOTEBOOK.read_text())
    code_cells = [cell for cell in notebook["cells"] if cell["cell_type"] == "code"]
    unexecuted = [index for index, cell in enumerate(code_cells) if cell.get("execution_count") is None]
    errors = [
        {"code_cell_index": index, "ename": output.get("ename"), "evalue": output.get("evalue")}
        for index, cell in enumerate(code_cells)
        for output in cell.get("outputs", [])
        if output.get("output_type") == "error"
    ]
    code_source = "\n".join("".join(cell["source"]) for cell in code_cells)
    if unexecuted:
        raise RuntimeError(f"Notebook contains unexecuted code cells: {unexecuted}")
    if errors:
        raise RuntimeError(f"Notebook contains execution errors: {errors}")
    if "genome_size_pg" in code_source:
        raise RuntimeError("Notebook code reads or references the blocked historical picogram field")
    if "results/data/corrected" not in code_source:
        raise RuntimeError("Notebook does not consume the corrected release branch")

    manifest = {
        "analysis_id": "Desmognathus_publication_audit_analysis18_v1",
        "audited_at_utc": datetime.now(timezone.utc).isoformat(),
        "notebook": {
            "path": str(NOTEBOOK.relative_to(ROOT)),
            "sha256": sha256(NOTEBOOK),
            "n_cells": len(notebook["cells"]),
            "n_code_cells": len(code_cells),
            "all_code_cells_executed": True,
            "n_error_outputs": 0,
        },
        "required_artifacts": [
            {"path": relative, "sha256": sha256(ROOT / relative)}
            for relative in REQUIRED_ARTIFACTS
        ],
        "absolute_genome_size_used": False,
        "causal_claim_allowed": False,
        "historical_outputs_overwritten": False,
    }
    OUTPUT.write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"Audited executed notebook: {manifest['notebook']['sha256']}")


if __name__ == "__main__":
    main()
