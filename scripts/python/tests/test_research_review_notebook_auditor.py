from __future__ import annotations

import importlib.util
import json
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
MODULE_PATH = ROOT / "scripts/processing/audit_research_review_notebooks.py"
SPEC = importlib.util.spec_from_file_location("audit_research_review_notebooks", MODULE_PATH)
assert SPEC is not None and SPEC.loader is not None
AUDIT = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(AUDIT)


def notebook(cells: list[dict[str, object]]) -> dict[str, object]:
    return {
        "cells": cells,
        "metadata": {"kernelspec": {"name": "python3"}},
        "nbformat": 4,
        "nbformat_minor": 5,
    }


def code_cell(
    source: str,
    *,
    execution_count: int | None = 1,
    outputs: list[dict[str, object]] | None = None,
) -> dict[str, object]:
    return {
        "cell_type": "code",
        "execution_count": execution_count,
        "metadata": {},
        "outputs": outputs or [],
        "source": source.splitlines(keepends=True),
    }


class ResearchReviewNotebookAuditorTests(unittest.TestCase):
    def test_frozen_product_names_are_not_mistaken_for_execution(self) -> None:
        source = "\n".join(
            [
                "repeatmasker = read_csv('corrected/repeatmasker_hits.csv')",
                "dna = read_csv('dnapipete_mass_accounting.csv')",
                "cellpose_qc = {'source': 'external_Cellpose_cpsam_pretraining'}",
            ]
        )
        self.assertEqual(AUDIT.detect_forbidden_execution(source), [])

    def test_process_and_cellpose_execution_surfaces_are_rejected(self) -> None:
        examples = [
            "!RepeatMasker assembly.fa",
            "import subprocess\nsubprocess.run(['dnaPipeTE', '--genome_size', '10'])",
            "import os\nos.system('RepeatModeler -database salamanders')",
            "from os import system as launch\nlaunch('RepeatMasker assembly.fa')",
            "from cellpose import models\nmodel = models.CellposeModel()\nmodel.eval(image)",
            "get_ipython().system('cellpose --train')",
        ]
        for source in examples:
            with self.subTest(source=source):
                self.assertTrue(AUDIT.detect_forbidden_execution(source))

    def test_markdown_is_not_scanned_but_code_must_be_executed_and_clean(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "good.ipynb"
            payload = notebook(
                [
                    {
                        "cell_type": "markdown",
                        "metadata": {},
                        "source": ["Methods used RepeatMasker, RepeatModeler, dnaPipeTE, and Cellpose."],
                    },
                    code_cell("value = 1"),
                ]
            )
            path.write_text(json.dumps(payload))
            record, failures = AUDIT.audit_notebook(path)
            self.assertEqual(failures, [])
            self.assertEqual(record["status"], "pass")

            payload["cells"][1] = code_cell("raise ValueError('bad')", execution_count=None)
            path.write_text(json.dumps(payload))
            record, failures = AUDIT.audit_notebook(path)
            self.assertEqual(record["status"], "fail")
            self.assertIn("all_code_cells_executed", {failure["gate"] for failure in failures})

    def test_error_outputs_fail_the_gate(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "error.ipynb"
            payload = notebook(
                [
                    code_cell(
                        "1 / 0",
                        outputs=[
                            {
                                "output_type": "error",
                                "ename": "ZeroDivisionError",
                                "evalue": "division by zero",
                                "traceback": [],
                            }
                        ],
                    )
                ]
            )
            path.write_text(json.dumps(payload))
            record, failures = AUDIT.audit_notebook(path)
            self.assertEqual(record["n_error_outputs"], 1)
            self.assertIn("no_error_outputs", {failure["gate"] for failure in failures})


if __name__ == "__main__":
    unittest.main()
