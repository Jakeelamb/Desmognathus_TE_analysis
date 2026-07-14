from __future__ import annotations

import json
import hashlib
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
NOTEBOOK = ROOT / "notebooks/Desmognathus_publication_audit_analysis18_v1.ipynb"
MANIFEST = ROOT / "notebooks/Desmognathus_publication_audit_analysis18_v1.manifest.json"


class PublicationAuditNotebookTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.notebook = json.loads(NOTEBOOK.read_text())

    def test_notebook_is_executed_without_error_outputs(self) -> None:
        code_cells = [cell for cell in self.notebook["cells"] if cell["cell_type"] == "code"]
        self.assertGreaterEqual(len(code_cells), 14)
        self.assertTrue(all(cell.get("execution_count") is not None for cell in code_cells))
        errors = [
            output
            for cell in code_cells
            for output in cell.get("outputs", [])
            if output.get("output_type") == "error"
        ]
        self.assertEqual(errors, [])

    def test_code_uses_corrected_inputs_and_never_reads_historical_picograms(self) -> None:
        code_source = "\n".join(
            "".join(cell["source"])
            for cell in self.notebook["cells"]
            if cell["cell_type"] == "code"
        )
        self.assertIn("results/data/corrected", code_source)
        self.assertIn("publication_release_gate_matrix_analysis18_v1.csv", code_source)
        self.assertIn("corrected_path_simulation_calibration_analysis18_v1.csv", code_source)
        self.assertNotIn("genome_size_pg", code_source)

    def test_notebook_contains_final_panel_and_weird_result_queue(self) -> None:
        all_source = "\n".join("".join(cell["source"]) for cell in self.notebook["cells"])
        for species in ["amphileucus", "fuscus", "gvnigeusgwotli", "valtos", "welteri"]:
            self.assertIn(species, all_source)
        self.assertIn("Automatically generated weird-result queue", all_source)
        self.assertIn("no terminal:internal assembly proxy", all_source)
        self.assertIn("anchor edge interval includes zero", all_source)

    def test_release_manifest_matches_executed_notebook(self) -> None:
        manifest = json.loads(MANIFEST.read_text())
        digest = hashlib.sha256(NOTEBOOK.read_bytes()).hexdigest()
        self.assertEqual(manifest["notebook"]["sha256"], digest)
        self.assertTrue(manifest["notebook"]["all_code_cells_executed"])
        self.assertEqual(manifest["notebook"]["n_error_outputs"], 0)
        self.assertFalse(manifest["absolute_genome_size_used"])
        self.assertFalse(manifest["causal_claim_allowed"])


if __name__ == "__main__":
    unittest.main()
