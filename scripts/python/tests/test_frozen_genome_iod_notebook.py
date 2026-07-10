from __future__ import annotations

import hashlib
import json
import sys
import unittest
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "path_analysis" / "scripts"))
import build_frozen_genome_iod_notebook as report  # noqa: E402


class FrozenGenomeIodNotebookTests(unittest.TestCase):
    def test_species_summary_weights_images_equally(self) -> None:
        frame = pd.DataFrame(
            {
                "species": ["D. alpha"] * 4 + ["D. beta"] * 2,
                "filename": ["a1", "a1", "a2", "a2", "b1", "b1"],
                "specimen_group": [1, 1, 2, 2, 3, 3],
                "review_key": [f"r{i}" for i in range(6)],
                "nuc_iod": [1.0, 3.0, 5.0, 7.0, 2.0, 4.0],
                "nuc_area_um2": [10.0, 12.0, 14.0, 16.0, 11.0, 13.0],
                "nuc_mean_od": [0.2, 0.2, 0.3, 0.3, 0.25, 0.25],
            }
        )
        species, images = report.summarize_species(frame, n_bootstrap=100, seed=11)
        indexed = species.set_index("species")

        self.assertEqual(len(images), 3)
        self.assertAlmostEqual(indexed.loc["D. alpha", "iod_equal_image_estimate"], 4.0)
        self.assertAlmostEqual(indexed.loc["D. beta", "iod_equal_image_estimate"], 3.0)
        self.assertAlmostEqual(indexed.loc["D. alpha", "relative_iod_index"], 4.0 / 3.5)
        self.assertEqual(indexed.loc["D. alpha", "n_images"], 2)

    def test_hierarchical_bootstrap_is_reproducible(self) -> None:
        frame = pd.DataFrame(
            {
                "filename": ["a", "a", "b", "b"],
                "nuc_iod": [10.0, 12.0, 20.0, 22.0],
            }
        )
        first = report.hierarchical_bootstrap_equal_image_estimate(
            frame, n_bootstrap=200, seed=42
        )
        second = report.hierarchical_bootstrap_equal_image_estimate(
            frame, n_bootstrap=200, seed=42
        )
        self.assertEqual(first, second)
        self.assertLessEqual(first[0], 16.0)
        self.assertGreaterEqual(first[1], 16.0)

    def test_rendered_notebook_and_report_are_complete(self) -> None:
        notebook = json.loads(report.EXECUTED_NOTEBOOK_PATH.read_text())
        code_cells = [cell for cell in notebook["cells"] if cell["cell_type"] == "code"]
        errors = [
            output
            for cell in code_cells
            for output in cell.get("outputs", [])
            if output.get("output_type") == "error"
        ]
        source = "\n".join("".join(cell["source"]) for cell in notebook["cells"])

        self.assertGreaterEqual(len(code_cells), 10)
        self.assertTrue(all(cell.get("execution_count") is not None for cell in code_cells))
        self.assertEqual(errors, [])
        self.assertIn("relative nuclear-IOD", source)
        self.assertIn("not an absolute genome-size estimate", source)
        self.assertNotIn("genome_size_pg", source)
        self.assertTrue(report.HTML_PATH.exists())
        self.assertGreater(report.HTML_PATH.stat().st_size, 500_000)

    def test_analysis_manifest_and_tables_match_frozen_source(self) -> None:
        manifest = json.loads(report.ANALYSIS_MANIFEST_PATH.read_text())
        frozen = pd.read_csv(report.FROZEN_PATH, low_memory=False)
        species = pd.read_csv(report.SPECIES_SUMMARY_PATH)
        source_sha = hashlib.sha256(report.FROZEN_PATH.read_bytes()).hexdigest()

        self.assertEqual(manifest["frozen_source_sha256"], source_sha)
        self.assertEqual(manifest["n_frozen_nuclei"], 721)
        self.assertEqual(manifest["n_primary_species"], 20)
        self.assertEqual(len(frozen), manifest["n_frozen_nuclei"])
        self.assertEqual(len(species), manifest["n_primary_species"])
        self.assertEqual(manifest["absolute_genome_size_claimed"], False)


if __name__ == "__main__":
    unittest.main()
