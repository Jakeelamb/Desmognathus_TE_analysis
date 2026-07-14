from __future__ import annotations

import json
import unittest
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[3]
PANEL = ROOT / "path_analysis/data/derived/panels/study_cell_linked_panel21_v1.csv"
DATA = ROOT / "results/data/corrected/cell21"


class Cell21DescriptiveBundleTests(unittest.TestCase):
    def test_cell21_output_exactly_matches_linked_cell_panel(self) -> None:
        expected = set(pd.read_csv(PANEL)["species"])
        observed = pd.read_csv(DATA / "cell_linked_traits_cell21_v1.csv")
        self.assertEqual(set(observed["species"]), expected)
        self.assertEqual(len(observed), 21)
        self.assertTrue(observed[["cell_area_um2", "nucleus_area_um2", "relative_nuclear_iod_raw_value"]].notna().all().all())

    def test_manifest_keeps_relative_iod_out_of_absolute_genome_size_claims(self) -> None:
        manifest = json.loads((DATA / "cell21_descriptive_bundle_v1.manifest.json").read_text())
        self.assertEqual(manifest["n_linked_cell_species"], 21)
        self.assertFalse(manifest["relative_nuclear_iod_is_absolute_genome_size"])
        self.assertTrue(manifest["analysis18_outputs_preserved"])
        self.assertFalse(manifest["eligible_for_integrated_phylogenetic_path_analysis"])


if __name__ == "__main__":
    unittest.main()
