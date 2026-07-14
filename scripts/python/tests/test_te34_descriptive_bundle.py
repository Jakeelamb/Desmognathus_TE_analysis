from __future__ import annotations

import json
import unittest
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[3]
DATA = ROOT / "results/data/corrected/te34"
PANEL = ROOT / "path_analysis/data/derived/panels/study_te_resource_panel34_v1.csv"


class Te34DescriptiveBundleTests(unittest.TestCase):
    def test_all_frozen_te_products_cover_the_declared_te_resource_panel(self) -> None:
        expected = set(pd.read_csv(PANEL)["species"])
        self.assertEqual(len(expected), 34)
        for filename in [
            "dnapipete_mass_accounting_te34_v1.csv",
            "te_diversity_mass_sensitivity_te34_v1.csv",
            "te_pca_scores_te34_v1.csv",
            "repeatmasker_hit_inventory_te34_v1.csv",
            "repeatmasker_divergence_landscape_te34_v1.csv",
        ]:
            observed = set(pd.read_csv(DATA / filename)["species"])
            self.assertEqual(observed, expected, filename)

    def test_repeatmasker_landscape_conserves_the_all_resource_hit_table(self) -> None:
        inventory = pd.read_csv(DATA / "repeatmasker_hit_inventory_te34_v1.csv")
        landscape = pd.read_csv(DATA / "repeatmasker_divergence_landscape_te34_v1.csv")
        self.assertEqual(int(inventory["hit_count"].sum()), int(landscape["hit_count"].sum()))
        self.assertEqual(int(inventory["aligned_hit_bp"].sum()), int(landscape["hit_bp"].sum()))
        fractions = landscape.groupby("species")["fraction_species_hit_bp"].sum()
        self.assertTrue(np.allclose(fractions, 1.0, atol=1e-12))

    def test_manifest_declares_descriptive_only_scope_and_preserves_path18(self) -> None:
        manifest = json.loads((DATA / "te34_descriptive_bundle_v1.manifest.json").read_text())
        self.assertEqual(manifest["n_te_species"], 34)
        self.assertEqual(manifest["analysis_scope"], "active_vetted_te_resource_panel34_descriptive")
        self.assertFalse(manifest["expensive_upstream_tools_executed"])
        self.assertTrue(manifest["analysis18_outputs_preserved"])
        self.assertFalse(manifest["eligible_for_integrated_phylogenetic_path_analysis"])


if __name__ == "__main__":
    unittest.main()
