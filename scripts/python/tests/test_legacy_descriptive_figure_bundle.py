from __future__ import annotations

import json
import unittest
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[3]
DATA = ROOT / "results/data/legacy_descriptive"
FIGURES = ROOT / "results/figures/legacy_descriptive"


class LegacyDescriptiveFigureBundleTests(unittest.TestCase):
    def test_te34_quality_pca_diversity_and_ectopic_products_exist(self) -> None:
        quality = pd.read_csv(DATA / "dnapipete_run_quality_te34_legacy_descriptive_v1.csv")
        species_quality = pd.read_csv(DATA / "dnapipete_species_quality_te34_legacy_descriptive_v1.csv")
        scores = pd.read_csv(DATA / "te_superfamily_pca_scores_te34_legacy_descriptive_v1.csv")
        clusters = pd.read_csv(DATA / "te_superfamily_cluster_metrics_te34_legacy_descriptive_v1.csv")
        self.assertEqual(species_quality.species.nunique(), 34)
        self.assertGreaterEqual(quality.run_label.nunique(), 35)
        self.assertEqual(scores.species.nunique(), 34)
        self.assertEqual(set(clusters.k), set(range(1, 11)))
        self.assertTrue((DATA / "ectopic_anova_te34_legacy_descriptive_v1.txt").exists())

    def test_manifest_and_key_slide_quality_figures_are_present(self) -> None:
        manifest = json.loads((DATA / "legacy_descriptive_figure_bundle_v1.manifest.json").read_text())
        self.assertEqual(manifest["n_te_species"], 34)
        self.assertFalse(manifest["expensive_upstream_tools_executed"])
        for stem in [
            "dnapipete_input_quality_te34_legacy_descriptive_v1",
            "te_mean_composition_te34_legacy_descriptive_v1",
            "te_superfamily_pca_scree_te34_legacy_descriptive_v1",
            "te_superfamily_pca_elbow_te34_legacy_descriptive_v1",
            "te_superfamily_pca_clusters_te34_legacy_descriptive_v1",
            "te_diversity_indices_te34_legacy_descriptive_v1",
            "ectopic_ratio_violin_te34_legacy_descriptive_v1",
        ]:
            self.assertGreater((FIGURES / f"{stem}.png").stat().st_size, 30_000)
            self.assertGreater((FIGURES / f"{stem}.pdf").stat().st_size, 5_000)


if __name__ == "__main__":
    unittest.main()
