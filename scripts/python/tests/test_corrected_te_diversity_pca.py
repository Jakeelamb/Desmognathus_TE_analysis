#!/usr/bin/env python3
"""Tests for mass-aware TE diversity and compositional PCA helpers."""

import unittest

import numpy as np
import pandas as pd

from scripts.processing.build_corrected_te_diversity_pca import (
    build_mass_aware_composition,
    calculate_diversity_indices,
    replace_zeros_and_clr,
)


class TestCorrectedTEDiversityPCA(unittest.TestCase):
    def test_diversity_indices_have_explicit_standard_definitions(self):
        metrics = calculate_diversity_indices(np.array([0.5, 0.5]))

        self.assertAlmostEqual(metrics["shannon_entropy"], np.log(2))
        self.assertAlmostEqual(metrics["gini_simpson"], 0.5)
        self.assertAlmostEqual(metrics["simpson_dominance"], 0.5)
        self.assertAlmostEqual(metrics["hill_q1"], 2.0)
        self.assertAlmostEqual(metrics["hill_q2"], 2.0)
        self.assertAlmostEqual(metrics["pielou_evenness"], 1.0)
        self.assertEqual(metrics["observed_richness"], 2)

    def test_mass_aware_composition_preserves_unresolved_mass(self):
        breakdown = pd.DataFrame(
            {"LTR": [75.0], "LINE": [25.0]}, index=["fuscus"]
        )
        retained = pd.Series([0.8], index=["fuscus"])
        unresolved = pd.Series([0.2], index=["fuscus"])

        result = build_mass_aware_composition(breakdown, retained, unresolved)

        self.assertAlmostEqual(result.loc["fuscus", "LTR"], 0.6)
        self.assertAlmostEqual(result.loc["fuscus", "LINE"], 0.2)
        self.assertAlmostEqual(result.loc["fuscus", "Unresolved"], 0.2)
        self.assertAlmostEqual(result.loc["fuscus"].sum(), 1.0)

    def test_zero_replacement_recloses_and_clr_centers_rows(self):
        composition = pd.DataFrame(
            {"A": [0.5, 0.0], "B": [0.5, 0.75], "C": [0.0, 0.25]},
            index=["sp1", "sp2"],
        )

        replaced, clr, metadata = replace_zeros_and_clr(
            composition, method="half_global_min_positive"
        )

        self.assertTrue((replaced > 0).all().all())
        np.testing.assert_allclose(replaced.sum(axis=1), 1.0)
        np.testing.assert_allclose(clr.mean(axis=1), 0.0, atol=1e-12)
        self.assertGreater(metadata["n_zeros_replaced"], 0)


if __name__ == "__main__":
    unittest.main()
