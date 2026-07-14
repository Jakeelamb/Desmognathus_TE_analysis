#!/usr/bin/env python3
"""Tests for configured-denominator dnaPipeTE repeat-load reconstruction."""

import unittest

import pandas as pd

from scripts.processing.build_dnapipete_absolute_load_sensitivity import (
    build_absolute_load,
)


class TestDnaPipeTEAbsoluteLoadSensitivity(unittest.TestCase):
    def test_class_components_conserve_repeat_aligned_mass(self):
        mass = pd.DataFrame(
            {
                "species": ["fuscus"],
                "te_sra_accession": ["SRX20497025"],
                "te_assembly_accession": ["GCA_032353935.1"],
                "total_aligned_bases": [1000],
                "class_retained_aligned_bases": [900],
                "class_unresolved_aligned_bases": [100],
            }
        )
        class_breakdown = pd.DataFrame(
            {
                "DNAtransposons Subclass1": [20.0],
                "Retrotransposons Autonomous": [60.0],
                "Other": [10.0],
                "Unknown": [10.0],
            },
            index=["D.fuscus"],
        )

        result = build_absolute_load(mass, class_breakdown, 20_000, 0.1)
        row = result.iloc[0]

        self.assertEqual(row["repeat_aligned_bases"], 1000)
        self.assertAlmostEqual(row["repeat_aligned_fraction"], 0.5)
        self.assertAlmostEqual(row["te_classified_fraction"], 0.36)
        self.assertAlmostEqual(row["other_repeat_fraction"], 0.045)
        self.assertAlmostEqual(row["class_unknown_fraction"], 0.045)
        self.assertAlmostEqual(row["class_unresolved_fraction"], 0.05)
        self.assertAlmostEqual(row["class_component_sum_fraction"], 0.5)

    def test_breakdown_must_sum_to_one_hundred(self):
        mass = pd.DataFrame(
            {
                "species": ["fuscus"],
                "te_sra_accession": ["SRX20497025"],
                "te_assembly_accession": ["GCA_032353935.1"],
                "total_aligned_bases": [1000],
                "class_retained_aligned_bases": [900],
                "class_unresolved_aligned_bases": [100],
            }
        )
        bad_breakdown = pd.DataFrame(
            {"Retrotransposons Autonomous": [80.0]}, index=["D.fuscus"]
        )

        with self.assertRaisesRegex(ValueError, "sum to 100"):
            build_absolute_load(mass, bad_breakdown, 20_000, 0.1)


if __name__ == "__main__":
    unittest.main()
