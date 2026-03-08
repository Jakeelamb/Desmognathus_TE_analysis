#!/usr/bin/env python3
"""
Unit tests for the superfamily proportions writer.
"""

import os
import sys
import tempfile
import unittest
from pathlib import Path

import pandas as pd

parent_dir = str(Path(__file__).resolve().parent.parent)
sys.path.insert(0, parent_dir)

from preprocessing.generate_superfamily_proportions import create_superfamily_proportions


class TestSuperfamilyProportions(unittest.TestCase):
    """Test superfamily proportions across supported input layouts."""

    def setUp(self):
        self.temp_dir = tempfile.mkdtemp()
        self.temp_input = os.path.join(self.temp_dir, "superfamily_breakdown.csv")
        self.temp_output = os.path.join(self.temp_dir, "superfamily_proportions.csv")

    def tearDown(self):
        import shutil

        shutil.rmtree(self.temp_dir)

    def test_create_superfamily_proportions_from_species_rows(self):
        df = pd.DataFrame({
            "Species": ["Desmognathus_aeneus", "Desmognathus_ocoee", "Desmognathus_orestes"],
            "DNA/hAT": [142, 157, 148],
            "LINE/BovB": [103, 92, 98],
            "LINE/L1": [254, 301, 287],
            "LTR/Gypsy": [187, 162, 176],
        })
        df.to_csv(self.temp_input, index=False)

        result = create_superfamily_proportions(self.temp_input, self.temp_output)
        self.assertTrue(result)

        result_df = pd.read_csv(self.temp_output, index_col=0)
        self.assertListEqual(
            result_df.columns.tolist(),
            ["DNA/hAT", "LINE/BovB", "LINE/L1", "LTR/Gypsy"],
        )
        self.assertListEqual(
            result_df.index.tolist(),
            ["Desmognathus_aeneus", "Desmognathus_ocoee", "Desmognathus_orestes"],
        )
        self.assertAlmostEqual(result_df.loc["Desmognathus_aeneus", "DNA/hAT"], 142 / 686, places=6)
        self.assertAlmostEqual(result_df.loc["Desmognathus_aeneus", "LINE/BovB"], 103 / 686, places=6)
        for _, row in result_df.iterrows():
            self.assertAlmostEqual(row.sum(), 1.0, places=6)

    def test_create_superfamily_proportions_from_feature_rows(self):
        df = pd.DataFrame({
            "superfamily": ["DNA/hAT", "LINE/BovB", "LINE/L1", "LTR/Gypsy"],
            "Desmognathus_aeneus": [142, 103, 254, 187],
            "Desmognathus_ocoee": [157, 92, 301, 162],
            "Desmognathus_orestes": [148, 98, 287, 176],
        })
        df.to_csv(self.temp_input, index=False)

        result = create_superfamily_proportions(self.temp_input, self.temp_output)
        self.assertTrue(result)

        result_df = pd.read_csv(self.temp_output, index_col=0)
        self.assertAlmostEqual(result_df.loc["Desmognathus_aeneus", "DNA/hAT"], 142 / 686, places=6)
        self.assertAlmostEqual(result_df.loc["Desmognathus_ocoee", "LINE/L1"], 301 / 712, places=6)
        for _, row in result_df.iterrows():
            self.assertAlmostEqual(row.sum(), 1.0, places=6)


if __name__ == "__main__":
    unittest.main()
