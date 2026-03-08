#!/usr/bin/env python3
"""
Unit tests for diversity metrics on supported TE proportion table layouts.
"""

import sys
import unittest
from pathlib import Path

import pandas as pd

parent_dir = str(Path(__file__).resolve().parent.parent)
sys.path.insert(0, parent_dir)

from analysis.run_diversity_analysis import calculate_diversity_metrics


class TestRunDiversityAnalysis(unittest.TestCase):
    """Verify diversity metrics for canonical and legacy table shapes."""

    def test_species_row_input(self):
        df = pd.DataFrame({
            "Species": ["D.aeneus", "D.ocoee"],
            "Gypsy": [0.5, 0.25],
            "LINE": [0.25, 0.25],
            "DNA": [0.25, 0.5],
        })
        result = calculate_diversity_metrics(df)
        self.assertListEqual(result["species"].tolist(), ["D.aeneus", "D.ocoee"])
        self.assertTrue((result["richness"] == 3).all())

    def test_feature_row_input(self):
        df = pd.DataFrame({
            "superfamily": ["Gypsy", "LINE", "DNA"],
            "D.aeneus": [0.5, 0.25, 0.25],
            "D.ocoee": [0.25, 0.25, 0.5],
        })
        result = calculate_diversity_metrics(df)
        self.assertListEqual(result["species"].tolist(), ["D.aeneus", "D.ocoee"])
        self.assertTrue((result["richness"] == 3).all())


if __name__ == "__main__":
    unittest.main()
