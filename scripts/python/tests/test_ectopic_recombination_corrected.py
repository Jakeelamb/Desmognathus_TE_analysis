#!/usr/bin/env python3
"""Tests for corrected ectopic-recombination depth and join semantics."""

import unittest

import pandas as pd

from scripts.processing.build_corrected_ectopic_recombination import (
    calculate_zero_aware_depth_metrics,
    merge_elements_with_tesorter,
)


class TestCorrectedEctopicRecombination(unittest.TestCase):
    def test_zero_depth_positions_remain_in_regional_means(self):
        depth = pd.DataFrame(
            {
                "position": list(range(1, 11)),
                "depth": [2, 0, 1, 1, 0, 0, 1, 1, 0, 2],
            }
        )

        metrics = calculate_zero_aware_depth_metrics(
            depth=depth,
            element_length=10,
            left_ltr_length=2,
            right_ltr_length=2,
        )

        self.assertAlmostEqual(metrics["mean_depth_terminal_all_positions"], 1.0)
        self.assertAlmostEqual(metrics["mean_depth_internal_all_positions"], 2 / 3)
        self.assertAlmostEqual(metrics["ratio_terminal_internal_all_positions"], 1.5)
        self.assertAlmostEqual(metrics["ratio_terminal_internal_nonzero_only"], 2.0)
        self.assertAlmostEqual(metrics["terminal_positive_coverage_fraction"], 0.5)
        self.assertAlmostEqual(metrics["internal_positive_coverage_fraction"], 4 / 6)

    def test_tesorter_join_uses_full_element_coordinates(self):
        elements = pd.DataFrame(
            {
                "sequence": ["contig", "contig"],
                "element start": [10, 100],
                "element end": [50, 150],
            }
        )
        tesorter = pd.DataFrame(
            {
                "#TE": ["contig_10_50", "contig_100_150"],
                "Domains": ["GAG|RT", "GAG|RT|RH"],
            }
        )

        merged = merge_elements_with_tesorter(elements, tesorter)

        self.assertEqual(len(merged), 2)
        self.assertEqual(merged.loc[0, "#TE"], "contig_10_50")
        self.assertEqual(merged.loc[1, "#TE"], "contig_100_150")
        self.assertEqual(merged.loc[0, "Domains"], "GAG|RT")
        self.assertEqual(merged.loc[1, "Domains"], "GAG|RT|RH")


if __name__ == "__main__":
    unittest.main()
