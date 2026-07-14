#!/usr/bin/env python3
"""Tests for non-destructive corrected RepeatMasker summaries."""

import unittest

import pandas as pd

from scripts.processing.build_corrected_divergence_summary import build_summary


class TestCorrectedDivergenceSummary(unittest.TestCase):
    def test_threshold_summary_uses_hit_class_and_preserves_zero_count_rows(self):
        hits = pd.DataFrame(
            {
                "Class": ["Retro", "Retro", "Retro"],
                "Order": ["LTR", "LTR", "LTR"],
                "Superfamily": ["Gypsy", "Gypsy", "Gypsy"],
                "Desmognathus_Species": ["D.fuscus", "D.fuscus", "D.welteri"],
                "score": [100.0, 200.0, 50.0],
                "percent_divergence": [10.0, 20.0, 30.0],
                "percent_deletions": [1.0, 3.0, 5.0],
                "percent_insertions": [2.0, 4.0, 6.0],
                "hitlength_contiglength": [0.95, 0.80, 0.50],
            }
        )

        summary = build_summary(hits, thresholds=[0.0, 0.9])
        order = summary.loc[
            summary["group_level"].eq("order") & summary["group_name"].eq("LTR")
        ]

        all_hits = order.loc[
            order["Desmognathus_Species"].eq("D.welteri")
            & order["threshold"].eq(-1.0)
        ].iloc[0]
        self.assertEqual(all_hits["count"], 1)
        self.assertEqual(all_hits["threshold_basis"], "none_all_hits")

        fuscus_zero = order.loc[
            order["Desmognathus_Species"].eq("D.fuscus")
            & order["threshold"].eq(0.0)
        ].iloc[0]
        self.assertEqual(fuscus_zero["count"], 2)
        self.assertEqual(fuscus_zero["percent_divergence_median"], 15.0)

        fuscus_high = order.loc[
            order["Desmognathus_Species"].eq("D.fuscus")
            & order["threshold"].eq(0.9)
        ].iloc[0]
        self.assertEqual(fuscus_high["count"], 1)
        self.assertEqual(fuscus_high["percent_deletions_median"], 1.0)

        welteri_high = order.loc[
            order["Desmognathus_Species"].eq("D.welteri")
            & order["threshold"].eq(0.9)
        ].iloc[0]
        self.assertEqual(welteri_high["count"], 0)
        self.assertTrue(pd.isna(welteri_high["percent_divergence_median"]))
        self.assertEqual(
            welteri_high["threshold_basis"],
            "dnapipete_contig_hitlength_ratio",
        )


if __name__ == "__main__":
    unittest.main()
