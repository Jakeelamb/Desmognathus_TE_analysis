#!/usr/bin/env python3
"""Tests for fail-closed phylogenetic path-model release gates."""

import unittest

import pandas as pd

from scripts.processing.audit_path_model_release_gates import audit_ranking


def ranking(rows):
    return pd.DataFrame(
        rows,
        columns=["model", "k", "q", "C", "p", "CICc", "delta_CICc", "l", "w"],
    )


class TestPathModelReleaseGates(unittest.TestCase):
    def test_supported_unique_top_model_passes(self):
        result = audit_ranking(
            ranking(
                [
                    ["top", 2, 3, 4, 0.20, 10, 0, 1, 0.90],
                    ["other", 2, 3, 5, 0.10, 15, 5, 0.08, 0.10],
                ]
            ),
            "family",
            12,
        )
        self.assertEqual(result["release_status"], "supported_unique_top_model")
        self.assertTrue(result["winner_claim_allowed_by_ranking_gate"])
        self.assertFalse(result["publication_claim_allowed"])

    def test_significant_top_model_blocks_winner(self):
        result = audit_ranking(
            ranking(
                [
                    ["rejected", 2, 3, 10, 0.01, 10, 0, 1, 0.80],
                    ["supported", 2, 3, 3, 0.20, 14, 4, 0.14, 0.20],
                ]
            ),
            "family",
            12,
        )
        self.assertEqual(result["release_status"], "blocked_top_model_rejected")
        self.assertEqual(result["best_supported_model"], "supported")
        self.assertFalse(result["winner_claim_allowed_by_ranking_gate"])

    def test_nonfinite_ranking_blocks_release(self):
        result = audit_ranking(
            ranking(
                [
                    ["top", 2, 3, 4, 0.20, 10, 0, 1, float("nan")],
                    ["other", 2, 3, 5, 0.10, float("inf"), float("nan"), 0, float("nan")],
                ]
            ),
            "family",
            8,
        )
        self.assertEqual(result["release_status"], "blocked_nonfinite_ranking")
        self.assertFalse(result["winner_claim_allowed_by_ranking_gate"])

    def test_delta_under_two_requires_competitive_set_language(self):
        result = audit_ranking(
            ranking(
                [
                    ["top", 2, 3, 4, 0.20, 10, 0, 1, 0.60],
                    ["competitive", 2, 3, 5, 0.10, 11, 1, 0.61, 0.40],
                ]
            ),
            "family",
            12,
        )
        self.assertEqual(result["release_status"], "supported_competitive_model_set")
        self.assertFalse(result["winner_claim_allowed_by_ranking_gate"])
        self.assertTrue(result["model_averaging_allowed_by_ranking_gate"])


if __name__ == "__main__":
    unittest.main()
