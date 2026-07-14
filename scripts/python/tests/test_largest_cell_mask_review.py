#!/usr/bin/env python3
"""Tests for literal-largest cell mask review and freeze gates."""

from __future__ import annotations

import sys
import unittest
from pathlib import Path

import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(PROJECT_ROOT / "path_analysis" / "scripts"))
import build_largest_cell_mask_review as review  # noqa: E402


def candidate(key: str, area: float) -> dict[str, object]:
    return {"species": "D. testus", "review_key": key, "cell_area_um2": area, "nuc_area_um2": 10.0, "nuc_iod": 100.0}


class TestLargestCellMaskReview(unittest.TestCase):
    def test_queue_ranks_by_cell_area(self):
        queue = review.build_review_queue(pd.DataFrame([candidate("a", 30), candidate("b", 50), candidate("c", 40)]), target_per_species=2, review_depth=3)
        self.assertEqual(queue["review_key"].tolist(), ["b", "c", "a"])

    def test_freeze_requires_higher_ranks_to_be_reviewed(self):
        queue = review.build_review_queue(pd.DataFrame([candidate("a", 40), candidate("b", 30), candidate("c", 20)]), target_per_species=2, review_depth=3)
        decisions = pd.DataFrame([{"species": "D. testus", "review_key": "b", "decision": "keep"}, {"species": "D. testus", "review_key": "c", "decision": "keep"}])
        with self.assertRaisesRegex(ValueError, "not fully reviewed"):
            review.freeze_queue(queue, decisions, target_per_species=2)

    def test_freeze_selects_first_keeps_after_explicit_rejection(self):
        queue = review.build_review_queue(pd.DataFrame([candidate("a", 40), candidate("b", 30), candidate("c", 20)]), target_per_species=2, review_depth=3)
        decisions = pd.DataFrame([{"species": "D. testus", "review_key": "a", "decision": "problem"}, {"species": "D. testus", "review_key": "b", "decision": "keep"}, {"species": "D. testus", "review_key": "c", "decision": "keep"}])
        frozen, audit = review.freeze_queue(queue, decisions, target_per_species=2)
        self.assertEqual(frozen["review_key"].tolist(), ["b", "c"])
        self.assertEqual(int(audit["frozen_selected"].sum()), 2)

    def test_freeze_uses_human_decisions_when_queue_has_historical_decision(self):
        queue = review.build_review_queue(
            pd.DataFrame([candidate("a", 40), candidate("b", 30), candidate("c", 20)]),
            target_per_species=2,
            review_depth=3,
        )
        queue["decision"] = "historical_model_value"
        decisions = pd.DataFrame(
            [
                {"species": "D. testus", "review_key": "a", "decision": "problem"},
                {"species": "D. testus", "review_key": "b", "decision": "keep"},
                {"species": "D. testus", "review_key": "c", "decision": "keep"},
            ]
        )
        frozen, _audit = review.freeze_queue(queue, decisions, target_per_species=2)
        self.assertEqual(frozen["review_key"].tolist(), ["b", "c"])


if __name__ == "__main__":
    unittest.main()
