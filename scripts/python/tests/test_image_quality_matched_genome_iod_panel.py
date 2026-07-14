#!/usr/bin/env python3
"""Tests for cross-species technical image-quality matching."""

import sys
import tempfile
import unittest
from pathlib import Path

import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(PROJECT_ROOT / "path_analysis" / "scripts"))
import build_image_quality_matched_genome_iod_panel as matching  # noqa: E402


FEATURES = ["q1", "q2"]


def fixture() -> pd.DataFrame:
    rows = []
    for species, offset in [("D. alpha", 0.0), ("D. beta", 0.15), ("D. gamma", -0.10)]:
        for image in ["i1", "i2"]:
            for i in range(8):
                rows.append(
                    {
                        "species": species,
                        "filename": f"{species}-{image}",
                        "specimen_group": image,
                        "review_key": f"{species}-{image}-{i}",
                        "q1": i / 8 + offset,
                        "q2": (7 - i) / 8 - offset,
                        "nuc_iod": 1000 + i * 100,
                        "nuc_area_um2": 10 + i,
                    }
                )
    return pd.DataFrame(rows)


class TestImageQualityMatchedGenomeIodPanel(unittest.TestCase):
    def test_matches_equal_counts_and_represents_images(self):
        panel, balance, _details = matching.match_panel(fixture(), target_per_species=8, feature_columns=FEATURES)
        self.assertTrue(panel.groupby("species").size().eq(8).all())
        self.assertEqual(len(balance), 3)
        self.assertTrue(set(balance["quality_match_status"]).issubset({"common_support", "limited_overlap"}))

    def test_changing_iod_does_not_change_membership(self):
        original, _balance, _details = matching.match_panel(fixture(), target_per_species=8, feature_columns=FEATURES)
        changed = fixture()
        changed["nuc_iod"] = list(reversed(changed["nuc_iod"].tolist()))
        rematched, _balance, _details = matching.match_panel(changed, target_per_species=8, feature_columns=FEATURES)
        self.assertEqual(original["review_key"].tolist(), rematched["review_key"].tolist())

    def test_target_is_capped_by_smallest_species_pool(self):
        frame = fixture()
        frame = frame.loc[~frame["species"].eq("D. gamma") | frame.groupby("species").cumcount().lt(5)].copy()
        panel, _balance, details = matching.match_panel(frame, target_per_species=8, feature_columns=FEATURES)
        self.assertEqual(details["metadata"]["target_per_species"], 5)
        self.assertTrue(panel.groupby("species").size().eq(5).all())

    def test_review_registry_uses_latest_decision(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            first = Path(tmpdir) / "first.csv"
            second = Path(tmpdir) / "second.csv"
            pd.DataFrame([{"species": "D. alpha", "review_key": "a", "decision": "problem"}]).to_csv(first, index=False)
            pd.DataFrame([{"species": "D. alpha", "review_key": "a", "decision": "keep"}]).to_csv(second, index=False)
            registry = matching.load_review_registry([first, second])
        self.assertEqual(registry.iloc[0]["decision"], "keep")

    def test_reviewed_panel_freezes_with_all_common_species_above_minimum(self):
        balance = pd.DataFrame(
            {
                "species": ["D. alpha", "D. beta", "D. gamma"],
                "quality_match_status": ["common_support", "common_support", "limited_overlap"],
            }
        )
        frozen = pd.DataFrame({"species": ["D. alpha"] * 3 + ["D. beta"] * 3})
        self.assertTrue(
            matching.reviewed_panel_freeze_ready(frozen, balance, minimum_per_species=3)
        )

    def test_reviewed_panel_does_not_freeze_with_missing_or_thin_species(self):
        balance = pd.DataFrame(
            {
                "species": ["D. alpha", "D. beta"],
                "quality_match_status": ["common_support", "common_support"],
            }
        )
        missing = pd.DataFrame({"species": ["D. alpha"] * 3})
        thin = pd.DataFrame({"species": ["D. alpha"] * 3 + ["D. beta"] * 2})
        self.assertFalse(
            matching.reviewed_panel_freeze_ready(missing, balance, minimum_per_species=3)
        )
        self.assertFalse(
            matching.reviewed_panel_freeze_ready(thin, balance, minimum_per_species=3)
        )


if __name__ == "__main__":
    unittest.main()
