#!/usr/bin/env python3

from __future__ import annotations

import sys
import unittest
from pathlib import Path

import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(PROJECT_ROOT / "path_analysis/scripts"))

import extend_finalized_microscopy_reviewers as extension  # noqa: E402


class TestFinalizedMicroscopyReviewerExtension(unittest.TestCase):
    def test_morphology_queue_preserves_shortfall(self):
        rows = []
        for species, count in [("D. aeneus", 3), ("D. orestes", 5), ("D. wrighti", 1)]:
            rows.extend(
                {"species": species, "review_key": f"{species}-{idx}", "cell_area_um2": float(count - idx)}
                for idx in range(count)
            )
        queue, support = extension.build_new_morphology_queue(
            pd.DataFrame(rows), target_per_species=3, review_depth=4
        )

        self.assertEqual(queue.groupby("species").size().to_dict(), {
            "D. aeneus": 3,
            "D. orestes": 4,
            "D. wrighti": 1,
        })
        self.assertEqual(support["D. wrighti"]["shortfall_to_target"], 2)
        self.assertFalse(support["D. wrighti"]["target_available"])

    def test_frozen_scaling_uses_saved_clips_and_centers(self):
        frame = pd.DataFrame({
            "match_log_edge_sharpness": [-5.0, 2.0, 20.0],
            "match_log_relative_ring_noise": [0.0, 5.0, 10.0],
        })
        scaling = {
            "match_log_edge_sharpness": {"clip_low": 0.0, "clip_high": 10.0, "median": 2.0, "scale_iqr": 2.0},
            "match_log_relative_ring_noise": {"clip_low": 1.0, "clip_high": 9.0, "median": 5.0, "scale_iqr": 4.0},
        }

        out = extension.apply_frozen_scaling(frame, scaling)

        self.assertEqual(out["z_match_log_edge_sharpness"].tolist(), [-1.0, 0.0, 4.0])
        self.assertEqual(out["z_match_log_relative_ring_noise"].tolist(), [-1.0, 0.0, 1.0])

    def test_tile_global_centroids_are_localized_for_pixel_qc(self):
        frame = pd.DataFrame({
            "nucleus_source_image_path": [
                "/tmp/tiles/sample/tile_y008192_x016384.tiff",
                "/tmp/tiles/sample/tile_y000000_x000000.tiff",
            ],
            "nuc_centroid_x": [17_930.5, 512.0],
            "nuc_centroid_y": [9_216.25, 256.0],
            "tile_x0": [16_384, 0],
            "tile_y0": [8_192, 0],
        })

        out = extension.localize_tile_coordinates_for_pixel_qc(frame)

        self.assertEqual(out["nuc_centroid_x"].tolist(), [1546.5, 512.0])
        self.assertEqual(out["nuc_centroid_y"].tolist(), [1024.25, 256.0])
        self.assertEqual(out["tile_x0"].tolist(), [0, 0])
        self.assertEqual(out["tile_y0"].tolist(), [0, 0])
        self.assertEqual(out["source_global_nuc_centroid_x"].tolist(), [17_930.5, 512.0])

    def test_partial_template_match_returns_available_candidates_only(self):
        template = pd.DataFrame({
            "quality_template_id": [1, 2, 3],
            "z_match_log_edge_sharpness": [0.0, 1.0, 2.0],
            "z_match_log_relative_ring_noise": [0.0, 1.0, 2.0],
        })
        candidates = pd.DataFrame({
            "review_key": ["near-1", "near-3"],
            "z_match_log_edge_sharpness": [0.1, 1.9],
            "z_match_log_relative_ring_noise": [0.1, 1.9],
        })

        matched, metrics = extension.match_species_to_template(template, candidates)

        self.assertEqual(len(matched), 2)
        self.assertEqual(matched["quality_template_id"].tolist(), [1, 3])
        self.assertEqual(metrics["n_matched"], 2)

    def test_iod_review_composition_recovers_legacy_limited_overlap_species(self):
        frozen = pd.DataFrame({
            "species": ["D. fuscus"],
            "review_key": ["frozen-1"],
        })
        original_matched = pd.DataFrame({
            "species": ["D. fuscus", "D. ochrophaeus"],
            "review_key": ["frozen-1", "ochro-1"],
            "quality_match_status": ["common_support", "limited_overlap"],
        })
        new_matched = pd.DataFrame({
            "species": ["D. aeneus", "D. orestes", "D. wrighti"],
            "review_key": ["aeneus-1", "orestes-1", "wrighti-1"],
        })

        combined, recovered = extension.compose_iod_review_extension(
            frozen,
            original_matched,
            new_matched,
        )

        self.assertEqual(
            set(combined["species"]),
            {"D. fuscus", "D. ochrophaeus", "D. aeneus", "D. orestes", "D. wrighti"},
        )
        self.assertEqual(recovered["review_key"].tolist(), ["ochro-1"])
        self.assertEqual(recovered["viewer_default_decision"].tolist(), [""])
        self.assertEqual(recovered["review_status"].tolist(), ["unreviewed"])
        self.assertEqual(recovered["review_cohort"].tolist(), ["legacy_limited_overlap_unreviewed"])
        self.assertEqual(combined["review_key"].nunique(), len(combined))

    def test_existing_new_iod_matches_require_all_three_species(self):
        complete = pd.DataFrame({
            "species": ["D. aeneus", "D. orestes", "D. wrighti", "D. fuscus"],
            "review_key": ["a-1", "o-1", "w-1", "f-1"],
        })

        recovered = extension.extract_existing_new_iod_matches(complete)

        self.assertEqual(set(recovered["species"]), set(extension.NEW_SPECIES))
        with self.assertRaisesRegex(ValueError, "missing new species"):
            extension.extract_existing_new_iod_matches(
                complete.loc[~complete["species"].eq("D. wrighti")]
            )


if __name__ == "__main__":
    unittest.main()
