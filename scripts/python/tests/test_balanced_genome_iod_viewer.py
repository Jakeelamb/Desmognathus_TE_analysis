#!/usr/bin/env python3
"""Tests for balanced genome-IOD viewer helpers."""

import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np


PROJECT_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(PROJECT_ROOT / "path_analysis" / "scripts"))

import build_balanced_genome_iod_viewer as viewer  # noqa: E402


class TestBalancedGenomeIodViewer(unittest.TestCase):
    def test_bbox_for_masks_spans_cell_and_nucleus(self):
        cell = np.zeros((20, 20), dtype=bool)
        nucleus = np.zeros((20, 20), dtype=bool)
        cell[4:8, 3:9] = True
        nucleus[11:14, 12:16] = True

        self.assertEqual(viewer.bbox_for_masks([cell, nucleus]), (4, 14, 3, 16))

    def test_square_crop_bounds_stays_inside_image(self):
        bounds = viewer.square_crop_bounds((0, 5, 0, 8), (30, 40), margin=10)

        top, bottom, left, right = bounds
        self.assertGreaterEqual(top, 0)
        self.assertGreaterEqual(left, 0)
        self.assertLessEqual(bottom, 30)
        self.assertLessEqual(right, 40)
        self.assertGreater(bottom, top)
        self.assertGreater(right, left)

    def test_species_slug_is_filesystem_safe(self):
        self.assertEqual(viewer.species_slug("D. gvnigeus/gwotli"), "gvnigeus_gwotli")

    def test_index_supports_species_query_links(self):
        records = [
            {
                "id": "fuscus-0001",
                "index": 1,
                "panel": "balanced_qc_only",
                "species": "D. fuscus",
                "species_slug": "fuscus",
                "rank": 1,
                "review_key": "rk1",
                "filename": "img_raw_green.ome.tiff",
                "tile_name": "tile_a",
                "cell_label": 1,
                "nucleus_label": 2,
                "cell_area_um2": 100.0,
                "nucleus_area_um2": 25.0,
                "nuc_iod": 700.0,
                "nuc_mean_od": 0.5,
                "qc_score": 0.8,
                "clarity": 0.7,
                "raw_src": "assets/fuscus/0001_raw.jpg",
                "cell_overlay_src": "assets/fuscus/0001_cell.png",
                "nucleus_overlay_src": "assets/fuscus/0001_nucleus.png",
                "missing_mask": False,
            }
        ]
        with tempfile.TemporaryDirectory() as tmpdir:
            index = viewer.write_index(Path(tmpdir), records, panel="balanced_qc_only")
            html = index.read_text(encoding="utf-8")

        self.assertIn("new URLSearchParams(window.location.search)", html)
        self.assertIn('params.get("species")', html)
        self.assertIn('params.set("species", state.species)', html)

    def test_index_accepts_a_queue_specific_title_and_storage_key(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            index = viewer.write_index(
                Path(tmpdir),
                [],
                panel="literal_largest_cell_mask_review",
                title="Largest Cell Mask Review",
                storage_key="largest-cell-mask-review-v1",
                default_decision="keep",
                blind_target_metrics=True,
                common_support_only=True,
                sort_by_match_distance=True,
            )
            html = index.read_text(encoding="utf-8")

        self.assertIn("Largest Cell Mask Review", html)
        self.assertIn('const STORAGE_KEY = "largest-cell-mask-review-v1";', html)
        self.assertIn('const EXPORT_FILENAME = "literal_largest_cell_mask_review_decisions.csv";', html)
        self.assertIn('const DEFAULT_DECISION = "keep";', html)
        self.assertIn('const EXPORT_ALL_RECORDS = Boolean(DEFAULT_DECISION);', html)
        self.assertIn('const BLIND_TARGET_METRICS = true;', html)
        self.assertIn('const COMMON_SUPPORT_ONLY = true;', html)
        self.assertIn('const REVIEW_RECORDS = COMMON_SUPPORT_ONLY', html)

    def test_index_supports_per_record_initial_decisions(self):
        records = [
            {
                "id": "old-0001",
                "species": "D. old",
                "review_key": "old-key",
                "filename": "old.tif",
                "tile_name": "tile",
                "initial_decision": "keep",
                "quality_match_status": "common_support",
            },
            {
                "id": "new-0002",
                "species": "D. new",
                "review_key": "new-key",
                "filename": "new.tif",
                "tile_name": "tile",
                "initial_decision": "",
                "quality_match_status": "common_support",
            },
        ]
        with tempfile.TemporaryDirectory() as tmpdir:
            index = viewer.write_index(Path(tmpdir), records, panel="extended_review")
            html = index.read_text(encoding="utf-8")

        self.assertIn('record.initial_decision || DEFAULT_DECISION || ""', html)
        self.assertIn('filter(record => resolvedDecision(record))', html)
        self.assertIn('state.decisions[record.id] ? "explicit" : "initial_default"', html)

    def test_index_can_auto_fill_a_ranked_species_target(self):
        records = [
            {
                "id": f"new-{rank}",
                "species": "D. new",
                "rank": rank,
                "review_key": f"new-key-{rank}",
                "filename": "new.tif",
                "tile_name": f"tile-{rank}",
                "initial_decision": "",
                "review_cohort": "new_species_unreviewed",
                "quality_match_status": "common_support",
            }
            for rank in range(1, 4)
        ]
        with tempfile.TemporaryDirectory() as tmpdir:
            index = viewer.write_index(
                Path(tmpdir),
                records,
                panel="auto_fill_review",
                auto_fill_target_per_species=2,
            )
            html = index.read_text(encoding="utf-8")

        self.assertIn("const AUTO_FILL_TARGET_PER_SPECIES = 2;", html)
        self.assertIn("function computeAutoSelectedIds()", html)
        self.assertIn('return autoSelectedIds.has(record.id) ? "keep" : "";', html)
        self.assertIn('return "auto_fill";', html)


if __name__ == "__main__":
    unittest.main()
