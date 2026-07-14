#!/usr/bin/env python3
"""Tests for balanced genome-IOD sensitivity helpers."""

import sys
import tempfile
import unittest
from pathlib import Path

import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(PROJECT_ROOT / "path_analysis" / "scripts"))

import build_balanced_genome_iod_sensitivity as balanced_iod  # noqa: E402


def row(species="D. testus", filename="img1", review_key="r1", **overrides):
    base = {
        "species": species,
        "filename": filename,
        "specimen_id": filename,
        "review_key": review_key,
        "has_cell_match": True,
        "physical_pair_ok": True,
        "one_to_one_cell": True,
        "keep_mask_pair": True,
        "keep_strict_core": True,
        "cell_edge_touch": False,
        "cell_area_um2": 100.0,
        "nuc_area_um2": 10.0,
        "nuc_iod": 100.0,
        "nuc_mean_od": 0.2,
        "nucleus_mask_overlap_label_count": 1,
        "nucleus_mask_best_overlap_fraction": 1.0,
        "distance_over_cell_radius": 0.1,
        "raw_focus_lap_var": 1.0,
        "raw_focus_grad_mean": 1.0,
        "raw_focus_intensity_iqr": 1.0,
        "cell_shape_smoothness": 1.0,
        "cell_shape_ellipse_iou": 1.0,
        "cell_shape_solidity": 1.0,
        "nucleus_shape_smoothness": 1.0,
        "nucleus_shape_ellipse_iou": 1.0,
        "nucleus_shape_solidity": 1.0,
        "shape_is_suspect": 0.0,
        "decision": "unlabeled",
    }
    base.update(overrides)
    return base


class TestBalancedGenomeIodSensitivity(unittest.TestCase):
    def test_qc_score_does_not_reward_darkness_or_large_cell_area(self):
        df = pd.DataFrame(
            [
                row(
                    review_key="dark_large_blurry",
                    cell_area_um2=500.0,
                    nuc_iod=1000.0,
                    nuc_mean_od=0.9,
                    raw_focus_lap_var=1.0,
                    raw_focus_grad_mean=1.0,
                    raw_focus_intensity_iqr=1.0,
                ),
                row(
                    review_key="clear_small_light",
                    cell_area_um2=100.0,
                    nuc_iod=100.0,
                    nuc_mean_od=0.1,
                    raw_focus_lap_var=10.0,
                    raw_focus_grad_mean=10.0,
                    raw_focus_intensity_iqr=10.0,
                ),
            ]
        )
        scored = balanced_iod.add_qc_scores(df)
        selected = balanced_iod.select_balanced_by_image(scored, target_per_species=1, panel="balanced_qc_only")

        self.assertEqual(selected.iloc[0]["review_key"], "clear_small_light")

    def test_balanced_selection_round_robins_across_images(self):
        rows = []
        for image in ["img1", "img2"]:
            for i in range(3):
                rows.append(
                    row(
                        filename=image,
                        specimen_id=image,
                        review_key=f"{image}_{i}",
                        raw_focus_lap_var=10 - i,
                        raw_focus_grad_mean=10 - i,
                        raw_focus_intensity_iqr=10 - i,
                    )
                )
        scored = balanced_iod.add_qc_scores(pd.DataFrame(rows))
        selected = balanced_iod.select_balanced_by_image(scored, target_per_species=4, panel="balanced_qc_only")

        self.assertEqual(selected["filename"].value_counts().to_dict(), {"img1": 2, "img2": 2})

    def test_panel_calibration_uses_reference_species_same_sampling_rule(self):
        df = pd.DataFrame(
            [
                row(species="D. fuscus", filename="f1", review_key="f1", nuc_iod=1000.0),
                row(species="D. fuscus", filename="f2", review_key="f2", nuc_iod=1000.0),
                row(species="D. other", filename="o1", review_key="o1", nuc_iod=500.0),
                row(species="D. other", filename="o2", review_key="o2", nuc_iod=500.0),
            ]
        )
        df["panel"] = "balanced_qc_only"
        df["specimen_group"] = df["filename"]
        df["image_group"] = df["filename"]
        summary = balanced_iod.summarize_panel(df, reference_species="D. fuscus", reference_genome_pg=16.36)
        other = summary.loc[summary["species"].eq("D. other")].iloc[0]

        self.assertAlmostEqual(other["estimated_genome_pg"], 8.18)
        self.assertAlmostEqual(other["iod_ratio_to_reference"], 0.5)

    def test_problem_viewer_decisions_are_loaded_and_normalized(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "decisions.csv"
            pd.DataFrame(
                [
                    {"species": "valtos", "review_key": "rk1", "decision": "problem"},
                    {"species": "D. valtos", "review_key": "rk2", "decision": "keep"},
                    {"species": "Desmognathus intermedius", "review_key": "rk3", "decision": "Problem"},
                ]
            ).to_csv(path, index=False)

            problems = balanced_iod.load_problem_review_keys(path)

        self.assertEqual(set(problems["review_key"]), {"rk1", "rk3"})
        self.assertEqual(set(problems["species"]), {"D. valtos", "D. intermedius"})

    def test_problem_review_keys_are_excluded_before_replacement_selection(self):
        candidates = pd.DataFrame(
            [
                row(species="D. valtos", filename="img1", review_key="bad", raw_focus_lap_var=10),
                row(species="D. valtos", filename="img2", review_key="replacement", raw_focus_lap_var=9),
                row(species="D. valtos", filename="img3", review_key="keep", raw_focus_lap_var=8),
            ]
        )
        scored = balanced_iod.add_qc_scores(candidates)
        problems = pd.DataFrame([{"species": "D. valtos", "review_key": "bad", "decision": "problem"}])
        filtered = balanced_iod.exclude_problem_review_keys(scored, problems)
        selected = balanced_iod.select_balanced_by_image(filtered, target_per_species=2, panel="balanced_qc_curated")

        self.assertNotIn("bad", set(selected["review_key"]))
        self.assertEqual(set(selected["review_key"]), {"replacement", "keep"})


if __name__ == "__main__":
    unittest.main()
