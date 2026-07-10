from __future__ import annotations

import json
import sys
import unittest
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "path_analysis" / "scripts"))
import build_frozen_genome_iod_phylogeny_figure as figure  # noqa: E402


class GenomeIodPhylogenyFigureTests(unittest.TestCase):
    def test_paired_size_bootstrap_is_deterministic_and_preserves_draw_count(self) -> None:
        frame = pd.DataFrame(
            {
                "cell_area_um2": [10.0, 20.0, 30.0, 40.0],
                "nuc_area_um2": [1.0, 2.0, 3.0, 4.0],
            }
        )
        first = figure.bootstrap_paired_size_medians(frame, n_bootstrap=100, seed=7)
        second = figure.bootstrap_paired_size_medians(frame, n_bootstrap=100, seed=7)
        np.testing.assert_array_equal(first["cell_area_um2"], second["cell_area_um2"])
        np.testing.assert_array_equal(first["nuc_area_um2"], second["nuc_area_um2"])
        self.assertEqual(len(first), 100)

    def test_figure_data_uses_all_twenty_one_frozen_cell_species(self) -> None:
        summary, draws, correlations = figure.build_figure_data(
            n_bootstrap=100, seed=20260710
        )
        cell_species = set(pd.read_csv(figure.CELL_PATH, low_memory=False)["species"])
        self.assertEqual(set(summary["species"]), cell_species)
        self.assertEqual(len(summary), 21)
        self.assertIn("D. ochrophaeus", set(summary["species"]))
        ochrophaeus = summary.loc[summary["species"].eq("D. ochrophaeus")].iloc[0]
        self.assertTrue(pd.isna(ochrophaeus["relative_iod_index"]))
        self.assertEqual(ochrophaeus["genome_panel_status"], "limited_overlap_no_frozen_primary_iod")
        self.assertEqual(set(draws["metric"]), set(figure.METRIC_ORDER))
        draw_counts = draws.groupby("metric").size()
        self.assertEqual(draw_counts["relative_iod_index"], 20 * 100)
        self.assertEqual(draw_counts["nucleus_area_um2"], 21 * 100)
        self.assertEqual(draw_counts["cell_area_um2"], 21 * 100)
        self.assertEqual(set(correlations["comparison"]), {
            "relative_iod_vs_nucleus_area",
            "relative_iod_vs_cell_area",
            "nucleus_area_vs_cell_area",
        })

    def test_cell_and_nucleus_points_are_medians_of_frozen_top50(self) -> None:
        summary, _draws, _correlations = figure.build_figure_data(
            n_bootstrap=50, seed=12
        )
        cells = pd.read_csv(figure.CELL_PATH, low_memory=False)
        expected = cells.groupby("species").agg(
            expected_cell=("cell_area_um2", "median"),
            expected_nucleus=("nuc_area_um2", "median"),
        )
        observed = summary.set_index("species")
        pd.testing.assert_series_equal(
            observed.loc[expected.index, "cell_area_um2"],
            expected["expected_cell"],
            check_names=False,
        )
        pd.testing.assert_series_equal(
            observed.loc[expected.index, "nucleus_area_um2"],
            expected["expected_nucleus"],
            check_names=False,
        )

    def test_rendered_artifacts_and_manifest_are_presentation_ready(self) -> None:
        manifest = json.loads(figure.MANIFEST_PATH.read_text())
        notebook = json.loads(figure.EXECUTED_NOTEBOOK_PATH.read_text())
        source = "\n".join("".join(cell["source"]) for cell in notebook["cells"])

        self.assertTrue(figure.PNG_PATH.exists())
        self.assertGreater(figure.PNG_PATH.stat().st_size, 200_000)
        self.assertTrue(figure.PDF_PATH.exists())
        self.assertGreater(figure.PDF_PATH.stat().st_size, 20_000)
        self.assertEqual(manifest["n_measured_species"], 21)
        self.assertEqual(manifest["n_primary_genome_species"], 20)
        self.assertEqual(manifest["missing_primary_genome_species"], ["D. ochrophaeus"])
        self.assertEqual(manifest["phylogenetic_fills_used"], False)
        self.assertEqual(manifest["absolute_genome_size_claimed"], False)
        self.assertIn("Measured-only time-calibrated phylogeny", source)
        self.assertIn(figure.PNG_PATH.name, source)


if __name__ == "__main__":
    unittest.main()
