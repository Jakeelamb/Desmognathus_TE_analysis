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

    def test_figure_data_uses_all_finalized_species_and_flags_support(self) -> None:
        summary, draws, correlations = figure.build_figure_data(
            n_bootstrap=100, seed=20260710
        )
        cell_species = set(pd.read_csv(figure.CELL_PATH, low_memory=False)["species"])
        self.assertEqual(set(summary["species"]), cell_species)
        self.assertEqual(len(summary), 24)
        self.assertIn("D. ochrophaeus", set(summary["species"]))
        ochrophaeus = summary.loc[summary["species"].eq("D. ochrophaeus")].iloc[0]
        self.assertTrue(pd.notna(ochrophaeus["genome_size_pg_fuscus_anchored"]))
        self.assertEqual(
            ochrophaeus["genome_panel_status"],
            "finalized_included_limited_overlap",
        )
        self.assertEqual(
            summary.set_index("species").loc["D. aeneus", "genome_panel_status"],
            "finalized_included_limited_overlap",
        )
        self.assertEqual(
            summary.set_index("species").loc["D. wrighti", "genome_panel_status"],
            "finalized_included_limited_overlap",
        )
        self.assertEqual(
            summary.set_index("species").loc["D. orestes", "genome_panel_status"],
            "finalized_included_common_support",
        )
        self.assertEqual(set(draws["metric"]), set(figure.METRIC_ORDER))
        draw_counts = draws.groupby("metric").size()
        self.assertEqual(draw_counts["genome_size_pg_fuscus_anchored"], 24 * 100)
        self.assertEqual(draw_counts["nucleus_area_um2"], 24 * 100)
        self.assertEqual(draw_counts["cell_area_um2"], 24 * 100)
        self.assertEqual(set(correlations["comparison"]), {
            "genome_size_vs_nucleus_area",
            "genome_size_vs_cell_area",
            "nucleus_area_vs_cell_area",
        })
        self.assertTrue(correlations["phylogenetically_corrected"].all())
        self.assertEqual(
            correlations.set_index("comparison")["n_species"].to_dict(),
            {
                "genome_size_vs_nucleus_area": 24,
                "genome_size_vs_cell_area": 24,
                "nucleus_area_vs_cell_area": 24,
            },
        )
        support = summary.set_index("species")
        self.assertEqual(int(support.loc["D. aeneus", "n_size_cells"]), 44)
        self.assertEqual(int(support.loc["D. orestes", "n_size_cells"]), 50)
        self.assertEqual(int(support.loc["D. wrighti", "n_size_cells"]), 8)
        fuscus = summary.loc[summary["species"].eq("D. fuscus")].iloc[0]
        self.assertAlmostEqual(
            fuscus["genome_size_pg_fuscus_anchored"],
            figure.genome_analysis.REFERENCE_GENOME_SIZE_PG,
        )
        fuscus_draws = draws.loc[
            draws["species"].eq("D. fuscus")
            & draws["metric"].eq("genome_size_pg_fuscus_anchored"),
            "value",
        ]
        np.testing.assert_allclose(
            fuscus_draws,
            figure.genome_analysis.REFERENCE_GENOME_SIZE_PG,
        )

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

    def test_pagel_lambda_pgls_matches_established_path_fit(self) -> None:
        summary, _draws, correlations = figure.build_figure_data(
            n_bootstrap=20, seed=22
        )
        observed = correlations.set_index("comparison")
        self.assertAlmostEqual(
            observed.loc[
                "genome_size_vs_nucleus_area", "pgls_standardized_beta"
            ],
            0.41901841,
            places=8,
        )
        self.assertAlmostEqual(
            observed.loc[
                "genome_size_vs_cell_area", "pgls_standardized_beta"
            ],
            0.71666442,
            places=8,
        )
        self.assertAlmostEqual(
            observed.loc[
                "nucleus_area_vs_cell_area", "pgls_standardized_beta"
            ],
            0.5793012834,
            places=8,
        )
        self.assertTrue((observed["pagel_lambda"] < 1e-6).all())
        self.assertEqual(len(summary), 24)

    def test_rendered_artifacts_and_manifest_are_presentation_ready(self) -> None:
        manifest = json.loads(figure.MANIFEST_PATH.read_text())
        notebook = json.loads(figure.EXECUTED_NOTEBOOK_PATH.read_text())
        source = "\n".join("".join(cell["source"]) for cell in notebook["cells"])
        figure_source = Path(figure.__file__).read_text()

        self.assertTrue(figure.PNG_PATH.exists())
        self.assertGreater(figure.PNG_PATH.stat().st_size, 200_000)
        self.assertTrue(figure.PDF_PATH.exists())
        self.assertGreater(figure.PDF_PATH.stat().st_size, 20_000)
        self.assertTrue(figure.PAIRWISE_PNG_PATH.exists())
        self.assertGreater(figure.PAIRWISE_PNG_PATH.stat().st_size, 150_000)
        self.assertTrue(figure.PAIRWISE_PDF_PATH.exists())
        self.assertGreater(figure.PAIRWISE_PDF_PATH.stat().st_size, 20_000)
        self.assertEqual(manifest["n_measured_species"], 24)
        self.assertEqual(manifest["n_primary_genome_species"], 24)
        self.assertEqual(manifest["n_genome_estimate_species"], 24)
        self.assertEqual(manifest["missing_primary_genome_species"], [])
        self.assertEqual(
            manifest["limited_overlap_included_species"],
            ["D. aeneus", "D. ochrophaeus", "D. wrighti"],
        )
        self.assertEqual(manifest["all_finalized_species_included_in_genome_models"], True)
        self.assertEqual(
            manifest["reduced_size_sample_species"],
            {"D. aeneus": 44, "D. wrighti": 8},
        )
        self.assertEqual(manifest["phylogenetic_fills_used"], False)
        self.assertEqual(manifest["absolute_genome_size_claimed"], False)
        self.assertEqual(manifest["conditional_genome_size_estimates_reported"], True)
        self.assertEqual(manifest["correlations_are_phylogenetically_corrected"], True)
        self.assertIn("Pagel-lambda", manifest["pairwise_phylogenetic_method"])
        self.assertEqual(manifest["genome_reference_species"], "D. fuscus")
        self.assertEqual(manifest["genome_reference_pg"], 16.36)
        self.assertIn("Measured-only time-calibrated phylogeny", source)
        self.assertIn("Fuscus-anchored genome-size estimates", source)
        self.assertIn(figure.PNG_PATH.name, source)
        self.assertIn(figure.PAIRWISE_PNG_PATH.name, source)
        self.assertNotIn("D. ochrophaeus has no finalized IOD estimate", figure_source)

        correlations = pd.read_csv(figure.CORRELATION_PATH)
        self.assertTrue(
            {
                "spearman_rho",
                "pearson_r",
                "pgls_standardized_beta",
                "pgls_standard_error",
                "pgls_ci_low",
                "pgls_ci_high",
                "pgls_p_value",
                "pagel_lambda",
                "n_species",
            }.issubset(correlations.columns)
        )
        self.assertEqual(
            set(correlations["comparison"]),
            {
                "genome_size_vs_nucleus_area",
                "genome_size_vs_cell_area",
                "nucleus_area_vs_cell_area",
            },
        )


if __name__ == "__main__":
    unittest.main()
