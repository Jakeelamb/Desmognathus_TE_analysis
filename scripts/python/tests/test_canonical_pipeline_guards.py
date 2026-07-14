#!/usr/bin/env python3
"""Focused guards for the canonical repo-level analysis helpers."""

import sys
import tempfile
import unittest
from pathlib import Path

import pandas as pd
import yaml


PROJECT_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(PROJECT_ROOT))
sys.path.insert(0, str(PROJECT_ROOT / "path_analysis" / "scripts"))

from scripts.processing.diversity_stats import calculate_simpson_diversity  # noqa: E402
from build_canonical_diversity_tables import (  # noqa: E402
    calculate_simpson_diversity as calculate_audit_simpson_diversity,
)
from build_master_dataset import load_genome_results, load_morphology_results, require_input_files  # noqa: E402
from verify_setup import path_is_relative_to  # noqa: E402


class TestCanonicalPipelineGuards(unittest.TestCase):
    def test_simpson_diversity_uses_proportions_not_count_combinatorics(self):
        row = pd.Series([50.0, 25.0, 25.0])
        self.assertAlmostEqual(calculate_simpson_diversity(row), 0.625)

    def test_simpson_diversity_is_scale_invariant(self):
        percent_row = pd.Series([50.0, 25.0, 25.0])
        proportion_row = pd.Series([0.5, 0.25, 0.25])
        self.assertAlmostEqual(
            calculate_simpson_diversity(percent_row),
            calculate_simpson_diversity(proportion_row),
        )

    def test_diversity_audit_uses_the_same_gini_simpson_definition(self):
        percent_row = pd.Series([50.0, 25.0, 25.0])
        proportion_row = pd.Series([0.5, 0.25, 0.25])
        self.assertAlmostEqual(calculate_audit_simpson_diversity(percent_row), 0.625)
        self.assertAlmostEqual(
            calculate_audit_simpson_diversity(percent_row),
            calculate_audit_simpson_diversity(proportion_row),
        )

    def test_path_config_has_canonical_schema(self):
        with (PROJECT_ROOT / "paths.yaml").open("r", encoding="utf-8") as handle:
            canonical = yaml.safe_load(handle)

        self.assertIn("genomes", canonical["input_data"])
        self.assertIn("reports", canonical["results"])

    def test_master_dataset_preflight_reports_missing_generated_inputs(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            missing_path = Path(tmpdir) / "results" / "data" / "dnaPipeTE_order_breakdown.csv"
            with self.assertRaises(FileNotFoundError) as context:
                require_input_files({"dnaPipeTE order breakdown": missing_path})

        message = str(context.exception)
        self.assertIn("Missing generated analysis inputs", message)
        self.assertIn("dnaPipeTE order breakdown", message)
        self.assertIn("scripts/processing/dnaPipe.py", message)

    def test_genome_loader_preserves_verified_uncertainty_columns(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "genome.csv"
            pd.DataFrame(
                [
                    {
                        "species": "D. testus",
                        "primary_genome_pg": 18.5,
                        "primary_genome_ci_low_pg": 17.2,
                        "primary_genome_ci_high_pg": 19.8,
                        "primary_support_tier": "high",
                        "primary_support_warnings": "auto_dominant",
                        "genome_primary_missing": False,
                        "result_status": "minor_caution",
                    }
                ]
            ).to_csv(path, index=False)

            out = load_genome_results(path)

        self.assertEqual(out.loc[0, "species"], "testus")
        self.assertEqual(out.loc[0, "genome_size_pg"], 18.5)
        self.assertEqual(out.loc[0, "genome_size_ci_low_pg"], 17.2)
        self.assertEqual(out.loc[0, "genome_support_tier"], "high")

    def test_morphology_loader_preserves_verified_cell_and_nucleus_intervals(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "morphology.csv"
            pd.DataFrame(
                [
                    {
                        "species": "D. testus",
                        "species_median_cell_area_um2": 250.0,
                        "cell_area_um2_ci_low": 240.0,
                        "cell_area_um2_ci_high": 260.0,
                        "species_median_nuc_area_um2": 42.0,
                        "nucleus_area_um2_ci_low": 40.0,
                        "nucleus_area_um2_ci_high": 44.0,
                        "linked_support_label": "high",
                    }
                ]
            ).to_csv(path, index=False)

            out = load_morphology_results(path)

        self.assertEqual(out.loc[0, "morph_cell_area_um2"], 250.0)
        self.assertEqual(out.loc[0, "morph_cell_area_ci_low_um2"], 240.0)
        self.assertEqual(out.loc[0, "morph_nucleus_area_ci_high_um2"], 44.0)
        self.assertEqual(out.loc[0, "morph_support_label"], "high")

    def test_setup_verifier_detects_paths_inside_conda_prefix(self):
        root = Path("/tmp/dusky-env")
        self.assertTrue(path_is_relative_to(root / "bin" / "python", root))
        self.assertFalse(path_is_relative_to(Path("/usr/bin/python"), root))


if __name__ == "__main__":
    unittest.main()
