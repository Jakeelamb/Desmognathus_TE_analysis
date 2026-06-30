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
from build_master_dataset import require_input_files  # noqa: E402
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

    def test_setup_verifier_detects_paths_inside_conda_prefix(self):
        root = Path("/tmp/dusky-env")
        self.assertTrue(path_is_relative_to(root / "bin" / "python", root))
        self.assertFalse(path_is_relative_to(Path("/usr/bin/python"), root))


if __name__ == "__main__":
    unittest.main()
