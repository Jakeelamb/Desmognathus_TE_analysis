#!/usr/bin/env python3
"""Tests for the RepeatMasker hit-classification audit metrics."""

import unittest
from pathlib import Path
import subprocess
import sys

import pandas as pd

from scripts.processing.audit_repeatmasker_hit_classification import summarize_chunk


PROJECT_ROOT = Path(__file__).resolve().parents[3]


class TestRepeatMaskerClassificationAudit(unittest.TestCase):
    def test_direct_entrypoint_resolves_repository_imports(self):
        script = (
            PROJECT_ROOT
            / "scripts"
            / "processing"
            / "audit_repeatmasker_hit_classification.py"
        )
        completed = subprocess.run(
            [sys.executable, str(script), "--help"],
            cwd=PROJECT_ROOT,
            capture_output=True,
            text=True,
            check=False,
        )

        self.assertEqual(completed.returncode, 0, completed.stderr)

    def test_summary_counts_rows_and_inclusive_bp_without_mixing_sources(self):
        current_combined = pd.DataFrame(
            {
                "repeat_class": ["LTR/Gypsy", "DNA/hAT", "Satellite/centromeric"],
                "Class": ["Retrotransposons Autonomous", "DNAtransposons Subclass1", "Other"],
                "Order": ["LINE", "TIR", pd.NA],
                "Superfamily": ["L2", "hAT", pd.NA],
                "query_start": [1, 11, 31],
                "query_end": [10, 30, 35],
            }
        )

        summary = summarize_chunk(current_combined)

        self.assertEqual(summary["n_hits"], 3)
        self.assertEqual(summary["aligned_bp"], 35)
        self.assertEqual(summary["repeatmasker_unmapped_hits"], 1)
        self.assertEqual(summary["repeatmasker_unmapped_bp"], 5)
        self.assertEqual(summary["order_comparable_hits"], 2)
        self.assertEqual(summary["order_mismatch_hits"], 1)
        self.assertEqual(summary["order_mismatch_bp"], 10)
        self.assertEqual(summary["superfamily_comparable_hits"], 2)
        self.assertEqual(summary["superfamily_mismatch_hits"], 1)


if __name__ == "__main__":
    unittest.main()
