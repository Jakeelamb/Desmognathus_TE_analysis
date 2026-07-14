#!/usr/bin/env python3
"""Tests for the accession-labeled fuscus resource benchmark."""

import unittest
from pathlib import Path

from scripts.processing.build_fuscus_resource_benchmark import build_benchmark


PROJECT_ROOT = Path(__file__).resolve().parents[3]


class TestFuscusResourceBenchmark(unittest.TestCase):
    def test_benchmark_keeps_denominators_and_accessions_explicit(self):
        benchmark = build_benchmark(PROJECT_ROOT).set_index("metric")

        assembly = benchmark.loc["assembly_span_gb"]
        self.assertEqual(assembly["current_assembly_accession"], "GCA_032353935.1")
        self.assertEqual(assembly["validation_assembly_accession"], "GCA_050004315.1")
        self.assertAlmostEqual(assembly["current_value"], 0.652813091)
        self.assertAlmostEqual(assembly["published_value"], 16.117580045)

        ltr = benchmark.loc["LTR"]
        self.assertAlmostEqual(ltr["current_value"], 64.696348, places=5)
        self.assertEqual(ltr["current_denominator"], "classified_dnaPipeTE_order_percent")
        self.assertEqual(ltr["published_denominator"], "whole_assembly_percent")
        self.assertFalse(bool(ltr["directly_comparable"]))

        repeat_load = benchmark.loc["repeat_aligned_sensitivity"]
        self.assertAlmostEqual(repeat_load["current_value"], 66.0099252, places=5)
        self.assertEqual(
            repeat_load["current_denominator"],
            "percent_of_configured_1.5Gb_dnaPipeTE_quantification_sample",
        )
        self.assertFalse(bool(repeat_load["directly_comparable"]))


if __name__ == "__main__":
    unittest.main()
