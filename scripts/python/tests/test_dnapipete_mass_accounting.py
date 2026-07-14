#!/usr/bin/env python3
"""Tests for final-panel dnaPipeTE mass accounting."""

import unittest

import pandas as pd

from scripts.processing.audit_dnapipete_mass_accounting import (
    summarize_mass,
    validate_analysis_resources,
)


class TestDnaPipeTEMassAccounting(unittest.TestCase):
    def test_unresolved_mass_is_retained_and_out_of_panel_species_are_excluded(self):
        merged = pd.DataFrame(
            {
                "Species": ["D.fuscus", "D.fuscus", "D.fuscus", "D.orestes"],
                "Source": ["SRX1", "SRX1", "SRX1", "SRX2"],
                "aligned_bases": [60, 30, 10, 100],
                "Class": ["Retro", "Retro", "Unknown", "Retro"],
                "Order": ["LTR", "LINE", pd.NA, "LTR"],
                "Superfamily": ["Gypsy", "L2", pd.NA, "Gypsy"],
            }
        )

        summary = summarize_mass(
            merged,
            analysis_species={"fuscus"},
            retained_categories={
                "Class": {"Retro", "Unknown"},
                "Order": {"LTR", "LINE"},
                "Superfamily": {"Gypsy", "L2"},
            },
        ).set_index("species")

        self.assertEqual(summary.index.tolist(), ["fuscus"])
        self.assertEqual(summary.loc["fuscus", "total_aligned_bases"], 100)
        self.assertEqual(summary.loc["fuscus", "order_retained_aligned_bases"], 90)
        self.assertEqual(summary.loc["fuscus", "order_unresolved_aligned_bases"], 10)
        self.assertEqual(summary.loc["fuscus", "order_retained_fraction"], 0.9)
        self.assertEqual(summary.loc["fuscus", "sra_accession"], "SRX1")

    def test_resource_validation_is_accession_specific(self):
        mass = pd.DataFrame(
            {
                "species": ["fuscus"],
                "sra_accession": ["SRX20497025"],
                "total_aligned_bases": [100],
            }
        )
        lookup = pd.DataFrame(
            {
                "Species": ["D.fuscus"],
                "SRA_Accension": ["SRX20497025"],
                "Genome_Accension": ["GCA_032353935.1"],
            }
        )

        validated = validate_analysis_resources(mass, lookup, {"fuscus"})

        self.assertEqual(validated.loc[0, "te_sra_accession"], "SRX20497025")
        self.assertEqual(
            validated.loc[0, "te_assembly_accession"], "GCA_032353935.1"
        )
        self.assertNotIn("sra_accession", validated.columns)

        wrong_lookup = lookup.copy()
        wrong_lookup.loc[0, "SRA_Accension"] = "SRX19953421"
        with self.assertRaisesRegex(ValueError, "do not match the active TE lookup"):
            validate_analysis_resources(mass, wrong_lookup, {"fuscus"})

    def test_retry_outputs_cannot_be_summed_as_independent_sources(self):
        merged = pd.DataFrame(
            {
                "Species": ["D.fuscus", "D.fuscus"],
                "Source": ["SRX20497025", "SRX20497025R2"],
                "aligned_bases": [100, 100],
                "Class": ["Retro", "Retro"],
                "Order": ["LTR", "LTR"],
                "Superfamily": ["Gypsy", "Gypsy"],
            }
        )

        with self.assertRaisesRegex(ValueError, "exactly one SRX"):
            summarize_mass(
                merged,
                analysis_species={"fuscus"},
                retained_categories={
                    "Class": {"Retro"},
                    "Order": {"LTR"},
                    "Superfamily": {"Gypsy"},
                },
            )


if __name__ == "__main__":
    unittest.main()
