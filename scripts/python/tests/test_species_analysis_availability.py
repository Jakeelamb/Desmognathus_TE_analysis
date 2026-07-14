from __future__ import annotations

import sys
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "scripts/processing"))
import build_species_analysis_availability as report  # noqa: E402


class SpeciesAnalysisAvailabilityTests(unittest.TestCase):
    def test_amphibiaweb_svl_audit_covers_te34_without_conflating_total_length(self) -> None:
        import pandas as pd

        svl = pd.read_csv(
            ROOT / "results/data/research_review/amphibiaweb_svl_te34_v1.csv"
        )
        te = pd.read_csv(
            ROOT
            / "path_analysis/data/derived/panels/study_te_resource_panel34_v1.csv"
        )

        self.assertEqual(len(svl), 34)
        self.assertTrue(svl["species"].is_unique)
        self.assertEqual(set(svl["species"]), set(te["species"]))
        self.assertEqual(int(svl["amphibiaweb_svl_reported"].sum()), 5)
        self.assertEqual(int(svl["analysis_ready_direct_svl"].sum()), 4)
        self.assertEqual(
            int(svl["account_result"].eq("length_only_not_svl").sum()), 8
        )
        self.assertTrue(
            svl.loc[
                svl["account_result"].eq("length_only_not_svl"),
                "direct_adult_svl_max_mm",
            ].isna().all()
        )
        amphileucus = svl.set_index("species").loc["amphileucus"]
        self.assertEqual(amphileucus["amphibiaweb_taxon"], "Desmognathus amphileucas")
        self.assertEqual(amphileucus["taxonomy_match"], "spelling_mismatch")

    def test_union_and_analysis_step_denominators_match_current_artifacts(self) -> None:
        table = report.build_availability()

        self.assertEqual(len(table), 37)
        self.assertTrue(table["species"].is_unique)
        self.assertEqual(int(table["panel_te_resource34"].sum()), 34)
        self.assertEqual(int(table["panel_cell_linked21"].sum()), 21)
        self.assertEqual(int(table["panel_integrated_path18"].sum()), 18)
        expected_step_counts = {
            "step01_phylogeny_tree_trimming": 18,
            "step02_panel_and_resource_provenance": 37,
            "step03_te_repeat_analysis": 34,
            "step04_ltr_deletion_footprint": 30,
            "step05_cell_modeling_and_measurement": 24,
            "step06_genome_size_estimation": 24,
            "step07_cell_nucleus_genome_path": 24,
            "step08_integrated_phylogenetic_path": 18,
        }
        self.assertEqual(
            {column: int(table[column].sum()) for column in expected_step_counts},
            expected_step_counts,
        )

    def test_known_cross_panel_gaps_remain_visible(self) -> None:
        table = report.build_availability().set_index("species")

        self.assertEqual(
            set(table.index[table["step07_cell_nucleus_genome_path"]])
            - set(table.index[table["step08_integrated_phylogenetic_path"]]),
            {"aeneus", "brimleyorum", "folkertsi", "ochrophaeus", "orestes", "wrighti"},
        )
        self.assertTrue(table.loc["ochrophaeus", "step05_cell_modeling_and_measurement"])
        self.assertTrue(table.loc["ochrophaeus", "step06_genome_size_estimation"])
        self.assertEqual(
            set(table.index[table["step08_integrated_phylogenetic_path"]])
            - set(table.index[table["step08_terminal_internal_proxy_available"].fillna(False)]),
            {"kanawha", "valtos"},
        )


if __name__ == "__main__":
    unittest.main()
