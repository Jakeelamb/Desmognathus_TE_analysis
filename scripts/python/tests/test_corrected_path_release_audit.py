from __future__ import annotations

import json
import unittest
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[3]
DATA = ROOT / "results/data/corrected/path_analysis"


class CorrectedPathReleaseAuditTests(unittest.TestCase):
    def test_rankings_are_complete_fail_closed_and_do_not_use_absolute_genome_size(self) -> None:
        expected_fits = {
            "data_sensitivity": 84,
            "tree_sensitivity": 804,
            "leave_one_out": 70,
        }
        for phase, n_fits in expected_fits.items():
            rankings = pd.read_csv(
                DATA / f"corrected_path_{phase}_rankings_analysis18_v1.csv"
            )
            self.assertEqual(rankings["fit_id"].nunique(), n_fits)
            self.assertTrue(rankings[["CICc", "delta_CICc", "w"]].notna().all().all())
            self.assertFalse(rankings["absolute_genome_size_used"].astype(bool).any())
            self.assertFalse(rankings["publication_claim_allowed"].astype(bool).any())
            for _, group in rankings.groupby("fit_id"):
                self.assertAlmostEqual(float(group["w"].sum()), 1.0, places=6)

    def test_each_family_has_a_null_and_terminal_proxy_is_rejected(self) -> None:
        rankings = pd.read_csv(
            DATA / "corrected_path_data_sensitivity_rankings_analysis18_v1.csv"
        )
        nulls = {
            "te_iod": "proxy_null",
            "iod_morphology": "morphology_null",
            "integrated": "integrated_null",
            "terminal_internal_iod": "proxy_null",
        }
        for family, model in nulls.items():
            self.assertIn(model, set(rankings.loc[rankings["family"].eq(family), "model"]))
        terminal_top = rankings[
            rankings["family"].eq("terminal_internal_iod") & rankings["rank"].eq(1)
        ]
        self.assertEqual(len(terminal_top), 6)
        self.assertTrue(
            terminal_top["release_gate_status"]
            .eq("blocked_no_globally_supported_model")
            .all()
        )

    def test_tree_and_leave_one_out_sensitivity_are_complete(self) -> None:
        tree = pd.read_csv(
            DATA / "corrected_path_tree_sensitivity_rankings_analysis18_v1.csv"
        )
        self.assertEqual(tree["tree_id"].nunique(), 201)
        self.assertEqual(tree["fit_id"].nunique(), 804)
        loo = pd.read_csv(
            DATA / "corrected_path_leave_one_out_rankings_analysis18_v1.csv"
        )
        expected = {"te_iod": 18, "iod_morphology": 18, "integrated": 18, "terminal_internal_iod": 16}
        observed = loo[loo["rank"].eq(1)].groupby("family")["fit_id"].nunique().to_dict()
        self.assertEqual(observed, expected)

    def test_simulation_calibration_has_full_design_and_edge_coverage(self) -> None:
        summary = pd.read_csv(DATA / "corrected_path_simulation_calibration_analysis18_v1.csv")
        self.assertEqual(len(summary), 12)
        null = summary[summary["scenario"].eq("independent_null")]
        signal = summary[summary["scenario"].eq("observed_chain")]
        self.assertTrue(null["n_attempted"].eq(200).all())
        self.assertTrue(signal["n_attempted"].eq(100).all())
        self.assertTrue(summary["n_failures"].eq(0).all())
        self.assertTrue(summary["primary_calibration_rate"].between(0, 1).all())
        edge = pd.read_csv(
            DATA / "corrected_path_simulation_edge_calibration_analysis18_v1.csv"
        )
        self.assertEqual(len(edge), 24)
        self.assertTrue(edge["n_successful"].eq(100).all())
        self.assertTrue(edge["interval_coverage_rate"].between(0, 1).all())

    def test_release_manifest_is_explicitly_exploratory(self) -> None:
        manifest = json.loads(
            (DATA / "corrected_path_release_audit_analysis18_v1.manifest.json").read_text()
        )
        self.assertEqual(manifest["n_path_fits"], 958)
        self.assertEqual(manifest["n_fit_failures"], 0)
        self.assertFalse(manifest["absolute_genome_size_used"])
        self.assertFalse(manifest["causal_claim_allowed"])


if __name__ == "__main__":
    unittest.main()
