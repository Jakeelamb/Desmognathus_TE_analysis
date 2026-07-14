from __future__ import annotations

import json
import unittest
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[3]
DATA = ROOT / "results/data/corrected/te34_replicate_averaged"


class Te34ReplicateAveragedBundleTests(unittest.TestCase):
    def test_orestes_retry_runs_are_one_equal_weight_species_estimate(self) -> None:
        mass = pd.read_csv(DATA / "dnapipete_mass_accounting_te34_replicate_averaged_v1.csv")
        row = mass.set_index("species").loc["orestes"]
        self.assertEqual(int(row["n_dnapipete_runs"]), 2)
        self.assertEqual(row["run_aggregation"], "equal_weight_mean_of_run_level_estimates")
        self.assertEqual(set(json.loads(row["dnapipete_source_labels"])), {"SRX19952890", "SRX19952890R2"})
        self.assertEqual(mass.species.nunique(), 34)

    def test_manifest_records_a_mean_not_double_count_policy(self) -> None:
        manifest = json.loads((DATA / "te34_replicate_averaged_v1.manifest.json").read_text())
        self.assertEqual(manifest["dnapipete_run_aggregation"], "equal_weight_mean_of_run_level_estimates")
        self.assertEqual(manifest["species_with_multiple_dnapipete_runs"], ["orestes"])
        self.assertFalse(manifest["retry_runs_counted_as_additional_species"])
        self.assertTrue(manifest["exact_srx_te34_v1_preserved"])


if __name__ == "__main__":
    unittest.main()
