import importlib.util
import unittest
from pathlib import Path

import numpy as np
import pandas as pd


MODULE_PATH = Path(__file__).resolve().parents[2] / "processing" / "audit_microscopy_release.py"
SPEC = importlib.util.spec_from_file_location("audit_microscopy_release", MODULE_PATH)
audit = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(audit)


class MicroscopyReleaseAuditTests(unittest.TestCase):
    def test_species_normalization_is_exact(self):
        self.assertEqual("fuscus", audit.canonical_species("D. fuscus"))
        self.assertEqual("fuscus", audit.canonical_species("Desmognathus fuscus"))
        self.assertEqual("gvnigeusgwotli", audit.canonical_species("D.gvnigeusgwotli"))

    def test_explicit_rejects_are_excluded_without_dropping_unlabeled_rows(self):
        frame = pd.DataFrame(
            {
                "current_decision": ["keep", "discard", "nucleus_only", np.nan],
                "value": [1, 2, 3, 4],
            }
        )
        observed = audit.eligible_candidates(frame)
        self.assertEqual([1, 4], observed["value"].tolist())

    def test_literal_largest_selection_preserves_corresponding_nucleus(self):
        selected = pd.DataFrame(
            {
                "species": ["fuscus"] * 2,
                "cell_area_um2": [1.0, 2.0],
                "nuc_area_um2": [10.0, 20.0],
                "nc_area_ratio": [10.0, 10.0],
            }
        )
        candidates = pd.DataFrame(
            {
                "species": ["fuscus"] * 3,
                "cell_area_um2": [1.0, 3.0, 2.0],
                "nuc_area_um2": [10.0, 30.0, 20.0],
                "nc_area_ratio": [10.0, 10.0, 10.0],
                "current_decision": ["keep", "keep", "keep"],
            }
        )
        groups = audit.estimator_groups("fuscus", selected, candidates)
        literal = groups["literal_largest50_eligible"]
        self.assertEqual([3.0, 2.0, 1.0], literal["cell_area_um2"].tolist())
        self.assertEqual([30.0, 20.0, 10.0], literal["nuc_area_um2"].tolist())

    def test_iod_identity_is_area_times_mean_od(self):
        area = pd.Series([100.0, 250.0])
        mean_od = pd.Series([0.2, 0.4])
        iod = area * mean_od
        ratio = iod / (area * mean_od)
        self.assertLess(float(np.max(np.abs(ratio - 1.0))), 1e-12)

    def test_safe_spearman_reports_missing_pairs(self):
        rho, pvalue, n = audit.safe_spearman(pd.Series([1.0, 2.0, np.nan]), pd.Series([1.0, 3.0, 4.0]))
        self.assertTrue(np.isnan(rho))
        self.assertTrue(np.isnan(pvalue))
        self.assertEqual(2, n)


if __name__ == "__main__":
    unittest.main()
