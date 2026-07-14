import importlib.util
import unittest
from pathlib import Path

import numpy as np
import pandas as pd


MODULE_PATH = Path(__file__).resolve().parents[2] / "processing" / "build_corrected_path_inputs.py"
SPEC = importlib.util.spec_from_file_location("build_corrected_path_inputs", MODULE_PATH)
module = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(module)


class CorrectedPathInputTests(unittest.TestCase):
    def test_order_features_use_exact_ltr_line_logratio(self):
        composition = pd.DataFrame(
            {
                "te_level": ["order", "order", "order", "order"],
                "composition_mode": ["classified_conditional"] * 4,
                "species": ["a", "a", "b", "b"],
                "feature": ["LTR", "LINE", "LTR", "LINE"],
                "proportion": [0.6, 0.3, 0.2, 0.4],
            }
        )
        diversity = pd.DataFrame(
            {
                "te_level": ["order", "order"],
                "composition_mode": ["classified_conditional"] * 2,
                "species": ["a", "b"],
                "observed_richness": [2, 2],
                "shannon_entropy": [0.5, 0.6],
                "gini_simpson": [0.4, 0.5],
                "hill_q1": [2, 2],
                "hill_q2": [2, 2],
                "pielou_evenness": [0.7, 0.8],
            }
        )
        observed = module.order_features(composition, diversity).set_index("species")
        self.assertAlmostEqual(np.log(2), observed.loc["a", "ltr_line_logratio"])
        self.assertAlmostEqual(np.log(0.5), observed.loc["b", "ltr_line_logratio"])

    def test_measurement_cross_preserves_all_specifications(self):
        morph = pd.DataFrame(
            {
                "species": ["a", "a", "b", "b"],
                "morphology_estimator": ["m1", "m2", "m1", "m2"],
            }
        )
        iod = pd.DataFrame(
            {
                "species": ["a", "a", "b", "b"],
                "iod_subset": ["i1", "i2", "i1", "i2"],
            }
        )
        observed = module.cross_measurement_specifications(morph, iod)
        self.assertEqual(8, len(observed))
        self.assertEqual(4, len(observed[observed["species"].eq("a")]))


if __name__ == "__main__":
    unittest.main()
