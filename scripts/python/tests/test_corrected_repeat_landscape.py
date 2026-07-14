from __future__ import annotations

import json
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

from scripts.processing.build_corrected_repeat_landscape import summarize_chunk, validate


ROOT = Path(__file__).resolve().parents[3]
OUTPUT_DIR = ROOT / "results/data/corrected/repeat_landscape"


class CorrectedRepeatLandscapeTests(unittest.TestCase):
    def test_chunk_summary_bins_and_retains_unclassified_hits(self) -> None:
        source = pd.DataFrame(
            {
                "Desmognathus_Species": ["D.fuscus", "D.fuscus", "D.fuscus"],
                "percent_divergence": [0.2, 12.9, 55.1],
                "query_start": [1, 11, 31],
                "query_end": [10, 30, 60],
                "Order": ["LTR", np.nan, "LINE"],
            }
        )
        result = summarize_chunk(source)
        self.assertEqual(int(result.hit_bp.sum()), 60)
        self.assertEqual(set(result.divergence_bin_start_pct), {0, 12, 50})
        self.assertIn("Unclassified", set(result.order))

    def test_frozen_landscape_conserves_hit_bp_and_fractions(self) -> None:
        frame = pd.read_csv(OUTPUT_DIR / "repeatmasker_divergence_landscape_analysis18_v1.csv")
        source_manifest = json.loads(
            (
                ROOT
                / "results/data/corrected/repeatmasker_detailed_classification_hit_level_analysis18_v1.csv.manifest.json"
            ).read_text()
        )
        result = validate(frame, int(source_manifest["inclusive_aligned_bp"]))
        self.assertEqual(result["n_species"], 18)
        self.assertEqual(result["hit_count"], int(source_manifest["n_hits"]))
        self.assertLess(result["maximum_fraction_sum_error"], 1e-12)

    def test_manifest_marks_downstream_cache_only(self) -> None:
        manifest = json.loads(
            (OUTPUT_DIR / "repeatmasker_divergence_landscape_analysis18_v1.manifest.json").read_text()
        )
        self.assertFalse(manifest["expensive_upstream_tools_executed"])
        self.assertEqual(
            manifest["abundance_denominator"],
            "total inclusive corrected RepeatMasker query-coordinate aligned bp within species; "
            "overlapping hits are not interval-deduplicated",
        )


if __name__ == "__main__":
    unittest.main()
