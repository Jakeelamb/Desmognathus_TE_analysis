#!/usr/bin/env python3

from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path

import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(PROJECT_ROOT / "path_analysis/scripts"))

import finalize_microscopy_review_extensions as finalizer  # noqa: E402


class FinalizeMicroscopyReviewExtensionsTests(unittest.TestCase):
    def test_default_iod_review_species_include_recovered_ochrophaeus(self):
        self.assertEqual(
            finalizer.IOD_REVIEW_SPECIES,
            ("D. aeneus", "D. orestes", "D. wrighti", "D. ochrophaeus"),
        )

    def test_legacy_frozen_rows_are_preserved_and_only_new_keeps_are_added(self):
        frozen = pd.DataFrame(
            {
                "species": ["D. old", "D. old"],
                "review_key": ["old-1", "old-2"],
                "measurement": [10.0, 20.0],
                "panel": ["legacy", "legacy"],
            }
        )
        extension = pd.DataFrame(
            {
                "species": ["D. old", "D. old", "D. new", "D. new"],
                "review_key": ["old-1", "old-2", "new-1", "new-2"],
                "measurement": [10.0, 20.0, 30.0, 40.0],
                "panel": ["extension"] * 4,
                "decision": ["stale"] * 4,
                "review_decision": ["stale"] * 4,
            }
        )
        decisions = pd.DataFrame(
            {
                "species": ["D. old", "D. old", "D. new", "D. new"],
                "review_key": ["old-1", "old-2", "new-1", "new-2"],
                "decision": ["keep", "keep", "keep", "problem"],
                "decision_source": [
                    "initial_default",
                    "initial_default",
                    "auto_fill",
                    "explicit",
                ],
            }
        )

        merged, audit = finalizer.finalize_review_table(
            frozen=frozen,
            extension=extension,
            decisions=decisions,
            new_species=("D. new",),
            panel="finalized_test",
            decision_file="decisions.csv",
        )

        self.assertEqual(merged["review_key"].tolist(), ["new-1", "old-1", "old-2"])
        self.assertEqual(merged.set_index("review_key").loc["old-1", "measurement"], 10.0)
        self.assertTrue(merged["review_status"].eq("reviewed_keep").all())
        self.assertEqual(audit.set_index("review_key").loc["new-1", "included_in_finalized"], True)
        self.assertEqual(audit.set_index("review_key").loc["new-2", "included_in_finalized"], False)

    def test_every_requested_new_species_must_have_exported_decisions(self):
        frozen = pd.DataFrame({"species": ["D. old"], "review_key": ["old-1"]})
        extension = pd.DataFrame(
            {
                "species": ["D. new", "D. missing"],
                "review_key": ["new-1", "missing-1"],
            }
        )
        decisions = pd.DataFrame(
            {
                "species": ["D. new"],
                "review_key": ["new-1"],
                "decision": ["keep"],
                "decision_source": ["explicit"],
            }
        )

        with self.assertRaisesRegex(ValueError, "missing exported decisions"):
            finalizer.finalize_review_table(
                frozen=frozen,
                extension=extension,
                decisions=decisions,
                new_species=("D. new", "D. missing"),
                panel="finalized_test",
                decision_file="decisions.csv",
            )

    def test_release_writer_emits_finalized_table_audit_and_hashes(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            frozen_path = root / "frozen.csv.gz"
            extension_path = root / "extension.csv.gz"
            decisions_path = root / "decisions.csv"
            output_path = root / "finalized.csv.gz"
            audit_path = root / "audit.csv"
            pd.DataFrame(
                {"species": ["D. old"], "review_key": ["old-1"], "value": [1.0]}
            ).to_csv(frozen_path, index=False)
            pd.DataFrame(
                {
                    "species": ["D. new", "D. new"],
                    "review_key": ["new-1", "new-2"],
                    "value": [2.0, 3.0],
                }
            ).to_csv(extension_path, index=False)
            pd.DataFrame(
                {
                    "species": ["D. new", "D. new"],
                    "review_key": ["new-1", "new-2"],
                    "decision": ["keep", "problem"],
                    "decision_source": ["auto_fill", "explicit"],
                }
            ).to_csv(decisions_path, index=False)

            summary = finalizer.write_finalized_review_release(
                frozen_path=frozen_path,
                extension_path=extension_path,
                decisions_path=decisions_path,
                output_path=output_path,
                audit_path=audit_path,
                new_species=("D. new",),
                panel="finalized_test",
            )

            self.assertTrue(output_path.exists())
            self.assertTrue(audit_path.exists())
            self.assertEqual(summary["finalized_rows"], 2)
            self.assertEqual(summary["new_species_selected_counts"], {"D. new": 1})
            self.assertEqual(summary["output_sha256"], finalizer.sha256_file(output_path))

    def test_unknown_review_decisions_are_rejected(self):
        frozen = pd.DataFrame({"species": ["D. old"], "review_key": ["old-1"]})
        extension = pd.DataFrame({"species": ["D. new"], "review_key": ["new-1"]})
        decisions = pd.DataFrame(
            {
                "species": ["D. new"],
                "review_key": ["new-1"],
                "decision": ["maybe"],
                "decision_source": ["explicit"],
            }
        )

        with self.assertRaisesRegex(ValueError, "Unsupported review decisions"):
            finalizer.finalize_review_table(
                frozen=frozen,
                extension=extension,
                decisions=decisions,
                new_species=("D. new",),
                panel="finalized_test",
                decision_file="decisions.csv",
            )

    def test_incremental_in_place_iod_finalization_preserves_prior_rows_and_hash(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            finalized_path = root / "finalized.csv.gz"
            extension_path = root / "extension.csv.gz"
            decisions_path = root / "decisions.csv"
            audit_path = root / "audit.csv"
            pd.DataFrame(
                {"species": ["D. old"], "review_key": ["old-1"], "value": [1.0]}
            ).to_csv(finalized_path, index=False)
            original_hash = finalizer.sha256_file(finalized_path)
            pd.DataFrame(
                {
                    "species": ["D. added", "D. added"],
                    "review_key": ["added-1", "added-2"],
                    "value": [2.0, 3.0],
                }
            ).to_csv(extension_path, index=False)
            pd.DataFrame(
                {
                    "species": ["D. added", "D. added"],
                    "review_key": ["added-1", "added-2"],
                    "decision": ["keep", "problem"],
                    "decision_source": ["explicit", "explicit"],
                }
            ).to_csv(decisions_path, index=False)

            summary = finalizer.write_finalized_review_release(
                frozen_path=finalized_path,
                extension_path=extension_path,
                decisions_path=decisions_path,
                output_path=finalized_path,
                audit_path=audit_path,
                new_species=("D. added",),
                panel="finalized_test",
            )

            merged = pd.read_csv(finalized_path)
            self.assertEqual(set(merged["species"]), {"D. old", "D. added"})
            self.assertEqual(summary["frozen_source_sha256"], original_hash)
            self.assertNotEqual(summary["output_sha256"], original_hash)


if __name__ == "__main__":
    unittest.main()
