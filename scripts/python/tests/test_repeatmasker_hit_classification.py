#!/usr/bin/env python3
"""Behavioral tests for source-explicit RepeatMasker hit classification."""

import unittest
from pathlib import Path
import subprocess
import sys

import pandas as pd

from scripts.processing.te_classification import (
    annotate_repeatmasker_hits,
    classify_repeatmasker_hits,
    prepare_dnapipete_context,
)
from scripts.processing.rebuild_repeatmasker_hit_classification import (
    filter_hits_to_analysis_species,
)


PROJECT_ROOT = Path(__file__).resolve().parents[3]


class TestRepeatMaskerHitClassification(unittest.TestCase):
    def test_corrected_rebuild_excludes_species_outside_final_panel(self):
        hits = pd.DataFrame(
            {
                "Desmognathus_Species": ["D.fuscus", "D.planiceps", "D.welteri"],
                "repeat_class": ["LTR/Gypsy", "LINE/L2", "DNA/hAT"],
            }
        )

        retained = filter_hits_to_analysis_species(hits, {"fuscus", "welteri"})

        self.assertEqual(
            retained["Desmognathus_Species"].tolist(), ["D.fuscus", "D.welteri"]
        )

    def test_non_destructive_rebuild_entrypoint_is_runnable(self):
        script = (
            PROJECT_ROOT
            / "scripts"
            / "processing"
            / "rebuild_repeatmasker_hit_classification.py"
        )
        completed = subprocess.run(
            [sys.executable, str(script), "--help"],
            cwd=PROJECT_ROOT,
            capture_output=True,
            text=True,
            check=False,
        )

        self.assertEqual(completed.returncode, 0, completed.stderr)

    def test_merge_keeps_each_hit_class_and_namespaces_contig_context(self):
        hits = pd.DataFrame(
            {
                "query_name": ["contig_1", "contig_1"],
                "SRX_ID": ["SRX1", "SRX1"],
                "repeat_class": ["LTR/Gypsy", "DNA/hAT"],
            }
        )
        dnapipete = pd.DataFrame(
            {
                "dnaPipeTE_contig_name": ["contig_1", "contig_1"],
                "Source": ["SRX1", "SRX1"],
                "RM_classification": ["LTR/Gypsy", "LINE/L2"],
                "Class": ["Retrotransposons Autonomous", "Retrotransposons Autonomous"],
                "Order": ["LTR", "LINE"],
                "Superfamily": ["Gypsy", "L2"],
                "hitlength_contiglength": [0.4, 0.9],
            }
        )

        context = prepare_dnapipete_context(dnapipete)
        annotated = annotate_repeatmasker_hits(hits, context)

        self.assertEqual(len(annotated), 2)
        self.assertEqual(annotated["repeatmasker_order"].tolist(), ["LTR", "TIR"])
        self.assertEqual(annotated["dnapipete_order"].tolist(), ["LINE", "LINE"])
        self.assertNotIn("Order_dnapipete", annotated.columns)

    def test_each_hit_keeps_its_own_classification_after_contig_merge(self):
        hits = pd.DataFrame(
            {
                "query_name": ["contig_1", "contig_1"],
                "repeat_class": ["LTR/Gypsy", "DNA/hAT"],
                "dnapipete_repeat_class": ["LINE/L2", "LINE/L2"],
                "dnapipete_order": ["LINE", "LINE"],
                "dnapipete_superfamily": ["L2", "L2"],
            }
        )

        classified = classify_repeatmasker_hits(hits)

        self.assertEqual(classified["repeatmasker_order"].tolist(), ["LTR", "TIR"])
        self.assertEqual(classified["repeatmasker_superfamily"].tolist(), ["Gypsy", "hAT"])
        self.assertEqual(classified["Order"].tolist(), ["LTR", "TIR"])
        self.assertEqual(classified["Superfamily"].tolist(), ["Gypsy", "hAT"])
        self.assertEqual(classified["dnapipete_order"].tolist(), ["LINE", "LINE"])
        self.assertEqual(classified["repeat_class"].tolist(), ["LTR/Gypsy", "DNA/hAT"])
        self.assertEqual(len(classified), len(hits))
        self.assertTrue(classified.index.equals(hits.index))

    def test_partial_and_unmapped_labels_remain_explicit(self):
        hits = pd.DataFrame(
            {"repeat_class": ["LTR/ERV1", "Satellite/centromeric"]}
        )

        classified = classify_repeatmasker_hits(hits)

        self.assertEqual(
            classified["repeatmasker_te_class"].tolist(),
            ["Retrotransposons Autonomous", "Unclassified"],
        )
        self.assertEqual(
            classified["repeatmasker_order"].tolist(),
            ["LTR", "Unclassified"],
        )
        self.assertEqual(
            classified["repeatmasker_superfamily"].tolist(),
            ["Unclassified", "Unclassified"],
        )
        self.assertEqual(
            classified["repeatmasker_classification_status"].tolist(),
            ["order_only", "unmapped"],
        )


if __name__ == "__main__":
    unittest.main()
