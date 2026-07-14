#!/usr/bin/env python3
"""Focused guards for the final-panel phylogeny release audit."""

import unittest
from io import StringIO

import numpy as np
from Bio import Phylo

from scripts.processing.audit_phylogeny_release import (
    canonical_species,
    internal_clade_signature,
    pad_terminal_rounding,
    prune_to_species,
    root_to_tip_depths,
)


class TestPhylogenyReleaseAudit(unittest.TestCase):
    def test_species_normalization_does_not_collapse_lineage_suffixes(self):
        self.assertEqual(canonical_species("D. fuscus"), "fuscus")
        self.assertEqual(canonical_species("fuscus_A"), "fuscus_a")

    def test_exact_pruning_preserves_requested_tips_without_substitution(self):
        tree = Phylo.read(StringIO("((fuscus:1,planiceps:1):1,welteri:2);"), "newick")
        pruned = prune_to_species(tree, ["fuscus", "welteri"])
        self.assertEqual({tip.name for tip in pruned.get_terminals()}, {"fuscus", "welteri"})

    def test_terminal_padding_makes_rounding_only_tree_ultrametric(self):
        tree = Phylo.read(StringIO("((a:1,b:0.999998):1,c:2);"), "newick")
        corrected, ledger = pad_terminal_rounding(tree)
        depths = np.asarray(list(root_to_tip_depths(corrected).values()))
        self.assertLessEqual(np.ptp(depths), 1e-12)
        self.assertAlmostEqual(ledger["terminal_padding_years"].max(), 2.0)

    def test_clade_signature_ignores_branch_scale(self):
        first = Phylo.read(StringIO("((a:1,b:1):1,c:2);"), "newick")
        second = Phylo.read(StringIO("((a:0.8,b:0.8):0.8,c:1.6);"), "newick")
        self.assertEqual(internal_clade_signature(first), internal_clade_signature(second))


if __name__ == "__main__":
    unittest.main()
