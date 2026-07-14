from __future__ import annotations

import json
import unittest
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[3]
DECLARATION = ROOT / "path_analysis/data/templates/analysis_species_panels.csv"
LOOKUP = ROOT / "input_data/lookup_table.txt"
PANEL_DIR = ROOT / "path_analysis/data/derived/panels"
MANIFEST = PANEL_DIR / "study_species_panels_v1.manifest.json"

PATH18 = {
    "amphileucus", "anicetus", "apalachicolae", "auriculatus", "bairdi",
    "campi", "fuscus", "gvnigeusgwotli", "intermedius", "kanawha",
    "marmoratus", "mavrokoilius", "monticola", "ocoee", "perlapsus",
    "tilleyi", "valtos", "welteri",
}
CELL_ONLY = {"brimleyorum", "folkertsi", "ochrophaeus"}
CELL21 = PATH18 | CELL_ONLY
TE_ONLY_WITH_LEGACY_NUCLEI = {"aeneus", "organi", "wrighti"}


def canonical(values: pd.Series) -> set[str]:
    return set(
        values.astype(str).str.strip().str.replace("Desmognathus ", "", regex=False)
        .str.replace("D.", "", regex=False).str.lower()
    )


class StudySpeciesPanelTests(unittest.TestCase):
    def test_declared_panels_match_vetted_resources_and_evidence_intersection(self) -> None:
        declared = pd.read_csv(DECLARATION)
        self.assertTrue(declared["species"].is_unique)
        te34 = set(declared.loc[declared["te_resource_panel34_v1"], "species"])
        cell21 = set(declared.loc[declared["cell_linked_panel21_v1"], "species"])
        path18 = set(declared.loc[declared["integrated_path_panel18_v1"], "species"])
        active_lookup = canonical(pd.read_csv(LOOKUP, sep="\t")["Species"])

        self.assertEqual(len(te34), 34)
        self.assertEqual(te34, active_lookup)
        self.assertEqual(cell21, CELL21)
        self.assertEqual(path18, PATH18)
        self.assertEqual(path18, te34 & cell21)
        self.assertEqual(cell21 - te34, CELL_ONLY)
        self.assertEqual(len(te34 - path18), 16)
        self.assertTrue(TE_ONLY_WITH_LEGACY_NUCLEI.issubset(te34 - path18))
        self.assertNotIn("planiceps", te34)
        self.assertNotIn("aeaneus", te34)

    def test_derived_panels_preserve_accessions_and_do_not_impute_microscopy(self) -> None:
        te = pd.read_csv(PANEL_DIR / "study_te_resource_panel34_v1.csv")
        cell = pd.read_csv(PANEL_DIR / "study_cell_linked_panel21_v1.csv")
        path = pd.read_csv(PANEL_DIR / "study_integrated_path_panel18_v1.csv")
        self.assertEqual(len(te), 34)
        self.assertEqual(set(cell.species), CELL21)
        self.assertEqual(set(path.species), PATH18)
        self.assertTrue(te[["te_sra_accession", "te_assembly_accession"]].notna().all().all())
        self.assertTrue(te[["has_tree_tip", "has_te"]].all().all())
        self.assertTrue(cell[["has_genome", "has_morphology"]].all().all())
        self.assertTrue(path[["has_genome", "has_morphology"]].all().all())

        extras = te.set_index("species").loc[sorted(TE_ONLY_WITH_LEGACY_NUCLEI)]
        self.assertFalse(extras["has_genome"].any())
        self.assertFalse(extras["has_morphology"].any())
        expected = {
            "aeneus": ("SRX19953415", "GCA_030264635.1"),
            "organi": ("SRX19952929", "GCA_030180145.1"),
            "wrighti": ("SRX19958881", "GCA_030265035.1"),
        }
        for species, accessions in expected.items():
            self.assertEqual(
                (extras.loc[species, "te_sra_accession"], extras.loc[species, "te_assembly_accession"]),
                accessions,
            )

    def test_manifest_records_three_separate_panels_without_imputation(self) -> None:
        manifest = json.loads(MANIFEST.read_text())
        self.assertEqual(manifest["te_resource_panel"]["n_species"], 34)
        self.assertEqual(manifest["cell_linked_panel"]["n_species"], 21)
        self.assertEqual(manifest["integrated_path_panel"]["n_species"], 18)
        self.assertEqual(manifest["cell_only_species_without_active_te_resource"], sorted(CELL_ONLY))
        self.assertEqual(len(manifest["te_only_species_without_current_linked_microscopy"]), 16)
        self.assertTrue(TE_ONLY_WITH_LEGACY_NUCLEI.issubset(
            manifest["te_only_species_without_current_linked_microscopy"]
        ))
        self.assertFalse(manifest["microscopy_traits_imputed_for_te_only_species"])
        self.assertTrue(manifest["analysis18_outputs_preserved"])


if __name__ == "__main__":
    unittest.main()
