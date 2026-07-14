from __future__ import annotations

import hashlib
import json
import unittest
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[3]
DATASETS = ROOT / "Publication" / "datasets"


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


class PublicationPackageTests(unittest.TestCase):
    def test_analysis_availability_flags_match_frozen_panels(self) -> None:
        frame = pd.read_csv(
            DATASETS / "Supplementary_Data_S01_species_analysis_availability.csv"
        )
        self.assertEqual(int(frame["Genomic_TE_resource34"].sum()), 34)
        self.assertEqual(int(frame["Microscopy_cell_linked21"].sum()), 21)
        self.assertEqual(int(frame["Conserved_integrated_path18"].sum()), 18)
        self.assertFalse("Microscopyl_cell_linked21" in frame.columns)

    def test_manifest_matches_every_exported_csv(self) -> None:
        manifest = pd.read_csv(DATASETS / "DATASET_MANIFEST.csv")
        self.assertEqual(len(manifest), 39)
        self.assertFalse(manifest["supplement_id"].duplicated().any())
        self.assertFalse(manifest["filename"].duplicated().any())
        for row in manifest.itertuples(index=False):
            path = DATASETS / row.filename
            self.assertTrue(path.exists(), path)
            self.assertEqual(sha256_file(path), row.output_sha256)
            frame = pd.read_csv(path, low_memory=False)
            self.assertEqual(len(frame), row.rows)
            self.assertEqual(len(frame.columns), row.columns)

    def test_frozen_morphology_release_is_exactly_fifty_per_species(self) -> None:
        frame = pd.read_csv(
            DATASETS / "Supplementary_Data_S18_frozen_cell_nucleus_objects.csv",
            low_memory=False,
        )
        self.assertEqual(len(frame), 1050)
        self.assertEqual(frame["species"].nunique(), 21)
        self.assertTrue(frame.groupby("species").size().eq(50).all())
        self.assertTrue(frame["frozen_selected"].astype(bool).all())
        self.assertTrue(frame["decision"].eq("keep").all())

    def test_main_morphology_table_does_not_relabel_iod_as_genome_size(self) -> None:
        frame = pd.read_csv(
            DATASETS / "Supplementary_Data_S19_cell_nucleus_species_estimates.csv"
        )
        self.assertEqual(len(frame), 21)
        self.assertFalse(any("genome" in column.lower() for column in frame.columns))
        self.assertFalse(any("iod" in column.lower() for column in frame.columns))

    def test_nuclear_iod_panel_contains_only_reviewed_keeps(self) -> None:
        frame = pd.read_csv(
            DATASETS / "Supplementary_Data_S21_frozen_nuclear_iod_objects.csv",
            low_memory=False,
        )
        self.assertEqual(len(frame), 721)
        self.assertTrue(frame["review_status"].eq("reviewed_keep").all())

    def test_tree_standard_formats_and_csv_views_are_present(self) -> None:
        non_csv = json.loads((DATASETS / "NON_CSV_FILES.json").read_text())
        self.assertEqual(len(non_csv), 3)
        for row in non_csv:
            path = DATASETS / row["filename"]
            self.assertTrue(path.exists())
            self.assertEqual(sha256_file(path), row["sha256"])
        edges = pd.read_csv(DATASETS / "Supplementary_Data_S38_phylogeny_edges.csv")
        nodes = pd.read_csv(DATASETS / "Supplementary_Data_S39_phylogeny_nodes.csv")
        self.assertEqual(int(nodes["is_tip"].sum()), 18)
        self.assertEqual(len(edges), len(nodes) - 1)


if __name__ == "__main__":
    unittest.main()
