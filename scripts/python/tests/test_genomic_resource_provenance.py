#!/usr/bin/env python3
"""Guards for accession-specific genomic-resource decisions."""

import unittest
from pathlib import Path
import sys

import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(PROJECT_ROOT / "path_analysis" / "scripts"))

from build_master_dataset import load_te_resource_mapping  # noqa: E402

SOURCE_MANIFEST = PROJECT_ROOT / "path_analysis" / "data" / "templates" / "source_manifest.csv"
LOOKUP_TABLE = PROJECT_ROOT / "input_data" / "lookup_table.txt"
CELLPROFILER_TABLES = (
    PROJECT_ROOT
    / "path_analysis"
    / "data"
    / "external"
    / "derived"
    / "cellprofiler_final_species_results.csv",
    PROJECT_ROOT
    / "path_analysis"
    / "data"
    / "external"
    / "derived"
    / "cellprofiler_species_morphology_summary.csv",
)


class TestGenomicResourceProvenance(unittest.TestCase):
    def setUp(self):
        self.manifest = pd.read_csv(SOURCE_MANIFEST).set_index("source_id")
        self.lookup = pd.read_csv(LOOKUP_TABLE, sep="\t")

    def test_fuscus_resource_roles_are_accession_explicit(self):
        current = self.manifest.loc[
            "genomic_fuscus_srx20497025_gca032353935_1"
        ]
        excluded = self.manifest.loc[
            "genomic_fuscus_srx19953421_gca030265095_1"
        ]
        validation = self.manifest.loc["genomic_fuscus_gca050004315_1"]

        self.assertEqual(current["analysis_taxon_name"], "Desmognathus fuscus")
        self.assertEqual(current["inclusion_decision"], "included_current")
        self.assertEqual(excluded["inclusion_decision"], "excluded")
        self.assertEqual(validation["inclusion_decision"], "validation_only")

        active_fuscus = self.lookup.loc[self.lookup["Species"].eq("D.fuscus")].iloc[0]
        self.assertEqual(active_fuscus["SRA_Accension"], current["sra_accessions"])
        self.assertEqual(
            active_fuscus["Genome_Accension"], current["assembly_accession"]
        )
        self.assertNotIn(
            excluded["assembly_accession"], self.lookup["Genome_Accension"].tolist()
        )
        self.assertNotIn(
            validation["assembly_accession"], self.lookup["Genome_Accension"].tolist()
        )

    def test_active_lookup_loads_as_te_namespaced_species_metadata(self):
        resources = load_te_resource_mapping(LOOKUP_TABLE).set_index("species")

        self.assertEqual(resources.loc["fuscus", "te_sra_accession"], "SRX20497025")
        self.assertEqual(
            resources.loc["fuscus", "te_assembly_accession"], "GCA_032353935.1"
        )
        self.assertNotIn("GCA_030265095.1", resources["te_assembly_accession"].tolist())
        self.assertNotIn("GCA_050004315.1", resources["te_assembly_accession"].tolist())

    def test_microscopy_tables_do_not_claim_genomic_accessions(self):
        prohibited = {
            "sra_accession",
            "sra_accension",
            "genome_accession",
            "genome_accension",
            "assembly_accession",
            "te_sra_accession",
            "te_assembly_accession",
        }

        for path in CELLPROFILER_TABLES:
            columns = {column.lower() for column in pd.read_csv(path, nrows=0).columns}
            self.assertTrue(
                prohibited.isdisjoint(columns),
                f"{path.name} must describe microscopy specimens, not genomic resources",
            )


if __name__ == "__main__":
    unittest.main()
