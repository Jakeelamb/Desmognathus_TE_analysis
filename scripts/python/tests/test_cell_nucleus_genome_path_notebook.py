from __future__ import annotations

import json
import sys
import unittest
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "path_analysis" / "scripts"))
import build_cell_nucleus_genome_path_notebook as report  # noqa: E402
import build_cell_nucleus_genome_path_presentation as presentation  # noqa: E402


class CellNucleusGenomePathNotebookTests(unittest.TestCase):
    def test_candidate_gallery_uses_biologically_retained_display_contract(self) -> None:
        rows = presentation.candidate_gallery_rows(report.enumerate_three_node_dags())

        self.assertEqual(len(rows), 10)
        self.assertNotIn("genome_cell_only", set(rows["equivalence_class"]))
        self.assertNotIn("saturated", set(rows["equivalence_class"]))
        self.assertIn(
            report.USER_MECHANISMS["cell_to_nucleus_to_genome"],
            set(rows["dag_id"]),
        )
        self.assertIn("Cell → nucleus → genome", set(rows["gallery_title"]))

    def test_all_three_node_dags_collapse_to_expected_equivalence_classes(self) -> None:
        dags = report.enumerate_three_node_dags()

        self.assertEqual(len(dags), 25)
        self.assertEqual(dags["dag_id"].nunique(), 25)
        self.assertEqual(dags["equivalence_class"].nunique(), 11)
        self.assertEqual(
            int(dags.groupby("equivalence_class")["class_member_count"].first().sum()),
            25,
        )
        saturated = dags.loc[dags["equivalence_class"].eq("saturated")]
        self.assertEqual(len(saturated), 6)
        self.assertFalse(saturated["testable_by_dsep"].any())
        self.assertEqual(dags.loc[~dags["testable_by_dsep"]].shape[0], 6)

    def test_user_hypotheses_are_same_markov_equivalence_class(self) -> None:
        dags = report.enumerate_three_node_dags().set_index("dag_id")
        observed_classes = {
            dags.loc[dag_id, "equivalence_class"]
            for dag_id in report.USER_MECHANISMS.values()
        }

        self.assertEqual(observed_classes, {"nucleus_bridge"})
        self.assertEqual(
            dags.loc[list(report.USER_MECHANISMS.values()), "basis_claim"].unique().tolist(),
            ["genome_size _||_ cell_size | nucleus_size"],
        )

    def test_primary_traits_use_all_finalized_species(self) -> None:
        traits = report.load_primary_traits()

        self.assertEqual(len(traits), 24)
        self.assertEqual(traits["species"].nunique(), 24)
        self.assertEqual(
            traits.columns.tolist(),
            [
                "species",
                "tree_tip",
                "genome_size_pg",
                "nucleus_area_um2",
                "cell_area_um2",
                "n_genome_nuclei",
                "n_genome_images",
                "n_size_cells",
                "n_size_images",
            ],
        )
        self.assertTrue(
            (traits[["genome_size_pg", "nucleus_area_um2", "cell_area_um2"]] > 0)
            .all()
            .all()
        )
        self.assertFalse(traits.isna().any().any())
        self.assertIn("D. ochrophaeus", set(traits["species"]))
        self.assertIn("D. aeneus", set(traits["species"]))
        self.assertIn("D. wrighti", set(traits["species"]))
        self.assertIn("D. orestes", set(traits["species"]))
        fuscus = traits.loc[traits["species"].eq("D. fuscus")].iloc[0]
        self.assertAlmostEqual(float(fuscus["genome_size_pg"]), 16.36)
        self.assertTrue(pd.api.types.is_integer_dtype(traits["n_size_cells"]))
        self.assertEqual(int(traits["n_size_cells"].min()), 8)

    def test_measurement_bootstrap_preserves_complete_species_panel(self) -> None:
        draws = report.build_measurement_bootstrap_traits(
            n_bootstrap=12, seed=20260710
        )

        self.assertEqual(len(draws), 12 * 24)
        self.assertEqual(
            draws.groupby("bootstrap_replicate")["species"].nunique().unique().tolist(),
            [24],
        )
        self.assertEqual(
            draws.columns.tolist(),
            [
                "bootstrap_replicate",
                "species",
                "tree_tip",
                "genome_size_pg",
                "nucleus_area_um2",
                "cell_area_um2",
            ],
        )
        self.assertFalse(draws.isna().any().any())
        self.assertTrue(
            (draws[["genome_size_pg", "nucleus_area_um2", "cell_area_um2"]] > 0)
            .all()
            .all()
        )
        fuscus = draws.loc[draws["species"].eq("D. fuscus"), "genome_size_pg"]
        self.assertTrue(fuscus.eq(16.36).all())

    def test_full_phylogenetic_path_release_passes_all_replication_gates(self) -> None:
        output = ROOT / "results/data/research_review/cell_nucleus_genome_path"
        manifest = json.loads((output / "analysis_manifest.json").read_text())
        primary = pd.read_csv(output / "primary_model_ranking.csv")
        mechanisms = pd.read_csv(output / "user_mechanism_equivalence_ranking.csv")
        measurement = pd.read_csv(output / "measurement_bootstrap_rankings.csv.gz")
        failures = pd.read_csv(output / "analysis_failures.csv")

        self.assertEqual(manifest["run_scope"], "full_release")
        self.assertEqual(manifest["primary_species"], 24)
        self.assertEqual(manifest["measurement_bootstrap_replicates"], 250)
        self.assertEqual(manifest["published_bootstrap_trees_analyzed"], 200)
        self.assertEqual(manifest["simulation_replicates_per_class_and_regime"], 100)
        self.assertEqual(manifest["failure_count"], 0)
        self.assertFalse(manifest["causal_direction_identified"])
        self.assertTrue(manifest["user_mechanisms_markov_equivalent"])
        self.assertTrue(failures.empty)

        primary = primary.sort_values("rank")
        self.assertEqual(
            primary.iloc[:2]["model"].tolist(),
            ["cell_bridge", "cell_collider"],
        )
        self.assertTrue(bool(primary.iloc[0]["admissible_competitive"]))
        self.assertFalse(bool(primary.iloc[1]["admissible_competitive"]))
        self.assertGreater(float(primary.iloc[1]["delta_CICc"]), 2)
        self.assertEqual(len(mechanisms), 3)
        self.assertLess(float(mechanisms["CICc"].max() - mechanisms["CICc"].min()), 1e-12)
        self.assertLess(float(mechanisms["p"].max() - mechanisms["p"].min()), 1e-12)

        winners = measurement.loc[measurement["rank"].eq(1), "model"].value_counts()
        self.assertEqual(int(winners.sum()), 250)
        self.assertEqual(int(winners["cell_bridge"]), 246)
        self.assertEqual(int(winners["cell_collider"]), 3)
        self.assertEqual(int(winners["genome_bridge"]), 1)

    def test_primary_ranking_footer_is_derived_from_current_global_fit(self) -> None:
        ranking = pd.read_csv(
            ROOT
            / "results/data/research_review/cell_nucleus_genome_path/primary_model_ranking.csv"
        )

        footer = presentation.primary_ranking_footer(ranking)

        self.assertIn(
            "Cell bridge is the sole globally supported class within ΔCICc ≤ 2.",
            footer,
        )
        self.assertIn("All other testable classes fail global fit.", footer)
        self.assertNotIn(
            "nucleus bridge and cell collider pass global fit",
            footer.lower(),
        )

    def test_single_canonical_notebook_is_executed_and_rendered(self) -> None:
        output = ROOT / "notebooks/research_review/cell_nucleus_genome_path_analysis"
        notebook_path = ROOT / "notebooks/research_review/07_cell_nucleus_genome_path_analysis.ipynb"
        html_path = output / "index.html"
        presentation_manifest = json.loads(
            (output / "presentation_manifest.json").read_text()
        )
        notebook = json.loads(notebook_path.read_text())
        code_cells = [
            cell for cell in notebook["cells"] if cell["cell_type"] == "code"
        ]
        errors = [
            item
            for cell in code_cells
            for item in cell.get("outputs", [])
            if item.get("output_type") == "error"
        ]
        source = "\n".join(
            "".join(cell["source"]) for cell in notebook["cells"]
        )
        html = html_path.read_text(encoding="utf-8")

        self.assertGreaterEqual(len(code_cells), 10)
        self.assertTrue(
            all(cell.get("execution_count") is not None for cell in code_cells)
        )
        self.assertEqual(errors, [])
        self.assertIn("calibrated genome size in pg", source)
        self.assertIn("Markov-equivalent", source)
        self.assertIn(
            "resolves association and bridge structure more strongly than",
            source,
        )
        self.assertIn("data:image/png;base64", html)
        self.assertIn("98.4%", html)
        self.assertIn("Cell bridge 200/200", html)
        self.assertEqual(presentation_manifest["status"], "pass")
        self.assertTrue(
            presentation_manifest["gates"]["all_code_cells_executed"]
        )
        self.assertTrue(presentation_manifest["gates"]["no_error_outputs"])
        self.assertFalse(
            presentation_manifest["gates"]["causal_direction_identified"]
        )


if __name__ == "__main__":
    unittest.main()
