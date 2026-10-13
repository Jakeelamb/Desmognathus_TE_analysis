from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
BUILDER_PATH = ROOT / "analyses/03_phylogenetic_path/build_path_release.py"
SPEC = importlib.util.spec_from_file_location("path_genome_release", BUILDER_PATH)
if SPEC is None or SPEC.loader is None:
    raise RuntimeError(f"Cannot import Path24 release builder from {BUILDER_PATH}")
BUILDER = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(BUILDER)

DATA = ROOT / "analyses/03_phylogenetic_path/data"
MORPHOLOGY = ROOT / "analyses/02_morphology"


def test_builder_reconstructs_canonical_path_inputs_and_validates_results() -> None:
    products = BUILDER.reconstruct_release()
    assert set(products) == {
        "candidate_dags.csv",
        "path24_traits.csv",
        "measurement_bootstrap_traits.csv",
        "te_genome_size_overlap.csv",
    }
    BUILDER.check_release(products)


def test_path24_uses_exact_conditional_genome_size_fields_from_s22() -> None:
    traits = pd.read_csv(DATA / "path24_traits.csv", float_precision="round_trip")
    source = pd.read_csv(
        MORPHOLOGY / "data/relative_nuclear_iod_species.csv",
        float_precision="round_trip",
    )
    assert list(traits.columns) == BUILDER.PATH24_SCHEMA
    assert len(traits) == 24
    assert traits["species"].is_unique

    genome_columns = [
        "genome_size_pg_fuscus_anchored",
        "genome_size_pg_ci_low",
        "genome_size_pg_ci_high",
    ]
    observed = traits[["species", *genome_columns]].sort_values("species")
    expected = source[["species", *genome_columns]].sort_values("species")
    pd.testing.assert_frame_equal(
        observed.reset_index(drop=True),
        expected.reset_index(drop=True),
        check_exact=True,
    )
    assert set(traits["genome_size_definition"]) == {BUILDER.GENOME_SIZE_DEFINITION}


def test_measurement_bootstrap_uses_exact_current_fuscus_anchor() -> None:
    draws = pd.read_csv(
        DATA / "measurement_bootstrap_traits.csv", float_precision="round_trip"
    )
    manifest = json.loads(
        (MORPHOLOGY / "provenance/iod_source_manifest.json").read_text(
            encoding="utf-8"
        )
    )
    anchor = manifest["conditional_genome_size_anchor"]
    anchor_pg = float(anchor["derived_anchor_pg_1c"])

    assert draws.shape == (6000, 6)
    assert draws.groupby("bootstrap_replicate")["species"].nunique().eq(24).all()
    assert float(anchor["reported_assembly_span_gbp"]) == 16.1
    assert float(anchor["gbp_per_pg"]) == 0.978
    fuscus = draws.loc[
        draws["species"].eq("D. fuscus"),
        "genome_size_pg_fuscus_anchored",
    ]
    np.testing.assert_array_equal(fuscus.to_numpy(), np.full(250, anchor_pg))


def test_all_path_result_identifiers_are_genome_size_native() -> None:
    comparison = pd.read_csv(DATA / "path_model_comparison.csv")
    regressions = pd.read_csv(DATA / "pairwise_phylogenetic_regressions.csv")
    edges = pd.read_csv(DATA / "path_edge_estimates.csv")
    dsep = pd.read_csv(DATA / "path_dsep_tests.csv")

    assert set(comparison["model"]) == BUILDER.MODEL_IDS
    assert set(regressions["relation_id"]) == BUILDER.RELATION_IDS
    assert set(edges["mechanism"]) == BUILDER.MECHANISM_IDS
    assert set(dsep["candidate_model"]) == BUILDER.MODEL_IDS
    nodes = set(BUILDER.NODES)
    assert set(regressions["outcome"]) <= nodes
    assert set(regressions["predictor"]) <= nodes
    assert set(edges["parent"]) <= nodes
    assert set(edges["child"]) <= nodes

    primary = comparison.sort_values("rank")
    assert primary.iloc[:3]["model"].tolist() == [
        "cell_bridge",
        "cell_collider",
        "genome_bridge",
    ]
    assert np.isclose(primary.iloc[0]["CICc"], 13.382923669911122)


def test_te_overlap_uses_exact_s22_genome_size_rows() -> None:
    overlap = pd.read_csv(
        DATA / "te_genome_size_overlap.csv", float_precision="round_trip"
    )
    source = pd.read_csv(
        MORPHOLOGY / "data/relative_nuclear_iod_species.csv",
        float_precision="round_trip",
    ).assign(species=lambda frame: frame["species"].str.removeprefix("D. "))
    genome_columns = [
        "genome_size_pg_fuscus_anchored",
        "genome_size_pg_ci_low",
        "genome_size_pg_ci_high",
    ]
    expected = source.loc[source["species"].isin(overlap["species"])]
    observed = overlap[["species", *genome_columns]].sort_values("species")
    expected = expected[["species", *genome_columns]].sort_values("species")
    pd.testing.assert_frame_equal(
        observed.reset_index(drop=True),
        expected.reset_index(drop=True),
        check_exact=True,
    )
    assert set(overlap["n_superfamilies_included"]) == {10}
    assert "observed_richness" not in overlap.columns
    assert set(overlap["te_diversity_definition"]) == {
        BUILDER.TE_DIVERSITY_DEFINITION
    }
    assert "ten most abundant" in BUILDER.TE_DIVERSITY_DEFINITION
    assert "reclosed" in BUILDER.TE_DIVERSITY_DEFINITION


def test_r_runner_has_no_relative_iod_analysis_interface() -> None:
    source = (
        ROOT / "analyses/03_phylogenetic_path/run_path_analysis.R"
    ).read_text(encoding="utf-8")
    assert "genome_size_pg_fuscus_anchored" in source
    assert 'node_labels <- c(g = "genome_size"' in source
    assert "genome_nucleus_only" in source
    assert "genome_bridge" in source
    assert "relative_iod_index" not in source
    assert "iod_nucleus_only" not in source
    assert "iod_bridge" not in source
    assert "16.36" not in source
    assert "16.1-Gbp" in source
    assert "16.462167689 pg/1C" in source
