from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import pandas as pd
from Bio import Phylo

ROOT = Path(__file__).resolve().parents[1]
COMPARISON_DIR = ROOT / "analyses" / "01_transposable_elements" / "bonsai_comparison"


def load_module(name: str, path: Path):  # type: ignore[no-untyped-def]
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


COMPARISON = load_module(
    "compare_bonsai_to_pca", COMPARISON_DIR / "compare_bonsai_to_pca.py"
)
VALIDATION = load_module(
    "validate_strategy_localization",
    COMPARISON_DIR / "validate_strategy_localization.py",
)


def test_fixed_tree_localizes_reproductive_strategy_labels() -> None:
    species = COMPARISON.canonical_species()
    tree = Phylo.read(VALIDATION.REFERENCE_TREE, "newick")
    codes, _ = VALIDATION.load_strategy_codes(species)
    metrics = VALIDATION.LocalizationContext(tree, species, codes).metrics(codes)

    assert metrics["global_fitch_parsimony_steps"] == 2
    assert metrics["global_one_nearest_neighbor_agreement"] == 33 / 34
    assert metrics["aquatic_eggs_smallest_containing_side_tips"] == 7
    assert metrics["direct_development_smallest_containing_side_tips"] == 3
    assert metrics[
        "terrestrial_eggs_aquatic_larvae_smallest_containing_side_tips"
    ] == 27


def test_tracked_strategy_validation_preserves_claim_boundaries() -> None:
    data_dir = VALIDATION.TRACKED_DATA
    tests = pd.read_csv(data_dir / VALIDATION.PERMUTATION_FILENAME).set_index("test_id")
    stability = pd.read_csv(data_dir / VALIDATION.BOOTSTRAP_STRATEGY_FILENAME).set_index(
        "strategy"
    )
    index = pd.read_csv(data_dir / VALIDATION.BOOTSTRAP_INDEX_FILENAME)
    manifest = json.loads((data_dir / VALIDATION.VALIDATION_MANIFEST_FILENAME).read_text())

    assert tests.loc["global_fitch_parsimony_steps", "p_value"] == 1e-5
    assert tests.loc[
        "global_fitch_parsimony_steps", "p_value_holm_all_nine_tests"
    ] == 9e-5
    assert stability.loc["aquatic_eggs", "bootstrap_edge_delimited_count"] == 55
    assert stability.loc[
        "direct_development", "bootstrap_edge_delimited_count"
    ] == 70
    assert stability.loc[
        "terrestrial_eggs_aquatic_larvae", "bootstrap_edge_delimited_count"
    ] == 2
    assert len(index) == 201
    assert set(index["replicate_id"]) == set(range(201))
    assert manifest["analysis_status"].startswith("exploratory_post_hoc")
    assert manifest["label_permutation"]["assigned_tips"] == 34
    assert manifest["label_permutation"]["unassigned_tips"] == []
    boundaries = " ".join(manifest["interpretation_boundaries"])
    assert "not causation" in boundaries
    assert "not a biological-specimen bootstrap" in boundaries
    assert "No per-species feature-error matrix" in boundaries
