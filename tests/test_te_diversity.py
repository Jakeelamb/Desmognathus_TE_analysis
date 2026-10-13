from __future__ import annotations

import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
REPLAY_PATH = ROOT / "analyses/01_transposable_elements/recompute_pca.py"
SPEC = importlib.util.spec_from_file_location("te_release_replay", REPLAY_PATH)
if SPEC is None or SPEC.loader is None:
    raise RuntimeError(f"Cannot import TE release replay from {REPLAY_PATH}")
REPLAY = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(REPLAY)


def test_compact_inputs_reproduce_canonical_diversity() -> None:
    products = REPLAY.reconstruct_release()

    diversity = products["diversity"]
    assert list(diversity.columns) == [
        "species",
        "n_superfamilies_included",
        "shannon_entropy",
        "gini_simpson",
    ]
    panel_species = REPLAY.read_panel()
    assert diversity.shape == (34, 4)
    assert diversity["species"].tolist() == panel_species
    assert diversity["n_superfamilies_included"].eq(10).all()
    REPLAY.check_release(products)


def test_te_genome_size_overlap_uses_exact_canonical_diversity_rows() -> None:
    overlap = pd.read_csv(
        ROOT / "analyses/03_phylogenetic_path/data/te_genome_size_overlap.csv",
        float_precision="round_trip",
    )
    te34 = pd.read_csv(ROOT / "data/identity/te34_panel.csv")
    path24 = pd.read_csv(ROOT / "data/identity/path24_panel.csv")
    diversity = pd.read_csv(
        ROOT / "analyses/01_transposable_elements/data/te_diversity.csv",
        float_precision="round_trip",
    )

    assert not any(
        unsupported in column.casefold()
        for column in overlap.columns
        for unsupported in ("hill", "pielou")
    )

    te34_species = set(te34.loc[te34["te_resource_panel34_v1"].eq(True), "species"])
    path24_species = set(path24.loc[path24["path_analysis_panel24_v1"].eq(True), "species"])
    overlap21_species = te34_species & path24_species

    assert len(te34_species) == 34
    assert len(path24_species) == 24
    assert len(overlap21_species) == 21
    assert overlap["species"].is_unique
    assert set(overlap["species"]) == overlap21_species
    assert "genome_size_pg_fuscus_anchored" in overlap.columns
    assert "relative_iod_index" not in overlap.columns

    assert list(diversity.columns) == [
        "species",
        "n_superfamilies_included",
        "shannon_entropy",
        "gini_simpson",
    ]
    assert diversity.shape == (34, 4)
    assert diversity["species"].is_unique
    assert diversity["n_superfamilies_included"].eq(10).all()

    metric_columns = [
        "n_superfamilies_included",
        "shannon_entropy",
        "gini_simpson",
    ]
    observed = overlap[["species", *metric_columns]].sort_values("species")
    expected = diversity.loc[
        diversity["species"].isin(overlap21_species),
        ["species", *metric_columns],
    ].sort_values("species")
    pd.testing.assert_frame_equal(
        observed.reset_index(drop=True),
        expected.reset_index(drop=True),
        check_exact=True,
    )


def test_top_superfamily_selection_recloses_and_breaks_cutoff_ties_by_name() -> None:
    proportions = pd.Series(
        {
            **{f"feature_{number:02d}": 0.1 for number in range(9, 0, -1)},
            "Zulu": 0.05,
            "Alpha": 0.05,
        }
    )

    selected = REPLAY.select_top_superfamilies(proportions)

    assert selected.index.tolist() == [
        *[f"feature_{number:02d}" for number in range(1, 10)],
        "Alpha",
    ]
    assert "Zulu" not in selected.index
    assert len(selected) == 10
    assert np.isclose(selected.sum(), 1, rtol=0, atol=1e-15)
    assert np.isclose(selected.loc["Alpha"], 0.05 / 0.95, rtol=0, atol=1e-15)


def test_top_superfamily_selection_requires_ten_positive_categories() -> None:
    proportions = pd.Series({f"feature_{number}": 0.1 for number in range(9)})

    with np.testing.assert_raises_regex(
        AssertionError,
        "at least 10 positive superfamilies",
    ):
        REPLAY.select_top_superfamilies(proportions)


def test_diversity_metrics_use_natural_log_shannon_and_gini_simpson() -> None:
    probabilities = np.repeat(0.1, 10)
    metrics = REPLAY.diversity_metrics(probabilities)

    assert np.isclose(metrics["shannon_entropy"], np.log(10), rtol=0, atol=1e-15)
    assert np.isclose(metrics["gini_simpson"], 0.9, rtol=0, atol=1e-15)
    assert set(metrics) == {
        "shannon_entropy",
        "gini_simpson",
    }


def test_mass_accounting_is_qc_only_and_closes_exactly() -> None:
    panel_species = REPLAY.read_panel()
    mass = REPLAY.read_mass_accounting(panel_species)

    assert mass["species"].tolist() == panel_species
    for level in ("order", "superfamily"):
        retained_fraction = mass[f"{level}_retained_fraction"].to_numpy(dtype=float)
        unresolved_fraction = mass[f"{level}_unresolved_fraction"].to_numpy(dtype=float)
        retained_bases = mass[f"{level}_retained_aligned_bases"].to_numpy(dtype=float)
        unresolved_bases = mass[f"{level}_unresolved_aligned_bases"].to_numpy(dtype=float)
        total_bases = mass["total_aligned_bases"].to_numpy(dtype=float)

        np.testing.assert_allclose(retained_fraction + unresolved_fraction, 1, rtol=0, atol=5e-12)
        np.testing.assert_allclose(
            retained_bases + unresolved_bases,
            total_bases,
            rtol=1e-12,
            atol=1e-3,
        )

    compositions = REPLAY.read_compositions(panel_species)
    diversity = REPLAY.build_diversity(compositions["superfamily"])
    assert diversity.shape == (34, 4)
    assert "composition_mode" not in diversity.columns
    assert not any("unresolved" in column.casefold() for column in diversity.columns)
    assert not hasattr(REPLAY, "add_unresolved_accounting_bin")
    assert not hasattr(REPLAY, "UNRESOLVED_ACCOUNTING_BIN")
