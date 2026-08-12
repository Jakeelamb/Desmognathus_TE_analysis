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
        "te_level",
        "observed_richness",
        "shannon_entropy",
        "gini_simpson",
    ]
    assert diversity.groupby("te_level", sort=False).size().to_dict() == {
        "order": 34,
        "superfamily": 34,
    }
    panel_species = REPLAY.read_panel()
    assert diversity.loc[diversity["te_level"].eq("order"), "species"].tolist() == panel_species
    assert (
        diversity.loc[diversity["te_level"].eq("superfamily"), "species"].tolist() == panel_species
    )
    REPLAY.check_release(products)


def test_te_relative_iod_overlap_uses_exact_canonical_diversity_rows() -> None:
    overlap = pd.read_csv(ROOT / "analyses/03_phylogenetic_path/data/te_relative_iod_overlap.csv")
    te34 = pd.read_csv(ROOT / "data/identity/te34_panel.csv")
    path24 = pd.read_csv(ROOT / "data/identity/path24_panel.csv")
    diversity = pd.read_csv(ROOT / "analyses/01_transposable_elements/data/te_diversity.csv")

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

    source_rows = diversity.loc[diversity["te_level"].eq("superfamily")]
    assert len(source_rows) == 34
    assert source_rows["species"].is_unique

    metric_columns = [
        "observed_richness",
        "shannon_entropy",
        "gini_simpson",
    ]
    observed = overlap[["species", *metric_columns]].sort_values("species")
    expected = source_rows.loc[
        source_rows["species"].isin(overlap21_species),
        ["species", *metric_columns],
    ].sort_values("species")
    pd.testing.assert_frame_equal(
        observed.reset_index(drop=True),
        expected.reset_index(drop=True),
        check_exact=True,
    )


def test_diversity_metrics_use_exact_positive_categories() -> None:
    probabilities = np.array([0.5, 0.25, 0.25, 0.0])
    metrics = REPLAY.diversity_metrics(probabilities)
    expected_shannon = -np.sum(probabilities[:3] * np.log(probabilities[:3]))
    expected_sum_squared = np.sum(probabilities**2)

    assert metrics["observed_richness"] == 3
    assert np.isclose(metrics["shannon_entropy"], expected_shannon)
    assert np.isclose(metrics["gini_simpson"], 1 - expected_sum_squared)
    assert set(metrics) == {
        "observed_richness",
        "shannon_entropy",
        "gini_simpson",
    }

    singleton = REPLAY.diversity_metrics(np.array([1.0, 0.0]))
    assert singleton == {
        "observed_richness": 1,
        "shannon_entropy": 0.0,
        "gini_simpson": 0.0,
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

    diversity = REPLAY.build_diversity(REPLAY.read_compositions(panel_species))
    assert diversity.shape == (68, 5)
    assert "composition_mode" not in diversity.columns
    assert not any("unresolved" in column.casefold() for column in diversity.columns)
    assert not hasattr(REPLAY, "add_unresolved_accounting_bin")
    assert not hasattr(REPLAY, "UNRESOLVED_ACCOUNTING_BIN")
