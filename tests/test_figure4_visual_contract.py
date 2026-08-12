from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
FIGURE_SOURCE = ROOT / "scripts/publication/build_gbe_figures.R"
FIGURE_MANIFEST = ROOT / "Publication/figures/FIGURE_MANIFEST.csv"
LTR_ELEMENTS = (
    ROOT
    / "Publication/datasets/Supplementary_Data_S14_ltr_terminal_internal_elements.csv"
)
LTR_SUMMARIES = (
    ROOT
    / "Publication/datasets/Supplementary_Data_S15_ltr_terminal_internal_species_robustness.csv"
)


def test_figure4_uses_the_primary_element_level_release() -> None:
    source = FIGURE_SOURCE.read_text()
    block = source.split("## ---- figure-4", maxsplit=1)[1].split(
        "## ---- figure-5", maxsplit=1
    )[0]

    assert "Supplementary_Data_S14_ltr_terminal_internal_elements.csv" in block
    assert "filter(primary_ltr_element)" in block
    assert "log2_ratio_terminal_internal_all_positions" in block
    assert "nrow(ltr_elements) != 380L" in block
    assert 'n_distinct(ltr_elements$species) != 30L' in block
    assert "position_jitter" in block
    assert "geom_col" not in block
    assert "equal terminal and internal mean depth (1:1)" in block


def test_figure4_element_counts_and_summaries_are_exact() -> None:
    elements = pd.read_csv(LTR_ELEMENTS, low_memory=False)
    elements = elements.loc[elements["primary_ltr_element"]].copy()
    summaries = pd.read_csv(LTR_SUMMARIES)
    summaries = summaries.loc[
        summaries["analysis_branch"].eq(
            "primary_iqr_filtered_coverage_ge_80pct"
        )
    ].copy()

    assert (len(elements), elements["species"].nunique()) == (380, 30)
    assert (len(summaries), summaries["species"].nunique()) == (30, 30)

    observed = elements.groupby("species").agg(
        n_elements=("element_id", "size"),
        ratio_median=("ratio_terminal_internal_all_positions", "median"),
    )
    released = summaries.set_index("species")
    assert observed["n_elements"].equals(released["n_elements"])
    assert np.allclose(
        observed["ratio_median"],
        released["ratio_median"],
        rtol=0,
        atol=1e-12,
    )


def test_figure4_manifest_describes_the_visible_marks() -> None:
    manifest = pd.read_csv(FIGURE_MANIFEST).set_index("figure_id")
    row = manifest.loc["Figure_4"]

    assert row["release_status"] == "sensitivity_only_not_ectopic_rate"
    assert "element-level observations" in row["legend"]
    assert "species medians" in row["legend"]
    assert "equal terminal and internal mean depth" in row["legend"]
    assert "support count" in row["alt_text"]
