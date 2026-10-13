from __future__ import annotations

import re
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
    assert "aes(x = estimate, y = y)" in block
    assert "ltr_tip_y_order" in block
    assert "ltr_canonical_order" in block
    assert "!identical(ltr_tip_y_order, ltr_canonical_order)" in block
    assert "ltr_tree_left + ltr_distribution" in block
    assert 'tag = "a"' in block
    assert 'tag = "b"' in block
    assert "GBE_DOUBLE_COLUMN_MM, 190" in block
    assert "geom_col" not in block
    assert "geom_violin" in block
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


def test_figure4_root_release_manifest_preserves_its_original_marks() -> None:
    manifest = pd.read_csv(FIGURE_MANIFEST).set_index("figure_id")
    row = manifest.loc["Figure_4"]

    assert row["release_status"] == "sensitivity_only_not_ectopic_rate"
    assert "element-level observations" in row["legend"]
    assert "species medians" in row["legend"]
    assert "horizontal blue intervals" in row["legend"]
    assert "element-bootstrap intervals" in row["legend"]
    assert "plotted tip-y order" in row["legend"]
    assert "equal terminal and internal mean depth" in row["legend"]
    assert "same 30 phylogenetic rows" in row["alt_text"]
    assert "380 hollow element observations" in row["alt_text"]
    assert "horizontal bootstrap intervals" in row["alt_text"]
    assert "dotted vertical line" in row["alt_text"]


def test_figure4_candidate_export_text_describes_violin_and_median_bar_marks() -> None:
    source = FIGURE_SOURCE.read_text()
    block = source.split("## ---- figure-4", maxsplit=1)[1].split(
        "## ---- figure-5", maxsplit=1
    )[0]
    panel_source = (ROOT / "scripts/publication/figure4_panels.R").read_text()

    # Inspect the text supplied to the candidate exporters. The root release is
    # deliberately frozen, and existing candidate files may precede rendering.
    for export_source in (block, panel_source):
        descriptions = [
            text.lower()
            for text in re.findall(r'"([^"\n]+)"', export_source)
            if "380" in text and "bootstrap" in text
        ]
        assert len(descriptions) == 2
        for description in descriptions:
            assert "violin" in description
            assert "median" in description
            assert "bar" in description
            assert "hollow" in description
            assert "bootstrap" in description
            assert "dotted" in description
            assert "shading" not in description
            assert "estimated genome" in description
            assert "n/a" in description
        legend = next(text for text in descriptions if "2,000" in text)
        assert "95%" in legend
        assert "bandwidth" in legend
        assert "not biological replicates" in legend
        assert "0.978 gbp per pg" in legend
        assert "without imputation" in legend
        assert "not independently measured species c-values" in legend
