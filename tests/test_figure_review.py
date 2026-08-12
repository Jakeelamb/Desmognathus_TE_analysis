from __future__ import annotations

import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "Publication/Desmognathus_figure_review.Rmd"
RENDERER = ROOT / "scripts/publication/render_figure_review.R"
MANIFEST = ROOT / "Publication/figures/FIGURE_MANIFEST.csv"
OUTPUT = ROOT / "Publication/Desmognathus_figure_review.html"
AUXILIARY_PNG = ROOT / "Publication/figures/Path_DAG_equivalence_class_reference.png"


def read_manifest() -> list[dict[str, str]]:
    with MANIFEST.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def figure_one_has_phylogeny(rows: list[dict[str, str]]) -> bool:
    figure_one = [row for row in rows if row["figure_id"] == "Figure_1"]
    if len(figure_one) != 1:
        return False
    evidence = " ".join(
        figure_one[0].get(field, "")
        for field in ("datasets", "legend", "alt_text", "basename")
    ).lower()
    return any(
        token in evidence
        for token in ("s38", "s39", "phylogen", "dated tree", "dated-tree", "time-calibrated")
    )


def test_review_source_is_self_contained_and_manifest_driven() -> None:
    assert SOURCE.exists()
    source = SOURCE.read_text(encoding="utf-8")

    assert "self_contained: true" in source
    assert '"FIGURE_MANIFEST.csv"' in source
    assert "for (index in seq_len(nrow(figure_manifest)))" in source
    assert "figure_manifest$width_mm / maximum_width_mm" in source
    assert "row$legend" in source
    assert "row$alt_text" in source
    assert "row$png" in source
    assert "AUXILIARY_DAG_BASENAME <- \"Path_DAG_equivalence_class_reference.png\"" in source
    assert "write_csv(" not in source
    assert "ggsave(" not in source


def test_review_exposes_critique_and_tree_boundaries() -> None:
    source = SOURCE.read_text(encoding="utf-8")

    for phrase in [
        "Narrative and material",
        "Cross-figure coherence",
        "Typography",
        "Color",
        "Spacing",
        "Legend and alt text",
        "Phylogeny",
        "Unresolved tree provenance",
        "publication or archive identifier",
        "calibration method",
        "branch-length provenance",
        "Figure 1 includes the Path24 dated phylogeny",
        "Auxiliary Path DAG reference",
        "Not a numbered manuscript figure.",
        "10 testable equivalence classes",
        "saturated 6-member class is not shown",
        "remains unscored because it has no d-separation claim",
        "current Path24 relative IOD-only analysis",
    ]:
        assert phrase in source
    assert "figure_1_has_phylogeny" in source


def test_manifest_assets_are_complete_for_the_gallery() -> None:
    rows = read_manifest()
    assert rows
    assert len({row["figure_id"] for row in rows}) == len(rows)

    for row in rows:
        assert float(row["width_mm"]) > 0
        assert float(row["height_mm"]) > 0
        assert row["legend"].strip()
        assert row["alt_text"].strip()
        assert (ROOT / "Publication/figures" / row["png"]).is_file()


def test_renderer_targets_and_validates_the_review_html() -> None:
    assert RENDERER.exists()
    renderer = RENDERER.read_text(encoding="utf-8")

    assert 'output_file <- "Desmognathus_figure_review.html"' in renderer
    assert "knit_root_dir = root" in renderer
    assert 'count_fixed("data:image/png;base64", html)' in renderer
    assert "expected_embedded_png_count <- nrow(manifest) + 1L" in renderer
    assert 'stop("Missing auxiliary Path24 DAG reference: ", auxiliary_png)' in renderer
    assert "embedded_png_count != expected_embedded_png_count" in renderer
    assert "Unresolved tree provenance" in renderer
    assert "Figure 1 includes the Path24 dated phylogeny" in renderer
    assert "Auxiliary Path DAG reference" in renderer
    assert "Not a numbered manuscript figure." in renderer
    assert "10 testable equivalence classes" in renderer
    assert "saturated 6-member class is not shown" in renderer
    assert "remains unscored because it has no d-separation claim" in renderer
    assert "current Path24 relative IOD-only analysis" in renderer


def test_rendered_review_if_present() -> None:
    if not OUTPUT.exists() or not AUXILIARY_PNG.exists():
        return

    rows = read_manifest()
    rendered = OUTPUT.read_text(encoding="utf-8")
    if "Auxiliary Path DAG reference" not in rendered:
        return

    assert rendered.count("data:image/png;base64") == len(rows) + 1
    assert 'src="figures/' not in rendered
    assert "Unresolved tree provenance" in rendered
    for phrase in [
        "Auxiliary Path DAG reference",
        "Not a numbered manuscript figure.",
        "10 testable equivalence classes",
        "saturated 6-member class is not shown",
        "remains unscored because it has no d-separation claim",
        "current Path24 relative IOD-only analysis",
    ]:
        assert phrase in rendered

    positions = [rendered.index(f'data-figure-id="{row["figure_id"]}"') for row in rows]
    assert positions == sorted(positions)
    if figure_one_has_phylogeny(rows):
        assert "Figure 1 includes the Path24 dated phylogeny" in rendered
    assert "data-figure-id=\"Auxiliary_Path_DAG_reference\"" in rendered
