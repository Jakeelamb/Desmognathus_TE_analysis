from __future__ import annotations

import csv
import json
import re
from pathlib import Path

from test_auxiliary_figure_contract import (
    assert_report_lineage,
    normalize_report_text,
    rendered_report_text,
)

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "Publication/Desmognathus_figure_review.Rmd"
RENDERER = ROOT / "scripts/publication/render_figure_review.R"
MANIFEST = ROOT / "Publication/figures/FIGURE_MANIFEST.csv"
AUXILIARY_MANIFEST = ROOT / "Publication/figures/AUXILIARY_FIGURE_MANIFEST.csv"
CONTRACT = ROOT / "scripts/publication/auxiliary_figure_contract.R"
HISTORICAL_CHECKSUMS = (
    ROOT / "outputs/manuscript_story_20260924_01/qa/original_checksums.json"
)
OUTPUT = ROOT / "Publication/Desmognathus_figure_review.html"
EXPECTED_AUXILIARY_IDS = [
    "Auxiliary_Wang_style_integrated_te_history_candidate",
    "Auxiliary_Wang_style_diversity_genome_size_candidate",
    "Auxiliary_Path24_development_traits_reference",
    "Auxiliary_Path_DAG_reference",
    "Auxiliary_Wang2025_cross_study_te_diversity_reference",
]
EXPECTED_AUXILIARY_TITLES = [
    "Wang-style integrated TE content and proliferation candidate",
    "Wang-style TE diversity and genome-size candidate",
    "Path24 development and phenotype reference",
    "Path DAG equivalence-class reference",
    "Cross-study TE-diversity reference (Wang et al. 2025)",
]


CURRENT_AUXILIARY_TITLES = [*EXPECTED_AUXILIARY_TITLES]
CURRENT_AUXILIARY_TITLES[0] = "Wang-style integrated TE history across all 34 species"
CURRENT_AUXILIARY_TITLES[2] = "Path24 phylogeny and phenotype reference"


def read_manifest() -> list[dict[str, str]]:
    with MANIFEST.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def read_auxiliary_manifest() -> list[dict[str, str]]:
    with AUXILIARY_MANIFEST.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def figure_one_has_phylogeny(rows: list[dict[str, str]]) -> bool:
    figure_one = [row for row in rows if row["figure_id"] == "Figure_1"]
    if len(figure_one) != 1:
        return False
    evidence = " ".join(
        figure_one[0].get(field, "") for field in ("datasets", "legend", "alt_text", "basename")
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
    assert '"AUXILIARY_FIGURE_MANIFEST.csv"' in source
    assert "for (index in seq_len(nrow(figure_manifest)))" in source
    assert "for (index in seq_len(nrow(auxiliary_manifest)))" in source
    assert "figure_manifest$width_mm / maximum_width_mm" in source
    assert "auxiliary_manifest$width_mm / maximum_width_mm" in source
    assert "row$legend" in source
    assert "row$alt_text" in source
    assert "row$png" in source
    assert "validate_auxiliary_figure_manifest(auxiliary_manifest)" in source
    assert "write_csv(" not in source
    assert "ggsave(" not in source


def test_review_exposes_critique_and_tree_boundaries() -> None:
    source = SOURCE.read_text(encoding="utf-8")
    assert "Six representatives are selected" not in source
    assert "auxiliary_manifest$boundary[[1]]" in source
    assert "auxiliary_manifest$legend[[1]]" in source
    source += CONTRACT.read_text()

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
        "Figure 1 is the integrated 21-species overview",
        "Figure 2 Bonsai boundary",
        "TE-composition representation, not a species phylogeny",
        "feature resampling perturbs correlated compositional coordinates",
        "Wang-style integrated TE content and proliferation candidate",
        "Review candidate, not a numbered manuscript figure",
        "Six representatives are selected",
        "fractions of aligned RepeatMasker hit bp",
        "divergence bins are not insertion ages",
        "Wang-style TE diversity and genome-size candidate",
        "Reproductive-strategy panels are post-hoc descriptive strata",
        "without a regression because n = 2",
        "Pagel-lambda PGLS fits",
        "Path24 development and phenotype reference",
        "Descriptive review only",
        "D. organi is absent",
        "D. fuscus assembly-anchored pg/1C estimates",
        "not direct independently measured species C-values",
        "Path DAG equivalence-class reference",
        "Not a numbered manuscript figure.",
        "10 testable equivalence classes",
        "saturated 6-member class is not shown",
        "remains unscored because it has no d-separation claim",
        "Genome size",
        "D. fuscus assembly-anchored pg/1C analysis",
        "not a direct species C-value measurement",
        "same released genome-size scale is used consistently",
        "do not orient a unique causal arrow",
        "Cross-study TE-diversity reference (Wang et al. 2025)",
        "Descriptive cross-study context only",
        "no pooled regression or hypothesis test",
        "10.5061/dryad.zpc866tkv",
        "10.1038/s42003-025-08127-3",
        "CC0-1.0",
        "published top-10 TE-superfamily indices",
        "RepeatModeler/RepeatMasker",
        "reclosed calculations from dnaPipeTE compositions",
        "pipelines are different and not method-harmonized",
        "conditional D. fuscus assembly-anchored estimates",
        "not direct species C-values",
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

    genome_size_basenames = {
        "Figure_6": "Figure_6_genome_size_pairwise_pgls",
    }
    rows_by_id = {row["figure_id"]: row for row in rows}
    assert rows_by_id["Figure_2"]["basename"] == "Figure_2_te_diversity_and_bonsai"
    figure_two_copy = " ".join(
        rows_by_id["Figure_2"][field] for field in ("datasets", "legend", "alt_text")
    )
    for phrase in ("S42", "S43", "Bonsai", "not a species phylogeny", "feature resampling"):
        assert phrase.casefold() in figure_two_copy.casefold()
    for figure_id, basename in genome_size_basenames.items():
        assert rows_by_id[figure_id]["basename"] == basename
        public_copy = " ".join(
            rows_by_id[figure_id][field]
            for field in ("legend", "alt_text", "release_status")
        )
        assert "genome size" in public_copy.lower()
        assert not re.search(
            r"relative[- ]?(?:nuclear[- ]?)?iod", public_copy, re.IGNORECASE
        )

    auxiliary_rows = read_auxiliary_manifest()
    assert [row["auxiliary_id"] for row in auxiliary_rows] == EXPECTED_AUXILIARY_IDS
    assert [row["title"] for row in auxiliary_rows] in (
        EXPECTED_AUXILIARY_TITLES, CURRENT_AUXILIARY_TITLES
    )
    for row in auxiliary_rows:
        assert float(row["width_mm"]) > 0
        assert float(row["height_mm"]) > 0
        assert row["boundary"].strip()
        assert row["legend"].strip()
        assert row["alt_text"].strip()
        assert (ROOT / "Publication/figures" / row["png"]).is_file()


def test_renderer_targets_and_validates_the_review_html() -> None:
    assert RENDERER.exists()
    renderer = RENDERER.read_text(encoding="utf-8")
    assert "validate_auxiliary_figure_manifest(auxiliary_manifest)" in renderer
    assert "validate_auxiliary_rendered_html(auxiliary_contract, html)" in renderer
    assert "validate_report_source_fingerprint(source_path, html)" in renderer
    renderer += CONTRACT.read_text()

    assert 'output_file <- "Desmognathus_figure_review.html"' in renderer
    assert "knit_root_dir = root" in renderer
    assert 'count_fixed("data:image/png;base64", html)' in renderer
    assert "expected_embedded_png_count <- nrow(manifest) + nrow(auxiliary_manifest)" in renderer
    assert "expected_embedded_png_count != 12L" in renderer
    assert "nrow(manifest) != 7L" in renderer
    assert "five ordered IDs" in renderer
    assert 'stop("Missing auxiliary figure manifest: ", auxiliary_manifest_path)' in renderer
    assert "Missing auxiliary review PNG files" in renderer
    assert "embedded_png_count != expected_embedded_png_count" in renderer
    assert "Unresolved tree provenance" in renderer
    assert "Figure 1 is the integrated 21-species overview" in renderer
    assert "Figure_2_te_diversity_and_bonsai" in renderer
    assert "Figure_2 is missing its Bonsai boundary phrase" in renderer
    assert "Auxiliary_Wang_style_integrated_te_history_candidate" in renderer
    assert "Wang-style integrated TE content and proliferation candidate" in renderer
    assert "Six representatives are selected" in renderer
    assert "fractions of aligned RepeatMasker hit bp" in renderer
    assert "divergence bins are not insertion ages" in renderer
    assert "Auxiliary_Wang_style_diversity_genome_size_candidate" in renderer
    assert "Wang-style TE diversity and genome-size candidate" in renderer
    assert "Reproductive-strategy panels are post-hoc descriptive strata" in renderer
    assert "without a regression because n = 2" in renderer
    assert "Pagel-lambda PGLS fits" in renderer
    assert "Path24 development and phenotype reference" in renderer
    assert "Descriptive review only" in renderer
    assert "D. organi is absent" in renderer
    assert "D. fuscus assembly-anchored pg/1C estimates" in renderer
    assert "Path DAG equivalence-class reference" in renderer
    assert "Not a numbered manuscript figure." in renderer
    assert "10 testable equivalence classes" in renderer
    assert "saturated 6-member class is not shown" in renderer
    assert "remains unscored because it has no d-separation claim" in renderer
    assert "Genome size" in renderer
    assert "D. fuscus assembly-anchored pg/1C analysis" in renderer
    assert "not a direct species C-value measurement" in renderer
    assert "same released genome-size scale is used consistently" in renderer
    assert "do not orient a unique causal arrow" in renderer
    assert "Auxiliary_Wang2025_cross_study_te_diversity_reference" in renderer
    assert "Cross-study TE-diversity reference (Wang et al. 2025)" in renderer
    assert "Descriptive cross-study context only" in renderer
    assert "no pooled regression or hypothesis test" in renderer
    assert "10.5061/dryad.zpc866tkv" in renderer
    assert "10.1038/s42003-025-08127-3" in renderer
    assert "CC0-1.0" in renderer
    assert "published top-10 TE-superfamily indices" in renderer
    assert "RepeatModeler/RepeatMasker" in renderer
    assert "reclosed calculations from dnaPipeTE compositions" in renderer
    assert "pipelines are different and not method-harmonized" in renderer
    assert "conditional D. fuscus assembly-anchored estimates" in renderer
    assert "Figure_6_genome_size_pairwise_pgls" in renderer
    assert "Figure_S1" not in renderer
    assert "not direct species C-values" in renderer


def test_rendered_review_if_present() -> None:
    if not OUTPUT.exists() or not AUXILIARY_MANIFEST.exists():
        return

    rows = read_manifest()
    auxiliary_rows = read_auxiliary_manifest()
    if [row["auxiliary_id"] for row in auxiliary_rows] != EXPECTED_AUXILIARY_IDS:
        return
    if not all((ROOT / "Publication/figures" / row["png"]).exists() for row in auxiliary_rows):
        return
    rendered = OUTPUT.read_text(encoding="utf-8")
    assert_report_lineage(
        ROOT, SOURCE, OUTPUT, CONTRACT, json.loads(HISTORICAL_CHECKSUMS.read_text()),
        "a78148a279e4e77a4a409b45dccfe4a297f84896e1ee40dba8511347ba0bfa68",
        historical_source_record=False,
    )
    rendered_text = rendered_report_text(rendered)

    assert rendered.count("data:image/png;base64") == len(rows) + len(auxiliary_rows)
    assert 'src="figures/' not in rendered
    assert "Unresolved tree provenance" in rendered
    for row in auxiliary_rows:
        for field in ("title", "boundary", "legend"):
            assert normalize_report_text(row[field]) in rendered_text
    if auxiliary_rows[0]["title"] == CURRENT_AUXILIARY_TITLES[0]:
        assert "All 34 TE species are included" in rendered_text
        assert "Six representatives are selected" not in rendered_text
    else:
        assert "Six representatives are selected" in rendered_text

    positions = [rendered.index(f'data-figure-id="{row["figure_id"]}"') for row in rows]
    assert positions == sorted(positions)
    if figure_one_has_phylogeny(rows):
        assert "Figure 1 is the integrated 21-species overview" in rendered
    auxiliary_positions = [
        rendered.index(f'data-figure-id="{row["auxiliary_id"]}"') for row in auxiliary_rows
    ]
    assert auxiliary_positions == sorted(auxiliary_positions)
    assert positions[-1] < auxiliary_positions[0]
