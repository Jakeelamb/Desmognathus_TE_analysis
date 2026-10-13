from __future__ import annotations

import csv
import json
from pathlib import Path

from test_auxiliary_figure_contract import (
    assert_report_lineage,
    normalize_report_text,
    rendered_report_text,
)

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "Publication/Desmognathus_study_data.Rmd"
OUTPUT = ROOT / "Publication/Desmognathus_study_data.html"
PATH_README = ROOT / "analyses/03_phylogenetic_path/README.md"
FIGURE_MANIFEST = ROOT / "Publication/figures/FIGURE_MANIFEST.csv"
AUXILIARY_MANIFEST = ROOT / "Publication/figures/AUXILIARY_FIGURE_MANIFEST.csv"
CONTRACT = ROOT / "scripts/publication/auxiliary_figure_contract.R"
HISTORICAL_CHECKSUMS = (
    ROOT / "outputs/manuscript_story_20260924_01/qa/original_checksums.json"
)

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


def read_csv_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def test_study_report_replaces_notebook_surface() -> None:
    assert SOURCE.exists()
    assert OUTPUT.exists()
    assert not list((ROOT / "analyses").glob("*/explore.ipynb"))
    assert not (ROOT / "scripts/build_notebooks.py").exists()

    project = (ROOT / "pyproject.toml").read_text()
    assert "jupyterlab" not in project
    assert "nbformat" not in project

    environment = (ROOT / "environment.yml").read_text()
    assert "r-rmarkdown=2.31" in environment

    assert not any(line.rstrip("\r\n") != line.rstrip() for line in OUTPUT.open())


def test_study_report_exposes_current_release_contract() -> None:
    source = SOURCE.read_text()
    rendered = OUTPUT.read_text()

    required_source_phrases = {
        "nrow(te_diversity) == 34L",
        "n_superfamilies_included == 10L",
        "slice_max(proportion, n = 10, with_ties = FALSE)",
        "nrow(manifest) == 43L",
        'sprintf("S%02d", seq_len(43))',
        "S01–S43",
        "Unresolved aligned-base",
        "S07",
        "primary_iqr_filtered_coverage_ge_80pct",
        "genome_size_pg_fuscus_anchored",
        "relative_iod_index",
        "0.0144",
        "Q1 - 1.5",
        "Gini-Simpson",
        'cross_study_diversity <- read_supp("S41")',
        'bonsai_tips <- read_supp("S42")',
        'bonsai_validation <- read_supp("S43")',
        "nrow(bonsai_tips) == 34L",
        "nrow(bonsai_validation) == 3L",
        "feature_bootstrap_edge_delimited_count",
        "not a species phylogeny",
        "not a biological-specimen bootstrap",
        'figure_path("Figure_2_te_diversity_and_bonsai.png")',
        'figure_path("Figure_1_integrated_TE_landscape.png")',
        'figure_path("Figure_5_phylogeny_and_cellular_traits.png")',
        "nrow(cross_study_diversity) == 118L",
        'c(84L, 34L)',
        'c(11L, 73L, 21L, 13L)',
        "scripts/publication/build_gbe_figures.R",
    }
    assert all(phrase in source for phrase in required_source_phrases)
    assert "mass_aware_unresolved_bin" not in source
    assert "mass_aware_unresolved_bin" not in rendered
    assert "composition_mode" not in source
    assert "write_csv(" not in source
    assert "write.csv(" not in source
    assert "ggsave(" not in source

    assert_report_lineage(
        ROOT, SOURCE, OUTPUT, CONTRACT, json.loads(HISTORICAL_CHECKSUMS.read_text()),
        "4f69183df61bd6b8c5c4ba760bdef04d81c55501e2e75bbb54100a6cfdda5f48",
        historical_source_record=True,
    )
    for expected in [
        "TE34",
        "Path24",
        "overlap21",
        "380",
        "1,152",
        "805",
        "Figure 2. Bonsai TE-composition structure and reproductive-strategy localization.",
        "Figure 4. LTR terminal:internal mapping/deletion-footprint proxy.",
        "Figure 5. Phylogenetically aligned conditional genome size and reviewed cellular traits.",
        "Figure 6. Conditional genome size, nucleus area, and cell area sensitivities.",
    ]:
        assert expected in rendered
    assert rendered.count("data:image/png;base64") == 12

    for error in ["Quitting from lines", "Execution halted", "Error in "]:
        assert error not in rendered


def test_study_report_auxiliary_manifest_contract() -> None:
    source = SOURCE.read_text()
    normalized_source = " ".join((source + CONTRACT.read_text()).split())
    renderer = (ROOT / "scripts/publication/render_study_report.R").read_text()
    assert "validate_auxiliary_figure_manifest(auxiliary_manifest)" in renderer
    assert "validate_auxiliary_rendered_html(auxiliary_contract, html)" in renderer
    assert "validate_report_source_fingerprint(source_path, html)" in renderer
    assert "Six representatives are selected" not in source
    assert "auxiliary_manifest$boundary[[1]]" in source
    assert "auxiliary_manifest$legend[[1]]" in source
    renderer += CONTRACT.read_text()

    assert AUXILIARY_MANIFEST.exists()
    assert len(read_csv_rows(FIGURE_MANIFEST)) == 7
    rows = read_csv_rows(AUXILIARY_MANIFEST)
    assert [row["auxiliary_id"] for row in rows] == EXPECTED_AUXILIARY_IDS
    assert [row["title"] for row in rows] in (
        EXPECTED_AUXILIARY_TITLES, CURRENT_AUXILIARY_TITLES
    )

    required_source_phrases = {
        '"AUXILIARY_FIGURE_MANIFEST.csv"',
        "for (index in seq_len(nrow(auxiliary_manifest)))",
        "Wang-style integrated TE-history candidate boundary",
        "Review candidate, not a numbered manuscript figure",
        "Six representatives are selected",
        "fractions of aligned RepeatMasker hit bp",
        "divergence bins are not insertion ages",
        "Wang-style TE-diversity candidate boundary",
        "Reproductive-strategy panels are post-hoc descriptive strata",
        "without a regression because n = 2",
        "Pagel-lambda PGLS fits",
        "Descriptive review only",
        "2 direct-development species and 22 aquatic-larval species",
        "D. organi is absent from Path24",
        "No development mode or phenotype is imputed",
        "conditional D. fuscus assembly-anchored estimates",
        "not direct independently measured species C-values",
        "tree provenance remains unresolved",
        "Genome size node",
        "conditional D. fuscus assembly-anchored pg/1C estimates",
        "same genome-size scale is used across all scored classes",
        "not direct species C-value measurements",
        "representative only and do not orient unique causal arrows",
        "Descriptive cross-study context only",
        "no pooled regression or hypothesis test",
        "10.5061/dryad.zpc866tkv",
        "10.1038/s42003-025-08127-3",
        "CC0-1.0",
        "published top-10 TE-superfamily indices",
        "RepeatModeler/RepeatMasker",
        "reclosed calculations from dnaPipeTE compositions",
        "pipelines are different and not method-harmonized",
        "not direct species C-values",
    }
    assert all(phrase in normalized_source for phrase in required_source_phrases)

    required_renderer_phrases = {
        "nrow(numbered_manifest) != 7L",
        "nrow(auxiliary_manifest) != 5L",
        "expected_embedded_png_count <- 12L",
        "Auxiliary_Wang_style_integrated_te_history_candidate",
        "Auxiliary_Wang_style_diversity_genome_size_candidate",
        "Wang-style integrated TE content and proliferation candidate",
        "Wang-style TE diversity and genome-size candidate",
        "Six representatives are selected",
        "fractions of aligned RepeatMasker hit bp",
        "divergence bins are not insertion ages",
        "Reproductive-strategy panels are post-hoc descriptive strata",
        "without a regression because n = 2",
        "Pagel-lambda PGLS fits",
        "Auxiliary_Path24_development_traits_reference",
        "Auxiliary_Path_DAG_reference",
        "Auxiliary_Wang2025_cross_study_te_diversity_reference",
        "Path24 development and phenotype reference",
        "Path DAG equivalence-class reference",
        "Cross-study TE-diversity reference (Wang et al. 2025)",
        "Figure_6_genome_size_pairwise_pgls",
        "Figure_2_te_diversity_and_bonsai",
        "not a species phylogeny",
        "feature resampling",
        "Genome size",
        "same released genome-size scale is used consistently",
        "10.5061/dryad.zpc866tkv",
        "10.1038/s42003-025-08127-3",
        "CC0-1.0",
        "no pooled regression or hypothesis test",
    }
    assert all(phrase in renderer for phrase in required_renderer_phrases)


def test_study_report_sources_preserve_dag_equivalence_boundary() -> None:
    source = SOURCE.read_text()
    normalized_source = " ".join(source.split())
    readme = PATH_README.read_text()

    required_source_phrases = {
        "Path DAG equivalence-class reference",
        "10 testable Markov-equivalence classes",
        "The saturated 6-member class is omitted",
        "Genome size node",
        "not direct species C-value measurements",
        "representative only and do not orient unique causal arrows",
    }
    assert all(phrase in normalized_source for phrase in required_source_phrases)
    assert "dag_reference <- tibble::tribble(" not in source
    assert 'title = "Auxiliary DAG reference for Figure 7"' not in source
    assert "Representative of 3 DAGs; ΔCICc = 0.000; Passed" not in source

    required_readme_phrases = {
        "25 possible three-node DAG members",
        "10 testable",
        "6-member saturated class",
        "one representative DAG per scored",
    }
    assert all(phrase in readme for phrase in required_readme_phrases)


def test_study_report_keeps_iod_inside_the_calibration_proof() -> None:
    source = SOURCE.read_text()
    proof_heading = "### Calibration proof: optical density to conditional genome size"
    assert source.count(proof_heading) == 1
    before, remainder = source.split(proof_heading, maxsplit=1)
    proof, after = remainder.split("## Integrated species views", maxsplit=1)

    assert "relative_iod_index" in proof
    assert "iod_equal_image_estimate" in proof
    assert "genome_size_pg_fuscus_anchored" in proof
    outside_proof = before + after
    assert "relative_iod_" not in outside_proof
    assert "relative nuclear IOD" not in outside_proof
    assert "relative-IOD" not in outside_proof

    for expected in [
        'figure_path("Figure_6_genome_size_pairwise_pgls.png")',
        "genome_size_vs_nucleus_area",
        "genome_size_vs_cell_area",
    ]:
        assert expected in source or expected in (
            ROOT / "scripts/publication/build_gbe_figures.R"
        ).read_text()
    assert "Figure_S1" not in source
    assert "Figure_S1" not in (
        ROOT / "scripts/publication/build_gbe_figures.R"
    ).read_text()
    for obsolete in [
        "Figure_6_relative_iod_pairwise_pgls",
        "Figure_S1_te_relative_iod_overlap21",
    ]:
        assert obsolete not in source


def test_rendered_study_report_has_ordered_auxiliary_figures_if_present() -> None:
    if not OUTPUT.exists() or not AUXILIARY_MANIFEST.exists():
        return

    auxiliary_rows = read_csv_rows(AUXILIARY_MANIFEST)
    if [row["auxiliary_id"] for row in auxiliary_rows] != EXPECTED_AUXILIARY_IDS:
        return
    if not all((ROOT / "Publication/figures" / row["png"]).exists() for row in auxiliary_rows):
        return

    rendered = OUTPUT.read_text(encoding="utf-8")
    if EXPECTED_AUXILIARY_IDS[0] not in rendered:
        return

    assert rendered.count("data:image/png;base64") == 12
    positions = [rendered.index(auxiliary_id) for auxiliary_id in EXPECTED_AUXILIARY_IDS]
    assert positions == sorted(positions)
    rendered_text = rendered_report_text(rendered)
    for row in auxiliary_rows:
        for field in ("title", "boundary", "legend"):
            assert normalize_report_text(row[field]) in rendered_text
    assert 'src="figures/' not in rendered
