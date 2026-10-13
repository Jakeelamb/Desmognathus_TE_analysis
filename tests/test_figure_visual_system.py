from __future__ import annotations

import csv
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
THEME_SOURCE = ROOT / "scripts/publication/gbe_theme.R"
FIGURE_SOURCE = ROOT / "scripts/publication/build_gbe_figures.R"
FIGURE_MANIFEST = ROOT / "Publication/figures/FIGURE_MANIFEST.csv"
ENVIRONMENT = ROOT / "environment.yml"

PATH24_TREE = "trees/desmognathus_time_tree_path24_v1.nwk"
SOURCE_TREE = "trees/desmo900dated_source46.tre"
FIGURE_1_DATASETS = {"S04", "S12", "S24", "S42"}
REFERENCE_PALETTE = {
    "#0072B2",  # accessible blue
    "#E69F00",  # accessible orange
    "#009E73",  # accessible green
    "#D55E00",  # accessible vermillion
    "#56B4E9",  # accessible sky blue
    "#6A3D9A",  # accessible purple
    "#222222",  # reference ink
}


def _extract_calls(source: str, function_name: str) -> list[str]:
    """Return balanced R calls, ignoring parentheses in strings and comments."""
    calls: list[str] = []
    pattern = re.compile(rf"\b{re.escape(function_name)}\s*\(")
    for match in pattern.finditer(source):
        opening = source.find("(", match.start())
        depth = 0
        quote: str | None = None
        escaped = False
        in_comment = False
        for index in range(opening, len(source)):
            character = source[index]
            if in_comment:
                if character == "\n":
                    in_comment = False
                continue
            if quote is not None:
                if escaped:
                    escaped = False
                elif character == "\\":
                    escaped = True
                elif character == quote:
                    quote = None
                continue
            if character in {'"', "'", "`"}:
                quote = character
            elif character == "#":
                in_comment = True
            elif character == "(":
                depth += 1
            elif character == ")":
                depth -= 1
                if depth == 0:
                    calls.append(source[match.start() : index + 1])
                    break
        else:
            raise AssertionError(f"Unbalanced {function_name}() call in figure source")
    return calls


def _top_level_arguments(call: str) -> list[str]:
    body = call[call.find("(") + 1 : call.rfind(")")]
    arguments: list[str] = []
    start = 0
    stack: list[str] = []
    quote: str | None = None
    escaped = False
    in_comment = False
    matching = {")": "(", "]": "[", "}": "{"}

    for index, character in enumerate(body):
        if in_comment:
            if character == "\n":
                in_comment = False
            continue
        if quote is not None:
            if escaped:
                escaped = False
            elif character == "\\":
                escaped = True
            elif character == quote:
                quote = None
            continue
        if character in {'"', "'", "`"}:
            quote = character
        elif character == "#":
            in_comment = True
        elif character in "([{":
            stack.append(character)
        elif character in ")]}":
            assert stack and stack.pop() == matching[character]
        elif character == "," and not stack:
            arguments.append(body[start:index].strip())
            start = index + 1
    arguments.append(body[start:].strip())
    return [argument for argument in arguments if argument]


def _named_arguments(call: str) -> dict[str, str]:
    named: dict[str, str] = {}
    for argument in _top_level_arguments(call):
        match = re.match(
            r"^([A-Za-z.][A-Za-z0-9._]*)\s*=\s*(.+)$", argument, re.DOTALL
        )
        if match:
            named[match.group(1)] = match.group(2).strip()
    return named


def _figure_1_block(source: str) -> str:
    assert "## ---- figure-1" in source
    assert "## ---- figure-2" in source
    return source.split("## ---- figure-1", maxsplit=1)[1].split(
        "## ---- figure-2", maxsplit=1
    )[0]


def _figure_2_block(source: str) -> str:
    assert "## ---- figure-2" in source
    assert "## ---- figure-3" in source
    return source.split("## ---- figure-2", maxsplit=1)[1].split(
        "## ---- figure-3", maxsplit=1
    )[0]


def _figure_block(source: str, chunk: str, next_chunk: str) -> str:
    start = f"## ---- {chunk}"
    end = f"## ---- {next_chunk}"
    assert start in source
    assert end in source
    return source.split(start, maxsplit=1)[1].split(end, maxsplit=1)[0]


def _references_path(source: str, relative_path: str) -> bool:
    if relative_path in source:
        return True
    directory, filename = relative_path.split("/", maxsplit=1)
    separated_path = re.compile(
        rf"file\.path\(\s*DATA\s*,\s*['\"]{re.escape(directory)}['\"]\s*,"
        rf"\s*['\"]{re.escape(filename)}['\"]\s*\)"
    )
    return separated_path.search(source) is not None


def test_theme_uses_arial_and_minimal_gridless_panels() -> None:
    theme = THEME_SOURCE.read_text()

    assert re.search(r"GBE_FONT\s*<-\s*['\"]Arial['\"]", theme)
    assert "GBE_FALLBACK_FONT" not in theme
    assert re.search(r"GBE_SINGLE_COLUMN_MM\s*<-\s*88", theme)
    assert re.search(r"GBE_DOUBLE_COLUMN_MM\s*<-\s*180", theme)
    assert re.search(r"GBE_LINE_MM\s*<-\s*1\s*/\s*ggplot2::\.pt", theme)
    assert re.search(r"panel\.border\s*=\s*element_blank\(\)", theme)
    assert "theme_gbe_framed_facets" in theme
    assert re.search(r"panel\.border\s*=\s*element_rect\(", theme)
    assert re.search(r"panel\.grid\.major\s*=\s*element_blank\(\)", theme)
    assert re.search(r"panel\.grid\.minor\s*=\s*element_blank\(\)", theme)
    assert not re.search(
        r"panel\.grid(?:\.[A-Za-z]+)*\s*=\s*element_(?:line|rect)\(", theme
    )


def test_theme_contains_the_reference_inspired_palette() -> None:
    theme = THEME_SOURCE.read_text().upper()
    missing = sorted(REFERENCE_PALETTE.difference(set(re.findall(r"#[0-9A-F]{6}", theme))))
    assert not missing, f"Reference palette colors missing from GBE_COLORS: {missing}"


def test_all_patchwork_panel_tags_are_lowercase_and_unparenthesized() -> None:
    source = FIGURE_SOURCE.read_text()
    tagged_calls = [
        call
        for call in _extract_calls(source, "plot_annotation")
        if "tag_levels" in _named_arguments(call)
    ]

    assert tagged_calls, "Expected patchwork plot_annotation() calls with tag_levels"
    for call in tagged_calls:
        arguments = _named_arguments(call)
        assert arguments.get("tag_levels", "").strip("'\"") == "a", call
        assert "tag_prefix" not in arguments, call
        assert "tag_suffix" not in arguments, call


def test_figure1_tree_uses_ape_and_native_ggplot_segments() -> None:
    source = FIGURE_SOURCE.read_text()
    figure_1 = _figure_1_block(source)
    dependency_text = source + "\n" + ENVIRONMENT.read_text()

    assert re.search(r"(?:library\(\s*ape\s*\)|ape::read\.tree\s*\()", source)
    assert "prune_source_tree(overview_traits$species)" in figure_1
    assert _references_path(source, SOURCE_TREE)
    assert re.search(r"\bggplot\s*\(", figure_1)
    assert re.search(r"\bgeom_segment\s*\(", figure_1)
    assert "overview_top_superfamilies" in figure_1
    assert "overview_landscape" in figure_1
    assert "overview_tip_count != 21L" in figure_1
    assert "c(6L, 13L, 2L)" in figure_1

    assert not re.search(
        r"(?:library|require)\(\s*['\"]?ggtree|ggtree::|\bggtree\s*\(",
        dependency_text,
        re.IGNORECASE,
    )


def test_figure1_source_and_manifest_declare_the_integrated_topology_inputs() -> None:
    source = FIGURE_SOURCE.read_text()
    figure_1 = _figure_1_block(source)
    record_calls = [
        call
        for call in _extract_calls(figure_1, "record_figure")
        if re.search(r"['\"]Figure_1['\"]", call)
    ]
    assert len(record_calls) == 1
    source_record = record_calls[0]

    source_ids = set(re.findall(r"\bS\d+\b", source_record))
    assert FIGURE_1_DATASETS <= source_ids
    source_language = source_record.lower()
    assert "topology-only" in source_language
    assert "branch lengths are omitted" in source_language
    assert "dnaPipeTE".lower() in source_language
    assert "repeatmasker" in source_language

    with FIGURE_MANIFEST.open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    figure_rows = [row for row in rows if row["figure_id"] == "Figure_1"]
    assert len(figure_rows) == 1
    row = figure_rows[0]
    manifest_ids = set(row["datasets"].split(";"))
    assert FIGURE_1_DATASETS <= manifest_ids
    manifest_language = f"{row['legend']} {row['alt_text']}".lower()
    assert row["basename"] == "Figure_1_integrated_TE_landscape"
    assert "topology-only" in manifest_language
    assert "conditional genome" in manifest_language
    assert "dnaPipeTE".lower() in manifest_language
    assert "repeatmasker" in manifest_language


def test_figure2_uses_the_validated_bonsai_and_strategy_release() -> None:
    source = FIGURE_SOURCE.read_text()
    figure_2 = _figure_2_block(source)
    diversity_input = figure_2.split("te_annotations <-", maxsplit=1)[0]

    assert "nrow(te_diversity) != 68L" in diversity_input
    assert "n_superfamilies_included != 10L" in diversity_input
    assert "mass_aware_unresolved_bin" not in figure_2
    assert "nrow(diversity_counts) != 2L" in diversity_input
    assert "!all(diversity_counts$n == 34L)" in diversity_input

    record_calls = [
        call
        for call in _extract_calls(figure_2, "record_figure")
        if re.search(r"['\"]Figure_2['\"]", call)
    ]
    assert len(record_calls) == 1
    record_language = record_calls[0].lower()
    assert "24-superfamily clr matrix" in record_language
    assert "display root has no phylogenetic meaning" in record_language
    assert "pca outputs remain supplementary audits" in record_language
    assert "Supplementary_Data_S42_te_bonsai_tip_annotations.csv" in figure_2
    assert "Supplementary_Data_S43_te_bonsai_strategy_validation.csv" in figure_2
    assert "te34_bonsai_composition_tree_with_feature_support.nwk" in figure_2
    assert "S09;S10;S11" not in record_calls[0]
    assert "not a species phylogeny" in record_language
    assert "99,999 permutations" in record_calls[0]
    assert "200 bounded feature-resampling refits" in record_language
    assert "ordinary label permutations are not phylogenetically corrected" in record_language
    assert "not evidence for discrete natural clusters" in record_language
    assert "unassigned" not in figure_2.lower()
    assert "sort(c(7L, 24L, 3L))" in figure_2
    assert "nrow(bonsai_topology_nodes) != 2L" in figure_2
    assert "global_fitch_observed_steps) != 2L" in figure_2
    assert "te_strategy_permanova$variance_explained_r2" in figure_2
    assert "te_strategy_permdisp$p_value" in figure_2
    assert "55/200" in record_calls[0]
    assert "70/200" in record_calls[0]
    assert "2/200" in record_calls[0]


def test_species_axes_share_the_supplied_phylogenetic_order() -> None:
    source = FIGURE_SOURCE.read_text()
    figure_3 = _figure_block(source, "figure-3", "figure-4")
    figure_4 = _figure_block(source, "figure-4", "figure-5")
    figure_5 = _figure_block(source, "figure-5", "figure-6")

    assert _references_path(source, SOURCE_TREE)
    assert "canonical_phylogenetic_order" in source
    assert "phylogenetic_species_order" in source
    assert "phylogenetic_species_levels" in source
    assert "rotateConstr" in source

    assert "phylogenetic_species_levels(species)" in figure_3
    assert "reorder(species_label, young_repeat_percent)" not in figure_3
    assert "shared supplied-tree order" in figure_3
    assert "facet_grid(" in figure_3
    assert "profile_group" in figure_3
    assert "sort(c(3L, 7L, 24L, 34L))" in figure_3
    assert "geom_line" in figure_3
    assert "nrow(young_repeat_matrix) != 34L * 3L" in figure_3

    assert "prune_source_tree(ltr$species)" in figure_4
    assert "ltr_tip_count != 30L" in figure_4
    assert "trees/desmo900dated_source46.tre" in figure_4
    assert "ltr_tip_y_order" in figure_4
    assert "ltr_canonical_order" in figure_4
    assert "not a phylogenetic analysis" in figure_4
    assert "arrange(estimate, species)" not in figure_4
    assert "geom_violin" in figure_4
    assert 'orientation = "y"' in figure_4
    assert "aes(group = species)" in figure_4
    assert "geom_crossbar" not in figure_4
    assert "nrow(ltr_elements) != 380L" in figure_4
    assert "position_jitter" in figure_4
    assert "median_bootstrap_ci_low" in figure_4
    assert "ratio_terminal_internal_all_positions" in figure_4
    assert "equal terminal and internal mean depth" in figure_4
    assert "direction expected under internal-sequence depletion" in figure_4
    assert "Species phylogeny (topology only)" in figure_4

    assert "phylogenetic_species_order(path_tree$tip.label)" in figure_5
    assert "Figure 5 tracks must align to all 24 Path24 tree tips" in figure_5
    assert "path_topology_depth" in figure_5
    assert "trait_forest_panel" in figure_5
    assert "Supplementary_Data_S18_frozen_cell_nucleus_objects.csv" in figure_5
    assert "nrow(morphology_objects) != 1152L" in figure_5
    assert "position_jitter" not in figure_5


def test_wang_style_integrated_candidate_is_data_bound_and_phylogenetically_aligned(
) -> None:
    source = FIGURE_SOURCE.read_text()
    candidate = _figure_block(
        source,
        "wang-style-integrated-te-history-candidate",
        "wang-style-diversity-genome-size-candidate",
    )

    for phrase in [
        "prune_source_tree(integrated_species)",
        "integrated_composition_panel",
        "integrated_fraction_panel",
        'c("LTR", "LINE", "TIR", "DIRS", "SINE", "Maverick", "Other")',
        "aligned RepeatMasker hit bp",
        "divergence bins are not insertion ages",
        "same phylogenetic order",
    ]:
        assert phrase in candidate
    assert "nearest_genome_quantile_species" not in candidate
    assert "geom_arc_bar" not in candidate


def test_wang_style_diversity_candidate_uses_pgls_and_refuses_the_n2_fit() -> None:
    source = FIGURE_SOURCE.read_text()
    candidate = _figure_block(
        source,
        "wang-style-diversity-genome-size-candidate",
        "path24-development-traits-reference",
    )

    for phrase in [
        "all = 21L",
        "aquatic_eggs = 6L",
        "terrestrial_eggs_aquatic_larvae = 13L",
        "direct_development = 2L",
        "if (n_species < 3L)",
        'model = "lambda"',
        "PGLS not fit; n = %d",
        "six PGLS fits and two n=2 no-fit panels",
        "Reproductive-strategy panels are post-hoc descriptive strata",
        "without a regression because n = 2",
    ]:
        assert phrase in candidate


def test_cross_study_group_panels_do_not_single_out_a_published_species() -> None:
    source = FIGURE_SOURCE.read_text()
    cross_study = _figure_block(
        source,
        "wang2025-cross-study-diversity-reference",
        "figure-manifest",
    )

    assert 'aes(label = "D. ochrophaeus")' not in cross_study
    assert "is_published_desmognathus" not in cross_study
    assert "focal_published_desmognathus" not in cross_study


def test_genome_size_is_the_default_trait_in_figure_6() -> None:
    source = FIGURE_SOURCE.read_text()
    figure_6 = _figure_block(source, "figure-6", "figure-7").replace("\\n", " ")

    for phrase in [
        "genome_size_pg_fuscus_anchored",
        "genome_size_pg_ci_low",
        "genome_size_pg_ci_high",
        "genome_size_vs_nucleus_area",
        "genome_size_vs_cell_area",
        "Figure_6_genome_size_pairwise_pgls",
        "Estimated genome size (pg/1C*)",
    ]:
        assert phrase in figure_6
    assert "relative_iod" not in figure_6.lower()
    assert "Figure_S1" not in source
    assert "Supplementary Figure S1" not in source

    for obsolete in [
        "Figure_6_relative_iod_pairwise_pgls",
        "Figure_S1_te_relative_iod_overlap21",
        "relative_iod_vs_nucleus_area",
        "relative_iod_vs_cell_area",
    ]:
        assert obsolete not in source


def test_direct_text_geoms_use_the_seven_point_label_constant() -> None:
    theme = THEME_SOURCE.read_text()
    source = FIGURE_SOURCE.read_text()

    minimum_match = re.search(
        r"GBE_MIN_TEXT_PT\s*<-\s*([0-9]+(?:\.[0-9]+)?)", theme
    )
    assert minimum_match is not None
    assert float(minimum_match.group(1)) >= 7
    assert re.search(
        r"GBE_LABEL_SIZE_MM\s*<-\s*GBE_MIN_TEXT_PT\s*/\s*(?:ggplot2::)?\.pt",
        theme,
    )

    text_calls: list[str] = []
    for function_name in ("geom_text", "geom_text_repel", "annotate"):
        for call in _extract_calls(source, function_name):
            if function_name != "annotate" or "label" in _named_arguments(call):
                text_calls.append(call)

    assert text_calls, "Expected direct text geoms in the publication builder"
    offenders: list[str] = []
    for call in text_calls:
        size = _named_arguments(call).get("size")
        if size != "GBE_LABEL_SIZE_MM":
            offenders.append(" ".join(call.split())[:180])
    assert not offenders, (
        "Every direct geom_text/geom_text_repel/text annotate call must use "
        f"GBE_LABEL_SIZE_MM (>=7 pt); offending calls: {offenders}"
    )
