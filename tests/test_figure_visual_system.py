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
FIGURE_1_DATASETS = {"S04", "S38", "S39"}
REFERENCE_PALETTE = {
    "#00BFC4",  # cyan paired-group accent
    "#F8766D",  # coral paired-group accent
    "#440154",  # viridis purple
    "#31688E",  # viridis blue
    "#35B779",  # viridis green
    "#FDE725",  # viridis yellow
    "#1A1A1A",  # reference ink
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


def _references_path(source: str, relative_path: str) -> bool:
    if relative_path in source:
        return True
    directory, filename = relative_path.split("/", maxsplit=1)
    separated_path = re.compile(
        rf"file\.path\(\s*DATA\s*,\s*['\"]{re.escape(directory)}['\"]\s*,"
        rf"\s*['\"]{re.escape(filename)}['\"]\s*\)"
    )
    return separated_path.search(source) is not None


def test_theme_keeps_nimbus_sans_with_boxed_gridless_panels() -> None:
    theme = THEME_SOURCE.read_text()

    assert re.search(
        r"GBE_FALLBACK_FONT\s*<-\s*['\"]Nimbus Sans['\"]", theme
    )
    assert '"Liberation Sans"' in theme
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


def test_all_patchwork_panel_tags_are_parenthesized() -> None:
    source = FIGURE_SOURCE.read_text()
    tagged_calls = [
        call
        for call in _extract_calls(source, "plot_annotation")
        if "tag_levels" in _named_arguments(call)
    ]

    assert tagged_calls, "Expected patchwork plot_annotation() calls with tag_levels"
    for call in tagged_calls:
        arguments = _named_arguments(call)
        assert arguments.get("tag_prefix", "").strip("'\"") == "(", call
        assert arguments.get("tag_suffix", "").strip("'\"") == ")", call


def test_figure1_tree_uses_ape_and_native_ggplot_segments() -> None:
    source = FIGURE_SOURCE.read_text()
    figure_1 = _figure_1_block(source)
    dependency_text = source + "\n" + ENVIRONMENT.read_text()

    assert re.search(r"(?:library\(\s*ape\s*\)|ape::read\.tree\s*\()", source)
    assert re.search(r"(?:ape::)?read\.tree\s*\(", figure_1)
    assert _references_path(figure_1, PATH24_TREE)
    assert re.search(r"\bggplot\s*\(", figure_1)
    assert re.search(r"\bgeom_segment\s*\(", figure_1)

    assert not re.search(
        r"(?:library|require)\(\s*['\"]?ggtree|ggtree::|\bggtree\s*\(",
        dependency_text,
        re.IGNORECASE,
    )


def test_figure1_source_and_manifest_declare_the_time_tree_inputs() -> None:
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
    assert "time-calibrated" in source_language
    assert "phylogen" in source_language or "dated tree" in source_language

    with FIGURE_MANIFEST.open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    figure_rows = [row for row in rows if row["figure_id"] == "Figure_1"]
    assert len(figure_rows) == 1
    row = figure_rows[0]
    manifest_ids = set(row["datasets"].split(";"))
    assert FIGURE_1_DATASETS <= manifest_ids
    manifest_language = f"{row['legend']} {row['alt_text']}".lower()
    assert "time-calibrated" in manifest_language
    assert "phylogen" in manifest_language or "dated tree" in manifest_language


def test_figure2_uses_one_classified_only_diversity_view() -> None:
    source = FIGURE_SOURCE.read_text()
    figure_2 = _figure_2_block(source)
    diversity_input = figure_2.split("te_scores <-", maxsplit=1)[0]
    diversity_plot = figure_2.split("te_a <-", maxsplit=1)[1].split(
        "te_b <-", maxsplit=1
    )[0]

    assert 'filter(te_level == "order")' in diversity_input
    assert "composition_mode" not in diversity_input
    assert "mass_aware_unresolved_bin" not in figure_2
    assert "Order-composition denominator" not in figure_2
    assert re.search(
        r"facet_wrap\(vars\(metric\),\s*ncol\s*=\s*1,\s*scales\s*=\s*['\"]free_y['\"]\)",
        diversity_plot,
    )
    assert 'fill = GBE_COLORS[["blue"]]' in diversity_plot
    assert "nrow(diversity_counts) != 2L" in diversity_input
    assert "!all(diversity_counts$n == 34L)" in diversity_input

    record_calls = [
        call
        for call in _extract_calls(figure_2, "record_figure")
        if re.search(r"['\"]Figure_2['\"]", call)
    ]
    assert len(record_calls) == 1
    record_language = record_calls[0].lower()
    assert "classified-only" in record_language
    assert "s07" in record_language
    assert "not a diversity category" in record_language
    assert "not included in diversity" in record_language


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
