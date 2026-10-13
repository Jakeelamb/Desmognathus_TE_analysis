from __future__ import annotations

import csv
import hashlib
import json
import os
import re
import shutil
import subprocess
from html.parser import HTMLParser
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
CONTRACT = ROOT / "scripts/publication/auxiliary_figure_contract.R"
# Exact retained release bytes from baseline 0819d226; independent of fresh producers.
RETAINED_MANIFEST = ROOT / "tests/fixtures/auxiliary_figure_manifest_retained.csv"
RSCRIPT = Path(os.environ["CONDA_PREFIX"]) / "bin/Rscript"


@pytest.fixture(scope="module")
def retained_manifest() -> Path:
    assert hashlib.sha256(RETAINED_MANIFEST.read_bytes()).hexdigest() == (
        "e1fe7a42c29eb0ecd5ccb48ab04b5f8133a0cc0e14f57d93d14684e6e91b4399"
    )
    return RETAINED_MANIFEST


def normalize_report_text(value: str) -> str:
    return " ".join(value.translate(str.maketrans({"‘": "'", "’": "'", "“": '"', "”": '"'})).split())


class RenderedNarrativeText(HTMLParser):
    def __init__(self) -> None:
        super().__init__()
        self.parts: list[str] = []
        self.hidden_element: str | None = None

    def handle_starttag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        if tag in {"script", "style"}:
            self.hidden_element = tag
        elif tag in {"p", "div", "section", "br", "tr", "td", "th", "h1", "h2", "h3"}:
            self.parts.append(" ")

    def handle_endtag(self, tag: str) -> None:
        if tag == self.hidden_element:
            self.hidden_element = None
        elif tag in {"p", "div", "section", "tr", "td", "th", "h1", "h2", "h3"}:
            self.parts.append(" ")

    def handle_data(self, data: str) -> None:
        if self.hidden_element is None:
            self.parts.append(data)


def rendered_report_text(rendered: str) -> str:
    parser = RenderedNarrativeText()
    parser.feed(rendered)
    return normalize_report_text("".join(parser.parts))


def test_report_text_preserves_markup_content_and_literal_inequalities() -> None:
    rendered = (
        '<h3>Cross-study reference (Wang <em>et al.</em> 2025)</h3>'
        '<p>Those orders’ fractions are &lt;1% and &gt;0%; not genome occupancy.</p>'
        '<script>hidden scientific caveat</script><style>hidden title</style>'
    )
    assert rendered_report_text(rendered) == (
        "Cross-study reference (Wang et al. 2025) "
        "Those orders' fractions are <1% and >0%; not genome occupancy."
    )
    assert normalize_report_text("fractions <1%") != normalize_report_text("fractions >1%")


class SourceFingerprintRecords(HTMLParser):
    def __init__(self) -> None:
        super().__init__()
        self.records: list[tuple[str, str | None]] = []

    def handle_starttag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        attributes = dict(attrs)
        if tag == "span" and "data-report-source" in attributes:
            self.records.append((attributes["data-report-source"], attributes.get("data-source-sha256")))


def assert_report_lineage(
    root: Path, source: Path, output: Path, contract: Path,
    historical: dict[str, str], historical_source_sha: str, *, historical_source_record: bool,
) -> None:
    rendered = output.read_text()
    source_artifact = str(source.relative_to(root))
    records = SourceFingerprintRecords()
    records.feed(rendered)
    if records.records:
        for artifact in (source, contract):
            recorded = [digest for name, digest in records.records if name == str(artifact.relative_to(root))]
            assert recorded == [hashlib.sha256(artifact.read_bytes()).hexdigest()]
    else:
        assert historical[source_artifact] == historical_source_sha
        assert hashlib.sha256(output.read_bytes()).hexdigest() == historical[str(output.relative_to(root))]
        if historical_source_record:
            matching = []
            for table in re.findall(
                r'<script\b[^>]*data-pagedtable-source[^>]*>(.*?)</script>', rendered, re.DOTALL,
            ):
                payload = json.loads(table)
                columns = {column["label"][0]: str(column["name"][0]) for column in payload["columns"]}
                if {"artifact", "sha256"} <= columns.keys():
                    matching.extend(
                        row[columns["sha256"]] for row in payload["data"]
                        if row[columns["artifact"]] == source_artifact
                    )
            assert matching == [historical_source_sha]


@pytest.mark.parametrize("historical_record", [True, False])
def test_report_lineage_accepts_only_registered_history_or_current_records(
    tmp_path: Path, historical_record: bool,
) -> None:
    source = tmp_path / "Publication/report.Rmd"
    source.parent.mkdir()
    source.write_text("current template")
    contract = tmp_path / "contract.R"
    contract.write_text("current contract")
    output = source.with_suffix(".html")
    old_digest = hashlib.sha256(b"old template").hexdigest()
    historical_table = '<script data-pagedtable-source type="application/json">' + json.dumps({
        "columns": [{"label": ["artifact"], "name": [1]}, {"label": ["sha256"], "name": [2]}],
        "data": [{"1": "Publication/report.Rmd", "2": old_digest}],
    }) + "</script>"
    output.write_text(historical_table if historical_record else "preserved gallery")
    registry = {"Publication/report.Rmd": old_digest, "Publication/report.html": hashlib.sha256(output.read_bytes()).hexdigest()}

    def check() -> None:
        assert_report_lineage(tmp_path, source, output, contract, registry, old_digest,
                              historical_source_record=historical_record)

    check()
    current = "\n".join(
        f'<span data-report-source="{path.relative_to(tmp_path)}" '
        f'data-source-sha256="{hashlib.sha256(path.read_bytes()).hexdigest()}"></span>'
        for path in (source, contract)
    )
    output.write_text(current)
    check()
    bad_lineages = [
        current.replace(hashlib.sha256(source.read_bytes()).hexdigest(), old_digest),
        current.replace("Publication/report.Rmd", "Publication/another.Rmd"),
        current + current,
        hashlib.sha256(source.read_bytes()).hexdigest(),
        historical_table.replace(old_digest, hashlib.sha256(source.read_bytes()).hexdigest()),
    ]
    for rendered in bad_lineages:
        output.write_text(rendered)
        with pytest.raises(AssertionError):
            check()


LOAD_VARIANT = r'''
args <- commandArgs(trailingOnly = TRUE)
source(args[[1]])
manifest <- read.csv(args[[2]], stringsAsFactors = FALSE, check.names = FALSE)
if (args[[3]] == "current") {
  literal <- function(expression) {
    if (is.character(expression)) return(expression)
    if (is.call(expression) && as.character(expression[[1]]) %in% c("paste", "paste0") &&
        all(vapply(as.list(expression)[-1], is.character, logical(1)))) {
      return(eval(expression, baseenv()))
    }
    stop("Producer metadata must be a string literal or paste of string literals")
  }
  calls <- Filter(function(expression) is.call(expression) &&
    identical(expression[[1]], as.name("record_auxiliary_figure")), parse(file = args[[4]]))
  stopifnot(length(calls) == 5L)
  for (expression in calls) {
    index <- match(literal(expression[[2]]), manifest$auxiliary_id)
    fields <- c(title = 3L, boundary = 9L, legend = 10L, alt_text = 11L)
    for (field in names(fields)) manifest[[field]][index] <- literal(expression[[fields[[field]]]])
  }
}
contract <- validate_auxiliary_figure_manifest(manifest)
stopifnot(identical(contract$variant, args[[3]]))
'''


def run_r(code: str, *arguments: str, cwd: Path = ROOT) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [str(RSCRIPT), "--vanilla", "-", *arguments],
        input=code,
        cwd=cwd,
        text=True,
        capture_output=True,
        check=False,
        timeout=30,
    )


@pytest.mark.parametrize("variant", ["retained", "current"])
def test_manifest_variants_drive_the_actual_integrated_narrative(
    tmp_path: Path, variant: str, retained_manifest: Path,
):
    arguments = [
        str(CONTRACT), str(retained_manifest), variant,
        str(ROOT / "scripts/publication/build_gbe_figures.R"),
    ]
    functions = []
    fingerprint_snippets = []
    for name in ("study_data", "figure_review"):
        source = (ROOT / f"Publication/Desmognathus_{name}.Rmd").read_text()
        for function in ("escape_html", "format_inline"):
            functions.append(re.search(rf"{function} <- function.*?\n\}}", source, re.DOTALL)[0])
        if name == "study_data":
            start = source.index("**Wang-style integrated TE-history candidate boundary.**")
            end = source.index("**Wang-style TE-diversity candidate boundary.**", start)
        else:
            start = source.index("```{r integrated-history-boundary")
            end = source.index("```", start + 3) + 3
        snippet = tmp_path / f"{name}.Rmd"
        snippet.write_text(source[start:end])
        arguments.append(str(snippet))
        start = source.index("```{r source-fingerprint-records")
        end = source.index("```", start + 3) + 3
        fingerprint_snippet = tmp_path / f"{name}-fingerprint.Rmd"
        fingerprint_snippet.write_text(source[start:end])
        fingerprint_snippets.append(str(fingerprint_snippet))
    arguments.extend(fingerprint_snippets)
    code = LOAD_VARIANT + "\n".join(functions) + r'''
auxiliary_manifest <- manifest
ROOT <- normalizePath(file.path(dirname(args[[4]]), "..", ".."))
for (path in args[5:6]) {
  knitted <- knitr::knit(path, output = paste0(path, ".md"), quiet = TRUE)
  text <- paste(readLines(knitted, warn = FALSE), collapse = " ")
  plain <- gsub("[[:space:]]+", " ", gsub("<[^>]+>", "", text))
  stopifnot(grepl(manifest$boundary[[1]], plain, fixed = TRUE))
  stopifnot(grepl(manifest$legend[[1]], plain, fixed = TRUE))
  if (contract$variant == "current") {
    stopifnot(!grepl("Six representatives are selected", text, fixed = TRUE))
  }
}
for (index in 1:2) {
  path <- args[[6 + index]]
  knitted <- knitr::knit(path, output = paste0(path, ".md"), quiet = TRUE)
  html <- paste(readLines(knitted, warn = FALSE), collapse = "\n")
  source_path <- file.path(ROOT, "Publication", paste0("Desmognathus_",
                           c("study_data", "figure_review")[[index]], ".Rmd"))
  validate_report_source_fingerprint(source_path, html)
}
html <- paste(paste0("<section id='", manifest$auxiliary_id, "'><h2>", manifest$title,
                     "</h2>", manifest$boundary, " ", manifest$legend, " ",
                     manifest$alt_text, "</section>"), collapse = "\n")
validate_auxiliary_rendered_html(contract, html)
broken <- gsub("divergence bins are not insertion ages", "", html, fixed = TRUE)
stopifnot(inherits(tryCatch(validate_auxiliary_rendered_html(contract, broken),
                           error = identity), "error"))
path <- args[[1]]
artifact <- paste0("Publication/", basename(path))
record <- report_source_fingerprint_html(path, artifact)
validate_report_source_fingerprint(path, record)
for (invalid in c("stale", figure_source_sha256(path), paste0(record, record),
                  sub(figure_source_sha256(path), paste(rep("0", 64), collapse = ""), record))) {
  stopifnot(inherits(tryCatch(validate_report_source_fingerprint(path, invalid),
                             error = identity), "error"))
}
cat("PASS: ", contract$variant, " metadata, both actual narrative snippets and rejection controls\n")
'''
    result = run_r(code, *arguments, cwd=tmp_path)
    assert result.returncode == 0, result.stdout + result.stderr
    assert f"PASS:  {variant}" in result.stdout


@pytest.mark.parametrize("renderer", ["render_study_report", "render_figure_review"])
@pytest.mark.parametrize(
    "problem,diagnostic",
    [
        ("id", "five ordered IDs"),
        ("order", "five ordered IDs"),
        ("title", "known ordered variant"),
        ("mixed-titles", "known ordered variant"),
        ("schema", "required columns"),
        ("dimensions", "positive dimensions"),
        ("boundary", "missing required review-boundary text"),
    ],
)
def test_renderers_reject_bad_manifests_before_writing(
    tmp_path: Path, renderer: str, problem: str, diagnostic: str, retained_manifest: Path,
):
    figures = tmp_path / "Publication/figures"
    figures.mkdir(parents=True)
    helper = tmp_path / "scripts/publication/auxiliary_figure_contract.R"
    helper.parent.mkdir(parents=True)
    shutil.copy2(CONTRACT, helper)
    shutil.copy2(ROOT / "Publication/figures/FIGURE_MANIFEST.csv", figures)
    for name in ("study_data", "figure_review"):
        (tmp_path / f"Publication/Desmognathus_{name}.Rmd").write_text("fixture")
    with retained_manifest.open(newline="") as handle:
        reader = csv.DictReader(handle)
        fields, rows = reader.fieldnames, list(reader)
    if problem == "id":
        rows[0]["auxiliary_id"] = "unknown"
    elif problem == "order":
        rows[0], rows[1] = rows[1], rows[0]
    elif problem == "title":
        rows[0]["title"] = "unknown"
    elif problem == "mixed-titles":
        rows[0]["title"] = "Wang-style integrated TE history across all 34 species"
    elif problem == "schema":
        fields.remove("boundary")
        for row in rows:
            del row["boundary"]
    elif problem == "dimensions":
        rows[0]["width_mm"] = "0"
    else:
        rows[0]["boundary"] = rows[0]["boundary"].replace(
            "divergence bins are not insertion ages", ""
        )
    with (figures / "AUXILIARY_FIGURE_MANIFEST.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    result = subprocess.run(
        [str(RSCRIPT), "--vanilla", str(ROOT / f"scripts/publication/{renderer}.R")],
        cwd=tmp_path,
        text=True,
        capture_output=True,
        check=False,
        timeout=30,
    )
    assert result.returncode == 1, result.stdout + result.stderr
    assert diagnostic in result.stderr
    assert not list((tmp_path / "Publication").glob("*.html"))
