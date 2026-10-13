#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(rmarkdown))

root <- normalizePath(getwd(), mustWork = TRUE)
source(file.path(root, "scripts/publication/auxiliary_figure_contract.R"))
source_path <- file.path(root, "Publication", "Desmognathus_study_data.Rmd")
output_dir <- file.path(root, "Publication")
output_file <- "Desmognathus_study_data.html"
numbered_manifest_path <- file.path(
  root, "Publication", "figures", "FIGURE_MANIFEST.csv"
)
auxiliary_manifest_path <- file.path(
  root, "Publication", "figures", "AUXILIARY_FIGURE_MANIFEST.csv"
)
figure_dir <- file.path(root, "Publication", "figures")

if (!file.exists(source_path)) {
  stop("Missing study-data report source: ", source_path)
}
if (!file.exists(numbered_manifest_path)) {
  stop("Missing numbered figure manifest: ", numbered_manifest_path)
}
if (!file.exists(auxiliary_manifest_path)) {
  stop("Missing auxiliary figure manifest: ", auxiliary_manifest_path)
}

numbered_manifest <- utils::read.csv(
  numbered_manifest_path,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
auxiliary_manifest <- utils::read.csv(
  auxiliary_manifest_path,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

if (nrow(numbered_manifest) != 7L) {
  stop("Expected exactly 7 numbered figures")
}
if (nrow(auxiliary_manifest) != 5L) {
  stop("Expected exactly 5 auxiliary review figures")
}
auxiliary_contract <- validate_auxiliary_figure_manifest(auxiliary_manifest)

expected_genome_size_figures <- c(
  Figure_6 = "Figure_6_genome_size_pairwise_pgls"
)
for (figure_id in names(expected_genome_size_figures)) {
  row <- numbered_manifest[numbered_manifest$figure_id == figure_id, , drop = FALSE]
  if (
    nrow(row) != 1L ||
      !identical(row$basename[[1]], expected_genome_size_figures[[figure_id]])
  ) {
    stop(figure_id, " does not use the canonical genome-size figure basename")
  }
  public_copy <- paste(row$legend, row$alt_text, row$release_status, collapse = " ")
  if (!grepl("genome size", public_copy, ignore.case = TRUE)) {
    stop(figure_id, " does not expose genome size as its paper-facing trait")
  }
  if (grepl("relative[- ]?(nuclear[- ]?)?IOD", public_copy, ignore.case = TRUE)) {
    stop(figure_id, " still exposes relative IOD as its paper-facing trait")
  }
}

figure_two <- numbered_manifest[
  numbered_manifest$figure_id == "Figure_2",
  ,
  drop = FALSE
]
if (
  nrow(figure_two) != 1L ||
    !identical(figure_two$basename[[1]], "Figure_2_te_diversity_and_bonsai")
) {
  stop("Figure_2 does not use the canonical Bonsai figure basename")
}
figure_two_copy <- paste(figure_two$legend, figure_two$alt_text, collapse = " ")
for (phrase in c("Bonsai", "not a species phylogeny", "feature resampling")) {
  if (!grepl(tolower(phrase), tolower(figure_two_copy), fixed = TRUE)) {
    stop("Figure_2 is missing its Bonsai boundary phrase: ", phrase)
  }
}

numbered_png_paths <- file.path(figure_dir, numbered_manifest$png)
auxiliary_png_paths <- file.path(figure_dir, auxiliary_manifest$png)
missing_png_paths <- c(
  numbered_png_paths[!file.exists(numbered_png_paths)],
  auxiliary_png_paths[!file.exists(auxiliary_png_paths)]
)
if (length(missing_png_paths) > 0L) {
  stop(
    "Missing report PNG figures: ",
    paste(basename(missing_png_paths), collapse = ", ")
  )
}

rendered <- render(
  input = source_path,
  output_file = output_file,
  output_dir = output_dir,
  knit_root_dir = root,
  clean = TRUE,
  envir = new.env(parent = globalenv()),
  quiet = FALSE
)

expected <- file.path(output_dir, output_file)
if (!identical(normalizePath(rendered, mustWork = TRUE), normalizePath(expected, mustWork = TRUE))) {
  stop("R Markdown rendered to an unexpected path: ", rendered)
}

# Pandoc preserves trailing spaces in generated syntax-highlighting CSS and
# sessionInfo output. Normalize them so the tracked self-contained report
# remains diff-clean without changing its rendered content.
report_lines <- readLines(expected, warn = FALSE, encoding = "UTF-8")
report_lines <- sub("[[:blank:]]+$", "", report_lines)
writeLines(report_lines, expected, useBytes = TRUE)

html <- paste(report_lines, collapse = "\n")
normalized_html <- gsub("<[^>]+>", "", html)
normalized_html <- gsub("[[:space:]]+", " ", normalized_html)

count_fixed <- function(pattern, text) {
  matches <- gregexpr(pattern, text, fixed = TRUE)[[1]]
  if (identical(matches[[1]], -1L)) 0L else length(matches)
}

expected_embedded_png_count <- 12L
embedded_png_count <- count_fixed("data:image/png;base64", html)
if (embedded_png_count != expected_embedded_png_count) {
  stop(
    "Expected ", expected_embedded_png_count,
    " embedded PNG figures (7 numbered plus 5 auxiliary), found ",
    embedded_png_count
  )
}
if (grepl('src="figures/', html, fixed = TRUE)) {
  stop("Study report contains non-embedded figure references")
}

validate_auxiliary_rendered_html(auxiliary_contract, html)
validate_report_source_fingerprint(source_path, html)

required_boundary_phrases <- c(
  "Review candidate, not a numbered manuscript figure",
  "fractions of aligned RepeatMasker hit bp",
  "divergence bins are not insertion ages",
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
  "10 testable Markov-equivalence classes",
  "saturated 6-member class is omitted",
  "Genome size node",
  "conditional D. fuscus assembly-anchored pg/1C estimates",
  "same genome-size scale is used across all scored classes",
  "not direct species C-value measurements",
  "representative only and do not orient unique causal arrows",
  "Descriptive cross-study context only; no pooled regression or hypothesis test",
  "10.5061/dryad.zpc866tkv",
  "10.1038/s42003-025-08127-3",
  "CC0-1.0",
  "published top-10 TE-superfamily indices",
  "RepeatModeler/RepeatMasker",
  "reclosed calculations from dnaPipeTE compositions",
  "pipelines are different and not method-harmonized",
  "conditional D. fuscus assembly-anchored estimates",
  "not direct species C-values"
)
for (phrase in required_boundary_phrases) {
  if (!grepl(phrase, normalized_html, fixed = TRUE)) {
    stop("Rendered report is missing an auxiliary boundary phrase: ", phrase)
  }
}

message(
  "Wrote researcher-facing study report with 7 numbered and 5 ordered ",
  "auxiliary figures to ", expected
)
