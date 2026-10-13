#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(rmarkdown))

root <- normalizePath(getwd(), mustWork = TRUE)
source_path <- file.path(root, "Publication", "Desmognathus_figure_review.Rmd")
output_dir <- file.path(root, "Publication")
output_file <- "Desmognathus_figure_review.html"
manifest_path <- file.path(root, "Publication", "figures", "FIGURE_MANIFEST.csv")
auxiliary_manifest_path <- file.path(
  root, "Publication", "figures", "AUXILIARY_FIGURE_MANIFEST.csv"
)

if (!file.exists(source_path)) {
  stop("Missing figure-review source: ", source_path)
}
if (!file.exists(manifest_path)) {
  stop("Missing figure manifest: ", manifest_path)
}
if (!file.exists(auxiliary_manifest_path)) {
  stop("Missing auxiliary figure manifest: ", auxiliary_manifest_path)
}

manifest <- utils::read.csv(
  manifest_path,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
auxiliary_manifest <- utils::read.csv(
  auxiliary_manifest_path,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
expected_auxiliary_ids <- c(
  "Auxiliary_Wang_style_integrated_te_history_candidate",
  "Auxiliary_Wang_style_diversity_genome_size_candidate",
  "Auxiliary_Path24_development_traits_reference",
  "Auxiliary_Path_DAG_reference",
  "Auxiliary_Wang2025_cross_study_te_diversity_reference"
)
expected_auxiliary_titles <- c(
  "Wang-style integrated TE content and proliferation candidate",
  "Wang-style TE diversity and genome-size candidate",
  "Path24 development and phenotype reference",
  "Path DAG equivalence-class reference",
  "Cross-study TE-diversity reference (Wang et al. 2025)"
)
if (nrow(manifest) != 7L) {
  stop("Expected exactly 7 numbered figures")
}
if (nrow(auxiliary_manifest) != 5L) {
  stop("Expected exactly 5 auxiliary review figures")
}
if (
  !identical(auxiliary_manifest$auxiliary_id, expected_auxiliary_ids) ||
    !identical(auxiliary_manifest$title, expected_auxiliary_titles)
) {
  stop(
    "Auxiliary figure manifest must contain exactly these ordered IDs and titles: ",
    paste(expected_auxiliary_ids, collapse = ", ")
  )
}

expected_genome_size_figures <- c(
  Figure_6 = "Figure_6_genome_size_pairwise_pgls"
)
for (figure_id in names(expected_genome_size_figures)) {
  row <- manifest[manifest$figure_id == figure_id, , drop = FALSE]
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

figure_two <- manifest[manifest$figure_id == "Figure_2", , drop = FALSE]
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
auxiliary_pngs <- file.path(
  root,
  "Publication",
  "figures",
  auxiliary_manifest$png
)
if (any(!file.exists(auxiliary_pngs))) {
  stop(
    "Missing auxiliary review PNG files: ",
    paste(basename(auxiliary_pngs[!file.exists(auxiliary_pngs)]), collapse = ", ")
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

html <- paste(readLines(expected, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
plain_html <- gsub("<[^>]+>", "", html)
plain_html <- gsub("[[:space:]]+", " ", plain_html)
plain_html <- chartr("\u00A0", " ", plain_html)

count_fixed <- function(pattern, text) {
  matches <- gregexpr(pattern, text, fixed = TRUE)[[1]]
  if (identical(matches[[1]], -1L)) 0L else length(matches)
}

embedded_png_count <- count_fixed("data:image/png;base64", html)
expected_embedded_png_count <- nrow(manifest) + nrow(auxiliary_manifest)
if (expected_embedded_png_count != 12L) {
  stop("Expected exactly 12 review PNGs (7 numbered plus 5 auxiliary)")
}
if (embedded_png_count != expected_embedded_png_count) {
  stop(
    "Expected ", expected_embedded_png_count, " embedded PNG figures, found ",
    embedded_png_count
  )
}
if (!all(vapply(manifest$figure_id, grepl, logical(1), x = html, fixed = TRUE))) {
  stop("Rendered review is missing one or more manifest figure identifiers")
}
if (!grepl("Unresolved tree provenance", html, fixed = TRUE)) {
  stop("Rendered review is missing the tree-provenance warning")
}
for (phrase in c(
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
  "conditional D. fuscus assembly-anchored pg/1C estimates",
  "not direct independently measured species C-values",
  "Path DAG equivalence-class reference",
  "Not a numbered manuscript figure.",
  "10 testable equivalence classes",
  "saturated 6-member class is not shown",
  "remains unscored because it has no d-separation claim",
  "Genome size",
  "conditional D. fuscus assembly-anchored pg/1C analysis",
  "not a direct species C-value measurement",
  "same released genome-size scale is used consistently",
  "do not orient a unique causal arrow",
  "Cross-study TE-diversity reference (Wang et al. 2025)",
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
)) {
  if (!grepl(phrase, plain_html, fixed = TRUE)) {
    stop("Rendered review is missing an auxiliary review boundary phrase: ", phrase)
  }
}

figure_one <- manifest[manifest$figure_id == "Figure_1", , drop = FALSE]
figure_one_evidence <- if (nrow(figure_one) == 1L) {
  paste(
    figure_one$datasets,
    figure_one$legend,
    figure_one$alt_text,
    figure_one$basename,
    collapse = " "
  )
} else {
  ""
}
figure_1_has_phylogeny <- nrow(figure_one) == 1L && grepl(
  "S38|S39|phylogen|dated[ -]tree|time-calibrated",
  figure_one_evidence,
  ignore.case = TRUE,
  perl = TRUE
)
if (
  figure_1_has_phylogeny &&
    !grepl("Figure 1 is the integrated 21-species overview", html, fixed = TRUE)
) {
  stop("Rendered review does not acknowledge the integrated phylogeny in Figure 1")
}

message(
  "Wrote self-contained figure review with ", nrow(manifest),
  " manifest figures plus ", nrow(auxiliary_manifest),
  " auxiliary review figures to ", expected
)
