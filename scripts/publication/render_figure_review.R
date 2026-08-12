#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(rmarkdown))

root <- normalizePath(getwd(), mustWork = TRUE)
source_path <- file.path(root, "Publication", "Desmognathus_figure_review.Rmd")
output_dir <- file.path(root, "Publication")
output_file <- "Desmognathus_figure_review.html"
manifest_path <- file.path(root, "Publication", "figures", "FIGURE_MANIFEST.csv")
auxiliary_png <- file.path(
  root, "Publication", "figures", "Path_DAG_equivalence_class_reference.png"
)

if (!file.exists(source_path)) {
  stop("Missing figure-review source: ", source_path)
}
if (!file.exists(manifest_path)) {
  stop("Missing figure manifest: ", manifest_path)
}
if (!file.exists(auxiliary_png)) {
  stop("Missing auxiliary Path24 DAG reference: ", auxiliary_png)
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

manifest <- utils::read.csv(
  manifest_path,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
html <- paste(readLines(expected, warn = FALSE, encoding = "UTF-8"), collapse = "\n")

count_fixed <- function(pattern, text) {
  matches <- gregexpr(pattern, text, fixed = TRUE)[[1]]
  if (identical(matches[[1]], -1L)) 0L else length(matches)
}

embedded_png_count <- count_fixed("data:image/png;base64", html)
expected_embedded_png_count <- nrow(manifest) + 1L
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
  "Auxiliary Path DAG reference",
  "Not a numbered manuscript figure.",
  "10 testable equivalence classes",
  "saturated 6-member class is not shown",
  "remains unscored because it has no d-separation claim",
  "current Path24 relative IOD-only analysis"
)) {
  if (!grepl(phrase, html, fixed = TRUE)) {
    stop("Rendered review is missing the auxiliary DAG boundary phrase: ", phrase)
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
    !grepl("Figure 1 includes the Path24 dated phylogeny", html, fixed = TRUE)
) {
  stop("Rendered review does not acknowledge the dated phylogeny in Figure 1")
}

message(
  "Wrote self-contained figure review with ", nrow(manifest),
  " manifest figures plus one auxiliary DAG review section to ", expected
)
