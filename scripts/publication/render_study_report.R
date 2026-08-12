#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(rmarkdown))

root <- normalizePath(getwd(), mustWork = TRUE)
source_path <- file.path(root, "Publication", "Desmognathus_study_data.Rmd")
output_dir <- file.path(root, "Publication")
output_file <- "Desmognathus_study_data.html"

if (!file.exists(source_path)) {
  stop("Missing study-data report source: ", source_path)
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

message("Wrote researcher-facing study report to ", expected)
