#!/usr/bin/env Rscript

# Compatibility wrapper for the canonical clade assignment workflow.
script_dir <- tryCatch({
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = FALSE))
  } else {
    "scripts/R/analysis"
  }
}, error = function(...) {
  "scripts/R/analysis"
})

project_root <- normalizePath(file.path(script_dir, "..", "..", ".."), mustWork = FALSE)
source(file.path(project_root, "scripts", "processing", "process_clade_assignments.R"))
