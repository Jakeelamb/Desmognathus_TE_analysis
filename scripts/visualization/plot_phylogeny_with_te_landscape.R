#!/usr/bin/env Rscript

# Compatibility wrapper for the canonical TE-on-phylogeny visualization workflow.
script_dir <- tryCatch({
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = FALSE))
  } else {
    "scripts/visualization"
  }
}, error = function(...) {
  "scripts/visualization"
})

project_root <- normalizePath(file.path(script_dir, "..", ".."), mustWork = FALSE)
source(file.path(project_root, "scripts", "R", "visualization", "te_phylo_landscape.R"))
