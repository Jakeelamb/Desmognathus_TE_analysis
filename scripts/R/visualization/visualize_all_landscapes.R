#!/usr/bin/env Rscript
# Compatibility wrapper for the canonical multi-sample landscape plotter.

script_args <- commandArgs(trailingOnly = FALSE)
script_path_arg <- grep("^--file=", script_args, value = TRUE)
if (length(script_path_arg) > 0) {
  script_path <- normalizePath(sub("^--file=", "", script_path_arg[1]), mustWork = FALSE)
  script_dir <- dirname(script_path)
} else {
  script_dir <- "scripts/R/visualization"
}

project_root <- normalizePath(file.path(script_dir, "..", "..", ".."), mustWork = FALSE)
canonical_script <- file.path(project_root, "scripts", "visualization", "plot_all_te_landscapes.R")

if (!file.exists(canonical_script)) {
  stop("Canonical landscape plotting script not found: ", canonical_script, call. = FALSE)
}

message("Delegating to canonical workflow: ", canonical_script)
status <- system2("Rscript", canonical_script)
if (!identical(status, 0L)) {
  stop("Canonical landscape plotting workflow failed with exit status ", status, call. = FALSE)
}
