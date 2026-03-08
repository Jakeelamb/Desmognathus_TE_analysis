#!/usr/bin/env Rscript
# Compatibility wrapper for the canonical PERMANOVA workflow.

script_path <- normalizePath(commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))], mustWork = FALSE)
if (length(script_path) == 0) {
  script_dir <- normalizePath("scripts/processing", mustWork = FALSE)
} else {
  script_dir <- dirname(sub("^--file=", "", script_path[[1]]))
}

canonical_script <- file.path(script_dir, "permanova_analysis.R")
if (!file.exists(canonical_script)) {
  stop("Canonical PERMANOVA script not found: ", canonical_script, call. = FALSE)
}

message("Delegating to canonical workflow: ", canonical_script)
status <- system2("Rscript", canonical_script)
if (!identical(status, 0L)) {
  stop("Canonical PERMANOVA workflow failed with exit status ", status, call. = FALSE)
}
