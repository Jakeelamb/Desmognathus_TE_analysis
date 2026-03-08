#!/usr/bin/env Rscript

args_all <- commandArgs(trailingOnly = FALSE)
script_flag <- grep("^--file=", args_all, value = TRUE)
script_path <- if (length(script_flag) > 0) sub("^--file=", "", script_flag[[1]]) else "scripts/visualization/ectopic_recomb_plots.R"
script_dir <- dirname(normalizePath(script_path, mustWork = FALSE))
canonical_script <- file.path(script_dir, "plot_ectopic_recombination.R")

args <- commandArgs(trailingOnly = TRUE)
status <- system2("Rscript", c(canonical_script, args))
quit(save = "no", status = status)
