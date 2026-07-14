#!/usr/bin/env Rscript

# Execute the original Research Talk IV PCA/clustering plotting program from
# the repository history. Only package bootstrap, command-line parsing, and
# current input/output paths are adapted; the data preparation, PCA, elbow,
# silhouette, cluster, and ggplot code are the historical source verbatim.

args <- commandArgs(trailingOnly = FALSE)
flag <- grep("^--file=", args, value = TRUE)
script_dir <- if (length(flag)) dirname(normalizePath(sub("^--file=", "", flag[[1]]), mustWork = FALSE)) else "scripts/processing"
root <- normalizePath(file.path(script_dir, "..", ".."), mustWork = TRUE)
historical_ref <- "1ad5e8c:te_pca_analysis.R"
output_dir <- file.path(root, "results", "figures", "talk_iv_restored_pca")
input_source <- file.path(root, "results", "data", "legacy_descriptive", "te_superfamily_composition_te34_legacy_descriptive_v1.csv")
staged_input <- file.path(root, "results", "data", "legacy_descriptive", "talk_iv_superfamily_pca_input.csv")

if (!file.exists(input_source)) stop("Run build_legacy_descriptive_te_figures.R first")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
input <- read.csv(input_source, check.names = FALSE)
names(input)[names(input) == "species"] <- "Species"
input$Species <- paste0("D.", input$Species)
write.csv(input, staged_input, row.names = FALSE)

code <- system2("git", c("show", historical_ref), stdout = TRUE, stderr = TRUE)
if (!length(code) || any(grepl("fatal:", code, fixed = TRUE))) stop("Could not recover historical PCA source: ", historical_ref)

# The original imports include optional packages that its plotting body never
# calls. The active R environment has the packages actually used below.
library_start <- grep("^suppressPackageStartupMessages\\(\\{", code)
library_end <- which(seq_along(code) > library_start & code == "})")[1]
code <- c(
  code[seq_len(library_start - 1L)],
  "suppressPackageStartupMessages({ library(ggrepel); library(ggplot2); library(cluster); library(viridis); library(dplyr) })",
  code[(library_end + 1L):length(code)]
)

# Replace only the optparse setup with fixed current paths. The original
# plotting/analysis code after this block is retained unchanged.
parse_start <- grep("^# Parse command line arguments", code)
create_dir <- grep("^# Create output directory", code)
code <- c(
  code[seq_len(parse_start - 1L)],
  sprintf("opt <- list(input = %s, output_dir = %s, filter_species = '', optimal_k = 5L)", deparse(staged_input), deparse(output_dir)),
  code[create_dir:length(code)]
)

# The current input is named Species already, but leave this guard directly
# after the historical reader so the recovered program remains robust.
reader <- grep('df <- read.csv\\(opt\\$input\\)', code)
code <- append(code, 'if (!("Species" %in% names(df))) names(df)[1] <- "Species"', after = reader)

eval(parse(text = paste(code, collapse = "\n")), envir = new.env(parent = globalenv()))
cat("Historical Talk IV PCA figures restored in ", output_dir, "\n", sep = "")
