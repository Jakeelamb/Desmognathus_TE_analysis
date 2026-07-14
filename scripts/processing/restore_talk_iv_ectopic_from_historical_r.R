#!/usr/bin/env Rscript

# Run the original Talk IV LTR violin/ANOVA source from Git history against a
# staged current TE34 element table. Only input and output paths are adapted.

args <- commandArgs(trailingOnly = FALSE)
flag <- grep("^--file=", args, value = TRUE)
script_dir <- if (length(flag)) dirname(normalizePath(sub("^--file=", "", flag[[1]]), mustWork = FALSE)) else "scripts/processing"
root <- normalizePath(file.path(script_dir, "..", ".."), mustWork = TRUE)
data_dir <- file.path(root, "results", "data", "legacy_descriptive")
figure_dir <- file.path(root, "results", "figures", "talk_iv_restored_ectopic")
input_source <- file.path(root, "results", "data", "ectopic_recombination_filtered_3000bp_5+domains_no_unknown_species.csv")
staged_input <- file.path(data_dir, "talk_iv_ectopic_five_domains.csv")
stats_dir <- file.path(data_dir, "talk_iv_ectopic_stats")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(stats_dir, recursive = TRUE, showWarnings = FALSE)

lookup <- read.delim(file.path(root, "input_data", "lookup_table.txt"))
active <- tolower(sub("^D\\.", "", lookup$Species))
source <- read.delim(input_source, check.names = FALSE)
species <- tolower(sub("^D\\.\\s*", "", source$species))
staged <- data.frame(Species = species, Ratio.of.LTR.to.Internal = as.numeric(source$ratio_terminal_internal))
staged <- staged[staged$Species %in% active & is.finite(staged$Ratio.of.LTR.to.Internal) & staged$Ratio.of.LTR.to.Internal > 0, ]
write.csv(staged, staged_input, row.names = FALSE)

code <- system2("git", c("show", "1ad5e8c:scripts/visualization/ectopic_recomb_plots.R"), stdout = TRUE, stderr = TRUE)
if (!length(code) || any(grepl("fatal:", code, fixed = TRUE))) stop("Could not recover historical ectopic source")
code <- sub('five_domain_path <- ".*"', paste0("five_domain_path <- ", deparse(staged_input)), code)
code <- sub('output_dir <- ".*"', paste0("output_dir <- ", deparse(figure_dir)), code)
code <- sub('dir.create\\(file.path\\("Results", "Stats"\\), recursive = TRUE, showWarnings = FALSE\\)', paste0("dir.create(", deparse(stats_dir), ", recursive = TRUE, showWarnings = FALSE)"), code)
code <- sub('sink\\(file.path\\("Results", "Stats", "anova_results.txt"\\)\\)', paste0("sink(file.path(", deparse(stats_dir), ", \"anova_results.txt\"))"), code)
eval(parse(text = paste(code, collapse = "\n")), envir = new.env(parent = globalenv()))
cat("Historical Talk IV ectopic figure restored in ", figure_dir, "\n", sep = "")
