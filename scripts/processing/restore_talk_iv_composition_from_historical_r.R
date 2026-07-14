#!/usr/bin/env Rscript

# Execute the original order/superfamily mean-composition ggplot sections from
# the Talk IV R source. Input paths are the only replacement.

args <- commandArgs(trailingOnly = FALSE)
flag <- grep("^--file=", args, value = TRUE)
script_dir <- if (length(flag)) dirname(normalizePath(sub("^--file=", "", flag[[1]]), mustWork = FALSE)) else "scripts/processing"
root <- normalizePath(file.path(script_dir, "..", ".."), mustWork = TRUE)
data_dir <- file.path(root, "results", "data", "legacy_descriptive")
plot_dir <- file.path(root, "results", "figures", "talk_iv_restored_composition")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(ggplot2) })
Order <- read.csv(file.path(data_dir, "te_order_composition_te34_legacy_descriptive_v1.csv"), check.names = FALSE)
Superfamily <- read.csv(file.path(data_dir, "te_superfamily_composition_te34_legacy_descriptive_v1.csv"), check.names = FALSE)
names(Order)[names(Order) == "species"] <- "Species"
names(Superfamily)[names(Superfamily) == "species"] <- "Species"

code <- system2("git", c("show", "1ad5e8c:scripts/visualization/te_landscape_plots.R"), stdout = TRUE, stderr = TRUE)
start <- grep("^### Order Breakdown Plot mean", code)
end <- grep("^### Classified vs Unclassified", code)
if (length(start) != 1L || length(end) != 1L) stop("Could not isolate historical composition plot sections")
section <- code[start:(end - 1L)]

# These declarations are copied directly from the historical script's setup.
Retro_orders <- c("DIRS", "LINE", "LTR", "SINE")
DNA_orders <- c("Helitron", "Maverick", "PLE", "TIR", "YR")
Retro_superfamilies <- c("Gypsy", "Jockey", "L1", "Penelope", "tRNA", "DIRS")
DNA_superfamilies <- c("Academ", "CACTA", "Chapaev", "Cyrypton", "Dada", "EnSpm", "Ginger", "Helitron", "MULE", "Maverick", "Mutator", "P", "PIF.Harbinger", "PiggyBac", "Tc1.mariner", "Unknown.TIR", "hAT")
eval(parse(text = paste(section, collapse = "\n")), envir = environment())
cat("Historical Talk IV composition figures restored in ", plot_dir, "\n", sep = "")
