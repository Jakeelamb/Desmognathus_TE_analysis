#!/usr/bin/env Rscript

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

source(file.path(dirname(script_dir), "R", "path_config_utils.R"))
prefer_active_conda_r_library()

suppressPackageStartupMessages({
  library(ape)
  library(phytools)
})

project_root <- find_project_root(script_dir)
config <- load_project_config(project_root)

phylo_dir <- resolve_config_path(project_root, config$input_data$phylogeny %||% config$data$phylogeny, "input_data/phylogeny")
results_data_dir <- resolve_config_path(project_root, config$results$data, "results/data")
results_phylo_dir <- resolve_config_path(project_root, config$results$phylogeny, "results/phylogeny")
figure_dir <- resolve_config_path(project_root, config$results$figures$phylogeny, "results/figures/phylogeny")

dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

pick_first_existing <- function(paths) {
  existing <- paths[file.exists(paths)]
  if (length(existing) == 0) {
    stop("No tree file found. Checked:\n", paste(paths, collapse = "\n"), call. = FALSE)
  }
  existing[[1]]
}

tree_file <- pick_first_existing(c(
  file.path(results_data_dir, "desmo900dated_test_cleaned_phylo.tre"),
  file.path(phylo_dir, "desmo900dated_test.tre"),
  file.path(results_phylo_dir, "processed_phylogeny.nwk")
))

message("Loading phylogenetic tree from: ", tree_file)
tree <- read.tree(tree_file)

create_phylogeny_plot <- function(phy, filename, title) {
  max_age <- max(phytools::nodeHeights(phy)[, 2])
  age_breaks <- seq(0, ceiling(max_age), by = 5)

  png(filename, width = 12, height = 12, units = "in", res = 300)
  par(mar = c(5, 4, 4, 8))
  plot.phylo(
    phy,
    type = "phylogram",
    direction = "right",
    show.tip.label = TRUE,
    cex = 0.8,
    label.offset = 0.5,
    main = title,
    edge.width = 2
  )
  axis(1, at = age_breaks, labels = paste(age_breaks, "MYA"), las = 1)
  abline(v = age_breaks, lty = 2, col = "gray90")
  dev.off()
}

output_file <- file.path(figure_dir, "rectangular_phylogeny.png")
message("Creating phylogenetic tree visualization")
create_phylogeny_plot(tree, output_file, "Desmognathus Phylogeny (Rectangular)")
message("Saved simple phylogeny plot to: ", output_file)
