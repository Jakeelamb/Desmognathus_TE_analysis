#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ape)
  library(yaml)
  library(phytools)
  library(plotrix)
})

# === Configuration ===
tree_file <- file.path("results/data/desmo900dated_test_cleaned.tre")
cat("Loading phylogenetic tree from:", tree_file, "\n")
tree <- read.tree(tree_file)

create_phylogeny_plot <- function(tree, type = "rectangular", filename) {
  # mono_nodes <- identify_monophyletic_groups(tree) # Removed
  # edge_colors <- create_node_colors(tree, mono_nodes) # Removed
  max_age <- max(nodeHeights(tree)[,2])
  age_breaks <- seq(0, ceiling(max_age), by = 5)
  png(filename, width = 12, height = 12, units = "in", res = 300)
  if (type == "rectangular") {
    par(mar = c(5, 4, 4, 8))
    plot.phylo(tree, type = "phylogram", direction = "right", show.tip.label = TRUE,
               cex = 0.8, label.offset = 0.5, # Removed edge.color = edge_colors
               main = "Desmognathus Phylogeny (Rectangular)", edge.width = 2)
    axis(1, at = age_breaks, labels = paste(age_breaks, "MYA"), las = 1)
    abline(v = age_breaks, lty = 2, col = "gray90")
  } else {
    par(mar = c(5, 4, 4, 8))
    plot.phylo(tree, type = "phylogram", direction = "right", show.tip.label = TRUE,
               cex = 0.8, label.offset = 0.5, # Removed edge.color = edge_colors
               main = "Desmognathus Phylogeny", edge.width = 2)
    axis(1, at = age_breaks, labels = paste(age_breaks, "MYA"), las = 1)
    abline(v = age_breaks, lty = 2, col = "gray90")
  }
  # Removed the entire if block for the legend
  # if (length(mono_nodes) > 0) {
  #   ... legend code ...
  # }
  dev.off()
}

# === Generate Plots ===
cat("Creating phylogenetic tree visualizations...\n")
create_phylogeny_plot(tree, type = "rectangular",
                      filename = file.path("results/figures/", "rectangular_phylogeny.png"))
cat("Done!\n")
