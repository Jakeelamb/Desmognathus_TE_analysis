#!/usr/bin/env Rscript

script_dir <- tryCatch({
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = FALSE))
  } else {
    "scripts/R/visualization"
  }
}, error = function(...) {
  "scripts/R/visualization"
})

source(file.path(dirname(dirname(script_dir)), "R", "path_config_utils.R"))
prefer_active_conda_r_library()

suppressPackageStartupMessages({
  library(ape)
  library(yaml)
  library(phytools)
  library(plotrix)
})

# === Configuration ===
project_root <- find_project_root(script_dir)
config <- load_project_config(project_root)

plot_dir <- resolve_config_path(project_root, config$results$figures$phylogeny, "results/figures/phylogeny")
phylo_dir <- resolve_config_path(project_root, config$results$phylogeny, "results/phylogeny")
lookup_file <- resolve_config_path(project_root, config$input_data$lookup_table %||% config$data$lookup, "input_data/lookup_table.txt")

tree_candidates <- c(
  file.path(resolve_config_path(project_root, config$input_data$phylogeny %||% config$data$phylogeny, "input_data/phylogeny"), "desmo900dated_test.tre"),
  file.path(project_root, "results", "data", "desmo900dated_test_cleaned_phylo.tre")
)
tree_file <- tree_candidates[file.exists(tree_candidates)][1]
if (is.na(tree_file)) {
  stop("No phylogeny file found for simple phylogeny workflow.", call. = FALSE)
}

dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(phylo_dir, recursive = TRUE, showWarnings = FALSE)

# === Load Data ===
cat("Loading phylogenetic tree from:", tree_file, "\n")
tree <- read.tree(tree_file)
cat("Loading species lookup table from:", lookup_file, "\n")
lookup_table <- read.delim(lookup_file)
lookup_table$TreeName <- tolower(gsub("^D\\.", "", lookup_table$Species))

# === Collapse Suffixes in Tip Labels ===
original_tip_labels <- tree$tip.label
base_names <- sapply(original_tip_labels, function(x) {
  tolower(sub("_.*$", "", x))
})

# Group tips by base name
tip_groups <- split(original_tip_labels, base_names)

# Collapse groups with >1 tip
tips_to_keep <- character()
tips_to_drop <- character()

for (group in tip_groups) {
  if (length(group) > 1) {
    tips_to_keep <- c(tips_to_keep, group[1])  # Keep first
    tips_to_drop <- c(tips_to_drop, group[-1]) # Drop others
  } else {
    tips_to_keep <- c(tips_to_keep, group)
  }
}

if (length(tips_to_drop) > 0) {
  cat("Collapsing duplicated species entries:\n")
  print(tips_to_drop)
  tree <- drop.tip(tree, tips_to_drop)
}

# Re-label remaining tips with base names
new_labels <- sapply(tree$tip.label, function(x) {
  tolower(sub("_.*$", "", x))
})
tree$tip.label <- new_labels

# === Match to Lookup Table ===
species_to_keep <- tree$tip.label %in% lookup_table$TreeName
if (!all(species_to_keep)) {
  missing <- tree$tip.label[!species_to_keep]
  cat("Dropping", length(missing), "tips not found in lookup table:\n")
  print(missing)
  tree <- drop.tip(tree, missing)
}

# === Final Display Names ===
display_names <- sapply(tree$tip.label, function(label) {
  match_row <- which(lookup_table$TreeName == label)
  if (length(match_row) > 0) {
    lookup_table$Species[match_row[1]]
  } else {
    paste0("D. ", label)
  }
})
tree$tip.label <- display_names
scaled_tree <- tree
scaled_tree$edge.length <- scaled_tree$edge.length * 0.8

# === Save Processed Tree ===
output_newick_file <- file.path(phylo_dir, "processed_phylogeny.nwk")
cat("Saving processed Newick tree to:", output_newick_file, "\n")
write.tree(scaled_tree, file = output_newick_file)

# === Monophyly Functions ===
identify_monophyletic_groups <- function(tree) {
  n_tips <- length(tree$tip.label)
  n_nodes <- tree$Nnode
  nodes <- (n_tips + 1):(n_tips + n_nodes)
  monophyletic_nodes <- c()
  for (node in nodes) {
    desc <- tree$tip.label[getDescendants(tree, node)]
    if (length(desc) >= 2) {
      stripped <- gsub("^D\\. ", "", desc)
      stripped <- gsub(" .*$", "", stripped)
      if (length(unique(stripped)) == 1) {
        monophyletic_nodes <- c(monophyletic_nodes, node)
      }
    }
  }
  return(monophyletic_nodes)
}

create_node_colors <- function(tree, monophyletic_nodes) {
  edge_colors <- rep("black", nrow(tree$edge))
  colors <- colorRampPalette(c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",
                               "#FFFF33", "#A65628", "#F781BF", "#999999"))(length(monophyletic_nodes))
  for (i in seq_along(monophyletic_nodes)) {
    node <- monophyletic_nodes[i]
    clade_edges <- which(tree$edge[,2] %in% getDescendants(tree, node))
    parent_edge <- which(tree$edge[,2] == node)
    edge_colors[c(clade_edges, parent_edge)] <- colors[i]
  }
  return(edge_colors)
}

create_phylogeny_plot <- function(tree, type = "rectangular", filename) {
  mono_nodes <- identify_monophyletic_groups(tree)
  edge_colors <- create_node_colors(tree, mono_nodes)
  max_age <- max(nodeHeights(tree)[,2])
  age_breaks <- seq(0, ceiling(max_age), by = 5)
  png(filename, width = 12, height = 12, units = "in", res = 300)
  if (type == "circular") {
    par(mar = c(2, 2, 2, 2))
    plot.phylo(tree, type = "fan", show.tip.label = TRUE, cex = 0.7,
               edge.color = edge_colors, main = "Desmognathus Phylogeny (Circular)",
               edge.width = 2, label.offset = 1)
    for (i in seq_along(age_breaks)) {
      draw.circle(0, 0, radius = age_breaks[i], 
                 lty = 2,)
    }
    for (i in seq_along(age_breaks)) {
      age <- age_breaks[i]
      text(age * cos(pi/4), age * sin(pi/4), paste(age, "MYA"), cex = 0.8,
           adj = c(-0.2, 0.5))
    }
  } else {
    par(mar = c(5, 4, 4, 8))
    plot.phylo(tree, type = "phylogram", direction = "right", show.tip.label = TRUE,
               cex = 0.8, edge.color = edge_colors, label.offset = 0.5,
               main = "Desmognathus Phylogeny", edge.width = 2)
    axis(1, at = age_breaks, labels = paste(age_breaks, "MYA"), las = 1)
    abline(v = age_breaks, lty = 2, col = "gray90")
  }
  if (length(mono_nodes) > 0) {
    group_labels <- sapply(mono_nodes, function(node) {
      tips <- tree$tip.label[getDescendants(tree, node)]
      species <- unique(gsub("^D\\. |\\s.*$", "", tips))[1]
      paste("D.", species, "group")
    })
    legend_pos <- if(type == "circular") "topright" else "right"
    legend(legend_pos, legend = group_labels,
           col = unique(edge_colors[edge_colors != "black"]), lwd = 2,
           cex = 0.7, bg = "white", box.col = "gray80", title = "Monophyletic Groups")
  }
  dev.off()
}

# === Generate Plots ===
cat("Creating phylogenetic tree visualizations...\n")
create_phylogeny_plot(scaled_tree, type = "rectangular",
                      filename = file.path(plot_dir, "rectangular_phylogeny.png"))
create_phylogeny_plot(scaled_tree, type = "circular",
                      filename = file.path(plot_dir, "circular_phylogeny.png"))
cat("Done!\n")
