#!/usr/bin/env Rscript

script_dir <- tryCatch({
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = FALSE))
  } else {
    "scripts/processing"
  }
}, error = function(...) {
  "scripts/processing"
})

source(file.path(dirname(script_dir), "R", "path_config_utils.R"))
prefer_active_conda_r_library()

suppressPackageStartupMessages({
  library(ape)
  library(dplyr)
  library(ggtree)
  library(ggplot2)
  library(RColorBrewer)
  library(readr)
  library(yaml)
})

project_root <- find_project_root(script_dir)
config <- load_project_config(project_root)

phylo_dir <- resolve_config_path(project_root, config$input_data$phylogeny %||% config$data$phylogeny, "input_data/phylogeny")
results_data_dir <- resolve_config_path(project_root, config$results$data, "results/data")
results_phylo_dir <- resolve_config_path(project_root, config$results$phylogeny, "results/phylogeny")
output_dir <- resolve_config_path(project_root, config$results$tables$phylogeny, "results/tables/phylogeny")
figure_dir <- resolve_config_path(project_root, config$results$figures$phylogeny, "results/figures/phylogeny")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

pick_first_existing <- function(paths) {
  existing <- paths[file.exists(paths)]
  if (length(existing) == 0) {
    stop(
      "No tree file found. Checked:\n",
      paste(paths, collapse = "\n"),
      call. = FALSE
    )
  }
  existing[[1]]
}

tree_candidates <- c(
  file.path(phylo_dir, "desmo900dated_test.tre"),
  file.path(results_data_dir, "desmo900dated_test_cleaned_phylo.tre"),
  file.path(results_phylo_dir, "processed_phylogeny.nwk")
)
tree_path <- pick_first_existing(tree_candidates)

message("Project root: ", project_root)
message("Loading tree from: ", tree_path)
tree <- read.tree(tree_path)

threshold_mya <- 5
message("Defining clades based on ", threshold_mya, " MYA threshold")

internal_nodes <- (length(tree$tip.label) + 1):(length(tree$tip.label) + tree$Nnode)
node_depths <- node.depth.edgelength(tree)
root_height <- max(node_depths)
node_heights <- root_height - node_depths

get_descendants <- function(node, phy) {
  if (node <= length(phy$tip.label)) {
    return(phy$tip.label[[node]])
  }
  extract.clade(phy, node)$tip.label
}

get_node_height <- function(node, heights, phy) {
  if (node <= length(phy$tip.label)) {
    return(0)
  }
  heights[[node]]
}

all_species <- tibble::tibble(
  species = tree$tip.label,
  clade = NA_character_
)

clade_counter <- 1
node_clade_map <- tibble::tibble(node = integer(), clade = character(), n_tips = integer())

edge_df <- tibble::tibble(
  parent = tree$edge[, 1],
  child = tree$edge[, 2]
) %>%
  mutate(
    parent_height = vapply(parent, get_node_height, numeric(1), heights = node_heights, phy = tree),
    child_height = vapply(child, get_node_height, numeric(1), heights = node_heights, phy = tree)
  ) %>%
  filter(parent_height > threshold_mya, child_height <= threshold_mya)

if (nrow(edge_df) == 0) {
  all_species$clade <- "Clade_1"
} else {
  for (i in seq_len(nrow(edge_df))) {
    child_node <- edge_df$child[[i]]
    tips <- get_descendants(child_node, tree)
    clade_name <- paste0("Clade_", clade_counter)
    all_species$clade[all_species$species %in% tips] <- clade_name

    if (child_node > length(tree$tip.label)) {
      node_clade_map <- dplyr::bind_rows(
        node_clade_map,
        tibble::tibble(node = child_node, clade = clade_name, n_tips = length(tips))
      )
    }

    clade_counter <- clade_counter + 1
  }
}

unassigned <- all_species$species[is.na(all_species$clade)]
if (length(unassigned) > 0) {
  for (sp in unassigned) {
    all_species$clade[all_species$species == sp] <- paste0("Clade_", clade_counter)
    clade_counter <- clade_counter + 1
  }
}

clade_assignments <- all_species %>%
  arrange(species)

message("\nClade assignments:")
print(clade_assignments)

n_clades <- dplyr::n_distinct(clade_assignments$clade)
if (n_clades == 1) {
  clade_colors <- c("Clade_1" = "#1b9e77")
} else {
  palette_size <- max(3, min(n_clades, 8))
  palette_values <- brewer.pal(palette_size, "Dark2")
  clade_colors <- palette_values[seq_len(n_clades)]
  names(clade_colors) <- sort(unique(clade_assignments$clade))
}

tip_metadata <- clade_assignments %>%
  mutate(clade_num = sub("^Clade_", "", clade)) %>%
  transmute(
    label = species,
    clade,
    label_with_clade = paste0(species, " (", clade_num, ")")
  )

p <- ggtree(tree, layout = "rectangular")
tree_x_max <- max(p$data$x, na.rm = TRUE)
label_offset <- max(1, tree_x_max * 0.25)

p$data <- p$data %>%
  left_join(tip_metadata, by = "label") %>%
  mutate(label_with_clade = dplyr::coalesce(label_with_clade, label))

p <- p +
  geom_tippoint(aes(color = clade), size = 2, na.rm = TRUE) +
  geom_tiplab(aes(label = label_with_clade), size = 3, hjust = -0.1) +
  scale_color_manual(values = clade_colors, guide = "none", na.value = "grey60") +
  theme_tree2() +
  theme(
    axis.text.x = element_text(size = 10),
    axis.title.x = element_text(size = 12, face = "bold"),
    plot.margin = margin(20, 180, 20, 20)
  ) +
  xlim(NA, tree_x_max + label_offset) +
  xlab("Time (MYA)")

if (nrow(node_clade_map) > 0) {
  for (i in seq_len(nrow(node_clade_map))) {
    clade_name <- node_clade_map$clade[[i]]
    p <- p + geom_cladelabel(
      node = node_clade_map$node[[i]],
      label = gsub("_", " ", clade_name),
      color = clade_colors[[clade_name]],
      offset.text = label_offset * 0.85,
      barsize = 0.8,
      fontsize = 3,
      hjust = 0,
      angle = 0,
      align = TRUE,
      offset = 0.3
    )
  }
}

output_plot_file <- file.path(figure_dir, "clade_phylogeny.png")
ggsave(output_plot_file, p, width = 15, height = 8, dpi = 300, bg = "white")
message("Saved clade figure to: ", output_plot_file)

output_csv_file <- file.path(output_dir, "clade_assignments.csv")
readr::write_csv(clade_assignments, output_csv_file)
message("Saved clade assignments to: ", output_csv_file)
