#!/usr/bin/env Rscript
#
# TE_PHYLO_LANDSCAPE.R
# 
# Description: Creates visualizations that integrate Transposable Element landscape data
# with a phylogenetic tree of Desmognathus species.
#
# Input: 
#   - Phylogenetic tree file: data/raw/Phylogeny/desmo900dated_test.tre
#   - Lookup table: data/raw/lookup/lookup_table.txt
#   - TE landscape files: results/landscapes/*.csv
#
# Output: Plot files (PNG) saved to results/figures/phylo_landscape/
#

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

# Load necessary libraries (assumes these are already installed)
suppressPackageStartupMessages({
  library(ape)
  library(ggtree)
  library(tidyverse)
  library(ggplot2)
  library(phytools)
  library(gridExtra)
  library(viridis)
  library(yaml)
  library(scales)
})

if (!exists("is.waive", mode = "function")) {
  is.waive <- function(x) inherits(x, "waiver")
}

project_root <- find_project_root(script_dir)
config <- load_project_config(project_root)

landscape_dir <- resolve_config_path(project_root, config$results$landscapes, "results/landscapes")
plot_dir <- resolve_config_path(project_root, config$results$figures$phylo_landscape, "results/figures/phylo_landscape")
phylo_dir <- resolve_config_path(project_root, config$input_data$phylogeny %||% config$data$phylogeny, "input_data/phylogeny")
results_data_dir <- resolve_config_path(project_root, config$results$data, "results/data")
lookup_file <- resolve_config_path(project_root, config$input_data$lookup_table %||% config$data$lookup, "input_data/lookup_table.txt")

tree_candidates <- c(
  file.path(phylo_dir, "desmo900dated_test.tre"),
  file.path(results_data_dir, "desmo900dated_test_cleaned_phylo.tre")
)
tree_file <- tree_candidates[file.exists(tree_candidates)][1]
if (is.na(tree_file)) {
  stop("No phylogeny file found for TE phylo landscape workflow.", call. = FALSE)
}

# Create output directory for plots
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

cat("Loading phylogenetic tree from:", tree_file, "\n")

# Read the tree file
tree <- read.tree(tree_file)

# Read the lookup table
cat("Loading species lookup table from:", lookup_file, "\n")
lookup_table <- read.delim(lookup_file)

# Process species names in the lookup table to match the tree format
lookup_table$TreeName <- tolower(gsub("^D\\.", "", lookup_table$Species))

# Process tree tips to remove suffixes (e.g., "fuscus_A" -> "fuscus")
original_tip_labels <- tree$tip.label
processed_tip_labels <- sapply(tree$tip.label, function(x) {
  # Remove everything after underscore or use the full name if no underscore
  ifelse(grepl("_", x), sub("_.*$", "", x), x)
})

# Create a mapping between original and processed tip labels
tip_mapping <- data.frame(
  original = original_tip_labels,
  processed = processed_tip_labels,
  stringsAsFactors = FALSE
)

# Update tree tip labels
tree$tip.label <- processed_tip_labels

# Identify species to keep (those in the lookup table)
species_to_keep <- tree$tip.label %in% lookup_table$TreeName

# If there are species to prune based on lookup table, do so
if (!all(species_to_keep)) {
  tips_to_drop_lookup <- original_tip_labels[!species_to_keep]
  cat("Pruning", length(tips_to_drop_lookup), "species not found in lookup table:\n")
  cat(paste(" -", tips_to_drop_lookup), sep = "\n")
  
  # Prune the tree using original labels before they were processed
  tree <- drop.tip(tree, tips_to_drop_lookup)
}

# === Collapse Redundant Species Tips ===
cat("\nChecking for and collapsing redundant species tips...\n")

# Get current tip labels after lookup-based pruning
current_tips <- tree$tip.label

# Find duplicated labels (TRUE for second, third, etc. occurrences)
duplicated_mask <- duplicated(current_tips)

# If there are duplicates
if (any(duplicated_mask)) {
  duplicated_labels <- unique(current_tips[duplicated_mask])
  cat("Found redundant tips for species:", paste(duplicated_labels, collapse=", "), "\n")
  
  # Get the indices of tips to drop (all duplicated occurrences)
  indices_to_drop <- which(duplicated_mask)
  
  cat("Dropping", length(indices_to_drop), "redundant tips.\n")
  
  # Drop the redundant tips
  tree <- drop.tip(tree, indices_to_drop)
  
  # Update current tips list
  current_tips <- tree$tip.label
} else {
  cat("No redundant species tips found after initial pruning.\n")
}
# === End Collapsing Redundant Tips ===

# Rebuild the mapping based on the final tree
cat("\nRebuilding species mapping for the final tree...\n")
tree_to_srx_map <- data.frame(
  tree_label = tree$tip.label, # Use final tip labels
  stringsAsFactors = FALSE
)

# Add SRX ID by matching with lookup table
tree_to_srx_map$srx_id <- sapply(tree_to_srx_map$tree_label, function(label) {
  match_row <- which(lookup_table$TreeName == label)
  if (length(match_row) > 0) {
    return(lookup_table$SRA_Accension[match_row[1]])
  } else {
    return(NA)
  }
})

# Create full species name for display
tree_to_srx_map$display_name <- sapply(tree_to_srx_map$tree_label, function(label) {
  match_row <- which(lookup_table$TreeName == label)
  if (length(match_row) > 0) {
    return(lookup_table$Species[match_row[1]])
  } else {
    # Fallback if somehow still not in lookup (shouldn't happen after pruning)
    return(paste0("D. ", label)) 
  }
})

# Remove rows with NA SRX IDs (shouldn't be any after pruning, but good practice)
tree_to_srx_map <- tree_to_srx_map[!is.na(tree_to_srx_map$srx_id), ]

# Ensure tree only contains tips present in the final map
final_tips_to_keep <- tree$tip.label %in% tree_to_srx_map$tree_label
if (!all(final_tips_to_keep)) {
  final_tips_to_drop <- tree$tip.label[!final_tips_to_keep]
  cat("Performing final prune for tips missing SRX info:", paste(final_tips_to_drop, collapse=", "), "\n")
  tree <- drop.tip(tree, final_tips_to_drop)
}

# Load landscape data for species in the final tree
cat("\nLoading landscape data for", nrow(tree_to_srx_map), "final species...\n")

landscape_data <- list()
landscape_files <- list.files(landscape_dir, pattern = "^repeat_landscape.*\\.csv$", full.names = TRUE)

for (i in 1:nrow(tree_to_srx_map)) {
  srx_id <- tree_to_srx_map$srx_id[i]
  species_name <- tree_to_srx_map$tree_label[i]
  
  # Find the corresponding landscape file
  landscape_file <- landscape_files[grep(srx_id, landscape_files)]
  
  if (length(landscape_file) > 0 && file.exists(landscape_file)) {
    # Read landscape data
    temp_data <- read.csv(landscape_file)
    temp_data$Species <- species_name
    landscape_data[[i]] <- temp_data
    cat("  Loaded", species_name, "from", basename(landscape_file), "\n")
  } else {
    cat("  No landscape file found for", species_name, "(", srx_id, ")\n")
  }
}

# Combine all landscape data
if (length(landscape_data) > 0) {
  combined_data <- do.call(rbind, landscape_data)
  
  # Create a mapping of species to distinct colors
  species_colors <- setNames(
    viridis(length(unique(combined_data$Species)), option = "D"),
    unique(combined_data$Species)
  )
  
  # Calculate summary statistics for heatmap
  heatmap_data <- combined_data %>%
    # Filter out bins with very low counts to reduce noise
    group_by(Kimura_bin) %>%
    mutate(total_bp = sum(aligned_bp)) %>%
    ungroup() %>%
    filter(total_bp > quantile(total_bp, 0.05, na.rm = TRUE)) %>%  # Remove lowest 5% of bins, handle NAs
    # Recalculate the binned data
    group_by(Species, Kimura_bin) %>%
    summarize(aligned_bp = sum(aligned_bp), .groups = "drop") %>%
    ungroup()

  # Create matrix for heatmap
  heatmap_matrix_unordered <- heatmap_data %>%
    pivot_wider(
      id_cols = Species,
      names_from = Kimura_bin,
      values_from = aligned_bp,
      values_fill = 0
    ) %>%
    column_to_rownames("Species") %>%
    as.matrix()

  # Get species in the tree in phylogenetic order
  tree_order <- tree$tip.label

  # Ensure matrix only contains species present in the final tree order
  heatmap_matrix_filtered <- heatmap_matrix_unordered[rownames(heatmap_matrix_unordered) %in% tree_order, , drop = FALSE]

  # Reorder the matrix rows to match the tree tip order
  # Use match() to handle potential missing species gracefully
  final_row_order <- match(tree_order, rownames(heatmap_matrix_filtered))
  heatmap_matrix <- heatmap_matrix_filtered[final_row_order, , drop = FALSE]

  # Check if matrix is valid before proceeding
  if (is.null(heatmap_matrix) || nrow(heatmap_matrix) == 0 || ncol(heatmap_matrix) == 0) {
      stop("Heatmap matrix is invalid after ordering. Check species matching.")
  }

  # Create display name mapping
  display_names <- sapply(tree_order, function(label) {
    idx <- which(tree_to_srx_map$tree_label == label)
    if (length(idx) > 0) {
      return(tree_to_srx_map$display_name[idx])
    } else {
      return(paste0("D. ", label))
    }
  })
  
  # Enhance tree for better visualization
  # Scale tree to make it more readable (shorter branches)
  scaled_tree <- tree
  scaled_tree$edge.length <- scaled_tree$edge.length * 0.8  # Scale branches to 80% of original length
  
  #----------------------------------------
  # Create a basic tree visualization with time scale
  #----------------------------------------
  cat("Creating basic phylogenetic tree visualization...\n")
  
  # Calculate time scale: assuming tree is in millions of years
  max_age <- max(nodeHeights(tree)[,2])
  age_breaks <- seq(0, ceiling(max_age), by = 5)
  
  # Plot the basic tree with species names and time scale
  # Create a mapping for the basic plot as well
  basic_label_map <- setNames(
    sapply(scaled_tree$tip.label, function(l) {
      idx <- which(tree_to_srx_map$tree_label == l)
      if (length(idx) > 0) {
        return(tree_to_srx_map$display_name[idx])
      } else {
        return(paste0("D. ", l))
      }
    }),
    scaled_tree$tip.label
  )

  ggtree_obj <- ggtree(scaled_tree) +
    geom_tiplab(aes(label = ifelse(isTip, basic_label_map[label], NA_character_)), # Conditional label mapping
                size = 3, hjust = -0.1, na.rm = TRUE) + # Add na.rm=TRUE
    geom_vline(xintercept = age_breaks, linetype = "dashed", alpha = 0.3) +
    theme_tree2() +
    labs(title = "Desmognathus Phylogeny", 
         subtitle = "With time scale (MYA)") +
    theme(plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5)) +
    xlim(0, max(nodeHeights(tree)[,2]) * 1.5)  # Extend x-axis for labels
  
  # Add time scale
  for (age in age_breaks) {
    ggtree_obj <- ggtree_obj + 
      annotate("text", x = age, y = 0, label = age, size = 3, vjust = -0.5)
  }

  tip_positions <- ggtree(scaled_tree)$data %>%
    filter(isTip) %>%
    select(label, y)

  add_tip_heatmap <- function(base_plot,
                              matrix_data,
                              fill_name,
                              fill_option = "viridis",
                              log_fill = FALSE,
                              label_angle = 90,
                              total_width = max_age * 0.8) {
    if (is.null(matrix_data) || nrow(matrix_data) == 0 || ncol(matrix_data) == 0) {
      stop("Heatmap matrix is empty.", call. = FALSE)
    }

    heatmap_start <- max_age * 1.25
    cell_width <- total_width / max(ncol(matrix_data), 1)
    column_lookup <- tibble(
      column = colnames(matrix_data),
      x = heatmap_start + (seq_along(colnames(matrix_data)) - 0.5) * cell_width
    )

    heatmap_long <- as.data.frame(matrix_data, check.names = FALSE) %>%
      tibble::rownames_to_column("label") %>%
      pivot_longer(cols = -label, names_to = "column", values_to = "value") %>%
      left_join(tip_positions, by = "label") %>%
      left_join(column_lookup, by = "column") %>%
      mutate(plot_value = if (log_fill) pmax(value, 1) else value)

    heatmap_end <- heatmap_start + total_width + cell_width
    label_y <- max(tip_positions$y, na.rm = TRUE) + 1

    plot <- base_plot +
      geom_tile(
        data = heatmap_long,
        aes(x = x, y = y, fill = plot_value),
        inherit.aes = FALSE,
        width = cell_width * 0.95,
        height = 0.9
      ) +
      geom_text(
        data = column_lookup,
        aes(x = x, y = label_y, label = column),
        inherit.aes = FALSE,
        angle = label_angle,
        size = 2.5,
        hjust = 1
      ) +
      theme(
        legend.position = "right",
        legend.key.size = unit(1, "cm")
      ) +
      coord_cartesian(xlim = c(0, heatmap_end), clip = "off")

    if (log_fill) {
      plot <- plot +
        scale_fill_viridis_c(
          name = fill_name,
          option = fill_option,
          trans = "log10",
          labels = scales::label_number(accuracy = 1, scale_cut = scales::cut_short_scale()),
          limits = c(1, NA),
          oob = scales::squish
        )
    } else {
      plot <- plot +
        scale_fill_viridis_c(name = fill_name, option = fill_option)
    }

    plot
  }
  
  #----------------------------------------
  # Create tree with TE landscape heatmap
  #----------------------------------------
  cat("Creating phylogeny with TE landscape heatmap...\n")
  
  # Create the phylogeny with heatmap
  phylo_landscape_base <- ggtree(scaled_tree) +
    geom_tiplab(aes(label = ifelse(isTip, basic_label_map[label], NA_character_)), 
                size = 3, hjust = -0.1, na.rm = TRUE) +
    theme_tree2() +
    xlim(0, max_age * 2.3)

  phylo_landscape <- add_tip_heatmap(
    phylo_landscape_base,
    heatmap_matrix,
    fill_name = "Aligned Base Pairs",
    fill_option = "magma",
    log_fill = TRUE,
    label_angle = 90,
    total_width = max_age * 0.95
  )
  
  #----------------------------------------
  # Create tree with TE class distribution
  #----------------------------------------
  cat("Creating phylogeny with TE class distribution...\n")
  
  # Calculate class distribution
  class_data <- combined_data %>%
    filter(!is.na(Class) & Class != "") %>%
    group_by(Species, Class) %>%
    summarize(aligned_bp = sum(aligned_bp), .groups = "drop") %>%
    group_by(Species) %>%
    mutate(Percentage = aligned_bp / sum(aligned_bp) * 100) %>%
    ungroup() %>%
    filter(Species %in% tree_order)

  # Create matrix for class distribution
  class_matrix <- class_data %>%
    pivot_wider(
      id_cols = Species,
      names_from = Class,
      values_from = Percentage,
      values_fill = 0
    ) %>%
    column_to_rownames("Species") %>%
    as.matrix()
  
  # Create the phylogeny with class distribution
  phylo_class_base <- ggtree(scaled_tree) +
    geom_tiplab(aes(label = ifelse(isTip, basic_label_map[label], NA_character_)), 
                size = 3, hjust = -0.1, na.rm = TRUE) +
    theme_tree2() +
    xlim(0, max_age * 1.9)

  phylo_class <- add_tip_heatmap(
    phylo_class_base,
    class_matrix,
    fill_name = "Percentage",
    fill_option = "viridis",
    log_fill = FALSE,
    label_angle = 45,
    total_width = max_age * 0.55
  )
  
  #----------------------------------------
  # Create tree with circular layout and TE landscape
  #----------------------------------------
  # cat("Creating circular phylogeny with TE landscape data...\n")
  # 
  # # Create a mapping from original label to display name
  # label_map <- setNames(
  #   sapply(tree$tip.label, function(l) {
  #     idx <- which(tree_to_srx_map$tree_label == l)
  #     if (length(idx) > 0) {
  #       return(tree_to_srx_map$display_name[idx])
  #     } else {
  #       return(paste0("D. ", l))
  #     }
  #   }),
  #   tree$tip.label # Names of the vector are the original tip labels
  # )
  # 
  # # Create the base circular ggtree object
  # circular_phylo_base <- ggtree(scaled_tree, layout = "circular") + 
  #   theme_tree()
  # 
  # # Extract data, filter for tips, and add label/angle/hjust info
  # tip_plot_data <- circular_phylo_base$data %>% 
  #   filter(isTip) %>%
  #   mutate(
  #     display_name = label_map[label],
  #     angle = (atan2(y, x) * 180/pi) + 90,
  #     angle = ifelse(angle > 90 & angle < 270, angle + 180, angle),
  #     hjust = ifelse(angle > 90 & angle < 270, 1, 0)
  #   )
  # 
  # # Add labels using geom_text with the explicitly prepared tip data
  # circular_phylo <- circular_phylo_base + 
  #   geom_text(data = tip_plot_data, 
  #             aes(x = x * 1.1, y = y, label = display_name, angle = angle, hjust = hjust),
  #             size = 3,
  #             inherit.aes = FALSE) + # Explicitly prevent inheritance
  #   labs(title = "Desmognathus Phylogeny", 
  #        subtitle = "Circular layout showing evolutionary relationships") +
  #   theme(plot.title = element_text(hjust = 0.5),
  #         plot.subtitle = element_text(hjust = 0.5))
  # 
  # # Calculate maximum tree radius from the prepared tip data
  # if (nrow(tip_plot_data) > 0) {
  #   max_radius <- max(tip_plot_data$x, na.rm = TRUE)
  # } else {
  #   max_radius <- 1 # Default fallback
  # }
  # 
  # # Add time scales as background circles
  # time_rings <- age_breaks[age_breaks > 0]
  # for (time in time_rings) {
  #   if (is.finite(max_radius) && time <= max_radius) {
  #     circular_phylo <- circular_phylo +
  #       annotate("path", 
  #               x = time * cos(seq(0, 2*pi, length.out = 100)),
  #               y = time * sin(seq(0, 2*pi, length.out = 100)),
  #               colour = "gray80", linetype = "dashed", linewidth = 0.3)
  #   }
  # }
  # 
  # # Add time labels at selected positions
  # for (time in time_rings) {
  #   if (is.finite(max_radius) && time <= max_radius) {
  #     circular_phylo <- circular_phylo +
  #       annotate("text", x = time, y = 0, 
  #               label = paste(time, "MYA"), size = 3,
  #               color = "gray40")
  #   }
  # }
  # 
  # # Expand plot margins to make room for labels
  # circular_phylo <- circular_phylo +
  #   coord_fixed(clip = "off") +
  #   expand_limits(x = c(-max_radius * 1.2, max_radius * 1.2), 
  #                 y = c(-max_radius * 1.2, max_radius * 1.2))
  
  # Save basic phylogenetic tree
  ggsave(
    filename = file.path(plot_dir, "basic_phylogeny.png"),
    plot = ggtree_obj,
    width = 12,
    height = 10,
    dpi = 300
  )

  # Save phylogeny with TE landscape heatmap
  ggsave(
    filename = file.path(plot_dir, "phylogeny_with_landscape.png"),
    plot = phylo_landscape,
    width = 16,
    height = 12,
    dpi = 300
  )

  # Save phylogeny with class distribution
  ggsave(
    filename = file.path(plot_dir, "phylogeny_with_classes.png"),
    plot = phylo_class,
    width = 14,
    height = 12,
    dpi = 300
  )

  # # Save circular phylogeny
  # ggsave(
  #   filename = file.path(plot_dir, "circular_phylogeny.png"),
  #   plot = circular_phylo,
  #   width = 16,
  #   height = 16,
  #   dpi = 300
  # )

  cat("Visualizations complete. Results saved to", plot_dir, "\n")
  cat("Created 3 visualization files.\n")
} else {
  cat("No landscape data was loaded. Unable to create visualizations.\n")
} 
