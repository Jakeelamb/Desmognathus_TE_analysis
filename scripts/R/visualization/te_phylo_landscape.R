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

# Read the configuration file
config_file <- "config/paths.yaml"
if (file.exists(config_file)) {
  config <- yaml::read_yaml(config_file)
  
  # Define paths from configuration
  landscape_dir <- config$results$landscapes
  plot_dir <- config$results$figures$phylo_landscape
  tree_file <- file.path(config$data$phylogeny, "desmo900dated_test.tre")
  lookup_file <- file.path(config$data$lookup, "lookup_table.txt")
} else {
  # Fallback paths if config file is not found
  landscape_dir <- "results/landscapes"
  plot_dir <- "results/figures/phylo_landscape"
  tree_file <- "data/raw/Phylogeny/desmo900dated_test.tre"
  lookup_file <- "data/raw/lookup/lookup_table.txt"
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
  
  #----------------------------------------
  # Create tree with TE landscape heatmap
  #----------------------------------------
  cat("Creating phylogeny with TE landscape heatmap...\n")
  
  # Create the phylogeny with heatmap
  phylo_landscape <- ggtree(scaled_tree) +
    geom_tiplab(aes(label = ifelse(isTip, basic_label_map[label], NA_character_)), 
                size = 3, hjust = -0.1, na.rm = TRUE) +
    theme_tree2() +
    xlim(0, max(nodeHeights(tree)[,2]) * 1.2)  # Add space for heatmap
  
  # Add the heatmap with improved visibility
  phylo_landscape <- phylo_landscape %>%
    gheatmap(heatmap_matrix, # Use the reordered matrix
             offset = max(nodeHeights(tree)[,2]) * 0.3,  # Position heatmap further away
             width = 8,  # Make the heatmap wider
             colnames_angle = 90,  # Make column names vertical for better readability
             colnames_offset_y = -0.5,  # Adjust label position
             font.size = 4) +  # Increase font size
    scale_fill_viridis_c(
      name = "Aligned Base Pairs",
      option = "magma",
      trans = "log10", 
      # Add pseudo-count (1) before log transform to handle zeros
      labels = scales::label_number(accuracy = 1, scale_cut = scales::cut_short_scale()), # Nicer labels for log scale
      limits = c(1, NA), # Ensure scale starts at 1 (log10(1)=0)
      oob = scales::squish # Handle values slightly out of bounds
    ) +
    theme(
      axis.text.x = element_text(size = 8, angle = 90, hjust = 1),  # Improve x-axis label readability
      legend.position = "right",
      legend.key.size = unit(1, "cm")
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
  phylo_class <- ggtree(scaled_tree) +
    geom_tiplab(aes(label = ifelse(isTip, basic_label_map[label], NA_character_)), 
                size = 3, hjust = -0.1, na.rm = TRUE) +
    theme_tree2() +
    xlim(0, max(nodeHeights(tree)[,2]) * 1.2)  # Add space for heatmap
  
  # Add the class distribution heatmap
  phylo_class <- phylo_class %>%
    gheatmap(class_matrix, 
             offset = max(nodeHeights(tree)[,2]) * 0.3,  # Position heatmap further away
             width = 6,  # Make the heatmap wider
             colnames_angle = 45,
             colnames_offset_y = 0,
             font.size = 3) +
    scale_fill_viridis_c(name = "Percentage",
                         option = "viridis")
  
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