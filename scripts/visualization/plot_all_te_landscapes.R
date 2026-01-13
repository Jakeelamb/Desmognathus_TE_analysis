#!/usr/bin/env Rscript
#
# visualize_all_landscapes.R
#
# Description: Generates visualizations of TE order distribution from multiple samples.
# The combined line plot and heatmap have been removed as requested.
#
# Usage: Rscript visualize_all_landscapes.R
#

# Load required libraries
suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(viridis)
  library(reshape2)
  library(gridExtra)
  library(RColorBrewer)
})

# Define paths
landscapes_dir <- "results/landscapes"
output_dir <- "results/figures"
lookup_table_path <- "data/raw/lookup/lookup_table.txt"

# Create output directory if it doesn't exist
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Load the species lookup table
load_species_lookup <- function(lookup_path = lookup_table_path) {
  if (!file.exists(lookup_path)) {
    cat(sprintf("Warning: Lookup table not found at %s\n", lookup_path))
    return(NULL)
  }
  
  # Read the lookup table (tab-separated)
  lookup_df <- read.table(lookup_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  
  # Create a named vector for easy lookup
  species_map <- lookup_df$Species
  names(species_map) <- lookup_df$SRA_Accension
  
  cat(sprintf("Loaded lookup table with %d species mappings\n", length(species_map)))
  return(species_map)
}

# Function to get species name from SRX ID
get_species_name <- function(srx_id, species_map) {
  if (is.null(species_map) || !(srx_id %in% names(species_map))) {
    return(srx_id)  # Return original ID if no mapping found
  }
  return(species_map[srx_id])
}

# Function to load all landscape files
load_all_landscapes <- function(dir = landscapes_dir, species_map = NULL) {
  cat("Loading landscape data...\n")
  
  # List all landscape files
  landscape_files <- list.files(dir, pattern = "repeat_landscape_.+\\.csv$", full.names = TRUE)
  
  if (length(landscape_files) == 0) {
    stop("No landscape files found in directory: ", dir)
  }
  
  # Extract sample IDs from filenames
  sample_ids <- basename(landscape_files) %>%
    gsub("repeat_landscape_", "", .) %>%
    gsub("\\.csv$", "", .)
  
  # Initialize lists for both original IDs and species names
  landscapes <- list()
  srx_to_species <- list()
  
  for (i in seq_along(landscape_files)) {
    srx_id <- sample_ids[i]
    tryCatch({
      df <- read.csv(landscape_files[i])
      landscapes[[srx_id]] <- df
      
      # Map SRX ID to species name if available
      species_name <- get_species_name(srx_id, species_map)
      srx_to_species[[srx_id]] <- species_name
      
      cat(sprintf("  Loaded %s (%s): %d rows\n", srx_id, species_name, nrow(df)))
    }, error = function(e) {
      cat(sprintf("  Error loading %s: %s\n", landscape_files[i], e$message))
    })
  }
  
  cat(sprintf("Loaded %d landscape files\n", length(landscapes)))
  return(list(landscapes = landscapes, srx_to_species = srx_to_species))
}

# Function to extract Kimura bin as numeric
extract_kimura_bin <- function(kimura_bin) {
  # Handle different formats:
  # - "0-1" or "0–1" (en dash)
  # - "0" (single number)
  # - "0%" (percentage)
  
  # Remove any percentage signs
  kimura_clean <- gsub("%", "", kimura_bin)
  
  # Try to extract the first number from ranges (handles both normal and en dashes)
  first_num <- gsub("[-–].*$", "", kimura_clean)
  
  # Convert to numeric, with error handling
  result <- suppressWarnings(as.numeric(first_num))
  
  # Return 0 for NA values (this helps prevent errors in plotting)
  ifelse(is.na(result), 0, result)
}

# Function to generate TE order distribution visualization
plot_order_distribution <- function(data_list, dir = output_dir) {
  cat("Generating TE order distribution visualization...\n")
  
  # Extract the components
  landscapes <- data_list$landscapes
  srx_to_species <- data_list$srx_to_species
  
  # Create combined class data
  class_data <- data.frame()
  
  for (sample_id in names(landscapes)) {
    df <- landscapes[[sample_id]]
    
    # Find the order column
    order_col <- NULL
    possible_order_cols <- c("Order", "ORDER", "order")
    for (col in possible_order_cols) {
      if (col %in% colnames(df)) {
        order_col <- col
        break
      }
    }
    
    if (is.null(order_col)) {
      cat(sprintf("  Warning: No Order column found in %s\n", sample_id))
      next
    }
    
    # Find aligned bases column
    aligned_col <- NULL
    if ("aligned_bp" %in% colnames(df)) {
      aligned_col <- "aligned_bp"
    } else if ("aligned_bases" %in% colnames(df)) {
      aligned_col <- "aligned_bases"
    } else {
      cat(sprintf("  Warning: No aligned bases column found in %s\n", sample_id))
      next
    }
    
    # Group by order and sum
    summed <- df %>%
      group_by(!!sym(order_col)) %>%
      summarize(total_aligned = sum(!!sym(aligned_col)), .groups = "drop")
    
    # Add sample ID and species name
    summed$sample_id <- sample_id
    summed$species_name <- srx_to_species[[sample_id]]
    
    # Append to class data
    class_data <- rbind(class_data, summed)
  }
  
  if (nrow(class_data) > 0) {
    # Calculate percentages
    class_data <- class_data %>%
      group_by(sample_id, species_name) %>%
      mutate(percentage = total_aligned / sum(total_aligned) * 100) %>%
      ungroup()
    
    # Rename order column to a standard name for plotting
    names(class_data)[names(class_data) == order_col] <- "Order"
    
    # Create the stacked bar chart
    order_plot <- ggplot(class_data, aes(x = species_name, y = percentage, fill = Order)) +
      geom_bar(stat = "identity") +
      scale_fill_brewer(palette = "Set3") +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        legend.position = "right"
      ) +
      labs(
        title = "TE Order Distribution Across Desmognathus Species",
        x = "Species",
        y = "Percentage of Aligned Bases (%)"
      )
    
    # Save the plot
    ggsave(file.path(dir, "te_order_distribution.png"),
           order_plot, width = 12, height = 10, dpi = 300)
    
    cat("  Created TE order distribution plot\n")
  } else {
    cat("  No order data found in any sample, skipping order distribution plot\n")
  }
  
  cat("Visualization complete. Results saved to", dir, "\n")
}

# Main function
main <- function() {
  # Load species lookup table
  species_map <- load_species_lookup()
  
  # Load all landscape data
  data_list <- load_all_landscapes(species_map = species_map)
  
  # Generate order distribution plot only
  # Combined landscape and heatmap plots have been removed as they weren't useful
  plot_order_distribution(data_list)
  
  cat("Visualization completed successfully\n")
}

# Run the main function
main() 