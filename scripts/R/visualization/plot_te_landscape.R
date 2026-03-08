#!/usr/bin/env Rscript
#
# plot_te_landscape.R
# 
# Description: Simple script to visualize TE landscape data from a single sample
#
# Usage: Rscript plot_te_landscape.R [SRX_ID]
#
# Example: Rscript plot_te_landscape.R SRX19953421
#

# Check command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  cat("Error: No SRX ID provided\n")
  cat("Usage: Rscript plot_te_landscape.R [SRX_ID]\n")
  quit(status = 1)
}

srx_id <- args[1]
cat(paste("Generating plots for sample:", srx_id, "\n"))

script_args <- commandArgs(trailingOnly = FALSE)
script_path_arg <- grep("^--file=", script_args, value = TRUE)
if (length(script_path_arg) > 0) {
  script_path <- normalizePath(sub("^--file=", "", script_path_arg[1]), mustWork = FALSE)
  script_dir <- dirname(script_path)
} else {
  script_dir <- "scripts/R/visualization"
}

source(file.path(dirname(dirname(script_dir)), "R", "path_config_utils.R"))
prefer_active_conda_r_library()

project_root <- find_project_root(script_dir)
config <- load_project_config(project_root)
landscape_dir <- resolve_config_path(project_root, config$results$landscapes, "results/landscapes")
output_dir <- resolve_config_path(project_root, config$results$figures$landscape, "results/figures/landscape")

# Load required libraries after activating the conda R library path.
suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(viridis)
})

input_file <- file.path(landscape_dir, paste0("repeat_landscape_", srx_id, ".csv"))

# Check if input file exists
if (!file.exists(input_file)) {
  cat(paste("Error: Input file not found:", input_file, "\n"))
  cat("Please run parse_repeatmasker_landscape.py first\n")
  quit(status = 1)
}

# Create output directory if it doesn't exist
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Load the data
cat("Loading landscape data...\n")
landscape <- read_csv(input_file, show_col_types = FALSE)

# List available columns
cat("Available columns in data:\n")
print(colnames(landscape))

# Check if required columns exist
required_cols <- c("Species", "Kimura_bin", "Order", "Superfamily", "aligned_bp")
missing_cols <- setdiff(required_cols, colnames(landscape))

# Fix column names if needed (case insensitive matching)
if (length(missing_cols) > 0) {
  cat("Missing required columns:", paste(missing_cols, collapse=", "), "\n")
  cat("Attempting to match column names case-insensitively...\n")
  
  # Map existing columns to required columns
  for (req_col in missing_cols) {
    # Try to find a match ignoring case
    matches <- grep(tolower(req_col), tolower(colnames(landscape)), value = TRUE)
    if (length(matches) > 0) {
      cat(paste("Mapping", matches[1], "to", req_col, "\n"))
      colnames(landscape)[colnames(landscape) == matches[1]] <- req_col
    } else {
      cat(paste("Warning: Cannot find a match for", req_col, "\n"))
      # Add placeholder column if necessary
      if (req_col == "aligned_bp") {
        cat("Adding placeholder aligned_bp column\n")
        landscape$aligned_bp <- 1000 # Default value
      } else {
        landscape[[req_col]] <- "Unknown"
      }
    }
  }
}

# Add a bin_number column for proper x-axis ordering
landscape <- landscape %>%
  mutate(bin_number = as.numeric(gsub("–.*", "", Kimura_bin)))

# Define key TE groups for filtering
dna_orders <- c("TIR", "Helitron", "Maverick", "Crypton")
retro_orders <- c("LINE", "LTR", "SINE", "DIRS")

# Create overall landscape plot
cat("Creating overall landscape plot...\n")
overall <- landscape %>%
  group_by(Species, Kimura_bin, bin_number) %>%
  summarize(Aligned_bp = sum(aligned_bp), .groups = "drop") %>%
  arrange(Species, bin_number)

p1 <- ggplot(overall, aes(x = reorder(Kimura_bin, bin_number), y = Aligned_bp/1e6)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5)
  ) +
  labs(
    title = paste("TE Landscape -", srx_id),
    subtitle = paste("Species:", unique(landscape$Species)),
    x = "Kimura Distance (% - older →)",
    y = "Aligned Base Pairs (Mb)"
  )

# Save plot
ggsave(file.path(output_dir, paste0("overall_landscape_", srx_id, ".png")), 
       p1, width = 10, height = 6, dpi = 300)

# Create order-level plot if there are multiple orders
cat("Creating order-level landscape plot...\n")
order_counts <- landscape %>% 
  group_by(Order) %>% 
  summarize(Count = n())

if (nrow(order_counts) > 1) {
  order_data <- landscape %>%
    group_by(Species, Order, Kimura_bin, bin_number) %>%
    summarize(Aligned_bp = sum(aligned_bp), .groups = "drop") %>%
    arrange(Species, Order, bin_number)
  
  # Identify top orders by total aligned bp
  top_orders <- order_data %>%
    group_by(Order) %>%
    summarize(Total_bp = sum(Aligned_bp)) %>%
    arrange(desc(Total_bp)) %>%
    head(8) %>%
    pull(Order)
  
  # Filter to top orders for better visualization
  order_data_filtered <- order_data %>%
    filter(Order %in% top_orders)
  
  p2 <- ggplot(order_data_filtered, 
              aes(x = reorder(Kimura_bin, bin_number), 
                  y = Aligned_bp/1e6, 
                  fill = Order)) +
    geom_bar(stat = "identity") +
    scale_fill_viridis_d() +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(hjust = 0.5)
    ) +
    labs(
      title = paste("TE Orders Landscape -", srx_id),
      subtitle = paste("Species:", unique(landscape$Species)),
      x = "Kimura Distance (% - older →)",
      y = "Aligned Base Pairs (Mb)"
    )
  
  # Save plot
  ggsave(file.path(output_dir, paste0("order_landscape_", srx_id, ".png")), 
         p2, width = 10, height = 6, dpi = 300)
} else {
  cat("Skipping order-level plot as only one order is available\n")
}

# Create a smooth landscape plot with lines
cat("Creating smooth landscape plot...\n")
smooth_data <- landscape %>%
  group_by(Species, Kimura_bin, bin_number) %>%
  summarize(Aligned_bp = sum(aligned_bp), .groups = "drop") %>%
  arrange(Species, bin_number)

p3 <- ggplot(smooth_data, aes(x = bin_number, y = Aligned_bp/1e6)) +
  geom_line(linewidth = 1, color = "darkblue") +
  geom_point(color = "darkblue") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5)
  ) +
  labs(
    title = paste("Smoothed TE Age Landscape -", srx_id),
    subtitle = paste("Species:", unique(landscape$Species)),
    x = "Kimura Distance (%)",
    y = "Aligned Base Pairs (Mb)"
  )

# Save plot
ggsave(file.path(output_dir, paste0("smooth_landscape_", srx_id, ".png")), 
       p3, width = 10, height = 6, dpi = 300)

cat(paste("All plots saved to:", output_dir, "\n"))
cat("Completed successfully\n") 
