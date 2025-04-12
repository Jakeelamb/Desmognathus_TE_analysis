#!/usr/bin/env Rscript

# TE PCA Analysis Script
# This script performs PCA analysis on the TE superfamily proportions data
# and generates plots for visualizing the results.

# Load required libraries
library(tidyverse)
library(FactoMineR)
library(factoextra)
library(yaml)

# Function to parse command-line arguments
parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  
  # Initialize default values
  params <- list(
    input_file = NULL,
    config_file = "config/paths.yaml",
    output_dir = NULL,
    verbose = FALSE,
    min_species_presence = 3
  )
  
  i <- 1
  while (i <= length(args)) {
    if (args[i] == "--input" || args[i] == "-i") {
      params$input_file <- args[i + 1]
      i <- i + 1
    } else if (args[i] == "--output" || args[i] == "-o") {
      params$output_dir <- args[i + 1]
      i <- i + 1
    } else if (args[i] == "--config" || args[i] == "-c") {
      params$config_file <- args[i + 1]
      i <- i + 1
    } else if (args[i] == "--verbose" || args[i] == "-v") {
      params$verbose <- TRUE
    } else if (args[i] == "--min-species-presence") {
      params$min_species_presence <- as.numeric(args[i + 1])
      i <- i + 1
    } else {
      stop(paste("Unknown argument:", args[i]))
    }
    i <- i + 1
  }
  
  return(params)
}

# Function to load paths from configuration file
load_paths <- function(config_file) {
  # Read the YAML configuration file
  if (!file.exists(config_file)) {
    stop(paste("Config file not found:", config_file))
  }
  
  tryCatch({
    config <- yaml::read_yaml(config_file)
    
    # Define paths
    paths <- list(
      input = NULL,
      output_dir = NULL,
      tables_dir = NULL
    )
    
    # Get the data paths
    data_paths <- config$data
    results_paths <- config$results
    
    # Construct the paths
    paths$input <- file.path(data_paths$processed$diversity, "superfamily_proportions.csv")
    paths$output_dir <- results_paths$figures$pca
    paths$tables_dir <- results_paths$tables$pca
    
    return(paths)
  }, error = function(e) {
    stop(paste("Error loading config:", e$message))
  })
}

# Function to clean species names
clean_species_names <- function(names) {
  # Clean up the species names for nicer plotting
  cleaned <- names %>%
    gsub("\\.", " ", .) %>%
    gsub("_", " ", .) %>%
    # Make genus name italic
    gsub("^(\\w+) (\\w+)$", "italic('\\1') \\2", .) %>%
    # Make genus name italic for subspecies too
    gsub("^(\\w+) (\\w+) (\\w+)$", "italic('\\1') \\2 \\3", .)
  
  return(cleaned)
}

# Function to filter TEs based on species presence
filter_tes_by_species_presence <- function(data, min_species = 3) {
  # Get species columns (everything except the first column)
  species_cols <- colnames(data)[-1]
  
  # Count how many species have each TE
  te_species_counts <- rowSums(data[, species_cols] > 0, na.rm = TRUE)
  
  # Filter TEs that appear in at least the minimum number of species
  filtered_data <- data[te_species_counts >= min_species, ]
  
  # Return filtered data
  return(filtered_data)
}

# Main function to run the analysis
run_pca_analysis <- function() {
  # Parse command-line arguments
  params <- parse_args()
  
  # Set up verbose output
  verbose <- params$verbose
  
  # Load configuration
  paths <- load_paths(params$config_file)
  
  # Determine input and output paths
  input_file <- ifelse(is.null(params$input_file), paths$input, params$input_file)
  output_dir <- ifelse(is.null(params$output_dir), paths$output_dir, params$output_dir)
  tables_dir <- paths$tables_dir
  
  # Create output directories if they don't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  if (!dir.exists(tables_dir)) {
    dir.create(tables_dir, recursive = TRUE)
  }
  
  # Log paths
  if (verbose) {
    cat("Input file:", input_file, "\n")
    cat("Output directory:", output_dir, "\n")
    cat("Tables directory:", tables_dir, "\n")
    cat("Minimum species presence:", params$min_species_presence, "\n")
  }
  
  # Read the data
  if (!file.exists(input_file)) {
    stop(paste("Input file not found:", input_file))
  }
  
  data <- read.csv(input_file, row.names = 1)
  
  if (verbose) {
    cat("Data loaded with dimensions:", nrow(data), "x", ncol(data), "\n")
  }
  
  # Clean up the data
  # 1. Remove Unknown.TIR column (if present)
  unknown_col <- which(colnames(data) == "Unknown TIR")
  if (length(unknown_col) > 0) {
    if (verbose) {
      cat("Removing Unknown.TIR category\n")
    }
    data <- data[, -unknown_col]
  }
  
  # 2. Filter TE columns based on species presence
  original_cols <- ncol(data)
  # Count how many species have each TE (with non-zero values)
  te_presence_counts <- colSums(data > 0, na.rm = TRUE)
  # Keep only TEs present in at least min_species species
  data <- data[, te_presence_counts >= params$min_species_presence, drop = FALSE]
  
  if (verbose) {
    cat("Filtered TEs from", original_cols, "to", ncol(data), 
        "based on minimum species presence of", params$min_species_presence, "\n")
  }
  
  # Replace NA values with zeros
  data[is.na(data)] <- 0
  
  # No need to transpose - data is already in correct format for PCA
  # (rows = observations/species, columns = variables/TEs)
  
  # Perform PCA
  pca_result <- PCA(data, graph = FALSE)
  
  # Save the PCA results
  saveRDS(pca_result, file.path(tables_dir, "pca_results.rds"))
  
  # Extract and save eigenvalues
  eigenvalues <- get_eigenvalue(pca_result)
  write.csv(eigenvalues, file.path(tables_dir, "eigenvalues.csv"), row.names = TRUE)
  
  # Save variable contributions
  var_contrib <- get_pca_var(pca_result)$contrib
  write.csv(var_contrib, file.path(tables_dir, "variable_contributions.csv"), row.names = TRUE)
  
  # Prepare species names for plotting
  species_names <- rownames(data)
  species_labels <- clean_species_names(species_names)
  
  # Extract genus names for coloring
  genera <- gsub("D\\.", "", species_names)  # Remove "D." prefix
  genera <- gsub("\\..*", "", genera)       # Keep only genus portion
  genus_colors <- setNames(scales::hue_pal()(length(unique(genera))), unique(genera))
  
  # Generate plots
  # 1. Scree plot (variance explained)
  pdf(file.path(output_dir, "scree_plot.pdf"), width = 10, height = 8)
  fviz_eig(pca_result, addlabels = TRUE, ylim = c(0, 50))
  dev.off()
  
  # 2. Biplot of PC1 and PC2
  pdf(file.path(output_dir, "biplot_pc1_pc2.pdf"), width = 12, height = 10)
  p <- fviz_pca_biplot(pca_result,
                   label = "var",
                   col.ind = genera,
                   palette = genus_colors,
                   repel = TRUE,
                   title = "PCA - Biplot (PC1 & PC2)")
  print(p)
  dev.off()
  
  # 3. Individuals plot (species)
  pdf(file.path(output_dir, "individuals_plot.pdf"), width = 12, height = 10)
  p <- fviz_pca_ind(pca_result,
                col.ind = genera,
                palette = genus_colors,
                repel = TRUE,
                title = "PCA - Species")
  print(p)
  dev.off()
  
  # 4. Variables plot (TE superfamilies)
  pdf(file.path(output_dir, "variables_plot.pdf"), width = 14, height = 12)
  p <- fviz_pca_var(pca_result,
                col.var = "contrib",
                gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                repel = TRUE,
                title = "PCA - TE Superfamilies")
  print(p)
  dev.off()
  
  # 5. Contribution of variables to PC1
  pdf(file.path(output_dir, "contrib_pc1.pdf"), width = 12, height = 16)
  p <- fviz_contrib(pca_result, choice = "var", axes = 1, top = 30)
  print(p)
  dev.off()
  
  # 6. Contribution of variables to PC2
  pdf(file.path(output_dir, "contrib_pc2.pdf"), width = 12, height = 16)
  p <- fviz_contrib(pca_result, choice = "var", axes = 2, top = 30)
  print(p)
  dev.off()
  
  # 7. PCA Individual plot with PC3 and PC4
  pdf(file.path(output_dir, "individuals_plot_pc3_pc4.pdf"), width = 12, height = 10)
  p <- fviz_pca_ind(pca_result,
                axes = c(3, 4),
                col.ind = genera,
                palette = genus_colors,
                repel = TRUE,
                title = "PCA - Species (PC3 & PC4)")
  print(p)
  dev.off()
  
  # Log success
  if (verbose) {
    cat("PCA analysis completed successfully. Results saved to:", output_dir, "\n")
  }
}

# Run the analysis
tryCatch({
  run_pca_analysis()
  quit(status = 0)
}, error = function(e) {
  cat("Error:", e$message, "\n")
  quit(status = 1)
}) 