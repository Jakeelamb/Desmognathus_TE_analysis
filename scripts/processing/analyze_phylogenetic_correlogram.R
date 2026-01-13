#!/usr/bin/env Rscript
# Description: Calculates and plots phylogenetic correlograms (Moran's I vs. distance)
#              for EACH TE superfamily proportion trait individually
#              using the phylosignal package and ggplot2.

# Load required libraries
suppressPackageStartupMessages({
  library(phylosignal)
  library(ape)
  library(phylobase) # Provides the phylo4d function/class
  library(readr)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(ggplot2)
  library(here)
  library(yaml)
  library(forcats)
  library(stringr)
  library(RColorBrewer)
  library(phytools)
})

# === Configuration ===
project_root <- here::here()
config_file <- file.path(project_root, "config/paths.yaml")

# Set default paths first
default_tree_path <- "results/phylogeny/processed_phylogeny.nwk"
default_traits_path <- "data/processed/diversity/superfamily_proportions.csv"
default_output_subdir <- "results/figures/phylo_signal"

# Initialize with defaults
tree_file <- file.path(project_root, default_tree_path)
traits_file <- file.path(project_root, default_traits_path)
output_dir <- file.path(project_root, default_output_subdir)
plot_output_file <- file.path(output_dir, "phylogenetic_correlogram_moran_per_trait.png")

if (file.exists(config_file)) {
  tryCatch({
    config <- yaml::read_yaml(config_file)
    cat("Config file found. Reading paths.\n")

    tree_path_from_config <- config$results$phylogeny
    if (!is.null(tree_path_from_config) && is.character(tree_path_from_config) && nzchar(tree_path_from_config)) {
      tree_file <- file.path(project_root, tree_path_from_config, "processed_phylogeny.nwk")
      cat("  Using tree path from config:", tree_path_from_config, "\n")
    } else {
      warning("  Config missing or invalid 'results$phylogeny'. Using default tree path.")
    }

    traits_path_from_config <- config$data$processed$diversity
    if (!is.null(traits_path_from_config) && is.character(traits_path_from_config) && nzchar(traits_path_from_config)) {
      traits_file <- file.path(project_root, traits_path_from_config, "superfamily_proportions.csv")
      cat("  Using traits path from config:", traits_path_from_config, "\n")
    } else {
      warning("  Config missing or invalid 'data$processed$diversity'. Using default traits path.")
    }

    output_dir_from_config <- config$results$figures$phylo_signal
    if (!is.null(output_dir_from_config) && is.character(output_dir_from_config) && nzchar(output_dir_from_config)) {
      output_dir <- file.path(project_root, output_dir_from_config)
      plot_output_file <- file.path(output_dir, "phylogenetic_correlogram_moran_per_trait.png")
      cat("  Using output directory from config:", output_dir_from_config, "\n")
    } else {
      warning("  Config missing or invalid 'results$figures$phylo_signal'. Using default output dir.")
    }

  }, error = function(e) {
    warning("Error reading config file '", config_file, "': ", conditionMessage(e), ". Using default paths.")
  })
} else {
  warning("config/paths.yaml not found. Using default paths.")
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# === Load Data ===
cat("Loading phylogenetic tree from:", tree_file, "\n")
if (!file.exists(tree_file)) stop("Phylogeny file not found: ", tree_file)
tree <- read.tree(tree_file)

cat("Loading trait data from:", traits_file, "\n")
if (!file.exists(traits_file)) stop("Traits file not found: ", traits_file)
traits_raw <- readr::read_csv(traits_file, show_col_types = FALSE)

if (ncol(traits_raw) > 0 && tolower(colnames(traits_raw)[1]) == "species") {
  traits <- traits_raw %>% tibble::column_to_rownames(var = colnames(traits_raw)[1])
  cat("Using column '", colnames(traits_raw)[1], "' as species identifiers.\n")
} else if (ncol(traits_raw) > 0) {
  warning("First column is not named 'Species'. Assuming it contains species identifiers and setting as row names.")
  traits <- traits_raw %>% tibble::column_to_rownames(var = colnames(traits_raw)[1])
} else {
  stop("Trait data file is empty or has no columns.")
}

traits <- traits %>% mutate(across(everything(), function(x) as.numeric(as.character(x))))

# === Data Matching ===
cat("Matching tree tips and trait data...\n")
species_to_keep <- intersect(tree$tip.label, rownames(traits))
if (length(species_to_keep) < length(tree$tip.label)) {
  cat("Pruning tree to", length(species_to_keep), "species found in trait data.\n")
  tree <- drop.tip(tree, setdiff(tree$tip.label, species_to_keep))
}
if (length(species_to_keep) < 4) {
  stop("Less than 4 matching species found. Cannot calculate correlogram reliably.")
}
traits <- traits[tree$tip.label, , drop = FALSE]
cat("Proceeding with", length(tree$tip.label), "matching species.\n")

traits_filtered <- traits[, sapply(traits, function(col) {
  valid_vals <- col[!is.na(col)]
  length(valid_vals) >= 4 && var(valid_vals, na.rm=TRUE) > 1e-8
}), drop = FALSE]

if (ncol(traits_filtered) == 0) stop("No traits remaining after filtering for variance and NAs.")
cat("Analyzing", ncol(traits_filtered), "traits with sufficient data and variance.\n")

# === Calculate Phylogenetic Correlograms (Per Trait) ===
cat("Calculating Moran's I phylogenetic correlograms for each trait...\n")

all_trait_correlograms <- list()

for (trait_name in colnames(traits_filtered)) {
  cat("  Processing trait:", trait_name, "\n")

  # Create single-trait data frame
  trait_data_single <- traits_filtered[, trait_name, drop = FALSE]

  # Create phylo4d object for this single trait
  tryCatch({
      p4d_single <- phylo4d(tree, trait_data_single)

      # Calculate correlogram (should compute Moran's I for single trait)
      correlogram_result <- phyloCorrelogram(p4d_single)

      # --- Process the result for this trait ---
      # Check if the result is a list and contains $res
      if (is.list(correlogram_result) && !is.null(correlogram_result$res)) {
          res_trait <- as.data.frame(correlogram_result$res)

          # Ensure res_trait is valid (not NULL and has rows and at least 2 columns)
          if (!is.null(res_trait) && nrow(res_trait) > 0 && ncol(res_trait) >= 2) {
              
              # --- Assume Column 2 is Moran's I --- 
              MoranI_values <- res_trait[[2]] # Select the second column
              # --- End Assumption --- 

              # Parse/Calculate distance midpoints
              distance_mid <- numeric(nrow(res_trait))
              parsed_ok <- FALSE
              if (!is.null(rownames(res_trait))) {
                  # Try parsing from rownames
                  matches <- regmatches(rownames(res_trait), regexec("\\[([0-9.]+), *([0-9.]+)\\]", rownames(res_trait)))
                  if(length(matches) == nrow(res_trait) && all(sapply(matches, length) == 3)) { # Check if parsing worked for all rows
                      lower_bounds <- suppressWarnings(as.numeric(sapply(matches, `[`, 2)))
                      upper_bounds <- suppressWarnings(as.numeric(sapply(matches, `[`, 3)))
                       if (!any(is.na(lower_bounds)) && !any(is.na(upper_bounds))) {
                           distance_mid <- (lower_bounds + upper_bounds) / 2
                           parsed_ok <- TRUE
                       }
                  }
              }
              if (!parsed_ok) {
                  # Fallback: Approximate midpoints
                  warning("Could not parse distance from rownames for trait: ", trait_name, ". Approximating.")
                  max_dist <- max(cophenetic(tree))
                  n_bins <- nrow(res_trait)
                  breaks_approx <- seq(0, max_dist, length.out = n_bins + 1)
                  distance_mid <- (breaks_approx[-1] + breaks_approx[-length(breaks_approx)]) / 2
              }

              # Store results
              all_trait_correlograms[[trait_name]] <- tibble(
                  distance_mid = distance_mid,
                  MoranI = MoranI_values, # Use the extracted values
                  trait = trait_name
              )

          } else {
               warning("phyloCorrelogram result $res was NULL, empty, or had fewer than 2 columns for trait: ", trait_name)
          }

      } else {
          warning("phyloCorrelogram did not return expected list structure with $res for trait: ", trait_name)
      }

  }, error = function(e) {
    warning("Error processing trait ", trait_name, ": ", conditionMessage(e))
  })
}

if (length(all_trait_correlograms) > 0) {
  plot_data <- bind_rows(all_trait_correlograms)
} else {
  stop("Failed to calculate correlograms for any trait.")
}

cat("Generating combined Moran's I correlogram plot...\n")

num_traits <- length(unique(plot_data$trait))
color_palette <- RColorBrewer::brewer.pal(min(num_traits, 9), "Set1")

plot_data <- plot_data %>%
  mutate(trait_display = str_replace_all(trait, "[._]", " ")) %>%
  mutate(trait_display = str_to_title(trait_display)) %>%
  mutate(trait_display = fct_reorder(trait_display, trait))

correlogram_plot <- ggplot(plot_data, aes(x = distance_mid, y = MoranI, color = trait_display, group = trait_display)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_line(linewidth = 1) +
  geom_point(size = 2.5, shape = 19) +
  scale_color_manual(values = rep(color_palette, length.out = num_traits)) +
  scale_y_continuous(limits = c(min(plot_data$MoranI, -0.2, na.rm = TRUE), max(plot_data$MoranI, 0.8, na.rm = TRUE))) +
  labs(
    title = "Phylogenetic Correlogram (Moran's I) per TE Superfamily",
    x = "Phylogenetic Distance Class Midpoint",
    y = "Moran's I",
    color = "TE Superfamily"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "right",
    legend.text = element_text(size = 9),
    legend.title = element_text(size = 10)
  ) +
  guides(color = guide_legend(ncol = 1))

ggsave(plot_output_file, plot = correlogram_plot, width = 10, height = 7, dpi = 300, bg = "white")
cat("Saved Moran's I correlogram plot to:", plot_output_file, "\n")
cat("Script finished.\n")
