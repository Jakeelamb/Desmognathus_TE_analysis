#!/usr/bin/env Rscript
# Description: Calculates and plots phylogenetic correlograms (Moran's I vs. distance)
#              for each TE superfamily trait using distance-binned Moran's I
#              computed directly from the tree cophenetic matrix.

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

# Load required libraries
suppressPackageStartupMessages({
  library(ape)
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(stringr)
  library(RColorBrewer)
})

# === Configuration ===
project_root <- find_project_root(script_dir)
config <- load_project_config(project_root)

phylo_dir <- resolve_config_path(project_root, config$input_data$phylogeny %||% config$data$phylogeny, "input_data/phylogeny")
results_data_dir <- resolve_config_path(project_root, config$results$data, "results/data")
output_dir <- resolve_config_path(project_root, config$results$figures$phylo_signal, "results/figures/phylo_signal")
output_table_dir <- resolve_config_path(project_root, config$results$tables$phylo_signal, "results/tables/phylogenetic_signal")

tree_file <- file.path(phylo_dir, "desmo900dated_test.tre")
if (!file.exists(tree_file)) {
  tree_file <- file.path(project_root, "results", "phylogeny", "processed_phylogeny.nwk")
}

traits_file <- file.path(results_data_dir, "dnaPipeTE_superfamily_breakdown.csv")

plot_output_file <- file.path(output_dir, "phylogenetic_correlogram_moran_per_trait.png")
table_output_file <- file.path(output_table_dir, "phylogenetic_correlogram_moran_per_trait.csv")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(output_table_dir, recursive = TRUE, showWarnings = FALSE)

# === Load Data ===
cat("Loading phylogenetic tree from:", tree_file, "\n")
if (!file.exists(tree_file)) stop("Phylogeny file not found: ", tree_file)
tree <- read.tree(tree_file)

cat("Loading trait data from:", traits_file, "\n")
if (!file.exists(traits_file)) stop("Traits file not found: ", traits_file)
traits_raw <- readr::read_csv(traits_file, show_col_types = FALSE)

if (ncol(traits_raw) > 0 && !nzchar(colnames(traits_raw)[1])) {
  colnames(traits_raw)[1] <- "species"
}

if (ncol(traits_raw) > 0 && tolower(colnames(traits_raw)[1]) == "species") {
  traits <- as.data.frame(traits_raw, check.names = FALSE) %>% tibble::column_to_rownames(var = colnames(traits_raw)[1])
  cat("Using column '", colnames(traits_raw)[1], "' as species identifiers.\n")
} else if (ncol(traits_raw) > 0) {
  warning("First column is not named 'Species'. Assuming it contains species identifiers and setting as row names.")
  traits <- as.data.frame(traits_raw, check.names = FALSE) %>% tibble::column_to_rownames(var = colnames(traits_raw)[1])
} else {
  stop("Trait data file is empty or has no columns.")
}

rownames(traits) <- gsub("^D\\.", "", rownames(traits))
trait_species <- rownames(traits)
traits <- as.data.frame(lapply(traits, function(x) as.numeric(as.character(x))), check.names = FALSE)
rownames(traits) <- trait_species

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

distance_matrix <- cophenetic(tree)

build_distance_bins <- function(distance_matrix, n_bins = 8) {
  pairwise_distances <- distance_matrix[lower.tri(distance_matrix)]
  pairwise_distances <- pairwise_distances[is.finite(pairwise_distances) & pairwise_distances > 0]

  if (length(pairwise_distances) < n_bins) {
    stop("Not enough pairwise distances to build correlogram bins.")
  }

  breaks <- unique(as.numeric(
    stats::quantile(pairwise_distances, probs = seq(0, 1, length.out = n_bins + 1), type = 8, na.rm = TRUE)
  ))

  if (length(breaks) < 3) {
    rng <- range(pairwise_distances, na.rm = TRUE)
    breaks <- seq(rng[1], rng[2], length.out = min(5, length(unique(pairwise_distances))) + 1)
  }

  tibble(
    bin_id = seq_len(length(breaks) - 1),
    lower = breaks[-length(breaks)],
    upper = breaks[-1]
  ) %>%
    filter(upper > lower) %>%
    mutate(distance_mid = (lower + upper) / 2)
}

compute_binned_morans_i <- function(trait_vector, distance_matrix, distance_bins, min_pairs = 5) {
  valid <- !is.na(trait_vector)
  values <- as.numeric(trait_vector[valid])

  if (length(values) < 4 || stats::var(values) <= 1e-8) {
    return(NULL)
  }

  dist_sub <- distance_matrix[valid, valid, drop = FALSE]
  z <- values - mean(values)
  denom <- sum(z ^ 2)

  if (!is.finite(denom) || denom <= 0) {
    return(NULL)
  }

  z_outer <- tcrossprod(z)
  lower_mask <- lower.tri(dist_sub)

  binned_results <- lapply(seq_len(nrow(distance_bins)), function(i) {
    lower <- distance_bins$lower[[i]]
    upper <- distance_bins$upper[[i]]

    pair_mask_lower <- lower_mask & dist_sub > lower & dist_sub <= upper
    pair_count <- sum(pair_mask_lower)
    if (pair_count < min_pairs) {
      return(NULL)
    }

    weight_matrix <- (dist_sub > lower & dist_sub <= upper) * 1
    diag(weight_matrix) <- 0
    s0 <- sum(weight_matrix)
    if (s0 <= 0) {
      return(NULL)
    }

    morans_i <- (length(values) / s0) * sum(weight_matrix * z_outer) / denom

    tibble(
      bin_id = distance_bins$bin_id[[i]],
      lower = lower,
      upper = upper,
      distance_mid = distance_bins$distance_mid[[i]],
      pair_count = pair_count,
      MoranI = morans_i
    )
  })

  bind_rows(binned_results)
}

distance_bins <- build_distance_bins(distance_matrix, n_bins = 8)
cat("Using", nrow(distance_bins), "quantile-based phylogenetic distance bins.\n")

all_trait_correlograms <- list()

for (trait_name in colnames(traits_filtered)) {
  cat("  Processing trait:", trait_name, "\n")

  tryCatch({
      trait_result <- compute_binned_morans_i(traits_filtered[[trait_name]], distance_matrix, distance_bins)
      if (!is.null(trait_result) && nrow(trait_result) > 0) {
        all_trait_correlograms[[trait_name]] <- trait_result %>%
          mutate(trait = trait_name, .before = 1)
      } else {
        warning("No valid correlogram bins retained for trait: ", trait_name)
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
if (num_traits > 9) {
  color_palette <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(8, "Dark2"))(num_traits)
}

plot_data <- plot_data %>%
  mutate(trait_display = str_replace_all(trait, "[._]", " ")) %>%
  mutate(trait_display = str_to_title(trait_display))

trait_levels <- plot_data %>%
  distinct(trait, trait_display) %>%
  arrange(trait) %>%
  pull(trait_display)

plot_data <- plot_data %>%
  mutate(trait_display = factor(trait_display, levels = trait_levels))

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

readr::write_csv(plot_data, table_output_file)
ggsave(plot_output_file, plot = correlogram_plot, width = 10, height = 7, dpi = 300, bg = "white")
cat("Saved correlogram data to:", table_output_file, "\n")
cat("Saved Moran's I correlogram plot to:", plot_output_file, "\n")
cat("Script finished.\n")
