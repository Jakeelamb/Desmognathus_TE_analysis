#!/usr/bin/env Rscript
#' PERMANOVA Analysis for TE Composition
#'
#' Performs permutational multivariate analysis of variance (PERMANOVA) to test
#' for compositional differences between phylogenetic clades. Also includes
#' beta dispersion tests and pairwise distance analysis.
#'
#' @usage Rscript scripts/processing/permanova_analysis.R

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

# --- Load Dependencies ---
required_packages <- c(
  "vegan",        # PERMANOVA (adonis2), betadisper
  "ape",          # Phylogenetic tree handling
  "dplyr",        # Data manipulation
  "tidyr",        # Data reshaping
  "readr",        # CSV reading
  "tibble",       # Modern data frames
  "ggplot2",      # Plotting
  "pheatmap",     # Heatmaps
  "RColorBrewer", # Color palettes
  "yaml"          # Config loading
)

missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop(
    "Missing required R packages: ",
    paste(missing_packages, collapse = ", "),
    ". Update the Dusky environment instead of installing packages at runtime."
  )
}
invisible(lapply(required_packages, library, character.only = TRUE))

RNG_SEED <- 20260310
set.seed(RNG_SEED)

# --- Configuration ---
project_root <- find_project_root(script_dir)
config <- load_project_config(project_root)

# Paths
data_dir <- resolve_config_path(project_root, config$results$data, "results/data")
figures_dir <- resolve_config_path(project_root, config$results$figures, "results/figures")
phylo_dir <- resolve_config_path(project_root, config$input_data$phylogeny %||% config$data$phylogeny, "input_data/phylogeny")
output_dir <- file.path(figures_dir, "permanova")
output_data_dir <- file.path(data_dir, "permanova")

# Create output directories
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(output_data_dir, showWarnings = FALSE, recursive = TRUE)

message("Project root: ", project_root)
message("Output directory: ", output_dir)

# --- Load Data ---
message("\n=== Loading Data ===")

# Load phylogenetic tree
tree_file <- file.path(phylo_dir, "desmo900dated_test.tre")
if (!file.exists(tree_file)) {
  tree_file <- file.path(data_dir, "desmo900dated_test_cleaned_phylo.tre")
}
if (!file.exists(tree_file)) {
  stop("Phylogenetic tree not found")
}
tree <- read.tree(tree_file)
message("Loaded tree with ", length(tree$tip.label), " tips")

# Load TE composition data
superfamily_file <- file.path(data_dir, "dnaPipeTE_superfamily_breakdown.csv")
order_file <- file.path(data_dir, "dnaPipeTE_order_breakdown.csv")
class_file <- file.path(data_dir, "dnaPipeTE_class_breakdown.csv")

superfamily_df <- read_csv(superfamily_file, show_col_types = FALSE)
order_df <- read_csv(order_file, show_col_types = FALSE)
class_df <- read_csv(class_file, show_col_types = FALSE)

message("Loaded TE composition data")

# --- Define Phylogenetic Clades ---
message("\n=== Defining Phylogenetic Clades ===")

define_curated_clades <- function(tree, min_named_clade_size = 4) {
  #' Define biologically interpretable clades and collapse undersized groups.

  tips <- tree$tip.label
  curated <- list(
    quadramaculatus_group = c("marmoratus", "intermedius", "quadramaculatus"),
    mountain_clade = c("orestes", "carolinensis", "anicetus", "ocoee",
                       "imitator", "santeetlah", "conanti"),
    small_clade = c("wrighti", "aeneus", "aureatus", "gvnigeusgwotli"),
    coastal_clade = c("auriculatus", "brimleyorum", "fuscus",
                      "monticola", "apalachicolae"),
    southern_appalachian = c("organi", "valentinei", "welteri",
                             "folkertsi", "pascagoula")
  )

  species_clade <- tibble::tibble(
    species = tips,
    clade = "other",
    clade_source = "other"
  )

  for (clade_name in names(curated)) {
    matched <- intersect(curated[[clade_name]], tips)
    if (length(matched) == 0) {
      next
    }
    species_clade <- species_clade %>%
      mutate(
        clade = ifelse(species %in% matched, clade_name, clade),
        clade_source = ifelse(species %in% matched, "curated_manual", clade_source)
      )
  }

  clade_sizes <- table(species_clade$clade)
  small_named <- names(clade_sizes)[names(clade_sizes) != "other" & clade_sizes < min_named_clade_size]
  if (length(small_named) > 0) {
    species_clade <- species_clade %>%
      mutate(
        clade_source = ifelse(clade %in% small_named, "collapsed_small_named", clade_source),
        clade = ifelse(clade %in% small_named, "other", clade)
      )
  }

  species_clade %>%
    mutate(pairwise_clade = ifelse(clade == "other", NA_character_, clade)) %>%
    arrange(clade, species)
}

species_clades <- define_curated_clades(tree, min_named_clade_size = 4)
message("Defined curated clades for PERMANOVA:")
print(table(species_clades$clade))

# --- Prepare Data for PERMANOVA ---
message("\n=== Preparing Data ===")

prepare_permanova_data <- function(te_df, species_clades, te_level = "order") {
  #' Prepare TE composition matrix and grouping factor for PERMANOVA

  # Get species column
  species_col <- names(te_df)[1]

  # Clean species names
  te_df[[species_col]] <- gsub("^D\\.", "", te_df[[species_col]])

  # Merge with clade assignments
  te_df <- te_df %>%
    rename(species = !!species_col) %>%
    inner_join(species_clades, by = "species")

  if (nrow(te_df) < 5) {
    stop("Too few species with clade assignments")
  }

  # Extract numeric columns for composition matrix
  numeric_cols <- names(te_df)[sapply(te_df, is.numeric)]
  composition_matrix <- as.matrix(te_df[, numeric_cols])
  rownames(composition_matrix) <- te_df$species

  # Remove columns with zero variance
  col_var <- apply(composition_matrix, 2, var, na.rm = TRUE)
  composition_matrix <- composition_matrix[, col_var > 0, drop = FALSE]

  # Handle any remaining NAs
  composition_matrix[is.na(composition_matrix)] <- 0

  list(
    matrix = composition_matrix,
    groups = te_df$clade,
    pairwise_groups = te_df$pairwise_clade,
    species = te_df$species,
    metadata = te_df
  )
}

order_data <- prepare_permanova_data(order_df, species_clades, "order")
superfamily_data <- prepare_permanova_data(superfamily_df, species_clades, "superfamily")

message("Order matrix: ", nrow(order_data$matrix), " species x ", ncol(order_data$matrix), " TEs")
message("Superfamily matrix: ", nrow(superfamily_data$matrix), " species x ", ncol(superfamily_data$matrix), " TEs")

# --- Calculate Distance Matrices ---
message("\n=== Calculating Distance Matrices ===")

# Bray-Curtis dissimilarity (most common for compositional data)
order_dist_bray <- vegdist(order_data$matrix, method = "bray")
superfamily_dist_bray <- vegdist(superfamily_data$matrix, method = "bray")

# Euclidean distance (for CLR-transformed data comparison)
# CLR transform for Euclidean-based metrics
clr_transform <- function(mat) {
  # Add small constant to avoid log(0)
  mat[mat == 0] <- min(mat[mat > 0]) / 10
  log_mat <- log(mat)
  row_means <- rowMeans(log_mat)
  sweep(log_mat, 1, row_means, "-")
}

order_clr <- clr_transform(order_data$matrix)
superfamily_clr <- clr_transform(superfamily_data$matrix)

order_dist_eucl <- dist(order_clr)
superfamily_dist_eucl <- dist(superfamily_clr)

# --- PERMANOVA Analysis ---
message("\n=== Running PERMANOVA ===")

run_permanova <- function(dist_matrix, groups, nperm = 999) {
  #' Run PERMANOVA using adonis2

  # Create data frame for formula
  perm_df <- data.frame(group = groups)

  # Run PERMANOVA
  result <- adonis2(
    dist_matrix ~ group,
    data = perm_df,
    permutations = nperm
  )

  return(result)
}

# Order-level PERMANOVA
message("\n--- Order-level PERMANOVA (Bray-Curtis) ---")
order_permanova_bray <- run_permanova(order_dist_bray, order_data$groups)
print(order_permanova_bray)

message("\n--- Order-level PERMANOVA (Euclidean on CLR) ---")
order_permanova_eucl <- run_permanova(order_dist_eucl, order_data$groups)
print(order_permanova_eucl)

# Superfamily-level PERMANOVA
message("\n--- Superfamily-level PERMANOVA (Bray-Curtis) ---")
superfamily_permanova_bray <- run_permanova(superfamily_dist_bray, superfamily_data$groups)
print(superfamily_permanova_bray)

message("\n--- Superfamily-level PERMANOVA (Euclidean on CLR) ---")
superfamily_permanova_eucl <- run_permanova(superfamily_dist_eucl, superfamily_data$groups)
print(superfamily_permanova_eucl)

# --- Beta Dispersion Test ---
message("\n=== Beta Dispersion Analysis ===")

test_beta_dispersion <- function(dist_matrix, groups, method_name = "") {
  #' Test for homogeneity of multivariate dispersions

  # Calculate beta dispersion
  betadisp_result <- betadisper(dist_matrix, groups)

  # Permutation test
  permutest_result <- permutest(betadisp_result, permutations = 999)

  message("\nBeta dispersion test (", method_name, "):")
  message("  F-value: ", round(permutest_result$tab$F[1], 3))
  message("  p-value: ", round(permutest_result$tab$`Pr(>F)`[1], 4))

  if (permutest_result$tab$`Pr(>F)`[1] < 0.05) {
    message("  WARNING: Significant heterogeneity in dispersions - interpret PERMANOVA with caution")
  }

  list(
    betadisper = betadisp_result,
    permutest = permutest_result
  )
}

order_betadisp <- test_beta_dispersion(order_dist_bray, order_data$groups, "Order Bray-Curtis")
superfamily_betadisp <- test_beta_dispersion(superfamily_dist_bray, superfamily_data$groups, "Superfamily Bray-Curtis")

# --- Pairwise Comparisons ---
run_pairwise_permanova <- function(dist_matrix, groups, min_group_n = 4, nperm = 999, exclude_groups = "other") {
  group_vec <- as.character(groups)
  group_vec[group_vec %in% exclude_groups] <- NA_character_
  group_counts <- table(group_vec)
  valid_groups <- sort(unique(stats::na.omit(group_vec)))
  group_pairs <- combn(valid_groups, 2, simplify = FALSE)

  pairwise_results <- lapply(group_pairs, function(pair) {
    pair_counts <- group_counts[pair]
    if (any(pair_counts < min_group_n)) {
      return(NULL)
    }

    keep <- group_vec %in% pair
    pair_groups <- factor(group_vec[keep], levels = pair)
    pair_dist <- as.dist(as.matrix(dist_matrix)[keep, keep, drop = FALSE])
    pair_df <- data.frame(group = pair_groups)
    fit <- adonis2(pair_dist ~ group, data = pair_df, permutations = nperm)

    tibble::tibble(
      group_1 = pair[[1]],
      group_2 = pair[[2]],
      n_group_1 = unname(pair_counts[[1]]),
      n_group_2 = unname(pair_counts[[2]]),
      r_squared = fit$R2[[1]],
      f_value = fit$F[[1]],
      p_value = fit$`Pr(>F)`[[1]]
    )
  })

  results <- dplyr::bind_rows(pairwise_results)
  if (nrow(results) == 0) {
    return(results)
  }

  results %>%
    mutate(
      distance_metric = "bray",
      clade_scope = "named_curated_clades_only",
      p_adjusted = p.adjust(p_value, method = "BH"),
      significant = p_adjusted < 0.05
    ) %>%
    arrange(p_adjusted, p_value)
}

message("\n=== Pairwise PERMANOVA Comparisons ===")

message("\n--- Order-level Pairwise Comparisons ---")
order_pairwise <- run_pairwise_permanova(order_dist_bray, order_data$pairwise_groups)
if (nrow(order_pairwise) == 0) {
  message("No order-level pairwise comparisons met the minimum group size requirement.")
} else {
  print(order_pairwise)
}

message("\n--- Superfamily-level Pairwise Comparisons ---")
superfamily_pairwise <- run_pairwise_permanova(superfamily_dist_bray, superfamily_data$pairwise_groups)
if (nrow(superfamily_pairwise) == 0) {
  message("No superfamily-level pairwise comparisons met the minimum group size requirement.")
} else {
  print(superfamily_pairwise)
}

# --- Distance Matrix Visualization ---
message("\n=== Generating Distance Heatmaps ===")

plot_distance_heatmap <- function(dist_matrix, groups, species_names, title, output_file) {
  #' Create heatmap of pairwise distances

  # Convert to matrix
  dist_mat <- as.matrix(dist_matrix)
  rownames(dist_mat) <- species_names
  colnames(dist_mat) <- species_names

  # Create annotation
  annotation_row <- data.frame(
    Clade = groups,
    row.names = species_names
  )

  # Order by clade
  clade_order <- order(groups)
  dist_mat <- dist_mat[clade_order, clade_order]
  annotation_row <- annotation_row[clade_order, , drop = FALSE]

  # Generate heatmap
  png(output_file, width = 1200, height = 1000, res = 150)
  pheatmap(
    dist_mat,
    annotation_row = annotation_row,
    annotation_col = annotation_row,
    main = title,
    fontsize_row = 7,
    fontsize_col = 7,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    color = colorRampPalette(c("white", "steelblue", "darkblue"))(100)
  )
  dev.off()

  message("Saved: ", output_file)
}

# Order distance heatmap
plot_distance_heatmap(
  order_dist_bray,
  order_data$groups,
  order_data$species,
  "Bray-Curtis Dissimilarity (TE Order Composition)",
  file.path(output_dir, "distance_heatmap_order_bray.png")
)

# Superfamily distance heatmap
plot_distance_heatmap(
  superfamily_dist_bray,
  superfamily_data$groups,
  superfamily_data$species,
  "Bray-Curtis Dissimilarity (TE Superfamily Composition)",
  file.path(output_dir, "distance_heatmap_superfamily_bray.png")
)

# --- PCoA Visualization ---
message("\n=== PCoA Visualization ===")

plot_pcoa <- function(dist_matrix, groups, species_names, title, output_file) {
  #' Create PCoA plot colored by clade

  # Run PCoA
  pcoa_result <- cmdscale(dist_matrix, k = 2, eig = TRUE)

  # Calculate variance explained
  var_explained <- round(100 * pcoa_result$eig[1:2] / sum(pcoa_result$eig[pcoa_result$eig > 0]), 1)

  # Create plot data
  plot_df <- data.frame(
    PCo1 = pcoa_result$points[, 1],
    PCo2 = pcoa_result$points[, 2],
    Species = species_names,
    Clade = groups
  )

  ellipse_df <- plot_df %>%
    dplyr::add_count(Clade, name = "clade_n") %>%
    dplyr::filter(clade_n >= 3)

  # Plot
  p <- ggplot(plot_df, aes(x = PCo1, y = PCo2, color = Clade, label = Species)) +
    geom_point(size = 3, alpha = 0.8) +
    ggrepel::geom_text_repel(size = 2.5, max.overlaps = 15) +
    stat_ellipse(
      data = ellipse_df,
      mapping = aes(x = PCo1, y = PCo2, color = Clade, group = Clade),
      inherit.aes = FALSE,
      level = 0.95,
      linetype = "dashed",
      alpha = 0.5
    ) +
    scale_color_brewer(palette = "Set1") +
    labs(
      title = title,
      x = paste0("PCo1 (", var_explained[1], "%)"),
      y = paste0("PCo2 (", var_explained[2], "%)")
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      legend.position = "right"
    ) +
    coord_fixed()

  ggsave(output_file, p, width = 12, height = 10, dpi = 300)
  message("Saved: ", output_file)
}

plot_pcoa(
  order_dist_bray,
  order_data$groups,
  order_data$species,
  "PCoA of TE Order Composition (Bray-Curtis)",
  file.path(output_dir, "pcoa_order_bray.png")
)

plot_pcoa(
  superfamily_dist_bray,
  superfamily_data$groups,
  superfamily_data$species,
  "PCoA of TE Superfamily Composition (Bray-Curtis)",
  file.path(output_dir, "pcoa_superfamily_bray.png")
)

# --- Beta Dispersion Plots ---
message("\n=== Beta Dispersion Plots ===")

png(file.path(output_dir, "betadispersion_order.png"), width = 1000, height = 800, res = 150)
plot(order_betadisp$betadisper, main = "Beta Dispersion: TE Order Composition",
     hull = FALSE, ellipse = TRUE)
dev.off()

png(file.path(output_dir, "betadispersion_superfamily.png"), width = 1000, height = 800, res = 150)
plot(superfamily_betadisp$betadisper, main = "Beta Dispersion: TE Superfamily Composition",
     hull = FALSE, ellipse = TRUE)
dev.off()

# --- Save Results ---
message("\n=== Saving Results ===")

# Create summary data frame
permanova_summary <- data.frame(
  analysis = c("Order_BrayCurtis", "Order_Euclidean_CLR",
               "Superfamily_BrayCurtis", "Superfamily_Euclidean_CLR"),
  distance_space = c("bray", "clr_euclidean", "bray", "clr_euclidean"),
  clade_scheme = "curated_manual_with_small_named_collapsed_to_other",
  rng_seed = RNG_SEED,
  permutations = 999,
  R2 = c(
    order_permanova_bray$R2[1],
    order_permanova_eucl$R2[1],
    superfamily_permanova_bray$R2[1],
    superfamily_permanova_eucl$R2[1]
  ),
  F_value = c(
    order_permanova_bray$F[1],
    order_permanova_eucl$F[1],
    superfamily_permanova_bray$F[1],
    superfamily_permanova_eucl$F[1]
  ),
  p_value = c(
    order_permanova_bray$`Pr(>F)`[1],
    order_permanova_eucl$`Pr(>F)`[1],
    superfamily_permanova_bray$`Pr(>F)`[1],
    superfamily_permanova_eucl$`Pr(>F)`[1]
  ),
  beta_disp_p = c(
    order_betadisp$permutest$tab$`Pr(>F)`[1],
    NA,  # Only have beta dispersion for Bray-Curtis
    superfamily_betadisp$permutest$tab$`Pr(>F)`[1],
    NA  # Only have beta dispersion for Bray-Curtis
  ),
  significant = c(
    order_permanova_bray$`Pr(>F)`[1] < 0.05,
    order_permanova_eucl$`Pr(>F)`[1] < 0.05,
    superfamily_permanova_bray$`Pr(>F)`[1] < 0.05,
    superfamily_permanova_eucl$`Pr(>F)`[1] < 0.05
  )
)

write_csv(permanova_summary, file.path(output_data_dir, "permanova_summary.csv"))
message("Saved PERMANOVA summary")

# Save distance matrices
write.csv(as.matrix(order_dist_bray),
          file.path(output_data_dir, "distance_matrix_order_bray.csv"))
write.csv(as.matrix(superfamily_dist_bray),
          file.path(output_data_dir, "distance_matrix_superfamily_bray.csv"))
message("Saved distance matrices")

# Save clade assignments
write_csv(species_clades, file.path(output_data_dir, "species_clade_assignments.csv"))
message("Saved clade assignments")

if (exists("order_pairwise") && nrow(order_pairwise) > 0) {
  order_pairwise <- order_pairwise %>% mutate(rng_seed = RNG_SEED, permutations = 999)
  write_csv(order_pairwise, file.path(output_data_dir, "permanova_order_pairwise.csv"))
  message("Saved order pairwise PERMANOVA results")
}

if (exists("superfamily_pairwise") && nrow(superfamily_pairwise) > 0) {
  superfamily_pairwise <- superfamily_pairwise %>% mutate(rng_seed = RNG_SEED, permutations = 999)
  write_csv(superfamily_pairwise, file.path(output_data_dir, "permanova_superfamily_pairwise.csv"))
  message("Saved superfamily pairwise PERMANOVA results")
}

# --- Summary ---
message("\n=== PERMANOVA Analysis Summary ===")
message("\nOrder-level PERMANOVA (Bray-Curtis):")
message("  R2 = ", round(order_permanova_bray$R2[1], 4))
message("  F = ", round(order_permanova_bray$F[1], 2))
message("  p = ", format(order_permanova_bray$`Pr(>F)`[1], scientific = FALSE))

message("\nSuperfamily-level PERMANOVA (Bray-Curtis):")
message("  R2 = ", round(superfamily_permanova_bray$R2[1], 4))
message("  F = ", round(superfamily_permanova_bray$F[1], 2))
message("  p = ", format(superfamily_permanova_bray$`Pr(>F)`[1], scientific = FALSE))

message("\nInterpretation:")
if (order_permanova_bray$`Pr(>F)`[1] < 0.05 || order_permanova_eucl$`Pr(>F)`[1] < 0.05) {
  message("  - Order-level signal is present in at least one distance space; inspect Bray and CLR results separately")
} else {
  message("  - No robust order-level clade difference was detected across the tested distance spaces")
}

if (superfamily_permanova_bray$`Pr(>F)`[1] < 0.05 && superfamily_permanova_eucl$`Pr(>F)`[1] < 0.05) {
  message("  - Superfamily-level TE composition differs between clades in both Bray and CLR spaces")
} else if (superfamily_permanova_bray$`Pr(>F)`[1] < 0.05 || superfamily_permanova_eucl$`Pr(>F)`[1] < 0.05) {
  message("  - Superfamily-level clade structure is distance-sensitive: significant in one space but not the other")
} else {
  message("  - No superfamily-level clade difference was detected in either distance space")
}

message("\n=== PERMANOVA Analysis Complete ===")
message("Results saved to: ", output_data_dir)
message("Figures saved to: ", output_dir)
