#!/usr/bin/env Rscript
#' Trait Evolution Modeling for Desmognathus TE Composition
#'
#' Analyzes phylogenetic signal and evolutionary patterns in TE traits.
#' Uses available packages (ape, caper) without requiring geiger/phytools.
#'
#' @usage Rscript scripts/processing/trait_evolution.R

# --- Load Dependencies ---
required_packages <- c(
  "ape",          # Phylogenetic tree handling, pic, ace
  "caper",        # Comparative analysis
  "nlme",         # GLS for phylogenetic models
  "dplyr",        # Data manipulation
  "tidyr",        # Data reshaping
  "readr",        # CSV reading
  "tibble",       # Modern data frames
  "ggplot2",      # Plotting
  "yaml"          # Config loading
)

# Load packages
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org", quiet = TRUE)
  }
  library(pkg, character.only = TRUE)
}

# --- Configuration ---
find_project_root <- function() {
  current <- getwd()
  for (i in 1:10) {
    if (file.exists(file.path(current, "paths.yaml"))) {
      return(current)
    }
    parent <- dirname(current)
    if (parent == current) break
    current <- parent
  }
  candidates <- c("~/Projects/Desmognathus_TE", getwd())
  for (candidate in candidates) {
    if (file.exists(file.path(candidate, "paths.yaml"))) {
      return(normalizePath(candidate))
    }
  }
  stop("Could not find project root")
}

project_root <- find_project_root()
config <- yaml::read_yaml(file.path(project_root, "paths.yaml"))

# Paths
data_dir <- file.path(project_root, config$results$data)
figures_dir <- file.path(project_root, config$results$figures)
phylo_dir <- file.path(project_root, config$input_data$phylogeny)
output_dir <- file.path(figures_dir, "trait_evolution")
output_data_dir <- file.path(data_dir, "trait_evolution")

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

message("Loaded superfamily data: ", nrow(superfamily_df), " species, ", ncol(superfamily_df) - 1, " traits")

# --- Data Preparation ---
message("\n=== Preparing Data ===")

prepare_trait_data <- function(df, tree) {
  species_col <- names(df)[1]
  df[[species_col]] <- gsub("^D\\.", "", df[[species_col]])
  matching_species <- intersect(tree$tip.label, df[[species_col]])

  if (length(matching_species) < 4) {
    stop("Fewer than 4 matching species between tree and data")
  }

  message("  Matching species: ", length(matching_species), " of ", length(tree$tip.label))
  pruned_tree <- drop.tip(tree, setdiff(tree$tip.label, matching_species))

  df_filtered <- df %>% filter(.data[[species_col]] %in% matching_species)
  rownames(df_filtered) <- df_filtered[[species_col]]
  df_filtered <- df_filtered[pruned_tree$tip.label, ]

  trait_matrix <- df_filtered %>% select(where(is.numeric)) %>% as.data.frame()
  rownames(trait_matrix) <- pruned_tree$tip.label

  list(tree = pruned_tree, traits = trait_matrix, species = matching_species)
}

superfamily_data <- prepare_trait_data(superfamily_df, tree)
order_data <- prepare_trait_data(order_df, tree)
class_data <- prepare_trait_data(class_df, tree)

message("Prepared data for ", length(superfamily_data$species), " species")

# --- Phylogenetic Signal Analysis using Pagel's Lambda ---
message("\n=== Analyzing Phylogenetic Signal ===")

analyze_trait <- function(tree, trait_vector, trait_name) {
  # Remove NA values
  valid_idx <- !is.na(trait_vector)
  if (sum(valid_idx) < 4 || var(trait_vector[valid_idx]) == 0) {
    return(tibble(
      trait = trait_name, n_species = sum(valid_idx),
      lambda = NA_real_, lambda_p = NA_real_,
      K = NA_real_, K_p = NA_real_,
      pic_var = NA_real_, ancestral_root = NA_real_
    ))
  }

  if (sum(!valid_idx) > 0) {
    tree_pruned <- drop.tip(tree, names(trait_vector)[!valid_idx])
    trait_vector <- trait_vector[valid_idx]
  } else {
    tree_pruned <- tree
  }

  result <- tibble(
    trait = trait_name, n_species = length(trait_vector),
    lambda = NA_real_, lambda_p = NA_real_,
    K = NA_real_, K_p = NA_real_,
    pic_var = NA_real_, ancestral_root = NA_real_
  )

  # Pagel's Lambda using caper
  tryCatch({
    # Create comparative data
    trait_df <- data.frame(species = names(trait_vector), trait = trait_vector)
    comp_data <- comparative.data(tree_pruned, trait_df, names.col = species, vcv = TRUE)

    # Fit PGLS with ML lambda
    pgls_ml <- pgls(trait ~ 1, data = comp_data, lambda = "ML")
    result$lambda <- pgls_ml$param["lambda"]

    # Test significance by comparing to lambda=0
    pgls_0 <- pgls(trait ~ 1, data = comp_data, lambda = 0)
    lrt <- 2 * (logLik(pgls_ml) - logLik(pgls_0))
    result$lambda_p <- pchisq(as.numeric(lrt), df = 1, lower.tail = FALSE)
  }, error = function(e) {
    message("    Lambda estimation failed for ", trait_name)
  })

  # Blomberg's K
  tryCatch({
    # Calculate phylogenetically independent contrasts
    pics <- pic(trait_vector, tree_pruned)
    result$pic_var <- var(pics)

    # Calculate K
    n <- length(trait_vector)
    V <- vcv.phylo(tree_pruned)

    # Expected variance under BM
    C_mean <- mean(V[lower.tri(V)])
    V_obs <- var(trait_vector)
    V_bm <- sum(diag(V)) / n - C_mean

    # MSE from PICs
    mse_pics <- mean(pics^2)
    mse_bm <- V_bm

    result$K <- (V_obs / C_mean) / (mse_pics / mse_bm)

    # Randomization test for K
    n_perm <- 999
    K_null <- numeric(n_perm)
    for (i in 1:n_perm) {
      trait_rand <- sample(trait_vector)
      names(trait_rand) <- names(trait_vector)
      pics_rand <- pic(trait_rand, tree_pruned)
      mse_rand <- mean(pics_rand^2)
      V_obs_rand <- var(trait_rand)
      K_null[i] <- (V_obs_rand / C_mean) / (mse_rand / mse_bm)
    }
    result$K_p <- mean(K_null >= result$K)
  }, error = function(e) {
    message("    K estimation failed for ", trait_name)
  })

  # Ancestral state at root
  tryCatch({
    anc <- ace(trait_vector, tree_pruned, type = "continuous", method = "ML")
    result$ancestral_root <- anc$ace[1]
  }, error = function(e) {
    message("    Ancestral reconstruction failed for ", trait_name)
  })

  result
}

# Analyze all traits
message("\nAnalyzing Superfamily traits...")
superfamily_results <- lapply(names(superfamily_data$traits), function(trait_name) {
  message("  Processing: ", trait_name)
  trait_vec <- superfamily_data$traits[[trait_name]]
  names(trait_vec) <- rownames(superfamily_data$traits)
  analyze_trait(superfamily_data$tree, trait_vec, trait_name)
})
superfamily_results_df <- bind_rows(superfamily_results) %>% mutate(level = "superfamily")

message("\nAnalyzing Order traits...")
order_results <- lapply(names(order_data$traits), function(trait_name) {
  message("  Processing: ", trait_name)
  trait_vec <- order_data$traits[[trait_name]]
  names(trait_vec) <- rownames(order_data$traits)
  analyze_trait(order_data$tree, trait_vec, trait_name)
})
order_results_df <- bind_rows(order_results) %>% mutate(level = "order")

message("\nAnalyzing Class traits...")
class_results <- lapply(names(class_data$traits), function(trait_name) {
  message("  Processing: ", trait_name)
  trait_vec <- class_data$traits[[trait_name]]
  names(trait_vec) <- rownames(class_data$traits)
  analyze_trait(class_data$tree, trait_vec, trait_name)
})
class_results_df <- bind_rows(class_results) %>% mutate(level = "class")

# Combine all results
all_results <- bind_rows(superfamily_results_df, order_results_df, class_results_df)

# Classify evolutionary pattern
all_results <- all_results %>%
  mutate(
    phylo_signal = case_when(
      is.na(lambda) ~ "unknown",
      lambda > 0.8 & lambda_p < 0.05 ~ "strong (BM-like)",
      lambda > 0.4 & lambda_p < 0.05 ~ "moderate",
      lambda_p < 0.05 ~ "weak",
      TRUE ~ "none detected"
    )
  )

# --- Save Results ---
message("\n=== Saving Results ===")
output_file <- file.path(output_data_dir, "evolutionary_model_comparison.csv")
write_csv(all_results, output_file)
message("Saved results to: ", output_file)

# --- Generate Plots ---
message("\n=== Generating Plots ===")

# Plot 1: Lambda distribution
p_lambda <- ggplot(all_results %>% filter(!is.na(lambda)),
                   aes(x = lambda, fill = level)) +
  geom_histogram(bins = 20, alpha = 0.7, position = "identity") +
  geom_vline(xintercept = c(0, 1), linetype = "dashed", color = c("red", "blue")) +
  scale_fill_brewer(palette = "Set1", name = "Level") +
  labs(
    title = "Distribution of Pagel's Lambda Across TE Categories",
    subtitle = "Red = no signal (lambda=0), Blue = Brownian motion (lambda=1)",
    x = expression("Pagel's " * lambda),
    y = "Count"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave(file.path(output_dir, "lambda_distribution.png"), p_lambda,
       width = 10, height = 6, dpi = 300)

# Plot 2: Lambda by trait (for Order level)
order_plot_data <- all_results %>%
  filter(level == "order", !is.na(lambda)) %>%
  arrange(desc(lambda))

p_order_lambda <- ggplot(order_plot_data,
                         aes(x = reorder(trait, lambda), y = lambda, fill = phylo_signal)) +
  geom_bar(stat = "identity") +
  geom_hline(yintercept = c(0, 1), linetype = "dashed", color = "gray50") +
  coord_flip() +
  scale_fill_brewer(palette = "Set2", name = "Signal") +
  labs(
    title = "Phylogenetic Signal in TE Order Composition",
    x = "TE Order",
    y = expression("Pagel's " * lambda)
  ) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave(file.path(output_dir, "lambda_by_order.png"), p_order_lambda,
       width = 10, height = 6, dpi = 300)

# Plot 3: K values
p_K <- all_results %>%
  filter(!is.na(K), level == "superfamily") %>%
  ggplot(aes(x = reorder(trait, K), y = K, fill = K_p < 0.05)) +
  geom_bar(stat = "identity") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "blue") +
  coord_flip() +
  scale_fill_manual(values = c("TRUE" = "steelblue", "FALSE" = "gray70"),
                   name = "Significant", labels = c("No", "Yes")) +
  labs(
    title = "Blomberg's K for TE Superfamily Composition",
    subtitle = "K=1 indicates Brownian motion; K<1 less signal than expected",
    x = "TE Superfamily",
    y = "Blomberg's K"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.y = element_text(size = 7)
  )

ggsave(file.path(output_dir, "K_by_superfamily.png"), p_K,
       width = 12, height = 10, dpi = 300)

# Plot 4: Summary pie chart
signal_summary <- all_results %>%
  filter(!is.na(phylo_signal)) %>%
  count(phylo_signal) %>%
  mutate(pct = round(100 * n / sum(n), 1))

p_summary <- ggplot(signal_summary, aes(x = "", y = n, fill = phylo_signal)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y") +
  scale_fill_brewer(palette = "Set2", name = "Signal Strength") +
  labs(title = "Phylogenetic Signal Distribution Across All TE Traits") +
  theme_void() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave(file.path(output_dir, "signal_summary_pie.png"), p_summary,
       width = 8, height = 6, dpi = 300)

# --- Print Summary ---
message("\n=== Summary ===")
message("Total traits analyzed: ", nrow(all_results))

signal_table <- all_results %>%
  filter(!is.na(phylo_signal)) %>%
  count(phylo_signal)
message("\nPhylogenetic signal distribution:")
print(signal_table)

strong_signal <- all_results %>%
  filter(phylo_signal == "strong (BM-like)") %>%
  select(trait, level, lambda, K)

if (nrow(strong_signal) > 0) {
  message("\nTraits with strong phylogenetic signal (BM-like evolution):")
  print(strong_signal)
}

message("\n=== Trait Evolution Analysis Complete ===")
