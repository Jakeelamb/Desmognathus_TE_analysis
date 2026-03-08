#!/usr/bin/env Rscript
# Description: Tests for phylogenetic signal (Pagel's Lambda and Blomberg's K)
#              in TE superfamily proportions using a given phylogeny.

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
  library(phytools)
  library(dplyr)
  library(tibble)
  library(yaml)
})

# === Configuration ===
project_root <- find_project_root(script_dir)
config <- load_project_config(project_root)

phylo_dir <- resolve_config_path(project_root, config$input_data$phylogeny %||% config$data$phylogeny, "input_data/phylogeny")
results_data_dir <- resolve_config_path(project_root, config$results$data, "results/data")
output_dir <- resolve_config_path(project_root, config$results$tables$phylo_signal, "results/tables/phylogenetic_signal")

tree_file <- file.path(phylo_dir, "desmo900dated_test.tre")
if (!file.exists(tree_file)) {
  tree_file <- file.path(project_root, "results", "phylogeny", "processed_phylogeny.nwk")
}

traits_file <- file.path(results_data_dir, "dnaPipeTE_superfamily_breakdown.csv")
if (!file.exists(traits_file)) {
  legacy_traits_dir <- resolve_config_path(project_root, config$data$processed$diversity, "data/processed/diversity")
  traits_file <- file.path(legacy_traits_dir, "superfamily_proportions.csv")
}

output_file <- file.path(output_dir, "phylogenetic_signal_results.csv")

cat("Project root:", project_root, "\n")
cat("Ensuring output directory exists:", output_dir, "\n")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# === Load Data ===
cat("Loading phylogenetic tree from:", tree_file, "\n")
if (!file.exists(tree_file)) stop("Phylogeny file not found: ", tree_file)
tree <- read.tree(tree_file)

cat("Loading trait data from:", traits_file, "\n")
if (!file.exists(traits_file)) stop("Traits file not found: ", traits_file)
# Assuming the first column is the species name/identifier
traits_raw <- read.csv(traits_file, header = TRUE, check.names = FALSE)

if (ncol(traits_raw) > 0 && !nzchar(colnames(traits_raw)[1])) {
  colnames(traits_raw)[1] <- "species"
}

# Check if the first column is 'Species' or similar and set it as row names
if (ncol(traits_raw) > 0 && tolower(colnames(traits_raw)[1]) == "species") {
  traits <- as.data.frame(traits_raw, check.names = FALSE) %>%
    tibble::column_to_rownames(var = colnames(traits_raw)[1])
  cat("Using column '", colnames(traits_raw)[1], "' as species identifiers.\n")
} else if (ncol(traits_raw) > 0) {
   warning("First column is not named 'Species'. Assuming it contains species identifiers and setting as row names.")
   traits <- as.data.frame(traits_raw, check.names = FALSE) %>%
     tibble::column_to_rownames(var = colnames(traits_raw)[1])
} else {
  stop("Trait data file is empty or has no columns.")
}

rownames(traits) <- gsub("^D\\.", "", rownames(traits))


# Convert trait data to numeric, coercing errors to NA
# Store original column names before manipulation
original_colnames <- colnames(traits)
numeric_traits <- as.data.frame(lapply(traits, function(x) as.numeric(as.character(x))), check.names = FALSE)
rownames(numeric_traits) <- rownames(traits)

# Check which columns became all NA after conversion (potential non-numeric columns)
all_na_cols <- names(which(sapply(numeric_traits, function(col) all(is.na(col)))) & !sapply(traits, function(col) all(is.na(col))))
if (length(all_na_cols) > 0) {
  warning("Columns possibly contained non-numeric data and were coerced to NA: ", paste(all_na_cols, collapse=", "))
}
traits <- numeric_traits

# === Data Matching ===
# Check for species present in the tree but not in the trait data
missing_in_traits <- tree$tip.label[!tree$tip.label %in% rownames(traits)]
if (length(missing_in_traits) > 0) {
  cat("Warning: Species in tree but not in traits file:\n")
  print(missing_in_traits)
}

# Check for species present in the trait data but not in the tree
missing_in_tree <- rownames(traits)[!rownames(traits) %in% tree$tip.label]
if (length(missing_in_tree) > 0) {
  cat("Warning: Species in traits file but not in tree:\n")
  print(missing_in_tree)
}

# Prune tree to match species present in both datasets
species_to_keep <- intersect(tree$tip.label, rownames(traits))
if (length(species_to_keep) < length(tree$tip.label)) {
  cat("Pruning tree to match species in trait data...\n")
  tree <- drop.tip(tree, setdiff(tree$tip.label, species_to_keep))
}
if (length(species_to_keep) < 2) { # Need at least 2 species for phylosig
  stop("Less than 2 matching species found between the tree and trait data. Cannot perform analysis.")
}
cat("Proceeding with", length(tree$tip.label), "matching species.\n")

# Reorder trait data rows to match the order of tips in the pruned tree
traits <- traits[tree$tip.label, , drop = FALSE]

# === Phylogenetic Signal Analysis ===
cat("Running phylogenetic signal tests for each trait...\n")

# Use lapply to iterate through each column (trait) in the traits dataframe
results <- lapply(colnames(traits), function(trait_name) {
  cat("  Processing trait:", trait_name, "\n")
  # Extract the trait data as a named vector
  trait_vector <- traits[[trait_name]]
  names(trait_vector) <- rownames(traits) # Ensure names are assigned

  # Check for variance in the trait data, ignoring NAs
  valid_values <- trait_vector[!is.na(trait_vector)]
  if (length(valid_values) < 2 || var(valid_values) == 0) {
      cat("    Skipping trait", trait_name, "- zero variance, all NAs, or less than 2 valid values.\n")
      return(tibble(
          trait = trait_name,
          lambda = NA_real_,
          lambda_p = NA_real_,
          K = NA_real_,
          K_p = NA_real_,
          error = "Zero variance, all NA, or < 2 valid values"
      ))
  }

  # Initialize result tibble for this trait
  result_tibble <- tibble(
      trait = trait_name,
      lambda = NA_real_,
      lambda_p = NA_real_,
      K = NA_real_,
      K_p = NA_real_,
      error = NA_character_
  )

  # Calculate Pagel's Lambda
  lambda_result <- tryCatch({
      phylosig(tree, trait_vector, method = "lambda", test = TRUE)
  }, error = function(e) {
      cat("    Error calculating Lambda for", trait_name, ":", conditionMessage(e), "\n")
      # Return a list matching the expected structure even on error
      list(lambda = NA_real_, P = NA_real_, error = conditionMessage(e))
  })
  # Ensure list elements exist even if tryCatch returned NULL or unexpected structure
   if(is.null(lambda_result$lambda)) lambda_result$lambda <- NA_real_
   if(is.null(lambda_result$P)) lambda_result$P <- NA_real_


  # Calculate Blomberg's K
  k_result <- tryCatch({
      phylosig(tree, trait_vector, method = "K", test = TRUE)
  }, error = function(e) {
      cat("    Error calculating K for", trait_name, ":", conditionMessage(e), "\n")
      list(K = NA_real_, P = NA_real_, error = conditionMessage(e))
  })
   if(is.null(k_result$K)) k_result$K <- NA_real_
   if(is.null(k_result$P)) k_result$P <- NA_real_


  # Store results
  result_tibble$lambda <- lambda_result$lambda
  result_tibble$lambda_p <- lambda_result$P
  result_tibble$K <- k_result$K
  result_tibble$K_p <- k_result$P
  # Combine error messages if errors occurred
  error_msgs <- na.omit(c(lambda_result$error, k_result$error))
  result_tibble$error <- if(length(error_msgs) > 0) paste(error_msgs, collapse="; ") else NA_character_


  return(result_tibble)
})

# Combine results from all traits into a single dataframe
signal_df <- bind_rows(results)

# === Output Results ===
cat("Phylogenetic signal analysis complete.\n")
print(signal_df)

# Save results to CSV
cat("Saving results to:", output_file, "\n")
tryCatch({
  write.csv(signal_df, output_file, row.names = FALSE)
}, error = function(e){
  cat("Error saving results file:", conditionMessage(e), "\n")
})


cat("Script finished.\n")
