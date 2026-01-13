#' PCA Utility Functions for Desmognathus TE Analysis
#'
#' This module provides shared functions for PCA analysis including:
#' - Configuration loading
#' - Data filtering (NA, zeros)
#' - CLR transformation
#' - Scree plot generation
#' - Clustering utilities
#'
#' @usage source("scripts/processing/pca_utils.R")

# Required libraries - load once
required_packages <- c(
  "stats",       # For prcomp, kmeans
  "dplyr",       # For data manipulation
  "ggplot2",     # For plotting
  "cluster",     # For silhouette
  "factoextra",  # For PCA visualization
  "readr",       # For read_csv
  "tidyr",       # For data tidying
  "compositions", # For CLR transformation
  "ggrepel",     # For text labels
  "yaml"         # For config loading
)

# Load required packages
invisible(lapply(required_packages, function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(paste("Package", pkg, "is required but not installed."))
  }
  library(pkg, character.only = TRUE)
}))


#' Find project root directory
#'
#' Searches upward from current directory for paths.yaml
#'
#' @return Character string of project root path
find_project_root <- function() {
  # Try to find paths.yaml starting from script location
  script_dir <- tryCatch({
    dirname(sys.frame(1)$ofile)
  }, error = function(e) {
    getwd()
  })

  current <- normalizePath(script_dir, mustWork = FALSE)

  # Search upward for paths.yaml
  for (i in 1:10) {  # Limit search depth
    if (file.exists(file.path(current, "paths.yaml"))) {
      return(current)
    }
    parent <- dirname(current)
    if (parent == current) break  # Reached root
    current <- parent
  }

  # Fallback: check common locations
  candidates <- c(
    Sys.getenv("DESMOGNATHUS_PROJECT_ROOT"),
    "~/Projects/Desmognathus_TE",
    getwd()
  )

  for (candidate in candidates) {
    if (nchar(candidate) > 0 && file.exists(file.path(candidate, "paths.yaml"))) {
      return(normalizePath(candidate))
    }
  }

  stop("Could not find project root (no paths.yaml found)")
}


#' Load project configuration from paths.yaml
#'
#' @param project_root Optional project root path
#' @return List with configuration paths
load_config <- function(project_root = NULL) {
  if (is.null(project_root)) {
    project_root <- find_project_root()
  }

  config_file <- file.path(project_root, "paths.yaml")
  if (!file.exists(config_file)) {
    stop(paste("Configuration file not found:", config_file))
  }

  config <- yaml::read_yaml(config_file)

  # Convert relative paths to absolute
  list(
    project_root = project_root,
    input_data = list(
      root = file.path(project_root, config$input_data$root),
      dnaPipeTE = file.path(project_root, config$input_data$dnaPipeTE),
      repeatmasker = file.path(project_root, config$input_data$repeatmasker),
      phylogeny = file.path(project_root, config$input_data$phylogeny),
      lookup_table = file.path(project_root, config$input_data$lookup_table)
    ),
    results = list(
      root = file.path(project_root, config$results$root),
      data = file.path(project_root, config$results$data),
      figures = file.path(project_root, config$results$figures)
    ),
    interim = list(
      root = file.path(project_root, config$interim$root)
    )
  )
}


#' Filter columns with NA values
#'
#' @param df Data frame to filter
#' @param verbose Print removed column names
#' @return Filtered data frame
filter_na_columns <- function(df, verbose = TRUE) {
  cols_before <- colnames(df)
  df_filtered <- df[, colSums(is.na(df)) == 0, drop = FALSE]
  cols_after <- colnames(df_filtered)

  if (verbose) {
    removed <- setdiff(cols_before, cols_after)
    if (length(removed) > 0) {
      message("Removed columns due to NA values: ", paste(removed, collapse = ", "))
    }
  }

  df_filtered
}


#' Filter columns containing zero values
#'
#' @param df Data frame to filter
#' @param verbose Print removed column names
#' @return Filtered data frame
filter_zero_columns <- function(df, verbose = TRUE) {
  if (ncol(df) == 0) return(df)

  cols_before <- colnames(df)
  contains_zero <- sapply(df, function(col) any(col == 0, na.rm = TRUE))
  df_filtered <- df[, !contains_zero, drop = FALSE]
  cols_after <- colnames(df_filtered)

  if (verbose) {
    removed <- setdiff(cols_before, cols_after)
    if (length(removed) > 0) {
      message("Removed columns containing zeros: ", paste(removed, collapse = ", "))
    }
  }

  df_filtered
}


#' Filter rows with NA values
#'
#' @param df Data frame to filter
#' @param id_vector Optional ID vector to filter in parallel
#' @param verbose Print number of removed rows
#' @return List with filtered df and id_vector
filter_na_rows <- function(df, id_vector = NULL, verbose = TRUE) {
  if (ncol(df) == 0) {
    return(list(df = df, ids = id_vector))
  }

  rows_before <- nrow(df)
  complete_mask <- complete.cases(df)
  df_filtered <- df[complete_mask, , drop = FALSE]

  if (!is.null(id_vector)) {
    id_vector <- id_vector[complete_mask]
  }

  if (verbose && rows_before > nrow(df_filtered)) {
    message("Removed ", rows_before - nrow(df_filtered), " rows containing NA values")
  }

  list(df = df_filtered, ids = id_vector)
}


#' Apply Centered Log-Ratio (CLR) transformation
#'
#' @param df Numeric data frame (all positive values required)
#' @return CLR-transformed data frame
apply_clr_transform <- function(df) {
  if (ncol(df) == 0 || nrow(df) == 0) {
    stop("Cannot apply CLR to empty data frame")
  }

  # Check for zeros or negative values
  if (any(df <= 0, na.rm = TRUE)) {
    stop("CLR transformation requires all positive values. Filter zeros first.")
  }

  # Apply CLR using compositions package
  clr_result <- compositions::clr(as.matrix(df))

  # Convert back to data frame
  as.data.frame(clr_result)
}


#' Prepare data for PCA
#'
#' Complete pipeline: extract numeric, filter NAs, filter zeros, apply CLR
#'
#' @param df Data frame with identifier column and numeric columns
#' @param id_col Name of identifier column
#' @param apply_clr Whether to apply CLR transformation
#' @param verbose Print progress messages
#' @return List with prepared numeric data and IDs
prepare_pca_data <- function(df, id_col = names(df)[1], apply_clr = TRUE, verbose = TRUE) {
  # Validate ID column
  if (!id_col %in% names(df)) {
    stop(paste("ID column", id_col, "not found in data frame"))
  }

  # Extract IDs
  ids <- df[[id_col]]

  # Select numeric columns only
  numeric_df <- df %>% select(where(is.numeric))

  if (verbose) {
    message("Starting with ", ncol(numeric_df), " numeric columns and ", nrow(numeric_df), " rows")
  }

  # Filter NA columns
  numeric_df <- filter_na_columns(numeric_df, verbose)

  # Filter zero columns (required for CLR)
  if (apply_clr) {
    numeric_df <- filter_zero_columns(numeric_df, verbose)
  }

  # Filter NA rows
  result <- filter_na_rows(numeric_df, ids, verbose)
  numeric_df <- result$df
  ids <- result$ids

  # Apply CLR transformation
  if (apply_clr && ncol(numeric_df) > 0) {
    if (verbose) message("Applying CLR transformation...")
    numeric_df <- apply_clr_transform(numeric_df)
  }

  if (verbose) {
    message("Final data: ", ncol(numeric_df), " columns, ", nrow(numeric_df), " rows")
  }

  list(
    data = numeric_df,
    ids = ids,
    clr_applied = apply_clr
  )
}


#' Determine number of components explaining threshold variance
#'
#' @param pca_result PCA result from prcomp
#' @param threshold Cumulative variance threshold (default 0.80)
#' @return Number of components
get_n_components <- function(pca_result, threshold = 0.80) {
  var_explained <- summary(pca_result)$importance[3, ]  # Cumulative proportion
  which(var_explained >= threshold)[1]
}


#' Create scree plot
#'
#' @param pca_result PCA result from prcomp
#' @param title Plot title
#' @param output_path Optional path to save plot
#' @return ggplot object
create_scree_plot <- function(pca_result, title = "Scree Plot", output_path = NULL) {
  var_explained <- summary(pca_result)$importance[2, ] * 100  # Proportion as percentage
  cumulative <- cumsum(var_explained)

  plot_data <- data.frame(
    PC = factor(1:length(var_explained)),
    Variance = var_explained,
    Cumulative = cumulative
  )

  # Limit to first 15 components for readability
  n_show <- min(15, nrow(plot_data))
  plot_data <- plot_data[1:n_show, ]

  p <- ggplot(plot_data, aes(x = PC)) +
    geom_bar(aes(y = Variance), stat = "identity", fill = "steelblue", alpha = 0.7) +
    geom_line(aes(y = Cumulative, group = 1), color = "red", size = 1) +
    geom_point(aes(y = Cumulative), color = "red", size = 2) +
    geom_hline(yintercept = 80, linetype = "dashed", color = "gray50") +
    scale_y_continuous(
      name = "Variance Explained (%)",
      sec.axis = sec_axis(~., name = "Cumulative (%)")
    ) +
    labs(title = title, x = "Principal Component") +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )

  if (!is.null(output_path)) {
    ggsave(output_path, p, width = 10, height = 6, dpi = 300)
    message("Saved scree plot to: ", output_path)
  }

  p
}


#' Load standard TE breakdown data
#'
#' @param config Configuration list from load_config()
#' @param level Either "order" or "superfamily"
#' @return Data frame
load_te_data <- function(config, level = c("order", "superfamily")) {
  level <- match.arg(level)

  filename <- switch(level,
    order = "dnaPipeTE_order_breakdown.csv",
    superfamily = "dnaPipeTE_superfamily_breakdown.csv"
  )

  filepath <- file.path(config$results$data, filename)

  if (!file.exists(filepath)) {
    stop(paste("Data file not found:", filepath))
  }

  read_csv(filepath, show_col_types = FALSE)
}


#' Clean species identifiers by removing 'D.' prefix
#'
#' @param df Data frame
#' @param id_col Name of ID column (default: first column)
#' @return Data frame with cleaned IDs
clean_species_ids <- function(df, id_col = names(df)[1]) {
  df[[id_col]] <- sub("^D\\.", "", df[[id_col]])
  df
}


# Print confirmation when sourced
message("PCA utilities loaded successfully")
