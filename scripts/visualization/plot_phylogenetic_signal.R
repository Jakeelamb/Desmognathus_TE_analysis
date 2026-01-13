#!/usr/bin/env Rscript
# Description: Plots phylogenetic signal results (Blomberg's K and Pagel's Lambda)
#              from the analysis output CSV.

# Load required libraries
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr) # Using readr for robust CSV reading
  library(here)
  library(yaml)
  library(forcats) # For reordering factor levels
  library(stringr) # For text manipulation if needed
})

# === Configuration ===
project_root <- here::here()
config_file <- file.path(project_root, "config/paths.yaml")

# Set default paths first
default_input_subdir <- "results/tables/phylogenetic_signal"
default_output_subdir <- "results/figures/phylo_signal" # Save alongside other phylogeny figs

# Initialize with defaults
input_dir <- file.path(project_root, default_input_subdir)
output_dir <- file.path(project_root, default_output_subdir)
input_file <- file.path(input_dir, "phylogenetic_signal_results.csv")

if (file.exists(config_file)) {
  tryCatch({
    config <- yaml::read_yaml(config_file)
    cat("Config file found. Reading paths.\n")

    # Safely get paths from config
    input_dir_from_config <- config$results$tables$phylo_signal
    if (!is.null(input_dir_from_config) && is.character(input_dir_from_config) && nzchar(input_dir_from_config)) {
      input_dir <- file.path(project_root, input_dir_from_config)
      input_file <- file.path(input_dir, "phylogenetic_signal_results.csv") # Reconstruct full input path
      cat("  Using input directory from config:", input_dir_from_config, "\n")
    } else {
      warning("  Config missing or invalid 'results$tables$phylo_signal'. Using default input dir.")
      # Default input_file is already set
    }

    output_dir_from_config <- config$results$figures$phylo_signal
    if (!is.null(output_dir_from_config) && is.character(output_dir_from_config) && nzchar(output_dir_from_config)) {
      output_dir <- file.path(project_root, output_dir_from_config)
      cat("  Using output directory from config:", output_dir_from_config, "\n")
    } else {
      warning("  Config missing or invalid 'results$figures$phylo_signal'. Using default output dir.")
      # Default output_dir is already set
    }

  }, error = function(e) {
    warning("Error reading config file '", config_file, "': ", conditionMessage(e), ". Using default paths.")
    # Ensure defaults are set
    input_dir <- file.path(project_root, default_input_subdir)
    output_dir <- file.path(project_root, default_output_subdir)
    input_file <- file.path(input_dir, "phylogenetic_signal_results.csv")
  })
} else {
  warning("config/paths.yaml not found. Using default paths.")
  # Defaults are already set
}

# Ensure output directory exists
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# === Load Data ===
cat("Loading phylogenetic signal results from:", input_file, "\n")
if (!file.exists(input_file)) stop("Input file not found: ", input_file)
signal_results <- readr::read_csv(input_file, show_col_types = FALSE)

# === Data Processing ===
# Filter out rows where K or Lambda calculation failed (NA) or had errors reported
# For now, just filter NAs in the value itself and where an error was logged
k_data <- signal_results %>%
  filter(!is.na(K), is.na(error)) %>% # Keep only rows with valid K and no error
  select(trait, K, K_p) %>%
  mutate(K_significant = ifelse(K_p < 0.05, TRUE, FALSE))


lambda_data <- signal_results %>%
  filter(!is.na(lambda), is.na(error)) %>% # Keep only rows with valid lambda and no error
  select(trait, lambda, lambda_p) %>%
  mutate(lambda_significant = ifelse(lambda_p < 0.05, TRUE, FALSE))


# === Plotting Function ===
create_signal_plot <- function(data, value_col, p_col, title, y_label, h_lines = NULL, filename) {
  # Reorder traits based on the signal value
  data <- data %>%
    mutate(trait = fct_reorder(trait, {{ value_col }}, .desc = TRUE)) %>%
    # Clean up trait names for display if needed (e.g., replace underscores)
    mutate(trait_display = str_replace_all(trait, "_", " ")) %>%
    mutate(trait_display = fct_reorder(trait_display, {{ value_col }}, .desc = TRUE))

  # Optional: Determine significance status for fill color
  # This requires the significance column (e.g., K_significant) to exist
  significance_col_name <- paste0(deparse(substitute(value_col)), "_significant") # e.g., "K_significant"
  use_fill <- significance_col_name %in% colnames(data)


  plot <- ggplot(data, aes(x = trait_display, y = {{ value_col }}))

  if (use_fill) {
      plot <- plot + geom_col(aes(fill = .data[[significance_col_name]]), color = "black", width = 0.7) +
                     scale_fill_manual(values = c("TRUE" = "firebrick", "FALSE" = "steelblue"),
                                       name = "Significant (p < 0.05)",
                                       labels = c("TRUE" = "Yes", "FALSE" = "No"))
  } else {
      plot <- plot + geom_col(fill = "steelblue", color = "black", width = 0.7)
  }

  plot <- plot +
    labs(
      title = title,
      x = "TE Superfamily",
      y = y_label
    ) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 60, hjust = 1, size = 10), # Rotate labels
      plot.title = element_text(hjust = 0.5, face = "bold"),
      panel.grid.major.x = element_blank(), # Remove vertical grid lines
      panel.grid.minor.y = element_blank(),
      legend.position = "bottom" # Move legend if using fill
    )

  # Add horizontal lines if specified
  if (!is.null(h_lines)) {
    for (line_val in h_lines) {
      plot <- plot + geom_hline(yintercept = line_val, linetype = "dashed", color = "black", size = 0.8)
    }
  }

  # Save the plot
  ggsave(filename, plot = plot, width = 14, height = 8, dpi = 300, bg = "white")
  cat("Saved plot:", filename, "\n")

  return(plot) # Return the ggplot object
}

# === Generate Plots ===

# Plot Blomberg's K
if (nrow(k_data) > 0) {
  cat("Generating Blomberg's K plot...\n")
  k_plot_file <- file.path(output_dir, "blombergs_k_signal.png")
  create_signal_plot(
    data = k_data,
    value_col = K,
    p_col = K_p, # Pass p-value col name for potential fill logic
    title = "Phylogenetic Signal (Blomberg's K) in TE Superfamily Proportions",
    y_label = "Blomberg's K Value",
    h_lines = 1, # Dashed line at K=1
    filename = k_plot_file
  )
} else {
  cat("Skipping Blomberg's K plot - no valid data points.\n")
}


# Plot Pagel's Lambda
if (nrow(lambda_data) > 0) {
  cat("Generating Pagel's Lambda plot...\n")
  lambda_plot_file <- file.path(output_dir, "pagels_lambda_signal.png")
  create_signal_plot(
    data = lambda_data,
    value_col = lambda,
    p_col = lambda_p, # Pass p-value col name for potential fill logic
    title = "Phylogenetic Signal (Pagel's Lambda) in TE Superfamily Proportions",
    y_label = "Pagel's Lambda Value",
    h_lines = c(0, 1), # Dashed lines at Lambda=0 and Lambda=1
    filename = lambda_plot_file
  )
} else {
  cat("Skipping Pagel's Lambda plot - no valid data points.\n")
}


cat("Plotting script finished.\n")
