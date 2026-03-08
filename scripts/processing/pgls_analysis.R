#!/usr/bin/env Rscript
#' Phylogenetic Generalized Least Squares (PGLS) Analysis
#'
#' Performs phylogenetically-corrected regression analysis on TE composition data.
#' Uses caper package for PGLS with lambda estimation.
#'
#' @usage Rscript scripts/processing/pgls_analysis.R

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
  "ape",          # Phylogenetic tree handling
  "caper",        # PGLS regression
  "nlme",         # GLS functions
  "dplyr",        # Data manipulation
  "tidyr",        # Data reshaping
  "readr",        # CSV reading
  "tibble",       # Modern data frames
  "ggplot2",      # Plotting
  "broom",        # Tidy model outputs
  "yaml"          # Config loading
)

# Install/load packages
invisible(lapply(required_packages, function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
  library(pkg, character.only = TRUE)
}))

# --- Configuration ---
project_root <- find_project_root(script_dir)
config <- load_project_config(project_root)

# Paths
data_dir <- resolve_config_path(project_root, config$results$data, "results/data")
figures_dir <- resolve_config_path(project_root, config$results$figures, "results/figures")
phylo_dir <- resolve_config_path(project_root, config$input_data$phylogeny %||% config$data$phylogeny, "input_data/phylogeny")
output_dir <- file.path(figures_dir, "pgls")
output_data_dir <- file.path(data_dir, "pgls")

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

# Load diversity metrics
diversity_order_file <- file.path(data_dir, "diversity_order_stats.csv")
diversity_superfamily_file <- file.path(data_dir, "diversity_superfamily_stats.csv")

superfamily_df <- read_csv(superfamily_file, show_col_types = FALSE)
order_df <- read_csv(order_file, show_col_types = FALSE)
class_df <- read_csv(class_file, show_col_types = FALSE)

# Load diversity if available
diversity_order_df <- NULL
diversity_superfamily_df <- NULL
if (file.exists(diversity_order_file)) {
  diversity_order_df <- read_csv(diversity_order_file, show_col_types = FALSE)
  message("Loaded order diversity metrics")
}
if (file.exists(diversity_superfamily_file)) {
  diversity_superfamily_df <- read_csv(diversity_superfamily_file, show_col_types = FALSE)
  message("Loaded superfamily diversity metrics")
}

message("Loaded TE composition data for ", nrow(superfamily_df), " species")

# --- Data Preparation ---
message("\n=== Preparing Data ===")

prepare_comparative_data <- function(tree, trait_df, trait_df_name = "traits") {
  #' Prepare data for caper::comparative.data
  #'
  #' Matches species between tree and trait data, creates comparative.data object

  # Get species column
  species_col <- names(trait_df)[1]

  # Clean species names
  trait_df[[species_col]] <- gsub("^D\\.", "", trait_df[[species_col]])

  # Rename species column to 'species' for caper
  trait_df <- trait_df %>%
    rename(species = !!species_col)

  # Find matching species
  matching <- intersect(tree$tip.label, trait_df$species)

  if (length(matching) < 4) {
    stop("Fewer than 4 matching species")
  }

  message("  Matching species: ", length(matching))

  # Filter to matching species
  trait_df <- trait_df %>%
    filter(species %in% matching)

  # Prune tree
  pruned_tree <- drop.tip(tree, setdiff(tree$tip.label, matching))

  # Create comparative.data object
  comp_data <- comparative.data(
    phy = pruned_tree,
    data = as.data.frame(trait_df),
    names.col = species,
    vcv = TRUE,
    na.omit = FALSE,
    warn.dropped = TRUE
  )

  return(comp_data)
}

# --- PGLS Analysis Functions ---
message("\n=== PGLS Analysis Functions ===")

run_pgls <- function(comp_data, response_var, predictor_vars, model_name = NULL) {
  #' Run PGLS regression with lambda estimation
  #'
  #' @param comp_data comparative.data object from caper
  #' @param response_var name of response variable
  #' @param predictor_vars vector of predictor variable names
  #' @param model_name optional name for the model
  #'
  #' @return list with model, summary, and diagnostics

  quote_formula_var <- function(x) {
    paste0("`", x, "`")
  }

  # Build a formula that remains valid for names like `Tc1-mariner`.
  formula_label <- paste(response_var, "~", paste(predictor_vars, collapse = " + "))
  formula_str <- paste(
    quote_formula_var(response_var),
    "~",
    paste(vapply(predictor_vars, quote_formula_var, character(1)), collapse = " + ")
  )
  formula_obj <- as.formula(formula_str)

  if (is.null(model_name)) {
    model_name <- formula_label
  }

  # Check that variables exist
  all_vars <- c(response_var, predictor_vars)
  missing_vars <- setdiff(all_vars, names(comp_data$data))
  if (length(missing_vars) > 0) {
    warning("Variables not found: ", paste(missing_vars, collapse = ", "))
    return(NULL)
  }

  # Check for zero variance
  for (var in all_vars) {
    if (var(comp_data$data[[var]], na.rm = TRUE) == 0) {
      warning("Zero variance in: ", var)
      return(NULL)
    }
  }

  result <- tryCatch({
    # Fit PGLS with maximum likelihood estimation of lambda
    pgls_model <- pgls(
      formula_obj,
      data = comp_data,
      lambda = "ML"  # Estimate lambda via maximum likelihood
    )

    # Get summary
    model_summary <- summary(pgls_model)

    # Extract key statistics
    list(
      model_name = model_name,
      formula = formula_label,
      model = pgls_model,
      summary = model_summary,
      coefficients = as.data.frame(model_summary$coefficients),
      lambda = pgls_model$param["lambda"],
      lambda_lower = pgls_model$param.CI$lambda$ci.val[1],
      lambda_upper = pgls_model$param.CI$lambda$ci.val[2],
      r_squared = model_summary$r.squared,
      adj_r_squared = model_summary$adj.r.squared,
      aic = AIC(pgls_model),
      n = nrow(comp_data$data),
      success = TRUE
    )
  }, error = function(e) {
    warning("PGLS failed for ", model_name, ": ", e$message)
    list(
      model_name = model_name,
      formula = formula_label,
      success = FALSE,
      error = e$message
    )
  })

  return(result)
}

run_pairwise_pgls <- function(comp_data, variables, response_prefix = NULL) {
  #' Run PGLS for all pairwise combinations of variables
  #'
  #' Tests correlations between all pairs of TE categories

  results <- list()

  for (i in seq_along(variables)) {
    for (j in seq_along(variables)) {
      if (i >= j) next  # Avoid duplicates and self-comparison

      response <- variables[i]
      predictor <- variables[j]

      model_name <- paste(response, "~", predictor)
      message("  Testing: ", model_name)

      result <- run_pgls(comp_data, response, predictor, model_name)
      if (!is.null(result) && result$success) {
        results[[model_name]] <- result
      }
    }
  }

  return(results)
}

summarize_pgls_results <- function(results_list) {
  #' Create summary data frame from PGLS results

  if (length(results_list) == 0) {
    return(data.frame())
  }

  scalar_numeric <- function(x) {
    if (is.null(x) || length(x) == 0) {
      return(NA_real_)
    }
    as.numeric(x[[1]])
  }

  summaries <- lapply(results_list, function(res) {
    if (!res$success) return(NULL)

    # Get coefficient for main predictor (not intercept)
    coef_df <- res$coefficients
    predictor_row <- coef_df[rownames(coef_df) != "(Intercept)", , drop = FALSE]

    if (nrow(predictor_row) == 0) return(NULL)

    tibble::tibble(
      model = as.character(res$model_name),
      formula = as.character(res$formula),
      lambda = scalar_numeric(res$lambda),
      lambda_lower = scalar_numeric(res$lambda_lower),
      lambda_upper = scalar_numeric(res$lambda_upper),
      r_squared = scalar_numeric(res$r_squared),
      adj_r_squared = scalar_numeric(res$adj_r_squared),
      aic = scalar_numeric(res$aic),
      n = as.integer(res$n),
      predictor = gsub("^`|`$", "", rownames(predictor_row)[1]),
      estimate = as.numeric(predictor_row$Estimate[1]),
      std_error = as.numeric(predictor_row$`Std. Error`[1]),
      t_value = as.numeric(predictor_row$`t value`[1]),
      p_value = as.numeric(predictor_row$`Pr(>|t|)`[1])
    )
  })

  dplyr::bind_rows(summaries)
}

# --- Main PGLS Analyses ---
message("\n=== Running PGLS Analyses ===")

# Prepare comparative data
message("\nPreparing Order-level comparative data...")
order_comp <- prepare_comparative_data(tree, order_df, "order")

message("\nPreparing Superfamily-level comparative data...")
superfamily_comp <- prepare_comparative_data(tree, superfamily_df, "superfamily")

# Get numeric variable names
order_vars <- names(order_comp$data)[sapply(order_comp$data, is.numeric)]
order_vars <- order_vars[order_vars != "species"]
superfamily_vars <- names(superfamily_comp$data)[sapply(superfamily_comp$data, is.numeric)]
superfamily_vars <- superfamily_vars[superfamily_vars != "species"]

message("\nOrder variables: ", paste(order_vars, collapse = ", "))
message("Superfamily variables: ", paste(head(superfamily_vars, 10), collapse = ", "), "...")

# --- Analysis 1: Pairwise correlations between TE Orders ---
message("\n--- Pairwise PGLS: TE Order Correlations ---")
order_pairwise <- run_pairwise_pgls(order_comp, order_vars)
order_pairwise_summary <- summarize_pgls_results(order_pairwise)

if (nrow(order_pairwise_summary) > 0) {
  order_pairwise_summary <- order_pairwise_summary %>%
    arrange(p_value) %>%
    mutate(
      p_adjusted = p.adjust(p_value, method = "BH"),
      significant = p_adjusted < 0.05
    )

  message("\nSignificant Order correlations (p.adj < 0.05):")
  sig_results <- order_pairwise_summary %>% filter(significant)
  if (nrow(sig_results) > 0) {
    print(sig_results %>% select(model, estimate, p_value, p_adjusted, lambda, r_squared))
  } else {
    message("  None found")
  }
}

# --- Analysis 2: Key hypothesis tests ---
message("\n--- Specific Hypothesis Tests ---")

# Test: Do LINE elements correlate with Gypsy LTRs?
# (Common pattern in genomes - autonomous and non-autonomous element dynamics)
if ("LINE" %in% order_vars && "LTR" %in% order_vars) {
  message("\nTesting: LINE ~ LTR")
  line_ltr <- run_pgls(order_comp, "LINE", "LTR", "LINE_vs_LTR")
  if (!is.null(line_ltr) && line_ltr$success) {
    message("  Lambda: ", round(line_ltr$lambda, 3))
    message("  R-squared: ", round(line_ltr$r_squared, 3))
    message("  P-value: ", format(line_ltr$summary$coefficients["LTR", "Pr(>|t|)"], scientific = TRUE))
  }
}

# Test: DNA transposons (TIR) vs Retrotransposons (LINE + LTR)
if (all(c("TIR", "LINE", "LTR") %in% order_vars)) {
  message("\nTesting: TIR ~ LINE + LTR")
  # Add combined retrotransposon variable
  order_comp$data$Retrotransposons <- order_comp$data$LINE + order_comp$data$LTR

  tir_retro <- run_pgls(order_comp, "TIR", "Retrotransposons", "TIR_vs_Retrotransposons")
  if (!is.null(tir_retro) && tir_retro$success) {
    message("  Lambda: ", round(tir_retro$lambda, 3))
    message("  R-squared: ", round(tir_retro$r_squared, 3))
    coef_row <- tir_retro$summary$coefficients["Retrotransposons", ]
    message("  P-value: ", format(coef_row["Pr(>|t|)"], scientific = TRUE))
  }
}

# --- Analysis 3: Superfamily correlations within LTR ---
message("\n--- Superfamily-level Correlations ---")

# Focus on major superfamilies
major_superfamilies <- c("Gypsy", "L1", "hAT", "Tc1-mariner", "Helitron", "DIRS")
available_major <- intersect(major_superfamilies, superfamily_vars)

if (length(available_major) >= 2) {
  message("Testing correlations among: ", paste(available_major, collapse = ", "))
  superfamily_pairwise <- run_pairwise_pgls(superfamily_comp, available_major)
  superfamily_pairwise_summary <- summarize_pgls_results(superfamily_pairwise)

  if (nrow(superfamily_pairwise_summary) > 0) {
    superfamily_pairwise_summary <- superfamily_pairwise_summary %>%
      arrange(p_value) %>%
      mutate(
        p_adjusted = p.adjust(p_value, method = "BH"),
        significant = p_adjusted < 0.05
      )
  }
}

# --- Analysis 4: PGLS with diversity metrics (if available) ---
if (!is.null(diversity_order_df)) {
  message("\n--- Diversity vs TE Composition ---")
  tryCatch({
    # Merge diversity with order composition
    diversity_order_df_clean <- diversity_order_df %>%
      rename(species = 1) %>%
      mutate(species = gsub("^D\\.", "", species))

    order_with_diversity <- order_comp$data %>%
      as.data.frame() %>%
      tibble::rownames_to_column("species") %>%
      left_join(diversity_order_df_clean, by = "species")

    # Test if Simpson diversity correlates with LINE content.
    if ("Simpson" %in% names(order_with_diversity) && "LINE" %in% names(order_with_diversity)) {
      # Need to rebuild comparative.data with diversity
      diversity_comp <- comparative.data(
        phy = order_comp$phy,
        data = as.data.frame(order_with_diversity),
        names.col = species,
        vcv = TRUE
      )

      message("\nTesting: Simpson ~ LINE")
      simpson_line <- run_pgls(diversity_comp, "Simpson", "LINE", "Simpson_vs_LINE")
      if (!is.null(simpson_line) && simpson_line$success) {
        message("  Lambda: ", round(simpson_line$lambda, 3))
        message("  R-squared: ", round(simpson_line$r_squared, 3))
      }
    }
  }, error = function(e) {
    message("Could not complete diversity vs composition PGLS block: ", e$message)
  })
}

# --- Save Results ---
message("\n=== Saving Results ===")

# Save order pairwise results
if (exists("order_pairwise_summary") && nrow(order_pairwise_summary) > 0) {
  write_csv(order_pairwise_summary, file.path(output_data_dir, "pgls_order_pairwise.csv"))
  message("Saved order pairwise results")
}

# Save superfamily pairwise results
if (exists("superfamily_pairwise_summary") && nrow(superfamily_pairwise_summary) > 0) {
  write_csv(superfamily_pairwise_summary, file.path(output_data_dir, "pgls_superfamily_pairwise.csv"))
  message("Saved superfamily pairwise results")
}

# --- Diagnostic Plots ---
message("\n=== Generating Diagnostic Plots ===")

# Lambda distribution plot
if (exists("order_pairwise_summary") && nrow(order_pairwise_summary) > 0) {
  p_lambda <- ggplot(order_pairwise_summary, aes(x = lambda)) +
    geom_histogram(bins = 20, fill = "steelblue", color = "black", alpha = 0.7) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
    geom_vline(xintercept = 1, linetype = "dashed", color = "blue") +
    labs(
      title = "Distribution of Pagel's Lambda in PGLS Models",
      subtitle = "Red = no phylogenetic signal, Blue = Brownian motion",
      x = expression("Pagel's " * lambda),
      y = "Count"
    ) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))

  ggsave(file.path(output_dir, "pgls_lambda_distribution.png"), p_lambda,
         width = 8, height = 6, dpi = 300)
}

# Volcano plot of PGLS results
if (exists("order_pairwise_summary") && nrow(order_pairwise_summary) > 0) {
  p_volcano <- ggplot(order_pairwise_summary,
                      aes(x = estimate, y = -log10(p_value), color = significant)) +
    geom_point(size = 3, alpha = 0.7) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray50") +
    scale_color_manual(values = c("FALSE" = "gray", "TRUE" = "red"),
                      name = "Significant\n(p.adj < 0.05)") +
    labs(
      title = "PGLS Results: TE Order Correlations",
      x = "Regression Coefficient",
      y = "-log10(p-value)"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "right"
    )

  ggsave(file.path(output_dir, "pgls_volcano_plot.png"), p_volcano,
         width = 10, height = 8, dpi = 300)
}

# R-squared vs Lambda plot
if (exists("order_pairwise_summary") && nrow(order_pairwise_summary) > 0) {
  p_rsq_lambda <- ggplot(order_pairwise_summary,
                         aes(x = lambda, y = r_squared, color = significant)) +
    geom_point(size = 3, alpha = 0.7) +
    scale_color_manual(values = c("FALSE" = "gray", "TRUE" = "red")) +
    labs(
      title = "Model Fit vs Phylogenetic Signal",
      x = expression("Pagel's " * lambda),
      y = expression(R^2)
    ) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))

  ggsave(file.path(output_dir, "pgls_rsq_vs_lambda.png"), p_rsq_lambda,
         width = 8, height = 6, dpi = 300)
}

# --- Summary ---
message("\n=== PGLS Analysis Summary ===")

if (exists("order_pairwise_summary") && nrow(order_pairwise_summary) > 0) {
  message("\nOrder-level pairwise comparisons: ", nrow(order_pairwise_summary))
  message("  Significant (p.adj < 0.05): ", sum(order_pairwise_summary$significant))
  message("  Mean lambda: ", round(mean(order_pairwise_summary$lambda, na.rm = TRUE), 3))
  message("  Lambda range: ", round(min(order_pairwise_summary$lambda, na.rm = TRUE), 3),
          " - ", round(max(order_pairwise_summary$lambda, na.rm = TRUE), 3))
}

if (exists("superfamily_pairwise_summary") && nrow(superfamily_pairwise_summary) > 0) {
  message("\nSuperfamily-level pairwise comparisons: ", nrow(superfamily_pairwise_summary))
  message("  Significant (p.adj < 0.05): ", sum(superfamily_pairwise_summary$significant))
}

message("\n=== PGLS Analysis Complete ===")
message("Results saved to: ", output_data_dir)
message("Figures saved to: ", output_dir)
