# Downstream exploratory comparison of released species-level summaries.
# No individual element is treated as an independent regression observation.

prepare_ltr_genome_data <- function(dataset_dir) {
  paths <- file.path(dataset_dir, c(
    "Supplementary_Data_S15_ltr_terminal_internal_species_robustness.csv",
    "Supplementary_Data_S04_path_analysis_traits24.csv",
    "trees/desmo900dated_source46.tre"
  ))
  ltr <- readr::read_csv(paths[[1]], show_col_types = FALSE)
  ltr <- ltr[ltr$analysis_branch == "primary_iqr_filtered_coverage_ge_80pct", ]
  traits <- readr::read_csv(paths[[2]], show_col_types = FALSE)
  tree <- ape::ladderize(ape::read.tree(paths[[3]]), right = FALSE)
  if (nrow(ltr) != 30L || sum(ltr$n_elements) != 380L ||
      anyDuplicated(ltr$species) || anyDuplicated(traits$tree_tip) ||
      ape::Ntip(tree) != 46L || anyDuplicated(tree$tip.label) ||
      !ape::is.rooted(tree) || !ape::is.ultrametric(tree, tol = 1e-6) ||
      any(!is.finite(tree$edge.length)) || any(tree$edge.length < 0)) {
    stop("Expected the released 30-species / 380-element branch and dated source tree")
  }
  canonical_order <- tree$tip.label[order(ape::node.height(tree)[seq_len(ape::Ntip(tree))])]
  current <- ltr[, c("species", "n_elements", "ratio_median",
                    "median_bootstrap_ci_low", "median_bootstrap_ci_high")]
  trait_index <- match(current$species, traits$tree_tip)
  # Shared conversion changes units only; no new genome-size calibration is fit.
  current$genome_size_gbp <- traits$genome_size_pg_fuscus_anchored[trait_index] * 0.978
  current$genome_size_ci_low_gbp <- traits$genome_size_pg_ci_low[trait_index] * 0.978
  current$genome_size_ci_high_gbp <- traits$genome_size_pg_ci_high[trait_index] * 0.978
  excluded <- current[is.na(trait_index), c("species", "n_elements")]
  excluded$reason <- "No released genome-size estimate in S04"
  current <- current[!is.na(trait_index), ]
  current <- current[order(match(current$species, canonical_order)), ]
  if (nrow(current) != 19L || nrow(excluded) != 11L ||
      anyNA(current) || any(!is.finite(as.matrix(current[, -1]))) ||
      any(as.matrix(current[, -1]) <= 0) ||
      !all(current$species %in% tree$tip.label)) {
    stop("Expected 19 complete positive species summaries; no imputation is permitted")
  }
  current$log10_ratio <- log10(current$ratio_median)
  current$log10_genome <- log10(current$genome_size_gbp)
  current <- as.data.frame(current)
  rownames(current) <- current$species
  list(data = current, excluded = excluded,
       tree = ape::keep.tip(tree, current$species),
       canonical_order = canonical_order, input_paths = paths)
}

fit_ltr_genome_model <- function(data, tree, covariance = "lambda", raw = FALSE,
                               model_id = "log_log_lambda") {
  formula <- if (raw) ratio_median ~ genome_size_gbp else log10_ratio ~ log10_genome
  fit_warnings <- character()
  fit <- withCallingHandlers({
    if (covariance == "OLS") {
      stats::lm(formula, data = data)
    } else if (covariance == "lambda") {
      phylolm::phylolm(formula, data = data, phy = tree, model = "lambda",
                      lower.bound = 0, upper.bound = 1, REML = FALSE)
    } else {
      phylolm::phylolm(formula, data = data, phy = tree, model = covariance, REML = FALSE)
    }
  }, warning = function(w) {
    fit_warnings <<- c(fit_warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  coefficients <- stats::coef(fit)
  slope <- unname(coefficients[[2]])
  standard_error <- sqrt(diag(stats::vcov(fit)))[[2]]
  df <- nrow(data) - length(coefficients)
  critical <- stats::qt(0.975, df)
  lambda <- if (covariance == "lambda") fit$optpar else if (covariance == "BM") 1 else 0
  row <- tibble::tibble(
    model = model_id, transformation = if (raw) "untransformed" else "log10-log10",
    n_species = nrow(data), intercept = unname(coefficients[[1]]),
    slope = slope, standard_error = standard_error,
    ci_low = slope - critical * standard_error,
    ci_high = slope + critical * standard_error,
    t_value = slope / standard_error,
    p_value = 2 * stats::pt(-abs(slope / standard_error), df = df),
    lambda = as.numeric(lambda),
    lambda_boundary = covariance == "lambda" && (lambda < 1e-6 || lambda > 1 - 1e-6),
    warnings = paste(unique(fit_warnings), collapse = " | ")
  )
  if (any(!is.finite(unlist(row[c("slope", "standard_error", "ci_low", "ci_high", "p_value", "lambda")])))) {
    stop("Non-finite PGLS coefficient or uncertainty")
  }
  list(fit = fit, summary = row)
}

fit_ltr_genome_pgls <- function(prepared) {
  primary <- fit_ltr_genome_model(prepared$data, prepared$tree)
  bm <- fit_ltr_genome_model(prepared$data, prepared$tree, "BM", model_id = "log_log_BM")
  ols <- fit_ltr_genome_model(prepared$data, prepared$tree, "OLS", model_id = "log_log_OLS")
  raw <- fit_ltr_genome_model(prepared$data, prepared$tree, raw = TRUE,
                            model_id = "raw_raw_lambda")
  loo <- lapply(prepared$data$species, function(omitted) {
    current <- prepared$data[prepared$data$species != omitted, , drop = FALSE]
    fit <- fit_ltr_genome_model(current, ape::keep.tip(prepared$tree, current$species))
    dplyr::mutate(fit$summary, omitted_species = omitted, .before = 1)
  })
  x <- seq(min(prepared$data$log10_genome), max(prepared$data$log10_genome), length.out = 200)
  design <- cbind(1, x)
  estimate <- as.numeric(design %*% stats::coef(primary$fit))
  se <- sqrt(rowSums((design %*% stats::vcov(primary$fit)) * design))
  critical <- stats::qt(.975, nrow(prepared$data) - 2L)
  # Back transformation is the fitted geometric mean; no arithmetic-mean
  # bias correction or tree-conditioned new-species prediction is implied.
  prediction <- tibble::tibble(
    genome_size_gbp = 10^x, ratio_fit = 10^estimate,
    ratio_ci_low = 10^(estimate - critical * se),
    ratio_ci_high = 10^(estimate + critical * se)
  )
  list(primary = primary$fit,
       summary = dplyr::bind_rows(primary$summary, bm$summary, ols$summary, raw$summary),
       loo = dplyr::bind_rows(loo), prediction = prediction)
}
