#!/usr/bin/env Rscript

# Exhaustive three-trait phylogenetic path analysis for the frozen microscopy panel.
#
# The analysis ranks all ten testable Markov-equivalence classes on relative
# nuclear IOD, nucleus size, and cell size. The eleventh (saturated) class has no d-separation
# claim and therefore cannot be scored by confirmatory phylogenetic path analysis.

prefer_active_conda_r_library <- function() {
  conda_prefix <- Sys.getenv("CONDA_PREFIX", "")
  if (nzchar(conda_prefix)) {
    conda_lib <- normalizePath(
      file.path(conda_prefix, "lib", "R", "library"),
      mustWork = FALSE
    )
    if (dir.exists(conda_lib)) .libPaths(conda_lib)
  }
}
prefer_active_conda_r_library()

suppressPackageStartupMessages({
  library(ape)
  library(dplyr)
  library(jsonlite)
  library(phylolm)
  library(phylopath)
  library(readr)
  library(tibble)
})

parse_args <- function(args) {
  options <- list(
    measurement_limit = Inf,
    tree_limit = 200L,
    simulation_replicates = 100L,
    simulation_regimes = c("iid", "brownian"),
    seed = 20260710L
  )
  index <- 1L
  while (index <= length(args)) {
    argument <- args[[index]]
    if (argument == "--quick") {
      options$measurement_limit <- 12L
      options$tree_limit <- 12L
      options$simulation_replicates <- 5L
    } else if (argument %in% c(
      "--measurement-limit", "--tree-limit", "--simulation-replicates", "--seed"
    )) {
      if (index == length(args)) stop(argument, " requires a value")
      index <- index + 1L
      value <- as.integer(args[[index]])
      if (!is.finite(value) || value < 1L) stop(argument, " must be a positive integer")
      name <- switch(
        argument,
        "--measurement-limit" = "measurement_limit",
        "--tree-limit" = "tree_limit",
        "--simulation-replicates" = "simulation_replicates",
        "--seed" = "seed"
      )
      options[[name]] <- value
    } else if (argument == "--simulation-regimes") {
      if (index == length(args)) stop(argument, " requires a value")
      index <- index + 1L
      values <- strsplit(args[[index]], ",", fixed = TRUE)[[1]]
      values <- trimws(values)
      if (!length(values) || any(!values %in% c("iid", "brownian"))) {
        stop("--simulation-regimes must contain iid and/or brownian")
      }
      options$simulation_regimes <- unique(values)
    } else {
      stop("Unknown argument: ", argument)
    }
    index <- index + 1L
  }
  options
}

raw_args <- commandArgs(trailingOnly = TRUE)
options <- parse_args(raw_args)
quick_mode <- "--quick" %in% raw_args
project_root <- normalizePath(getwd(), mustWork = TRUE)
input_dir <- file.path(
  project_root,
  "analyses/03_phylogenetic_path/data"
)
output_dir <- file.path(project_root, "analyses/03_phylogenetic_path/output")

trait_path <- file.path(input_dir, "path24_traits.csv")
measurement_path <- file.path(input_dir, "measurement_bootstrap_traits.csv")
focal_tree_path <- file.path(
  project_root, "analyses/03_phylogenetic_path/trees/source_time_tree_46.tre"
)
published_main_path <- file.path(
  project_root, "analyses/03_phylogenetic_path/trees/published_main_analysis18.nwk"
)
published_bootstrap_path <- file.path(
  project_root,
  "analyses/03_phylogenetic_path/trees/published_bootstrap_200_analysis18.nex"
)

required_paths <- c(
  trait_path,
  measurement_path,
  focal_tree_path,
  published_main_path,
  published_bootstrap_path
)
missing_paths <- required_paths[!file.exists(required_paths)]
if (length(missing_paths)) {
  stop("Missing required path-analysis input(s): ", paste(missing_paths, collapse = ", "))
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

alpha <- 0.05
competitive_delta <- 2
node_labels <- c(i = "relative_iod", n = "nucleus_size", c = "cell_size")
model_order <- c(
  "independent",
  "iod_nucleus_only",
  "nucleus_cell_only",
  "iod_cell_only",
  "nucleus_bridge",
  "nucleus_collider",
  "cell_bridge",
  "cell_collider",
  "iod_bridge",
  "iod_collider"
)

zscore <- function(values) {
  values <- as.numeric(values)
  spread <- stats::sd(values)
  if (!is.finite(spread) || spread <= 0) stop("Cannot standardize a zero-variance trait")
  as.numeric((values - mean(values)) / spread)
}

normalize_tree_labels <- function(tree) {
  tree$tip.label <- trimws(sub("^D\\.\\s*", "", tree$tip.label))
  tree
}

prune_tree <- function(tree, tips) {
  tree <- normalize_tree_labels(tree)
  if (anyDuplicated(tree$tip.label)) stop("Tree contains duplicate tip labels")
  missing <- setdiff(tips, tree$tip.label)
  if (length(missing)) stop("Tree is missing species: ", paste(missing, collapse = ", "))
  tree <- drop.tip(tree, setdiff(tree$tip.label, tips))
  if (length(tree$tip.label) != length(tips)) stop("Pruned tree does not match analysis data")
  if (is.null(tree$edge.length) || any(!is.finite(tree$edge.length))) {
    stop("Analysis tree has missing or non-finite branch lengths")
  }
  # Zero root-adjacent lengths are tolerated by ape/phylolm; negative lengths are not.
  if (any(tree$edge.length < 0)) stop("Analysis tree has negative branch lengths")
  tree
}

validate_traits <- function(source, expected_species = NULL) {
  required <- c(
    "species", "tree_tip", "relative_iod_index", "nucleus_area_um2", "cell_area_um2"
  )
  missing <- setdiff(required, colnames(source))
  if (length(missing)) stop("Trait input is missing columns: ", paste(missing, collapse = ", "))
  if (anyDuplicated(source$tree_tip)) stop("Trait input contains duplicate species")
  numeric <- source[, c("relative_iod_index", "nucleus_area_um2", "cell_area_um2")]
  if (any(!is.finite(as.matrix(numeric))) || any(as.matrix(numeric) <= 0)) {
    stop("All modeled trait values must be finite and positive")
  }
  if (!is.null(expected_species) && !setequal(source$tree_tip, expected_species)) {
    stop("Trait input species do not match the expected frozen panel")
  }
  invisible(source)
}

prepare_model_data <- function(source, tree) {
  source <- as.data.frame(source, stringsAsFactors = FALSE)
  validate_traits(source)
  current_tree <- prune_tree(tree, source$tree_tip)
  data <- data.frame(
    i = zscore(log10(source$relative_iod_index)),
    n = zscore(log10(source$nucleus_area_um2)),
    c = zscore(log10(source$cell_area_um2)),
    row.names = source$tree_tip,
    check.names = FALSE
  )
  data <- data[current_tree$tip.label, , drop = FALSE]
  list(data = data, tree = current_tree)
}

exhaustive_model_set <- function() {
  define_model_set(
    independent = c(),
    iod_nucleus_only = c(n ~ i),
    nucleus_cell_only = c(c ~ n),
    iod_cell_only = c(c ~ i),
    nucleus_bridge = c(n ~ i, c ~ n),
    nucleus_collider = c(n ~ i + c),
    cell_bridge = c(c ~ i, n ~ c),
    cell_collider = c(c ~ i + n),
    iod_bridge = c(i ~ n, c ~ i),
    iod_collider = c(i ~ n + c)
  )
}

user_mechanism_model_set <- function() {
  define_model_set(
    iod_to_nucleus_to_cell = c(n ~ i, c ~ n),
    cell_to_nucleus_to_iod = c(n ~ c, i ~ n),
    nucleus_to_iod_and_cell = c(i ~ n, c ~ n)
  )
}

single_mechanism_model_set <- function(mechanism) {
  switch(
    mechanism,
    iod_to_nucleus_to_cell = define_model_set(
      iod_to_nucleus_to_cell = c(n ~ i, c ~ n)
    ),
    cell_to_nucleus_to_iod = define_model_set(
      cell_to_nucleus_to_iod = c(n ~ c, i ~ n)
    ),
    nucleus_to_iod_and_cell = define_model_set(
      nucleus_to_iod_and_cell = c(i ~ n, c ~ n)
    ),
    stop("Unknown user mechanism: ", mechanism)
  )
}

empty_failures <- function() {
  tibble(
    analysis_type = character(),
    replicate_id = character(),
    tree_id = character(),
    evolutionary_model = character(),
    generating_class = character(),
    residual_regime = character(),
    error = character()
  )
}

metadata_defaults <- function(metadata = list()) {
  defaults <- list(
    analysis_type = "unspecified",
    replicate_id = NA_character_,
    tree_id = "focal_primary",
    evolutionary_model = "lambda",
    generating_class = NA_character_,
    residual_regime = NA_character_
  )
  utils::modifyList(defaults, metadata)
}

fit_path <- function(
  source,
  tree,
  evolutionary_model = "lambda",
  model_set = exhaustive_model_set(),
  metadata = list(),
  retain_fit = FALSE
) {
  meta <- metadata_defaults(metadata)
  meta$evolutionary_model <- evolutionary_model
  runtime_warnings <- character()
  result <- tryCatch(
    {
      prepared <- prepare_model_data(source, tree)
      fit <- withCallingHandlers(
        phylo_path(
          model_set,
          data = prepared$data,
          tree = prepared$tree,
          model = evolutionary_model
        ),
        warning = function(warning) {
          runtime_warnings <<- c(runtime_warnings, conditionMessage(warning))
          invokeRestart("muffleWarning")
        }
      )
      ranking <- as_tibble(as.data.frame(summary(fit), stringsAsFactors = FALSE)) %>%
        mutate(
          analysis_type = meta$analysis_type,
          replicate_id = meta$replicate_id,
          tree_id = meta$tree_id,
          evolutionary_model = meta$evolutionary_model,
          generating_class = meta$generating_class,
          residual_regime = meta$residual_regime,
          n_species = nrow(prepared$data),
          rank = row_number(),
          globally_supported = is.finite(.data$p) & .data$p >= alpha,
          delta_competitive = is.finite(.data$delta_CICc) & .data$delta_CICc <= competitive_delta,
          admissible_competitive = .data$globally_supported & .data$delta_competitive,
          stored_warning_count = length(fit$warnings),
          runtime_warning_count = length(runtime_warnings)
        ) %>%
        select(
          analysis_type,
          replicate_id,
          tree_id,
          evolutionary_model,
          generating_class,
          residual_regime,
          n_species,
          rank,
          everything()
        )

      dsep_rows <- list()
      for (candidate in names(fit$d_sep)) {
        component <- as.data.frame(fit$d_sep[[candidate]], stringsAsFactors = FALSE)
        if (nrow(component)) {
          dsep_rows[[length(dsep_rows) + 1L]] <- tibble(
            analysis_type = meta$analysis_type,
            replicate_id = meta$replicate_id,
            tree_id = meta$tree_id,
            evolutionary_model = meta$evolutionary_model,
            candidate_model = candidate,
            independence_claim = component$d_sep,
            component_p = component$p,
            fitted_phylogenetic_parameter = component$phylo_par
          )
        }
      }
      warnings <- tibble(
        analysis_type = meta$analysis_type,
        replicate_id = meta$replicate_id,
        tree_id = meta$tree_id,
        evolutionary_model = meta$evolutionary_model,
        warning_source = c(
          rep("phylopath_stored", length(fit$warnings)),
          rep("runtime", length(runtime_warnings))
        ),
        warning = c(unlist(fit$warnings, use.names = FALSE), runtime_warnings)
      )
      list(
        ranking = ranking,
        dsep = bind_rows(dsep_rows),
        warnings = warnings,
        failure = empty_failures(),
        fit = if (retain_fit) fit else NULL,
        prepared = if (retain_fit) prepared else NULL
      )
    },
    error = function(error) {
      list(
        ranking = tibble(),
        dsep = tibble(),
        warnings = tibble(),
        failure = tibble(
          analysis_type = meta$analysis_type,
          replicate_id = meta$replicate_id,
          tree_id = meta$tree_id,
          evolutionary_model = meta$evolutionary_model,
          generating_class = meta$generating_class,
          residual_regime = meta$residual_regime,
          error = conditionMessage(error)
        ),
        fit = NULL,
        prepared = NULL
      )
    }
  )
  result
}

matrix_edges <- function(fitted_dag, mechanism) {
  coefficient <- fitted_dag$coef
  standard_error <- fitted_dag$se
  rows <- list()
  for (parent in rownames(coefficient)) {
    for (child in colnames(coefficient)) {
      estimate <- coefficient[parent, child]
      if (is.finite(estimate) && estimate != 0) {
        se <- standard_error[parent, child]
        rows[[length(rows) + 1L]] <- tibble(
          mechanism = mechanism,
          equivalence_class = "nucleus_bridge",
          parent = unname(node_labels[[parent]]),
          child = unname(node_labels[[child]]),
          standardized_coefficient = as.numeric(estimate),
          standard_error = as.numeric(se),
          ci_low = as.numeric(estimate - 1.96 * se),
          ci_high = as.numeric(estimate + 1.96 * se),
          interval_excludes_zero = (estimate - 1.96 * se) * (estimate + 1.96 * se) > 0
        )
      }
    }
  }
  bind_rows(rows)
}

fit_phylogenetic_regressions <- function(source, tree, evolutionary_model = "lambda") {
  prepared <- prepare_model_data(source, tree)
  specifications <- list(
    iod_to_nucleus = n ~ i,
    nucleus_to_iod = i ~ n,
    iod_to_cell = c ~ i,
    cell_to_iod = i ~ c,
    nucleus_to_cell = c ~ n,
    cell_to_nucleus = n ~ c,
    iod_and_nucleus_to_cell = c ~ i + n,
    iod_and_cell_to_nucleus = n ~ i + c,
    nucleus_and_cell_to_iod = i ~ n + c
  )
  rows <- list()
  for (relation_id in names(specifications)) {
    formula <- specifications[[relation_id]]
    fit <- suppressWarnings(phylolm(
      formula,
      data = prepared$data,
      phy = prepared$tree,
      model = evolutionary_model
    ))
    coefficient_table <- as.data.frame(summary(fit)$coefficients)
    coefficient_table$term <- rownames(coefficient_table)
    coefficient_table <- coefficient_table[coefficient_table$term != "(Intercept)", , drop = FALSE]
    response <- all.vars(formula)[[1]]
    for (row_index in seq_len(nrow(coefficient_table))) {
      row <- coefficient_table[row_index, , drop = FALSE]
      rows[[length(rows) + 1L]] <- tibble(
        relation_id = relation_id,
        evolutionary_model = evolutionary_model,
        outcome = unname(node_labels[[response]]),
        predictor = unname(node_labels[[row$term[[1]]]]),
        standardized_coefficient = as.numeric(row$Estimate),
        standard_error = as.numeric(row$StdErr),
        ci_low = as.numeric(row$Estimate - 1.96 * row$StdErr),
        ci_high = as.numeric(row$Estimate + 1.96 * row$StdErr),
        t_value = as.numeric(row$t.value),
        p_value = as.numeric(row$p.value),
        fitted_phylogenetic_parameter = if (length(fit$optpar)) as.numeric(fit$optpar) else NA_real_,
        r_squared = as.numeric(fit$r.squared),
        n_species = nrow(prepared$data)
      )
    }
  }
  bind_rows(rows)
}

raw_correlations <- function(source) {
  specifications <- list(
    iod_nucleus = c("relative_iod_index", "nucleus_area_um2"),
    iod_cell = c("relative_iod_index", "cell_area_um2"),
    nucleus_cell = c("nucleus_area_um2", "cell_area_um2")
  )
  rows <- list()
  for (comparison in names(specifications)) {
    columns <- specifications[[comparison]]
    left <- source[[columns[[1]]]]
    right <- source[[columns[[2]]]]
    spearman <- suppressWarnings(cor.test(left, right, method = "spearman", exact = FALSE))
    pearson_log <- cor.test(log10(left), log10(right), method = "pearson")
    rows[[length(rows) + 1L]] <- tibble(
      comparison = comparison,
      left_trait = columns[[1]],
      right_trait = columns[[2]],
      spearman_rho = unname(spearman$estimate),
      spearman_p = spearman$p.value,
      log10_pearson_r = unname(pearson_log$estimate),
      log10_pearson_p = pearson_log$p.value,
      n_species = length(left),
      phylogenetically_corrected = FALSE
    )
  }
  bind_rows(rows)
}

fit_many <- function(specifications, progress_label) {
  rankings <- list()
  failures <- list()
  for (index in seq_along(specifications)) {
    specification <- specifications[[index]]
    result <- fit_path(
      specification$source,
      specification$tree,
      evolutionary_model = specification$evolutionary_model,
      metadata = specification$metadata
    )
    if (nrow(result$ranking)) rankings[[length(rankings) + 1L]] <- result$ranking
    if (nrow(result$failure)) failures[[length(failures) + 1L]] <- result$failure
    if (index %% 25L == 0L || index == length(specifications)) {
      cat(sprintf("%s: %d/%d fits complete\n", progress_label, index, length(specifications)))
    }
  }
  list(ranking = bind_rows(rankings), failures = bind_rows(failures))
}

summarize_simulation <- function(ranking) {
  if (!nrow(ranking)) return(tibble())
  per_fit <- ranking %>%
    group_by(.data$generating_class, .data$residual_regime, .data$replicate_id) %>%
    summarize(
      top_model = .data$model[which.min(.data$rank)],
      top_global_fit_pass = .data$globally_supported[which.min(.data$rank)],
      true_delta = ifelse(
        .data$generating_class[[1]] %in% .data$model,
        .data$delta_CICc[match(.data$generating_class[[1]], .data$model)],
        NA_real_
      ),
      true_global_fit_pass = ifelse(
        .data$generating_class[[1]] %in% .data$model,
        .data$globally_supported[match(.data$generating_class[[1]], .data$model)],
        NA
      ),
      n_admissible_competitive = sum(.data$admissible_competitive),
      .groups = "drop"
    )
  per_fit %>%
    group_by(.data$generating_class, .data$residual_regime) %>%
    summarize(
      n_simulations = n(),
      true_class_top_rate = mean(.data$top_model == .data$generating_class),
      true_class_delta_le_2_rate = mean(.data$true_delta <= competitive_delta, na.rm = TRUE),
      true_class_global_fit_rate = mean(.data$true_global_fit_pass, na.rm = TRUE),
      unique_true_class_rate = mean(
        .data$top_model == .data$generating_class & .data$n_admissible_competitive == 1
      ),
      any_bridge_top_rate = mean(grepl("_bridge$", .data$top_model)),
      median_admissible_competitive_models = median(.data$n_admissible_competitive),
      .groups = "drop"
    )
}

phylogenetic_draw <- function(tree, regime) {
  if (regime == "iid") {
    values <- stats::rnorm(length(tree$tip.label))
    names(values) <- tree$tip.label
    return(zscore(values))
  }
  if (regime == "brownian") {
    values <- ape::rTraitCont(tree, model = "BM", sigma = 1)
    return(zscore(values[tree$tip.label]))
  }
  stop("Unknown simulation residual regime: ", regime)
}

simulate_traits <- function(tree, generating_class, residual_regime, coefficients) {
  first <- phylogenetic_draw(tree, residual_regime)
  second_error <- phylogenetic_draw(tree, residual_regime)
  third_error <- phylogenetic_draw(tree, residual_regime)
  combine <- function(parent, error, beta) {
    zscore(beta * parent + sqrt(max(1 - beta^2, 0.05)) * error)
  }
  if (generating_class == "independent") {
    iod <- first
    n <- second_error
    c <- third_error
  } else if (generating_class == "nucleus_bridge") {
    iod <- first
    n <- combine(iod, second_error, coefficients[["iod_to_nucleus"]])
    c <- combine(n, third_error, coefficients[["nucleus_to_cell"]])
  } else if (generating_class == "cell_bridge") {
    iod <- first
    c <- combine(iod, second_error, coefficients[["iod_to_cell"]])
    n <- combine(c, third_error, coefficients[["cell_to_nucleus"]])
  } else if (generating_class == "iod_bridge") {
    n <- first
    iod <- combine(n, second_error, coefficients[["nucleus_to_iod"]])
    c <- combine(iod, third_error, coefficients[["iod_to_cell"]])
  } else {
    stop("Unknown simulation generating class: ", generating_class)
  }
  tibble(
    species = paste0("D. ", tree$tip.label),
    tree_tip = tree$tip.label,
    relative_iod_index = 10^iod,
    nucleus_area_um2 = 10^n,
    cell_area_um2 = 10^c
  )
}

cat("Loading and validating frozen path-analysis inputs...\n")
traits <- read_csv(trait_path, show_col_types = FALSE)
measurement_draws <- read_csv(measurement_path, show_col_types = FALSE)
focal_tree <- read.tree(focal_tree_path)
published_main_tree <- read.tree(published_main_path)
published_bootstrap_trees <- read.nexus(published_bootstrap_path)
if (nrow(traits) < 3L) stop("Primary path panel must contain at least three species")
if (n_distinct(traits$species) != nrow(traits)) {
  stop("Primary path panel must contain one row per species")
}
validate_traits(traits)
if (length(published_bootstrap_trees) != 200L) {
  stop("Expected exactly 200 published bootstrap trees")
}

cat("Running primary exhaustive equivalence-class comparison...\n")
primary <- fit_path(
  traits,
  focal_tree,
  evolutionary_model = "lambda",
  metadata = list(analysis_type = "primary", replicate_id = "point_estimate"),
  retain_fit = TRUE
)
if (nrow(primary$failure)) stop("Primary path fit failed: ", primary$failure$error[[1]])
if (!identical(sort(primary$ranking$model), sort(model_order))) {
  stop("Primary ranking did not contain the exact ten testable equivalence classes")
}

standardized_traits <- primary$prepared$data %>%
  as.data.frame() %>%
  rownames_to_column("tree_tip") %>%
  rename(
    relative_iod_log10_z = i,
    nucleus_area_log10_z = n,
    cell_area_log10_z = c
  ) %>%
  left_join(select(traits, tree_tip, species), by = "tree_tip") %>%
  select(species, tree_tip, everything())

pairwise <- fit_phylogenetic_regressions(traits, focal_tree, "lambda")
correlations <- raw_correlations(traits)

cat("Estimating all three user-proposed orientations...\n")
mechanism_comparison <- fit_path(
  traits,
  focal_tree,
  evolutionary_model = "lambda",
  model_set = user_mechanism_model_set(),
  metadata = list(analysis_type = "user_mechanisms", replicate_id = "point_estimate")
)
if (nrow(mechanism_comparison$failure)) {
  stop("User-mechanism comparison failed: ", mechanism_comparison$failure$error[[1]])
}
if (max(mechanism_comparison$ranking$CICc) - min(mechanism_comparison$ranking$CICc) > 1e-8 ||
    max(mechanism_comparison$ranking$p) - min(mechanism_comparison$ranking$p) > 1e-8) {
  stop("Markov-equivalent user mechanisms unexpectedly received different fit scores")
}

mechanism_edges <- list()
mechanism_failures <- list()
for (mechanism in c(
  "iod_to_nucleus_to_cell",
  "cell_to_nucleus_to_iod",
  "nucleus_to_iod_and_cell"
)) {
  result <- fit_path(
    traits,
    focal_tree,
    evolutionary_model = "lambda",
    model_set = single_mechanism_model_set(mechanism),
    metadata = list(analysis_type = "user_mechanism_edges", replicate_id = mechanism),
    retain_fit = TRUE
  )
  if (nrow(result$failure)) {
    mechanism_failures[[length(mechanism_failures) + 1L]] <- result$failure
  } else {
    mechanism_edges[[length(mechanism_edges) + 1L]] <- matrix_edges(best(result$fit), mechanism)
  }
}
mechanism_edges <- bind_rows(mechanism_edges)

cat("Running evolutionary-residual-model sensitivity...\n")
evolutionary_specs <- lapply(
  c("lambda", "BM", "OUfixedRoot", "OUrandomRoot"),
  function(model) list(
    source = traits,
    tree = focal_tree,
    evolutionary_model = model,
    metadata = list(
      analysis_type = "evolutionary_model_sensitivity",
      replicate_id = model,
      tree_id = "focal_primary"
    )
  )
)
evolutionary <- fit_many(evolutionary_specs, "Evolutionary-model sensitivity")

cat("Running leave-one-species-out sensitivity...\n")
loo_specs <- lapply(sort(traits$tree_tip), function(omitted_tip) {
  list(
    source = filter(traits, .data$tree_tip != omitted_tip),
    tree = focal_tree,
    evolutionary_model = "lambda",
    metadata = list(
      analysis_type = "leave_one_species_out",
      replicate_id = omitted_tip,
      tree_id = "focal_primary_minus_one"
    )
  )
})
loo <- fit_many(loo_specs, "Leave-one-species-out")

cat("Running measurement-bootstrap sensitivity...\n")
measurement_ids <- sort(unique(measurement_draws$bootstrap_replicate))
if (is.finite(options$measurement_limit)) {
  measurement_ids <- head(measurement_ids, options$measurement_limit)
}
measurement_specs <- lapply(measurement_ids, function(replicate) {
  list(
    source = filter(measurement_draws, .data$bootstrap_replicate == replicate),
    tree = focal_tree,
    evolutionary_model = "lambda",
    metadata = list(
      analysis_type = "measurement_bootstrap",
      replicate_id = sprintf("measurement_%04d", replicate),
      tree_id = "focal_primary"
    )
  )
})
measurement <- fit_many(measurement_specs, "Measurement bootstrap")

cat("Running focal-versus-published tree/species-set sensitivity...\n")
published_main_tree <- normalize_tree_labels(published_main_tree)
published_tips <- published_main_tree$tip.label
published_traits <- filter(traits, .data$tree_tip %in% published_tips)
if (nrow(published_traits) != 18L) stop("Published tree overlap must contain exactly 18 species")
tree_specs <- list(
  list(
    source = traits,
    tree = focal_tree,
    evolutionary_model = "lambda",
    metadata = list(
      analysis_type = "tree_species_set_sensitivity",
      replicate_id = "focal_primary",
      tree_id = "focal_primary"
    )
  ),
  list(
    source = published_traits,
    tree = focal_tree,
    evolutionary_model = "lambda",
    metadata = list(
      analysis_type = "tree_species_set_sensitivity",
      replicate_id = "focal_tree_18",
      tree_id = "focal_tree_18"
    )
  ),
  list(
    source = published_traits,
    tree = published_main_tree,
    evolutionary_model = "lambda",
    metadata = list(
      analysis_type = "tree_species_set_sensitivity",
      replicate_id = "published_main_18",
      tree_id = "published_main_18"
    )
  )
)
tree_limit <- min(length(published_bootstrap_trees), options$tree_limit)
for (tree_index in seq_len(tree_limit)) {
  tree_specs[[length(tree_specs) + 1L]] <- list(
    source = published_traits,
    tree = published_bootstrap_trees[[tree_index]],
    evolutionary_model = "lambda",
    metadata = list(
      analysis_type = "published_bootstrap_tree_sensitivity",
      replicate_id = sprintf("published_bootstrap_%03d", tree_index),
      tree_id = sprintf("published_bootstrap_%03d", tree_index)
    )
  )
}
tree_sensitivity <- fit_many(tree_specs, "Tree/species-set sensitivity")

cat("Running actual-tree small-sample simulation calibration...\n")
coefficient_lookup <- pairwise %>%
  filter(.data$relation_id %in% c(
    "iod_to_nucleus",
    "nucleus_to_iod",
    "iod_to_cell",
    "nucleus_to_cell",
    "cell_to_nucleus"
  )) %>%
  group_by(.data$relation_id) %>%
  summarize(value = .data$standardized_coefficient[[1]], .groups = "drop")
coefficients <- setNames(coefficient_lookup$value, coefficient_lookup$relation_id)
required_coefficients <- c(
  "iod_to_nucleus", "nucleus_to_iod", "iod_to_cell",
  "nucleus_to_cell", "cell_to_nucleus"
)
if (!all(required_coefficients %in% names(coefficients))) {
  stop("Missing observed coefficients required for simulation")
}

set.seed(options$seed)
simulation_specs <- list()
simulation_classes <- c("independent", "nucleus_bridge", "cell_bridge", "iod_bridge")
simulation_tree <- prune_tree(focal_tree, traits$tree_tip)
for (regime in options$simulation_regimes) {
  for (generating_class in simulation_classes) {
    for (replicate in seq_len(options$simulation_replicates)) {
      simulated <- simulate_traits(
        simulation_tree,
        generating_class,
        regime,
        coefficients
      )
      simulation_specs[[length(simulation_specs) + 1L]] <- list(
        source = simulated,
        tree = simulation_tree,
        evolutionary_model = "lambda",
        metadata = list(
          analysis_type = "simulation_calibration",
          replicate_id = sprintf("%s_%s_%04d", regime, generating_class, replicate),
          tree_id = "focal_primary",
          generating_class = generating_class,
          residual_regime = regime
        )
      )
    }
  }
}
simulation <- fit_many(simulation_specs, "Simulation calibration")
simulation_summary <- summarize_simulation(simulation$ranking)

all_failures <- bind_rows(
  primary$failure,
  mechanism_comparison$failure,
  bind_rows(mechanism_failures),
  evolutionary$failures,
  loo$failures,
  measurement$failures,
  tree_sensitivity$failures,
  simulation$failures
)

write_csv(primary$ranking, file.path(output_dir, "primary_model_ranking.csv"))
write_csv(primary$dsep, file.path(output_dir, "primary_dsep_tests.csv"))
write_csv(primary$warnings, file.path(output_dir, "primary_fit_warnings.csv"))
write_csv(standardized_traits, file.path(output_dir, "primary_standardized_traits.csv"))
write_csv(pairwise, file.path(output_dir, "pairwise_phylogenetic_regressions.csv"))
write_csv(correlations, file.path(output_dir, "raw_trait_correlations.csv"))
write_csv(
  mechanism_comparison$ranking,
  file.path(output_dir, "user_mechanism_equivalence_ranking.csv")
)
write_csv(mechanism_edges, file.path(output_dir, "user_mechanism_edges.csv"))
write_csv(
  evolutionary$ranking,
  file.path(output_dir, "evolutionary_model_sensitivity_rankings.csv")
)
write_csv(loo$ranking, file.path(output_dir, "leave_one_species_out_rankings.csv.gz"))
write_csv(
  measurement$ranking,
  file.path(output_dir, "measurement_bootstrap_rankings.csv.gz")
)
write_csv(tree_sensitivity$ranking, file.path(output_dir, "tree_sensitivity_rankings.csv.gz"))
write_csv(simulation$ranking, file.path(output_dir, "simulation_rankings.csv.gz"))
write_csv(simulation_summary, file.path(output_dir, "simulation_summary.csv"))
write_csv(all_failures, file.path(output_dir, "analysis_failures.csv"))

manifest <- list(
  analysis_id = "cell_nucleus_iod_phylogenetic_path_analysis_v1",
  generated_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  run_scope = if (quick_mode) "quick_validation" else "full_release",
  primary_species = nrow(traits),
  published_tree_species = nrow(published_traits),
  tested_equivalence_classes = length(model_order),
  labeled_dags = 25,
  markov_equivalence_classes = 11,
  saturated_class_testable = FALSE,
  primary_evolutionary_model = "lambda",
  alpha = alpha,
  competitive_delta_cicc = competitive_delta,
  trait_transformation = "log10 then sample z-score within every fit",
  relative_iod_units = "relative image-derived nuclear IOD index",
  iod_calibration_note = paste(
    "A shared positive scale factor does not change standardized log-scale",
    "path results; image-derived IOD measurement uncertainty does."
  ),
  measurement_bootstrap_replicates = length(measurement_ids),
  leave_one_out_fits = length(loo_specs),
  evolutionary_models = c("lambda", "BM", "OUfixedRoot", "OUrandomRoot"),
  published_bootstrap_trees_analyzed = tree_limit,
  simulation_replicates_per_class_and_regime = options$simulation_replicates,
  simulation_generating_classes = simulation_classes,
  simulation_residual_regimes = options$simulation_regimes,
  simulation_seed = options$seed,
  phylopath_version = as.character(packageVersion("phylopath")),
  phylolm_version = as.character(packageVersion("phylolm")),
  ape_version = as.character(packageVersion("ape")),
  r_version = R.version.string,
  failure_count = nrow(all_failures),
  primary_lambda_boundary_warning_count = sum(
    primary$warnings$warning_source == "phylopath_stored"
  ),
  primary_top_model = primary$ranking$model[[1]],
  primary_second_model = primary$ranking$model[[2]],
  primary_top_two_delta_cicc = primary$ranking$delta_CICc[[2]],
  path_results_are_exploratory = TRUE,
  causal_direction_identified = FALSE,
  user_mechanisms_markov_equivalent = TRUE,
  common_user_mechanism_basis_claim = "relative_iod independent of cell_size conditional on nucleus_size"
)
write_json(
  manifest,
  file.path(output_dir, "analysis_manifest.json"),
  pretty = TRUE,
  auto_unbox = TRUE,
  digits = 12
)

cat("Completed exhaustive cell–nucleus–iod path analysis.\n")
cat("Primary top models:\n")
print(select(head(primary$ranking, 3), model, p, CICc, delta_CICc, w))
cat("Failures:", nrow(all_failures), "\n")
