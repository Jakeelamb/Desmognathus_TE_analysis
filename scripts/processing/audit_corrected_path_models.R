#!/usr/bin/env Rscript

# Corrected final-18 phylogenetic path-model sensitivity audit.
#
# This script never reads the historical genome_size_pg column. The focal
# image-derived node is a relative nuclear-IOD proxy and every output is marked
# exploratory. Historical path-analysis outputs are preserved.

prefer_active_conda_r_library <- function() {
  conda_prefix <- Sys.getenv("CONDA_PREFIX", "")
  if (nzchar(conda_prefix)) {
    conda_lib <- normalizePath(file.path(conda_prefix, "lib", "R", "library"), mustWork = FALSE)
    if (dir.exists(conda_lib)) .libPaths(conda_lib)
  }
}
prefer_active_conda_r_library()

suppressPackageStartupMessages({
  library(ape)
  library(dplyr)
  library(jsonlite)
  library(phylopath)
  library(readr)
  library(tibble)
})

args <- commandArgs(trailingOnly = TRUE)
phase <- "all"
if ("--phase" %in% args) {
  position <- match("--phase", args)
  if (position == length(args)) stop("--phase requires a value")
  phase <- args[[position + 1]]
}
valid_phases <- c("data", "trees", "loo", "all")
if (!phase %in% valid_phases) stop("Unknown phase: ", phase)

project_root <- normalizePath(getwd())
input_path <- file.path(
  project_root,
  "results/data/corrected/path_analysis/corrected_path_input_sensitivity_analysis18_v1.csv"
)
focal_tree_path <- file.path(
  project_root,
  "results/data/corrected/phylogeny/desmognathus_time_tree_analysis18_v1.nwk"
)
published_tree_path <- file.path(
  project_root,
  "results/data/corrected/phylogeny/stewart_wiens_2025_main_analysis18_v1.nwk"
)
bootstrap_tree_path <- file.path(
  project_root,
  "results/data/corrected/phylogeny/stewart_wiens_2025_bootstrap_analysis18_v1.nex"
)
output_dir <- file.path(project_root, "results/data/corrected/path_analysis")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

for (path in c(input_path, focal_tree_path, published_tree_path, bootstrap_tree_path)) {
  if (!file.exists(path)) stop("Required corrected path input missing: ", path)
}

alpha <- 0.05
competitive_delta <- 2
anchor_morphology <- "image_balanced_selected50"
anchor_iod <- "image_qc_pass"
families <- c("te_iod", "iod_morphology", "integrated", "terminal_internal_iod")

zscore <- function(x) {
  x <- as.numeric(x)
  value <- stats::sd(x, na.rm = TRUE)
  if (!is.finite(value) || value == 0) stop("Cannot standardize zero-variance variable")
  as.numeric((x - mean(x, na.rm = TRUE)) / value)
}

model_set_for_family <- function(family) {
  switch(
    family,
    te_iod = define_model_set(
      proxy_null = c(),
      ltr_only = c(iod ~ ltr),
      evenness_only = c(iod ~ even),
      additive_composition = c(iod ~ ltr + even),
      mediated_evenness = c(even ~ ltr, iod ~ even)
    ),
    iod_morphology = define_model_set(
      morphology_null = c(),
      iod_to_nucleus = c(ns ~ iod),
      iod_to_cell = c(cs ~ iod),
      iod_nucleus_cell_chain = c(ns ~ iod, cs ~ ns)
    ),
    integrated = define_model_set(
      integrated_null = c(),
      iod_morphology_only = c(ns ~ iod, cs ~ ns),
      te_evenness_chain = c(even ~ ltr, iod ~ even, ns ~ iod, cs ~ ns),
      te_additive_chain = c(iod ~ ltr + even, ns ~ iod, cs ~ ns),
      te_nucleus_bypass = c(iod ~ even, ns ~ iod + even, cs ~ ns),
      te_cell_bypass = c(iod ~ even, ns ~ iod, cs ~ ns + even)
    ),
    terminal_internal_iod = define_model_set(
      proxy_null = c(),
      terminal_internal_only = c(iod ~ ti),
      ltr_terminal_internal_chain = c(ti ~ ltr, iod ~ ti),
      terminal_internal_evenness_additive = c(iod ~ ti + even),
      te_terminal_internal_state = c(ti ~ ltr + even, iod ~ ti)
    ),
    stop("Unknown corrected path family: ", family)
  )
}

prepare_data <- function(source, family, morphology_estimator, iod_subset, omit_species = NA_character_) {
  selected <- source %>%
    filter(
      .data$morphology_estimator == !!morphology_estimator,
      .data$iod_subset == !!iod_subset
    )
  if (!is.na(omit_species)) selected <- selected %>% filter(.data$species != !!omit_species)

  transformed <- selected %>%
    mutate(
      iod = zscore(log(.data$relative_nuclear_iod_proxy)),
      ns = zscore(log(.data$nuc_area_um2)),
      cs = zscore(log(.data$cell_area_um2)),
      ltr = zscore(.data$ltr_line_logratio),
      even = zscore(.data$pielou_evenness),
      ti = ifelse(
        is.finite(.data$terminal_internal_ratio_median) & .data$terminal_internal_ratio_median > 0,
        log(.data$terminal_internal_ratio_median),
        NA_real_
      )
    )
  if (family == "terminal_internal_iod") {
    transformed <- transformed %>% filter(is.finite(.data$ti)) %>% mutate(ti = zscore(.data$ti))
  }

  columns <- switch(
    family,
    te_iod = c("species", "iod", "ltr", "even"),
    iod_morphology = c("species", "iod", "ns", "cs"),
    integrated = c("species", "iod", "ns", "cs", "ltr", "even"),
    terminal_internal_iod = c("species", "iod", "ti", "ltr", "even")
  )
  out <- transformed[, columns, drop = FALSE]
  out <- out[stats::complete.cases(out), , drop = FALSE]
  if (anyDuplicated(out$species)) stop("Duplicate species in corrected path input")
  if (nrow(out) < 8) stop("Fewer than 8 complete species for family ", family)
  out
}

prune_tree <- function(tree, species) {
  tree$tip.label <- sub("^D\\.\\s*", "", tree$tip.label)
  missing <- setdiff(species, tree$tip.label)
  if (length(missing)) stop("Tree missing path species: ", paste(missing, collapse = ", "))
  tree <- drop.tip(tree, setdiff(tree$tip.label, species))
  if (any(tree$edge.length <= 0) || any(!is.finite(tree$edge.length))) {
    stop("Path tree contains non-positive or non-finite branches")
  }
  tree
}

matrix_edges <- function(coef, se, fit_id, best_model, metadata) {
  rows <- list()
  for (child in colnames(coef)) {
    for (parent in rownames(coef)) {
      value <- coef[parent, child]
      if (is.finite(value) && value != 0) {
        standard_error <- se[parent, child]
        rows[[length(rows) + 1]] <- tibble(
          fit_id = fit_id,
          family = metadata$family,
          phase = metadata$phase,
          tree_id = metadata$tree_id,
          morphology_estimator = metadata$morphology_estimator,
          iod_subset = metadata$iod_subset,
          omitted_species = metadata$omitted_species,
          n_species = metadata$n_species,
          best_model = best_model,
          parent = parent,
          child = child,
          coefficient = as.numeric(value),
          std_error = as.numeric(standard_error),
          approx_ci_low = as.numeric(value - 1.96 * standard_error),
          approx_ci_high = as.numeric(value + 1.96 * standard_error),
          interval_excludes_zero = as.logical((value - 1.96 * standard_error) * (value + 1.96 * standard_error) > 0),
          absolute_genome_size_used = FALSE,
          publication_claim_allowed = FALSE
        )
      }
    }
  }
  bind_rows(rows)
}

release_gate <- function(ranking) {
  finite <- is.finite(ranking$CICc) & is.finite(ranking$delta_CICc) & is.finite(ranking$w)
  supported <- finite & is.finite(ranking$p) & ranking$p >= alpha
  weights_valid <- all(finite) && abs(sum(ranking$w) - 1) <= 1e-6
  competitive <- supported & ranking$delta_CICc <= competitive_delta
  if (!all(finite) || !weights_valid) {
    "blocked_nonfinite_ranking"
  } else if (!any(supported)) {
    "blocked_no_globally_supported_model"
  } else if (!supported[[1]]) {
    "blocked_top_model_rejected"
  } else if (sum(competitive) > 1) {
    "supported_competitive_model_set"
  } else {
    "supported_unique_top_model"
  }
}

empty_result <- function() {
  list(ranking = tibble(), edges = tibble(), dsep = tibble(), failures = tibble())
}

fit_once <- function(
  source,
  tree,
  family,
  phase_name,
  tree_id,
  morphology_estimator,
  iod_subset,
  omitted_species = NA_character_
) {
  fit_id <- paste(
    phase_name, tree_id, family, morphology_estimator, iod_subset,
    ifelse(is.na(omitted_species), "none", omitted_species),
    sep = "__"
  )
  result <- empty_result()
  attempt <- tryCatch(
    {
      data <- prepare_data(source, family, morphology_estimator, iod_subset, omitted_species)
      current_tree <- prune_tree(tree, data$species)
      data <- as.data.frame(data)
      rownames(data) <- data$species
      data$species <- NULL
      data <- data[current_tree$tip.label, , drop = FALSE]
      model_set <- model_set_for_family(family)
      fit <- suppressWarnings(phylo_path(model_set, data = data, tree = current_tree))
      ranking <- as.data.frame(summary(fit), stringsAsFactors = FALSE)
      ranking$model <- as.character(ranking$model)
      gate <- release_gate(ranking)
      top_supported <- is.finite(ranking$p[[1]]) && ranking$p[[1]] >= alpha
      metadata <- list(
        family = family,
        phase = phase_name,
        tree_id = tree_id,
        morphology_estimator = morphology_estimator,
        iod_subset = iod_subset,
        omitted_species = omitted_species,
        n_species = nrow(data)
      )
      result$ranking <- as_tibble(ranking) %>%
        mutate(
          fit_id = fit_id,
          family = family,
          phase = phase_name,
          tree_id = tree_id,
          morphology_estimator = morphology_estimator,
          iod_subset = iod_subset,
          omitted_species = omitted_species,
          n_species = nrow(data),
          rank = row_number(),
          global_fit_pass = is.finite(.data$p) & .data$p >= alpha,
          finite_ranking_row = is.finite(.data$CICc) & is.finite(.data$delta_CICc) & is.finite(.data$w),
          release_gate_status = gate,
          winner_claim_allowed_by_ranking_gate = gate == "supported_unique_top_model",
          absolute_genome_size_used = FALSE,
          publication_claim_allowed = FALSE,
          warning_count = length(fit$warnings)
        ) %>%
        select(
          fit_id, family, phase, tree_id, morphology_estimator, iod_subset,
          omitted_species, n_species, rank, everything()
        )

      best_model <- ranking$model[[1]]
      best_dag <- tryCatch(best(fit), error = function(e) NULL)
      if (!is.null(best_dag)) {
        result$edges <- matrix_edges(best_dag$coef, best_dag$se, fit_id, best_model, metadata)
      }

      dsep_rows <- list()
      for (model_name in names(fit$d_sep)) {
        component <- as.data.frame(fit$d_sep[[model_name]])
        if (nrow(component)) {
          dsep_rows[[length(dsep_rows) + 1]] <- as_tibble(component) %>%
            transmute(
              fit_id = fit_id,
              family = family,
              phase = phase_name,
              tree_id = tree_id,
              morphology_estimator = morphology_estimator,
              iod_subset = iod_subset,
              omitted_species = omitted_species,
              n_species = nrow(data),
              candidate_model = model_name,
              independence_claim = .data$d_sep,
              component_p = .data$p,
              fitted_phylogenetic_parameter = .data$phylo_par,
              component_fit_pass = is.finite(.data$p) & .data$p >= alpha
            )
        }
      }
      result$dsep <- bind_rows(dsep_rows)
      result
    },
    error = function(e) {
      result$failures <- tibble(
        fit_id = fit_id,
        family = family,
        phase = phase_name,
        tree_id = tree_id,
        morphology_estimator = morphology_estimator,
        iod_subset = iod_subset,
        omitted_species = omitted_species,
        error = conditionMessage(e)
      )
      result
    }
  )
  attempt
}

combine_results <- function(results) {
  list(
    ranking = bind_rows(lapply(results, `[[`, "ranking")),
    edges = bind_rows(lapply(results, `[[`, "edges")),
    dsep = bind_rows(lapply(results, `[[`, "dsep")),
    failures = bind_rows(lapply(results, `[[`, "failures"))
  )
}

write_phase <- function(name, result) {
  paths <- list(
    ranking = file.path(output_dir, paste0("corrected_path_", name, "_rankings_analysis18_v1.csv")),
    edges = file.path(output_dir, paste0("corrected_path_", name, "_best_edges_analysis18_v1.csv")),
    dsep = file.path(output_dir, paste0("corrected_path_", name, "_basis_sets_analysis18_v1.csv")),
    failures = file.path(output_dir, paste0("corrected_path_", name, "_failures_analysis18_v1.csv"))
  )
  write_csv(result$ranking, paths$ranking)
  write_csv(result$edges, paths$edges)
  write_csv(result$dsep, paths$dsep)
  write_csv(result$failures, paths$failures)
  paths
}

run_data_phase <- function(source, focal_tree, published_tree) {
  results <- list()
  tree_specs <- list(focal = focal_tree, published_main = published_tree)
  morphology_values <- sort(unique(source$morphology_estimator))
  iod_values <- sort(unique(source$iod_subset))
  for (tree_id in names(tree_specs)) {
    for (iod_subset in iod_values) {
      for (family in c("te_iod", "terminal_internal_iod")) {
        results[[length(results) + 1]] <- fit_once(
          source, tree_specs[[tree_id]], family, "data", tree_id,
          anchor_morphology, iod_subset
        )
      }
      for (morphology_estimator in morphology_values) {
        for (family in c("iod_morphology", "integrated")) {
          results[[length(results) + 1]] <- fit_once(
            source, tree_specs[[tree_id]], family, "data", tree_id,
            morphology_estimator, iod_subset
          )
        }
      }
    }
  }
  combine_results(results)
}

run_tree_phase <- function(source, published_tree, bootstrap_trees) {
  results <- list()
  tree_specs <- c(list(published_main = published_tree), unclass(bootstrap_trees))
  names(tree_specs)[-1] <- sprintf("bootstrap_%03d", seq_along(bootstrap_trees))
  for (tree_id in names(tree_specs)) {
    for (family in families) {
      results[[length(results) + 1]] <- fit_once(
        source, tree_specs[[tree_id]], family, "trees", tree_id,
        anchor_morphology, anchor_iod
      )
    }
  }
  combine_results(results)
}

run_loo_phase <- function(source, published_tree) {
  results <- list()
  for (family in families) {
    base <- prepare_data(source, family, anchor_morphology, anchor_iod)
    for (omitted in sort(base$species)) {
      results[[length(results) + 1]] <- fit_once(
        source, published_tree, family, "loo", "published_main",
        anchor_morphology, anchor_iod, omitted_species = omitted
      )
    }
  }
  combine_results(results)
}

if (Sys.getenv("DESMO_PATH_LIBRARY_ONLY", "0") != "1") {
  source <- read_csv(input_path, show_col_types = FALSE)
  focal_tree <- read.tree(focal_tree_path)
  published_tree <- read.tree(published_tree_path)
  bootstrap_trees <- read.nexus(bootstrap_tree_path)
  if (length(bootstrap_trees) != 200) stop("Expected 200 published bootstrap time trees")

  phase_outputs <- list()
  if (phase %in% c("data", "all")) {
    cat("Running corrected path data-specification phase...\n")
    phase_outputs$data <- write_phase("data_sensitivity", run_data_phase(source, focal_tree, published_tree))
  }
  if (phase %in% c("trees", "all")) {
    cat("Running corrected path 200-tree phase...\n")
    phase_outputs$trees <- write_phase("tree_sensitivity", run_tree_phase(source, published_tree, bootstrap_trees))
  }
  if (phase %in% c("loo", "all")) {
    cat("Running corrected path leave-one-species-out phase...\n")
    phase_outputs$loo <- write_phase("leave_one_out", run_loo_phase(source, published_tree))
  }

  manifest <- list(
    analysis_id = "corrected_path_models_analysis18_v1",
    phase = phase,
    phylopath_version = as.character(packageVersion("phylopath")),
    alpha = alpha,
    competitive_delta_cicc = competitive_delta,
    anchor_morphology = anchor_morphology,
    anchor_iod_subset = anchor_iod,
    n_published_bootstrap_trees = length(bootstrap_trees),
    families = families,
    absolute_genome_size_used = FALSE,
    publication_claim_allowed = FALSE,
    historical_outputs_overwritten = FALSE,
    outputs = lapply(phase_outputs, function(group) lapply(group, function(path) sub(paste0(project_root, "/"), "", path)))
  )
  write_json(
    manifest,
    file.path(output_dir, paste0("corrected_path_models_", phase, "_analysis18_v1.manifest.json")),
    pretty = TRUE,
    auto_unbox = TRUE
  )
  cat("Completed corrected path phase:", phase, "\n")
}
