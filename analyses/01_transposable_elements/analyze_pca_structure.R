#!/usr/bin/env Rscript

# Rebuild and formally audit structure in the released TE34 superfamily CLR-PCA.

suppressPackageStartupMessages({
  library(ape)
  library(cluster)
  library(dplyr)
  library(jsonlite)
  library(purrr)
  library(readr)
  library(tidyr)
  library(vegan)
})

args <- commandArgs(trailingOnly = TRUE)
write_outputs <- "--write" %in% args
unknown_args <- setdiff(args, "--write")
if (length(unknown_args) > 0L) {
  stop("Unknown argument(s): ", paste(unknown_args, collapse = ", "))
}

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_arg[[1L]]))
analysis_dir <- dirname(script_path)
project_root <- normalizePath(file.path(analysis_dir, "..", ".."))
data_dir <- file.path(analysis_dir, "data")

base_seed <- 20260810L
permutations <- 99999L
gap_replicates <- 1000L
subsample_replicates <- 1000L
k_values <- 1:8
kmeans_nstart <- 1000L
stability_nstart <- 1000L
gap_nstart <- 200L
minimum_cluster_size <- 3L
downstream_minimum_cluster_size <- 5L
minimum_mean_silhouette <- 0.25
minimum_stability_median_ari <- 0.80
minimum_stability_q10_ari <- 0.60
minimum_dimension_ari <- 0.80

paths <- list(
  te34_panel = file.path(project_root, "data", "identity", "te34_panel.csv"),
  order_composition = file.path(data_dir, "te_order_composition.csv"),
  superfamily_composition = file.path(data_dir, "te_superfamily_composition.csv"),
  scores = file.path(data_dir, "te_pca_scores.csv"),
  variance = file.path(data_dir, "te_pca_variance.csv"),
  loadings = file.path(data_dir, "te_pca_loadings.csv"),
  clr = file.path(data_dir, "te_pca_clr_matrix.csv"),
  prevalence = file.path(data_dir, "te_feature_prevalence.csv"),
  traits = file.path(
    project_root,
    "analyses",
    "03_phylogenetic_path",
    "data",
    "organismal_traits.csv"
  ),
  tree = file.path(
    project_root,
    "analyses",
    "03_phylogenetic_path",
    "trees",
    "source_time_tree_46.tre"
  ),
  pca_recompute_script = file.path(analysis_dir, "recompute_pca.py"),
  analysis_script = script_path,
  environment = file.path(project_root, "environment.yml"),
  python_project = file.path(project_root, "pyproject.toml"),
  python_lock = file.path(project_root, "uv.lock")
)

missing_inputs <- unlist(paths)[!file.exists(unlist(paths))]
if (length(missing_inputs) > 0L) {
  stop("Missing manifest input(s): ", paste(missing_inputs, collapse = ", "))
}

output_paths <- list(
  diagnostics = file.path(data_dir, "te_pca_clustering_diagnostics.csv"),
  stability = file.path(data_dir, "te_pca_clustering_stability.csv"),
  assignments = file.path(data_dir, "te_pca_candidate_cluster_assignments.csv"),
  traits = file.path(data_dir, "te_pca_trait_association_tests.csv"),
  phylogeny = file.path(data_dir, "te_pca_phylogenetic_signal_tests.csv"),
  manifest = file.path(data_dir, "te_pca_clustering_analysis_manifest.json")
)

sha256 <- function(path) {
  line <- system2("sha256sum", path, stdout = TRUE)
  strsplit(line[[1L]], "[[:space:]]+")[[1L]][[1L]]
}

format_species <- function(values) {
  values <- sort(unique(as.character(values)))
  if (length(values) == 0L) "none" else paste(values, collapse = ";")
}

format_group_sizes <- function(groups) {
  counts <- sort(table(as.character(groups)))
  paste(paste(names(counts), as.integer(counts), sep = "="), collapse = ";")
}

fit_kmeans <- function(x, k, seed, nstart = kmeans_nstart) {
  set.seed(seed)
  stats::kmeans(
    x,
    centers = k,
    nstart = nstart,
    iter.max = 1000L,
    algorithm = "Hartigan-Wong"
  )
}

canonicalize_clusters <- function(labels, scores) {
  label_values <- sort(unique(labels))
  centroids <- lapply(label_values, function(label) {
    members <- labels == label
    c(mean(scores[members, 1L]), mean(scores[members, 2L]))
  })
  centroid_matrix <- do.call(rbind, centroids)
  ordering <- order(centroid_matrix[, 1L], centroid_matrix[, 2L], label_values)
  mapping <- setNames(seq_along(ordering), label_values[ordering])
  unname(mapping[as.character(labels)])
}

adjusted_rand_index <- function(left, right) {
  contingency <- table(left, right)
  choose_two <- function(values) sum(values * (values - 1) / 2)
  observed <- choose_two(contingency)
  row_pairs <- choose_two(rowSums(contingency))
  column_pairs <- choose_two(colSums(contingency))
  total_pairs <- choose_two(sum(contingency))
  expected <- row_pairs * column_pairs / total_pairs
  maximum <- 0.5 * (row_pairs + column_pairs)
  if (isTRUE(all.equal(maximum, expected))) return(1)
  (observed - expected) / (maximum - expected)
}

silhouette_statistics <- function(labels, distances) {
  values <- cluster::silhouette(labels, distances)[, "sil_width"]
  c(mean = mean(values), median = median(values), minimum = min(values))
}

calinski_harabasz <- function(fit, n, k) {
  if (k == 1L) return(NA_real_)
  between <- fit$totss - fit$tot.withinss
  (between / (k - 1L)) / (fit$tot.withinss / (n - k))
}

davies_bouldin <- function(x, fit) {
  k <- nrow(fit$centers)
  if (k == 1L) return(NA_real_)
  scatter <- vapply(seq_len(k), function(cluster_id) {
    members <- x[fit$cluster == cluster_id, , drop = FALSE]
    centered <- sweep(members, 2L, fit$centers[cluster_id, ], "-")
    mean(sqrt(rowSums(centered^2)))
  }, numeric(1L))
  center_distances <- as.matrix(dist(fit$centers))
  ratios <- outer(scatter, scatter, "+") / center_distances
  diag(ratios) <- -Inf
  mean(apply(ratios, 1L, max))
}

summarize_ari <- function(values) {
  tibble(
    n_replicates = length(values),
    ari_mean = mean(values),
    ari_sd = if (length(values) > 1L) sd(values) else 0,
    ari_minimum = min(values),
    ari_q10 = unname(quantile(values, 0.10, names = FALSE, type = 7)),
    ari_median = median(values),
    ari_q90 = unname(quantile(values, 0.90, names = FALSE, type = 7)),
    ari_maximum = max(values),
    fraction_ari_at_least_0_80 = mean(values >= 0.80)
  )
}

permutation_p <- function(observed, randomized, alternative = c("greater", "less")) {
  alternative <- match.arg(alternative)
  exceedances <- if (alternative == "greater") {
    sum(randomized >= observed - sqrt(.Machine$double.eps))
  } else {
    sum(randomized <= observed + sqrt(.Machine$double.eps))
  }
  (exceedances + 1) / (length(randomized) + 1)
}

scores <- readr::read_csv(paths$scores, show_col_types = FALSE)
variance <- readr::read_csv(paths$variance, show_col_types = FALSE)
clr_long <- readr::read_csv(paths$clr, show_col_types = FALSE) |>
  filter(te_level == "superfamily")
prevalence <- readr::read_csv(paths$prevalence, show_col_types = FALSE)
traits <- readr::read_csv(paths$traits, show_col_types = FALSE)

pc_columns <- grep("^PC[0-9]+$", names(scores), value = TRUE)
if (nrow(scores) != 34L || length(pc_columns) != 23L) {
  stop("Released PCA must contain 34 species and 23 nonzero axes")
}
if (anyDuplicated(scores$species)) stop("PCA species are not unique")

clr_wide <- clr_long |>
  select(species, feature, clr_value) |>
  pivot_wider(names_from = feature, values_from = clr_value) |>
  arrange(match(species, scores$species))
if (!identical(clr_wide$species, scores$species)) stop("CLR/PCA species order mismatch")
clr_features <- prevalence |>
  filter(te_level == "superfamily", pca_included) |>
  pull(feature)
if (!setequal(clr_features, setdiff(names(clr_wide), "species"))) {
  stop("CLR features do not equal the prevalence-declared PCA features")
}
clr_matrix <- as.matrix(clr_wide[, clr_features])
rownames(clr_matrix) <- clr_wide$species
all_pc_matrix <- as.matrix(scores[, pc_columns])
rownames(all_pc_matrix) <- scores$species

if (ncol(clr_matrix) != 24L || qr(clr_matrix)$rank != 23L) {
  stop("Released superfamily CLR matrix must contain 24 features with rank 23")
}
distance_difference <- max(abs(as.matrix(dist(clr_matrix)) - as.matrix(dist(all_pc_matrix))))
if (distance_difference > 1e-10) stop("All released PCs do not preserve CLR distances")

representations <- list(
  all_23_pcs = all_pc_matrix,
  first_6_pcs = all_pc_matrix[, 1:6, drop = FALSE],
  first_2_pcs = all_pc_matrix[, 1:2, drop = FALSE]
)
representation_pc_counts <- c(all_23_pcs = 23L, first_6_pcs = 6L, first_2_pcs = 2L)
representation_variance <- c(
  all_23_pcs = 1,
  first_6_pcs = variance$cumulative_variance_explained[[6L]],
  first_2_pcs = variance$cumulative_variance_explained[[2L]]
)

fits <- list()
diagnostic_rows <- list()
for (representation_index in seq_along(representations)) {
  representation_name <- names(representations)[[representation_index]]
  x <- representations[[representation_name]]
  distances <- dist(x)
  for (k in k_values) {
    seed <- if (representation_index == 1L) {
      base_seed + k
    } else {
      base_seed + representation_index * 100L + k
    }
    fit <- fit_kmeans(x, k, seed)
    fits[[paste(representation_name, k, sep = "::")]] <- fit
    sizes <- sort(as.integer(table(fit$cluster)))
    silhouette <- if (k == 1L) {
      c(mean = NA_real_, median = NA_real_, minimum = NA_real_)
    } else {
      silhouette_statistics(fit$cluster, distances)
    }
    diagnostic_rows[[length(diagnostic_rows) + 1L]] <- tibble(
      representation = representation_name,
      n_pcs = representation_pc_counts[[representation_name]],
      cumulative_variance_explained = representation_variance[[representation_name]],
      k = k,
      kmeans_seed = seed,
      kmeans_nstart = kmeans_nstart,
      within_cluster_sum_squares = fit$tot.withinss,
      mean_silhouette = silhouette[["mean"]],
      median_silhouette = silhouette[["median"]],
      minimum_silhouette = silhouette[["minimum"]],
      calinski_harabasz = calinski_harabasz(fit, nrow(x), k),
      davies_bouldin = davies_bouldin(x, fit),
      minimum_cluster_size = min(sizes),
      singleton_count = sum(sizes == 1L),
      cluster_sizes_ascending = paste(sizes, collapse = ";")
    )
  }
}
diagnostics <- bind_rows(diagnostic_rows)

gap_kmeans <- function(x, k) {
  stats::kmeans(
    x,
    centers = k,
    nstart = gap_nstart,
    iter.max = 1000L,
    algorithm = "Hartigan-Wong"
  )
}

run_gap <- function(d_power, seed) {
  set.seed(seed)
  result <- cluster::clusGap(
    all_pc_matrix,
    FUNcluster = gap_kmeans,
    K.max = max(k_values),
    B = gap_replicates,
    d.power = d_power,
    spaceH0 = "scaledPCA"
  )
  table <- as.data.frame(result$Tab)
  tibble(
    k = seq_len(nrow(table)),
    gap = table$gap,
    gap_standard_error = table$SE.sim,
    distance_power = d_power,
    seed = seed
  )
}

gap_d1 <- run_gap(1, base_seed + 1L)
gap_d2 <- run_gap(2, base_seed)

gap_selection <- function(gap_table, method) {
  cluster::maxSE(
    gap_table$gap,
    gap_table$gap_standard_error,
    method = method,
    SE.factor = 1
  )
}

gap_methods <- c("firstSEmax", "Tibs2001SEmax", "globalSEmax")
gap_selected <- list(
  d1 = setNames(vapply(gap_methods, function(method) {
    gap_selection(gap_d1, method)
  }, integer(1L)), gap_methods),
  d2 = setNames(vapply(gap_methods, function(method) {
    gap_selection(gap_d2, method)
  }, integer(1L)), gap_methods)
)

diagnostics <- diagnostics |>
  left_join(
    gap_d1 |> select(k, gap_d1 = gap, gap_d1_se = gap_standard_error),
    by = "k"
  ) |>
  left_join(
    gap_d2 |> select(k, gap_d2 = gap, gap_d2_se = gap_standard_error),
    by = "k"
  ) |>
  mutate(
    gap_d1 = if_else(representation == "all_23_pcs", gap_d1, NA_real_),
    gap_d1_se = if_else(representation == "all_23_pcs", gap_d1_se, NA_real_),
    gap_d2 = if_else(representation == "all_23_pcs", gap_d2, NA_real_),
    gap_d2_se = if_else(representation == "all_23_pcs", gap_d2_se, NA_real_)
  )

full_diagnostics <- diagnostics |> filter(representation == "all_23_pcs")
unconstrained_best_k <- full_diagnostics |>
  filter(k >= 2L) |>
  slice_max(mean_silhouette, n = 1L, with_ties = FALSE) |>
  pull(k)
singleton_free_best_k <- full_diagnostics |>
  filter(k >= 2L, singleton_count == 0L) |>
  slice_max(mean_silhouette, n = 1L, with_ties = FALSE) |>
  pull(k)

# Stability is only evaluated for the strongest singleton-free candidate.
# All other k already fail the gap, silhouette, or minimum-size gates.
stability_k_values <- singleton_free_best_k
reference_labels <- lapply(stability_k_values, function(k) {
  fits[[paste("all_23_pcs", k, sep = "::")]]$cluster
})
names(reference_labels) <- as.character(stability_k_values)

stability_rows <- list()

# Exhaustive leave-one-superfamily-out subcomposition stability.
for (k in stability_k_values) {
  ari_values <- numeric(length(clr_features))
  for (feature_index in seq_along(clr_features)) {
    retained <- setdiff(seq_along(clr_features), feature_index)
    subcomposition_clr <- clr_matrix[, retained, drop = FALSE]
    subcomposition_clr <- subcomposition_clr - rowMeans(subcomposition_clr)
    fit <- fit_kmeans(
      subcomposition_clr,
      k,
      if (k == 2L) {
        20261000L + feature_index
      } else {
        20261000L + k * 100L + feature_index
      },
      nstart = stability_nstart
    )
    ari_values[[feature_index]] <- adjusted_rand_index(reference_labels[[as.character(k)]], fit$cluster)
  }
  stability_rows[[length(stability_rows) + 1L]] <- summarize_ari(ari_values) |>
    mutate(
      k = k,
      stability_type = "leave_one_superfamily_out",
      perturbation = "24 exhaustive subcompositions",
      seed = NA_integer_,
      .before = 1L
    )
}

# Refit each k to 80% species subsamples; compare labels only on sampled species.
set.seed(20260820L)
subsample_size <- ceiling(0.80 * nrow(all_pc_matrix))
species_subsamples <- replicate(
  subsample_replicates,
  sample.int(nrow(all_pc_matrix), subsample_size, replace = FALSE),
  simplify = FALSE
)
for (k in stability_k_values) {
  ari_values <- numeric(subsample_replicates)
  for (replicate_index in seq_len(subsample_replicates)) {
    retained <- species_subsamples[[replicate_index]]
    fit <- fit_kmeans(
      all_pc_matrix[retained, , drop = FALSE],
      k,
      20262000L + (k - 2L) * 2000L + replicate_index,
      nstart = stability_nstart
    )
    ari_values[[replicate_index]] <- adjusted_rand_index(
      reference_labels[[as.character(k)]][retained],
      fit$cluster
    )
  }
  stability_rows[[length(stability_rows) + 1L]] <- summarize_ari(ari_values) |>
    mutate(
      k = k,
      stability_type = "species_subsample_28_of_34",
      perturbation = paste0(subsample_size, " of 34 species without replacement"),
      seed = 20260820L,
      .before = 1L
    )
}

# Fixed-k dimensionality sensitivities.
for (representation_name in c("first_6_pcs", "first_2_pcs")) {
  for (k in stability_k_values) {
    alternate <- fits[[paste(representation_name, k, sep = "::")]]$cluster
    ari <- adjusted_rand_index(reference_labels[[as.character(k)]], alternate)
    stability_rows[[length(stability_rows) + 1L]] <- summarize_ari(ari) |>
      mutate(
        k = k,
        stability_type = "dimension_sensitivity",
        perturbation = representation_name,
        seed = diagnostics |>
          filter(representation == representation_name, .data$k == .env$k) |>
          pull(kmeans_seed),
        .before = 1L
      )
  }
}

stability <- bind_rows(stability_rows) |>
  select(
    k,
    stability_type,
    perturbation,
    seed,
    n_replicates,
    ari_mean,
    ari_sd,
    ari_minimum,
    ari_q10,
    ari_median,
    ari_q90,
    ari_maximum,
    fraction_ari_at_least_0_80
  ) |>
  arrange(k, stability_type, perturbation)

stability_gate <- function(k, type) {
  row <- stability |>
    filter(.data$k == .env$k, stability_type == type)
  if (nrow(row) != 1L) return(FALSE)
  row$ari_median >= minimum_stability_median_ari &&
    row$ari_q10 >= minimum_stability_q10_ari
}

dimension_gate <- function(k) {
  row <- stability |>
    filter(
      .data$k == .env$k,
      stability_type == "dimension_sensitivity",
      perturbation == "first_6_pcs"
    )
  if (nrow(row) != 1L) return(FALSE)
  row$ari_median >= minimum_dimension_ari
}

gap_supported_k <- unname(gap_selected$d2[["Tibs2001SEmax"]])
diagnostics <- diagnostics |>
  mutate(
    unconstrained_silhouette_best = representation == "all_23_pcs" & k == unconstrained_best_k,
    singleton_free_silhouette_best = representation == "all_23_pcs" & k == singleton_free_best_k,
    passes_minimum_cluster_size = minimum_cluster_size >= !!minimum_cluster_size,
    passes_downstream_cluster_size = minimum_cluster_size >= !!downstream_minimum_cluster_size,
    passes_mean_silhouette = !is.na(mean_silhouette) & mean_silhouette >= !!minimum_mean_silhouette,
    gap_supported = representation == "all_23_pcs" & k %in% gap_supported_k & k > 1L,
    passes_feature_stability = if_else(
      representation == "all_23_pcs" & k > 1L,
      vapply(k, function(value) stability_gate(value, "leave_one_superfamily_out"), logical(1L)),
      FALSE
    ),
    passes_species_stability = if_else(
      representation == "all_23_pcs" & k > 1L,
      vapply(k, function(value) stability_gate(value, "species_subsample_28_of_34"), logical(1L)),
      FALSE
    ),
    passes_six_pc_sensitivity = if_else(
      representation == "all_23_pcs" & k > 1L,
      vapply(k, dimension_gate, logical(1L)),
      FALSE
    ),
    accepted_cluster_solution = representation == "all_23_pcs" & k > 1L &
      passes_minimum_cluster_size & passes_downstream_cluster_size &
      passes_mean_silhouette & gap_supported & passes_feature_stability &
      passes_species_stability & passes_six_pc_sensitivity
  ) |>
  arrange(factor(representation, levels = names(representations)), k)

accepted_k_values <- diagnostics |>
  filter(accepted_cluster_solution) |>
  pull(k)
selected_k <- if (length(accepted_k_values) == 0L) {
  NULL
} else {
  diagnostics |>
    filter(accepted_cluster_solution) |>
    slice_max(mean_silhouette, n = 1L, with_ties = FALSE) |>
    pull(k)
}
analysis_status <- if (is.null(selected_k)) {
  "exploratory_post_hoc_no_discrete_cluster_accepted"
} else {
  "exploratory_post_hoc_discrete_cluster_candidate_accepted"
}

trait_columns <- c(
  "species",
  "reproductive_strategy",
  "max_svl_mm",
  "aquaticity_index"
)
joined_traits <- scores |>
  select(species, PC1, PC2) |>
  left_join(traits |> select(all_of(trait_columns)), by = "species")

rejection_reasons_for <- function(k) {
  row <- diagnostics |>
    filter(representation == "all_23_pcs", .data$k == .env$k)
  has_feature_stability <- any(
    stability$k == k & stability$stability_type == "leave_one_superfamily_out"
  )
  has_species_stability <- any(
    stability$k == k & stability$stability_type == "species_subsample_28_of_34"
  )
  has_dimension_stability <- any(
    stability$k == k & stability$stability_type == "dimension_sensitivity" &
      stability$perturbation == "first_6_pcs"
  )
  reasons <- c(
    if (!row$gap_supported) "gap_not_supported" else NULL,
    if (!row$passes_mean_silhouette) "mean_silhouette_below_0.25" else NULL,
    if (!row$passes_minimum_cluster_size) "minimum_cluster_size_below_3" else NULL,
    if (!row$passes_downstream_cluster_size) "minimum_cluster_size_below_5" else NULL,
    if (!has_feature_stability) {
      "feature_stability_not_evaluated_after_prior_gate_failure"
    } else if (!row$passes_feature_stability) {
      "feature_instability"
    } else {
      NULL
    },
    if (!has_species_stability) {
      "species_stability_not_evaluated_after_prior_gate_failure"
    } else if (!row$passes_species_stability) {
      "species_subsample_instability"
    } else {
      NULL
    },
    if (!has_dimension_stability) {
      "six_pc_sensitivity_not_evaluated_after_prior_gate_failure"
    } else if (!row$passes_six_pc_sensitivity) {
      "six_pc_dimension_instability"
    } else {
      NULL
    }
  )
  if (length(reasons) == 0L) "none" else paste(reasons, collapse = ";")
}

singleton_candidate_id <- paste0("singleton_free_k", singleton_free_best_k)
unconstrained_candidate_id <- paste0(
  "unconstrained_silhouette_k",
  unconstrained_best_k
)
reported_candidate_k_values <- unique(c(
  singleton_free_best_k,
  unconstrained_best_k,
  if (is.null(selected_k)) integer(0L) else selected_k
))
reported_candidate_ids <- vapply(reported_candidate_k_values, function(k) {
  if (k == singleton_free_best_k) {
    singleton_candidate_id
  } else if (k == unconstrained_best_k) {
    unconstrained_candidate_id
  } else {
    paste0("accepted_k", k)
  }
}, character(1L))
candidate_status <- if (is.null(selected_k)) {
  rep("rejected", length(reported_candidate_k_values))
} else {
  ifelse(reported_candidate_k_values == selected_k, "accepted", "rejected")
}
candidate_specs <- tibble(
  candidate_id = reported_candidate_ids,
  k = reported_candidate_k_values,
  candidate_status = candidate_status,
  rejection_reasons = vapply(
    reported_candidate_k_values,
    rejection_reasons_for,
    character(1L)
  )
)
singleton_candidate_status <- candidate_specs |>
  filter(candidate_id == singleton_candidate_id) |>
  pull(candidate_status)

assignment_rows <- lapply(seq_len(nrow(candidate_specs)), function(index) {
  spec <- candidate_specs[index, ]
  fit <- fits[[paste("all_23_pcs", spec$k, sep = "::")]]
  canonical <- canonicalize_clusters(fit$cluster, all_pc_matrix[, 1:2, drop = FALSE])
  cluster_sizes <- table(canonical)
  diagnostic <- diagnostics |>
    filter(representation == "all_23_pcs", k == spec$k)
  joined_traits |>
    mutate(
      candidate_id = spec$candidate_id,
      representation = "all_23_pcs",
      k = spec$k,
      cluster = paste0("C", canonical),
      cluster_size = as.integer(cluster_sizes[as.character(canonical)]),
      mean_silhouette = diagnostic$mean_silhouette,
      candidate_status = spec$candidate_status,
      rejection_reasons = spec$rejection_reasons,
      .after = species
    )
})
assignments <- bind_rows(assignment_rows) |>
  arrange(candidate_id, cluster, species)
singleton_cluster_sizes <- assignments |>
  filter(candidate_id == singleton_candidate_id) |>
  distinct(cluster, cluster_size) |>
  pull(cluster_size) |>
  sort()
singleton_candidate_phrase <- paste(
  singleton_candidate_status,
  paste0("k=", singleton_free_best_k),
  "candidate"
)

make_test_row <- function(
  test_id,
  analysis_role,
  test_family,
  predictor,
  n_species,
  excluded_species,
  group_sizes,
  statistic_name,
  statistic,
  r_squared = NA_real_,
  permutation_count = permutations,
  seed,
  p_value,
  multiplicity_family = "none",
  adjusted_method = "none",
  adjusted_p_value = p_value,
  interpretation
) {
  tibble(
    test_id = test_id,
    analysis_role = analysis_role,
    test_family = test_family,
    response = "24-superfamily CLR/Aitchison composition",
    predictor = predictor,
    n_species = n_species,
    excluded_species = excluded_species,
    group_sizes = group_sizes,
    statistic_name = statistic_name,
    statistic = statistic,
    variance_explained_r2 = r_squared,
    permutations = permutation_count,
    seed = seed,
    p_value = p_value,
    multiplicity_family = multiplicity_family,
    adjusted_method = adjusted_method,
    adjusted_p_value = adjusted_p_value,
    interpretation_boundary = interpretation
  )
}

run_factor_tests <- function(groups, test_prefix, role, seed, interpretation) {
  valid <- !is.na(groups) & groups != ""
  x <- clr_matrix[valid, , drop = FALSE]
  group <- droplevels(factor(groups[valid]))
  distances <- dist(x)
  metadata <- data.frame(group = group)
  set.seed(seed)
  permanova <- vegan::adonis2(distances ~ group, data = metadata, permutations = permutations)
  set.seed(seed + 1L)
  dispersion <- vegan::betadisper(
    distances,
    group,
    type = "median",
    bias.adjust = TRUE
  )
  dispersion_test <- vegan::permutest(dispersion, permutations = permutations)
  excluded <- rownames(clr_matrix)[!valid]
  list(
    permanova = make_test_row(
      paste0(test_prefix, "_permanova"),
      role,
      "PERMANOVA",
      test_prefix,
      sum(valid),
      format_species(excluded),
      format_group_sizes(group),
      "pseudo_F",
      permanova$F[[1L]],
      permanova$R2[[1L]],
      seed = seed,
      p_value = permanova$`Pr(>F)`[[1L]],
      interpretation = interpretation
    ),
    dispersion = make_test_row(
      paste0(test_prefix, "_permdisp"),
      "dispersion_diagnostic",
      "PERMDISP_spatial_median_bias_adjusted",
      test_prefix,
      sum(valid),
      format_species(excluded),
      format_group_sizes(group),
      "F",
      dispersion_test$tab$F[[1L]],
      seed = seed + 1L,
      p_value = dispersion_test$tab$`Pr(>F)`[[1L]],
      interpretation = "Dispersion diagnostic only; not evidence for or against a location effect."
    )
  )
}

species_index <- match(rownames(clr_matrix), joined_traits$species)
reproductive_strategy <- joined_traits$reproductive_strategy[species_index]
max_svl_mm <- joined_traits$max_svl_mm[species_index]
aquaticity <- joined_traits$aquaticity_index[species_index]

trait_test_rows <- list()
primary_tests <- run_factor_tests(
  reproductive_strategy,
  "reproductive_strategy",
  "primary_post_hoc_exploratory",
  base_seed + 10000L,
  paste(
    "Species-label permutation ignores phylogenetic non-independence;",
    "the expert categories were supplied after review of the initial PCA."
  )
)
trait_test_rows <- c(trait_test_rows, primary_tests)

direct_vs_larval <- ifelse(
  is.na(reproductive_strategy) | reproductive_strategy == "",
  NA_character_,
  ifelse(reproductive_strategy == "direct_development", "direct_development", "aquatic_larvae")
)
direct_tests <- run_factor_tests(
  direct_vs_larval,
  "direct_development_vs_aquatic_larvae",
  "planned_post_hoc_contrast",
  base_seed + 10010L,
  "Post-hoc contrast; ordinary species-label permutation is not phylogenetically corrected."
)
trait_test_rows <- c(trait_test_rows, direct_tests)

egg_environment <- ifelse(
  reproductive_strategy == "aquatic_eggs",
  "aquatic_eggs",
  ifelse(
    reproductive_strategy == "terrestrial_eggs_aquatic_larvae",
    "terrestrial_eggs",
    NA_character_
  )
)
egg_tests <- run_factor_tests(
  egg_environment,
  "egg_environment_among_aquatic_larvae",
  "planned_post_hoc_contrast",
  base_seed + 10020L,
  "Post-hoc contrast excludes direct developers; ordinary species-label permutation is not phylogenetically corrected."
)
trait_test_rows <- c(trait_test_rows, egg_tests)

valid_svl <- is.finite(max_svl_mm) & max_svl_mm > 0
svl_distance <- dist(clr_matrix[valid_svl, , drop = FALSE])
svl_metadata <- data.frame(log10_max_svl_mm = log10(max_svl_mm[valid_svl]))
set.seed(base_seed + 10030L)
svl_fit <- vegan::adonis2(
  svl_distance ~ log10_max_svl_mm,
  data = svl_metadata,
  permutations = permutations
)
trait_test_rows[[length(trait_test_rows) + 1L]] <- make_test_row(
  "max_svl_log10_permanova",
  "size_sensitivity",
  "PERMANOVA",
  "log10(max_svl_mm)",
  sum(valid_svl),
  format_species(rownames(clr_matrix)[!valid_svl]),
  "continuous",
  "pseudo_F",
  svl_fit$F[[1L]],
  svl_fit$R2[[1L]],
  seed = base_seed + 10030L,
  p_value = svl_fit$`Pr(>F)`[[1L]],
  interpretation = "Species-level size association; ordinary permutations are not phylogenetically corrected."
)

valid_adjusted <- valid_svl & !is.na(reproductive_strategy) & reproductive_strategy != ""
adjusted_distance <- dist(clr_matrix[valid_adjusted, , drop = FALSE])
adjusted_metadata <- data.frame(
  log10_max_svl_mm = log10(max_svl_mm[valid_adjusted]),
  reproductive_strategy = droplevels(factor(reproductive_strategy[valid_adjusted]))
)
set.seed(base_seed + 10040L)
adjusted_fit <- vegan::adonis2(
  adjusted_distance ~ log10_max_svl_mm + reproductive_strategy,
  data = adjusted_metadata,
  permutations = permutations,
  by = "margin"
)
for (term in c("log10_max_svl_mm", "reproductive_strategy")) {
  row <- adjusted_fit[term, ]
  trait_test_rows[[length(trait_test_rows) + 1L]] <- make_test_row(
    paste0("size_adjusted_model_", term),
    "size_adjusted_sensitivity",
    "marginal_PERMANOVA",
    paste0(term, " adjusted for the other model term"),
    sum(valid_adjusted),
    format_species(rownames(clr_matrix)[!valid_adjusted]),
    if (term == "reproductive_strategy") {
      format_group_sizes(adjusted_metadata$reproductive_strategy)
    } else {
      "continuous"
    },
    "pseudo_F",
    row$F[[1L]],
    row$R2[[1L]],
    seed = base_seed + 10040L,
    p_value = row$`Pr(>F)`[[1L]],
    interpretation = paste(
      "Marginal species-level association after the other recorded term;",
      "not a phylogenetically corrected or causal decomposition."
    )
  )
}

aquaticity_tests <- run_factor_tests(
  as.character(aquaticity),
  "adult_aquaticity_factor",
  "adult_lifestyle_sensitivity",
  base_seed + 10050L,
  paste(
    "Adult aquaticity is distinct from egg environment; factor coding does not assume equal ordinal spacing.",
    "Ordinary species-label permutation is not phylogenetically corrected."
  )
)
trait_test_rows <- c(trait_test_rows, aquaticity_tests)

# Conditional association of the rejected singleton-free k=2 candidate with expert strategy.
k2_labels <- assignments |>
  filter(candidate_id == singleton_candidate_id) |>
  arrange(match(species, rownames(clr_matrix))) |>
  pull(cluster)
valid_ari <- !is.na(reproductive_strategy) & reproductive_strategy != ""
observed_ari <- adjusted_rand_index(k2_labels[valid_ari], reproductive_strategy[valid_ari])
set.seed(base_seed + 10060L)
random_ari <- replicate(
  permutations,
  adjusted_rand_index(k2_labels[valid_ari], sample(reproductive_strategy[valid_ari]))
)
contingency <- table(k2_labels[valid_ari], reproductive_strategy[valid_ari])
contingency_frame <- as.data.frame(contingency, stringsAsFactors = FALSE)
contingency_text <- paste0(
  contingency_frame$Var1,
  ":",
  contingency_frame$Var2,
  "=",
  contingency_frame$Freq,
  collapse = ";"
)
trait_test_rows[[length(trait_test_rows) + 1L]] <- make_test_row(
  paste0(singleton_candidate_status, "_k", singleton_free_best_k, "_reproductive_strategy_ari"),
  paste0(singleton_candidate_status, "_cluster_candidate_sensitivity"),
  "fixed_partition_label_permutation",
  "reproductive_strategy",
  sum(valid_ari),
  format_species(rownames(clr_matrix)[!valid_ari]),
  contingency_text,
  "adjusted_Rand_index",
  observed_ari,
  seed = base_seed + 10060L,
  p_value = permutation_p(observed_ari, random_ari, "greater"),
  interpretation = paste(
    "Conditional on the", singleton_candidate_phrase,
    "; this conditional test is not independent validation of clusters and is not phylogenetically corrected."
  )
)

trait_tests <- bind_rows(trait_test_rows)
contrast_location <- trait_tests$analysis_role == "planned_post_hoc_contrast" &
  trait_tests$test_family == "PERMANOVA"
trait_tests$adjusted_p_value[contrast_location] <- p.adjust(
  trait_tests$p_value[contrast_location],
  method = "holm"
)
trait_tests$multiplicity_family[contrast_location] <- "two_strategy_contrasts"
trait_tests$adjusted_method[contrast_location] <- "Holm"
contrast_dispersion <- trait_tests$analysis_role == "dispersion_diagnostic" &
  grepl("direct_development|egg_environment", trait_tests$test_id)
trait_tests$adjusted_p_value[contrast_dispersion] <- p.adjust(
  trait_tests$p_value[contrast_dispersion],
  method = "holm"
)
trait_tests$multiplicity_family[contrast_dispersion] <- "two_contrast_dispersion_diagnostics"
trait_tests$adjusted_method[contrast_dispersion] <- "Holm"

# Prune the unreproducible collaborator source tree directly to the exact TE34 tips.
source_tree <- ape::read.tree(paths$tree)
missing_tree_species <- setdiff(rownames(clr_matrix), source_tree$tip.label)
if (length(missing_tree_species) > 0L) {
  stop("Source tree is missing TE34 species: ", paste(missing_tree_species, collapse = ", "))
}
te_tree <- ape::drop.tip(source_tree, setdiff(source_tree$tip.label, rownames(clr_matrix)))
if (length(te_tree$tip.label) != 34L || any(te_tree$edge.length <= 0)) {
  stop("Pruned TE34 tree must have 34 tips and positive edge lengths")
}

k_mult_components <- function(y, covariance_inverse, covariance, one, q_value) {
  phylogenetic_mean <- drop(crossprod(one, covariance_inverse %*% y)) / q_value
  residuals <- sweep(y, 2L, phylogenetic_mean, "-")
  expected_ratio <- (sum(diag(covariance)) - nrow(y) / q_value) / (nrow(y) - 1L)
  (sum(residuals * residuals) / sum(residuals * (covariance_inverse %*% residuals))) /
    expected_ratio
}

tree_order <- te_tree$tip.label
tree_indices <- match(tree_order, rownames(clr_matrix))
covariance <- ape::vcv.phylo(te_tree)[tree_order, tree_order]
covariance_inverse <- solve(covariance)
one <- matrix(1, nrow = length(tree_order), ncol = 1L)
q_value <- drop(crossprod(one, covariance_inverse %*% one))
k_representations <- list(
  full_24_feature_clr = clr_matrix[tree_indices, , drop = FALSE],
  first_6_pcs = all_pc_matrix[tree_indices, 1:6, drop = FALSE],
  first_2_pcs = all_pc_matrix[tree_indices, 1:2, drop = FALSE]
)
observed_k <- vapply(
  k_representations,
  k_mult_components,
  numeric(1L),
  covariance_inverse = covariance_inverse,
  covariance = covariance,
  one = one,
  q_value = q_value
)

set.seed(base_seed)
k_null <- matrix(
  NA_real_,
  nrow = permutations,
  ncol = length(k_representations),
  dimnames = list(NULL, names(k_representations))
)
for (replicate_index in seq_len(permutations)) {
  permutation <- sample.int(length(tree_order))
  for (representation_index in seq_along(k_representations)) {
    randomized <- k_representations[[representation_index]][permutation, , drop = FALSE]
    k_null[replicate_index, representation_index] <- k_mult_components(
      randomized,
      covariance_inverse,
      covariance,
      one,
      q_value
    )
  }
}

make_phylogeny_row <- function(
  test_id,
  role,
  representation,
  statistic_name,
  statistic,
  randomized,
  alternative,
  seed,
  branch_length_use,
  details,
  interpretation
) {
  tibble(
    test_id = test_id,
    analysis_role = role,
    representation = representation,
    response = "TE composition",
    n_species = 34L,
    tree_path = "analyses/03_phylogenetic_path/trees/source_time_tree_46.tre",
    statistic_name = statistic_name,
    statistic = statistic,
    null_mean = mean(randomized),
    null_sd = sd(randomized),
    null_q025 = unname(quantile(randomized, 0.025, names = FALSE, type = 8)),
    null_median = median(randomized),
    null_q975 = unname(quantile(randomized, 0.975, names = FALSE, type = 8)),
    alternative = alternative,
    permutations = length(randomized),
    seed = seed,
    p_value = permutation_p(statistic, randomized, alternative),
    branch_length_use = branch_length_use,
    details = details,
    interpretation_boundary = interpretation
  )
}

phylogeny_rows <- lapply(seq_along(k_representations), function(index) {
  representation <- names(k_representations)[[index]]
  make_phylogeny_row(
    paste0("k_mult_", representation),
    if (index == 1L) "primary_phylogenetic_signal" else "dimension_sensitivity",
    representation,
    "K_mult",
    observed_k[[index]],
    k_null[, index],
    "greater",
    base_seed,
    "dated_branch_lengths",
    "Adams 2014 Equation 7; complete multivariate species vectors permuted among fixed tips.",
    paste(
      "Tree topology, calibration, and citation provenance are unresolved;",
      "K_mult below 1 is weaker than Brownian expectation even when above the randomized null."
    )
  )
})

patristic <- ape::cophenetic.phylo(te_tree)[tree_order, tree_order]
upper <- upper.tri(patristic)

cluster_distance_statistic <- function(labels) {
  same <- outer(labels, labels, "==")[upper]
  distances <- patristic[upper]
  within <- mean(distances[same])
  between <- mean(distances[!same])
  list(
    statistic = between - within,
    within = within,
    between = between,
    cluster_means = vapply(sort(unique(labels)), function(label) {
      members <- labels == label
      block <- patristic[members, members, drop = FALSE]
      mean(block[upper.tri(block)])
    }, numeric(1L))
  )
}

tree_labels_for_candidate <- function(representation_name, k) {
  labels <- fits[[paste(representation_name, k, sep = "::")]]$cluster
  labels[match(tree_order, rownames(all_pc_matrix))]
}

k2_tree_labels <- tree_labels_for_candidate("all_23_pcs", singleton_free_best_k)
k2_distance <- cluster_distance_statistic(k2_tree_labels)
set.seed(base_seed + 2L)
k2_distance_null <- replicate(
  permutations,
  cluster_distance_statistic(sample(k2_tree_labels))$statistic
)
phylogeny_rows[[length(phylogeny_rows) + 1L]] <- make_phylogeny_row(
  paste0(singleton_candidate_status, "_k", singleton_free_best_k, "_patristic_separation"),
  paste0(singleton_candidate_status, "_cluster_candidate_sensitivity"),
  paste0("all_23_pcs_k", singleton_free_best_k),
  "mean_between_minus_mean_within_patristic_distance",
  k2_distance$statistic,
  k2_distance_null,
  "greater",
  base_seed + 2L,
  "dated_branch_lengths",
  paste0(
    "pooled_within=", signif(k2_distance$within, 10L),
    ";between=", signif(k2_distance$between, 10L),
    ";cluster_specific_within=", paste(signif(k2_distance$cluster_means, 10L), collapse = ";")
  ),
  paste(
    "Conditional on the", singleton_candidate_phrase, "; neither group is monophyletic,",
    "and the pooled within distance is dominated by the larger group of",
    max(singleton_cluster_sizes), "species."
  )
)

fitch_tree <- ape::reorder.phylo(te_tree, order = "postorder")
fitch_n_tip <- length(fitch_tree$tip.label)
fitch_n_node <- fitch_tree$Nnode
fitch_parents <- unique(fitch_tree$edge[, 1L])
fitch_children <- split(fitch_tree$edge[, 2L], fitch_tree$edge[, 1L])
fitch_root <- setdiff(fitch_tree$edge[, 1L], fitch_tree$edge[, 2L])[[1L]]

fitch_steps <- function(tip_states) {
  n_tip <- fitch_n_tip
  n_node <- fitch_n_node
  state_levels <- sort(unique(tip_states))
  n_states <- length(state_levels)
  costs <- matrix(Inf, nrow = n_tip + n_node, ncol = n_states)
  tip_values <- as.integer(factor(tip_states, levels = state_levels))
  costs[cbind(seq_len(n_tip), tip_values)] <- 0
  for (node in fitch_parents) {
    children <- fitch_children[[as.character(node)]]
    for (state in seq_len(n_states)) {
      costs[node, state] <- sum(vapply(children, function(child) {
        transition_cost <- as.integer(seq_len(n_states) != state)
        min(costs[child, ] + transition_cost)
      }, numeric(1L)))
    }
  }
  min(costs[fitch_root, ])
}

observed_fitch <- fitch_steps(k2_tree_labels)
set.seed(base_seed + 3L)
fitch_null <- replicate(permutations, fitch_steps(sample(k2_tree_labels)))
phylogeny_rows[[length(phylogeny_rows) + 1L]] <- make_phylogeny_row(
  paste0(singleton_candidate_status, "_k", singleton_free_best_k, "_fitch_transitions"),
  paste0(singleton_candidate_status, "_cluster_candidate_topology_sensitivity"),
  paste0("all_23_pcs_k", singleton_free_best_k),
  "minimum_label_transitions",
  observed_fitch,
  fitch_null,
  "less",
  base_seed + 3L,
  "topology_only",
  paste0(
    "Sankoff dynamic program with equal transition costs; exact label counts ",
    paste(singleton_cluster_sizes, collapse = "/"),
    " retained."
  ),
  paste(
    "Conditional on the", singleton_candidate_phrase,
    "; few transitions mean labels occupy few tree regions,",
    "not that either candidate group is monophyletic."
  )
)

first2_labels <- tree_labels_for_candidate("first_2_pcs", singleton_free_best_k)
first2_distance <- cluster_distance_statistic(first2_labels)
set.seed(base_seed + 4L)
first2_distance_null <- replicate(
  permutations,
  cluster_distance_statistic(sample(first2_labels))$statistic
)
phylogeny_rows[[length(phylogeny_rows) + 1L]] <- make_phylogeny_row(
  paste0("first_2_pc_k", singleton_free_best_k, "_patristic_separation"),
  "dimension_and_rejected_cluster_sensitivity",
  paste0("first_2_pcs_k", singleton_free_best_k),
  "mean_between_minus_mean_within_patristic_distance",
  first2_distance$statistic,
  first2_distance_null,
  "greater",
  base_seed + 4L,
  "dated_branch_lengths",
  paste0(
    "pooled_within=", signif(first2_distance$within, 10L),
    ";between=", signif(first2_distance$between, 10L),
    ";ARI_vs_all_23_PCs=",
    stability |>
      filter(k == singleton_free_best_k, stability_type == "dimension_sensitivity", perturbation == "first_2_pcs") |>
      pull(ari_median) |>
      signif(10L)
  ),
  "The visually cleaner two-PC partition is not the full-composition partition and is not phylogenetically compact."
)

phylogeny_tests <- bind_rows(phylogeny_rows)

root_to_tip <- ape::node.depth.edgelength(te_tree)[seq_along(te_tree$tip.label)]
manifest <- list(
  analysis_id = "te34_final_clr_pca_structure_audit_v1",
  analysis_date = "2026-08-10",
  analysis_status = analysis_status,
  selected_k = if (is.null(selected_k)) NA_integer_ else selected_k,
  inputs = lapply(paths, function(path) {
    list(
      path = sub(paste0(project_root, "/"), "", path, fixed = TRUE),
      sha256 = sha256(path)
    )
  }),
  pca_geometry = list(
    species = nrow(clr_matrix),
    ubiquitous_superfamilies = ncol(clr_matrix),
    clr_rank = qr(clr_matrix)$rank,
    nonzero_pcs = length(pc_columns),
    maximum_clr_to_all_pc_distance_difference = distance_difference,
    pc1_pc2_variance = representation_variance[["first_2_pcs"]],
    pc1_pc6_variance = representation_variance[["first_6_pcs"]],
    zero_replacement = FALSE,
    sparse_features_excluded_from_primary_pca = prevalence |>
      filter(te_level == "superfamily", !pca_included) |>
      select(feature, n_species_positive, prevalence_fraction) |>
      as.data.frame()
  ),
  clustering = list(
    primary_representation = "all 23 unscaled released PC scores; exactly preserves CLR/Aitchison distances",
    k_grid = k_values,
    kmeans_nstart = kmeans_nstart,
    gap_replicates = gap_replicates,
    gap_reference = "scaledPCA uniform reference",
    gap_primary_rule = list(
      distance_power = 2L,
      selector = "Tibs2001SEmax",
      standard_error_multiplier = 1,
      selected_k = gap_selected$d2[["Tibs2001SEmax"]]
    ),
    gap_selector_sensitivities = list(
      distance_power_1 = as.list(gap_selected$d1),
      distance_power_2 = as.list(gap_selected$d2)
    ),
    stability = list(
      leave_one_superfamily_out = "recompute the CLR of each 23-part subcomposition",
      species_subsampling = paste(subsample_replicates, "replicates of", subsample_size, "of 34 species"),
      sensitivity_representations = c("first_6_pcs", "first_2_pcs")
    ),
    conservative_release_gates = list(
      note = "Audit gates adopted for this rebuild; not preregistered or confirmatory.",
      requires_primary_gap_rule_to_select_more_than_one_cluster = TRUE,
      minimum_cluster_size = minimum_cluster_size,
      downstream_minimum_cluster_size = downstream_minimum_cluster_size,
      minimum_mean_silhouette = minimum_mean_silhouette,
      minimum_feature_stability_median_ari = minimum_stability_median_ari,
      minimum_feature_stability_q10_ari = minimum_stability_q10_ari,
      minimum_species_stability_median_ari = minimum_stability_median_ari,
      minimum_species_stability_q10_ari = minimum_stability_q10_ari,
      minimum_first_6_pc_ari = minimum_dimension_ari
    ),
    unconstrained_silhouette_best_k = unconstrained_best_k,
    singleton_free_silhouette_best_k = singleton_free_best_k,
    accepted_k = if (is.null(selected_k)) NA_integer_ else selected_k,
    conclusion = if (is.null(selected_k)) {
      "No k-means partition passes all release gates; candidate labels are audit-only."
    } else {
      paste0("k=", selected_k, " passes the rebuild release gates.")
    }
  ),
  trait_association = list(
    permutations = permutations,
    primary_endpoint = "reproductive_strategy PERMANOVA on full Aitchison distances",
    primary_status = "post hoc exploratory; expert categories followed review of initial PCA",
    dispersion_diagnostic = "PERMDISP using spatial medians with small-sample bias adjustment",
    planned_post_hoc_contrasts = c(
      "direct development versus aquatic larvae",
      "aquatic versus terrestrial eggs among aquatic-larval species"
    ),
    contrast_adjustment = "Holm within the two contrast tests",
    size_sensitivity = "marginal PERMANOVA of strategy and log10 maximum SVL",
    deliberately_not_duplicated = c(
      "development_mode is a deterministic collapse of reproductive_strategy on common complete cases",
      "microhabitat_class is the same partition as aquaticity_index"
    ),
    phylogenetic_correction = FALSE
  ),
  phylogenetic_signal = list(
    permutations = permutations,
    primary_test = "Adams 2014 Equation 7 K_mult",
    primary_k_mult = observed_k[["full_24_feature_clr"]],
    clr_all_pc_k_mult_absolute_difference = abs(
      observed_k[["full_24_feature_clr"]] -
        k_mult_components(
          all_pc_matrix[tree_indices, , drop = FALSE],
          covariance_inverse,
          covariance,
          one,
          q_value
        )
    ),
    source_tree_tips = length(source_tree$tip.label),
    retained_te34_tips = length(te_tree$tip.label),
    maximum_root_to_tip_spread = max(root_to_tip) - min(root_to_tip),
    all_edges_positive = all(te_tree$edge.length > 0),
    provenance_status = "provisional: source-tree citation, calibration, and branch-length provenance unresolved"
  ),
  random_seeds = list(
    base = base_seed,
    gap_d1 = base_seed + 1L,
    gap_d2 = base_seed,
    species_subsampling = 20260820L,
    trait_tests = paste0(base_seed + 10000L, " through ", base_seed + 10060L),
    k_mult = base_seed,
    candidate_patristic = base_seed + 2L,
    candidate_fitch = base_seed + 3L,
    first_2_pc_patristic = base_seed + 4L
  ),
  software = list(
    R = R.version.string,
    ape = as.character(packageVersion("ape")),
    cluster = as.character(packageVersion("cluster")),
    vegan = as.character(packageVersion("vegan"))
  ),
  interpretation_boundaries = c(
    "Continuous composition can have phylogenetic or trait association without forming stable discrete clusters.",
    if (is.null(selected_k)) {
      "No cluster labels are accepted for publication as biological groups."
    } else {
      "A release-gate partition remains exploratory and does not itself establish biological groups."
    },
    "Trait-label permutations are species-level and do not correct for phylogenetic non-independence.",
    "The reproductive-strategy hypothesis is post hoc rather than confirmatory.",
    "Tree-based results remain provisional until the source tree citation and calibration are supplied."
  )
)

frames <- list(
  diagnostics = diagnostics,
  stability = stability,
  assignments = assignments,
  traits = trait_tests,
  phylogeny = phylogeny_tests
)

if (write_outputs) {
  walk2(frames, output_paths[names(frames)], function(frame, path) {
    readr::write_csv(frame, path, na = "")
  })
  jsonlite::write_json(
    manifest,
    output_paths$manifest,
    auto_unbox = TRUE,
    pretty = TRUE,
    na = "null",
    digits = 16
  )
  cat(
    "WROTE: final PCA structure audit; selected_k=",
    if (is.null(selected_k)) "none" else selected_k,
    "\n",
    sep = ""
  )
} else {
  missing_outputs <- unlist(output_paths)[!file.exists(unlist(output_paths))]
  if (length(missing_outputs) > 0L) {
    stop("Missing canonical outputs; rerun with --write: ", paste(missing_outputs, collapse = ", "))
  }
  for (name in names(frames)) {
    expected <- readr::read_csv(output_paths[[name]], show_col_types = FALSE)
    comparison <- all.equal(
      as.data.frame(frames[[name]]),
      as.data.frame(expected),
      tolerance = 1e-10,
      check.attributes = FALSE
    )
    if (!isTRUE(comparison)) stop(name, " output mismatch: ", comparison)
  }
  expected_manifest_text <- paste(readLines(output_paths$manifest, warn = FALSE), collapse = "\n")
  generated_manifest_text <- as.character(jsonlite::toJSON(
    manifest,
    auto_unbox = TRUE,
    pretty = TRUE,
    na = "null",
    digits = 16
  ))
  if (!identical(trimws(expected_manifest_text), trimws(generated_manifest_text))) {
    stop("Manifest content mismatch")
  }
  cat(
    "PASS: final PCA structure audit reproduces canonical results; selected_k=",
    if (is.null(selected_k)) "none" else selected_k,
    "\n",
    sep = ""
  )
}

primary_trait <- trait_tests |> filter(test_id == "reproductive_strategy_permanova")
cat(
  sprintf(
    "K_mult=%.6f (p=%.5g); reproductive-strategy PERMANOVA R2=%.4f (p=%.5g)\n",
    observed_k[["full_24_feature_clr"]],
    phylogeny_tests |> filter(test_id == "k_mult_full_24_feature_clr") |> pull(p_value),
    primary_trait$variance_explained_r2,
    primary_trait$p_value
  )
)
