#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(ape))

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
if (length(file_arg) == 0) {
  stop("Could not determine script path from commandArgs().")
}

script_path <- normalizePath(sub("^--file=", "", file_arg[[1]]))
project_root <- dirname(dirname(dirname(script_path)))
derived_dir <- file.path(project_root, "path_analysis", "data", "derived")

tree_file <- file.path(project_root, "input_data", "phylogeny", "desmo900dated_test.tre")
traits_file <- file.path(derived_dir, "organismal_traits_curated.csv")
output_long <- file.path(derived_dir, "phylogenetic_trait_imputation_long.csv")
output_wide <- file.path(derived_dir, "organismal_traits_phylo_inference.csv")
output_summary <- file.path(derived_dir, "phylogenetic_trait_imputation_summary.csv")

if (!file.exists(tree_file)) {
  stop(sprintf("Tree file not found: %s", tree_file))
}
if (!file.exists(traits_file)) {
  stop(sprintf("Trait file not found: %s", traits_file))
}

confidence_to_score <- function(x) {
  scores <- c(low = 1, medium = 2, high = 3)
  unname(scores[as.character(x)])
}

safe_collapse <- function(x, digits = 4) {
  x <- x[!is.na(x)]
  if (length(x) == 0) {
    return(NA_character_)
  }
  if (is.numeric(x)) {
    return(paste(format(round(x, digits), trim = TRUE, scientific = FALSE), collapse = ";"))
  }
  paste(as.character(x), collapse = ";")
}

safe_numeric <- function(x) {
  out <- suppressWarnings(as.numeric(x))
  ifelse(is.na(out), NA_real_, out)
}

normalize_missing_strings <- function(df) {
  for (col in names(df)) {
    if (is.character(df[[col]])) {
      trimmed <- trimws(df[[col]])
      df[[col]][trimmed == ""] <- NA_character_
    }
  }
  df
}

normalize_tree_tip_labels <- function(labels) {
  labels <- sub("^D\\.\\s*", "", labels)
  labels <- sub("_.*$", "", labels)
  trimws(tolower(labels))
}

predict_from_neighbors <- function(target_species, train_species, train_values, dist_mat, kind,
                                   k = 5, log_transform = FALSE) {
  if (!(target_species %in% rownames(dist_mat))) {
    return(NULL)
  }
  if (length(train_species) == 0) {
    return(NULL)
  }

  distances <- dist_mat[target_species, train_species]
  distances <- distances[is.finite(distances)]
  if (length(distances) == 0) {
    return(NULL)
  }

  ordered_species <- names(sort(distances, method = "radix"))
  chosen_species <- ordered_species[seq_len(min(k, length(ordered_species)))]
  chosen_distances <- as.numeric(distances[chosen_species])
  weights <- 1 / (chosen_distances + 1e-6)^2
  weights <- weights / sum(weights)
  chosen_values <- train_values[chosen_species]

  if (kind == "continuous") {
    raw_values <- as.numeric(chosen_values)
    transformed_values <- if (log_transform) log(raw_values) else raw_values
    pred_t <- sum(weights * transformed_values)
    pred_raw <- if (log_transform) exp(pred_t) else pred_t
    weighted_sd <- sqrt(sum(weights * (raw_values - pred_raw)^2))
    relative_sd <- weighted_sd / max(pred_raw, 1e-6)
    consistency <- 1 / (1 + relative_sd)
    return(list(
      predicted = pred_raw,
      support = consistency,
      support_label = "neighbor_consistency_score",
      top_donor = chosen_species[[which.max(weights)]],
      top_distance = min(chosen_distances),
      donor_species = chosen_species,
      donor_distances = chosen_distances,
      donor_weights = weights
    ))
  }

  states <- as.character(chosen_values)
  weight_by_state <- tapply(weights, states, sum)
  weight_by_state <- sort(weight_by_state, decreasing = TRUE)
  predicted_state <- names(weight_by_state)[[1]]
  support <- as.numeric(weight_by_state[[1]])

  list(
    predicted = predicted_state,
    support = support,
    support_label = "weighted_vote_share",
    top_donor = chosen_species[[which.max(weights)]],
    top_distance = min(chosen_distances),
    donor_species = chosen_species,
    donor_distances = chosen_distances,
    donor_weights = weights
  )
}

compute_cv_summary <- function(train_species, train_values, dist_mat, kind,
                               k = 5, log_transform = FALSE) {
  if (length(train_species) < 3) {
    return(list(
      n_cv = 0L,
      metric_primary_type = if (kind == "continuous") "mae" else "accuracy",
      metric_primary_value = NA_real_,
      metric_secondary_type = if (kind == "continuous") "rmse" else "mean_vote_share",
      metric_secondary_value = NA_real_
    ))
  }

  predictions <- list()
  truths <- c()
  supports <- c()

  for (species in train_species) {
    other_species <- setdiff(train_species, species)
    other_values <- train_values[other_species]
    pred <- predict_from_neighbors(
      target_species = species,
      train_species = other_species,
      train_values = other_values,
      dist_mat = dist_mat,
      kind = kind,
      k = k,
      log_transform = log_transform
    )
    if (is.null(pred)) {
      next
    }
    predictions[[species]] <- pred$predicted
    truths[[species]] <- train_values[[species]]
    supports[[species]] <- pred$support
  }

  if (length(predictions) == 0) {
    return(list(
      n_cv = 0L,
      metric_primary_type = if (kind == "continuous") "mae" else "accuracy",
      metric_primary_value = NA_real_,
      metric_secondary_type = if (kind == "continuous") "rmse" else "mean_vote_share",
      metric_secondary_value = NA_real_
    ))
  }

  pred_vec <- unlist(predictions, use.names = TRUE)
  truth_vec <- unlist(truths, use.names = TRUE)
  support_vec <- unlist(supports, use.names = TRUE)

  if (kind == "continuous") {
    errors <- abs(as.numeric(pred_vec) - as.numeric(truth_vec))
    rmse <- sqrt(mean((as.numeric(pred_vec) - as.numeric(truth_vec))^2))
    return(list(
      n_cv = length(errors),
      metric_primary_type = "mae",
      metric_primary_value = mean(errors),
      metric_secondary_type = "rmse",
      metric_secondary_value = rmse
    ))
  }

  return(list(
    n_cv = length(pred_vec),
    metric_primary_type = "accuracy",
    metric_primary_value = mean(as.character(pred_vec) == as.character(truth_vec)),
    metric_secondary_type = "mean_vote_share",
    metric_secondary_value = mean(as.numeric(support_vec))
  ))
}

build_trait_outputs <- function(df, dist_mat, trait_col, confidence_col, kind,
                                k = 5, log_transform = FALSE) {
  df <- df
  observed_values <- df[[trait_col]]
  confidence_scores <- confidence_to_score(df[[confidence_col]])
  in_tree <- df$species %in% rownames(dist_mat)

  train_mask <- in_tree & !is.na(observed_values) & !is.na(confidence_scores) & confidence_scores >= 2
  train_species <- df$species[train_mask]
  train_values <- observed_values[train_mask]
  names(train_values) <- train_species

  cv_summary <- compute_cv_summary(
    train_species = train_species,
    train_values = train_values,
    dist_mat = dist_mat,
    kind = kind,
    k = k,
    log_transform = log_transform
  )

  long_rows <- vector("list", length = nrow(df))
  wide_rows <- vector("list", length = nrow(df))

  for (i in seq_len(nrow(df))) {
    species <- df$species[[i]]
    observed_value <- observed_values[[i]]
    observed_conf <- df[[confidence_col]][[i]]

    status <- "observed_source_backed"
    recommendation <- "keep_observed_source_backed_value"
    inferred_value <- NA
    support_value <- NA_real_
    support_label <- NA_character_
    top_donor <- NA_character_
    top_distance <- NA_real_
    neighbor_species <- NA_character_
    neighbor_distances <- NA_character_
    neighbor_weights <- NA_character_
    method <- "phylogenetic_inverse_distance_knn"

    if (!in_tree[[i]]) {
      status <- "not_in_tree"
      recommendation <- "no_phylogenetic_inference_available"
      method <- NA_character_
    } else if (is.na(observed_value) || is.na(observed_conf) || observed_conf == "low") {
      pred <- predict_from_neighbors(
        target_species = species,
        train_species = train_species,
        train_values = train_values,
        dist_mat = dist_mat,
        kind = kind,
        k = k,
        log_transform = log_transform
      )

      if (is.null(pred)) {
        status <- "insufficient_training"
        recommendation <- "no_phylogenetic_inference_available"
      } else {
        inferred_value <- pred$predicted
        support_value <- pred$support
        support_label <- pred$support_label
        top_donor <- pred$top_donor
        top_distance <- pred$top_distance
        neighbor_species <- safe_collapse(pred$donor_species)
        neighbor_distances <- safe_collapse(pred$donor_distances)
        neighbor_weights <- safe_collapse(pred$donor_weights)

        if (is.na(observed_value)) {
          status <- "phylo_inferred_for_missing"
          recommendation <- "sensitivity_only_fill_missing_with_explicit_phylo_flag"
        } else if (kind == "continuous") {
          allowed_error <- cv_summary$metric_primary_value
          if (is.na(allowed_error)) {
            allowed_error <- 0
          }
          if (abs(as.numeric(observed_value) - as.numeric(inferred_value)) <= allowed_error) {
            status <- "low_confidence_observed_supported_by_phylogeny"
            recommendation <- "keep_observed_but_mark_phylo_support"
          } else {
            status <- "low_confidence_observed_discordant_with_phylogeny"
            recommendation <- "manual_review_priority_do_not_overwrite"
          }
        } else if (as.character(observed_value) == as.character(inferred_value)) {
          status <- "low_confidence_observed_supported_by_phylogeny"
          recommendation <- "keep_observed_but_mark_phylo_support"
        } else {
          status <- "low_confidence_observed_discordant_with_phylogeny"
          recommendation <- "manual_review_priority_do_not_overwrite"
        }
      }
    }

    observed_state <- if (!in_tree[[i]]) {
      "not_in_tree"
    } else if (is.na(observed_value)) {
      "missing"
    } else if (!is.na(observed_conf) && observed_conf == "low") {
      "low_confidence"
    } else {
      "source_backed"
    }

    long_rows[[i]] <- data.frame(
      species = species,
      trait_name = trait_col,
      trait_type = kind,
      observed_value = if (is.na(observed_value)) NA_character_ else as.character(observed_value),
      observed_confidence = if (is.na(observed_conf)) NA_character_ else as.character(observed_conf),
      observed_state = observed_state,
      phylo_inferred_value = if (is.na(inferred_value)) NA_character_ else as.character(inferred_value),
      phylo_inference_status = status,
      phylo_recommendation = recommendation,
      phylo_method = method,
      phylo_training_n = length(train_species),
      phylo_top_donor_species = top_donor,
      phylo_top_donor_distance = top_distance,
      phylo_neighbor_species = neighbor_species,
      phylo_neighbor_distances = neighbor_distances,
      phylo_neighbor_weights = neighbor_weights,
      phylo_support_label = support_label,
      phylo_support_value = support_value,
      cv_primary_metric_type = cv_summary$metric_primary_type,
      cv_primary_metric_value = cv_summary$metric_primary_value,
      cv_secondary_metric_type = cv_summary$metric_secondary_type,
      cv_secondary_metric_value = cv_summary$metric_secondary_value,
      stringsAsFactors = FALSE
    )

    wide_rows[[i]] <- data.frame(
      species = species,
      observed_value = observed_value,
      observed_confidence = observed_conf,
      inferred_value = if (is.na(inferred_value)) NA else inferred_value,
      status = status,
      recommendation = recommendation,
      method = method,
      training_n = length(train_species),
      top_donor_species = top_donor,
      top_donor_distance = top_distance,
      neighbor_species = neighbor_species,
      neighbor_support = support_value,
      neighbor_support_label = support_label,
      cv_primary_metric_type = cv_summary$metric_primary_type,
      cv_primary_metric_value = cv_summary$metric_primary_value,
      stringsAsFactors = FALSE
    )
  }

  long_df <- do.call(rbind, long_rows)
  wide_df <- do.call(rbind, wide_rows)

  rename_map <- c(
    observed_value = sprintf("phylo_%s_observed", trait_col),
    observed_confidence = sprintf("phylo_%s_observed_confidence", trait_col),
    inferred_value = sprintf("phylo_%s", trait_col),
    status = sprintf("phylo_%s_status", trait_col),
    recommendation = sprintf("phylo_%s_recommendation", trait_col),
    method = sprintf("phylo_%s_method", trait_col),
    training_n = sprintf("phylo_%s_training_n", trait_col),
    top_donor_species = sprintf("phylo_%s_top_donor_species", trait_col),
    top_donor_distance = sprintf("phylo_%s_top_donor_distance", trait_col),
    neighbor_species = sprintf("phylo_%s_neighbor_species", trait_col),
    neighbor_support = sprintf("phylo_%s_neighbor_support", trait_col),
    neighbor_support_label = sprintf("phylo_%s_neighbor_support_label", trait_col),
    cv_primary_metric_type = sprintf("phylo_%s_cv_primary_metric_type", trait_col),
    cv_primary_metric_value = sprintf("phylo_%s_cv_primary_metric_value", trait_col)
  )

  names(wide_df)[match(names(rename_map), names(wide_df))] <- unname(rename_map)

  summary_df <- data.frame(
    trait_name = trait_col,
    trait_type = kind,
    training_species_n = length(train_species),
    cv_n = cv_summary$n_cv,
    cv_primary_metric_type = cv_summary$metric_primary_type,
    cv_primary_metric_value = cv_summary$metric_primary_value,
    cv_secondary_metric_type = cv_summary$metric_secondary_type,
    cv_secondary_metric_value = cv_summary$metric_secondary_value,
    n_missing = sum(in_tree & is.na(observed_values)),
    n_low_confidence = sum(in_tree & !is.na(observed_values) & !is.na(confidence_scores) & confidence_scores < 2),
    stringsAsFactors = FALSE
  )

  list(long = long_df, wide = wide_df, summary = summary_df)
}

message("Reading tree and curated trait table...")
tree <- read.tree(tree_file)
clean_tip_labels <- normalize_tree_tip_labels(tree$tip.label)
duplicated_tips <- duplicated(clean_tip_labels)
if (any(duplicated_tips)) {
  message(sprintf(
    "Dropping %d duplicate phylogeny tips after species-label normalization.",
    sum(duplicated_tips)
  ))
  tree <- drop.tip(tree, tree$tip.label[duplicated_tips])
  clean_tip_labels <- clean_tip_labels[!duplicated_tips]
}
tree$tip.label <- clean_tip_labels
dist_mat <- cophenetic.phylo(tree)

traits <- read.csv(
  traits_file,
  stringsAsFactors = FALSE,
  check.names = FALSE,
  na.strings = c("", "NA")
)
traits <- normalize_missing_strings(traits)
traits <- traits[order(traits$species), ]

trait_specs <- list(
  list(trait_col = "body_size_proxy_mm", confidence_col = "body_size_proxy_confidence", kind = "continuous", log_transform = TRUE),
  list(trait_col = "development_mode", confidence_col = "development_confidence", kind = "discrete", log_transform = FALSE),
  list(trait_col = "aquaticity_index", confidence_col = "lifestyle_confidence", kind = "discrete", log_transform = FALSE),
  list(trait_col = "microhabitat_class", confidence_col = "lifestyle_confidence", kind = "discrete", log_transform = FALSE)
)

long_blocks <- list()
wide_blocks <- list(data.frame(species = traits$species, stringsAsFactors = FALSE))
summary_blocks <- list()

for (spec in trait_specs) {
  message(sprintf("Imputing %s...", spec$trait_col))
  result <- build_trait_outputs(
    df = traits,
    dist_mat = dist_mat,
    trait_col = spec$trait_col,
    confidence_col = spec$confidence_col,
    kind = spec$kind,
    log_transform = spec$log_transform
  )
  long_blocks[[spec$trait_col]] <- result$long
  wide_blocks[[length(wide_blocks) + 1]] <- result$wide
  summary_blocks[[spec$trait_col]] <- result$summary
}

long_df <- do.call(rbind, long_blocks)
wide_df <- Reduce(function(x, y) merge(x, y, by = "species", all = TRUE), wide_blocks)
summary_df <- do.call(rbind, summary_blocks)

write.csv(long_df, output_long, row.names = FALSE, na = "")
write.csv(wide_df, output_wide, row.names = FALSE, na = "")
write.csv(summary_df, output_summary, row.names = FALSE, na = "")

message(sprintf("Wrote %s", output_long))
message(sprintf("Wrote %s", output_wide))
message(sprintf("Wrote %s", output_summary))
