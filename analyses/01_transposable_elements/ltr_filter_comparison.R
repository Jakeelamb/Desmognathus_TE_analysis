# Downstream, source-only sensitivity analysis of the released LTR audit table.
# No canonical release files are written or upstream analyses rerun.
# The existing five-domain branch is copied exactly from S14/S15. The broader
# branch uses the same eligibility and filtering sequence with >=4 domains.
prepare_ltr_filter_comparison <- function(dataset_dir) {
  element_path <- file.path(
    dataset_dir, "Supplementary_Data_S14_ltr_terminal_internal_elements.csv"
  )
  summary_path <- file.path(
    dataset_dir, "Supplementary_Data_S15_ltr_terminal_internal_species_robustness.csv"
  )
  elements <- utils::read.csv(element_path, check.names = FALSE,
                             stringsAsFactors = FALSE)
  released <- utils::read.csv(summary_path, check.names = FALSE,
                             stringsAsFactors = FALSE)
  ratio_column <- "ratio_terminal_internal_all_positions"
  required <- c(
    "species", "element_id", "Domains", "Order", "Superfamily",
    "element length", "tesorter_annotation_found", "domain_count",
    "domain_gate_keep", "iqr_filter_keep", "coverage_ge_80pct",
    "primary_ltr_element", ratio_column,
    "left_ltr_positive_coverage_fraction",
    "right_ltr_positive_coverage_fraction",
    "internal_positive_coverage_fraction"
  )
  missing <- setdiff(required, names(elements))
  if (length(missing)) stop("LTR table lacks columns: ", paste(missing, collapse = ", "))
  summary_columns <- c("species", "n_elements", "ratio_median",
                       "median_bootstrap_ci_low", "median_bootstrap_ci_high")
  missing <- setdiff(c("analysis_branch", summary_columns), names(released))
  if (length(missing)) stop("LTR summaries lack columns: ", paste(missing, collapse = ", "))
  if (nrow(elements) != 1086L || anyDuplicated(elements$element_id)) {
    stop("Expected 1,086 unique usable historical LTR elements")
  }
  if (anyNA(elements$species) || any(!nzchar(elements$species))) {
    stop("LTR species identifiers must be present")
  }
  ratios <- elements[[ratio_column]]
  if (!is.numeric(ratios) || any(!is.finite(ratios) | ratios <= 0)) {
    stop("Every LTR ratio must be finite and positive")
  }
  as_flag <- function(values, label) {
    flags <- as.logical(values)
    if (anyNA(flags)) stop("Missing or invalid flag: ", label)
    flags
  }
  count_domains <- function(value) {
    if (is.na(value) || !nzchar(trimws(value))) stop("Missing TEsorter Domains value")
    value <- trimws(value)
    if (tolower(value) == "none") return(0L)
    tokens <- strsplit(value, "[[:space:]]+")[[1L]]
    if (any(!grepl("^[^|]+[|][^|]+$", tokens))) {
      stop("Malformed TEsorter DOMAIN|MODEL annotation")
    }
    length(tokens)
  }
  domain_count <- vapply(elements$Domains, count_domains, integer(1), USE.NAMES = FALSE)
  if (!identical(as.integer(elements$domain_count), domain_count)) {
    stop("Released domain counts disagree with whitespace-delimited annotations")
  }
  annotation_found <- as_flag(elements$tesorter_annotation_found,
                              "tesorter_annotation_found")
  common_eligible <- elements[["element length"]] >= 3000 & annotation_found &
    elements$Order == "LTR" & elements$Superfamily == "Gypsy"
  if (anyNA(common_eligible)) stop("Missing LTR eligibility fields")
  five_eligible <- common_eligible & domain_count >= 5L
  four_eligible <- common_eligible & domain_count >= 4L
  if (!identical(five_eligible, as_flag(elements$domain_gate_keep, "domain_gate_keep"))) {
    stop("Released five-domain eligibility does not reproduce")
  }
  coverage_columns <- c("left_ltr_positive_coverage_fraction",
                         "right_ltr_positive_coverage_fraction",
                         "internal_positive_coverage_fraction")
  coverage_values <- as.matrix(elements[, coverage_columns])
  if (any(!is.finite(coverage_values) | coverage_values < 0 | coverage_values > 1)) {
    stop("Regional positive-depth coverage must be finite and between zero and one")
  }
  coverage_keep <- apply(coverage_values >= 0.8, 1L, all)
  if (!identical(coverage_keep, as_flag(elements$coverage_ge_80pct, "coverage_ge_80pct"))) {
    stop("Released regional coverage flags do not reproduce")
  }

  # Fences use all eligible ratios BEFORE intersecting with regional coverage.
  # R type 7 is the linear-interpolation convention used by the released pandas
  # quantiles; equality with either fence is retained.
  q1 <- q3 <- rep(NA_real_, nrow(elements))
  for (species in sort(unique(elements$species[four_eligible]))) {
    rows <- which(four_eligible & elements$species == species)
    quartiles <- stats::quantile(ratios[rows], c(0.25, 0.75), type = 7,
                                names = FALSE)
    q1[rows] <- quartiles[1L]
    q3[rows] <- quartiles[2L]
  }
  iqr <- q3 - q1
  lower <- q1 - 1.5 * iqr
  upper <- q3 + 1.5 * iqr
  four_iqr_keep <- four_eligible & ratios >= lower & ratios <= upper
  four_keep <- four_iqr_keep & coverage_keep
  five_iqr_keep <- as_flag(elements$iqr_filter_keep, "iqr_filter_keep")
  five_keep <- as_flag(elements$primary_ltr_element, "primary_ltr_element")
  if (!identical(five_keep, five_iqr_keep & coverage_keep)) {
    stop("Released primary branch is not the IQR/coverage intersection")
  }
  observed <- c(sum(five_eligible), sum(five_eligible & coverage_keep),
                sum(five_iqr_keep), sum(five_keep), sum(four_eligible),
                sum(four_eligible & coverage_keep), sum(four_keep))
  expected <- c(408L, 407L, 381L, 380L, 1083L, 1080L, 989L)
  if (!identical(as.integer(observed), expected)) {
    stop("Frozen comparison cohort counts changed: ", paste(observed, collapse = ", "))
  }

  decisions <- elements
  decisions$four_domain_count <- domain_count
  decisions$four_domain_eligible <- four_eligible
  decisions$four_domain_q1 <- q1
  decisions$four_domain_q3 <- q3
  decisions$four_domain_iqr <- iqr
  decisions$four_domain_lower_bound <- lower
  decisions$four_domain_upper_bound <- upper
  decisions$four_domain_iqr_keep <- four_iqr_keep
  decisions$four_domain_coverage_keep <- four_eligible & coverage_keep
  decisions$four_domain_keep <- four_keep
  decisions$four_domain_iqr_direction <- ifelse(
    !four_eligible, "not_eligible",
    ifelse(ratios < lower, "lower_outlier",
           ifelse(ratios > upper, "upper_outlier", "retained"))
  )
  decisions$four_domain_eligibility_reason <- ifelse(
    elements[["element length"]] < 3000, "element_shorter_than_3000bp",
    ifelse(!annotation_found, "missing_exact_tesorter_annotation",
           ifelse(elements$Order != "LTR", "not_ltr_order",
                  ifelse(elements$Superfamily != "Gypsy", "not_gypsy_superfamily",
                         ifelse(domain_count < 4L, "fewer_than_4_tesorter_domain_annotations",
                                "eligible"))))
  )
  branches <- list(five_domains = five_keep, four_domains = four_keep)
  selected <- do.call(rbind, lapply(names(branches), function(branch_id) {
    rows <- elements[branches[[branch_id]], , drop = FALSE]
    rows$branch_id <- branch_id
    # Original S14 flags above remain the RELEASED decisions even in the broader
    # branch. Use branch_retained or decisions$four_domain_keep for this analysis.
    rows$branch_retained <- TRUE
    rows
  }))
  rownames(selected) <- NULL

  # Preserve released NumPy PCG64 bootstrap intervals exactly for the reference.
  five_summary <- released[
    released$analysis_branch == "primary_iqr_filtered_coverage_ge_80pct",
    summary_columns, drop = FALSE
  ]
  five_summary <- five_summary[order(five_summary$species), , drop = FALSE]
  if (nrow(five_summary) != 30L || anyDuplicated(five_summary$species)) {
    stop("Expected exactly 30 released primary species summaries")
  }
  for (row in seq_len(nrow(five_summary))) {
    values <- ratios[five_keep & elements$species == five_summary$species[row]]
    if (length(values) != five_summary$n_elements[row] ||
        !isTRUE(all.equal(stats::median(values), five_summary$ratio_median[row],
                          tolerance = 1e-13))) {
      stop("Released primary species count/median does not reproduce: ",
           five_summary$species[row])
    }
  }
  five_summary$branch_id <- "five_domains"
  five_summary$bootstrap_replicates <- 2000L
  five_summary$bootstrap_method <- "released_S15_NumPy_PCG64_percentile_linear"
  five_summary$bootstrap_seed <- NA_integer_
  five_summary$bootstrap_seed_text <- paste0(
    five_summary$species, ":primary_iqr_filtered_coverage_ge_80pct:median"
  )

  # New branch: same percentile-bootstrap estimand and replicate count, with a
  # documented R RNG rather than claiming byte-equivalence to NumPy PCG64.
  # UTF-8 polynomial hash modulo 2^31-1 provides a reproducible species seed.
  species_seed <- function(species) {
    seed <- 20260924
    for (byte in utf8ToInt(paste0(species, ":four_domains:median"))) {
      seed <- (seed * 131 + byte) %% 2147483647
    }
    as.integer(if (seed == 0) 1 else seed)
  }
  bootstrap_interval <- function(values, seed) {
    old_kind <- RNGkind()
    had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    on.exit({
      do.call(RNGkind, as.list(old_kind))
      if (had_seed) {
        assign(".Random.seed", old_seed, envir = .GlobalEnv)
      } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    }, add = TRUE)
    RNGkind("Mersenne-Twister", "Inversion", "Rejection")
    set.seed(seed)
    estimates <- replicate(2000L, stats::median(
      values[sample.int(length(values), length(values), replace = TRUE)]
    ))
    stats::quantile(estimates, c(0.025, 0.975), type = 7, names = FALSE)
  }
  four_summary <- do.call(rbind, lapply(
    sort(unique(elements$species[four_keep])), function(species) {
      rows <- which(four_keep & elements$species == species)
      # Stable ordering makes the new intervals independent of input row order.
      rows <- rows[order(elements$element_id[rows])]
      values <- ratios[rows]
      seed <- species_seed(species)
      interval <- bootstrap_interval(values, seed)
      data.frame(
        species = species, n_elements = length(values),
        ratio_median = stats::median(values),
        median_bootstrap_ci_low = interval[1L],
        median_bootstrap_ci_high = interval[2L],
        branch_id = "four_domains", bootstrap_replicates = 2000L,
        bootstrap_method = "R_Mersenne_Twister_Rejection_percentile_type7",
        bootstrap_seed = seed,
        bootstrap_seed_text = paste0(species, ":four_domains:median"),
        stringsAsFactors = FALSE
      )
    }
  ))
  if (nrow(four_summary) != 30L ||
      !identical(five_summary$species, four_summary$species) ||
      anyDuplicated(four_summary$bootstrap_seed)) {
    stop("Comparison must retain all 30 species with unique deterministic seeds")
  }
  summaries <- rbind(five_summary, four_summary)
  summaries <- summaries[, c("branch_id", setdiff(names(summaries), "branch_id"))]
  rownames(summaries) <- NULL
  summary_values <- as.matrix(summaries[, c("ratio_median", "median_bootstrap_ci_low",
                                           "median_bootstrap_ci_high")])
  if (any(!is.finite(summary_values) | summary_values <= 0) ||
      any(summaries$median_bootstrap_ci_low > summaries$ratio_median) ||
      any(summaries$median_bootstrap_ci_high < summaries$ratio_median)) {
    stop("Species summaries must have finite positive medians and enclosing intervals")
  }
  counts <- data.frame(
    branch_id = names(branches), historical_usable = nrow(elements),
    eligible = c(sum(five_eligible), sum(four_eligible)),
    coverage_qualified = c(sum(five_eligible & coverage_keep),
                           sum(four_eligible & coverage_keep)),
    iqr_retained = c(sum(five_iqr_keep), sum(four_iqr_keep)),
    retained = c(sum(five_keep), sum(four_keep)),
    n_species = c(length(unique(elements$species[five_keep])),
                  length(unique(elements$species[four_keep])))
  )
  paired <- merge(five_summary[, summary_columns], four_summary[, summary_columns],
                  by = "species", suffixes = c("_five_domains", "_four_domains"),
                  sort = TRUE)
  paired$n_elements_change <- paired$n_elements_four_domains - paired$n_elements_five_domains
  paired$ratio_median_change <- paired$ratio_median_four_domains - paired$ratio_median_five_domains
  paired$ratio_median_change_pct <- 100 * paired$ratio_median_change / paired$ratio_median_five_domains
  list(elements = selected, summaries = summaries, decisions = decisions,
       counts = counts, paired_species_summary = paired)
}
