# Source-only plots comparing the current and relaxed domain eligibility rules.
# The caller loads gbe_theme.R and analyses/01_transposable_elements/ltr_genome_pgls.R.
# Species, genome estimates, and actual phylogenetic branch lengths are identical
# in both fits; only the retained element summaries change.

make_ltr_filter_pgls_panels <- function(comparison, dataset_dir) {
  prepared <- prepare_ltr_genome_data(dataset_dir)
  branch_ids <- c("five_domains", "four_domains")
  branch_labels <- c(five_domains = "Current: >=5 domains",
                     four_domains = "Sensitivity: >=4 domains")
  summary_columns <- c("species", "n_elements", "ratio_median",
                       "median_bootstrap_ci_low", "median_bootstrap_ci_high")
  if (!all(c("branch_id", summary_columns) %in% names(comparison$summaries))) {
    stop("Comparison must contain species medians, counts, and bootstrap intervals for both branches")
  }

  species_data <- list()
  model_summary <- list()
  predictions <- list()
  for (branch in branch_ids) {
    summary <- comparison$summaries[
      comparison$summaries$branch_id == branch, summary_columns, drop = FALSE
    ]
    if (nrow(summary) != 30L || anyDuplicated(summary$species) ||
        !setequal(summary$species, c(prepared$data$species, prepared$excluded$species))) {
      stop("Each branch must retain the same 30 species before the genome-size join")
    }
    data <- prepared$data
    joined <- match(data$species, summary$species)
    data[, summary_columns[-1]] <- summary[joined, summary_columns[-1], drop = FALSE]
    if (nrow(data) != 19L || anyNA(data) ||
        any(!is.finite(as.matrix(data[, summary_columns[-1]]))) ||
        any(as.matrix(data[, summary_columns[-1]]) <= 0) ||
        any(data$median_bootstrap_ci_low > data$ratio_median) ||
        any(data$median_bootstrap_ci_high < data$ratio_median)) {
      stop("Expected 19 positive species medians contained within valid bootstrap intervals")
    }
    if (branch == "five_domains" &&
        !isTRUE(all.equal(data[, summary_columns], prepared$data[, summary_columns],
                          tolerance = 1e-12, check.attributes = FALSE))) {
      stop("The current branch must exactly preserve the released primary species summaries")
    }
    data$log10_ratio <- log10(data$ratio_median)
    data$log10_genome <- log10(data$genome_size_gbp)
    rownames(data) <- data$species
    model <- fit_ltr_genome_model(data, prepared$tree, model_id = "log_log_lambda")

    x <- seq(min(data$log10_genome), max(data$log10_genome), length.out = 200L)
    design <- cbind(1, x)
    estimate <- as.numeric(design %*% stats::coef(model$fit))
    standard_error <- sqrt(rowSums((design %*% stats::vcov(model$fit)) * design))
    critical <- stats::qt(.975, nrow(data) - 2L)
    # This is a coefficient interval for the fitted geometric mean, not an
    # individual-species prediction interval or measurement-error propagation.
    predictions[[branch]] <- data.frame(
      branch_id = branch, genome_size_gbp = 10^x, ratio_fit = 10^estimate,
      ratio_ci_low = 10^(estimate - critical * standard_error),
      ratio_ci_high = 10^(estimate + critical * standard_error)
    )
    model_summary[[branch]] <- dplyr::mutate(model$summary, branch_id = branch, .before = 1)
    data$species_label <- paste("D.", data$species)
    species_data[[branch]] <- dplyr::mutate(data, branch_id = branch, .before = 1)
  }

  all_data <- dplyr::bind_rows(species_data)
  all_predictions <- dplyr::bind_rows(predictions)
  # Use the previous plot's limits when possible and expand both panels
  # together whenever any data or interval would otherwise be clipped.
  x_limits <- range(c(7.7, 23, all_data$genome_size_ci_low_gbp,
                       all_data$genome_size_ci_high_gbp))
  y_limits <- range(c(.23, 1.8, all_data$median_bootstrap_ci_low,
                       all_data$median_bootstrap_ci_high,
                       all_predictions$ratio_ci_low, all_predictions$ratio_ci_high))
  blue <- GBE_COLORS[["blue"]]
  # One shared destination map keeps every species label outside the dense
  # point cloud in both branches and at both standalone/collage dimensions.
  # Positions are presentation-only; tethers always start at observed medians.
  label_layout <- data.frame(
    species = c("campi", "orestes", "welteri", "perlapsus", "anicetus", "aeneus",
      "wrighti", "fuscus", "amphileucus", "auriculatus", "gvnigeusgwotli",
      "mavrokoilius", "apalachicolae", "tilleyi", "marmoratus", "intermedius",
      "monticola", "ocoee", "bairdi"),
    label_x = c(11.6, 11.6, 11.6, 11.6, 11.6, 9.4, 13.8, 15.7, 17.9, 20.9,
      13.3, 11.6, 16.2, 16.4, 19.3, 19.4, 19.4, 19.4, 19.4),
    label_y = c(.98, .84, .49, .41, .34, .70, 1.12, 1, 1.19, 1.05,
      .275, .58, .30, .40, .34, .83, .70, .59, .47),
    label_hjust = c(rep(1, 5), rep(.5, 6), 1, rep(.5, 3), rep(0, 4))
  )
  stopifnot(setequal(label_layout$species, prepared$data$species))
  plots <- list()
  for (i in seq_along(branch_ids)) {
    branch <- branch_ids[[i]]
    data <- species_data[[branch]]
    prediction <- predictions[[branch]]
    primary <- model_summary[[branch]]
    lambda_label <- if (primary$lambda < .001) "< 0.001" else sprintf("= %.3f", primary$lambda)
    annotation <- sprintf(
      "Pagel-lambda PGLS; n = %d species\nSlope = %.3f (95%% CI %.3f to %.3f)\nP = %.3f; lambda %s",
      primary$n_species, primary$slope, primary$ci_low, primary$ci_high,
      primary$p_value, lambda_label
    )
    annotation_x <- 10^(log10(x_limits[[1]]) + .055 * diff(log10(x_limits)))
    annotation_y <- 10^(log10(y_limits[[1]]) + .974 * diff(log10(y_limits)))
    label_data <- cbind(data, label_layout[match(data$species, label_layout$species), -1])
    # End lines just outside the label edge rather than crossing the text.
    label_data$line_x <- ifelse(label_data$label_hjust == 1, label_data$label_x * 1.015,
      ifelse(label_data$label_hjust == 0, label_data$label_x / 1.015, label_data$label_x))
    label_data$line_y <- ifelse(label_data$label_hjust == .5,
      label_data$label_y * ifelse(label_data$label_y > label_data$ratio_median, 1 / 1.035, 1.035),
      label_data$label_y)
    plots[[branch]] <- ggplot2::ggplot(data, ggplot2::aes(genome_size_gbp, ratio_median)) +
      ggplot2::geom_ribbon(data = prediction,
        ggplot2::aes(x = genome_size_gbp, ymin = ratio_ci_low, ymax = ratio_ci_high),
        inherit.aes = FALSE, fill = blue, alpha = .12) +
      ggplot2::geom_hline(yintercept = 1, linetype = "dotted", linewidth = .35, colour = "#888888") +
      ggplot2::geom_segment(
        ggplot2::aes(x = genome_size_ci_low_gbp, xend = genome_size_ci_high_gbp, yend = ratio_median),
        colour = "#BDBDBD", linewidth = .3) +
      ggplot2::geom_segment(
        ggplot2::aes(y = median_bootstrap_ci_low, yend = median_bootstrap_ci_high, xend = genome_size_gbp),
        colour = "#BDBDBD", linewidth = .3) +
      ggplot2::geom_line(data = prediction, ggplot2::aes(x = genome_size_gbp, y = ratio_fit),
        inherit.aes = FALSE, colour = blue, linewidth = .7) +
      ggplot2::geom_point(shape = 21, fill = "white", colour = GBE_COLORS[["charcoal"]],
        size = 1.9, stroke = .45) +
      ggplot2::geom_segment(data = label_data,
        ggplot2::aes(x = genome_size_gbp, y = ratio_median, xend = line_x, yend = line_y),
        inherit.aes = FALSE, colour = "#999999", linewidth = .2) +
      ggplot2::geom_text(data = label_data,
        ggplot2::aes(x = label_x, y = label_y, label = species_label, hjust = label_hjust),
        inherit.aes = FALSE, family = GBE_FONT, fontface = "italic", size = GBE_LABEL_SIZE_MM) +
      ggplot2::annotate("text", x = annotation_x, y = annotation_y,
        label = annotation, hjust = 0, vjust = 1,
        family = GBE_FONT, size = GBE_LABEL_SIZE_MM, lineheight = 1.1) +
      ggplot2::scale_x_log10(breaks = c(8, 10, 12, 15, 18, 22), limits = x_limits,
        expand = ggplot2::expansion(mult = .02)) +
      ggplot2::scale_y_log10(breaks = c(.2, .25, .4, .6, .8, 1, 1.5, 2, 3), limits = y_limits,
        labels = scales::label_number(), expand = ggplot2::expansion(mult = .02)) +
      ggplot2::labs(
        subtitle = sprintf("%s | %d elements across %d species",
                           branch_labels[[branch]], sum(data$n_elements), nrow(data)),
        x = "Estimated genome size (Gbp; log scale)",
        y = "Species median LTR terminal-to-internal\ndepth ratio (log scale)",
        tag = letters[[i]],
        caption = paste(
          "Blue line/band: PGLS fit and 95% coefficient interval; dotted line: equal depth.",
          "Gray bars: estimated genome-size intervals and element-bootstrap median intervals.",
          "Both branches retain regional coverage >=80% and within-species IQR filtering.",
          sep = "\n"
        )
      ) +
      theme_gbe() + theme_gbe_framed_facets() +
      ggplot2::theme(
        plot.subtitle = ggplot2::element_text(family = GBE_FONT, size = GBE_BASE_SIZE,
          face = "bold", margin = ggplot2::margin(b = 5, unit = "pt")),
        plot.tag.position = "topleft",
        plot.margin = ggplot2::margin(6, 8, 6, 6, "pt")
      )
  }
  list(plots = plots, model_summary = dplyr::bind_rows(model_summary),
       species_data = all_data, predictions = all_predictions, tree = prepared$tree,
       excluded_species = prepared$excluded,
       axes_common = list(x_limits = x_limits, y_limits = y_limits))
}
