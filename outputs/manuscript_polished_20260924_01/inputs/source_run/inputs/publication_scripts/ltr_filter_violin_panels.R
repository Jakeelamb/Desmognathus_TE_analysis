# Matched views of the current and four-domain sensitivity LTR distributions.
make_ltr_filter_violin_panels <- function(comparison, dataset_dir) {
  source_tree <- ape::ladderize(ape::read.tree(file.path(
    dataset_dir, "trees/desmo900dated_source46.tre"
  )), right = FALSE)
  canonical <- source_tree$tip.label[order(ape::node.height(source_tree)[
    seq_len(ape::Ntip(source_tree))
  ])]
  species <- canonical[canonical %in% comparison$summaries$species]
  if (length(species) != 30L) stop("Comparison requires the same 30 species")
  rows <- data.frame(species = species, y = seq_along(species),
                     species_label = paste("D.", species))
  traits <- readr::read_csv(file.path(dataset_dir,
    "Supplementary_Data_S04_path_analysis_traits24.csv"), show_col_types = FALSE)
  rows$genome_size_gbp <- traits$genome_size_pg_fuscus_anchored[
    match(rows$species, traits$tree_tip)
  ] * .978
  elements <- comparison$elements
  elements$y <- rows$y[match(elements$species, rows$species)]
  elements$estimate <- log2(elements$ratio_terminal_internal_all_positions)
  summaries <- comparison$summaries
  summaries$y <- rows$y[match(summaries$species, rows$species)]
  summaries$estimate <- log2(summaries$ratio_median)
  summaries$ci_low <- log2(summaries$median_bootstrap_ci_low)
  summaries$ci_high <- log2(summaries$median_bootstrap_ci_high)
  bandwidth <- stats::bw.nrd0(elements$estimate[elements$branch_id == "five_domains"])
  limits <- log2(c(.025, 2.5))
  if (any(!is.finite(elements$estimate)) || any(elements$estimate < limits[1]) ||
      any(elements$estimate > limits[2]) || anyNA(elements$y) ||
      any(summaries$ci_low < limits[1]) || any(summaries$ci_high > limits[2])) {
    stop("Shared comparison axis must include all elements and intervals")
  }
  branch_ids <- c("five_domains", "four_domains")
  make_panel <- function(branch_id, show_species, include_genome) {
    current <- elements[elements$branch_id == branch_id, ]
    summary <- summaries[summaries$branch_id == branch_id, ]
    count <- summary$n_elements[match(rows$species, summary$species)]
    right_label <- if (include_genome) paste0(
      ifelse(is.na(rows$genome_size_gbp), "n/a", sprintf("%.1f", rows$genome_size_gbp)),
      " / ", count
    ) else as.character(count)
    right_title <- if (include_genome) "Estimated genome size (Gbp*) / elements (n)" else "n"
    subtitle <- sprintf("%s: >=%d domains\n%d elements; 30 species",
      if (branch_id == "five_domains") "Current" else "Sensitivity",
      if (branch_id == "five_domains") 5 else 4, nrow(current))
    p <- ggplot2::ggplot(current, ggplot2::aes(x = estimate, y = y)) +
      ggplot2::geom_violin(ggplot2::aes(group = species), orientation = "y",
        position = "identity", bw = bandwidth, trim = TRUE, scale = "width", width = .8,
        quantiles = NULL, fill = "white", colour = "#B3B3B3", linewidth = .35) +
      ggplot2::geom_vline(xintercept = 0, linetype = 3, linewidth = .5,
                         colour = GBE_COLORS[["charcoal"]]) +
      ggplot2::geom_point(position = ggplot2::position_jitter(
        width = 0, height = .15, seed = 20260828), shape = 21, size = .95,
        stroke = .25, fill = "white", colour = GBE_COLORS[["charcoal"]], alpha = .75) +
      ggplot2::geom_segment(data = summary,
        ggplot2::aes(x = ci_low, xend = ci_high, y = y, yend = y),
        inherit.aes = FALSE, linewidth = .35, colour = GBE_COLORS[["blue"]]) +
      ggplot2::geom_segment(data = summary,
        ggplot2::aes(x = estimate, xend = estimate, y = y - .25, yend = y + .25),
        inherit.aes = FALSE, linewidth = .9, colour = GBE_COLORS[["blue"]]) +
      ggplot2::scale_x_continuous(limits = limits,
        breaks = log2(c(.025, .1, .5, 1, 2.5)),
        labels = c("0.025", "0.1", "0.5", "1", "2.5"), expand = c(0, 0)) +
      ggplot2::scale_y_continuous(limits = c(.5, 30.5), breaks = rows$y,
        labels = rows$species_label, expand = c(0, 0),
        sec.axis = ggplot2::dup_axis(breaks = rows$y, labels = right_label, name = right_title)) +
      ggplot2::labs(x = "Terminal-to-internal depth ratio\n(log2 scale; 1 = equal depth)",
        y = NULL, subtitle = subtitle,
        tag = if (branch_id == "five_domains") "a" else "b") +
      theme_gbe() + ggplot2::theme(
        plot.subtitle = ggplot2::element_text(family = GBE_FONT, size = GBE_BASE_SIZE,
                                              face = "bold", margin = ggplot2::margin(b = 6)),
        panel.border = ggplot2::element_rect(fill = NA, colour = "#B3B3B3", linewidth = .35),
        axis.line.y = ggplot2::element_blank(), axis.ticks.y = ggplot2::element_blank(),
        axis.text.y.left = if (show_species) ggplot2::element_text(
          family = GBE_FONT, face = "italic", size = GBE_MIN_TEXT_PT,
          margin = ggplot2::margin(r = 4)) else ggplot2::element_blank(),
        axis.text.y.right = ggplot2::element_text(family = GBE_FONT, size = GBE_MIN_TEXT_PT),
        axis.title.y.right = ggplot2::element_text(family = GBE_FONT, size = GBE_MIN_TEXT_PT),
        plot.tag.position = "topleft", plot.margin = ggplot2::margin(6, 6, 6, 6, "pt"))
    if (include_genome) p <- p + ggplot2::labs(caption = paste(
      "Points: elements; blue bar: median; thin line: 95% element-bootstrap interval.",
      "Gbp*: assembly-anchored estimate; n/a: genome size unavailable.", sep = "\n"))
    p
  }
  standalone <- stats::setNames(lapply(branch_ids, make_panel,
                                      show_species = TRUE, include_genome = TRUE), branch_ids)
  collage_panels <- stats::setNames(lapply(branch_ids, function(branch_id) {
    make_panel(branch_id, show_species = branch_id == "five_domains", include_genome = FALSE)
  }), branch_ids)
  list(plots = standalone, collage_panels = collage_panels, elements = elements,
       summaries = summaries, rows = rows, bandwidth = bandwidth, limits = limits)
}
