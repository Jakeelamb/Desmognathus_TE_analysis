# Preserve the shared species rows and raw observations in standalone exports.
figure4_review_panels <- function(panels, tips) {
  if (!identical(names(panels), letters[1:2]) || nrow(tips) != 30L ||
      !all(c("y", "species_label") %in% names(tips)) ||
      !identical(tips$y, as.numeric(1:30)) || anyNA(tips$species_label) ||
      any(!nzchar(tips$species_label)) || anyDuplicated(tips$species_label)) {
    stop("Figure 4 inspection requires panels a-b and the 30 original species rows")
  }
  result <- lapply(names(panels), function(tag) {
    plot <- panels[[tag]]
    if (tag == "b") {
      source_axis <- plot$scales$get_scales("y")
      plot <- suppressMessages(plot + scale_y_continuous(
        limits = source_axis$limits, breaks = source_axis$breaks,
        labels = tips$species_label, expand = source_axis$expand,
        sec.axis = source_axis$secondary.axis
      )) + theme(axis.text.y.left = element_text(
        family = GBE_FONT, face = "italic", size = GBE_MIN_TEXT_PT,
        margin = margin(r = 4, unit = "pt")
      ))
    }
    plot + labs(tag = tag) + theme(
      plot.tag.position = "topleft", plot.margin = margin(8, 12, 8, 8, unit = "pt")
    )
  })
  names(result) <- names(panels)
  result
}

export_figure4_panels <- function(panels, tips, output_dir) {
  plots <- figure4_review_panels(panels, tips)
  panel_dir <- file.path(output_dir, "panels")
  titles <- c("Species phylogeny for the LTR panel", "LTR terminal-to-internal depth")
  basenames <- c("Figure_4a_LTR_species_phylogeny", "Figure_4b_terminal_internal_depth")
  legends <- c(
    "Topology-only 30-tip prune of the supplied 46-tip species tree, retaining the shared canonical phylogenetic row order. These are the 30 species in the primary LTR panel. Branch lengths are omitted because the source tree's publication or archive identifier, calibration method, and branch-length provenance remain unresolved. The tree provides context and ordering, not a new phylogenetic analysis.",
    "Zero-aware terminal:internal LTR mean-depth ratios for the exact 380-element primary branch across 30 species. Light gray horizontal violin outlines show smoothed log2-ratio distributions with a common bandwidth selected from the pooled log2 observations, equal maximum widths, and tails trimmed to each species' observed range. Hollow points show individual elements; short vertical blue bars mark released medians of untransformed ratios, displayed after log2 transformation. Thin horizontal blue intervals are 95% percentile intervals from 2,000 within-species element-bootstrap replicates. Vertical jitter separates overlapping elements without changing their ratios. The right-hand labels pair estimated genome size (Gbp*) with retained element count (3-50), separated by a slash. Genome sizes are released D. fuscus assembly-anchored IOD-derived estimates converted at 0.978 Gbp per pg; the asterisk denotes assembly-anchored estimates, not independently measured species C-values. The n/a entries mark unavailable genome sizes without imputation. The common log2 axis labels ratios directly as decimals; the dotted 1:1 line marks equal terminal and internal mean depth. Violin width does not encode sample size; density shapes are descriptive smoothing and should be read alongside the raw points, especially for species with few recovered elements. Primary elements satisfy within-species two-sided 1.5-IQR fences and at least 80% positive-depth coverage in each terminal and internal region. These intervals are conditional on recovered elements and sequencing resources; elements are not biological replicates. This mapping/deletion-footprint proxy is not an ectopic-recombination, solo-LTR, or DNA-loss rate."
  )
  alt_text <- c(
    "Topology-only species tree with 30 italic species labels in the same canonical order as the depth-ratio panel.",
    "Thirty labelled species rows contain light gray horizontal violins, 380 hollow element observations, short vertical blue median bars and thin horizontal bootstrap intervals, a dotted vertical line at a one-to-one ratio, and paired estimated genome-size and retained-element-count labels at right, with n/a for unavailable genome sizes."
  )
  rows <- lapply(seq_along(plots), function(i) {
    save_gbe_figure(plots[[i]], file.path(panel_dir, basenames[[i]]),
                    width_mm = 180, height_mm = 190, dpi = 300)
    tibble(
      parent_figure = "Figure_4", panel_id = names(plots)[[i]], title = titles[[i]],
      width_mm = 180, height_mm = 190,
      final_font_family = resolve_gbe_font(), minimum_target_text_pt = GBE_MIN_TEXT_PT,
      legend = legends[[i]], alt_text = alt_text[[i]],
      pdf = paste0(basenames[[i]], ".pdf"), png = paste0(basenames[[i]], ".png"),
      tiff = paste0(basenames[[i]], ".tif")
    )
  })
  manifest <- bind_rows(rows)
  write_csv(manifest, file.path(panel_dir, "PANEL_MANIFEST.csv"))
  write_lines(c(
    "# Figure 4 isolated panels", "",
    "Plots preserve all observations, intervals, axes, counts, and canonical species rows.", "",
    unlist(lapply(seq_len(nrow(manifest)), function(i) c(
      paste0("## ", manifest$panel_id[[i]], " - ", manifest$title[[i]]), "",
      manifest$legend[[i]], "", paste0("**Alt text:** ", manifest$alt_text[[i]]), ""
    )))
  ), file.path(panel_dir, "PANEL_LEGENDS_AND_ALT_TEXT.md"))
  invisible(manifest)
}
