# Isolated panels reuse the collage's data, axes, species order, and mappings.
figure3_review_panels <- function(panels) {
  if (!identical(names(panels), letters[1:2])) stop("Figure 3 requires panels a-b")
  result <- lapply(names(panels), function(tag) {
    panels[[tag]] + labs(tag = tag) + theme(
      plot.tag.position = "topleft", plot.margin = margin(8, 12, 8, 8, unit = "pt")
    )
  })
  names(result) <- names(panels)
  result
}

export_figure3_panels <- function(panels, output_dir) {
  plots <- figure3_review_panels(panels)
  panel_dir <- file.path(output_dir, "panels")
  titles <- c("Repeat divergence profiles", "Low-divergence fractions by species")
  basenames <- c("Figure_3a_repeat_divergence_profiles", "Figure_3b_low_divergence_fractions")
  heights <- c(85, 170)
  legends <- c(
    "LTR, LINE, and TIR profiles for all 34 species. Thin gray curves are individual species; coloured curves with distinct line types are equal-weight strategy means for 7 aquatic-egg species, 24 terrestrial-egg species with aquatic larvae, and 3 direct developers. Each facet uses 0-50% sequence divergence and 0-3% of all species RepeatMasker-aligned hit bp per 1-percentage-point bin. Absent order/bin combinations are zero-filled for bins 0-49 without subset renormalization. These alignments are from dnaPipeTE repeat-contig assemblies. Strategy strata are descriptive, not phylogenetically corrected comparisons; divergence is not insertion age and hit-bp fractions are not genome occupancy.",
    "Fractions of all species RepeatMasker-aligned hit bp below 5% divergence (bins 0-4), separately for LTR, LINE, and TIR. The 34 species retain the shared supplied-tree order; all three facets use the same 0-6% axis. Point colours and shapes identify reproductive strategy. Gray stems connect zero to the observed fraction and are not uncertainty intervals. Selected orders and bins are not renormalized. These are observed repeat-contig assay fractions, not genome occupancy, insertion ages, or biological confidence intervals."
  )
  alt_text <- c(
    "Three matched LTR, LINE, and TIR facets with thin gray curves for 34 species and coloured, differently styled mean curves for three reproductive strategies.",
    "Three matched dot-plot columns show 102 low-divergence fractions across 34 italic-labelled species in shared phylogenetic order, using common zero-to-six-percent axes and reproductive-strategy colours and shapes."
  )
  rows <- lapply(seq_along(plots), function(i) {
    save_gbe_figure(plots[[i]], file.path(panel_dir, basenames[[i]]),
                    width_mm = 180, height_mm = heights[[i]], dpi = 300)
    tibble(
      parent_figure = "Figure_3", panel_id = names(plots)[[i]], title = titles[[i]],
      width_mm = 180, height_mm = heights[[i]],
      final_font_family = resolve_gbe_font(), minimum_target_text_pt = GBE_MIN_TEXT_PT,
      legend = legends[[i]], alt_text = alt_text[[i]],
      pdf = paste0(basenames[[i]], ".pdf"), png = paste0(basenames[[i]], ".png"),
      tiff = paste0(basenames[[i]], ".tif")
    )
  })
  manifest <- bind_rows(rows)
  write_csv(manifest, file.path(panel_dir, "PANEL_MANIFEST.csv"))
  write_lines(c(
    "# Figure 3 isolated panels", "",
    "Plots reuse the collage data, fixed scales, colours, species order, and inference boundaries.", "",
    unlist(lapply(seq_len(nrow(manifest)), function(i) c(
      paste0("## ", manifest$panel_id[[i]], " - ", manifest$title[[i]]), "",
      manifest$legend[[i]], "", paste0("**Alt text:** ", manifest$alt_text[[i]]), ""
    )))
  ), file.path(panel_dir, "PANEL_LEGENDS_AND_ALT_TEXT.md"))
  invisible(manifest)
}
