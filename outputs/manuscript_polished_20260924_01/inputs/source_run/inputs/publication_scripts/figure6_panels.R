# Reuse released PGLS fits and their original measurement intervals.
figure6_review_panels <- function(panels) {
  if (!identical(names(panels), letters[1:3])) stop("Figure 6 requires panels a-c")
  titles <- c(a = "Genome size and nucleus area", b = "Genome size and cell area",
              c = "Nucleus area and cell area")
  plots <- lapply(names(panels), function(tag) {
    plot <- panels[[tag]]
    if (tag %in% c("a", "b")) plot <- plot + guides(shape = guide_legend(nrow = 1))
    plot + labs(tag = tag, title = titles[[tag]]) + theme(
      plot.title = element_text(family = GBE_FONT, size = GBE_BASE_SIZE,
                                face = "bold", hjust = 0.5),
      legend.position = "bottom", plot.tag.position = "topleft",
      plot.margin = margin(8, 12, 8, 8, unit = "pt")
    )
  })
  names(plots) <- names(panels)
  plots
}

export_figure6_panels <- function(panels, output_dir) {
  plots <- figure6_review_panels(panels)
  panel_dir <- file.path(output_dir, "panels")
  titles <- c("Genome size and nucleus area", "Genome size and cell area",
              "Nucleus area and cell area")
  basenames <- paste0("Figure_6", letters[1:3], "_", c(
    "genome_size_nucleus_area", "genome_size_cell_area", "nucleus_cell_area"
  ))
  common <- paste(
    "Released 24-species Pagel-lambda PGLS sensitivity fit, using log10 traits followed by sample standardization.",
    "The line and 95% model confidence ribbon are back-transformed onto measurement scales; bars are conditional measurement-bootstrap intervals.",
    "The annotation reports the standardized slope and 95% interval with an unadjusted P-value.",
    "The lambda estimate lies at the lower fitting boundary. Associations do not establish causal direction."
  )
  support <- paste(
    "Estimated genome size is derived from IOD and anchored to the D. fuscus assembly, not an independent species C-value; its intervals exclude anchor uncertainty.",
    "Circles show balanced multiple-image support, triangles show three one-image/one-animal estimates, and crosses flag four species failing the quality-balance gate."
  )
  legends <- c(
    paste("Estimated genome size versus nucleus area.", common, support),
    paste("Estimated genome size versus cell area.", common, support),
    paste("Reviewed upper-tail nucleus area versus cell area.", common,
          "A common circle symbol is used; measurements are conditional on the observed images and animals.")
  )
  alt_text <- paste(titles, "for 24 species, with measurement intervals, a fitted curve and confidence ribbon, and the released standardized slope and P-value.")
  rows <- lapply(seq_along(plots), function(i) {
    save_gbe_figure(plots[[i]], file.path(panel_dir, basenames[[i]]),
                    width_mm = 140, height_mm = 120, dpi = 300)
    tibble(
      parent_figure = "Figure_6", panel_id = names(plots)[[i]], title = titles[[i]],
      width_mm = 140, height_mm = 120,
      final_font_family = resolve_gbe_font(), minimum_target_text_pt = GBE_MIN_TEXT_PT,
      legend = legends[[i]], alt_text = alt_text[[i]],
      pdf = paste0(basenames[[i]], ".pdf"), png = paste0(basenames[[i]], ".png"),
      tiff = paste0(basenames[[i]], ".tif")
    )
  })
  manifest <- bind_rows(rows)
  write_csv(manifest, file.path(panel_dir, "PANEL_MANIFEST.csv"))
  write_lines(c(
    "# Figure 6 isolated panels", "",
    "Panels preserve the released data, PGLS fits, measurement intervals and scales.", "",
    unlist(lapply(seq_len(nrow(manifest)), function(i) c(
      paste0("## ", manifest$panel_id[[i]], " - ", manifest$title[[i]]), "",
      manifest$legend[[i]], "", paste0("**Alt text:** ", manifest$alt_text[[i]]), ""
    )))
  ), file.path(panel_dir, "PANEL_LEGENDS_AND_ALT_TEXT.md"))
  invisible(manifest)
}
