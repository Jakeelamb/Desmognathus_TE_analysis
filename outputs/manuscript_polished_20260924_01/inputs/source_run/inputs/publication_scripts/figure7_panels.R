# Preserve frozen model comparisons and simulation recovery rates.
figure7_review_panels <- function(panels) {
  if (!identical(names(panels), letters[1:2])) stop("Figure 7 requires panels a-b")
  titles <- c(a = "Path24 model comparison", b = "Simulation calibration")
  plots <- lapply(names(panels), function(tag) {
    panels[[tag]] + labs(tag = tag, title = titles[[tag]]) + theme(
      plot.title = element_text(family = GBE_FONT, size = GBE_BASE_SIZE,
                                face = "bold", hjust = 0.5),
      legend.position = "bottom", plot.tag.position = "topleft",
      plot.margin = margin(8, 12, 8, 8, unit = "pt")
    )
  })
  names(plots) <- names(panels)
  plots
}

export_figure7_panels <- function(panels, output_dir) {
  plots <- figure7_review_panels(panels)
  panel_dir <- file.path(output_dir, "panels")
  titles <- c("Path24 model comparison", "Simulation calibration")
  basenames <- c("Figure_7a_path_model_comparison", "Figure_7b_simulation_calibration")
  legends <- c(
    "All ten released primary Path24 candidate equivalence classes, ranked by delta CICc from the best model. The dashed line marks delta CICc = 2; blue circles pass the global d-separation fit gate and vermillion triangles fail it. The cell-bridge class ranks first. Estimated genome size is the D. fuscus assembly-anchored pg/1C estimate derived from IOD, not an independent species C-value. These exploratory model comparisons do not orient causal arrows within an equivalence class.",
    "True-class top-model recovery in 100 simulations per generating-class and residual-regime combination, retaining all eight released rates and the complete 0-100% scale. Colours and symbols distinguish Brownian and independent residuals. Unequal and incomplete recovery qualifies the path-model comparison. Estimated genome size refers to the IOD-derived assembly-anchored trait; these simulations do not establish unique causal directions."
  )
  alt_text <- c(
    "Ten ranked equivalence classes have delta CICc points; only the best cell-bridge class has a blue global-fit passing symbol.",
    "Paired points show recovery rates for four generating classes under Brownian and independent residuals on a zero-to-100-percent axis."
  )
  heights <- c(135, 110)
  rows <- lapply(seq_along(plots), function(i) {
    save_gbe_figure(plots[[i]], file.path(panel_dir, basenames[[i]]),
                    width_mm = 150, height_mm = heights[[i]], dpi = 300)
    tibble(
      parent_figure = "Figure_7", panel_id = names(plots)[[i]], title = titles[[i]],
      width_mm = 150, height_mm = heights[[i]],
      final_font_family = resolve_gbe_font(), minimum_target_text_pt = GBE_MIN_TEXT_PT,
      legend = legends[[i]], alt_text = alt_text[[i]],
      pdf = paste0(basenames[[i]], ".pdf"), png = paste0(basenames[[i]], ".png"),
      tiff = paste0(basenames[[i]], ".tif")
    )
  })
  manifest <- bind_rows(rows)
  write_csv(manifest, file.path(panel_dir, "PANEL_MANIFEST.csv"))
  write_lines(c(
    "# Figure 7 isolated panels", "",
    "Panels preserve the released model order, global-fit status, simulation rates and scales.", "",
    unlist(lapply(seq_len(nrow(manifest)), function(i) c(
      paste0("## ", manifest$panel_id[[i]], " - ", manifest$title[[i]]), "",
      manifest$legend[[i]], "", paste0("**Alt text:** ", manifest$alt_text[[i]]), ""
    )))
  ), file.path(panel_dir, "PANEL_LEGENDS_AND_ALT_TEXT.md"))
  invisible(manifest)
}
