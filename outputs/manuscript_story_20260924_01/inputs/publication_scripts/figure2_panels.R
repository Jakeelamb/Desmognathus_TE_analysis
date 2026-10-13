# Reuse the three collage ggplots without refitting or transforming source data.
figure2_review_panels <- function(panels) {
  if (!identical(names(panels), letters[1:3])) stop("Figure 2 requires panels a-c")
  result <- lapply(names(panels), function(tag) {
    panels[[tag]] + labs(tag = tag) + theme(
      plot.tag.position = "topleft", plot.margin = margin(8, 12, 8, 8, unit = "pt")
    )
  })
  names(result) <- names(panels)
  result
}

export_figure2_panels <- function(panels, output_dir) {
  plots <- figure2_review_panels(panels)
  panel_dir <- file.path(output_dir, "panels")
  titles <- c("Bonsai TE-composition topology", "Feature-resampling stability",
              "Reproductive-strategy association")
  basenames <- paste0("Figure_2", letters[1:3], "_", c(
    "TE_composition_topology", "feature_resampling_stability", "strategy_association"
  ))
  widths <- c(180, 140, 140)
  heights <- c(190, 100, 100)
  legends <- c(
    "The 34-tip Bonsai topology fitted to the retry-averaged 24-superfamily CLR matrix, retaining its composition-tree order. This is not a species phylogeny; fitted branch lengths are omitted and the display root has no phylogenetic meaning. Tip symbols show the post-hoc reproductive-strategy annotation: 7 aquatic-egg species, 24 terrestrial-egg species with aquatic larvae, and 3 direct developers. Only the exact aquatic-egg and direct-development edges carry feature-resampling split support labels. Feature resampling perturbs compositional coordinates, not biological specimens.",
    "Exact-edge recovery across 200 bounded feature-resampling refits, on a common 0-100% axis: aquatic eggs 55/200 (27.5%), terrestrial eggs with aquatic larvae 2/200 (1%), and direct development 70/200 (35%). Intervals are 95% Wilson binomial intervals describing the finite number of refits. They are not biological confidence intervals. The relatively low recovery qualifies the stability of exact group-defining edges.",
    "Global localization of the post-hoc reproductive-strategy labels: 2 observed Fitch steps versus a permutation-null mean of 9.41 (99,999 ordinary label permutations; P = 0.00001). The Aitchison-distance association is PERMANOVA R-squared = 0.260, P = 0.00001, with PERMDISP P = 0.707. Ordinary label permutations are not phylogenetically corrected. The result supports strategy-associated compositional structure, not discrete natural clusters, a species phylogeny, or causation."
  )
  alt_text <- c(
    "Topology-only Bonsai tree with 34 italic species labels, shared reproductive-strategy colours and symbols, and two feature-resampling support labels.",
    "Three recovery points and Wilson intervals on a zero-to-100-percent axis; the counts are 55, 2, and 70 successful recoveries out of 200 feature refits.",
    "A filled point at 2 observed Fitch steps and an open point at the null mean of 9.41, with the permutation and multivariate association statistics printed above."
  )
  rows <- lapply(seq_along(plots), function(i) {
    save_gbe_figure(plots[[i]], file.path(panel_dir, basenames[[i]]),
                    width_mm = widths[[i]], height_mm = heights[[i]], dpi = 300)
    tibble(
      parent_figure = "Figure_2", panel_id = names(plots)[[i]], title = titles[[i]],
      width_mm = widths[[i]], height_mm = heights[[i]],
      final_font_family = resolve_gbe_font(), minimum_target_text_pt = GBE_MIN_TEXT_PT,
      legend = legends[[i]], alt_text = alt_text[[i]],
      pdf = paste0(basenames[[i]], ".pdf"), png = paste0(basenames[[i]], ".png"),
      tiff = paste0(basenames[[i]], ".tif")
    )
  })
  manifest <- bind_rows(rows)
  write_csv(manifest, file.path(panel_dir, "PANEL_MANIFEST.csv"))
  write_lines(c(
    "# Figure 2 isolated panels", "",
    "Plots reuse the collage data, scales, colours, tip order, and inference boundaries.", "",
    unlist(lapply(seq_len(nrow(manifest)), function(i) c(
      paste0("## ", manifest$panel_id[[i]], " - ", manifest$title[[i]]), "",
      manifest$legend[[i]], "", paste0("**Alt text:** ", manifest$alt_text[[i]]), ""
    )))
  ), file.path(panel_dir, "PANEL_LEGENDS_AND_ALT_TEXT.md"))
  invisible(manifest)
}
