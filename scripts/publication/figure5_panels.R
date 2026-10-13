# Isolated panels retain source marks and restore labels needed outside the collage.
figure5_review_panels <- function(panels, tips) {
  if (!identical(names(panels), letters[1:4]) || nrow(tips) != 24L ||
      !all(c("y", "display_species") %in% names(tips)) ||
      !identical(tips$y, as.numeric(1:24)) || anyNA(tips$display_species) ||
      any(!nzchar(tips$display_species)) || anyDuplicated(tips$display_species)) {
    stop("Figure 5 inspection requires panels a-d and the 24 original species rows")
  }
  titles <- c(a = "Phylogeny and development", b = "Estimated genome size",
              c = "Reviewed nucleus area", d = "Reviewed cell area")
  plots <- lapply(names(panels), function(tag) {
    plot <- panels[[tag]]
    if (tag != "a") {
      source_axis <- plot$scales$get_scales("y")
      support_values <- if (tag == "b") {
        c(`Balanced multi-image support` = 21, `One image and animal` = 24,
          `Quality-balance caution` = 4)
      } else {
        c(`Two or more animals` = 21, `One animal` = 24)
      }
      plot <- suppressMessages(plot +
        scale_y_continuous(limits = source_axis$limits, breaks = tips$y,
                           labels = tips$display_species, expand = source_axis$expand) +
        scale_shape_manual(values = support_values, name = NULL,
                           guide = guide_legend(nrow = 1))) +
        theme(axis.text.y = element_text(family = GBE_FONT, face = "italic",
                                        size = GBE_MIN_TEXT_PT),
              legend.position = "bottom")
    }
    plot + labs(tag = tag, title = titles[[tag]]) + theme(
      plot.title = element_text(family = GBE_FONT, size = GBE_BASE_SIZE,
                                face = "bold", hjust = 0.5),
      plot.tag.position = "topleft", plot.margin = margin(8, 12, 8, 8, unit = "pt")
    )
  })
  names(plots) <- names(panels)
  plots
}

export_figure5_panels <- function(panels, tips, output_dir) {
  plots <- figure5_review_panels(panels, tips)
  panel_dir <- file.path(output_dir, "panels")
  titles <- c("Species phylogeny and development", "Estimated genome size",
              "Reviewed nucleus area", "Reviewed cell area")
  basenames <- c("Figure_5a_phylogeny_and_development", "Figure_5b_conditional_genome_size",
                 "Figure_5c_nucleus_area", "Figure_5d_cell_area")
  legends <- c(
    "Topology-only prune of the supplied species tree for all 24 Path24 species in the shared canonical phylogenetic row order. Light gray circles identify 22 species with aquatic larvae and purple triangles identify the two direct developers, D. aeneus and D. wrighti. D. organi is absent from Path24 and is not imputed. Development mode is a two-level trait here, distinct from the three-level reproductive strategy in the TE figures. Branch lengths are omitted because the source tree's publication, calibration, and time provenance remain unresolved.",
    "Estimated genome sizes for all 24 Path24 species, derived from IOD and anchored to the D. fuscus assembly, with the published conditional bootstrap intervals. These estimates are not independently measured species C-values; pg/1C* denotes the assembly anchor and intervals exclude anchor uncertainty. Circles indicate balanced multi-image support, triangles indicate one image and one animal, and crosses flag the prespecified image-quality balance caution. The axis is fixed to 5-30 pg/1C*, preserving every original interval and the shared species row order.",
    "Species nucleus-area medians and the published conditional bootstrap intervals for all 24 Path24 species, based on the 1,152 accepted reviewed cell-nucleus pairs in S18. Circles identify support from two or more animals and triangles identify the seven species represented by one animal. Measurements describe the reviewed upper-tail erythrocyte morphology; intervals are conditional on the observed images and animals, not population-level uncertainty. The axis is fixed to 30-60 square micrometres and rows follow the shared canonical phylogenetic order.",
    "Species cell-area medians and the published conditional bootstrap intervals for all 24 Path24 species, based on the same 1,152 accepted reviewed cell-nucleus pairs in S18. Circles identify support from two or more animals and triangles identify the seven species represented by one animal. Measurements describe the reviewed upper-tail erythrocyte morphology; intervals are conditional on the observed images and animals, not population-level uncertainty. The axis is fixed to 180-340 square micrometres and rows follow the shared canonical phylogenetic order."
  )
  alt_text <- c(
    "A topology-only 24-tip species tree with italic species names, gray aquatic-larva circles and purple direct-development triangles.",
    "Twenty-four labelled species rows show estimated genome sizes and horizontal intervals, with circle, triangle or cross symbols identifying image and animal support.",
    "Twenty-four labelled species rows show nucleus-area medians and horizontal intervals, with triangles identifying single-animal support.",
    "Twenty-four labelled species rows show cell-area medians and horizontal intervals, with triangles identifying single-animal support."
  )
  rows <- lapply(seq_along(plots), function(i) {
    save_gbe_figure(plots[[i]], file.path(panel_dir, basenames[[i]]),
                    width_mm = 180, height_mm = 190, dpi = 300)
    tibble(
      parent_figure = "Figure_5", panel_id = names(plots)[[i]], title = titles[[i]],
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
    "# Figure 5 isolated panels", "",
    "Panels preserve the original estimates, intervals, support symbols and species rows.", "",
    unlist(lapply(seq_len(nrow(manifest)), function(i) c(
      paste0("## ", manifest$panel_id[[i]], " - ", manifest$title[[i]]), "",
      manifest$legend[[i]], "", paste0("**Alt text:** ", manifest$alt_text[[i]]), ""
    )))
  ), file.path(panel_dir, "PANEL_LEGENDS_AND_ALT_TEXT.md"))
  invisible(manifest)
}
