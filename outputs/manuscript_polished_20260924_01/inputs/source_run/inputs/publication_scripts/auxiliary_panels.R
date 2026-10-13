# Inspection exports reuse the auxiliary collage objects and their source values.
auxiliary_label_rows <- function(plot, tips, label_column) {
  if (!all(c("y", label_column) %in% names(tips)) || anyDuplicated(tips$y) ||
      anyNA(tips[[label_column]])) {
    stop("Isolated aligned tracks require unique, labelled source species rows")
  }
  source_axis <- plot$scales$get_scales("y")
  suppressMessages(plot + scale_y_continuous(
    limits = source_axis$limits, breaks = tips$y, labels = tips[[label_column]],
    expand = source_axis$expand
  )) + theme(axis.text.y = element_text(
    family = GBE_FONT, face = "italic", size = GBE_MIN_TEXT_PT
  ))
}

auxiliary_dag_panel <- function(plot, model, title) {
  result <- plot
  result$layers <- lapply(plot$layers, function(layer) {
    if (is.data.frame(layer$data) && "model" %in% names(layer$data)) {
      ggproto(NULL, layer, data = layer$data[layer$data$model == model, , drop = FALSE])
    } else {
      layer
    }
  })
  result + facet_null() + labs(title = title) + theme(
    plot.title = element_text(family = GBE_FONT, face = "bold", size = 9, hjust = 0.5)
  )
}

auxiliary_review_spec <- function(id) {
  if (id == "Auxiliary_Wang_style_integrated_te_history_candidate") {
    plots <- list(
      integrated_tree_panel,
      auxiliary_label_rows(integrated_genome_panel, integrated_tree_tips, "display_species"),
      auxiliary_label_rows(integrated_composition_panel, integrated_tree_tips, "display_species"),
      auxiliary_label_rows(integrated_fraction_panel, integrated_tree_tips, "display_species"),
      auxiliary_label_rows(integrated_profile_panel, integrated_tree_tips, "display_species")
    )
    return(list(number = 8L, plots = plots, height = rep(240, 5), titles = c(
      "All 34 TE species", "Estimated genome size", "Repeat-hit composition",
      "Repeat-order fractions", "Repeat divergence profiles"
    ), slugs = c("representative_phylogeny", "conditional_genome_size", "repeat_hit_composition",
                 "repeat_order_fractions", "repeat_divergence_profiles"), legends = c(
      "Topology-only prune of the supplied species tree to all 34 TE species, retaining canonical species order. Branch lengths are omitted because source publication and calibration provenance remain unresolved.",
      "All 34 TE species retain their phylogenetic row order. Points and intervals show estimated genome sizes in pg/1C for 21 species, derived from IOD and anchored to the D. fuscus assembly; n/a marks the 13 species without estimates. These are not independent species C-values. Conditional measurement intervals omit anchor uncertainty; no missing estimates are imputed.",
      "RepeatMasker aligned-hit-bp composition for all 34 TE species, shown as stacked bars closed across LTR, LINE, TIR, DIRS, SINE, Maverick, and Other. Other includes Unclassified. Shared order colours and phylogenetic rows are retained. These fractions describe the aligned-hit assay, not genome occupancy.",
      "LTR, LINE, and TIR fractions for all 34 species on a common 0-60% scale. Each point and stem shows the order's share of all RepeatMasker-aligned hit bp; the denominator is not restricted to these three orders. Rarer categories remain visible in panel c. These are not genome-occupancy estimates.",
      "Zero-filled LTR, LINE, and TIR RepeatMasker profiles for all 34 species, with 0-50% sequence divergence and a common profile-height scale of 0-3% of species aligned-hit bp. Boxes and sparse rare-order columns are omitted; rarer categories remain in panel c. The full species denominator is retained; divergence is not insertion age, and hit-bp fraction is not genome occupancy."
    )))
  }
  if (id == "Auxiliary_Wang_style_diversity_genome_size_candidate") {
    return(list(number = 9L, plots = list(
      wang_diversity_shannon + labs(x = "Estimated genome size (pg/1C*)"), wang_diversity_gini
    ), height = c(90, 90), titles = c("Shannon diversity and estimated genome size",
                                     "Gini-Simpson diversity and estimated genome size"),
    slugs = c("Shannon_diversity", "Gini_Simpson_diversity"), legends = c(
      "Canonical top-ten-superfamily Shannon entropy against estimated genome size for overlap21, derived from IOD and anchored to the D. fuscus assembly. Columns show all 21 species and post-hoc reproductive strata of 6, 13, and 2 species. Horizontal intervals are conditional measurement intervals. Lines and coefficient confidence bands are Pagel-lambda PGLS fits; the two direct developers have no fitted regression. These estimates are not direct species C-values, and the strata are descriptive.",
      "Canonical top-ten-superfamily Gini-Simpson index against the same estimated genome sizes and intervals. All-species and reproductive-strategy columns retain 21/6/13/2 observations and shared x limits. Pagel-lambda PGLS fits and coefficient confidence bands are shown only where n is at least 3. The n = 2 direct-development panel intentionally has no fit; post-hoc strata do not establish causal effects."
    )))
  }
  if (id == "Auxiliary_Path24_development_traits_reference") {
    plots <- list(path24_tree_reference,
      auxiliary_label_rows(path24_genome_reference, path24_development_traits, "species"),
      auxiliary_label_rows(path24_nucleus_reference, path24_development_traits, "species"),
      auxiliary_label_rows(path24_cell_reference, path24_development_traits, "species"))
    return(list(number = 10L, plots = plots, height = rep(190, 4),
    panel_ids = c("a", "c", "d", "e"), titles = c(
      "Path24 species phylogeny", "Estimated genome size",
      "Reviewed nucleus area", "Reviewed cell area"
    ), slugs = c("species_phylogeny", "conditional_genome_size",
                 "nucleus_area", "cell_area"), legends = c(
      "Topology-only view of the supplied 24-tip Path24 tree, retaining canonical species rows. Source publication and calibration provenance are unresolved, so branch lengths and time axes are omitted. Development mode is described in the caption: 22 species have aquatic larvae and undergo metamorphosis; D. aeneus and D. wrighti are direct developers. D. organi is absent and not imputed. This is descriptive, not a development-mode test.",
      "Estimated genome sizes in pg/1C for Path24, derived from IOD and anchored to the D. fuscus assembly, with released conditional measurement intervals and original phylogenetic row order. Numbers label estimates to one decimal; the 8-27 pg/1C axis retains every interval endpoint. The estimates are not independent species C-values and intervals do not include uncertainty in the anchor.",
      "Reviewed species nucleus-area estimates in square micrometres, with released measurement intervals, for all 24 Path24 species in the original phylogenetic row order. Numbers label estimates to one decimal; the 35-55 square micrometre axis retains every interval endpoint. The plot is descriptive and preserves variation in the underlying measurement support.",
      "Reviewed species cell-area estimates in square micrometres, with released measurement intervals, for all 24 Path24 species in the original phylogenetic row order. Numbers label estimates to one decimal; the 190-335 square micrometre axis retains every interval endpoint. The plot is descriptive and preserves variation in the underlying measurement support."
    )))
  }
  if (id == "Auxiliary_Path_DAG_reference") {
    spec <- representative_audit[order(representative_audit$model_order), ]
    plots <- lapply(seq_len(nrow(spec)), function(i) auxiliary_dag_panel(
      Path_DAG_reference, spec$model[[i]], as.character(spec$panel_label[[i]])
    ))
    return(list(number = 11L, plots = plots, height = rep(110, 10),
      titles = spec$display_name, slugs = spec$model,
      legends = vapply(seq_len(nrow(spec)), function(i) sprintf(
        "%s: one representative of a testable Path24 equivalence class containing %d DAG%s. Estimated genome size is the IOD-derived, D. fuscus assembly-anchored pg/1C estimate, not an independent species C-value. Arrow orientation is not uniquely established by class support. The saturated six-member class has no d-separation claim and remains unscored.",
        spec$display_name[[i]], spec$class_members[[i]], if (spec$class_members[[i]] == 1L) "" else "s"
      ), character(1))))
  }
  if (id == "Auxiliary_Wang2025_cross_study_te_diversity_reference") {
    boundary <- paste(
      "Exact Wang et al. 2025 CC0 Dryad values use RepeatModeler/RepeatMasker; current indices use dnaPipeTE.",
      "Pipelines are not harmonized and no pooled fit or hypothesis test is shown."
    )
    return(list(number = 12L, plots = list(wang_a, wang_b, wang_c, wang_d),
      height = rep(125, 4), titles = c("Shannon diversity and genome size across studies",
        "Shannon diversity across vertebrate groups", "Gini-Simpson diversity and genome size across studies",
        "Gini-Simpson diversity across vertebrate groups"),
      slugs = c("Shannon_scatter", "Shannon_groups", "Gini_Simpson_scatter", "Gini_Simpson_groups"),
      legends = paste(c(
        "Top-ten-superfamily Shannon entropy for 11 published salamanders and 21 current Desmognathus species. Current estimated genome sizes and horizontal intervals are IOD-derived and anchored to the D. fuscus assembly in Gbp, not direct species C-values.",
        "Top-ten-superfamily Shannon distributions for all 84 published vertebrates and 34 current Desmognathus species across nine labelled groups. Boxes and individual points retain source values and group counts; only vertical jitter is applied.",
        "Top-ten-superfamily Gini-Simpson index for 11 published salamanders and 21 current Desmognathus species. Current estimated genome sizes and horizontal intervals are IOD-derived and anchored to the D. fuscus assembly in Gbp, not direct species C-values.",
        "Top-ten-superfamily Gini-Simpson distributions for all 84 published vertebrates and 34 current Desmognathus species across nine labelled groups. Boxes and individual points retain source values and group counts; only vertical jitter is applied."
      ), boundary)))
  }
  stop("Unknown auxiliary inspection figure: ", id)
}

export_auxiliary_panels <- function(id, output_dir) {
  spec <- auxiliary_review_spec(id)
  count <- length(spec$plots)
  expected <- c(`8` = 5L, `9` = 2L, `10` = 4L, `11` = 10L, `12` = 4L)
  if (count != expected[[as.character(spec$number)]] ||
      any(vapply(spec[c("height", "titles", "slugs", "legends")], length, integer(1)) != count)) {
    stop("Auxiliary panel metadata does not match its complete set of plots")
  }
  panel_ids <- if (is.null(spec$panel_ids)) letters[seq_len(count)] else spec$panel_ids
  if (length(panel_ids) != count || anyDuplicated(panel_ids) ||
      any(!panel_ids %in% letters)) {
    stop("Auxiliary panels require unique stable lowercase identifiers")
  }
  panel_dir <- file.path(output_dir, "panels")
  parent <- paste0("Figure_", spec$number)
  manifest <- bind_rows(lapply(seq_len(count), function(i) {
    tag <- letters[[i]]
    plot <- spec$plots[[i]] + labs(tag = tag) + theme(
      legend.position = "bottom", legend.box = "vertical",
      plot.tag = element_text(family = GBE_FONT, face = "bold", size = 9),
      plot.tag.position = "topleft", plot.margin = margin(8, 12, 8, 8, unit = "pt")
    )
    basename <- paste0(parent, panel_ids[[i]], "_", spec$slugs[[i]])
    save_gbe_figure(plot, file.path(panel_dir, basename),
                    width_mm = GBE_DOUBLE_COLUMN_MM, height_mm = spec$height[[i]], dpi = 300)
    entry <- tibble(
      parent_figure = parent, panel_id = panel_ids[[i]], title = spec$titles[[i]],
      width_mm = GBE_DOUBLE_COLUMN_MM, height_mm = spec$height[[i]],
      final_font_family = resolve_gbe_font(), minimum_target_text_pt = GBE_MIN_TEXT_PT,
      legend = spec$legends[[i]], alt_text = paste(spec$titles[[i]], spec$legends[[i]]),
      pdf = paste0(basename, ".pdf"), png = paste0(basename, ".png"), tiff = paste0(basename, ".tif")
    )
    if (!is.null(spec$panel_ids)) {
      entry <- mutate(entry, display_panel_id = tag, .after = panel_id)
    }
    entry
  }))
  write_csv(manifest, file.path(panel_dir, "PANEL_MANIFEST.csv"))
  write_lines(c(paste0("# ", parent, " isolated panels"), "",
    "Inspection exports preserve the collage's data, scales, colours, and species order.", "",
    unlist(lapply(seq_len(nrow(manifest)), function(i) c(
      paste0("## ", letters[[i]], " - ", manifest$title[[i]]), "",
      manifest$legend[[i]], "", paste0("**Alt text:** ", manifest$alt_text[[i]]), ""
    )))
  ), file.path(panel_dir, "PANEL_LEGENDS_AND_ALT_TEXT.md"))
  invisible(manifest)
}
