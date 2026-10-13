suppressPackageStartupMessages({
  library(ape); library(dplyr); library(tidyr); library(readr)
  library(ggplot2); library(ggrepel); library(patchwork); library(scales)
})
RUN <- normalizePath(getwd(), mustWork = TRUE)
stopifnot(basename(RUN) == "manuscript_story_20260924_01")
DATA <- file.path(RUN, "inputs", "datasets")
source(file.path(RUN, "inputs", "publication_scripts", "gbe_theme.R"))
GBE_BASE_SIZE <- 8.5
GBE_LABEL_SIZE_MM <- 7.5 / ggplot2::.pt
resolve_gbe_font()
read_supp <- function(filename) read_csv(file.path(DATA, filename), show_col_types = FALSE, progress = FALSE)
species_label <- function(x) ifelse(grepl("^D\\. ", x), x, paste("D.", x))
source_tree <- ladderize(read.tree(file.path(DATA, "trees", "desmo900dated_source46.tre")), right = FALSE)
canonical_phylogenetic_order <- source_tree$tip.label[order(node.height(source_tree)[seq_len(Ntip(source_tree))])]
phylogenetic_species_order <- function(species) {
  species <- unique(sub("^D\\. ", "", as.character(species)))
  stopifnot(all(species %in% canonical_phylogenetic_order))
  canonical_phylogenetic_order[canonical_phylogenetic_order %in% species]
}
phylogenetic_species_levels <- function(species) species_label(phylogenetic_species_order(species))
prune_source_tree <- function(species) {
  order <- phylogenetic_species_order(species)
  rotateConstr(keep.tip(source_tree, order), order)
}
te_strategy_levels <- c("aquatic_eggs", "terrestrial_eggs_aquatic_larvae", "direct_development")
te_strategy_colors <- c(aquatic_eggs = GBE_COLORS[["blue"]], terrestrial_eggs_aquatic_larvae = GBE_COLORS[["vermillion"]], direct_development = GBE_COLORS[["purple"]])
te_strategy_shapes <- c(aquatic_eggs = 21, terrestrial_eggs_aquatic_larvae = 22, direct_development = 24)
te_strategy_labels <- c(aquatic_eggs = "Aquatic eggs", terrestrial_eggs_aquatic_larvae = "Terrestrial eggs + aquatic larvae", direct_development = "Direct development")

save_candidate <- function(number, title, plot, panels, width_mm, height_mm, datasets, legend, alt_text) {
  stopifnot(number %in% 1:8, length(panels) >= 1, identical(names(panels), letters[seq_along(panels)]))
  id <- paste0("Figure_", number)
  out <- file.path(RUN, "figures", paste0("figure_", number))
  basename <- paste0(id, "_manuscript_candidate")
  save_gbe_figure(plot, file.path(out, basename), width_mm, height_mm)
  manifest <- tibble(figure_id = id, title = title, basename = basename,
    width_mm = width_mm, height_mm = height_mm, final_font_family = resolve_gbe_font(),
    minimum_target_text_pt = GBE_MIN_TEXT_PT, datasets = datasets,
    release_status = "isolated_manuscript_candidate", legend = legend, alt_text = alt_text,
    pdf = paste0(basename, ".pdf"), png = paste0(basename, ".png"), tiff = paste0(basename, ".tif"))
  write_csv(manifest, file.path(out, "FIGURE_MANIFEST.csv"))
  panel_rows <- lapply(names(panels), function(tag) {
    p <- panels[[tag]]
    panel_basename <- paste0(id, tag, "_isolated")
    panel_plot <- p$plot + labs(tag = tag) + theme(plot.tag.position = "topleft")
    save_gbe_figure(panel_plot, file.path(out, "panels", panel_basename), p$width_mm, p$height_mm)
    tibble(parent_figure = id, panel_id = tag, title = p$title, width_mm = p$width_mm,
      height_mm = p$height_mm, final_font_family = resolve_gbe_font(), minimum_target_text_pt = GBE_MIN_TEXT_PT,
      legend = p$legend, alt_text = p$alt_text,
      pdf = paste0(panel_basename, ".pdf"), png = paste0(panel_basename, ".png"), tiff = paste0(panel_basename, ".tif"))
  })
  panel_manifest <- bind_rows(panel_rows)
  write_csv(panel_manifest, file.path(out, "panels", "PANEL_MANIFEST.csv"))
  write_lines(c(paste0("# ", title), "", legend, "", paste0("Alt text: ", alt_text)), file.path(out, "FIGURE_LEGENDS_AND_ALT_TEXT.md"))
  write_lines(unlist(lapply(seq_len(nrow(panel_manifest)), function(i) c(paste0("# ", panel_manifest$panel_id[i], " ", panel_manifest$title[i]), "", panel_manifest$legend[i], "", paste0("Alt text: ", panel_manifest$alt_text[i]), ""))), file.path(out, "panels", "PANEL_LEGENDS_AND_ALT_TEXT.md"))
  capture.output(sessionInfo(), file = file.path(out, "session_info.txt"))
  invisible(manifest)
}
