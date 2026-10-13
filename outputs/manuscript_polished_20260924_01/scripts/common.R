suppressPackageStartupMessages({
  library(ape); library(dplyr); library(tidyr); library(readr)
  library(ggplot2); library(ggrepel); library(patchwork); library(scales)
})
RUN <- normalizePath(getwd(), mustWork = TRUE)
stopifnot(basename(RUN) == "manuscript_polished_20260924_01")
SOURCE <- file.path(RUN, "inputs", "source_run")
DATA <- file.path(SOURCE, "inputs", "datasets")
ORIGINAL <- normalizePath(file.path(RUN, "..", "manuscript_story_20260924_01"))
source(file.path(SOURCE, "inputs", "publication_scripts", "gbe_theme.R"))
GBE_BASE_SIZE <- 8.5
GBE_MIN_TEXT_PT <- 7.5
GBE_LABEL_SIZE_MM <- 7.5 / ggplot2::.pt
resolve_gbe_font()
base_theme_gbe <- theme_gbe
theme_gbe <- function(base_size = 8.5, base_family = resolve_gbe_font()) {
  base_theme_gbe(base_size, base_family) + theme(
    plot.title = element_text(size = 9, face = "bold", margin = margin(b = 5)),
    plot.subtitle = element_text(size = 7.5, margin = margin(b = 5)),
    plot.caption = element_text(size = 7.5, hjust = 0, lineheight = 1.05),
    axis.text = element_text(size = 8, colour = "#222222"),
    legend.text = element_text(size = 7.5), legend.title = element_text(size = 7.5),
    strip.text = element_text(size = 8), plot.tag = element_text(size = 11, face = "bold"),
    plot.margin = margin(5, 6, 5, 6))
}
theme_set(theme_gbe())
read_supp <- function(filename) read_csv(file.path(DATA, filename), show_col_types = FALSE, progress = FALSE)
read_analysis <- function(number, filename, folder = "analysis") read_csv(file.path(SOURCE, "figures", paste0("figure_", number), folder, filename), show_col_types = FALSE, progress = FALSE)
species_label <- function(x) ifelse(grepl("^D\\. ", x), x, paste("D.", x))
source_tree <- ladderize(read.tree(file.path(DATA, "trees", "desmo900dated_source46.tre")), right = FALSE)
canonical_phylogenetic_order <- source_tree$tip.label[order(node.height(source_tree)[seq_len(Ntip(source_tree))])]
phylogenetic_species_order <- function(species) {
  species <- unique(sub("^D\\. ", "", as.character(species)))
  stopifnot(all(species %in% canonical_phylogenetic_order))
  canonical_phylogenetic_order[canonical_phylogenetic_order %in% species]
}
phylogenetic_species_levels <- function(species) species_label(phylogenetic_species_order(species))
prune_source_tree <- function(species) {order <- phylogenetic_species_order(species); rotateConstr(keep.tip(source_tree, order), order)}
te_strategy_levels <- c("aquatic_eggs", "terrestrial_eggs_aquatic_larvae", "direct_development")
te_strategy_colors <- c(aquatic_eggs = "#0072B2", terrestrial_eggs_aquatic_larvae = "#D55E00", direct_development = "#6A3D9A", unassigned = "#777777")
te_strategy_shapes <- c(aquatic_eggs = 21, terrestrial_eggs_aquatic_larvae = 22, direct_development = 24, unassigned = 23)
te_strategy_labels <- c(aquatic_eggs = "Aquatic eggs", terrestrial_eggs_aquatic_larvae = "Terrestrial eggs + aquatic larvae", direct_development = "Direct development")
all_strategy_labels <- c(te_strategy_labels, unassigned = "Egg placement unassigned")
strategy_colours <- function(name = NULL) scale_colour_manual(name = name, values = te_strategy_colors, breaks = te_strategy_levels, labels = te_strategy_labels, drop = FALSE)
strategy_shapes <- function(name = NULL) scale_shape_manual(name = name, values = te_strategy_shapes, breaks = te_strategy_levels, labels = te_strategy_labels, drop = FALSE)
format_p <- function(p) ifelse(p < .001, formatC(p, format = "e", digits = 1), sprintf("%.3f", p))

save_candidate <- function(number, title, plot, panels, width_mm = 180, height_mm, datasets, legend, alt_text) {
  stopifnot(number %in% 1:8, height_mm <= 245, width_mm <= 180,
    length(panels) >= 1, identical(names(panels), letters[seq_along(panels)]))
  id <- paste0("Figure_", number)
  out <- file.path(RUN, "figures", paste0("figure_", number))
  dir.create(file.path(out, "panels"), recursive = TRUE, showWarnings = FALSE)
  basename <- paste0(id, "_polished")
  save_gbe_figure(plot, file.path(out, basename), width_mm, height_mm, dpi = 300)
  saveRDS(plot, file.path(out, paste0(basename, ".rds")))
  manifest <- tibble(figure_id = id, title = title, basename = basename,
    width_mm = width_mm, height_mm = height_mm, final_font_family = resolve_gbe_font(),
    minimum_target_text_pt = GBE_MIN_TEXT_PT, datasets = datasets,
    release_status = "polished_review_candidate", legend = legend, alt_text = alt_text,
    pdf = paste0(basename, ".pdf"), png = paste0(basename, ".png"), tiff = paste0(basename, ".tif"))
  write_csv(manifest, file.path(out, "FIGURE_MANIFEST.csv"))
  panel_rows <- lapply(names(panels), function(tag) {
    p <- panels[[tag]]; stem <- paste0(id, tag, "_isolated")
    pp <- p$plot + labs(tag = tag) + theme(plot.tag.position = "topleft")
    save_gbe_figure(pp, file.path(out, "panels", stem), p$width_mm, p$height_mm, dpi = 300)
    tibble(parent_figure = id, panel_id = tag, title = p$title, width_mm = p$width_mm,
      height_mm = p$height_mm, final_font_family = resolve_gbe_font(), minimum_target_text_pt = GBE_MIN_TEXT_PT,
      legend = p$legend, alt_text = p$alt_text,
      pdf = paste0(stem, ".pdf"), png = paste0(stem, ".png"), tiff = paste0(stem, ".tif"))
  })
  pm <- bind_rows(panel_rows)
  write_csv(pm, file.path(out, "panels", "PANEL_MANIFEST.csv"))
  write_lines(c(paste0("# Fig. ", number, ". ", title), "", legend, "", paste0("Alt text: ", alt_text)), file.path(out, "FIGURE_LEGENDS_AND_ALT_TEXT.md"))
  write_lines(unlist(lapply(seq_len(nrow(pm)), function(i) c(paste0("# ", pm$panel_id[i], " ", pm$title[i]), "", pm$legend[i], "", paste0("Alt text: ", pm$alt_text[i]), ""))), file.path(out, "panels", "PANEL_LEGENDS_AND_ALT_TEXT.md"))
  capture.output(sessionInfo(), file = file.path(out, "session_info.txt"))
  invisible(manifest)
}
