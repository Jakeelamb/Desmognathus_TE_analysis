#!/usr/bin/env Rscript
# Build review-only comparisons without promoting new element eligibility.
suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(readr)
})
source("scripts/publication/gbe_theme.R")
source("analyses/01_transposable_elements/ltr_genome_pgls.R")
source("analyses/01_transposable_elements/ltr_filter_comparison.R")
source("scripts/publication/ltr_filter_violin_panels.R")
source("scripts/publication/ltr_filter_pgls_panels.R")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0L) {
  output_root <- "Publication/figures/candidates"
} else if (length(args) == 2L && args[[1]] == "--output-root") {
  output_root <- args[[2]]
} else stop("Usage: build_ltr_filter_comparison.R [--output-root DIRECTORY]")
dataset_dir <- "Publication/datasets"
comparison <- prepare_ltr_filter_comparison(dataset_dir)
violins <- make_ltr_filter_violin_panels(comparison, dataset_dir)
pgls <- make_ltr_filter_pgls_panels(comparison, dataset_dir)

boundary <- paste(
  "Review-only comparison; the four-domain branch is a sensitivity analysis, not a replacement release.",
  "The same >=3,000-bp LTR/Gypsy selection, exact annotation join and >=80% positive-depth",
  "coverage in each of the two LTRs and internal region are required.",
  "Within-species two-sided 1.5-IQR fences are calculated on each branch's eligible raw ratios",
  "before intersecting with coverage. The corrected domain-token count is used throughout.",
  "The broader branch changes the reference-element selection; its structural validity remains to be assessed.",
  "The depth ratio is a mapping/deletion-footprint proxy, not a deletion rate or direct solo-LTR count.",
  "Genome sizes are D. fuscus assembly-anchored IOD-derived estimates, converted at 0.978 Gbp per pg."
)
bootstrap_note <- paste(
  "Primary medians and intervals are preserved exactly from S15.",
  "Four-domain medians are calculated on untransformed ratios, and new 95% percentile intervals",
  "use 2,000 within-species element-bootstrap replicates with deterministic per-species R seeds.",
  "The statistical bootstrap method matches the primary release, while the R random-number stream",
  "differs from the released NumPy PCG64 stream. Seeds and methods are recorded in species_summaries.csv.",
  "Elements are not biological replicates; these intervals describe the recovered-element sample."
)
violin_common <- paste(
  "The current >=5-domain branch contains 380 elements, and the >=4-domain sensitivity contains 989,",
  "across the same 30 species in canonical supplied-tree order.",
  "Both panels have identical 0.025-2.5 ratio axes in log2 space, equal-width trimmed violin outlines",
  "and the same density bandwidth estimated from the original 380 observations.",
  "Hollow points are every retained element; vertical blue bars are species medians and horizontal blue",
  "lines are bootstrap intervals. Vertical jitter does not alter ratios. The dotted line means equal depth.",
  "Counts appear at right. In standalone exports, right-hand labels also give estimated genome size in Gbp.",
  "Broader domain eligibility with recalculated IQR fences adds 611 elements and excludes two previously",
  "retained elements, for a net gain of 609.",
  bootstrap_note, boundary
)
pgls_common <- paste(
  "Matched PGLS panels compare the same 19 species with estimated genome size in both branches,",
  "using the original branch lengths of the shared supplied dated phylogeny.",
  "Each point represents one species median; element counts do not act as regression weights.",
  "Both panels use the same log axes and fit log10(median ratio) ~ log10(genome Gbp),",
  "with maximum-likelihood Pagel lambda constrained to [0,1]. Each branch is refitted independently.",
  "Blue lines are back-transformed geometric fits; ribbons are 95% t coefficient intervals (17 degrees of freedom),",
  "not individual-species prediction intervals. Gray bars are genome-size and element-bootstrap median intervals.",
  "The fitted line and its interval do not propagate measurement error or tree uncertainty.",
  "The supplied tree's publication and calibration provenance remain unresolved.",
  bootstrap_note, boundary
)

write_comparison_candidate <- function(directory, id, title, stem, collage, plots,
                                        height, panel_height, legend, alt, panel_titles) {
  dir.create(file.path(directory, "panels"), recursive = TRUE, showWarnings = FALSE)
  save_gbe_figure(collage, file.path(directory, stem), width_mm = 180, height_mm = height)
  panel_stems <- paste0(stem, c("_a_five_domains", "_b_four_domains"))
  for (i in 1:2) save_gbe_figure(plots[[i]], file.path(directory, "panels", panel_stems[[i]]),
                                width_mm = 180, height_mm = panel_height)
  manifest <- tibble(figure_id = id, title = title, basename = stem,
    width_mm = 180, height_mm = height, final_font_family = resolve_gbe_font(),
    minimum_target_text_pt = GBE_MIN_TEXT_PT, datasets = "S04;S14;S15;trees/desmo900dated_source46.tre",
    release_status = "exploratory_filter_comparison_only", legend = legend, alt_text = alt,
    pdf = paste0(stem, ".pdf"), png = paste0(stem, ".png"), tiff = paste0(stem, ".tif"))
  panels <- tibble(parent_figure = id, panel_id = c("a", "b"), title = panel_titles,
    width_mm = 180, height_mm = panel_height, final_font_family = resolve_gbe_font(),
    minimum_target_text_pt = GBE_MIN_TEXT_PT,
    legend = paste(c("Current branch: >=5 domains.", "Sensitivity branch: >=4 domains."), legend),
    alt_text = paste(panel_titles, alt), pdf = paste0(panel_stems, ".pdf"),
    png = paste0(panel_stems, ".png"), tiff = paste0(panel_stems, ".tif"))
  write_csv(manifest, file.path(directory, "FIGURE_MANIFEST.csv"))
  write_csv(panels, file.path(directory, "panels", "PANEL_MANIFEST.csv"))
  write_lines(c(paste("#", title), "", legend, "", paste("**Alt text:**", alt)),
              file.path(directory, "FIGURE_LEGENDS_AND_ALT_TEXT.md"))
  write_lines(c(paste("#", title, "isolated panels"), "",
    unlist(lapply(1:2, function(i) c(paste("##", panels$panel_id[[i]], panel_titles[[i]]),
      "", panels$legend[[i]], "", paste("**Alt text:**", panels$alt_text[[i]]), "")))),
    file.path(directory, "panels", "PANEL_LEGENDS_AND_ALT_TEXT.md"))
}

violin_dir <- file.path(output_root, "figure_4_filter_comparison")
violin_collage <- (violins$collage_panels[[1]] | violins$collage_panels[[2]]) +
  plot_layout(widths = c(1, 1)) + plot_annotation(caption = paste(
    "Matched species rows, axes and smoothing. Points: retained elements; blue bars: medians and 95% intervals.",
    "Both branches use coverage and within-species IQR filters; the four-domain branch is exploratory.", sep = "\n"
  ), theme = theme(plot.caption = element_text(family = GBE_FONT, size = GBE_MIN_TEXT_PT,
                                               hjust = 0, colour = GBE_COLORS[["gray"]])))
write_comparison_candidate(violin_dir, "Figure_4_Filter", "LTR violin filtering comparison",
  "Figure_4_filter_violin_comparison", violin_collage, violins$plots, 210, 190,
  violin_common, "Two aligned violin panels compare all 30 species under five-domain and four-domain eligibility, with 380 and 989 hollow element points, blue median bars and intervals, identical axes and smoothing, and counts at right.",
  c("Current five-domain violins", "Four-domain sensitivity violins"))

pgls_dir <- file.path(output_root, "figure_4_pgls_comparison")
pgls_collage <- (pgls$plots[[1]] / pgls$plots[[2]]) + plot_layout(heights = c(1, 1))
write_comparison_candidate(pgls_dir, "Figure_4_PGLS_Filter", "PGLS filtering comparison",
  "Figure_4_filter_PGLS_comparison", pgls_collage, pgls$plots, 265, 130,
  pgls_common, "Two matched scatterplots show the same 19 species, estimated genome size and median LTR ratio, with separate Pagel-lambda fits for five-domain and four-domain element selection. Each panel shows all species labels, measurement intervals and a blue confidence ribbon.",
  c("Current five-domain PGLS", "Four-domain sensitivity PGLS"))

analysis_dir <- file.path(violin_dir, "analysis")
dir.create(analysis_dir, showWarnings = FALSE)
write_csv(comparison$elements, file.path(analysis_dir, "retained_elements.csv"))
write_csv(comparison$summaries, file.path(analysis_dir, "species_summaries.csv"))
write_csv(comparison$decisions, file.path(analysis_dir, "element_decisions.csv"))
write_csv(comparison$counts, file.path(analysis_dir, "branch_counts.csv"))
write_csv(comparison$paired_species_summary, file.path(analysis_dir, "paired_species_summary.csv"))
write_csv(violins$rows, file.path(analysis_dir, "species_plot_order.csv"))
write_csv(pgls$model_summary, file.path(analysis_dir, "pgls_model_summary.csv"))
write_csv(pgls$species_data, file.path(analysis_dir, "pgls_species_data.csv"))
write_csv(pgls$predictions, file.path(analysis_dir, "pgls_predictions.csv"))
ape::write.tree(pgls$tree, file = file.path(analysis_dir, "pgls_19_species.tre"), digits = 15)
paths <- c(file.path(dataset_dir, c("Supplementary_Data_S14_ltr_terminal_internal_elements.csv",
  "Supplementary_Data_S15_ltr_terminal_internal_species_robustness.csv",
  "Supplementary_Data_S04_path_analysis_traits24.csv", "trees/desmo900dated_source46.tre")),
  "analyses/01_transposable_elements/ltr_filter_comparison.R",
  "analyses/01_transposable_elements/ltr_genome_pgls.R",
  "scripts/publication/build_ltr_filter_comparison.R",
  "scripts/publication/ltr_filter_violin_panels.R", "scripts/publication/ltr_filter_pgls_panels.R",
  "scripts/publication/gbe_theme.R")
write_csv(tibble(path = paths, md5 = unname(tools::md5sum(paths))),
          file.path(analysis_dir, "input_checksums.csv"))
capture.output(sessionInfo(), file = file.path(analysis_dir, "session_info.txt"))
write_lines(c("# LTR domain-filter visual comparison", "", violin_common, "", pgls_common,
  "", sprintf("Shared violin bandwidth: %.12f log2 units.", violins$bandwidth), "",
  "The analysis tables record both retained sets, every row's four-domain decision, original and new species summaries, seeds, exact membership changes and both PGLS fits. Canonical source data and prior figure exports are unchanged.",
  "No baseline calibration, alignment rerun or structural validation of the expanded cohort is implied.",
  "Rebuild: scripts/run_r.sh Rscript scripts/publication/build_ltr_filter_comparison.R"),
  file.path(analysis_dir, "README.md"))
write_lines(c("# PGLS comparison analysis records", "",
  "Shared data, model fits and provenance are in ../figure_4_filter_comparison/analysis/.", "", pgls_common),
  file.path(pgls_dir, "ANALYSIS_README.md"))
message("Wrote both LTR filtering comparisons under ", output_root)
