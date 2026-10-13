#!/usr/bin/env Rscript
# An exploratory companion to Figure 4, built only from released summaries.

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
  library(readr)
})
source("scripts/publication/gbe_theme.R")
source("analyses/01_transposable_elements/ltr_genome_pgls.R")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0L) {
  out <- "Publication/figures/candidates/figure_4_pgls"
} else if (length(args) == 2L && args[[1]] == "--output-dir") {
  out <- args[[2]]
} else {
  stop("Usage: build_ltr_genome_pgls.R [--output-dir DIRECTORY]")
}
dir.create(file.path(out, "analysis"), recursive = TRUE, showWarnings = FALSE)
prepared <- prepare_ltr_genome_data("Publication/datasets")
results <- fit_ltr_genome_pgls(prepared)
data <- prepared$data
data$species_label <- paste("D.", data$species)
primary <- results$summary[results$summary$model == "log_log_lambda", ]
blue <- GBE_COLORS[["blue"]]
lambda_label <- if (primary$lambda < .001) "< 0.001" else sprintf("= %.3f", primary$lambda)
association_description <- if (primary$p_value >= .05) {
  "No detectable association in these 19 species; the interval permits effects in either direction."
} else {
  "The fitted association is exploratory and should be interpreted with the effect interval and sensitivity results."
}

annotation <- sprintf(
  "Pagel-lambda PGLS; n = %d species\nSlope = %.3f (95%% CI %.3f to %.3f)\nP = %.3f; lambda %s",
  primary$n_species, primary$slope, primary$ci_low, primary$ci_high, primary$p_value, lambda_label
)
scatter <- ggplot(data, aes(genome_size_gbp, ratio_median)) +
  geom_ribbon(data = results$prediction,
              aes(x = genome_size_gbp, ymin = ratio_ci_low, ymax = ratio_ci_high),
              inherit.aes = FALSE, fill = blue, alpha = .12) +
  geom_hline(yintercept = 1, linetype = "dotted", linewidth = .35, colour = "#888888") +
  geom_segment(aes(x = genome_size_ci_low_gbp, xend = genome_size_ci_high_gbp,
                   yend = ratio_median), colour = "#BDBDBD", linewidth = .3) +
  geom_segment(aes(y = median_bootstrap_ci_low, yend = median_bootstrap_ci_high,
                   xend = genome_size_gbp), colour = "#BDBDBD", linewidth = .3) +
  geom_line(data = results$prediction, aes(x = genome_size_gbp, y = ratio_fit),
            inherit.aes = FALSE, colour = blue, linewidth = .7) +
  geom_point(shape = 21, fill = "white", colour = GBE_COLORS[["charcoal"]],
             size = 1.9, stroke = .45) +
  geom_text_repel(aes(label = species_label), family = GBE_FONT, fontface = "italic",
                  size = GBE_LABEL_SIZE_MM, seed = 20260923, max.overlaps = Inf,
                  box.padding = .6, point.padding = .25, min.segment.length = 0,
                  segment.colour = "#999999", segment.size = .2,
                  force = 3, force_pull = .1, max.time = 5, max.iter = 30000) +
  annotate("text", x = 8.2, y = 1.7, label = annotation, hjust = 0, vjust = 1,
           family = GBE_FONT, size = GBE_LABEL_SIZE_MM, lineheight = 1.1) +
  scale_x_log10(breaks = c(8, 10, 12, 15, 18, 22), limits = c(7.7, 23),
                expand = expansion(mult = .02)) +
  scale_y_log10(breaks = c(.25, .4, .6, .8, 1, 1.5), limits = c(.23, 1.8),
                labels = scales::label_number(), expand = expansion(mult = .02)) +
  labs(x = "Estimated genome size (Gbp; log scale)",
       y = "Species median LTR terminal-to-internal\ndepth ratio (log scale)", tag = "a",
       caption = paste(
         "Blue line/band: PGLS fit and 95% coefficient interval; dotted line: equal depth.",
         "Gray bars: estimated genome-size intervals and element-bootstrap median intervals.",
         sep = "\n"
       )) +
  theme_gbe() + theme_gbe_framed_facets() +
  theme(plot.tag.position = "topleft", plot.margin = margin(6, 8, 6, 6, "pt"))

# The author omitted the influence plot from review; numerical diagnostics
# remain in analysis/leave_one_species_out.csv. Keep panel a's note identity.
collage <- scatter
fig_id <- "Figure_4_PGLS"
basename <- "Figure_4_PGLS_LTR_ratio_genome_size"
save_gbe_figure(collage, file.path(out, basename), width_mm = 180, height_mm = 130)
panel_dir <- file.path(out, "panels")
panel_names <- "Figure_4_PGLS_a_ratio_genome_size"
plots <- list(scatter)
for (i in seq_along(plots)) {
  save_gbe_figure(plots[[i]], file.path(panel_dir, panel_names[[i]]),
                  width_mm = 180, height_mm = 130)
}

boundary <- paste(
  "Exploratory companion, not a manuscript figure assignment.",
  "Genome sizes are D. fuscus assembly-anchored IOD-derived estimates",
  "converted at 0.978 Gbp per pg, not independent species C-values.",
  "The LTR ratio is a mapping/deletion-footprint proxy, not a deletion rate.",
  "The supplied dated tree's publication and calibration provenance remain unresolved.",
  "Coefficient intervals and P values are conditional on fitted lambda and the supplied tree;",
  "the fit does not propagate measurement uncertainty in either variable."
)
scatter_legend <- paste(
  "Nineteen species have both released genome-size estimates and primary LTR ratios",
  "(239 retained elements); 11 LTR species without genome sizes are excluded without imputation.",
  "Each hollow point is the released median of untransformed element ratios for one species.",
  "Gray bars are estimated genome-size intervals and 95% within-species element-bootstrap median intervals.",
  "Pagel-lambda PGLS fits log10(median ratio) as the response to log10(genome Gbp), by maximum likelihood",
  "with lambda constrained to [0,1], on the original branch lengths of the source tree pruned to 19 species.",
  "Each species contributes one observation; elements are not biological replicates and do not act as regression weights.",
  "The blue line is the back-transformed fitted geometric mean; its band is a 95% t-based coefficient interval,",
  "not an individual-species prediction interval. The dotted ratio-one line denotes equal terminal and internal depth.",
  sprintf("Slope %.6f (95%% CI %.6f to %.6f), P = %.6f; lambda = %.3g%s.",
          primary$slope, primary$ci_low, primary$ci_high, primary$p_value, primary$lambda,
          if (primary$lambda_boundary) " at a parameter boundary" else ""),
  association_description
)
influence_audit <- paste(
  "Leave-one-species-out diagnostics are retained as analysis tables; the author omitted their plot from the review gallery.",
  "Each diagnostic fit omits one species, prunes the tree, and re-estimates lambda with the remaining 18 species.",
  sprintf("Leave-one-out slopes range from %.6f to %.6f; every P value is at least %.6f.",
          min(results$loo$slope), max(results$loo$slope), min(results$loo$p_value)),
  "These deletions assess influence, not independent replications or 19 additional primary hypothesis tests."
)
alt <- "Nineteen labelled species medians and gray measurement intervals on log-scaled axes, with a nearly horizontal blue PGLS fit and broad blue confidence band."
manifest <- tibble(
  figure_id = fig_id, title = "LTR ratio and genome size: exploratory PGLS", basename = basename,
  width_mm = 180, height_mm = 130, final_font_family = resolve_gbe_font(), minimum_target_text_pt = GBE_MIN_TEXT_PT,
  datasets = "S04;S15;trees/desmo900dated_source46.tre", release_status = "exploratory_review_companion",
  legend = paste(scatter_legend, boundary), alt_text = alt,
  pdf = paste0(basename, ".pdf"), png = paste0(basename, ".png"), tiff = paste0(basename, ".tif")
)
panels <- tibble(
  parent_figure = fig_id, panel_id = "a",
  title = "LTR ratio versus estimated genome size",
  width_mm = 180, height_mm = 130, final_font_family = resolve_gbe_font(),
  minimum_target_text_pt = GBE_MIN_TEXT_PT,
  legend = paste(scatter_legend, boundary), alt_text = alt,
  pdf = paste0(panel_names, ".pdf"), png = paste0(panel_names, ".png"), tiff = paste0(panel_names, ".tif")
)
write_csv(manifest, file.path(out, "FIGURE_MANIFEST.csv"))
write_csv(panels, file.path(panel_dir, "PANEL_MANIFEST.csv"))
write_lines(c("# LTR ratio and genome-size PGLS", "", manifest$legend,
              "", paste("**Alt text:**", manifest$alt_text)),
            file.path(out, "FIGURE_LEGENDS_AND_ALT_TEXT.md"))
write_lines(c("# PGLS companion isolated panels", "",
              unlist(lapply(seq_len(nrow(panels)), function(i) c(paste("##", panels$panel_id[[i]], panels$title[[i]]),
                                            "", panels$legend[[i]], "", paste("**Alt text:**", alt[[i]]), "")))),
            file.path(panel_dir, "PANEL_LEGENDS_AND_ALT_TEXT.md"))

analysis_dir <- file.path(out, "analysis")
write_csv(prepared$data, file.path(analysis_dir, "species_data.csv"))
write_csv(prepared$excluded, file.path(analysis_dir, "excluded_species.csv"))
write_csv(results$summary, file.path(analysis_dir, "model_summary.csv"))
write_csv(results$loo, file.path(analysis_dir, "leave_one_species_out.csv"))
write_csv(results$prediction, file.path(analysis_dir, "fitted_relationship.csv"))
ape::write.tree(prepared$tree, file = file.path(analysis_dir, "pgls_19_species.tre"), digits = 15)
input_paths <- c(prepared$input_paths, "analyses/01_transposable_elements/ltr_genome_pgls.R",
                 "scripts/publication/build_ltr_genome_pgls.R", "scripts/publication/gbe_theme.R")
write_csv(tibble(path = input_paths, md5 = unname(tools::md5sum(input_paths))),
          file.path(analysis_dir, "input_checksums.csv"))
capture.output(sessionInfo(), file = file.path(analysis_dir, "session_info.txt"))
write_lines(c(
  "# Exploratory LTR ratio / estimated genome-size association", "",
  scatter_legend, "", influence_audit, "", boundary, "",
  "## Model specification", "",
  "Primary model: log10(ratio_median) ~ log10(genome_size_gbp); phylolm, model=lambda, ML, bounds [0,1].",
  "Slope is a log-log regression coefficient, not a Pearson correlation coefficient. A 1% increase in genome size corresponds approximately to a slope-percent change in fitted ratio.",
  "The medians are calculated on the original ratio scale before taking logarithms. All data are positive; no pseudocount is added.",
  "Standard errors use phylolm's coefficient covariance. Two-sided t tests and 95% t intervals use n-2 residual degrees of freedom, conditional on estimated lambda. The fitted lambda reaches its lower boundary; this is reported rather than suppressed as a successful interior estimate.",
  "No element-count weighting, measurement-error model, causal interpretation, or uncertainty over tree topology/calibration is included.",
  "Model sensitivities (predefined before plotting results): Brownian covariance (lambda=1), ordinary regression, and an untransformed PGLS. They are recorded in model_summary.csv and are not selected by significance.",
  "The 19 leave-one-species-out fits each re-estimate lambda, retain source branch lengths, and use 18 species. All outputs are downstream derivatives; canonical publication tables remain unchanged.", "",
  "## Files", "",
  "species_data.csv: exact included species, original medians, counts, converted genome values, both uncertainty intervals, and transformed variables.",
  "excluded_species.csv: 11 LTR species without a genome estimate, with explicit exclusion reason.",
  "model_summary.csv: primary and sensitivity slopes, t intervals, P values, lambda boundaries, and captured model warnings.",
  "leave_one_species_out.csv: all 19 influence fits; fitted_relationship.csv: plotted geometric fit and confidence ribbon.",
  "pgls_19_species.tre: actual pruned tree with source branch lengths; input_checksums.csv and session_info.txt: provenance.", "",
  "Method reference: [phylolm package documentation](https://search.r-project.org/CRAN/refmans/phylolm/html/phylolm.html); [coefficient inference and prediction implementation](https://github.com/lamho86/phylolm/blob/master/R/phylolm.R)."
), file.path(analysis_dir, "README.md"))
message("Wrote exploratory Figure 4 PGLS companion to ", out)
