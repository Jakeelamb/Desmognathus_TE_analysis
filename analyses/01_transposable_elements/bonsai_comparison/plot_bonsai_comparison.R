#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(ape))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 1L) stop("Usage: plot_bonsai_comparison.R [comparison-data-dir]")

script_arg <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", script_arg[grepl("^--file=", script_arg)])
comparison_dir <- normalizePath(dirname(script_path), mustWork = TRUE)
analysis_dir <- normalizePath(file.path(comparison_dir, ".."), mustWork = TRUE)
root <- normalizePath(file.path(analysis_dir, "..", ".."), mustWork = TRUE)
data_dir <- if (length(args) == 1L) normalizePath(args[[1]], mustWork = TRUE) else file.path(comparison_dir, "data")
figure_dir <- if (length(args) == 1L) file.path(args[[1]], "figures") else file.path(comparison_dir, "figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

scores <- read.csv(file.path(analysis_dir, "data", "te_pca_scores.csv"), check.names = FALSE)
scores <- scores[scores$te_level == "superfamily" & scores$composition_mode == "classified_conditional", ]
traits <- read.csv(file.path(root, "analyses", "03_phylogenetic_path", "data", "organismal_traits.csv"), stringsAsFactors = FALSE)
fidelity <- read.csv(file.path(data_dir, "te34_bonsai_distance_fidelity_by_species.csv"), stringsAsFactors = FALSE)
tree <- read.tree(file.path(data_dir, "te34_bonsai_tree.nwk"))

if (nrow(scores) != 34L || length(tree$tip.label) != 34L || !setequal(scores$species, tree$tip.label)) {
  stop("Bonsai figure inputs do not match TE34")
}

strategy <- setNames(traits$reproductive_strategy, traits$species)[scores$species]
strategy[is.na(strategy) | strategy == ""] <- "unassigned"
strategy_labels <- c(
  aquatic_eggs = "Aquatic eggs",
  terrestrial_eggs_aquatic_larvae = "Terrestrial eggs + aquatic larvae",
  direct_development = "Direct development",
  unassigned = "Unassigned"
)
strategy_colors <- c(
  aquatic_eggs = "#31688E",
  terrestrial_eggs_aquatic_larvae = "#F8766D",
  direct_development = "#35B779",
  unassigned = "#767676"
)
point_colors <- unname(strategy_colors[strategy])
names(point_colors) <- scores$species

species_label <- function(species) {
  paste0("italic(D.)~italic(", species, ")")
}

draw_figure <- function() {
  layout(matrix(c(1, 2, 3, 4, 4, 4), nrow = 2, byrow = TRUE), widths = c(1.05, 1.45, 0.75), heights = c(1, 0.12))
  par(family = "sans", fg = "#1A1A1A", col.axis = "#1A1A1A", col.lab = "#1A1A1A")

  par(mar = c(4.2, 4.2, 2.2, 0.5))
  plot(
    scores$PC1,
    scores$PC2,
    pch = 21,
    bg = point_colors[scores$species],
    col = "#1A1A1A",
    cex = 1.15,
    xlab = "PC1 score",
    ylab = "PC2 score",
    main = "A  PCA (2-D view)",
    cex.main = 1,
    adj = 0
  )
  text(scores$PC1, scores$PC2, labels = scores$species, pos = 3, cex = 0.47, xpd = NA)
  box(col = "#1A1A1A")

  par(mar = c(4.2, 0.5, 2.2, 1.0))
  tree <- ladderize(tree, right = FALSE)
  tip_colors <- unname(point_colors[tree$tip.label])
  plot(
    tree,
    type = "phylogram",
    use.edge.length = TRUE,
    show.tip.label = TRUE,
    tip.color = tip_colors,
    edge.color = "#4D4D4D",
    edge.width = 0.8,
    cex = 0.48,
    label.offset = max(node.depth.edgelength(tree)) * 0.012,
    no.margin = FALSE
  )
  title(main = "B  Bonsai tree", adj = 0, cex.main = 1)
  axisPhylo(side = 1, backward = FALSE, las = 1, cex.axis = 0.7)
  mtext("Fitted diffusion-time branch length", side = 1, line = 2.6, cex = 0.8)

  par(mar = c(6.0, 3.8, 2.2, 0.8))
  shown <- fidelity[fidelity$representation_id %in% c("pca_2d", "bonsai_tree"), ]
  shown$representation_label <- factor(
    shown$representation_label,
    levels = c("PCA PC1-PC2", "Bonsai tree")
  )
  shown$display_label <- factor(
    ifelse(shown$representation_id == "pca_2d", "PCA 2-D", "Bonsai"),
    levels = c("PCA 2-D", "Bonsai")
  )
  boxplot(
    distance_pearson_r ~ display_label,
    data = shown,
    ylim = c(0.65, 1.01),
    ylab = "Per-species distance Pearson r",
    xlab = "",
    col = c("#D6D6D6", "#00BFC4"),
    border = "#1A1A1A",
    outline = FALSE,
    las = 1,
    cex.axis = 0.75,
    main = "C  Distance fidelity",
    cex.main = 1
  )
  stripchart(
    distance_pearson_r ~ display_label,
    data = shown,
    vertical = TRUE,
    method = "jitter",
    jitter = 0.12,
    pch = 16,
    cex = 0.45,
    col = adjustcolor("#1A1A1A", alpha.f = 0.55),
    add = TRUE
  )

  par(mar = c(0, 0, 0, 0))
  plot.new()
  legend(
    "center",
    legend = unname(strategy_labels),
    pt.bg = unname(strategy_colors[names(strategy_labels)]),
    pch = 21,
    pt.cex = 1.1,
    col = "#1A1A1A",
    bty = "n",
    horiz = TRUE,
    cex = 0.75,
    title = "Expert reproductive strategy (post hoc annotation)"
  )
}

pdf_path <- file.path(figure_dir, "te34_bonsai_vs_pca_exploratory.pdf")
png_path <- file.path(figure_dir, "te34_bonsai_vs_pca_exploratory.png")
cairo_pdf(pdf_path, width = 13.2, height = 8.5, family = "sans")
draw_figure()
dev.off()
png(png_path, width = 3960, height = 2550, res = 300, type = "cairo")
draw_figure()
dev.off()

cat("PASS: wrote exploratory Bonsai comparison figure\n")
