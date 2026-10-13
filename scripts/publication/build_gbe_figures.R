#!/usr/bin/env Rscript

# Build the Genome Biology and Evolution ggplot figure set from Publication/datasets.

## ---- figure-setup
parse_figure_options <- function(args = commandArgs(trailingOnly = TRUE)) {
  options <- list(figure = NULL, output_dir = NULL)
  while (length(args) > 0L) {
    option <- args[[1L]]
    if (!option %in% c("--figure", "--output-dir")) {
      stop("Unknown option: ", option, call. = FALSE)
    }
    if (length(args) < 2L || !nzchar(args[[2L]]) || startsWith(args[[2L]], "--")) {
      stop("Missing value for ", option, call. = FALSE)
    }
    key <- if (option == "--figure") "figure" else "output_dir"
    if (!is.null(options[[key]])) {
      stop("Duplicate option: ", option, call. = FALSE)
    }
    options[[key]] <- args[[2L]]
    args <- args[-c(1L, 2L)]
  }
  if (!is.null(options$figure) && !options$figure %in% paste0("Figure_", 1:12)) {
    stop("Unsupported --figure: ", options$figure,
         ". Supported figures: Figure_1 through Figure_12.", call. = FALSE)
  }
  if (!is.null(options$figure) && is.null(options$output_dir)) {
    stop("--figure requires --output-dir to preserve the full figure release.", call. = FALSE)
  }
  options
}

figure_options <- parse_figure_options()
ROOT <- normalizePath(getwd(), mustWork = TRUE)
default_figure_output <- file.path(ROOT, "Publication/figures")
OUT <- if (is.null(figure_options$output_dir)) {
  default_figure_output
} else {
  path.expand(figure_options$output_dir)
}
if (
  !is.null(figure_options$figure) &&
    identical(
      normalizePath(OUT, mustWork = FALSE),
      normalizePath(default_figure_output, mustWork = FALSE)
    )
) {
  stop("Single-figure builds require an output directory separate from Publication/figures.", call. = FALSE)
}
candidate_manifest <- file.path(OUT, "FIGURE_MANIFEST.csv")
if (!is.null(figure_options$figure) && file.exists(candidate_manifest)) {
  existing_candidate <- utils::read.csv(candidate_manifest, stringsAsFactors = FALSE)
  if (!identical(existing_candidate$figure_id, figure_options$figure)) {
    stop("Output directory belongs to a different figure; use a separate candidate directory.", call. = FALSE)
  }
}

suppressPackageStartupMessages({
  library(ape)
  library(dplyr)
  library(ggforce)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(phylolm)
  library(readr)
  library(scales)
  library(tidyr)
})

source(file.path(ROOT, "scripts/publication/gbe_theme.R"))
DATA <- file.path(ROOT, "Publication/datasets")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

read_supp <- function(filename) {
  read_csv(file.path(DATA, filename), show_col_types = FALSE, progress = FALSE)
}

species_label <- function(x) {
  text <- as.character(x)
  ifelse(grepl("^D\\. ", text), text, paste("D.", text))
}

source_tree_path <- file.path(DATA, "trees", "desmo900dated_source46.tre")
source_tree <- ladderize(read.tree(source_tree_path), right = FALSE)
source_tree_tip_y <- node.height(source_tree)[seq_len(Ntip(source_tree))]
canonical_phylogenetic_order <- source_tree$tip.label[order(source_tree_tip_y)]

if (
  Ntip(source_tree) != 46L ||
    Nnode(source_tree) != 45L ||
    anyDuplicated(source_tree$tip.label) ||
    !is.rooted(source_tree) ||
    !is.binary(source_tree) ||
    !is.ultrametric(source_tree, tol = 1e-6) ||
    any(source_tree$edge.length < 0)
) {
  stop("The shared species-axis order requires the exact supplied 46-tip dated tree")
}

phylogenetic_species_order <- function(species) {
  species <- unique(sub("^D\\. ", "", as.character(species)))
  missing_species <- setdiff(species, canonical_phylogenetic_order)
  if (length(missing_species) > 0L) {
    stop(
      "Species are missing from the supplied source tree: ",
      paste(sort(missing_species), collapse = ", ")
    )
  }
  ordered_species <- canonical_phylogenetic_order[
    canonical_phylogenetic_order %in% species
  ]
  if (length(ordered_species) != length(species)) {
    stop("Phylogenetic species ordering lost or duplicated panel members")
  }
  ordered_species
}

phylogenetic_species_levels <- function(species) {
  species_label(phylogenetic_species_order(species))
}

prune_source_tree <- function(species) {
  constraint <- phylogenetic_species_order(species)
  pruned <- keep.tip(source_tree, constraint)
  rotateConstr(pruned, constraint)
}

te_strategy_levels <- c(
  "aquatic_eggs",
  "terrestrial_eggs_aquatic_larvae",
  "direct_development"
)
te_strategy_colors <- c(
  aquatic_eggs = GBE_COLORS[["blue"]],
  terrestrial_eggs_aquatic_larvae = GBE_COLORS[["vermillion"]],
  direct_development = GBE_COLORS[["purple"]]
)
te_strategy_shapes <- c(
  aquatic_eggs = 21,
  terrestrial_eggs_aquatic_larvae = 22,
  direct_development = 24
)
te_strategy_labels <- c(
  aquatic_eggs = "Aquatic eggs",
  terrestrial_eggs_aquatic_larvae = "Terrestrial eggs + aquatic larvae",
  direct_development = "Direct development"
)

manifest_rows <- list()
auxiliary_manifest_rows <- list()

write_numbered_figure_metadata <- function() {
  figure_manifest <- bind_rows(manifest_rows)
  write_csv(figure_manifest, file.path(OUT, "FIGURE_MANIFEST.csv"))
  legend_lines <- c(
    "# Figure legends and alt text",
    "",
    unlist(lapply(seq_len(nrow(figure_manifest)), function(i) {
      row <- figure_manifest[i, ]
      c(
        paste0("## ", row$figure_id),
        "",
        row$legend,
        "",
        paste0("**Alt text:** ", row$alt_text),
        ""
      )
    }))
  )
  write_lines(
    head(legend_lines, -1),
    file.path(OUT, "FIGURE_LEGENDS_AND_ALT_TEXT.md")
  )
  invisible(figure_manifest)
}

record_figure <- function(id, basename, width_mm, height_mm, datasets, legend, alt_text, status) {
  if (!is.null(figure_options$figure) && !identical(id, figure_options$figure)) {
    return(invisible(NULL))
  }
  paths <- save_gbe_figure(
    get(id, envir = .GlobalEnv),
    file.path(OUT, basename),
    width_mm = width_mm,
    height_mm = height_mm,
    dpi = 300
  )
  manifest_rows[[length(manifest_rows) + 1L]] <<- tibble(
    figure_id = id,
    basename = basename,
    width_mm = width_mm,
    height_mm = height_mm,
    final_font_family = resolve_gbe_font(),
    minimum_target_text_pt = GBE_MIN_TEXT_PT,
    datasets = datasets,
    release_status = status,
    legend = legend,
    alt_text = alt_text,
    pdf = basename |> paste0(".pdf"),
    png = basename |> paste0(".png"),
    tiff = basename |> paste0(".tif")
  )
  invisible(paths)
}

record_auxiliary_figure <- function(
  id,
  title,
  plot,
  basename,
  width_mm,
  height_mm,
  datasets,
  boundary,
  legend,
  alt_text,
  status
) {
  review_ids <- c(
    Auxiliary_Wang_style_integrated_te_history_candidate = "Figure_8",
    Auxiliary_Wang_style_diversity_genome_size_candidate = "Figure_9",
    Auxiliary_Path24_development_traits_reference = "Figure_10",
    Auxiliary_Path_DAG_reference = "Figure_11",
    Auxiliary_Wang2025_cross_study_te_diversity_reference = "Figure_12"
  )
  review_id <- unname(review_ids[[id]])
  if (!is.null(figure_options$figure) && !identical(review_id, figure_options$figure)) {
    return(invisible(NULL))
  }
  paths <- save_gbe_figure(
    plot,
    file.path(OUT, basename),
    width_mm = width_mm,
    height_mm = height_mm,
    dpi = 300
  )
  if (!is.null(figure_options$figure)) {
    manifest_rows[[length(manifest_rows) + 1L]] <<- tibble(
      figure_id = review_id, title = title, basename = basename,
      width_mm = width_mm, height_mm = height_mm,
      final_font_family = resolve_gbe_font(), minimum_target_text_pt = GBE_MIN_TEXT_PT,
      datasets = datasets, release_status = status,
      legend = paste(legend, boundary), alt_text = alt_text,
      pdf = paste0(basename, ".pdf"), png = paste0(basename, ".png"),
      tiff = paste0(basename, ".tif")
    )
    write_numbered_figure_metadata()
    source(file.path(ROOT, "scripts/publication/auxiliary_panels.R"))
    export_auxiliary_panels(id, OUT)
    message("Wrote ", review_id, " auxiliary review candidate and its isolated panels to ", OUT)
    quit(save = "no", status = 0, runLast = FALSE)
  }
  auxiliary_manifest_rows[[length(auxiliary_manifest_rows) + 1L]] <<- tibble(
    auxiliary_id = id,
    title = title,
    basename = basename,
    width_mm = width_mm,
    height_mm = height_mm,
    final_font_family = resolve_gbe_font(),
    minimum_target_text_pt = GBE_MIN_TEXT_PT,
    datasets = datasets,
    release_status = status,
    boundary = boundary,
    legend = legend,
    alt_text = alt_text,
    pdf = basename |> paste0(".pdf"),
    png = basename |> paste0(".png"),
    tiff = basename |> paste0(".tif")
  )
  invisible(paths)
}

## ---- figure-1
# Figure 1: Wang-inspired aligned species tree, estimated genome size,
# dnaPipeTE composition, and RepeatMasker divergence profiles.
path_traits <- read_supp("Supplementary_Data_S04_path_analysis_traits24.csv")
tree_audit_edges <- read_supp("Supplementary_Data_S38_phylogeny_edges.csv")
tree_audit_nodes <- read_supp("Supplementary_Data_S39_phylogeny_nodes.csv")
tree_path <- file.path(DATA, "trees", "desmognathus_time_tree_path24_v1.nwk")
path_tree <- read.tree(tree_path)
path_tree <- rotateConstr(
  path_tree,
  phylogenetic_species_order(path_tree$tip.label)
)

tip_count <- Ntip(path_tree)
node_count <- Nnode(path_tree)
root_distance <- node.depth.edgelength(path_tree)
root_age <- max(root_distance[seq_len(tip_count)])
node_age <- root_age - root_distance
node_y <- node.height(path_tree)

tree_horizontal <- tibble(
  parent = path_tree$edge[, 1],
  child = path_tree$edge[, 2]
) %>%
  mutate(
    x = node_age[parent],
    xend = node_age[child],
    y = node_y[child],
    yend = node_y[child]
  )
tree_vertical <- tibble(parent = sort(unique(path_tree$edge[, 1]))) %>%
  rowwise() %>%
  mutate(
    x = node_age[parent],
    ymin = min(node_y[path_tree$edge[path_tree$edge[, 1] == parent, 2]]),
    ymax = max(node_y[path_tree$edge[path_tree$edge[, 1] == parent, 2]])
  ) %>%
  ungroup()
tree_tips <- tibble(
  node = seq_len(tip_count),
  tree_tip = path_tree$tip.label,
  x = 0,
  y = node_y[seq_len(tip_count)]
) %>%
  left_join(
    path_traits %>% select(species, tree_tip),
    by = "tree_tip",
    relationship = "one-to-one"
  )

audit_tip_labels <- tree_audit_nodes %>%
  filter(is_tip) %>%
  pull(tip_label)
if (
  tip_count != 24L ||
    node_count != 23L ||
    nrow(tree_horizontal) != 46L ||
    nrow(tree_vertical) != 23L ||
    nrow(tree_audit_edges) != 46L ||
    nrow(tree_audit_nodes) != 47L ||
    any(is.na(tree_tips$species)) ||
    !setequal(path_tree$tip.label, path_traits$tree_tip) ||
    !setequal(path_tree$tip.label, audit_tip_labels) ||
    !isTRUE(all.equal(
      tree_horizontal$x - tree_horizontal$xend,
      path_tree$edge.length,
      tolerance = 1e-12,
      check.attributes = FALSE
    )) ||
    abs(root_age - max(tree_audit_nodes$root_distance_myr)) > 1e-4
) {
  stop("Path24 figure references require the exact audited 24-tip time tree")
}
if (
  !is.rooted(path_tree) ||
    !is.binary(path_tree) ||
    !is.ultrametric(path_tree, tol = 1e-6) ||
    any(path_tree$edge.length < 0)
) {
  stop("Path24 figure references require a rooted, binary, time-scaled tree with nonnegative edges")
}

phylogeny_a <- ggplot() +
  geom_segment(
    data = tree_horizontal,
    aes(x = x, xend = xend, y = y, yend = yend),
    linewidth = 0.35,
    lineend = "square",
    colour = GBE_COLORS[["charcoal"]]
  ) +
  geom_segment(
    data = tree_vertical,
    aes(x = x, xend = x, y = ymin, yend = ymax),
    linewidth = 0.35,
    lineend = "square",
    colour = GBE_COLORS[["charcoal"]]
  ) +
  geom_text(
    data = tree_tips,
    aes(x = -0.6, y = y, label = species),
    hjust = 0,
    size = GBE_LABEL_SIZE_MM,
    fontface = "italic"
  ) +
  scale_x_reverse(
    breaks = c(30, 20, 10, 0),
    limits = c(root_age, -14),
    expand = c(0, 0)
  ) +
  scale_y_continuous(limits = c(0.5, tip_count + 0.5), expand = c(0, 0)) +
  labs(x = "Time before present (Ma; as supplied)", y = NULL) +
  theme_gbe(base_size = 7.5) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    plot.margin = margin(3, 3, 3, 4, unit = "pt")
  )
overview_annotations <- read_supp(
  "Supplementary_Data_S42_te_bonsai_tip_annotations.csv"
)
overview_quality <- read_supp(
  "Supplementary_Data_S24_nuclear_iod_quality_balance.csv"
)
overview_traits <- path_traits %>%
  transmute(
    species = sub("^D\\. ", "", species),
    genome_size_pg_fuscus_anchored,
    genome_size_pg_ci_low,
    genome_size_pg_ci_high,
    n_genome_size_images,
    n_genome_size_specimens
  ) %>%
  inner_join(
    overview_annotations %>%
      select(species, reproductive_strategy, reproductive_strategy_label),
    by = "species",
    relationship = "one-to-one"
  ) %>%
  left_join(
    overview_quality %>%
      transmute(
        species = sub("^D\\. ", "", species),
        passes_prespecified_balance_thresholds
      ),
    by = "species",
    relationship = "one-to-one"
  )
if (
  nrow(overview_traits) != 21L ||
    !identical(
      as.integer(table(factor(
        overview_traits$reproductive_strategy,
        levels = te_strategy_levels
      ))),
      c(6L, 13L, 2L)
    ) ||
    any(is.na(overview_traits$passes_prespecified_balance_thresholds))
) {
  stop("Figure 1 requires the exact 21-species, 6/13/2 integrated panel")
}

overview_tree <- prune_source_tree(overview_traits$species)
overview_tip_count <- Ntip(overview_tree)
overview_node_y <- node.height(overview_tree)
overview_edges <- reorder.phylo(overview_tree, "cladewise")$edge
overview_depth <- numeric(overview_tip_count + Nnode(overview_tree))
overview_root <- setdiff(overview_edges[, 1], overview_edges[, 2])
if (length(overview_root) != 1L) stop("Figure 1 overview tree requires one root")
for (edge_index in seq_len(nrow(overview_edges))) {
  parent <- overview_edges[edge_index, 1]
  child <- overview_edges[edge_index, 2]
  overview_depth[child] <- overview_depth[parent] + 1
}
overview_tip_x <- max(overview_depth) + 0.5
overview_tree_horizontal <- tibble(
  parent = overview_tree$edge[, 1],
  child = overview_tree$edge[, 2]
) %>%
  mutate(
    x = overview_depth[parent],
    xend = if_else(child <= overview_tip_count, overview_tip_x, overview_depth[child]),
    y = overview_node_y[child],
    yend = overview_node_y[child]
  )
overview_tree_vertical <- tibble(
  parent = sort(unique(overview_tree$edge[, 1]))
) %>%
  rowwise() %>%
  mutate(
    x = overview_depth[parent],
    ymin = min(overview_node_y[
      overview_tree$edge[overview_tree$edge[, 1] == parent, 2]
    ]),
    ymax = max(overview_node_y[
      overview_tree$edge[overview_tree$edge[, 1] == parent, 2]
    ])
  ) %>%
  ungroup()
overview_tips <- tibble(
  species = overview_tree$tip.label,
  x = overview_tip_x,
  y = overview_node_y[seq_len(overview_tip_count)]
) %>%
  left_join(overview_traits, by = "species", relationship = "one-to-one") %>%
  mutate(
    reproductive_strategy = factor(
      reproductive_strategy,
      levels = te_strategy_levels
    ),
    genome_support = factor(
      case_when(
        !passes_prespecified_balance_thresholds ~ "Quality-balance caution",
        n_genome_size_images == 1L & n_genome_size_specimens == 1L ~
          "One image and animal",
        TRUE ~ "Balanced multi-image support"
      ),
      levels = c(
        "Balanced multi-image support",
        "One image and animal",
        "Quality-balance caution"
      )
    ),
    display_species = species_label(species)
  )
overview_order <- phylogenetic_species_order(overview_traits$species)
if (
  overview_tip_count != 21L ||
    !identical(overview_tree$tip.label[order(overview_tips$y)], overview_order) ||
    any(is.na(overview_tips$genome_size_pg_fuscus_anchored))
) {
  stop("Figure 1 species tracks must share the exact overlap21 tree-tip order")
}

overview_order_composition <- read_csv(
  file.path(
    ROOT,
    "analyses",
    "01_transposable_elements",
    "data",
    "te_order_composition.csv"
  ),
  show_col_types = FALSE,
  progress = FALSE
) %>%
  filter(species %in% overview_traits$species) %>%
  pivot_longer(-species, names_to = "source_order", values_to = "fraction") %>%
  mutate(
    repeat_order = if_else(
      source_order %in% c("LTR", "LINE", "TIR"),
      source_order,
      "Other"
    )
  ) %>%
  group_by(species, repeat_order) %>%
  summarise(fraction = sum(fraction), .groups = "drop") %>%
  mutate(
    repeat_order = factor(repeat_order, levels = c("LTR", "LINE", "TIR", "Other"))
  ) %>%
  arrange(species, repeat_order) %>%
  group_by(species) %>%
  mutate(xmax = cumsum(fraction), xmin = lag(xmax, default = 0)) %>%
  ungroup() %>%
  left_join(overview_tips %>% select(species, y), by = "species")
if (
  nrow(overview_order_composition) != 21L * 4L ||
    any(abs(
      overview_order_composition %>%
        group_by(species) %>%
        summarise(total = sum(fraction), .groups = "drop") %>%
        pull(total) - 1
    ) > 1e-10)
) {
  stop("Figure 1 broad-order compositions must close to one for overlap21")
}
overview_order_colors <- GBE_TE_ORDER_COLORS

overview_superfamilies <- read_csv(
  file.path(
    ROOT,
    "analyses",
    "01_transposable_elements",
    "data",
    "te_superfamily_composition.csv"
  ),
  show_col_types = FALSE,
  progress = FALSE
) %>%
  filter(species %in% overview_traits$species) %>%
  pivot_longer(-species, names_to = "superfamily", values_to = "fraction")
overview_top_superfamilies <- overview_superfamilies %>%
  group_by(superfamily) %>%
  summarise(total_fraction = sum(fraction), .groups = "drop") %>%
  slice_max(total_fraction, n = 9, with_ties = FALSE) %>%
  arrange(desc(total_fraction)) %>%
  pull(superfamily)
expected_overview_superfamilies <- c(
  "Gypsy",
  "hAT",
  "L2",
  "DIRS",
  "L1",
  "PIF-Harbinger",
  "Mutator",
  "Tc1-mariner",
  "CACTA"
)
if (!setequal(overview_top_superfamilies, expected_overview_superfamilies)) {
  stop("Figure 1 top-nine superfamily set drifted")
}
overview_superfamily_bars <- overview_superfamilies %>%
  filter(superfamily %in% overview_top_superfamilies) %>%
  mutate(
    superfamily = factor(superfamily, levels = overview_top_superfamilies),
    superfamily_index = as.integer(superfamily)
  ) %>%
  left_join(overview_tips %>% select(species, y), by = "species")
overview_superfamily_max <- 0.70
if (any(overview_superfamily_bars$fraction > overview_superfamily_max)) {
  stop("Figure 1 superfamily values exceed the shared 0-70% display scale")
}
overview_superfamily_bars <- overview_superfamily_bars %>%
  mutate(
    xmin = superfamily_index - 0.38,
    xmax = superfamily_index + 0.38,
    ymin = y - 0.38,
    ymax = ymin + 0.76 * fraction / overview_superfamily_max
  )
overview_superfamily_frames <- overview_tips %>%
  mutate(
    xmin = 0.5,
    xmax = 9.5,
    ymin = y - 0.39,
    ymax = y + 0.39
  )

# Each profile uses the same 0-50% divergence and 0-3% hit-bp scales.
# Missing category/bin combinations are true zeroes within the released assay.
overview_landscape <- read_supp(
  "Supplementary_Data_S12_repeatmasker_divergence_landscape.csv"
) %>%
  filter(
    species %in% overview_traits$species,
    order %in% c("LTR", "LINE", "TIR")
  ) %>%
  group_by(species, order, divergence_bin_start_pct) %>%
  summarise(percent_species_hit_bp = sum(percent_species_hit_bp), .groups = "drop")
if (
  any(overview_landscape$divergence_bin_start_pct < 0 |
        overview_landscape$divergence_bin_start_pct >= 50) ||
    any(overview_landscape$percent_species_hit_bp > 3)
) {
  stop("Figure 1 profiles exceed the fixed divergence or abundance scales")
}
overview_landscape <- overview_landscape %>%
  complete(
    species = overview_traits$species,
    order = c("LTR", "LINE", "TIR"),
    divergence_bin_start_pct = 0:49,
    fill = list(percent_species_hit_bp = 0)
  ) %>%
  mutate(
    order = factor(order, levels = c("LTR", "LINE", "TIR")),
    order_index = as.integer(order)
  ) %>%
  left_join(overview_tips %>% select(species, y), by = "species") %>%
  mutate(
    xmin = order_index - 0.38 + 0.76 * divergence_bin_start_pct / 50,
    xmax = order_index - 0.38 + 0.76 * (divergence_bin_start_pct + 1) / 50,
    ymin = y - 0.38,
    ymax = ymin + 0.76 * percent_species_hit_bp / 3
  )
if (
  nrow(overview_landscape) != 21L * 3L * 50L ||
    any(is.na(overview_landscape$y)) ||
    any(!is.finite(overview_landscape$percent_species_hit_bp))
) {
  stop("Figure 1 profiles require exactly 21 species by three orders by 50 bins")
}
overview_profile_frames <- overview_landscape %>%
  distinct(species, order, order_index, y)

overview_tree_panel <- ggplot() +
  geom_segment(
    data = overview_tree_horizontal,
    aes(x = x / overview_tip_x, xend = xend / overview_tip_x, y = y, yend = yend),
    linewidth = 0.25,
    lineend = "square",
    colour = GBE_COLORS[["charcoal"]]
  ) +
  geom_segment(
    data = overview_tree_vertical,
    aes(x = x / overview_tip_x, xend = x / overview_tip_x, y = ymin, yend = ymax),
    linewidth = 0.25,
    lineend = "square",
    colour = GBE_COLORS[["charcoal"]]
  ) +
  geom_point(
    data = overview_tips,
    aes(x = 1, y = y, fill = reproductive_strategy, shape = reproductive_strategy),
    size = 1.75,
    stroke = 0.45,
    colour = GBE_COLORS[["charcoal"]]
  ) +
  geom_text(
    data = overview_tips,
    aes(x = 1.15, y = y, label = display_species),
    hjust = 0,
    size = GBE_LABEL_SIZE_MM,
    family = resolve_gbe_font(),
    fontface = "italic"
  ) +
  scale_fill_manual(
    values = te_strategy_colors,
    breaks = te_strategy_levels,
    labels = te_strategy_labels,
    name = "Reproduction"
  ) +
  scale_shape_manual(
    values = te_strategy_shapes,
    breaks = te_strategy_levels,
    labels = te_strategy_labels,
    name = "Reproduction"
  ) +
  scale_x_continuous(
    limits = c(-0.03, 2.8),
    breaks = NULL,
    expand = c(0, 0)
  ) +
  scale_y_continuous(limits = c(0.5, 21.5), expand = c(0, 0)) +
  labs(title = "Phylogeny +\nreproductive strategy", x = NULL, y = NULL) +
  coord_cartesian(clip = "off") +
  theme_gbe() +
  theme(
    axis.line = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    plot.title = element_text(size = 8, face = "bold", hjust = 0.5),
    plot.margin = margin(3, 3, 3, 3, unit = "pt")
  )

overview_genome_panel <- ggplot(
  overview_tips,
  aes(genome_size_pg_fuscus_anchored, y)
) +
  geom_segment(
    aes(x = genome_size_pg_ci_low, xend = genome_size_pg_ci_high, yend = y),
    linewidth = 0.45,
    colour = GBE_COLORS[["gray"]]
  ) +
  geom_point(
    shape = 21,
    size = 1.9,
    stroke = 0.55,
    fill = "white",
    colour = GBE_COLORS[["charcoal"]]
  ) +
  geom_text(
    aes(y = y + 0.38, label = sprintf("%.1f", genome_size_pg_fuscus_anchored)),
    family = GBE_FONT, size = GBE_LABEL_SIZE_MM, vjust = 0.5,
    colour = GBE_COLORS[["charcoal"]]
  ) +
  scale_x_continuous(limits = c(8, 27), breaks = c(10, 18, 26),
                     expand = expansion(mult = 0.08)) +
  scale_y_continuous(limits = c(0.5, 21.5), breaks = NULL, expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  labs(
    title = "Estimated\ngenome size",
    x = "pg/1C*",
    y = NULL
  ) +
  theme_gbe() +
  theme(
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    plot.title = element_text(size = 8, face = "bold", hjust = 0.5),
    plot.margin = margin(3, 7, 3, 3, unit = "pt")
  )

overview_order_panel <- ggplot(overview_order_composition) +
  geom_rect(
    aes(
      xmin = xmin,
      xmax = xmax,
      ymin = y - 0.36,
      ymax = y + 0.36,
      fill = repeat_order
    ),
    colour = NA
  ) +
  scale_fill_manual(values = overview_order_colors, name = "TE order") +
  scale_x_continuous(
    limits = c(0, 1),
    # Cumulative fractions can exceed one by floating-point roundoff; keep
    # the final segment instead of censoring its endpoint to NA.
    oob = scales::squish,
    breaks = c(0, 0.5, 1),
    labels = percent_format(accuracy = 1),
    expand = c(0, 0)
  ) +
  scale_y_continuous(limits = c(0.5, 21.5), breaks = NULL, expand = c(0, 0)) +
  labs(title = "Broad TE\ncomposition", x = "Classified\nbase share", y = NULL) +
  theme_gbe() +
  theme(
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    plot.title = element_text(size = 8, face = "bold", hjust = 0.5),
    plot.margin = margin(3, 6, 3, 6, unit = "pt")
  )

overview_superfamily_panel <- ggplot() +
  geom_rect(
    data = overview_superfamily_frames,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    fill = "white",
    colour = GBE_COLORS[["light_gray"]],
    linewidth = 0.22
  ) +
  geom_rect(
    data = overview_superfamily_bars,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = superfamily),
    colour = NA
  ) +
  scale_fill_manual(
    values = GBE_TE_SUPERFAMILY_COLORS,
    breaks = overview_top_superfamilies,
    labels = paste(seq_along(overview_top_superfamilies), overview_top_superfamilies),
    name = "TE superfamily",
    guide = guide_legend(nrow = 3, byrow = TRUE)
  ) +
  scale_x_continuous(
    limits = c(0.5, 9.5),
    breaks = 1:9,
    labels = seq_along(overview_top_superfamilies),
    expand = c(0, 0),
    position = "top"
  ) +
  scale_y_continuous(limits = c(0.5, 21.5), breaks = NULL, expand = c(0, 0)) +
  labs(title = "TE superfamilies\n0–70% classified bases", x = NULL, y = NULL) +
  theme_gbe() +
  theme(
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text.x.top = element_text(size = GBE_MIN_TEXT_PT),
    plot.title = element_text(size = 8, face = "bold", hjust = 0.5)
  )

overview_landscape_panel <- ggplot() +
  geom_rect(
    data = overview_profile_frames,
    aes(xmin = order_index - 0.38, xmax = order_index + 0.38,
        ymin = y - 0.39, ymax = y + 0.39),
    fill = "white", colour = GBE_COLORS[["light_gray"]], linewidth = 0.22
  ) +
  geom_rect(
    data = overview_landscape,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = order),
    colour = NA
  ) +
  scale_fill_manual(values = overview_order_colors, guide = "none") +
  scale_x_continuous(
    limits = c(0.5, 3.5),
    breaks = as.vector(outer(c(-0.38, 0, 0.38), 1:3, "+")),
    labels = rep(c("0", "25", "50"), 3),
    expand = c(0, 0),
    sec.axis = dup_axis(name = NULL, breaks = 1:3, labels = c("LTR", "LINE", "TIR"))
  ) +
  scale_y_continuous(limits = c(0.5, 21.5), breaks = NULL, expand = c(0, 0)) +
  labs(title = "Divergence profiles\n0–3% aligned hit bp / bin", x = "Sequence divergence (%)", y = NULL) +
  theme_gbe() +
  theme(
    axis.line = element_blank(),
    axis.ticks.y = element_blank(),
    axis.ticks.x.top = element_blank(),
    axis.text.x = element_text(size = GBE_MIN_TEXT_PT),
    axis.text.x.top = element_text(size = GBE_MIN_TEXT_PT),
    plot.title = element_text(size = 8, face = "bold", hjust = 0.5)
  )

Figure_1 <-
  overview_tree_panel +
  overview_genome_panel +
  overview_order_panel +
  overview_superfamily_panel +
  overview_landscape_panel +
  plot_layout(
    widths = c(2.5, 1.1, 1.3, 2.3, 3.3),
    guides = "collect"
  ) +
  plot_annotation(tag_levels = "a") &
  theme(
    legend.position = "bottom",
    legend.box = "vertical",
    plot.tag.position = "topleft"
  )

record_figure(
  "Figure_1", "Figure_1_integrated_TE_landscape", GBE_DOUBLE_COLUMN_MM, 190,
  "S04;S12;S24;S42;te_order_composition;te_superfamily_composition;trees/desmo900dated_source46.tre",
  "TE composition and divergence profiles vary among the 21 species with both genomic and cellular-trait data. Every panel shares one fixed phylogenetic row order, adapting the aligned-track layout of Wang et al. (2025; doi:10.1038/s42003-025-08127-3). Panel a is a topology-only prune of the supplied species tree; branch lengths are omitted because the tree's publication, calibration, and time provenance remain unresolved. Panel b shows D. fuscus assembly-anchored IOD-derived pg/1C estimates and conditional 95% measurement-bootstrap intervals. Uniform circles and one-decimal labels show each estimate on an 8-27 pg/1C axis. Image-quality and specimen-support cautions remain in S24; intervals are conditional on the observed support. These are not independently measured species C-values, and intervals exclude anchor uncertainty. Panel c shows broad-order shares of classified dnaPipeTE aligned bases. Panel d shows the nine most abundant superfamilies from the same assay: every miniature bar uses the common 0–70% height scale, and numbered positions match the legend. Panel e shows LTR, LINE, and TIR RepeatMasker divergence profiles in 1-percentage-point bins from 0 to 50%; all 63 profiles use the same 0–3% vertical scale, measured as each bin's share of all RepeatMasker-aligned hit bp for that species. Absent order/bin combinations are zero-filled without renormalizing the displayed subset. dnaPipeTE composition and RepeatMasker landscape denominators are distinct; neither is genome occupancy, and divergence is not insertion age.",
  "Twenty-one species share aligned rows across a topology-only tree, reproductive-strategy symbols, estimated genome-size intervals, broad TE composition, nine numbered superfamily bars, and three divergence profiles. All superfamily bars share a 0–70% scale; all profiles share 0–50% divergence and 0–3% aligned-hit-bp scales. Uniform genome-size circles carry numeric estimate labels; support cautions remain documented in S24.",
  "review_candidate_wang_aligned_profiles_overlap21"
)

if (identical(figure_options$figure, "Figure_1")) {
  figure_manifest <- write_numbered_figure_metadata()
  source(file.path(ROOT, "scripts/publication/figure1_panels.R"))
  export_figure1_panels(
    list(a = overview_tree_panel, b = overview_genome_panel, c = overview_order_panel,
         d = overview_superfamily_panel, e = overview_landscape_panel),
    overview_tips, OUT, figure_manifest[1, ]
  )
  message("Wrote ", nrow(figure_manifest), " numbered GBE figure to ", OUT)
  quit(save = "no", status = 0, runLast = FALSE)
}

## ---- figure-2
# Figure 2: TE diversity and the Bonsai TE-composition representation.
te_diversity <- read_supp("Supplementary_Data_S08_te_diversity.csv") %>%
  select(
    species,
    n_superfamilies_included,
    shannon_entropy,
    gini_simpson
  ) %>%
  pivot_longer(
    cols = c(shannon_entropy, gini_simpson),
    names_to = "metric",
    values_to = "diversity_value"
  ) %>%
  mutate(
    metric = factor(
      recode(
        metric,
        shannon_entropy = "Shannon entropy",
        gini_simpson = "Gini-Simpson index"
      ),
      levels = c("Shannon entropy", "Gini-Simpson index")
    )
  )
diversity_counts <- te_diversity %>% count(metric)
if (
  nrow(te_diversity) != 68L ||
  n_distinct(te_diversity$species) != 34L ||
    any(te_diversity$n_superfamilies_included != 10L) ||
    nrow(diversity_counts) != 2L ||
    !all(diversity_counts$n == 34L)
) {
  stop("The TE34 top-10-superfamily diversity audit requires 34 species per metric")
}
te_annotations <- read_supp("Supplementary_Data_S42_te_bonsai_tip_annotations.csv")
te_validation <- read_supp("Supplementary_Data_S43_te_bonsai_strategy_validation.csv")
bonsai_tree_path <- file.path(
  DATA,
  "trees",
  "te34_bonsai_composition_tree_with_feature_support.nwk"
)
bonsai_tree <- ladderize(read.tree(bonsai_tree_path), right = FALSE)

if (
  Ntip(bonsai_tree) != 34L ||
    Nnode(bonsai_tree) != 33L ||
    !is.rooted(bonsai_tree) ||
    !is.binary(bonsai_tree) ||
    any(bonsai_tree$edge.length < 0) ||
    nrow(te_annotations) != 34L ||
    !setequal(bonsai_tree$tip.label, te_annotations$species) ||
    !identical(
      unname(sort(as.integer(table(te_annotations$reproductive_strategy)))),
      sort(c(7L, 24L, 3L))
    ) ||
    nrow(te_validation) != 3L ||
    any(te_validation$feature_bootstrap_replicates != 200L) ||
    n_distinct(te_validation$global_fitch_observed_steps) != 1L ||
    unique(te_validation$global_fitch_observed_steps) != 2L
) {
  stop("Figure 2 requires the exact validated 34-tip Bonsai publication release")
}

bonsai_tip_count <- Ntip(bonsai_tree)
bonsai_node_x <- node.depth.edgelength(bonsai_tree)
bonsai_node_y <- node.height(bonsai_tree)
bonsai_horizontal <- tibble(
  parent = bonsai_tree$edge[, 1],
  child = bonsai_tree$edge[, 2]
) %>%
  mutate(
    x = bonsai_node_x[parent],
    xend = bonsai_node_x[child],
    y = bonsai_node_y[child],
    yend = bonsai_node_y[child]
  )
bonsai_vertical <- tibble(parent = sort(unique(bonsai_tree$edge[, 1]))) %>%
  rowwise() %>%
  mutate(
    x = bonsai_node_x[parent],
    ymin = min(bonsai_node_y[bonsai_tree$edge[bonsai_tree$edge[, 1] == parent, 2]]),
    ymax = max(bonsai_node_y[bonsai_tree$edge[bonsai_tree$edge[, 1] == parent, 2]])
  ) %>%
  ungroup()
bonsai_tips <- tibble(
  node = seq_len(bonsai_tip_count),
  species = bonsai_tree$tip.label,
  x = bonsai_node_x[seq_len(bonsai_tip_count)],
  y = bonsai_node_y[seq_len(bonsai_tip_count)]
) %>%
  left_join(te_annotations, by = "species", relationship = "one-to-one") %>%
  mutate(
    reproductive_strategy = factor(
      reproductive_strategy,
      levels = te_strategy_levels
    )
  )

descendant_tips <- function(node) {
  if (node <= bonsai_tip_count) return(bonsai_tree$tip.label[[node]])
  children <- bonsai_tree$edge[bonsai_tree$edge[, 1] == node, 2]
  unlist(lapply(children, descendant_tips), use.names = FALSE)
}
strategy_groups <- split(te_annotations$species, te_annotations$reproductive_strategy)
all_bonsai_tips <- bonsai_tree$tip.label
bonsai_nodes <- tibble(
  node = bonsai_tip_count + seq_len(Nnode(bonsai_tree)),
  x = bonsai_node_x[bonsai_tip_count + seq_len(Nnode(bonsai_tree))],
  y = bonsai_node_y[bonsai_tip_count + seq_len(Nnode(bonsai_tree))],
  support = suppressWarnings(as.numeric(bonsai_tree$node.label))
) %>%
  rowwise() %>%
  mutate(
    focal_strategy = {
      side <- descendant_tips(node)
      matches <- names(strategy_groups)[vapply(
        strategy_groups,
        function(group) {
          setequal(side, group) || setequal(setdiff(all_bonsai_tips, side), group)
        },
        logical(1)
      )]
      if (length(matches) == 1L) matches[[1]] else NA_character_
    },
    support_label = if_else(
      !is.na(support) & (support >= 50 | !is.na(focal_strategy)),
      paste0(round(support), "%"),
      NA_character_
    ),
    support_colour = if_else(
      is.na(focal_strategy),
      GBE_COLORS[["charcoal"]],
      te_strategy_colors[focal_strategy]
    )
  ) %>%
  ungroup()
bonsai_cladewise_edges <- reorder.phylo(bonsai_tree, "cladewise")$edge
bonsai_topology_depth <- numeric(bonsai_tip_count + Nnode(bonsai_tree))
bonsai_root <- setdiff(bonsai_cladewise_edges[, 1], bonsai_cladewise_edges[, 2])
if (length(bonsai_root) != 1L) stop("Bonsai topology display requires one root")
for (edge_index in seq_len(nrow(bonsai_cladewise_edges))) {
  parent <- bonsai_cladewise_edges[edge_index, 1]
  child <- bonsai_cladewise_edges[edge_index, 2]
  bonsai_topology_depth[child] <- bonsai_topology_depth[parent] + 1
}
bonsai_tip_x <- max(bonsai_topology_depth) + 0.5
bonsai_topology_horizontal <- tibble(
  parent = bonsai_tree$edge[, 1],
  child = bonsai_tree$edge[, 2]
) %>%
  mutate(
    x = bonsai_topology_depth[parent],
    xend = if_else(child <= bonsai_tip_count, bonsai_tip_x, bonsai_topology_depth[child]),
    y = bonsai_node_y[child],
    yend = bonsai_node_y[child]
  )
bonsai_topology_vertical <- bonsai_vertical %>%
  mutate(x = bonsai_topology_depth[parent])
bonsai_topology_tips <- bonsai_tips %>% mutate(x = bonsai_tip_x)
bonsai_topology_nodes <- bonsai_nodes %>%
  mutate(x = bonsai_topology_depth[node]) %>%
  filter(!is.na(focal_strategy)) %>%
  left_join(
    te_validation %>% transmute(
      focal_strategy = strategy,
      exact_support = 100 * feature_bootstrap_edge_delimited_count / feature_bootstrap_replicates
    ), by = "focal_strategy", relationship = "many-to-one"
  ) %>%
  mutate(support_label = sprintf("%.1f%%", exact_support))
if (
  length(bonsai_root) != 1L ||
    nrow(bonsai_topology_horizontal) != 66L ||
    nrow(bonsai_topology_vertical) != 33L ||
    nrow(bonsai_topology_nodes) != 2L ||
    any(!is.finite(bonsai_topology_nodes$exact_support)) ||
    any(abs(bonsai_topology_nodes$support - bonsai_topology_nodes$exact_support) > 0.5 + 1e-10)
) {
  stop("Figure 2 topology display must retain 34 tips and the two exact strategy edges")
}

te_association_tests <- read_csv(
  file.path(
    ROOT,
    "analyses",
    "01_transposable_elements",
    "data",
    "te_pca_trait_association_tests.csv"
  ),
  show_col_types = FALSE,
  progress = FALSE
)
te_strategy_permanova <- te_association_tests %>%
  filter(test_id == "reproductive_strategy_permanova")
te_strategy_permdisp <- te_association_tests %>%
  filter(test_id == "reproductive_strategy_permdisp")
if (
  nrow(te_strategy_permanova) != 1L ||
    nrow(te_strategy_permdisp) != 1L ||
    te_strategy_permanova$n_species != 34L ||
    te_strategy_permanova$p_value != 0.00001 ||
    te_strategy_permdisp$p_value <= 0.05
) {
  stop("Figure 2 requires the released post-hoc strategy association tests")
}

te_b <- ggplot() +
  geom_segment(
    data = bonsai_topology_horizontal,
    aes(x = x, xend = xend, y = y, yend = yend),
    linewidth = GBE_LINE_MM,
    lineend = "square",
    colour = GBE_COLORS[["charcoal"]]
  ) +
  geom_segment(
    data = bonsai_topology_vertical,
    aes(x = x, xend = x, y = ymin, yend = ymax),
    linewidth = GBE_LINE_MM,
    lineend = "square",
    colour = GBE_COLORS[["charcoal"]]
  ) +
  geom_point(
    data = bonsai_topology_tips,
    aes(x = x, y = y, fill = reproductive_strategy, shape = reproductive_strategy),
    size = 1.85,
    stroke = 0.45,
    colour = GBE_COLORS[["charcoal"]]
  ) +
  geom_text(
    data = bonsai_topology_tips,
    aes(x = x + 0.34, y = y, label = display_species),
    hjust = 0,
    size = GBE_LABEL_SIZE_MM,
    family = resolve_gbe_font(),
    fontface = "italic",
    colour = GBE_COLORS[["charcoal"]]
  ) +
  geom_label(
    data = bonsai_topology_nodes,
    aes(x = x - 0.12, y = y, label = support_label, colour = support_colour),
    hjust = 1,
    size = GBE_LABEL_SIZE_MM,
    family = resolve_gbe_font(),
    linewidth = 0,
    label.padding = unit(0.6, "pt"),
    fill = alpha("white", 0.9),
    fontface = "bold",
    show.legend = FALSE
  ) +
  scale_fill_manual(
    values = te_strategy_colors,
    breaks = te_strategy_levels,
    labels = te_strategy_labels,
    name = NULL
  ) +
  scale_shape_manual(
    values = te_strategy_shapes,
    breaks = te_strategy_levels,
    labels = te_strategy_labels,
    name = NULL
  ) +
  scale_colour_identity() +
  scale_x_continuous(
    limits = c(-0.4, bonsai_tip_x + 4.2),
    breaks = NULL,
    expand = c(0, 0)
  ) +
  scale_y_continuous(limits = c(0.5, bonsai_tip_count + 0.5), expand = c(0, 0)) +
  labs(
    title = "TE-composition topology",
    x = "Fitted branch lengths omitted;\nnot a species phylogeny",
    y = NULL
  ) +
  coord_cartesian(clip = "off") +
  theme_gbe() +
  theme(
    plot.title = element_text(size = 8, face = "bold", hjust = 0.5),
    axis.line = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    plot.margin = margin(3, 13, 3, 3, unit = "pt")
  )

z_value <- qnorm(0.975)
te_recovery <- te_validation %>%
  mutate(
    strategy = factor(strategy, levels = te_strategy_levels),
    strategy_short = factor(
      recode(
        as.character(strategy),
        aquatic_eggs = "Aquatic eggs",
        terrestrial_eggs_aquatic_larvae = "Terrestrial eggs +\naquatic larvae",
        direct_development = "Direct development"
      ),
      levels = rev(c(
        "Aquatic eggs",
        "Terrestrial eggs +\naquatic larvae",
        "Direct development"
      ))
    ),
    denominator = 1 + z_value^2 / feature_bootstrap_replicates,
    interval_center = (
      feature_bootstrap_edge_delimited_fraction +
        z_value^2 / (2 * feature_bootstrap_replicates)
    ) / denominator,
    interval_half = z_value * sqrt(
      feature_bootstrap_edge_delimited_fraction *
        (1 - feature_bootstrap_edge_delimited_fraction) /
        feature_bootstrap_replicates +
        z_value^2 / (4 * feature_bootstrap_replicates^2)
    ) / denominator,
    ci_low = pmax(interval_center - interval_half, 0),
    ci_high = pmin(interval_center + interval_half, 1),
    count_label = sprintf(
      "%d/%d",
      feature_bootstrap_edge_delimited_count,
      feature_bootstrap_replicates
    )
  )

te_c <- ggplot(te_recovery, aes(feature_bootstrap_edge_delimited_fraction, strategy_short)) +
  geom_segment(
    aes(x = ci_low, xend = ci_high, yend = strategy_short, colour = strategy),
    linewidth = 0.55
  ) +
  geom_point(
    aes(fill = strategy, shape = strategy),
    size = 2.7,
    stroke = 0.5,
    colour = GBE_COLORS[["charcoal"]]
  ) +
  geom_text(
    aes(x = 0.98, label = count_label, colour = strategy),
    hjust = 1,
    size = GBE_LABEL_SIZE_MM,
    family = resolve_gbe_font(),
    fontface = "bold",
    show.legend = FALSE
  ) +
  scale_colour_manual(values = te_strategy_colors, guide = "none") +
  scale_fill_manual(values = te_strategy_colors, guide = "none") +
  scale_shape_manual(values = te_strategy_shapes, guide = "none") +
  scale_x_continuous(
    limits = c(0, 1),
    breaks = c(0, 0.5, 1),
    labels = percent_format(accuracy = 1),
    expand = c(0, 0)
  ) +
  labs(title = "Feature-resampling stability",
       x = "Exact-edge recovery\n(200 feature refits; 95% interval)", y = NULL) +
  theme_gbe() +
  theme(axis.text.y = element_text(size = GBE_MIN_TEXT_PT), legend.position = "none",
        plot.title = element_text(size = 8, face = "bold", hjust = 0.5))

fitch_observed <- unique(te_validation$global_fitch_observed_steps)
fitch_null_mean <- unique(te_validation$global_fitch_null_mean_steps)
te_d <- ggplot() +
  geom_segment(
    aes(x = fitch_observed, xend = fitch_null_mean, y = 1, yend = 1),
    linewidth = 0.7,
    colour = GBE_COLORS[["light_gray"]]
  ) +
  geom_point(
    aes(x = fitch_null_mean, y = 1),
    shape = 21,
    size = 3,
    stroke = 0.55,
    fill = "white",
    colour = GBE_COLORS[["gray"]]
  ) +
  geom_point(
    aes(x = fitch_observed, y = 1),
    shape = 21,
    size = 3.2,
    stroke = 0.55,
    fill = GBE_COLORS[["charcoal"]],
    colour = GBE_COLORS[["charcoal"]]
  ) +
  annotate(
    "text",
    x = fitch_observed,
    y = 1.18,
    label = sprintf("Observed = %d", fitch_observed),
    hjust = 0,
    size = GBE_LABEL_SIZE_MM,
    fontface = "bold"
  ) +
  annotate(
    "text",
    x = fitch_null_mean,
    y = 0.82,
    label = sprintf("Null mean = %.2f", fitch_null_mean),
    hjust = 1,
    size = GBE_LABEL_SIZE_MM,
    colour = GBE_COLORS[["gray"]]
  ) +
  annotate(
    "text",
    x = 0.2,
    y = 1.56,
    hjust = 0,
    label = sprintf(
      "Reproductive strategy\nPERMANOVA R² = %.3f; P = 0.00001\nPERMDISP P = %.3f",
      te_strategy_permanova$variance_explained_r2,
      te_strategy_permdisp$p_value
    ),
    size = GBE_LABEL_SIZE_MM,
    lineheight = 0.95,
    fontface = "plain"
  ) +
  scale_x_continuous(limits = c(0, 11), breaks = c(0, 2, 4, 6, 8, 10), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0.58, 1.8), breaks = NULL, expand = c(0, 0)) +
  labs(title = "Strategy association",
       x = "Fitch steps (99,999 permutations;\nP = 0.00001)", y = NULL) +
  theme_gbe() +
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(),
        plot.title = element_text(size = 8, face = "bold", hjust = 0.5))

figure_2_layout <- "
AAB
AAB
AAC
"
Figure_2 <- te_b + te_c + te_d +
  plot_layout(
    design = figure_2_layout,
    widths = c(1, 1, 1.2),
    heights = c(1, 1, 1),
    guides = "collect"
  ) +
  plot_annotation(tag_levels = "a") &
  theme(
    legend.position = "bottom",
    plot.tag.position = c(0.005, 1)
  )

record_figure(
  "Figure_2", "Figure_2_te_diversity_and_bonsai", GBE_DOUBLE_COLUMN_MM, 190,
  "S42;S43;te_pca_trait_association_tests;trees/te34_bonsai_composition_tree_with_feature_support.nwk",
  "Panel a shows the 34-tip Bonsai topology fitted to the retry-averaged 24-superfamily CLR matrix. Branch lengths are deliberately omitted from this main display so the full topology remains legible; the display root has no phylogenetic meaning, and this is not a species phylogeny. Tip symbols show the post-hoc reproductive-strategy annotation, including D. fuscus as terrestrial eggs with aquatic larvae. Only the exact aquatic-egg and direct-development edges are labeled with feature-resampling split support. Panel b shows exact-edge recovery with 95% binomial intervals across 200 bounded feature-resampling refits: aquatic eggs 55/200, direct development 70/200, and terrestrial eggs with aquatic larvae 2/200. Panel c shows the global strategy-label localization (2 observed Fitch steps versus a 9.41-step permutation-null mean; 99,999 permutations, P = 0.00001) and the post-hoc Aitchison-distance association (PERMANOVA R² = 0.260, P = 0.00001; PERMDISP P = 0.707). Ordinary label permutations are not phylogenetically corrected, and feature resampling perturbs compositional coordinates rather than biological specimens; the result is strategy-associated compositional structure, not evidence for discrete natural clusters or causation. Exact fitted branch lengths and PCA outputs remain supplementary audits.",
  "Panel a is a topology-only TE-composition tree for 34 species with all tips colored and shaped by post-hoc reproductive strategy; aquatic-egg species and the three direct developers localize to exact edges. Panel b gives exact-edge recovery estimates and binomial intervals from 200 feature refits. Panel c contrasts two observed Fitch steps with a permutation-null mean of 9.41 and summarizes the post-hoc multivariate association and dispersion diagnostic.",
  "main_candidate_bonsai_with_pca_supplementary_audit"
)

if (identical(figure_options$figure, "Figure_2")) {
  manifest_rows[[1]]$title <- "TE composition and reproductive strategy"
  figure_manifest <- write_numbered_figure_metadata()
  source(file.path(ROOT, "scripts/publication/figure2_panels.R"))
  export_figure2_panels(list(a = te_b, b = te_c, c = te_d), OUT)
  message("Wrote Figure 2 and its three isolated panels to ", OUT)
  quit(save = "no", status = 0, runLast = FALSE)
}

## ---- figure-3
# Figure 3: fixed-scale profiles and species-level low-divergence fractions.
# The same all-hit-bp denominator is retained across species and TE orders.
landscape <- read_supp("Supplementary_Data_S12_repeatmasker_divergence_landscape.csv")
profile_orders <- c("LTR", "LINE", "TIR")
profile_group_levels <- c(
  "All species (n = 34)",
  "Aquatic eggs (n = 7)",
  "Terrestrial eggs + aquatic larvae (n = 24)",
  "Direct development (n = 3)"
)
profile_group_colors <- c(
  `All species (n = 34)` = GBE_COLORS[["charcoal"]],
  `Aquatic eggs (n = 7)` = te_strategy_colors[["aquatic_eggs"]],
  `Terrestrial eggs + aquatic larvae (n = 24)` =
    te_strategy_colors[["terrestrial_eggs_aquatic_larvae"]],
  `Direct development (n = 3)` = te_strategy_colors[["direct_development"]]
)

landscape_profiles <- landscape %>%
  filter(order %in% profile_orders, divergence_bin_start_pct < 50) %>%
  select(species, order, divergence_bin_start_pct, percent_species_hit_bp) %>%
  right_join(
    expand_grid(
      species = sort(unique(landscape$species)),
      order = profile_orders,
      divergence_bin_start_pct = 0:49
    ),
    by = c("species", "order", "divergence_bin_start_pct")
  ) %>%
  mutate(
    percent_species_hit_bp = replace_na(percent_species_hit_bp, 0),
    order = factor(order, levels = profile_orders)
  ) %>%
  left_join(
    te_annotations %>% select(species, reproductive_strategy),
    by = "species",
    relationship = "many-to-one"
  )
if (
  nrow(landscape_profiles) != 34L * 3L * 50L ||
    any(is.na(landscape_profiles$reproductive_strategy))
) {
  stop("Figure 3 requires the complete zero-filled TE34 LTR/LINE/TIR grid")
}

profile_display <- bind_rows(
  landscape_profiles %>% mutate(profile_group = "All species (n = 34)"),
  landscape_profiles %>%
    mutate(
      profile_group = recode(
        reproductive_strategy,
        aquatic_eggs = "Aquatic eggs (n = 7)",
        terrestrial_eggs_aquatic_larvae =
          "Terrestrial eggs + aquatic larvae (n = 24)",
        direct_development = "Direct development (n = 3)"
      )
    )
) %>%
  mutate(
    profile_group = factor(profile_group, levels = profile_group_levels),
    order = factor(order, levels = profile_orders)
  )
profile_means <- profile_display %>%
  group_by(profile_group, order, divergence_bin_start_pct) %>%
  summarise(
    mean_percent = mean(percent_species_hit_bp),
    n_species = n(),
    .groups = "drop"
  )
if (
  !identical(
    sort(unique(profile_means$n_species)),
    sort(c(3L, 7L, 24L, 34L))
  )
) {
  stop("Figure 3 strategy means must retain the exact 34/7/24/3 species counts")
}

# All-species means are retained as a check; only the three strategy means
# are drawn, avoiding duplication of the same 34 species in four facet rows.
profile_strategy_means <- profile_means %>%
  filter(profile_group != "All species (n = 34)")
profile_group_linetypes <- setNames(c("solid", "longdash", "dotdash"),
                                   profile_group_levels[-1])
if (any(!is.finite(landscape_profiles$percent_species_hit_bp)) ||
    any(landscape_profiles$percent_species_hit_bp < 0 |
        landscape_profiles$percent_species_hit_bp > 3)) {
  stop("Figure 3 profiles exceed the fixed 0-3% all-hit-bp axis")
}

repeat_a <- ggplot(
  landscape_profiles,
  aes(divergence_bin_start_pct, percent_species_hit_bp, group = species)
) +
  geom_line(colour = GBE_COLORS[["gray"]], linewidth = 0.2, alpha = 0.2) +
  geom_line(
    data = profile_strategy_means,
    aes(divergence_bin_start_pct, mean_percent,
        group = profile_group, colour = profile_group, linetype = profile_group),
    inherit.aes = FALSE, linewidth = 0.65
  ) +
  facet_grid(cols = vars(order)) +
  scale_colour_manual(
    values = profile_group_colors, breaks = profile_group_levels[-1],
    labels = c("Aquatic eggs (n = 7)", "Terrestrial eggs + aquatic larvae (n = 24)",
               "Direct development (n = 3)"), name = NULL
  ) +
  scale_linetype_manual(
    values = profile_group_linetypes, breaks = profile_group_levels[-1],
    labels = c("Aquatic eggs (n = 7)", "Terrestrial eggs + aquatic larvae (n = 24)",
               "Direct development (n = 3)"), name = NULL
  ) +
  scale_x_continuous(limits = c(0, 50), breaks = c(0, 25, 50), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 3), breaks = 0:3, expand = expansion(mult = c(0, .02))) +
  labs(
    x = "RepeatMasker sequence divergence (%)",
    y = "All aligned hit bp per bin (%)",
    caption = "Thin gray curves: individual species. Coloured curves: equal-weight strategy means."
  ) +
  theme_gbe() + theme_gbe_framed_facets() +
  theme(text = element_text(family = GBE_FONT),
        panel.spacing = unit(12, "pt"), legend.key.width = unit(8, "mm"),
        legend.text = element_text(size = GBE_MIN_TEXT_PT))

young_repeat_matrix <- landscape_profiles %>%
  filter(divergence_bin_start_pct < 5) %>%
  group_by(species, order, reproductive_strategy) %>%
  summarise(young_repeat_percent = sum(percent_species_hit_bp), .groups = "drop") %>%
  mutate(
    species_label = factor(
      species_label(species),
      levels = phylogenetic_species_levels(species)
    ),
    order_index = match(order, profile_orders),
    reproductive_strategy = factor(
      reproductive_strategy,
      levels = te_strategy_levels
    )
  )
if (nrow(young_repeat_matrix) != 34L * 3L ||
    any(!is.finite(young_repeat_matrix$young_repeat_percent)) ||
    any(young_repeat_matrix$young_repeat_percent < 0 |
        young_repeat_matrix$young_repeat_percent > 6)) {
  stop("Figure 3 requires 34 species by three orders within the fixed 0-6% axis")
}

repeat_b <- ggplot(young_repeat_matrix, aes(young_repeat_percent, species_label)) +
  geom_segment(aes(x = 0, xend = young_repeat_percent, yend = species_label),
               colour = GBE_COLORS[["light_gray"]], linewidth = 0.35) +
  geom_point(aes(colour = reproductive_strategy, shape = reproductive_strategy,
                 fill = reproductive_strategy), size = 1.8, stroke = 0.35) +
  facet_grid(cols = vars(order)) +
  scale_colour_manual(values = te_strategy_colors, breaks = te_strategy_levels,
                       labels = te_strategy_labels, name = NULL) +
  scale_fill_manual(values = te_strategy_colors, breaks = te_strategy_levels,
                     labels = te_strategy_labels, name = NULL) +
  scale_shape_manual(values = te_strategy_shapes, breaks = te_strategy_levels,
                      labels = te_strategy_labels, name = NULL) +
  scale_x_continuous(limits = c(0, 6), breaks = c(0, 2, 4, 6),
                     expand = expansion(mult = c(.025, .035))) +
  labs(x = "All aligned hit bp below 5% divergence (%)", y = NULL) +
  theme_gbe() + theme_gbe_framed_facets() +
  theme(
    text = element_text(family = GBE_FONT),
    axis.text.y = element_text(face = "italic", size = GBE_MIN_TEXT_PT),
    panel.spacing = unit(6, "pt"),
    legend.position = "bottom"
  )

Figure_3 <- repeat_a / repeat_b +
  plot_layout(heights = c(1, 2.4)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag.position = c(0.005, 1))

record_figure(
  "Figure_3", "Figure_3_repeat_divergence_landscape", GBE_DOUBLE_COLUMN_MM, 220,
  "S12;S02;trees/desmo900dated_source46.tre",
  "RepeatMasker divergence profiles vary among species and reproductive strategies. Panel a shows LTR, LINE, and TIR profiles for all 34 species as thin gray curves and equal-weight means for aquatic eggs (7 species), terrestrial eggs with aquatic larvae (24), and direct development (3) as coloured curves with distinct line types. All three panels share 0-50% divergence and 0-3% of all species RepeatMasker-aligned hit bp per 1-percentage-point bin; absent order/bin combinations are zero-filled for bins 0-49 without renormalizing the displayed subset. The strategy strata are descriptive, not phylogenetically corrected comparisons. Panel b uses three matched 0-6% dot-plot axes to show each species' fraction of all aligned hit bp below 5% divergence (bins 0-4), in the shared supplied-tree order. Colour and shape identify reproductive strategy. These are observed assay fractions without biological confidence intervals. Divergence is not insertion age, and aligned-hit-bp fractions are not genome occupancy.",
  "Panel a has three shared-scale facets for LTR, LINE, and TIR, showing 34 individual-species curves and three strategy means. Panel b has 34 labelled species in shared phylogenetic order and three matched dot-plot columns for fractions below five percent divergence, with strategy colours and shapes.",
  "main_candidate"
)

if (identical(figure_options$figure, "Figure_3")) {
  manifest_rows[[1]]$title <- "Repeat divergence across 34 species"
  figure_manifest <- write_numbered_figure_metadata()
  source(file.path(ROOT, "scripts/publication/figure3_panels.R"))
  export_figure3_panels(list(a = repeat_a, b = repeat_b), OUT)
  message("Wrote Figure 3 and its two isolated panels to ", OUT)
  quit(save = "no", status = 0, runLast = FALSE)
}

## ---- figure-4
# Figure 4: terminal:internal LTR deletion-footprint proxy.
ltr <- read_supp("Supplementary_Data_S15_ltr_terminal_internal_species_robustness.csv") %>%
  filter(analysis_branch == "primary_iqr_filtered_coverage_ge_80pct", ratio_median > 0) %>%
  mutate(
    species_label = species_label(species),
    estimate = log2(ratio_median),
    ci_low = log2(median_bootstrap_ci_low),
    ci_high = log2(median_bootstrap_ci_high)
  )

ltr_tree <- prune_source_tree(ltr$species)

ltr_tip_count <- Ntip(ltr_tree)
ltr_node_count <- Nnode(ltr_tree)
ltr_node_y <- node.height(ltr_tree)
ltr_edges <- reorder.phylo(ltr_tree, "cladewise")$edge
ltr_topology_depth <- numeric(ltr_tip_count + ltr_node_count)
ltr_root <- setdiff(ltr_edges[, 1], ltr_edges[, 2])
if (length(ltr_root) != 1L) stop("Figure 4 LTR tree requires one root")
for (edge_index in seq_len(nrow(ltr_edges))) {
  parent <- ltr_edges[edge_index, 1]
  child <- ltr_edges[edge_index, 2]
  ltr_topology_depth[child] <- ltr_topology_depth[parent] + 1
}
ltr_tip_x <- max(ltr_topology_depth) + 0.5

ltr_tree_horizontal <- tibble(
  parent = ltr_edges[, 1],
  child = ltr_edges[, 2]
) %>%
  mutate(
    x = ltr_topology_depth[parent],
    xend = if_else(child <= ltr_tip_count, ltr_tip_x, ltr_topology_depth[child]),
    y = ltr_node_y[child],
    yend = ltr_node_y[child]
  )
ltr_tree_vertical <- tibble(parent = sort(unique(ltr_edges[, 1]))) %>%
  rowwise() %>%
  mutate(
    x = ltr_topology_depth[parent],
    ymin = min(ltr_node_y[ltr_edges[ltr_edges[, 1] == parent, 2]]),
    ymax = max(ltr_node_y[ltr_edges[ltr_edges[, 1] == parent, 2]])
  ) %>%
  ungroup()
ltr_tree_tips <- tibble(
  node = seq_len(ltr_tip_count),
  species = ltr_tree$tip.label,
  species_label = species_label(ltr_tree$tip.label),
  x = ltr_tip_x,
  y = ltr_node_y[seq_len(ltr_tip_count)]
)
ltr_tree_tips_sorted <- ltr_tree_tips %>% arrange(y)
ltr_tip_y_order <- ltr_tree_tips_sorted$species
ltr_canonical_order <- phylogenetic_species_order(ltr$species)

if (
  ltr_tip_count != 30L ||
    ltr_node_count != 29L ||
    !is.rooted(ltr_tree) ||
    !is.binary(ltr_tree) ||
    !is.ultrametric(ltr_tree, tol = 1e-6) ||
    any(ltr_tree$edge.length < 0) ||
    !setequal(ltr_tree$tip.label, ltr$species) ||
    nrow(ltr_tree_horizontal) != 58L ||
    nrow(ltr_tree_vertical) != 29L ||
    !identical(ltr_tip_y_order, ltr_canonical_order) ||
    !identical(ltr_tree_tips_sorted$y, as.numeric(seq_len(ltr_tip_count)))
) {
  stop("Figure 4 requires the exact rooted 30-tip prune in canonical tip-y order")
}

ltr <- ltr %>%
  left_join(
    ltr_tree_tips %>% select(species, y),
    by = "species",
    relationship = "one-to-one"
  )

ltr_elements <- read_supp("Supplementary_Data_S14_ltr_terminal_internal_elements.csv") %>%
  filter(primary_ltr_element) %>%
  mutate(
    estimate = log2_ratio_terminal_internal_all_positions
  ) %>%
  left_join(
    ltr_tree_tips %>% select(species, y),
    by = "species",
    relationship = "many-to-one"
  )
if (nrow(ltr_elements) != 380L || n_distinct(ltr_elements$species) != 30L ||
    any(!is.finite(ltr_elements$estimate)) || any(is.na(ltr_elements$y)) ||
    any(is.na(ltr$y))) {
  stop("Figure 4 must expose the exact 380-element, 30-species primary LTR branch")
}

ltr_element_checks <- ltr_elements %>%
  group_by(species) %>%
  summarise(
    element_n = n(),
    element_median = median(ratio_terminal_internal_all_positions),
    .groups = "drop"
  ) %>%
  right_join(
    ltr %>% select(species, n_elements, ratio_median),
    by = "species"
  ) %>%
  left_join(
    read_supp("Supplementary_Data_S16_ltr_resource_coverage.csv") %>%
      select(species, n_primary_iqr_filtered_elements),
    by = "species"
  )
if (any(is.na(ltr_element_checks$element_n)) ||
    !all(ltr_element_checks$element_n == ltr_element_checks$n_elements) ||
    !all(ltr_element_checks$n_elements == ltr_element_checks$n_primary_iqr_filtered_elements) ||
    !isTRUE(all.equal(
      ltr_element_checks$element_median,
      ltr_element_checks$ratio_median,
      tolerance = 1e-12,
      check.attributes = FALSE
    ))) {
  stop("Figure 4 element observations, medians, and displayed support must agree")
}

# Draw the tree and raw element-level proxy on a shared continuous tip-y scale.
# This makes the supplied-tree geometry, species labels, observations, and
# species summaries use one canonical row order rather than relying on the
# phylo object's storage order.
ltr_tree_left <- ggplot() +
  geom_segment(
    data = ltr_tree_horizontal,
    aes(x = x, xend = xend, y = y, yend = yend),
    linewidth = 0.35,
    lineend = "square",
    colour = GBE_COLORS[["charcoal"]]
  ) +
  geom_segment(
    data = ltr_tree_vertical,
    aes(x = x, xend = x, y = ymin, yend = ymax),
    linewidth = 0.35,
    lineend = "square",
    colour = GBE_COLORS[["charcoal"]]
  ) +
  scale_x_continuous(
    limits = c(-0.2, ltr_tip_x),
    breaks = NULL,
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    limits = c(0.5, ltr_tip_count + 0.5),
    breaks = ltr_tree_tips_sorted$y,
    labels = ltr_tree_tips_sorted$species_label,
    expand = c(0, 0),
    position = "right"
  ) +
  labs(
    x = "Species phylogeny (topology only)",
    y = NULL,
    tag = "a"
  ) +
  theme_gbe() +
  theme(
    text = element_text(family = GBE_FONT),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_text(
      size = GBE_MIN_TEXT_PT,
      face = "italic",
      margin = margin(l = 2, unit = "pt")
    ),
    plot.margin = margin(4, 2, 4, 4, unit = "pt")
  )

# Estimate density in the displayed log2 space with one shared bandwidth.
# Equal maximum widths show distribution shape; the count column carries n.
ltr_violin_bandwidth <- stats::bw.nrd0(ltr_elements$estimate)
ltr_genome_labels <- ltr_tree_tips_sorted %>%
  select(species, y) %>%
  left_join(
    read_supp("Supplementary_Data_S04_path_analysis_traits24.csv") %>%
      transmute(species = tree_tip,
                genome_size_gbp = genome_size_pg_fuscus_anchored * 0.978),
    by = "species", relationship = "one-to-one"
  ) %>%
  left_join(ltr %>% select(species, n_elements),
            by = "species", relationship = "one-to-one") %>%
  mutate(label = paste0(if_else(is.na(genome_size_gbp), "n/a",
                                sprintf("%.1f", genome_size_gbp)),
                        " / ", n_elements))
ltr_distribution <- ggplot(
  ltr_elements,
  aes(x = estimate, y = y)
) +
  geom_violin(
    aes(group = species),
    orientation = "y",
    position = "identity",
    bw = ltr_violin_bandwidth,
    trim = TRUE,
    scale = "width",
    width = 0.8,
    quantiles = NULL,
    fill = "white",
    colour = "#B3B3B3",
    linewidth = 0.35
  ) +
  geom_vline(
    xintercept = 0,
    linetype = 3,
    linewidth = 0.5,
    colour = GBE_COLORS[["charcoal"]]
  ) +
  geom_point(
    position = position_jitter(width = 0, height = 0.15, seed = 20260828),
    shape = 21,
    size = 0.95,
    stroke = 0.25,
    fill = "white",
    colour = GBE_COLORS[["charcoal"]],
    alpha = 0.75
  ) +
  geom_segment(
    data = ltr,
    aes(
      x = ci_low,
      xend = ci_high,
      y = y,
      yend = y
    ),
    inherit.aes = FALSE,
    linewidth = 0.35,
    lineend = "butt",
    colour = GBE_COLORS[["blue"]]
  ) +
  geom_segment(
    data = ltr,
    aes(x = estimate, xend = estimate, y = y - 0.25, yend = y + 0.25),
    inherit.aes = FALSE,
    linewidth = 0.9,
    lineend = "butt",
    colour = GBE_COLORS[["blue"]]
  ) +
  scale_x_continuous(
    breaks = log2(c(0.025, 0.05, 0.1, 0.25, 0.5, 1, 1.5)),
    labels = c("0.025", "0.05", "0.1", "0.25", "0.5", "1", "1.5"),
    limits = c(log2(0.025), log2(1.6)),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    limits = c(0.5, ltr_tip_count + 0.5),
    breaks = ltr_tree_tips_sorted$y,
    expand = c(0, 0),
    sec.axis = dup_axis(
      breaks = ltr_tree_tips_sorted$y,
      labels = ltr_genome_labels$label,
      name = "Estimated genome size (Gbp*) / elements (n)"
    )
  ) +
  labs(
    x = "LTR terminal-to-internal depth ratio\n(log2 scale; 1 = equal depth)",
    y = NULL,
    tag = "b",
    caption = paste0("Outline: density; points: elements; bar: median.\n",
                     "Thin line: 95% interval. Gbp*: assembly-anchored estimate.\n",
                     "n/a: genome size unavailable.")
  ) +
  theme_gbe() +
  theme(
    text = element_text(family = GBE_FONT),
    panel.border = element_rect(fill = NA, colour = "#B3B3B3", linewidth = 0.35),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    axis.text.y.right = element_text(family = GBE_FONT, size = GBE_MIN_TEXT_PT),
    axis.title.y.right = element_text(family = GBE_FONT, size = GBE_MIN_TEXT_PT),
    plot.margin = margin(4, 4, 4, 2, unit = "pt")
  )

Figure_4 <- ltr_tree_left + ltr_distribution +
  plot_layout(widths = c(1.18, 1.82))

record_figure(
  "Figure_4", "Figure_4_ltr_terminal_internal_proxy", GBE_DOUBLE_COLUMN_MM, 190,
  "S04;S14;S15;S16;S17;trees/desmo900dated_source46.tre",
  "Panel a shows a topology-only 30-tip prune of the supplied 46-tip species tree; branch lengths are omitted because its publication or archive identifier, calibration method, and branch-length provenance remain unresolved. Species rows are derived from the plotted tip-y order and therefore follow the shared canonical phylogenetic order. Panel b shows the zero-aware terminal:internal LTR depth proxy across the exact 380-element primary branch. Light gray violin outlines show smoothed distributions of log2 element ratios, using one bandwidth selected from the pooled log2 observations, equal maximum widths, and tails trimmed to each species' observed range. Hollow points are element-level observations; short vertical blue bars are released species medians of untransformed ratios, displayed after log2 transformation. Thin horizontal blue intervals are 95% percentile intervals from 2,000 within-species element-bootstrap replicates. The right-hand labels give estimated genome size in Gbp followed by retained element count (3-50), separated by a slash. Genome sizes use released D. fuscus assembly-anchored IOD-derived estimates converted at 0.978 Gbp per pg; the asterisk denotes assembly-anchored estimates, not independently measured species C-values, and n/a marks unavailable values without imputation. Violin width does not encode sample size, and shapes for species with few recovered elements are descriptive smoothing only. Elements are not biological replicates, and intervals are conditional on the recovered elements and sequencing resources. Elements were retained inclusively within each species' two-sided 1.5-IQR fences and only when the left LTR, right LTR, and internal region each had at least 80% positive-depth coverage. The primary branch derives from 408 LTR/Gypsy candidates at least 3,000 bp long with at least five correctly parsed TEsorter domain annotations. The dotted line marks equal terminal and internal mean depth (1:1); values to its right point in the direction expected under internal-sequence depletion. The tree supplies context and ordering only; this mapping/deletion-footprint proxy is not a phylogenetic analysis or an estimate of an ectopic-recombination, solo-LTR, or DNA-loss rate.",
  "Panel a shows a topology-only species tree pruned to the 30 LTR species. Panel b aligns the same 30 phylogenetic rows with light gray horizontal violins, 380 hollow element observations, short vertical blue median bars, and thin horizontal bootstrap intervals on a log2 terminal:internal ratio axis labelled with decimal ratios. A dotted vertical line marks equal terminal and internal depth at one, and the right-hand column pairs available estimated genome sizes in Gbp with retained element counts, using n/a for unavailable genome sizes.",
  "sensitivity_only_not_ectopic_rate"
)

if (identical(figure_options$figure, "Figure_4")) {
  manifest_rows[[1]]$title <- "LTR terminal-to-internal depth across 30 species"
  figure_manifest <- write_numbered_figure_metadata()
  source(file.path(ROOT, "scripts/publication/figure4_panels.R"))
  export_figure4_panels(list(a = ltr_tree_left, b = ltr_distribution),
                        ltr_tree_tips_sorted, OUT)
  message("Wrote Figure 4 and its two isolated panels to ", OUT)
  quit(save = "no", status = 0, runLast = FALSE)
}

## ---- figure-5
# Figure 5: phylogenetically aligned estimated genome size and reviewed
# upper-tail cell/nucleus morphology.  Raw objects stay in S18 and are validated
# here, but the main panel emphasizes species estimates and biological context.
morphology <- read_supp("Supplementary_Data_S19_cell_nucleus_species_estimates.csv") %>%
  mutate(
    species_label = species_label(species),
    species_label = factor(
      species_label,
      levels = phylogenetic_species_levels(species)
    ),
    morphology_support = if_else(n_size_specimens == 1L, "One animal", "Two or more animals")
  )

morphology_objects <- read_supp("Supplementary_Data_S18_frozen_cell_nucleus_objects.csv") %>%
  transmute(
    species,
    species_label = factor(
      species_label(species),
      levels = phylogenetic_species_levels(species)
    ),
    cell_area_um2,
    nucleus_area_um2 = nuc_area_um2
  )
if (
  nrow(morphology_objects) != 1152L ||
    n_distinct(morphology_objects$species) != 24L ||
    any(!is.finite(morphology_objects$cell_area_um2)) ||
    any(!is.finite(morphology_objects$nucleus_area_um2)) ||
    !setequal(morphology_objects$species, as.character(morphology$species))
) {
  stop("Figure 5 requires the exact 1,152 reviewed cell-nucleus object pairs")
}

morph_organismal <- read_csv(
  file.path(
    ROOT,
    "analyses",
    "03_phylogenetic_path",
    "data",
    "organismal_traits.csv"
  ),
  show_col_types = FALSE,
  progress = FALSE
) %>%
  select(tree_tip = species, development_mode)
morph_quality <- read_supp(
  "Supplementary_Data_S24_nuclear_iod_quality_balance.csv"
) %>%
  select(species, passes_prespecified_balance_thresholds)
morph_context <- tree_tips %>%
  select(species, tree_tip, y) %>%
  left_join(path_traits, by = c("species", "tree_tip"), relationship = "one-to-one") %>%
  left_join(
    morphology %>%
      select(
        species,
        n_size_specimens,
        cell_area_um2,
        cell_area_um2_ci_low,
        cell_area_um2_ci_high,
        nucleus_area_um2,
        nucleus_area_um2_ci_low,
        nucleus_area_um2_ci_high
      ),
    by = "species",
    relationship = "one-to-one",
    suffix = c("_path", "")
  ) %>%
  left_join(morph_quality, by = "species", relationship = "one-to-one") %>%
  left_join(morph_organismal, by = "tree_tip", relationship = "one-to-one") %>%
  mutate(
    display_species = species_label(species),
    development_mode = factor(
      development_mode,
      levels = c("aquatic_larva", "direct_development")
    ),
    genome_support = factor(
      case_when(
        !passes_prespecified_balance_thresholds ~ "Quality-balance caution",
        n_genome_size_images == 1L & n_genome_size_specimens == 1L ~
          "One image and animal",
        TRUE ~ "Balanced multi-image support"
      ),
      levels = c(
        "Balanced multi-image support",
        "One image and animal",
        "Quality-balance caution"
      )
    ),
    morphology_support = factor(
      if_else(n_size_specimens == 1L, "One animal", "Two or more animals"),
      levels = c("Two or more animals", "One animal")
    )
  )
if (
  nrow(morph_context) != 24L ||
    any(!complete.cases(morph_context[, c(
      "genome_size_pg_fuscus_anchored",
      "cell_area_um2",
      "nucleus_area_um2",
      "development_mode"
    )])) ||
    !identical(
      path_tree$tip.label[order(morph_context$y)],
      phylogenetic_species_order(path_tree$tip.label)
    )
) {
  stop("Figure 5 tracks must align to all 24 Path24 tree tips")
}

path_edges <- reorder.phylo(path_tree, "cladewise")$edge
path_topology_depth <- numeric(tip_count + node_count)
path_root <- setdiff(path_edges[, 1], path_edges[, 2])
if (length(path_root) != 1L) stop("Figure 5 Path24 tree requires one root")
for (edge_index in seq_len(nrow(path_edges))) {
  parent <- path_edges[edge_index, 1]
  child <- path_edges[edge_index, 2]
  path_topology_depth[child] <- path_topology_depth[parent] + 1
}
path_tip_x <- max(path_topology_depth) + 0.5
path_topology_horizontal <- tibble(
  parent = path_tree$edge[, 1],
  child = path_tree$edge[, 2]
) %>%
  mutate(
    x = path_topology_depth[parent],
    xend = if_else(child <= tip_count, path_tip_x, path_topology_depth[child]),
    y = node_y[child],
    yend = node_y[child]
  )
path_topology_vertical <- tree_vertical %>%
  mutate(x = path_topology_depth[parent])

morph_tree_panel <- ggplot() +
  geom_segment(
    data = path_topology_horizontal,
    aes(x = x, xend = xend, y = y, yend = yend),
    linewidth = GBE_LINE_MM,
    lineend = "square",
    colour = GBE_COLORS[["charcoal"]]
  ) +
  geom_segment(
    data = path_topology_vertical,
    aes(x = x, xend = x, y = ymin, yend = ymax),
    linewidth = GBE_LINE_MM,
    lineend = "square",
    colour = GBE_COLORS[["charcoal"]]
  ) +
  geom_point(
    data = morph_context,
    aes(
      x = path_tip_x,
      y = y,
      fill = development_mode,
      shape = development_mode
    ),
    size = 1.7,
    stroke = 0.45,
    colour = GBE_COLORS[["charcoal"]]
  ) +
  geom_text(
    data = morph_context,
    aes(x = path_tip_x + 0.28, y = y, label = display_species),
    hjust = 0,
    size = GBE_LABEL_SIZE_MM,
    family = GBE_FONT,
    fontface = "italic"
  ) +
  scale_fill_manual(
    values = c(
      aquatic_larva = GBE_COLORS[["light_gray"]],
      direct_development = GBE_COLORS[["purple"]]
    ),
    labels = c("Aquatic larvae", "Direct development"),
    name = NULL
  ) +
  scale_shape_manual(
    values = c(aquatic_larva = 21, direct_development = 24),
    labels = c("Aquatic larvae", "Direct development"),
    name = NULL,
    guide = guide_legend(nrow = 1)
  ) +
  scale_x_continuous(limits = c(-0.2, path_tip_x + 4.6), breaks = NULL, expand = c(0, 0)) +
  scale_y_continuous(limits = c(0.5, 24.5), breaks = NULL, expand = c(0, 0)) +
  labs(title = "Phylogeny and development", x = NULL, y = NULL) +
  coord_cartesian(clip = "off") +
  theme_gbe() +
  theme(
    axis.line = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    plot.title = element_text(size = GBE_BASE_SIZE, face = "bold", hjust = 0.5),
    plot.margin = margin(3, 10, 3, 3, unit = "pt")
  )

trait_forest_panel <- function(
  data,
  estimate,
  low,
  high,
  x_label,
  support,
  support_values,
  x_limits,
  x_breaks,
  show_support_legend = TRUE
) {
  ggplot(data, aes(x = .data[[estimate]], y = y)) +
    geom_segment(
      aes(x = .data[[low]], xend = .data[[high]], yend = y),
      linewidth = 0.48,
      colour = GBE_COLORS[["light_gray"]]
    ) +
    geom_point(
      aes(shape = .data[[support]]),
      size = 1.9,
      stroke = 0.55,
      fill = "white",
      colour = GBE_COLORS[["charcoal"]]
    ) +
    scale_shape_manual(
      values = support_values,
      name = NULL,
      guide = if (show_support_legend) guide_legend(nrow = 1) else "none"
    ) +
    scale_x_continuous(limits = x_limits, breaks = x_breaks,
                       expand = expansion(mult = 0.04)) +
    scale_y_continuous(limits = c(0.5, 24.5), breaks = NULL, expand = c(0, 0)) +
    labs(x = x_label, y = NULL) +
    theme_gbe() +
    theme(
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank(),
      legend.position = if (show_support_legend) "bottom" else "none"
    )
}

morph_genome_panel <- trait_forest_panel(
  morph_context,
  "genome_size_pg_fuscus_anchored",
  "genome_size_pg_ci_low",
  "genome_size_pg_ci_high",
  "Estimated genome\nsize (pg/1C*)",
  "genome_support",
  c(
    `Balanced multi-image support` = 21,
    `One image and animal` = 24,
    `Quality-balance caution` = 4
  ),
  c(5, 30), c(5, 15, 25),
  FALSE
)
morph_nucleus_panel <- trait_forest_panel(
  morph_context,
  "nucleus_area_um2",
  "nucleus_area_um2_ci_low",
  "nucleus_area_um2_ci_high",
  "Nucleus area\n(µm²)",
  "morphology_support",
  c(`Two or more animals` = 21, `One animal` = 24),
  c(30, 60), c(30, 45, 60),
  FALSE
)
morph_cell_panel <- trait_forest_panel(
  morph_context,
  "cell_area_um2",
  "cell_area_um2_ci_low",
  "cell_area_um2_ci_high",
  "Cell area\n(µm²)",
  "morphology_support",
  c(`Two or more animals` = 21, `One animal` = 24),
  c(180, 340), c(200, 260, 320),
  FALSE
)

Figure_5 <-
  morph_tree_panel +
  morph_genome_panel +
  morph_nucleus_panel +
  morph_cell_panel +
  plot_layout(widths = c(2.2, 0.78, 0.72, 0.72), guides = "keep") +
  plot_annotation(tag_levels = "a") &
  theme(
    legend.position = "bottom",
    plot.tag.position = c(0.005, 1)
  )

record_figure(
  "Figure_5", "Figure_5_phylogeny_and_cellular_traits", GBE_DOUBLE_COLUMN_MM, 190,
  "S04;S18;S19;S20;S24;organismal_traits;trees/desmognathus_time_tree_path24_v1.nwk",
  "Twenty-four species are aligned across a topology-only prune of the supplied species tree, development mode, estimated genome size, and reviewed upper-tail erythrocyte nucleus and cell areas. Branch lengths are omitted because the species tree's publication, calibration, and time provenance remain unresolved. Genome-size points are D. fuscus assembly-anchored IOD-derived estimates; triangles identify one-image, one-animal support and crosses flag the prespecified image-quality balance caution. These are not independently measured species C-values, and intervals exclude anchor uncertainty. Nucleus and cell points are species medians from the 1,152 accepted reviewed pairs retained in S18; triangles identify the seven species represented by one animal. All intervals are conditional on observed images/animals and do not represent population-level uncertainty.",
  "A topology-only tree for 24 Desmognathus species aligns development mode with estimated genome-size, nucleus-area, and cell-area forest plots. All tracks share one phylogenetic row order, and point shapes expose single-animal and image-quality support cautions.",
  "main_candidate_integrated_phylogeny_and_cellular_traits"
)

if (identical(figure_options$figure, "Figure_5")) {
  manifest_rows[[1]]$title <- "Phylogeny and cellular traits across 24 species"
  figure_manifest <- write_numbered_figure_metadata()
  source(file.path(ROOT, "scripts/publication/figure5_panels.R"))
  export_figure5_panels(
    list(a = morph_tree_panel, b = morph_genome_panel,
         c = morph_nucleus_panel, d = morph_cell_panel),
    morph_context %>% arrange(y), OUT
  )
  message("Rendered Figure_5 and four isolated panels to ", OUT)
  quit(save = "no", status = 0)
}

## ---- figure-6
# Figure 6: estimated genome size and pairwise Pagel-lambda PGLS sensitivity.
genome_size <- path_traits %>%
  transmute(
    species,
    genome_size_pg_fuscus_anchored,
    genome_size_pg_ci_low,
    genome_size_pg_ci_high,
    n_images = n_genome_size_images,
    n_specimens = n_genome_size_specimens
  ) %>%
  left_join(
    read_supp("Supplementary_Data_S24_nuclear_iod_quality_balance.csv") %>%
      select(species, passes_prespecified_balance_thresholds),
    by = "species"
  ) %>%
  mutate(
    genome_size_support = case_when(
      !passes_prespecified_balance_thresholds ~ "Quality-balance caution",
      n_images == 1L & n_specimens == 1L ~ "One image and animal",
      TRUE ~ "Balanced multiple-image support"
    )
  )
pair_data <- inner_join(
  genome_size,
  morphology %>% mutate(species = as.character(species)) %>% select(-species_label),
  by = "species"
)
pgls <- read_supp("Supplementary_Data_S26_pairwise_pgls_sensitivity.csv")

pgls_panel <- function(data, comparison_id, x, y, x_low, x_high, y_low, y_high, x_label, y_label, colour, shape_variable = NULL) {
  stat <- pgls %>% filter(.data$comparison == .env$comparison_id) %>% slice(1)
  x_values <- data[[x]]
  y_values <- data[[y]]
  x_log <- log10(x_values)
  y_log <- log10(y_values)
  x_sequence <- seq(min(x_values), max(x_values), length.out = 200)
  x_z <- (log10(x_sequence) - mean(x_log)) / sd(x_log)
  predicted_z <- stat$pgls_intercept_z + stat$pgls_standardized_beta * x_z
  model_variance <- stat$pgls_var_intercept + 2 * x_z * stat$pgls_cov_intercept_slope + x_z^2 * stat$pgls_var_slope
  critical <- qt(0.975, stat$pgls_degrees_freedom)
  model <- tibble(
    x = x_sequence,
    fit = 10^(mean(y_log) + sd(y_log) * predicted_z),
    low = 10^(mean(y_log) + sd(y_log) * (predicted_z - critical * sqrt(pmax(model_variance, 0)))),
    high = 10^(mean(y_log) + sd(y_log) * (predicted_z + critical * sqrt(pmax(model_variance, 0))))
  )
  p_text <- ifelse(
    stat$pgls_p_value < 0.001,
    "<0.001",
    paste("=", sprintf("%.3f", stat$pgls_p_value))
  )
  annotation <- sprintf(
    "β = %.2f [%.2f, %.2f]\nP %s",
    stat$pgls_standardized_beta, stat$pgls_ci_low, stat$pgls_ci_high,
    p_text
  )
  panel <- ggplot(data, aes(x = .data[[x]], y = .data[[y]])) +
    geom_ribbon(data = model, aes(x = x, ymin = low, ymax = high), inherit.aes = FALSE, fill = colour, alpha = 0.12) +
    geom_line(data = model, aes(x = x, y = fit), inherit.aes = FALSE, colour = colour, linewidth = 0.65) +
    geom_errorbar(aes(ymin = .data[[y_low]], ymax = .data[[y_high]]), width = 0, linewidth = 0.3, colour = alpha(colour, 0.35)) +
    geom_errorbar(aes(xmin = .data[[x_low]], xmax = .data[[x_high]]), orientation = "y", width = 0, linewidth = 0.3, colour = alpha(colour, 0.35))
  if (is.null(shape_variable)) {
    panel <- panel + geom_point(
      shape = 21,
      size = 1.8,
      stroke = 0.55,
      fill = "white",
      colour = colour
    )
  } else {
    panel <- panel +
      geom_point(
        aes(shape = .data[[shape_variable]]),
        size = 2,
        stroke = 0.55,
        fill = "white",
        colour = colour
      ) +
      scale_shape_manual(
        name = NULL,
        values = c(
          "Balanced multiple-image support" = 21,
          "One image and animal" = 24,
          "Quality-balance caution" = 4
        ),
        labels = c("Balanced support", "One image/animal", "Quality caution")
      )
  }
  panel +
    annotate(
      "text",
      x = -Inf,
      y = Inf,
      label = annotation,
      hjust = -0.02,
      vjust = 1.25,
      size = GBE_LABEL_SIZE_MM,
      family = GBE_FONT
    ) +
    labs(x = x_label, y = y_label) +
    theme_gbe() +
    theme(plot.margin = margin(5, 6, 5, 6, unit = "pt"))
}

genome_a <- pgls_panel(
  pair_data, "genome_size_vs_nucleus_area", "genome_size_pg_fuscus_anchored", "nucleus_area_um2",
  "genome_size_pg_ci_low", "genome_size_pg_ci_high", "nucleus_area_um2_ci_low", "nucleus_area_um2_ci_high",
  "Estimated genome size\n(pg/1C*)", "Nucleus area (µm²)", GBE_COLORS[["blue"]], "genome_size_support"
)
genome_b <- pgls_panel(
  pair_data, "genome_size_vs_cell_area", "genome_size_pg_fuscus_anchored", "cell_area_um2",
  "genome_size_pg_ci_low", "genome_size_pg_ci_high", "cell_area_um2_ci_low", "cell_area_um2_ci_high",
  "Estimated genome size\n(pg/1C*)", "Cell area (µm²)", GBE_COLORS[["orange"]], "genome_size_support"
) + scale_y_continuous(limits = c(180, 370), breaks = seq(200, 350, 50)) +
  guides(shape = "none")
genome_c <- pgls_panel(
  morphology, "nucleus_area_vs_cell_area", "nucleus_area_um2", "cell_area_um2",
  "nucleus_area_um2_ci_low", "nucleus_area_um2_ci_high", "cell_area_um2_ci_low", "cell_area_um2_ci_high",
  "Nucleus area (µm²)", "Cell area (µm²)", GBE_COLORS[["purple"]]
) + scale_y_continuous(limits = c(180, 370), breaks = seq(200, 350, 50))

Figure_6 <- genome_a + genome_b + genome_c +
  plot_layout(widths = c(1, 1, 1), guides = "collect") +
  plot_annotation(tag_levels = "a") & theme(legend.position = "bottom")

record_figure(
  "Figure_6", "Figure_6_genome_size_pairwise_pgls", GBE_DOUBLE_COLUMN_MM, 105,
  "S04;S19;S24;S26",
  "Estimated genome size, nucleus area, and cell area are positively associated in the released 24-species Pagel-lambda PGLS sensitivity fits. Fits use log10 traits followed by sample standardization; lines and 95% model confidence ribbons are back-transformed onto the displayed measurement scales. Annotations give standardized slopes with 95% intervals and unadjusted P-values. In panels a and b, circles identify balanced multiple-image support, triangles identify three estimates from one image and animal, and crosses identify four species that fail the declared quality-balance gate. Panel c uses a common circle symbol. Bars are conditional measurement-bootstrap intervals. The asterisk marks D. fuscus assembly-anchored pg/1C estimates, not independently measured species C-values; their intervals exclude anchor uncertainty. All three released lambda estimates lie at the lower fitting boundary, so phylogenetic covariance has little influence in these fitted sensitivity models. These associations do not establish causal direction.",
  "Three scatterplots show estimated genome size versus nucleus area, estimated genome size versus cell area, and nucleus area versus cell area, each with measurement intervals and a phylogenetically corrected fitted line; point shape exposes single-image and quality-balance limitations in the genome-size panels.",
  "sensitivity_only_conditional_genome_size"
)

if (identical(figure_options$figure, "Figure_6")) {
  manifest_rows[[1]]$title <- "Estimated genome size and cellular associations"
  figure_manifest <- write_numbered_figure_metadata()
  source(file.path(ROOT, "scripts/publication/figure6_panels.R"))
  export_figure6_panels(list(a = genome_a, b = genome_b, c = genome_c), OUT)
  message("Wrote Figure 6 and its three isolated panels to ", OUT)
  quit(save = "no", status = 0, runLast = FALSE)
}

## ---- figure-7
# Figure 7: estimated genome-size display of the frozen path comparison and calibration.
path_models <- read_supp("Supplementary_Data_S30_path_model_comparison.csv") %>%
  filter(analysis_type == "primary") %>%
  arrange(delta_CICc) %>%
  slice_head(n = 10) %>%
  mutate(
    model_label = gsub("_", " ", model),
    model_label = sub("^genome", "Genome size", model_label),
    model_label = factor(model_label, levels = rev(model_label))
  )
simulation <- read_supp("Supplementary_Data_S36_path_simulation_calibration.csv") %>%
  mutate(
    generating_label = gsub("_", " ", generating_class),
    generating_label = sub("^genome", "Genome size", generating_label),
    residual_label = recode(residual_regime, iid = "Independent residuals", brownian = "Brownian residuals")
  )

path_a <- ggplot(
  path_models,
  aes(
    delta_CICc,
    model_label,
    colour = globally_supported,
    shape = globally_supported
  )
) +
  geom_vline(xintercept = 2, linetype = 2, linewidth = 0.35, colour = GBE_COLORS[["gray"]]) +
  geom_point(size = 1.7) +
  scale_color_manual(values = c(`TRUE` = GBE_COLORS[["blue"]], `FALSE` = GBE_COLORS[["vermillion"]]), name = "Global fit", labels = c("Rejected", "Passed")) +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 17), name = "Global fit", labels = c("Rejected", "Passed")) +
  scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 18)) +
  labs(x = "ΔCICc from best Path24 model", y = NULL) +
  theme_gbe() +
  theme(plot.margin = margin(6, 8, 6, 8, unit = "pt"))

path_b <- ggplot(
  simulation,
  aes(
    generating_label,
    true_class_top_rate,
    colour = residual_label,
    shape = residual_label,
    group = residual_label
  )
) +
  geom_point(size = 1.7, position = position_dodge(width = 0.45)) +
  scale_color_manual(values = unname(GBE_COLORS[c("blue", "orange")]), name = "Simulation regime") +
  scale_shape_manual(values = c(16, 17), name = "Simulation regime") +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1),
                     breaks = c(0, .25, .5, .75, 1)) +
  scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 12)) +
  coord_flip() +
  labs(x = "Generating equivalence class", y = "True-class top-model rate") +
  theme_gbe() +
  theme(plot.margin = margin(6, 8, 6, 8, unit = "pt"))

Figure_7 <- path_a + path_b +
  plot_layout(widths = c(1.2, 1), guides = "collect") +
  plot_annotation(tag_levels = "a") & theme(legend.position = "bottom", legend.box = "vertical")

record_figure(
  "Figure_7", "Figure_7_exploratory_phylogenetic_path_sensitivity", GBE_DOUBLE_COLUMN_MM, 135,
  "S04;S30-S36",
  "The released 24-species estimated genome-size path-model comparison favors the cell-bridge Markov-equivalence class. Panel a shows all ten primary candidate classes, ranked by delta CICc relative to the best model; color and shape record the global d-separation fit gate, and the dashed line marks delta CICc = 2. Panel b reports the frequency with which the generating class is the top-ranked model in 100 simulations for each generating-class and residual-regime combination. The full 0-100% axis makes incomplete and unequal recovery visible. Genome size is the released D. fuscus assembly-anchored pg/1C estimate, not an independently measured species C-value. These exploratory conditional-trait results do not orient causal arrows within an equivalence class.",
  "The left panel compares candidate estimated genome-size equivalence classes in the finalized 24-species analysis by delta CICc and global-fit status. The right panel shows true-class recovery rates under independent and Brownian residual simulations; genome-size labels denote the released D. fuscus assembly-anchored pg/1C scale.",
  "sensitivity_only_conditional_genome_size_no_causal_claim"
)

if (identical(figure_options$figure, "Figure_7")) {
  manifest_rows[[1]]$title <- "Exploratory path-model comparison and calibration"
  figure_manifest <- write_numbered_figure_metadata()
  source(file.path(ROOT, "scripts/publication/figure7_panels.R"))
  export_figure7_panels(list(a = path_a, b = path_b), OUT)
  message("Wrote Figure 7 and its two isolated panels to ", OUT)
  quit(save = "no", status = 0, runLast = FALSE)
}

## ---- wang-style-integrated-te-history-candidate
# Unnumbered review candidate: adapt Wang et al. Figure 1's aligned phylogeny,
# estimated genome-size, repeat-composition, and divergence-profile grammar
# to all 34 TE species, retaining missing genome estimates explicitly.
te_diversity_wide <- te_diversity %>%
  mutate(
    metric_key = recode(
      as.character(metric),
      `Shannon entropy` = "shannon_entropy",
      `Gini-Simpson index` = "gini_simpson"
    )
  ) %>%
  select(species, n_superfamilies_included, metric_key, diversity_value) %>%
  pivot_wider(names_from = metric_key, values_from = diversity_value)

te_genome_diversity <- path_traits %>%
  transmute(
    species = sub("^D\\. ", "", species),
    genome_size_pg_fuscus_anchored,
    genome_size_pg_ci_low,
    genome_size_pg_ci_high
  ) %>%
  inner_join(te_diversity_wide, by = "species", relationship = "one-to-one") %>%
  inner_join(
    te_annotations %>%
      select(
        species,
        reproductive_strategy,
        reproductive_strategy_label
      ),
    by = "species",
    relationship = "one-to-one"
  )

if (
  nrow(te_genome_diversity) != 21L ||
    any(te_genome_diversity$n_superfamilies_included != 10L) ||
    any(!is.finite(as.matrix(te_genome_diversity[, c(
      "genome_size_pg_fuscus_anchored",
      "genome_size_pg_ci_low",
      "genome_size_pg_ci_high",
      "shannon_entropy",
      "gini_simpson"
    )]))) ||
    !identical(
      as.integer(table(factor(
        te_genome_diversity$reproductive_strategy,
        levels = te_strategy_levels
      ))),
      c(6L, 13L, 2L)
    )
) {
  stop("Wang-style candidates require the exact 21-species TE-diversity/genome-size overlap")
}

# Every TE34 species remains visible; unavailable genome estimates stay missing.
integrated_species <- phylogenetic_species_order(te_annotations$species)
integrated_tree <- prune_source_tree(integrated_species)
integrated_tip_count <- Ntip(integrated_tree)
integrated_edges <- reorder.phylo(integrated_tree, "cladewise")$edge
integrated_topology_depth <- numeric(integrated_tip_count + Nnode(integrated_tree))
for (edge_index in seq_len(nrow(integrated_edges))) {
  parent <- integrated_edges[edge_index, 1]
  child <- integrated_edges[edge_index, 2]
  integrated_topology_depth[child] <- integrated_topology_depth[parent] + 1
}
integrated_topology_depth <- integrated_topology_depth / (max(integrated_topology_depth) + 0.5)
integrated_tip_x <- 1
integrated_node_y <- node.height(integrated_tree)
integrated_tree_horizontal <- tibble(
  parent = integrated_tree$edge[, 1], child = integrated_tree$edge[, 2]
) %>% mutate(
  x = integrated_topology_depth[parent],
  xend = if_else(child <= integrated_tip_count, integrated_tip_x, integrated_topology_depth[child]),
  y = integrated_node_y[child], yend = integrated_node_y[child]
)
integrated_tree_vertical <- tibble(parent = sort(unique(integrated_tree$edge[, 1]))) %>%
  rowwise() %>% mutate(
    x = integrated_topology_depth[parent],
    ymin = min(integrated_node_y[integrated_tree$edge[integrated_tree$edge[, 1] == parent, 2]]),
    ymax = max(integrated_node_y[integrated_tree$edge[integrated_tree$edge[, 1] == parent, 2]])
  ) %>% ungroup()
integrated_tree_tips <- tibble(
  species = integrated_tree$tip.label,
  y = integrated_node_y[seq_len(integrated_tip_count)]
) %>%
  left_join(te_genome_diversity, by = "species", relationship = "one-to-one") %>%
  mutate(display_species = species_label(species))
if (
  integrated_tip_count != 34L || Nnode(integrated_tree) != 33L ||
    !setequal(integrated_species, landscape$species) ||
    sum(is.finite(integrated_tree_tips$genome_size_pg_fuscus_anchored)) != 21L ||
    sum(is.na(integrated_tree_tips$genome_size_pg_fuscus_anchored)) != 13L ||
    !identical(integrated_tree_tips$species[order(integrated_tree_tips$y)], integrated_species) ||
    !is.rooted(integrated_tree) || !is.binary(integrated_tree)
) stop("Integrated TE-history candidate requires TE34 in canonical order and 21 observed genome estimates")

integrated_category_levels <- c("LTR", "LINE", "TIR", "DIRS", "SINE", "Maverick", "Other")
integrated_top_orders <- head(integrated_category_levels, -1L)
integrated_major_orders <- c("LTR", "LINE", "TIR")
integrated_category_colors <- c(
  GBE_TE_ORDER_COLORS[integrated_major_orders],
  DIRS = GBE_COLORS[["sky"]], SINE = GBE_COLORS[["yellow"]],
  Maverick = GBE_COLORS[["orange"]], Other = GBE_TE_ORDER_COLORS[["Other"]]
)[integrated_category_levels]
integrated_composition <- landscape %>%
  mutate(repeat_category = if_else(order %in% integrated_top_orders, order, "Other")) %>%
  group_by(species, repeat_category) %>%
  summarise(hit_bp = sum(hit_bp), .groups = "drop") %>%
  right_join(expand_grid(species = integrated_species, repeat_category = integrated_category_levels),
             by = c("species", "repeat_category")) %>%
  mutate(hit_bp = replace_na(hit_bp, 0)) %>%
  group_by(species) %>% mutate(repeat_fraction = hit_bp / sum(hit_bp)) %>% ungroup() %>%
  left_join(integrated_tree_tips %>% select(species, y), by = "species", relationship = "many-to-one") %>%
  mutate(repeat_category = factor(repeat_category, levels = integrated_category_levels)) %>%
  arrange(species, repeat_category) %>%
  group_by(species) %>% mutate(
    composition_end = 100 * cumsum(repeat_fraction),
    composition_start = lag(composition_end, default = 0)
  ) %>% ungroup()
if (
  nrow(integrated_composition) != 34L * 7L || anyNA(integrated_composition$y) ||
    any(abs((integrated_composition %>% group_by(species) %>%
      summarise(total = sum(repeat_fraction), .groups = "drop") %>% pull(total)) - 1) > 1e-12)
) stop("Integrated TE-history composition requires 34 closed seven-category rows")

# Select major orders only after calculating fractions of ALL species hit bp.
integrated_fraction_max <- 60
integrated_fraction_data <- integrated_composition %>%
  filter(repeat_category %in% integrated_major_orders) %>%
  mutate(repeat_category = factor(repeat_category, levels = integrated_major_orders),
         fraction_percent = 100 * repeat_fraction)
if (any(integrated_fraction_data$fraction_percent > integrated_fraction_max)) {
  stop("Integrated repeat fractions exceed the fixed 0-60% scale")
}
integrated_profile_max <- 3
integrated_profiles <- landscape %>%
  filter(order %in% integrated_major_orders) %>%
  transmute(species, repeat_category = order, divergence_bin = divergence_bin_start_pct,
            percent_species_hit_bp) %>%
  right_join(expand_grid(species = integrated_species,
                        repeat_category = integrated_major_orders, divergence_bin = 0:49),
             by = c("species", "repeat_category", "divergence_bin")) %>%
  mutate(percent_species_hit_bp = replace_na(percent_species_hit_bp, 0)) %>%
  left_join(integrated_tree_tips %>% select(species, y), by = "species", relationship = "many-to-one") %>%
  mutate(repeat_category = factor(repeat_category, levels = integrated_major_orders),
         profile_bottom = y - 0.34,
         profile_top = profile_bottom + 0.68 * percent_species_hit_bp / integrated_profile_max)
if (
  nrow(integrated_profiles) != 34L * 3L * 50L || anyNA(integrated_profiles$y) ||
    any(integrated_profiles$percent_species_hit_bp > integrated_profile_max) ||
    any(!is.finite(integrated_profiles$percent_species_hit_bp))
) stop("Integrated profiles require TE34 by three orders by 50 bins on the fixed 0-3% scale")

integrated_row_theme <- theme_gbe() + theme(
  axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank(),
  axis.text.x = element_text(size = GBE_MIN_TEXT_PT),
  axis.title = element_text(size = GBE_MIN_TEXT_PT),
  plot.title = element_text(size = GBE_MIN_TEXT_PT, face = "bold", hjust = 0.5,
                            lineheight = 1.1, margin = margin(b = 6)),
  plot.tag.position = "topleft", plot.tag = element_text(size = 9, face = "bold"),
  strip.background = element_blank(), strip.text = element_text(size = GBE_MIN_TEXT_PT),
  panel.spacing.x = unit(2.5, "mm"), plot.margin = margin(3, 3, 3, 3, unit = "pt")
)
integrated_row_scale <- function() scale_y_continuous(
  limits = c(0.5, integrated_tip_count + 0.5), expand = expansion(mult = 0)
)
integrated_tree_panel <- ggplot() +
  geom_segment(data = integrated_tree_horizontal,
               aes(x = x, xend = xend, y = y, yend = yend),
               linewidth = 0.30, colour = GBE_COLORS[["charcoal"]]) +
  geom_segment(data = integrated_tree_vertical,
               aes(x = x, xend = x, y = ymin, yend = ymax),
               linewidth = 0.30, colour = GBE_COLORS[["charcoal"]]) +
  geom_point(data = integrated_tree_tips, aes(x = integrated_tip_x, y = y),
             shape = 21, size = 1.25, stroke = 0.3, fill = "white") +
  geom_text(data = integrated_tree_tips, aes(x = 1.12, y = y, label = display_species),
            family = GBE_FONT, fontface = "italic", size = GBE_LABEL_SIZE_MM, hjust = 0) +
  scale_x_continuous(limits = c(-0.03, 2.65), breaks = NULL, expand = expansion(mult = 0)) +
  integrated_row_scale() + labs(title = "Species phylogeny\n34 TE species", x = NULL, y = NULL) +
  integrated_row_theme + theme(axis.line.x = element_blank(), axis.ticks.x = element_blank())
integrated_genome_panel <- ggplot() +
  geom_segment(data = filter(integrated_tree_tips, !is.na(genome_size_pg_fuscus_anchored)),
               aes(x = genome_size_pg_ci_low, xend = genome_size_pg_ci_high, y = y, yend = y),
               linewidth = 0.45, colour = GBE_COLORS[["charcoal"]]) +
  geom_point(data = filter(integrated_tree_tips, !is.na(genome_size_pg_fuscus_anchored)),
             aes(x = genome_size_pg_fuscus_anchored, y = y), size = 1.4,
             shape = 21, stroke = 0.3, fill = "white", colour = GBE_COLORS[["charcoal"]]) +
  geom_text(data = filter(integrated_tree_tips, is.na(genome_size_pg_fuscus_anchored)),
            aes(x = 17.5, y = y), label = "n/a", family = GBE_FONT,
            size = GBE_LABEL_SIZE_MM, colour = GBE_COLORS[["gray"]]) +
  scale_x_continuous(limits = c(5, 30), breaks = c(5, 30), expand = expansion(mult = 0.04)) +
  integrated_row_scale() + labs(title = "Estimated\ngenome size", x = "pg/1C*", y = NULL) +
  integrated_row_theme
integrated_composition_panel <- ggplot(integrated_composition) +
  geom_rect(aes(xmin = composition_start, xmax = composition_end, ymin = y - 0.29,
                ymax = y + 0.29, fill = repeat_category), colour = NA) +
  scale_fill_manual(values = integrated_category_colors, drop = FALSE, name = NULL) +
  scale_x_continuous(limits = c(0, 100), breaks = c(0, 100), expand = expansion(mult = 0.04)) +
  integrated_row_scale() + labs(title = "Repeat-hit\ncomposition", x = "% of hit bp", y = NULL) +
  integrated_row_theme + guides(fill = guide_legend(nrow = 1))
integrated_fraction_panel <- ggplot(integrated_fraction_data) +
  geom_segment(aes(x = 0, xend = fraction_percent, y = y, yend = y, colour = repeat_category),
               linewidth = 0.45) +
  geom_point(aes(x = fraction_percent, y = y, colour = repeat_category), size = 1.2) +
  facet_grid(cols = vars(repeat_category)) +
  scale_colour_manual(values = integrated_category_colors, guide = "none") +
  scale_x_continuous(limits = c(0, integrated_fraction_max), breaks = c(0, 60),
                     expand = expansion(mult = 0.12)) +
  integrated_row_scale() + labs(title = "Order fractions\n0-60% per column", x = "% of hit bp", y = NULL) +
  integrated_row_theme
integrated_profile_panel <- ggplot(integrated_profiles) +
  geom_segment(data = distinct(integrated_profiles, species, repeat_category, profile_bottom),
               aes(x = 0, xend = 50, y = profile_bottom, yend = profile_bottom),
               linewidth = 0.18, colour = "#D9D9D9") +
  geom_ribbon(aes(x = divergence_bin, ymin = profile_bottom, ymax = profile_top,
                  group = species, fill = repeat_category), colour = NA) +
  facet_grid(cols = vars(repeat_category)) +
  scale_fill_manual(values = integrated_category_colors, guide = "none") +
  scale_x_continuous(limits = c(0, 50), breaks = c(0, 50), expand = expansion(mult = 0.06)) +
  integrated_row_scale() + labs(title = "Divergence profiles\nHeight: 0-3% of hit bp", x = "Divergence (%)", y = NULL) +
  integrated_row_theme
Wang_style_integrated_TE_history_candidate <-
  (integrated_tree_panel + integrated_genome_panel + integrated_composition_panel +
     integrated_fraction_panel + integrated_profile_panel) +
  plot_layout(nrow = 1, widths = c(2.6, 0.9, 1.1, 1.8, 2.3), guides = "collect") +
  plot_annotation(tag_levels = "a", theme = theme(plot.margin = margin(5, 5, 5, 5))) &
  theme(legend.position = "bottom")
record_auxiliary_figure(
  "Auxiliary_Wang_style_integrated_te_history_candidate",
  "Wang-style integrated TE history across all 34 species",
  Wang_style_integrated_TE_history_candidate,
  "Wang_style_integrated_TE_history_candidate", GBE_DOUBLE_COLUMN_MM, 240,
  "S04;S08;S12;S42;trees/desmo900dated_source46.tre",
  paste(
    "Review candidate, not a numbered manuscript figure. All 34 TE species are included;",
    "estimated genome-size estimates and measurement intervals are available for 21, with 13 labelled n/a.",
    "Estimated genome size is IOD-derived, anchored to D. fuscus, and not a direct species C-value;",
    "intervals omit anchor uncertainty. Repeat proportions are fractions of aligned RepeatMasker hit bp,",
    "not genome occupancy; divergence bins are not insertion ages."
  ),
  paste(
    "Panel a shows the topology-only supplied-tree prune for all 34 TE species; branch lengths are omitted",
    "because publication and calibration provenance remain unresolved. Panel b shows 21 estimated",
    "D. fuscus-anchored genome-size estimates and measurement intervals; n/a marks 13 unavailable estimates.",
    "Panel c shows stacked aligned-hit composition across six named orders plus Other (including Unclassified).",
    "Panel d compares LTR, LINE, and TIR on a common 0-60% fraction scale. Panel e shows those three orders'",
    "zero-filled 0-50% divergence profiles with a common 0-3% height scale. Every repeat track retains the",
    "full species aligned-hit denominator; all tracks use the same phylogenetic order."
  ),
  paste(
    "Thirty-four Desmognathus species align across a topology-only phylogeny, estimated genome-size",
    "estimates with explicit missing values, stacked repeat composition, three repeat-order dot plots,",
    "and three columns of divergence profiles without boxes. All rows share the same phylogenetic order."
  ),
  "auxiliary_review_candidate_wang_style_all34"
)

## ---- wang-style-diversity-genome-size-candidate
# Unnumbered review candidate: reproduce Wang et al. Figure 2's two-row framed
# regression grammar for the exact overlap21, while retaining PGLS and an honest
# no-fit boundary for the two direct developers.
diversity_facet_spec <- tibble(
  facet_id = c(
    "all",
    "aquatic_eggs",
    "terrestrial_eggs_aquatic_larvae",
    "direct_development"
  ),
  facet_label = c(
    "All Desmognathus",
    "Aquatic eggs",
    "Terrestrial eggs + aquatic larvae",
    "Direct development"
  ),
  facet_colour = unname(GBE_COLORS[c("charcoal", "blue", "vermillion", "purple")])
)
diversity_facet_colors <- setNames(
  diversity_facet_spec$facet_colour,
  diversity_facet_spec$facet_label
)
diversity_facet_levels <- diversity_facet_spec$facet_label

diversity_by_facet <- bind_rows(
  te_genome_diversity %>%
    mutate(facet_id = "all", facet_label = "All Desmognathus"),
  te_genome_diversity %>%
    mutate(
      facet_id = reproductive_strategy,
      facet_label = reproductive_strategy_label
    )
) %>%
  mutate(
    facet_label = factor(facet_label, levels = diversity_facet_levels)
  )

expected_diversity_facet_counts <- c(
  all = 21L,
  aquatic_eggs = 6L,
  terrestrial_eggs_aquatic_larvae = 13L,
  direct_development = 2L
)
observed_diversity_facet_counts <- diversity_by_facet %>%
  count(facet_id) %>%
  arrange(match(facet_id, names(expected_diversity_facet_counts)))
if (
  !identical(
    setNames(observed_diversity_facet_counts$n, observed_diversity_facet_counts$facet_id),
    expected_diversity_facet_counts
  )
) {
  stop("Diversity-genome-size facets require exact 21/6/13/2 group sizes")
}

fit_diversity_pgls <- function(data, metric, facet_id, facet_label) {
  current <- data %>% filter(.data$facet_id == .env$facet_id)
  n_species <- nrow(current)
  if (n_species < 3L) {
    return(list(
      prediction = tibble(),
      annotation = tibble(
        metric = metric,
        facet_id = facet_id,
        facet_label = factor(facet_label, levels = diversity_facet_levels),
        label = sprintf("PGLS not fit; n = %d", n_species)
      )
    ))
  }

  current_tree <- prune_source_tree(current$species)
  model_data <- as.data.frame(current, stringsAsFactors = FALSE)
  rownames(model_data) <- model_data$species
  fit <- phylolm(
    reformulate("genome_size_pg_fuscus_anchored", response = metric),
    data = model_data,
    phy = current_tree,
    model = "lambda",
    lower.bound = 0,
    upper.bound = 1
  )
  coefficient_table <- summary(fit)$coefficients
  slope_t <- coefficient_table["genome_size_pg_fuscus_anchored", "t.value"]
  slope_p <- coefficient_table["genome_size_pg_fuscus_anchored", "p.value"]
  if (
    any(!is.finite(coefficient_table)) ||
      !is.finite(fit$optpar) ||
      fit$optpar < 0 ||
      fit$optpar > 1
  ) {
    stop("Diversity-genome-size PGLS produced a non-finite or invalid fit")
  }

  x_sequence <- seq(
    min(current$genome_size_pg_fuscus_anchored),
    max(current$genome_size_pg_fuscus_anchored),
    length.out = 160
  )
  design <- cbind(
    `(Intercept)` = 1,
    genome_size_pg_fuscus_anchored = x_sequence
  )
  fit_values <- as.numeric(design %*% stats::coef(fit))
  fit_variance <- rowSums((design %*% stats::vcov(fit)) * design)
  critical <- stats::qt(0.975, df = n_species - 2L)
  fit_se <- sqrt(pmax(fit_variance, 0))
  p_label <- if_else(
    slope_p < 0.001,
    "< 0.001",
    paste("=", sprintf("%.3f", slope_p))
  )

  list(
    prediction = tibble(
      metric = metric,
      facet_id = facet_id,
      facet_label = factor(facet_label, levels = diversity_facet_levels),
      genome_size_pg_fuscus_anchored = x_sequence,
      fit = fit_values,
      low = fit_values - critical * fit_se,
      high = fit_values + critical * fit_se
    ),
    annotation = tibble(
      metric = metric,
      facet_id = facet_id,
      facet_label = factor(facet_label, levels = diversity_facet_levels),
      label = sprintf(
        "PGLS F = %.2f, P %s\nlambda = %.2f; n = %d",
        slope_t^2,
        p_label,
        fit$optpar,
        n_species
      )
    )
  )
}

diversity_fit_results <- lapply(
  c("shannon_entropy", "gini_simpson"),
  function(metric) {
    lapply(seq_len(nrow(diversity_facet_spec)), function(index) {
      fit_diversity_pgls(
        diversity_by_facet,
        metric,
        diversity_facet_spec$facet_id[[index]],
        diversity_facet_spec$facet_label[[index]]
      )
    })
  }
) %>%
  unlist(recursive = FALSE)
diversity_fit_predictions <- bind_rows(lapply(
  diversity_fit_results,
  `[[`,
  "prediction"
))
diversity_fit_annotations <- bind_rows(lapply(
  diversity_fit_results,
  `[[`,
  "annotation"
))
if (
  nrow(diversity_fit_annotations) != 8L ||
    nrow(diversity_fit_predictions) != 6L * 160L ||
    sum(grepl("not fit", diversity_fit_annotations$label, fixed = TRUE)) != 2L
) {
  stop("Diversity-genome-size candidate requires six PGLS fits and two n=2 no-fit panels")
}

wang_diversity_panel <- function(metric, y_label, tag, show_x_title) {
  point_data <- diversity_by_facet %>%
    transmute(
      species,
      facet_id,
      facet_label,
      genome_size_pg_fuscus_anchored,
      genome_size_pg_ci_low,
      genome_size_pg_ci_high,
      diversity_value = .data[[metric]]
    )
  prediction_data <- diversity_fit_predictions %>%
    filter(.data$metric == .env$metric)
  annotation_data <- diversity_fit_annotations %>%
    filter(.data$metric == .env$metric)

  ggplot() +
    geom_ribbon(
      data = prediction_data,
      aes(
        x = genome_size_pg_fuscus_anchored,
        ymin = low,
        ymax = high,
        fill = facet_label
      ),
      alpha = 0.24,
      colour = NA,
      show.legend = FALSE
    ) +
    geom_line(
      data = prediction_data,
      aes(
        x = genome_size_pg_fuscus_anchored,
        y = fit,
        colour = facet_label
      ),
      linewidth = 0.65,
      show.legend = FALSE
    ) +
    geom_segment(
      data = point_data,
      aes(
        x = genome_size_pg_ci_low,
        xend = genome_size_pg_ci_high,
        y = diversity_value,
        yend = diversity_value,
        colour = facet_label
      ),
      linewidth = 0.3,
      alpha = 0.34,
      show.legend = FALSE
    ) +
    geom_point(
      data = point_data,
      aes(
        x = genome_size_pg_fuscus_anchored,
        y = diversity_value,
        colour = facet_label
      ),
      shape = 21,
      size = 1.65,
      stroke = 0.45,
      fill = "white",
      show.legend = FALSE
    ) +
    geom_text(
      data = annotation_data,
      aes(x = -Inf, y = -Inf, label = label, colour = facet_label),
      hjust = -0.04,
      vjust = -0.35,
      size = GBE_LABEL_SIZE_MM,
      lineheight = 0.95,
      family = GBE_FONT,
      show.legend = FALSE
    ) +
    facet_grid(cols = vars(facet_label), drop = FALSE, labeller = label_wrap_gen(21)) +
    scale_colour_manual(values = diversity_facet_colors) +
    scale_fill_manual(values = diversity_facet_colors) +
    scale_y_continuous(expand = expansion(mult = c(0.22, 0.05))) +
    scale_x_continuous(
      limits = range(te_genome_diversity$genome_size_pg_ci_low,
                     te_genome_diversity$genome_size_pg_ci_high),
      breaks = scales::breaks_pretty(n = 4),
      expand = expansion(mult = c(0.025, 0.025))
    ) +
    labs(
      x = if (show_x_title) "Estimated genome size (pg/1C*)" else NULL,
      y = y_label,
      tag = tag
    ) +
    coord_cartesian(clip = "off") +
    theme_gbe(base_size = 7.5) +
    theme(
      panel.spacing.x = grid::unit(4, "pt"),
      strip.text = element_text(size = GBE_MIN_TEXT_PT, face = "bold"),
      plot.margin = margin(4, 4, 4, 4, unit = "pt")
    )
}

wang_diversity_shannon <- wang_diversity_panel(
  "shannon_entropy",
  "Shannon entropy",
  "a",
  FALSE
)
wang_diversity_gini <- wang_diversity_panel(
  "gini_simpson",
  "Gini-Simpson index",
  "b",
  TRUE
)
Wang_style_diversity_genome_size_candidate <-
  wang_diversity_shannon / wang_diversity_gini +
    plot_layout(heights = c(1, 1))

record_auxiliary_figure(
  "Auxiliary_Wang_style_diversity_genome_size_candidate",
  "Wang-style TE diversity and genome-size candidate",
  Wang_style_diversity_genome_size_candidate,
  "Wang_style_TE_diversity_genome_size_candidate",
  GBE_DOUBLE_COLUMN_MM,
  145,
  "S04;S08;S42;trees/desmo900dated_source46.tre",
  paste(
    "Review candidate, not a numbered manuscript figure. All genome-size values are",
    "D. fuscus assembly-anchored IOD-derived estimates, not direct species",
    "C-values. Reproductive-strategy panels are post-hoc descriptive strata. The two",
    "direct developers are shown without a regression because n = 2."
  ),
  paste(
    "Rows show the canonical top-ten-superfamily Shannon entropy and Gini-Simpson index.",
    "Columns show all 21 overlap species, six aquatic-egg species, 13 terrestrial-egg",
    "species with aquatic larvae, and two direct developers. Hollow points are species",
    "and horizontal bars are estimated genome-size intervals. Lines and bands are",
    "Pagel-lambda PGLS fits and coefficient confidence bands; the n = 2 direct-development",
    "panels intentionally have no fit."
  ),
  paste(
    "Two aligned rows of four framed scatterplots show Shannon and Gini-Simpson diversity",
    "against estimated genome size for all Desmognathus and three reproductive strategies.",
    "PGLS statistics appear inside six panels; the two direct-development panels show only",
    "their two observations and state that a fit was not estimated."
  ),
  "auxiliary_review_candidate_wang_style_pgls_post_hoc_strata"
)

## ---- path24-development-traits-reference
# Descriptive Path24 reference: topology and three released phenotype estimates.
# Development mode is retained in the caption rather than a redundant track.
organismal_traits <- read_csv(
  file.path(
    ROOT, "analyses", "03_phylogenetic_path", "data", "organismal_traits.csv"
  ),
  show_col_types = FALSE,
  progress = FALSE
)
genome_size_path24 <- path_traits %>%
  select(
    species,
    genome_size_pg_fuscus_anchored,
    genome_size_pg_ci_low,
    genome_size_pg_ci_high
  )

path24_development_traits <- tree_tips %>%
  select(species, tree_tip, y) %>%
  left_join(
    morphology %>%
      transmute(
        species = as.character(species),
        cell_area_um2,
        cell_area_um2_ci_low,
        cell_area_um2_ci_high,
        nucleus_area_um2,
        nucleus_area_um2_ci_low,
        nucleus_area_um2_ci_high
      ),
    by = "species",
    relationship = "one-to-one"
  ) %>%
  left_join(
    genome_size_path24,
    by = "species",
    relationship = "one-to-one"
  ) %>%
  left_join(
    organismal_traits %>%
      select(tree_tip = species, development_mode),
    by = "tree_tip",
    relationship = "one-to-one"
  )

trait_value_columns <- c(
  "genome_size_pg_fuscus_anchored",
  "genome_size_pg_ci_low",
  "genome_size_pg_ci_high",
  "nucleus_area_um2",
  "nucleus_area_um2_ci_low",
  "nucleus_area_um2_ci_high",
  "cell_area_um2",
  "cell_area_um2_ci_low",
  "cell_area_um2_ci_high"
)
path24_direct_tips <- path24_development_traits %>%
  filter(development_mode == "direct_development") %>%
  pull(tree_tip)
all_direct_tips <- organismal_traits %>%
  filter(development_mode == "direct_development") %>%
  pull(species)
if (
  nrow(path24_development_traits) != 24L ||
    n_distinct(path24_development_traits$tree_tip) != 24L ||
    any(is.na(path24_development_traits$development_mode)) ||
    any(!is.finite(as.matrix(path24_development_traits[, trait_value_columns]))) ||
    any(path24_development_traits$genome_size_pg_ci_low > path24_development_traits$genome_size_pg_fuscus_anchored) ||
    any(path24_development_traits$genome_size_pg_ci_high < path24_development_traits$genome_size_pg_fuscus_anchored) ||
    any(path24_development_traits$nucleus_area_um2_ci_low > path24_development_traits$nucleus_area_um2) ||
    any(path24_development_traits$nucleus_area_um2_ci_high < path24_development_traits$nucleus_area_um2) ||
    any(path24_development_traits$cell_area_um2_ci_low > path24_development_traits$cell_area_um2) ||
    any(path24_development_traits$cell_area_um2_ci_high < path24_development_traits$cell_area_um2) ||
    !setequal(path24_direct_tips, c("aeneus", "wrighti")) ||
    !setequal(all_direct_tips, c("aeneus", "organi", "wrighti")) ||
    "organi" %in% path24_development_traits$tree_tip ||
    sum(path24_development_traits$development_mode == "aquatic_larva") != 22L
) {
  stop("Path24 development-trait reference requires the exact 24/22/2 audited contract")
}

path24_tree_reference <- morph_tree_panel +
  labs(title = "Phylogeny (topology only)", tag = "a") +
  guides(fill = "none", shape = "none") +
  theme(
    plot.margin = margin(3, 1, 3, 4, unit = "pt")
  )
path24_tree_reference$layers <- Filter(
  function(layer) !inherits(layer$geom, "GeomPoint"), path24_tree_reference$layers
)

aligned_trait_panel <- function(
  data, estimate, low, high, x_label, colour, x_limits, x_breaks
) {
  if (any(data[[low]] < x_limits[[1]]) || any(data[[high]] > x_limits[[2]])) {
    stop("The displayed trait range must retain every released interval endpoint")
  }
  ggplot(data, aes(x = .data[[estimate]], y = y)) +
    geom_errorbar(
      aes(xmin = .data[[low]], xmax = .data[[high]]),
      orientation = "y",
      width = 0,
      linewidth = 0.4,
      colour = colour,
      alpha = 0.72
    ) +
    geom_point(
      shape = 21,
      size = 1.8,
      stroke = 0.4,
      fill = colour,
      colour = GBE_COLORS[["charcoal"]]
    ) +
    geom_text(
      aes(
        label = sprintf("%.1f", .data[[estimate]]),
        hjust = ifelse(
          .data[[estimate]] > x_limits[[2]] - 0.15 * diff(x_limits), 1,
          ifelse(.data[[estimate]] < x_limits[[1]] + 0.15 * diff(x_limits), 0, 0.5)
        )
      ),
      nudge_y = 0.38,
      family = GBE_FONT,
      size = GBE_LABEL_SIZE_MM,
      colour = GBE_COLORS[["charcoal"]]
    ) +
    scale_x_continuous(limits = x_limits, breaks = x_breaks,
                       expand = expansion(mult = 0.04)) +
    scale_y_continuous(limits = c(0.5, tip_count + 0.5), expand = c(0, 0)) +
    coord_cartesian(clip = "off") +
    labs(x = x_label, y = NULL) +
    theme_gbe(base_size = 7.5) +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.line.y = element_blank(),
      plot.margin = margin(3, 6, 3, 6, unit = "pt")
    )
}

path24_genome_reference <- aligned_trait_panel(
  path24_development_traits,
  "genome_size_pg_fuscus_anchored",
  "genome_size_pg_ci_low",
  "genome_size_pg_ci_high",
  "Estimated genome size\n(pg/1C*)",
  GBE_COLORS[["blue"]],
  c(8, 27), c(8, 12, 16, 20, 24, 27)
) + labs(tag = "b")
path24_nucleus_reference <- aligned_trait_panel(
  path24_development_traits,
  "nucleus_area_um2",
  "nucleus_area_um2_ci_low",
  "nucleus_area_um2_ci_high",
  "Nucleus area\n(µm²)",
  GBE_COLORS[["orange"]],
  c(35, 55), seq(35, 55, 5)
) + labs(tag = "c")
path24_cell_reference <- aligned_trait_panel(
  path24_development_traits,
  "cell_area_um2",
  "cell_area_um2_ci_low",
  "cell_area_um2_ci_high",
  "Cell area\n(µm²)",
  GBE_COLORS[["purple"]],
  c(190, 335), seq(200, 325, 25)
) + labs(tag = "d")

Path24_development_traits_reference <-
  path24_tree_reference +
    path24_genome_reference +
    path24_nucleus_reference +
    path24_cell_reference +
    plot_layout(widths = c(2, 1, 1, 1), guides = "collect") +
    plot_annotation(
      caption = paste0(
        "Development: 22 species have aquatic larvae; D. aeneus and D. wrighti are direct developers.\n",
        "D. organi is absent. Numbers label species estimates; bars show released measurement intervals.\n",
        "Genome-size values are D. fuscus assembly-anchored pg/1C estimates.\n",
        "Tree branch lengths are omitted."
      )
    ) &
  theme(
    legend.position = "bottom",
    plot.caption = element_text(
      size = GBE_MIN_TEXT_PT,
      hjust = 0,
      colour = GBE_COLORS[["gray"]]
    )
  )

record_auxiliary_figure(
  "Auxiliary_Path24_development_traits_reference",
  "Path24 phylogeny and phenotype reference",
  Path24_development_traits_reference,
  "Path24_development_and_trait_phylogeny_reference",
  GBE_DOUBLE_COLUMN_MM,
  190,
  "S04;S19;S38;S39;organismal_traits",
  paste(
    "Descriptive review only. D. organi is absent from Path24 and is not imputed.",
    "The genome-size track uses D. fuscus assembly-anchored pg/1C estimates,",
    "not direct independently measured species C-values. The tree shows topology only;",
    "publication and calibration provenance remain unresolved."
  ),
  paste(
    "The supplied 24-tip Path24 tree (a) is shown without branch lengths, aligned to estimated",
    "D. fuscus-anchored genome-size (b), reviewed nucleus-area (c), and reviewed cell-area (d) estimates.",
    "Numbers label species estimates to one decimal; horizontal bars retain all released conditional",
    "measurement intervals. Axes are tightened to 8-27 pg/1C, 35-55 square micrometres, and",
    "190-335 square micrometres respectively. The two sampled direct",
    "developers are D. aeneus and D. wrighti; the other 22 Path24 species have aquatic larvae",
    "and undergo metamorphosis. This is a descriptive visualization, not a development-mode test."
  ),
  paste(
    "A topology-only phylogeny for 24 Desmognathus species is aligned with labelled",
    "forest-plot tracks for estimated genome size, nucleus area, and cell area.",
    "D. aeneus and D. wrighti are the two direct developers represented in Path24."
  ),
  "auxiliary_review_only_descriptive"
)

## ---- path-dag-reference
# Unnumbered reference: one representative DAG for each scored equivalence class.
candidate_dags <- read_csv(
  file.path(
    ROOT, "analyses", "03_phylogenetic_path", "data", "candidate_dags.csv"
  ),
  show_col_types = FALSE,
  progress = FALSE
)
path_reference_models <- read_supp(
  "Supplementary_Data_S30_path_model_comparison.csv"
) %>%
  filter(analysis_type == "primary")

dag_representatives <- tribble(
  ~model, ~display_name, ~representative_dag, ~equivalence_class, ~class_members,
  "independent", "Independent", "no_edges", "independent", 1L,
  "genome_nucleus_only", "Genome-nucleus only", "GS->NS", "genome_nucleus_only", 2L,
  "nucleus_cell_only", "Nucleus-cell only", "NS->CS", "nucleus_cell_only", 2L,
  "genome_cell_only", "Genome-cell only", "GS->CS", "genome_cell_only", 2L,
  "nucleus_bridge", "Nucleus bridge", "GS->NS;NS->CS", "nucleus_bridge", 3L,
  "nucleus_collider", "Nucleus collider", "GS->NS;CS->NS", "nucleus_collider", 1L,
  "cell_bridge", "Cell bridge", "GS->CS;CS->NS", "cell_bridge", 3L,
  "cell_collider", "Cell collider", "GS->CS;NS->CS", "cell_collider", 1L,
  "genome_bridge", "Genome bridge", "GS->CS;NS->GS", "genome_bridge", 3L,
  "genome_collider", "Genome collider", "NS->GS;CS->GS", "genome_collider", 1L
) %>%
  mutate(model_order = row_number())

if (
  nrow(candidate_dags) != 25L ||
    n_distinct(candidate_dags$equivalence_class) != 11L ||
    n_distinct(candidate_dags$equivalence_class[candidate_dags$testable_by_dsep]) != 10L ||
    nrow(path_reference_models) != 10L ||
    !setequal(path_reference_models$model, dag_representatives$model)
) {
  stop("Path-DAG reference requires the frozen 25-DAG, 11-class, 10-testable-class contract")
}

representative_audit <- dag_representatives %>%
  left_join(
    candidate_dags %>%
      select(
        representative_dag = dag_id,
        candidate_class = equivalence_class,
        candidate_members = class_member_count,
        testable_by_dsep
      ),
    by = "representative_dag",
    relationship = "one-to-one"
  )
if (
  any(is.na(representative_audit$candidate_class)) ||
    any(representative_audit$candidate_class != representative_audit$equivalence_class) ||
    any(representative_audit$candidate_members != representative_audit$class_members) ||
    !all(representative_audit$testable_by_dsep)
) {
  stop("Path-DAG representative mapping drifted from candidate_dags.csv")
}

representative_audit <- representative_audit %>%
  mutate(
    panel_label = sprintf(
      "%s\n(%d equivalent %s)",
      display_name,
      class_members,
      if_else(class_members == 1L, "DAG", "DAGs")
    ),
    panel_label = factor(panel_label, levels = panel_label[order(model_order)])
  )

dag_nodes <- crossing(
  representative_audit %>% select(model, panel_label),
  tibble(
    node = c("genome_size", "nucleus_size", "cell_size"),
    node_label = c("Genome\nsize", "Nucleus\narea", "Cell\narea"),
    x = c(0.15, 0.85, 0.50),
    y = c(0.23, 0.23, 0.78),
    text_colour = c(
      GBE_COLORS[["charcoal"]],
      GBE_COLORS[["charcoal"]],
      "white"
    )
  )
)

segment_templates <- tribble(
  ~edge_key, ~x, ~y, ~xend, ~yend,
  "genome_size->nucleus_size", 0.26, 0.23, 0.74, 0.23,
  "nucleus_size->genome_size", 0.74, 0.23, 0.26, 0.23,
  "genome_size->cell_size", 0.23, 0.34, 0.44, 0.67,
  "cell_size->genome_size", 0.44, 0.67, 0.23, 0.34,
  "nucleus_size->cell_size", 0.77, 0.34, 0.56, 0.67,
  "cell_size->nucleus_size", 0.56, 0.67, 0.77, 0.34
)

dag_edges <- representative_audit %>%
  filter(representative_dag != "no_edges") %>%
  separate_rows(representative_dag, sep = ";") %>%
  separate(representative_dag, into = c("parent_code", "child_code"), sep = "->") %>%
  mutate(
    parent = recode(parent_code, GS = "genome_size", NS = "nucleus_size", CS = "cell_size"),
    child = recode(child_code, GS = "genome_size", NS = "nucleus_size", CS = "cell_size"),
    edge_key = paste(parent, child, sep = "->")
  ) %>%
  left_join(segment_templates, by = "edge_key", relationship = "many-to-one")
if (any(!complete.cases(dag_edges[, c("x", "y", "xend", "yend")]))) {
  stop("Path-DAG reference contains an unsupported directed edge")
}

Path_DAG_reference <- ggplot() +
  geom_segment(
    data = dag_edges,
    aes(x = x, y = y, xend = xend, yend = yend),
    linewidth = 0.65,
    colour = GBE_COLORS[["charcoal"]],
    arrow = grid::arrow(type = "closed", angle = 22, length = grid::unit(2.2, "mm"))
  ) +
  geom_point(
    data = dag_nodes,
    aes(x = x, y = y, fill = node),
    shape = 21,
    size = 15.5,
    stroke = 0.65,
    colour = GBE_COLORS[["charcoal"]]
  ) +
  geom_text(
    data = dag_nodes,
    aes(x = x, y = y, label = node_label, colour = text_colour),
    size = GBE_LABEL_SIZE_MM,
    family = GBE_FONT,
    fontface = "bold"
  ) +
  facet_wrap(vars(panel_label), ncol = 2, drop = FALSE) +
  scale_fill_manual(
    values = c(
      genome_size = GBE_COLORS[["blue"]],
      nucleus_size = GBE_COLORS[["orange"]],
      cell_size = GBE_COLORS[["purple"]]
    ),
    guide = "none"
  ) +
  scale_colour_identity() +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  labs(
    caption = paste0(
      "Genome size is the D. fuscus assembly-anchored pg/1C estimate, not a direct species C-value.\n",
      "Arrows are not unique within classes; saturated class (6 members) unscored."
    )
  ) +
  theme_void(base_size = 8, base_family = resolve_gbe_font()) +
  theme(
    panel.border = element_rect(
      fill = NA,
      linewidth = 0.4,
      colour = GBE_COLORS[["charcoal"]]
    ),
    panel.spacing = grid::unit(5, "pt"),
    strip.background = element_rect(
      fill = "white",
      linewidth = 0.35,
      colour = GBE_COLORS[["charcoal"]]
    ),
    strip.text = element_text(
      size = 8,
      face = "bold",
      lineheight = 0.95,
      margin = margin(3, 2, 3, 2, unit = "pt")
    ),
    plot.caption = element_text(
      size = 7,
      hjust = 0,
      lineheight = 1.05,
      colour = GBE_COLORS[["gray"]],
      margin = margin(t = 5, unit = "pt")
    ),
    plot.margin = margin(4, 4, 4, 4, unit = "pt"),
    aspect.ratio = 0.47
  )

record_auxiliary_figure(
  "Auxiliary_Path_DAG_reference",
  "Path DAG equivalence-class reference",
  Path_DAG_reference,
  "Path_DAG_equivalence_class_reference",
  GBE_DOUBLE_COLUMN_MM,
  220,
  "S04;S30;candidate_dags",
  paste(
    "Not a numbered manuscript figure. The DAG node is Genome size, expressed on the",
    "D. fuscus assembly-anchored pg/1C scale, not a direct species C-value measurement.",
    "The same released genome-size scale is used consistently across all scored classes.",
    "The 10 testable equivalence classes do not orient a unique causal arrow; the saturated",
    "6-member class is not shown and remains unscored because it has no d-separation claim."
  ),
  paste(
    "One representative DAG is shown for each of the 10 testable Path24 equivalence classes.",
    "Genome size denotes the D. fuscus assembly-anchored pg/1C estimate, not a",
    "direct species C-value measurement. Equivalent arrow",
    "orientations are omitted, and the saturated 6-member class remains unscored."
  ),
  paste(
    "Ten boxed path diagrams connect Genome size, Nucleus area, and Cell area with arrows.",
    "Each is one representative of a testable equivalence class, not a uniquely oriented",
    "causal graph; the six-member saturated class is omitted as unscored."
  ),
  "auxiliary_review_only_conditional_genome_size_no_causal_claim"
)

## ---- wang2025-cross-study-diversity-reference
# Unnumbered descriptive reference: exact Wang et al. 2025 Dryad Figure 2
# values beside current top-10-superfamily Desmognathus summaries. The two
# studies remain method-labelled and are never pooled into a fitted model.
cross_study_diversity <- read_supp("Supplementary_Data_S41_wang2025_cross_study_te_diversity.csv")
expected_cross_study_groups <- c(
  "Birds",
  "Bony fishes",
  "Cartilaginous fishes",
  "Lobe-finned fishes",
  "Mammals",
  "Non-avian reptiles",
  "Frogs",
  "Published salamanders",
  "Current Desmognathus"
)
expected_cross_study_counts <- c(
  wang2025_published = 84L,
  desmognathus_current = 34L
)
expected_inclusion_counts <- c(
  published_group_only = 73L,
  published_scatter_and_group = 11L,
  current_group_only = 13L,
  current_scatter_and_group = 21L
)
if (
  nrow(cross_study_diversity) != 118L ||
    !identical(
      as.integer(
        table(factor(
          cross_study_diversity$data_origin,
          levels = names(expected_cross_study_counts)
        ))
      ),
      unname(expected_cross_study_counts)
    ) ||
    !identical(
      as.integer(
        table(factor(
          cross_study_diversity$figure_inclusion,
          levels = names(expected_inclusion_counts)
        ))
      ),
      unname(expected_inclusion_counts)
    ) ||
    !setequal(cross_study_diversity$comparison_group, expected_cross_study_groups) ||
    any(!is.finite(cross_study_diversity$shannon_top10_superfamilies)) ||
    any(!is.finite(cross_study_diversity$gini_simpson_top10_superfamilies))
) {
  stop("S41 must remain the exact 84-published plus 34-current descriptive release")
}

cross_study_scatter <- cross_study_diversity %>%
  filter(grepl("scatter_and_group$", figure_inclusion)) %>%
  mutate(
    study_label = factor(
      if_else(
        data_origin == "desmognathus_current",
        "Current Desmognathus",
        "Wang et al. 2025 salamanders"
      ),
      levels = c("Wang et al. 2025 salamanders", "Current Desmognathus")
    )
  )
if (
  nrow(cross_study_scatter) != 32L ||
    sum(cross_study_scatter$data_origin == "wang2025_published") != 11L ||
    sum(cross_study_scatter$data_origin == "desmognathus_current") != 21L ||
    any(!is.finite(cross_study_scatter$genome_size_gb)) ||
    any(!is.finite(
      cross_study_scatter %>%
        filter(data_origin == "desmognathus_current") %>%
        select(genome_size_ci_low_gb, genome_size_ci_high_gb) %>%
        unlist(use.names = FALSE)
    ))
) {
  stop("S41 scatter subset must remain 11 published salamanders plus overlap21")
}

group_counts <- cross_study_diversity %>%
  count(comparison_group, name = "group_n") %>%
  mutate(group_label = paste0(comparison_group, " (n = ", group_n, ")"))
cross_study_groups <- cross_study_diversity %>%
  left_join(group_counts, by = "comparison_group", relationship = "many-to-one") %>%
  mutate(
    group_label = factor(
      group_label,
      levels = rev(
        group_counts$group_label[
          match(expected_cross_study_groups, group_counts$comparison_group)
        ]
      )
    ),
    group_emphasis = factor(
      case_when(
        comparison_group == "Current Desmognathus" ~ "Current Desmognathus",
        comparison_group == "Published salamanders" ~ "Published salamanders",
        TRUE ~ "Other Wang et al. groups"
      ),
      levels = c(
        "Other Wang et al. groups",
        "Published salamanders",
        "Current Desmognathus"
      )
    )
  )

cross_study_scatter_panel <- function(data, metric, y_label) {
  ggplot(data, aes(x = genome_size_gb, y = .data[[metric]])) +
    geom_segment(
      data = filter(data, data_origin == "desmognathus_current"),
      aes(
        x = genome_size_ci_low_gb,
        xend = genome_size_ci_high_gb,
        y = .data[[metric]],
        yend = .data[[metric]],
        colour = study_label
      ),
      linewidth = 0.35,
      alpha = 0.48,
      show.legend = FALSE
    ) +
    geom_point(
      aes(fill = study_label, shape = study_label),
      size = 1.8,
      stroke = 0.45,
      colour = GBE_COLORS[["charcoal"]]
    ) +
    scale_fill_manual(
      values = c(
        "Wang et al. 2025 salamanders" = "white",
        "Current Desmognathus" = GBE_COLORS[["orange"]]
      ),
      name = NULL
    ) +
    scale_shape_manual(
      values = c(
        "Wang et al. 2025 salamanders" = 21,
        "Current Desmognathus" = 24
      ),
      name = NULL
    ) +
    scale_colour_manual(
      values = c(
        "Wang et al. 2025 salamanders" = GBE_COLORS[["sky"]],
        "Current Desmognathus" = GBE_COLORS[["orange"]]
      ),
      guide = "none"
    ) +
    scale_x_continuous(breaks = scales::breaks_pretty(n = 5)) +
    labs(
      x = "Estimated genome size (Gbp)",
      y = y_label
    ) +
    theme_gbe(base_size = 7.5)
}

cross_study_group_panel <- function(data, metric, x_label) {
  ggplot(data, aes(x = .data[[metric]], y = group_label)) +
    geom_boxplot(
      aes(fill = group_emphasis),
      width = 0.58,
      outlier.shape = NA,
      linewidth = 0.35,
      colour = GBE_COLORS[["charcoal"]],
      show.legend = FALSE
    ) +
    geom_point(
      aes(fill = group_emphasis),
      position = position_jitter(width = 0, height = 0.12, seed = 20260812),
      shape = 21,
      size = 0.8,
      stroke = 0.2,
      colour = GBE_COLORS[["charcoal"]],
      alpha = 0.68,
      show.legend = FALSE
    ) +
    scale_fill_manual(
      values = c(
        "Other Wang et al. groups" = GBE_COLORS[["light_gray"]],
        "Published salamanders" = GBE_COLORS[["sky"]],
        "Current Desmognathus" = GBE_COLORS[["orange"]]
      )
    ) +
    labs(x = x_label, y = NULL) +
    theme_gbe(base_size = 7.5) +
    theme(axis.text.y = element_text(size = GBE_MIN_TEXT_PT))
}

wang_a <- cross_study_scatter_panel(
  cross_study_scatter,
  "shannon_top10_superfamilies",
  "Shannon entropy\n(top 10 superfamilies)"
)
wang_b <- cross_study_group_panel(
  cross_study_groups,
  "shannon_top10_superfamilies",
  "Shannon entropy (top 10 superfamilies)"
)
wang_c <- cross_study_scatter_panel(
  cross_study_scatter,
  "gini_simpson_top10_superfamilies",
  "Gini-Simpson index\n(top 10 superfamilies)"
)
wang_d <- cross_study_group_panel(
  cross_study_groups,
  "gini_simpson_top10_superfamilies",
  "Gini-Simpson index (top 10 superfamilies)"
)

Wang2025_cross_study_TE_diversity_reference <-
  wang_a + wang_b + wang_c + wang_d +
    plot_layout(design = "AB\nCD", widths = c(1, 1.32), guides = "collect") +
    plot_annotation(
      tag_levels = "a",
      caption = paste0(
        "Exact Wang et al. (2025) CC0 Dryad Figure 2 values; current values were recomputed from dnaPipeTE.\n",
        "Pipelines are not harmonized and no pooled fit is shown."
      )
    ) &
  theme(
    legend.position = "bottom",
    plot.caption = element_text(
      size = GBE_MIN_TEXT_PT,
      hjust = 0,
      colour = GBE_COLORS[["gray"]]
    )
  )

record_auxiliary_figure(
  "Auxiliary_Wang2025_cross_study_te_diversity_reference",
  "Cross-study TE-diversity reference (Wang et al. 2025)",
  Wang2025_cross_study_TE_diversity_reference,
  "Wang2025_cross_study_TE_diversity_reference",
  GBE_DOUBLE_COLUMN_MM,
  158,
  "S41",
  paste(
    "Descriptive cross-study context only; no pooled regression or hypothesis test is shown.",
    "Wang et al. 2025 published top-10 TE-superfamily indices are exact CC0-1.0 Dryad",
    "values (dataset DOI 10.5061/dryad.zpc866tkv; article DOI 10.1038/s42003-025-08127-3)",
    "produced from RepeatModeler/RepeatMasker annotations. Current top-10 indices are reclosed",
    "calculations from dnaPipeTE compositions. The pipelines are different and not method-harmonized.",
    "Current genome sizes are D. fuscus assembly-anchored estimates, not direct species C-values."
  ),
  paste(
    "Panels a and c show genome size against top-ten-superfamily Shannon entropy and the",
    "Gini-Simpson index for 11 published salamanders and 21 current Desmognathus species.",
    "Panels b and d show the exact published vertebrate-group distributions beside all 34",
    "current Desmognathus summaries. Wang et al. values come from RepeatModeler/RepeatMasker;",
    "current values come from dnaPipeTE. The comparison is descriptive and has no pooled fit."
  ),
  paste(
    "Two scatterplots compare genome size with top-ten-superfamily Shannon and Gini-Simpson",
    "diversity for published salamanders and current Desmognathus. Two aligned box-and-point",
    "panels compare nine vertebrate groups, highlighting published salamanders in blue and",
    "current Desmognathus in orange. No regression line or p-value is shown."
  ),
  "auxiliary_review_only_cross_study_context_non_harmonized"
)

## ---- figure-manifest
figure_manifest <- write_numbered_figure_metadata()
auxiliary_figure_manifest <- bind_rows(auxiliary_manifest_rows)
if (
  nrow(auxiliary_figure_manifest) != 5L ||
    !identical(
      auxiliary_figure_manifest$auxiliary_id,
      c(
        "Auxiliary_Wang_style_integrated_te_history_candidate",
        "Auxiliary_Wang_style_diversity_genome_size_candidate",
        "Auxiliary_Path24_development_traits_reference",
        "Auxiliary_Path_DAG_reference",
        "Auxiliary_Wang2025_cross_study_te_diversity_reference"
      )
    )
) {
  stop("Auxiliary figure manifest requires two Wang-style candidates followed by the three reference figures")
}
write_csv(
  auxiliary_figure_manifest,
  file.path(OUT, "AUXILIARY_FIGURE_MANIFEST.csv")
)
message(
  "Wrote ", nrow(figure_manifest), " numbered GBE figures and ",
  nrow(auxiliary_figure_manifest), " auxiliary review figures to ", OUT
)
