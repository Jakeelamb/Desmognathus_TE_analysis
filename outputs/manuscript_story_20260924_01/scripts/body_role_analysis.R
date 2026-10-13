#!/usr/bin/env Rscript
# Focused extension of the body-size candidate set. Existing analyses are inputs only.
source("scripts/common.R")
suppressPackageStartupMessages({library(phylopath); library(phylolm); library(jsonlite)})
OUT <- file.path(RUN, "figures/figure_13")
ANA <- file.path(OUT, "analysis")
SNAP <- file.path(RUN, "inputs/body_roles")
dir.create(ANA, recursive = TRUE, showWarnings = FALSE)
dir.create(SNAP, recursive = TRUE, showWarnings = FALSE)
sha <- function(p) sub(" .*", "", system2("sha256sum", shQuote(p), stdout = TRUE))
sources <- c(traits = file.path(RUN, "figures/figure_8/analysis/traits_joined.csv"),
  tree = file.path(RUN, "inputs/body_size/path24_time_tree.nwk"))
dest <- file.path(SNAP, c("traits24.csv", "tree24.nwk")); names(dest) <- names(sources)
for (id in names(sources)) {
  if (!file.exists(dest[[id]])) stopifnot(file.copy(sources[[id]], dest[[id]]))
  stopifnot(identical(sha(sources[[id]]), sha(dest[[id]])))
}
write_csv(tibble(source = names(sources), path = unname(sources), snapshot = unname(dest),
  sha256 = vapply(dest, sha, character(1))), file.path(SNAP, "MANIFEST.csv"))
d <- read_csv(dest[["traits"]], show_col_types = FALSE)
tree <- read.tree(dest[["tree"]])
stopifnot(nrow(d) == 24L, !anyDuplicated(d$tree_tip), setequal(d$tree_tip, tree$tip.label))
models <- define_model_set(
  body_upstream_cell = c(g ~ b, c ~ b + g, n ~ b + c),
  body_upstream_nucleus = c(g ~ b, n ~ b + g, c ~ b + n),
  body_upstream_genome = c(g ~ b, c ~ b + g, n ~ b + g),
  body_only = c(g ~ b, c ~ b, n ~ b),
  cell_body_branch = c(c ~ g, n ~ c, b ~ c),
  nucleus_cell_body_chain = c(n ~ g, c ~ n, b ~ c),
  genome_body_branch = c(c ~ g, n ~ c, b ~ g),
  body_disconnected = c(c ~ g, n ~ c))
model_names <- c(body_upstream_cell = "SVL upstream: cell bridge",
  body_upstream_nucleus = "SVL upstream: nucleus bridge",
  body_upstream_genome = "SVL upstream: genome fork", body_only = "SVL-only fork",
  cell_body_branch = "Cell linked to SVL", nucleus_cell_body_chain = "Nucleus–cell–SVL chain",
  genome_body_branch = "Genome linked to SVL", body_disconnected = "SVL independent")
codes <- setNames(paste0("H", seq_along(models)), names(models))
edge_rows <- function(m) {
  ix <- which(m != 0, arr.ind = TRUE)
  tibble(parent = rownames(m)[ix[, 1]], child = colnames(m)[ix[, 2]])
}
dag_key <- function(m) {
  sk <- (m + t(m)) > 0
  skeleton <- apply(which(upper.tri(sk) & sk, arr.ind = TRUE), 1, function(i)
    paste(sort(rownames(m)[i]), collapse = "-"))
  colliders <- character()
  for (v in colnames(m)) {
    pa <- rownames(m)[m[, v] != 0]
    if (length(pa) >= 2) for (pair in combn(pa, 2, simplify = FALSE)) {
      if (!sk[pair[1], pair[2]]) colliders <- c(colliders,
        paste0(paste(sort(pair), collapse = "+"), "->", v))
    }
  }
  paste(paste(sort(skeleton), collapse = ";"), paste(sort(colliders), collapse = ";"), sep = " / ")
}
acyclic <- function(m) {
  left <- colnames(m)
  while (length(left)) {
    roots <- left[colSums(m[left, left, drop = FALSE]) == 0]
    if (!length(roots)) return(FALSE)
    left <- setdiff(left, roots)
  }
  TRUE
}
equivalent <- function(m) {
  sk <- m + t(m); ij <- which(upper.tri(sk) & sk > 0, arr.ind = TRUE)
  orientations <- expand.grid(rep(list(c(0L, 1L)), nrow(ij)))
  found <- list()
  for (k in seq_len(nrow(orientations))) {
    a <- m * 0
    for (e in seq_len(nrow(ij))) {
      u <- ij[e, 1]; v <- ij[e, 2]
      if (orientations[k, e] == 1L) a[u, v] <- 1 else a[v, u] <- 1
    }
    if (acyclic(a) && identical(dag_key(a), dag_key(m))) found[[length(found) + 1L]] <- a
  }
  found
}
keys <- vapply(models, dag_key, character(1))
stopifnot(!anyDuplicated(keys), all(vapply(models, acyclic, logical(1))))
equiv <- lapply(models, equivalent)
write_csv(bind_rows(lapply(names(models), function(id) edge_rows(models[[id]]) %>%
  mutate(model = id, code = codes[[id]], hypothesis = model_names[[id]]))), file.path(ANA, "candidate_edges.csv"))
write_csv(tibble(model = names(models), code = unname(codes), hypothesis = unname(model_names),
  equivalence_key = keys, equivalent_orientations = lengths(equiv)), file.path(ANA, "equivalence_classes.csv"))
write_csv(bind_rows(lapply(names(equiv), function(id) bind_rows(lapply(seq_along(equiv[[id]]), function(k)
  edge_rows(equiv[[id]][[k]]) %>% mutate(model = id, orientation = k))))), file.path(ANA, "equivalent_orientations.csv"))
# Fixed in advance. phylopath uses these to choose the test response only
# for pairs with no directed path between them; ancestral pairs keep their
# graph-implied testing direction. This does not reorient candidate arrows.
orders <- list(primary = c("b", "g", "c", "n"), reverse_test_order = c("n", "c", "g", "b"))
panels <- list(full24 = d, subset22 = d[!d$tree_tip %in% c("aeneus", "wrighti"), ])
stopifnot(nrow(panels$subset22) == 22L)
warn <- list()
capture <- function(expr, label) withCallingHandlers(expr, warning = function(w) {
  warn[[length(warn) + 1L]] <<- tibble(fit = label, warning = conditionMessage(w))
  invokeRestart("muffleWarning")
})
ranks <- claims <- coeffs <- list()
for (panel in names(panels)) {
  dat <- as.data.frame(panels[[panel]][, c("b", "g", "c", "n")]); rownames(dat) <- panels[[panel]]$tree_tip
  phy <- keep.tip(tree, rownames(dat)); dat <- dat[phy$tip.label, , drop = FALSE]
  for (order_id in names(orders)) {
    pp <- capture(phylo_path(models, data = dat, tree = phy, model = "lambda",
      order = orders[[order_id]], lower.bound = 1e-7, upper.bound = 1), paste(panel, order_id))
    if (length(pp$warnings)) warn[[length(warn) + 1L]] <- tibble(fit = paste(panel, order_id),
      warning = unlist(pp$warnings, use.names = FALSE))
    rank <- as_tibble(as.data.frame(summary(pp))) %>% mutate(panel = panel, test_order = order_id,
      n_species = nrow(dat), rank = row_number())
    stopifnot(abs(sum(rank$w) - 1) < 1e-12, all(rank$k > 0), all(rank$q + 1 < nrow(dat)),
      max(abs(rank$delta_CICc - rank$CICc + min(rank$CICc))) < 1e-10)
    ranks[[length(ranks) + 1L]] <- rank
    for (id in names(pp$d_sep)) {
      ds <- pp$d_sep[[id]]
      claims[[length(claims) + 1L]] <- tibble(panel = panel, test_order = order_id, model = id,
        independence_claim = ds$d_sep, p_value = ds$p, lambda = ds$phylo_par)
      r <- rank[rank$model == id, ]
      stopifnot(abs(r$C + 2 * sum(log(ds$p))) < 1e-9,
        abs(r$p - pchisq(r$C, 2 * r$k, lower.tail = FALSE)) < 1e-10)
    }
  }
  for (id in names(models)) for (child in colnames(models[[id]])) {
    pa <- rownames(models[[id]])[models[[id]][, child] != 0]
    if (!length(pa)) next
    f <- reformulate(pa, response = child)
    fit <- capture(phylolm(f, data = dat, phy = phy, model = "lambda", lower.bound = 1e-7,
      upper.bound = 1), paste(panel, id, child))
    sm <- summary(fit)$coefficients; crit <- qt(.975, fit$n - fit$d)
    coeffs[[length(coeffs) + 1L]] <- tibble(panel = panel, model = id, response = child,
      equation = paste(deparse(f), collapse = " "), term = rownames(sm), estimate = sm[, "Estimate"],
      standard_error = sm[, "StdErr"], ci_low = sm[, "Estimate"] - crit * sm[, "StdErr"],
      ci_high = sm[, "Estimate"] + crit * sm[, "StdErr"], p_value = sm[, "p.value"],
      lambda = fit$optpar, n_species = fit$n, residual_df = fit$n - fit$d,
      scaling = "Fixed full24 log10 means and sample SDs")
  }
}
ranking <- bind_rows(ranks); ds <- bind_rows(claims)
write_csv(ranking, file.path(ANA, "model_comparison.csv"))
write_csv(ds, file.path(ANA, "dsep_tests.csv"))
write_csv(bind_rows(coeffs), file.path(ANA, "structural_coefficients.csv"))
write_csv(bind_rows(warn), file.path(ANA, "warnings.csv"))

style <- theme_gbe() + theme(plot.title = element_text(size = 10, face = "bold"),
  plot.subtitle = element_text(size = 8), plot.caption = element_text(size = 7.5, hjust = 0),
  axis.text.y = element_text(size = 8), strip.text = element_text(size = 8),
  plot.margin = margin(6, 8, 6, 8), legend.position = "bottom")
blue <- GBE_COLORS[["blue"]]; gray <- "#777777"
primary <- ranking %>% filter(test_order == "primary")
ordering <- primary %>% filter(panel == "full24") %>% arrange(CICc) %>% pull(model)
ranking_plot_data <- ranking %>% mutate(y = 9 - match(model, ordering),
  y = y + ifelse(test_order == "primary", .12, -.12),
  order_label = ifelse(test_order == "primary", "Primary test order", "Reverse test order"),
  cohort = factor(ifelse(panel == "full24", "All species (n = 24)", "Without aeneus + wrighti (n = 22)"),
    levels = c("All species (n = 24)", "Without aeneus + wrighti (n = 22)")))
maxdelta <- max(ranking$delta_CICc)
pa <- ggplot(ranking_plot_data, aes(delta_CICc, y, colour = order_label, shape = order_label)) +
  geom_vline(xintercept = 2, linetype = "dashed", colour = gray, linewidth = .3) +
  geom_segment(aes(x = 0, xend = delta_CICc, yend = y), linewidth = .3) + geom_point(size = 2) +
  scale_colour_manual(values = c(`Primary test order` = blue, `Reverse test order` = gray), name = NULL) +
  scale_shape_manual(values = c(`Primary test order` = 16, `Reverse test order` = 1), name = NULL) +
  scale_y_continuous(breaks = 8:1, labels = paste(codes[ordering], model_names[ordering]), limits = c(.5, 8.5)) +
  scale_x_continuous(expand = expansion(mult = c(.02, .06))) + facet_wrap(~cohort, nrow = 1) +
  labs(title = "Does support depend on where body size enters?", x = "ΔCICc within species panel and test order", y = NULL,
    subtitle = "Eight hypotheses, one representative of each distinct equivalence class.",
    caption = "Dashed line: ΔCICc = 2. Test order changes response choice for pairs without a directed path.\nCandidate arrows are hypotheses; rankings do not identify the direction of causation.") + style +
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank())

# Consistent four-node geometry, with the genome–nucleus link routed below cell.
nodes <- tibble(node = c("b", "g", "c", "n"), x = c(2, .35, 2, 3.65), y = c(2.25, .55, .55, .55),
  label = c("SVL", "Genome", "Cell", "Nucleus"))
diagram <- function(m, title, subtitle = NULL, coefficients = NULL) {
  e <- edge_rows(m) %>% left_join(nodes %>% select(parent = node, x, y), by = "parent") %>%
    left_join(nodes %>% select(child = node, xend = x, yend = y), by = "child") %>%
    mutate(long = parent %in% c("g", "n") & child %in% c("g", "n"),
      dx = xend - x, dy = yend - y, shortening = pmin(.52 / pmax(abs(dx), 1e-8), .23 / pmax(abs(dy), 1e-8)),
      sx = x + dx * shortening, sy = y + dy * shortening,
      ex = xend - dx * shortening, ey = yend - dy * shortening)
  short <- e %>% filter(!long); long <- e %>% filter(long)
  p <- ggplot() + geom_segment(data = short, aes(x = sx, y = sy, xend = ex, yend = ey),
    linewidth = .5, colour = "#333333", arrow = grid::arrow(length = unit(1.8, "mm"), type = "closed"))
  if (nrow(long)) p <- p +
    geom_segment(data = long, aes(x = x, y = y - .23, xend = x, yend = -.08), linewidth = .5) +
    geom_segment(data = long, aes(x = x, y = -.08, xend = xend, yend = -.08), linewidth = .5) +
    geom_segment(data = long, aes(x = xend, y = -.08, xend = xend, yend = yend - .23), linewidth = .5,
      arrow = grid::arrow(length = unit(1.8, "mm"), type = "closed"))
  p + geom_rect(data = nodes, aes(xmin = x - .52, xmax = x + .52, ymin = y - .23, ymax = y + .23),
    fill = "white", colour = "#333333", linewidth = .4) +
    geom_text(data = nodes, aes(x, y, label = label), family = GBE_FONT, size = 7.5 / .pt) +
    coord_cartesian(xlim = c(-.25, 4.25), ylim = c(-.2, 2.65), expand = FALSE) +
    labs(title = title, subtitle = subtitle) + theme_void(base_family = GBE_FONT, base_size = 8) +
    theme(plot.title = element_text(size = 8, face = "bold", hjust = 0),
      plot.subtitle = element_text(size = 7.5), plot.margin = margin(5, 6, 5, 6))
}
atlas <- wrap_plots(lapply(names(models), function(id) diagram(models[[id]], paste(codes[[id]], model_names[[id]]),
  sprintf("%d equivalent orientations", length(equiv[[id]])))), ncol = 2) +
  plot_annotation(title = "Candidate structures", caption = "SVL = maximum snout–vent length. Each displayed direction is a representative hypothesis.",
    theme = theme(plot.title = element_text(family = GBE_FONT, size = 10, face = "bold"),
      plot.caption = element_text(family = GBE_FONT, size = 7.5, hjust = 0)))
winner <- primary %>% filter(panel == "full24", rank == 1)
we <- equiv[[winner$model]]
# Show the supplied representative, then alternatives that change the direction
# of body-size edges most. These are selected for exposition, not model refitting.
distances <- vapply(we, function(m) sum(abs(m["b", ] - models[[winner$model]]["b", ])) +
  sum(abs(m[, "b"] - models[[winner$model]][, "b"])), numeric(1))
selected <- list(models[[winner$model]])
for (i in order(distances, decreasing = TRUE)) {
  if (!any(vapply(selected, function(m) identical(m, we[[i]]), logical(1)))) selected[[length(selected) + 1L]] <- we[[i]]
  if (length(selected) >= min(3L, length(we))) break
}
stopifnot(all(vapply(selected, function(m) identical(dag_key(m), dag_key(models[[winner$model]])), logical(1))))
equiv_plot <- wrap_plots(lapply(seq_along(selected), function(i) diagram(selected[[i]],
  if (i == 1) "Displayed representative" else paste("Equivalent direction", i - 1))), ncol = length(selected)) +
  plot_annotation(title = sprintf("%s: the best-ranked class permits %d arrow orientations", codes[[winner$model]], length(we)),
    subtitle = sprintf("Full24 weight = %.3f; global-fit P = %.3f. Same skeleton and unshielded colliders.", winner$w, winner$p),
    caption = "Observational independence claims cannot select among these directions. Arrow widths are schematic.",
    theme = theme(plot.title = element_text(family = GBE_FONT, size = 10, face = "bold"),
      plot.subtitle = element_text(family = GBE_FONT, size = 8),
      plot.caption = element_text(family = GBE_FONT, size = 7.5, hjust = 0)))
legend_a <- paste("Eight biologically motivated candidates are compared on the same four traits using Pagel-lambda phylogenetic path analysis.",
  "All log10 traits retain the full24 means and sample SDs. Full24 and the targeted exclusion of aeneus and wrighti (n=22) are fitted separately.",
  "Primary test ordering is b,g,c,n; a prespecified reverse n,c,g,b checks sensitivity to regression direction with estimated lambda for node pairs without a directed path. Ancestral pairs retain the graph-implied test direction under phylopath.",
  "Delta CICc and normalized weights apply only within one cohort, order and eight-candidate set. The earlier four-candidate weights remain valid within that earlier set.",
  "No saturated graph or exhaustive search is ranked. Higher global-fit P indicates no detected contradiction of the tested claims, not proof of a causal model.")
legend_b <- paste("Candidate arrows define representatives of eight distinct Markov equivalence classes.",
  "Equivalence is checked by the exact undirected skeleton and unshielded colliders; all orientations of each small skeleton are enumerated and cyclic graphs rejected.",
  "Counts describe equivalent acyclic orientations, not additional fitted hypotheses. Body-size downstream models and alternative cellular chains broaden the earlier upstream-only assumption.")
legend_c <- paste(sprintf("The full24 primary winner is %s, weight %.6f and global-fit P %.6f.", winner$model, winner$w, winner$p),
  sprintf("Its class has %d acyclic orientations with identical observational conditional independences.", length(we)),
  "The representative and up to two reorientations are drawn to make this causal nonidentifiability visible. These schematic arrows have no effect-size encoding.",
  "Equivalence is a graphical statement; separately fitted reversed phylogenetic regressions may differ numerically because their residual covariance parameters are estimated.")
atlas_panel <- wrap_elements(panel = atlas) + theme(plot.margin = margin(20, 8, 6, 8))
equiv_panel <- wrap_elements(panel = equiv_plot) + theme(plot.margin = margin(20, 8, 6, 8))
save_candidate(13, "Alternative positions of body size in the path analysis",
  (pa / patchwork::free(atlas_panel) / patchwork::free(equiv_panel)) +
    plot_layout(heights = c(1.05, 1.5, .60)) + plot_annotation(tag_levels = "a"),
  list(a = list(plot = pa, title = "Eight body-size hypotheses and test-direction sensitivity", width_mm = 180, height_mm = 125,
    legend = legend_a, alt_text = "Eight candidate classes are ranked for 24 and 22 species under two fixed independence-test orders."),
    b = list(plot = atlas_panel, title = "Candidate path structures", width_mm = 180, height_mm = 190,
    legend = legend_b, alt_text = "Eight four-node path hypotheses include body size upstream, linked through cell or genome, and independent."),
    c = list(plot = equiv_panel, title = "Equivalent causal orientations of the best-ranked class", width_mm = 180, height_mm = 82,
    legend = legend_c, alt_text = "Different arrow orientations have the same skeleton and unshielded colliders, making them observationally equivalent.")),
  width_mm = 180, height_mm = 325, datasets = "body_roles/traits24.csv;body_roles/tree24.nwk",
  legend = paste("(a)", legend_a, "(b)", legend_b, "(c)", legend_c),
  alt_text = "Alternative body-size roles in a finite eight-class phylogenetic path comparison, with an atlas of hypotheses and explicit equivalent orientations.")
write_json(list(status = "PASS", candidates = length(models), unique_equivalence_classes = length(unique(keys)),
  cohorts = sapply(panels, nrow), test_orders = orders, winner = winner$model, equivalent_orientations = length(we),
  weights_sum_to_one = TRUE, fisher_C_reproduced = TRUE, positive_CICc_denominators = TRUE,
  scaling = "Frozen full24 log10 means and SDs", no_old_results_modified = TRUE), file.path(ANA, "verification.json"),
  pretty = TRUE, auto_unbox = TRUE)
write_lines(c("# Alternative body-size roles", "", legend_a, "", legend_b, "", legend_c,
  "", "Intervals in structural_coefficients.csv use estimate +/- qt(.975, n-p) times SE. They condition on fitted lambda.",
  "No measurement or tree bootstrap is added for this expanded candidate set. Figure 9 tests sampling sensitivity for the original four SVL-adjusted hypotheses."), file.path(ANA, "METHODS.md"))
message("Completed Figure 13 and eight-class path comparison.")
