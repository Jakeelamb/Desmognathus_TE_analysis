#!/usr/bin/env Rscript
# Quick PERMANOVA analysis

suppressPackageStartupMessages({
  library(ape); library(vegan); library(dplyr); library(readr); library(ggplot2); library(yaml)
})

project_root <- '/home/jake/Projects/Desmognathus_TE'
config <- yaml::read_yaml(file.path(project_root, 'paths.yaml'))
data_dir <- file.path(project_root, config$results$data)
figures_dir <- file.path(project_root, config$results$figures)
phylo_dir <- file.path(project_root, config$input_data$phylogeny)
output_dir <- file.path(figures_dir, 'permanova')
output_data_dir <- file.path(data_dir, 'permanova')
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(output_data_dir, showWarnings = FALSE, recursive = TRUE)

# Load data
tree <- read.tree(file.path(phylo_dir, 'desmo900dated_test.tre'))
order_df <- read_csv(file.path(data_dir, 'dnaPipeTE_order_breakdown.csv'), show_col_types = FALSE)
superfamily_df <- read_csv(file.path(data_dir, 'dnaPipeTE_superfamily_breakdown.csv'), show_col_types = FALSE)

# Prepare data
order_df[[1]] <- gsub('^D\\.', '', order_df[[1]])
order_df <- order_df %>% rename(species = 1) %>% filter(species %in% tree$tip.label)

superfamily_df[[1]] <- gsub('^D\\.', '', superfamily_df[[1]])
superfamily_df <- superfamily_df %>% rename(species = 1) %>% filter(species %in% tree$tip.label)

# Assign clades by cutting tree
pruned_tree <- drop.tip(tree, setdiff(tree$tip.label, order_df$species))
cophen <- cophenetic(pruned_tree)
hc <- hclust(as.dist(cophen))
clusters <- cutree(hc, k = 5)
species_clades <- data.frame(species = names(clusters), clade = paste0("clade_", clusters))

# Merge with data
order_data <- order_df %>% inner_join(species_clades, by = "species")
superfamily_data <- superfamily_df %>% inner_join(species_clades, by = "species")

# Create composition matrices
order_matrix <- order_data %>% select(where(is.numeric)) %>% as.matrix()
rownames(order_matrix) <- order_data$species
order_matrix[is.na(order_matrix)] <- 0

superfamily_matrix <- superfamily_data %>% select(where(is.numeric)) %>% as.matrix()
rownames(superfamily_matrix) <- superfamily_data$species
superfamily_matrix[is.na(superfamily_matrix)] <- 0

# Calculate distances
order_dist <- vegdist(order_matrix, method = "bray")
superfamily_dist <- vegdist(superfamily_matrix, method = "bray")

# Run PERMANOVA
cat("=== Order-level PERMANOVA ===\n")
order_perm <- adonis2(order_dist ~ clade, data = order_data, permutations = 999)
print(order_perm)

cat("\n=== Superfamily-level PERMANOVA ===\n")
superfamily_perm <- adonis2(superfamily_dist ~ clade, data = superfamily_data, permutations = 999)
print(superfamily_perm)

# Save results
results <- data.frame(
  level = c("Order", "Superfamily"),
  R2 = c(order_perm$R2[1], superfamily_perm$R2[1]),
  F_value = c(order_perm$F[1], superfamily_perm$F[1]),
  p_value = c(order_perm$`Pr(>F)`[1], superfamily_perm$`Pr(>F)`[1])
)
results$significant <- results$p_value < 0.05

write_csv(results, file.path(output_data_dir, 'permanova_summary.csv'))
write_csv(species_clades, file.path(output_data_dir, 'species_clade_assignments.csv'))

# PCoA plot
pcoa_result <- cmdscale(order_dist, k = 2, eig = TRUE)
var_explained <- round(100 * pcoa_result$eig[1:2] / sum(pcoa_result$eig[pcoa_result$eig > 0]), 1)

plot_df <- data.frame(
  PCo1 = pcoa_result$points[, 1],
  PCo2 = pcoa_result$points[, 2],
  Species = order_data$species,
  Clade = order_data$clade
)

p_pcoa <- ggplot(plot_df, aes(x = PCo1, y = PCo2, color = Clade, label = Species)) +
  geom_point(size = 3, alpha = 0.8) +
  ggrepel::geom_text_repel(size = 2.5, max.overlaps = 15) +
  stat_ellipse(aes(group = Clade), level = 0.95, linetype = "dashed", alpha = 0.5) +
  scale_color_brewer(palette = "Set1") +
  labs(title = "PCoA of TE Order Composition",
       subtitle = paste0("PERMANOVA: R2 = ", round(results$R2[1], 3), ", p = ", results$p_value[1]),
       x = paste0("PCo1 (", var_explained[1], "%)"),
       y = paste0("PCo2 (", var_explained[2], "%)")) +
  theme_minimal() +
  coord_fixed()

ggsave(file.path(output_dir, "pcoa_order.png"), p_pcoa, width = 12, height = 10, dpi = 300)

# Distance heatmap
png(file.path(output_dir, "distance_heatmap_order.png"), width = 1200, height = 1000, res = 150)
heatmap(as.matrix(order_dist), main = "Bray-Curtis Dissimilarity (TE Order)",
        col = colorRampPalette(c("white", "steelblue", "darkblue"))(100))
dev.off()

cat("\n=== Summary ===\n")
print(results)
cat("\nSaved to:", output_data_dir, "\n")
cat("Figures saved to:", output_dir, "\n")
