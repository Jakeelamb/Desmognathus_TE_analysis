#!/usr/bin/env Rscript
#
# DIVERSITY_INDICES_VISUALIZER.R
#
# Description: Creates visualizations for Shannon and Simpson diversity indices
# calculated by the diversity_indices_processor.py script.
#

# Load necessary libraries
suppressPackageStartupMessages({
  library(ggplot2)
  library(viridis)
  library(dplyr)
  library(tidyr)
  library(gridExtra)
  library(optparse)
})

# Parse command line arguments
option_list <- list(
  make_option(c("-d", "--data_dir"), type="character", default="data/results",
              help="Directory containing the diversity CSV files [default=%default]", metavar="character"),
  make_option(c("-o", "--output_dir"), type="character", default="data/results/plots",
              help="Output directory for plots [default=%default]", metavar="character"),
  make_option(c("-s", "--shannon_file"), type="character", default="shannon_diversity.csv",
              help="Shannon diversity filename [default=%default]", metavar="character"),
  make_option(c("-p", "--simpson_file"), type="character", default="simpson_diversity.csv",
              help="Simpson diversity filename [default=%default]", metavar="character"),
  make_option(c("-c", "--combined_file"), type="character", default="combined_diversity.csv",
              help="Combined diversity filename [default=%default]", metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Check required arguments
if (is.null(opt$data_dir)) {
  stop("Data directory (-d) is required")
}

if (is.null(opt$output_dir)) {
  stop("Output directory (-o) is required")
}

# Create output directory if it doesn't exist
dir.create(opt$output_dir, recursive = TRUE, showWarnings = FALSE)

# Define paths
shannon_path <- file.path(opt$data_dir, opt$shannon_file)
simpson_path <- file.path(opt$data_dir, opt$simpson_file)
combined_path <- file.path(opt$data_dir, opt$combined_file)

# Check if files exist
if (!file.exists(shannon_path) && !file.exists(combined_path)) {
  stop(paste("Shannon diversity file not found:", shannon_path))
}

if (!file.exists(simpson_path) && !file.exists(combined_path)) {
  stop(paste("Simpson diversity file not found:", simpson_path))
}

# Read data
if (file.exists(combined_path)) {
  cat("Reading combined diversity file...\n")
  diversity_df <- read.csv(combined_path)
} else {
  cat("Reading individual diversity files...\n")
  shannon_df <- read.csv(shannon_path)
  simpson_df <- read.csv(simpson_path)
  
  # Merge data
  diversity_df <- merge(shannon_df, simpson_df, by="Species")
}

# Rename columns for consistency if needed
if ("Shannon_Diversity_Index" %in% colnames(diversity_df)) {
  names(diversity_df)[names(diversity_df) == "Shannon_Diversity_Index"] <- "Shannon"
}
if ("Simpson_Diversity_Index" %in% colnames(diversity_df)) {
  names(diversity_df)[names(diversity_df) == "Simpson_Diversity_Index"] <- "Simpson"
}

# Sort data by Shannon index for better visualization
diversity_df <- diversity_df[order(diversity_df$Shannon),]

# Create a custom theme for plots
theme_publication <- theme_minimal() +
  theme(
    text = element_text(family = "serif", size = 12),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  )

# 1. Create Shannon Diversity barplot
shannon_plot <- ggplot(diversity_df, aes(x = reorder(Species, Shannon), y = Shannon)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  labs(
    title = "Shannon Diversity Index by Species",
    x = "Species",
    y = "Shannon Diversity Index"
  ) +
  theme_publication

# Save Shannon plot
shannon_plot_path <- file.path(opt$output_dir, "shannon_diversity_plot.png")
ggsave(shannon_plot_path, shannon_plot, width = 12, height = 7, dpi = 300)
cat(paste("Shannon diversity plot saved to:", shannon_plot_path, "\n"))

# 2. Create Simpson Diversity barplot
simpson_plot <- ggplot(diversity_df, aes(x = reorder(Species, Simpson), y = Simpson)) +
  geom_bar(stat = "identity", fill = "tomato") +
  labs(
    title = "Simpson Diversity Index by Species",
    x = "Species",
    y = "Simpson Diversity Index"
  ) +
  theme_publication

# Save Simpson plot
simpson_plot_path <- file.path(opt$output_dir, "simpson_diversity_plot.png")
ggsave(simpson_plot_path, simpson_plot, width = 12, height = 7, dpi = 300)
cat(paste("Simpson diversity plot saved to:", simpson_plot_path, "\n"))

# 3. Create combined plot
combined_plot <- grid.arrange(shannon_plot, simpson_plot, nrow = 2)

# Save combined plot
combined_plot_path <- file.path(opt$output_dir, "diversity_indices_plot.png")
ggsave(combined_plot_path, combined_plot, width = 12, height = 14, dpi = 300)
cat(paste("Combined diversity plot saved to:", combined_plot_path, "\n"))

# 4. Create correlation plot
correlation_plot <- ggplot(diversity_df, aes(x = Shannon, y = Simpson)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE, color = "blue", alpha = 0.2) +
  geom_text(aes(label = Species), hjust = -0.2, vjust = 0.5, size = 3) +
  labs(
    title = "Correlation between Shannon and Simpson Indices",
    x = "Shannon Diversity Index",
    y = "Simpson Diversity Index"
  ) +
  theme_publication

# Save correlation plot
correlation_plot_path <- file.path(opt$output_dir, "diversity_correlation_plot.png")
ggsave(correlation_plot_path, correlation_plot, width = 10, height = 8, dpi = 300)
cat(paste("Correlation plot saved to:", correlation_plot_path, "\n"))

# 5. Create heatmap of diversity indices
# First reshape data for heatmap
diversity_long <- pivot_longer(
  diversity_df,
  cols = c("Shannon", "Simpson"),
  names_to = "Index",
  values_to = "Value"
)

# Calculate Z-scores for proper comparison across indices
diversity_long <- diversity_long %>%
  group_by(Index) %>%
  mutate(Z_Score = scale(Value)[,1]) %>%
  ungroup()

# Create heatmap
heatmap_plot <- ggplot(diversity_long, aes(x = reorder(Species, Z_Score), y = Index, fill = Z_Score)) +
  geom_tile() +
  scale_fill_viridis() +
  labs(
    title = "Diversity Indices Heatmap",
    x = "Species",
    y = "Diversity Index",
    fill = "Z-Score"
  ) +
  theme_publication +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )

# Save heatmap
heatmap_plot_path <- file.path(opt$output_dir, "diversity_heatmap.png")
ggsave(heatmap_plot_path, heatmap_plot, width = 14, height = 6, dpi = 300)
cat(paste("Heatmap saved to:", heatmap_plot_path, "\n"))

# Calculate and print summary statistics
cat("\nSummary Statistics:\n")
cat("Shannon Diversity Index:\n")
print(summary(diversity_df$Shannon))
cat("\nSimpson Diversity Index:\n")
print(summary(diversity_df$Simpson))

# Calculate correlation between indices
correlation <- cor.test(diversity_df$Shannon, diversity_df$Simpson)
cat("\nCorrelation between Shannon and Simpson indices:\n")
cat(paste("Pearson's r:", round(correlation$estimate, 4), "\n"))
cat(paste("p-value:", format.pval(correlation$p.value), "\n"))

cat("\nAll visualizations completed!\n") 