# Set CRAN mirror
options(repos = c(CRAN = "https://cloud.r-project.org"))

# Load required libraries
library(ggplot2)
library(viridis)
library(gridExtra)
library(dplyr)
library(ggbeeswarm)

setwd("/home/jake/Projects/Desmognathus_RPC_Results")

# Define the publication theme (unchanged)
theme_publication <- theme_bw() +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
    panel.grid.major = element_line(colour = "grey92"),
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14),
    plot.title = element_text(size = 16, hjust = 0.5),
    strip.background = element_blank(),
    strip.text = element_text(size = 14)
  )

# Read the datasets
five_plus_domains <- read.csv('/home/jake/Projects/Ectopic_recombination/results/five_or_more_domains.csv')
all_domains <- read.csv('/home/jake/Projects/Ectopic_recombination/results/Gypsy_coverage_results_with_species.csv')

# Function to apply IQR filtering (unchanged)
apply_iqr_filter <- function(data) {
  data %>%
    group_by(Species) %>%
    filter(
      Ratio.of.LTR.to.Internal >=
        quantile(Ratio.of.LTR.to.Internal, 0.25) - 1.5 * IQR(Ratio.of.LTR.to.Internal) &
        Ratio.of.LTR.to.Internal <=
        quantile(Ratio.of.LTR.to.Internal, 0.75) + 1.5 * IQR(Ratio.of.LTR.to.Internal)
    ) %>%
    ungroup()
}

# Apply IQR filtering to both datasets
five_plus_domains_filtered <- apply_iqr_filter(five_plus_domains)
all_domains_filtered <- apply_iqr_filter(all_domains)

# Function to create violin plot
create_violin_plot <- function(data, title) {
  ggplot(data, aes(x = Species, y = Ratio.of.LTR.to.Internal, fill = Species)) +
    # Add horizontal threshold line
    geom_hline(yintercept = 1, linetype = "dashed", color = "red", linewidth = 0.8) +
    # Add violin plot
    geom_violin(trim = FALSE, alpha = 0.8) +
    # Add mean bar
    stat_summary(fun = mean, geom = "crossbar", width = 0.5, color = "black", linewidth = 0.6) +
    labs(title = title,
         x = "Species",
         y = "Ratio of LTR to Internal") +
    theme_publication +
    scale_fill_viridis_d() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          legend.position = "none")
}

# Function to create beeswarm plot with reference points
create_beeswarm_plot <- function(data, title) {
  ggplot(data, aes(x = Species, y = Ratio.of.LTR.to.Internal, color = Species)) +
    # Add horizontal threshold line
    geom_hline(yintercept = 1, linetype = "dashed", color = "red", linewidth = 0.8) +
    # Add beeswarm plot with more separation between points
    geom_beeswarm(alpha = 0.8, size = 3, cex = 2.5, priority = "density") +
    # Add mean point
    stat_summary(fun = mean, geom = "point", size = 5, color = "black", shape = 23, fill = "white") +
    # Add vertical guides to separate species
    geom_vline(xintercept = seq(1.5, length(unique(data$Species))-0.5, 1), 
               color = "gray90", linetype = "solid", linewidth = 0.5) +
    # Add horizontal grid lines for better reference
    geom_hline(yintercept = 1:12, color = "gray90", linetype = "solid", linewidth = 0.5) +
    labs(title = title,
         x = "Species",
         y = "Ratio of LTR to Internal") +
    theme_publication +
    scale_color_viridis_d() +
    # Improve x-axis labeling
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 10),
          legend.position = "none",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    # Set reasonable y-axis limits
    scale_y_continuous(breaks = seq(0, 12, 1), limits = c(0, 12.5))
}

# Create all four plots using the violin plot function
p1 <- create_violin_plot(five_plus_domains, "5+ Domains\nOriginal Distribution")
p2 <- create_violin_plot(five_plus_domains_filtered, "5+ Domains\nIQR Filtered")
p3 <- create_violin_plot(all_domains, "All 6 Domains\nOriginal Distribution")
p4 <- create_violin_plot(all_domains_filtered, "All 6 Domains\nIQR Filtered")

# Create beeswarm plot for five plus domains (no IQR filtering)
beeswarm_plot <- create_beeswarm_plot(five_plus_domains, "5+ Domains\nBeeswarm Plot")

# Save the beeswarm plot
ggsave("five_plus_domains_beeswarm.png", beeswarm_plot, width = 10, height = 8, dpi = 300)

# Run one-way ANOVA on five_plus_domains dataset
anova_result <- aov(Ratio.of.LTR.to.Internal ~ Species, data = five_plus_domains)
anova_summary <- summary(anova_result)
print("One-way ANOVA results for five_plus_domains dataset:")
print(anova_summary)

# Save ANOVA results to file
sink("five_plus_domains_anova_results.txt")
cat("One-way ANOVA results for five_plus_domains dataset:\n\n")
print(anova_summary)

# Add Tukey's HSD post-hoc test
tukey_result <- TukeyHSD(anova_result)
cat("\n\nTukey's HSD post-hoc test:\n\n")
print(tukey_result)
sink()

# Arrange plots in a 2x2 grid
combined_plot <- grid.arrange(p1, p2, p3, p4, ncol = 2,
                            padding = unit(1, "line"))

# Save individual plots
ggsave("five_plus_domains_original.png", p1, width = 8, height = 6, dpi = 300)
ggsave("five_plus_domains_filtered.png", p2, width = 8, height = 6, dpi = 300)
ggsave("all_domains_original.png", p3, width = 8, height = 6, dpi = 300)
ggsave("all_domains_filtered.png", p4, width = 8, height = 6, dpi = 300)

# Save combined plots
ggsave("combined_violins.png", combined_plot, width = 15, height = 12, dpi = 300)

IQR_combined <- grid.arrange(p2, p4, ncol = 2)
ggsave("IQR_combined_violins.png", IQR_combined, width = 32, height = 12, dpi = 300)

NoIQR_combined <- grid.arrange(p1, p3, ncol = 2)
ggsave("NoIQR_combined_violins.png", NoIQR_combined, width = 32, height = 12, dpi = 300)