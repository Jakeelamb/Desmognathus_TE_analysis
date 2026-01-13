#!/usr/bin/env Rscript

# Set CRAN mirror
options(repos = c(CRAN = "https://cloud.r-project.org"))

# Load required libraries
library(ggplot2)
library(viridis)
library(gridExtra)
library(dplyr)

# Set working directory - adjust as needed
# setwd("/home/jake/Projects/Desmognathus_TE")

# Define file paths - updated to reflect reorganized directory structure
five_domain_path <- "Data/ectopic_recombination/five_or_more_domains.csv" # Updated path
tesorter_output <- "Projects/04_ectopic_recombination/combined_sequences.fasta.rexdb-metazoa.cls.tsv" # TEsorter output
ltr_dir <- "Projects/04_ectopic_recombination/LTRs" # Directory containing LTR files
output_dir <- "Results/Plots/Ectopic_Recombination"

# Create output directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Print status information
cat("Using updated directory structure:\n")
cat("- Five domains file:", five_domain_path, "\n")
cat("- TEsorter output:", tesorter_output, "\n")
cat("- LTR directory:", ltr_dir, "\n")
cat("- Output directory:", output_dir, "\n\n")

# Check if required files exist
if (!file.exists(five_domain_path)) {
  stop("Error: Five or more domains file not found at: ", five_domain_path)
}

if (!file.exists(tesorter_output)) {
  cat("Warning: TEsorter output file not found at:", tesorter_output, "\n")
  cat("Continuing with existing five_or_more_domains.csv file only\n\n")
}

if (!dir.exists(ltr_dir)) {
  cat("Warning: LTR directory not found at:", ltr_dir, "\n")
  cat("Continuing with existing five_or_more_domains.csv file only\n\n")
}

# Define the publication theme
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

# Load the prepared five or more domains file
five_plus_domains <- read.csv(five_domain_path)
cat("Loaded data from", five_domain_path, "\n")

# Ensure the Ratio column is correctly named and converted to numeric
if ("Ratio.of.LTR.to.Internal" %in% colnames(five_plus_domains)) {
  # Column already has the right name
  five_plus_domains$Ratio.of.LTR.to.Internal <- as.numeric(five_plus_domains$Ratio.of.LTR.to.Internal)
} else if ("Ratio_of_LTR_to_Internal" %in% colnames(five_plus_domains)) {
  # Column name with underscores
  five_plus_domains$Ratio.of.LTR.to.Internal <- as.numeric(five_plus_domains$Ratio_of_LTR_to_Internal)
} else if ("Ratio of LTR to Internal" %in% colnames(five_plus_domains)) {
  # Original column name with spaces
  five_plus_domains$Ratio.of.LTR.to.Internal <- as.numeric(five_plus_domains$`Ratio of LTR to Internal`)
} else {
  # List all columns to help with debugging
  cat("Available columns in the file:\n")
  print(colnames(five_plus_domains))
  stop("Couldn't find the ratio column in the five_or_more_domains.csv file")
}

# Extract species column if exists, or stop
if (!"Species" %in% colnames(five_plus_domains)) {
  cat("Available columns in the file:\n")
  print(colnames(five_plus_domains))
  stop("Species column not found in the five_or_more_domains.csv file")
}

# Filter out rows with NA or 0 ratio
five_plus_domains <- five_plus_domains[!is.na(five_plus_domains$Ratio.of.LTR.to.Internal) & 
                                     five_plus_domains$Ratio.of.LTR.to.Internal > 0, ]

# Print species distribution for verification
cat("Species distribution in the dataset (5+ domains):\n")
print(table(five_plus_domains$Species))

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

# Create violin plot for five plus domains
p1 <- create_violin_plot(five_plus_domains, "LTR Elements with 5+ Domains\nDistribution")

# Run one-way ANOVA on five_plus_domains dataset
anova_result <- aov(Ratio.of.LTR.to.Internal ~ Species, data = five_plus_domains)
anova_summary <- summary(anova_result)
print("One-way ANOVA results:")
print(anova_summary)

# Save ANOVA results to file
dir.create(file.path("Results", "Stats"), recursive = TRUE, showWarnings = FALSE)
sink(file.path("Results", "Stats", "anova_results.txt"))
cat("One-way ANOVA results for LTR to Internal ratio (elements with 5+ domains):\n\n")
print(anova_summary)

# Add Tukey's HSD post-hoc test
tukey_result <- TukeyHSD(anova_result)
cat("\n\nTukey's HSD post-hoc test:\n\n")
print(tukey_result)
sink()

# Save the plot
ggsave(file.path(output_dir, "elements_5plus_domains.png"), p1, width = 10, height = 8, dpi = 300)

cat("Plot saved to", file.path(output_dir, "elements_5plus_domains.png"), "\n")
cat("Statistical results saved to Results/Stats/anova_results.txt\n")

# Note for future reference
cat("\nNOTE: If you need to regenerate the analysis data, follow these steps:\n")
cat("1. Run the filter_sequences_for3k.py script: python Projects/04_ectopic_recombination/filter_sequences_for3k.py\n")
cat("2. Run the TEsorter.sh script: bash Projects/04_ectopic_recombination/TEsorter.sh\n")
cat("3. The output will be in Projects/04_ectopic_recombination/combined_sequences.fasta.rexdb-metazoa.cls.tsv\n")
cat("4. Run this R script again to generate the plots\n") 