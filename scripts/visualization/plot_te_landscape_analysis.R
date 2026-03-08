#!/usr/bin/env Rscript
#
# TE_LANDSCAPE_PLOTS.R
# 
# Description: Creates visualizations for Transposable Element landscape analysis
# using the current frozen dnaPipeTE breakdown tables in results/data.
#
# Input:
#   - CSV files in results/data/
#
# Output: Plot files (PNG) saved to results/figures/te_landscape/
#

script_dir <- tryCatch({
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = FALSE))
  } else {
    "scripts/visualization"
  }
}, error = function(...) {
  "scripts/visualization"
})

source(file.path(dirname(script_dir), "R", "path_config_utils.R"))
prefer_active_conda_r_library()

# Load necessary libraries (assumes these are already installed)
suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
})

project_root <- find_project_root(script_dir)
config <- load_project_config(project_root)
data_dir <- resolve_config_path(project_root, config$results$data, "results/data")
figures_dir <- resolve_config_path(project_root, config$results$figures, "results/figures")
plot_dir <- file.path(figures_dir, "te_landscape")

# Create output directory for plots
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

# Check if frozen TE data directory exists
if (!dir.exists(data_dir)) {
  stop(paste("TE data directory", data_dir, "not found."))
}

# Load the CSV data & define variables
Class <- read.csv(file.path(data_dir, "dnaPipeTE_class_breakdown.csv"))
Order <- read.csv(file.path(data_dir, "dnaPipeTE_order_breakdown.csv"))
Superfamily <- read.csv(file.path(data_dir, "dnaPipeTE_superfamily_breakdown.csv"))

# Standardize the leading species column across the frozen tables.
colnames(Class)[1] <- "Species"
colnames(Order)[1] <- "Species"
colnames(Superfamily)[1] <- "Species"

# Derive legacy helper tables on the fly from the current class breakdown.
dna_vs_retro <- Class %>%
  transmute(
    Species = Species,
    DNA = DNAtransposons.Subclass1 + DNAtransposons.Subclass2 + DNAtransposons.Unknown,
    Retrotransposons = Retrotransposons.Autonomous + Retrotransposons.Non.autonomous + Retrotransposons.Unknown
  )

known_vs_unknown <- Class %>%
  transmute(
    Species = Species,
    Known = rowSums(across(-Species), na.rm = TRUE) - Unknown,
    Unknown = Unknown
  )

# Define key classification lists
Retro_classes <- c("Retrotransposons Autonomous", "Retrotransposons Non-Autonomous", "Retrotransposons Unknown")
DNA_classes <- c("DNAtransposons Subclass1", "DNAtransposons Subclass2", "DNAtransposons Unknown")
Retro_orders <- c("DIRS", "LINE", "LTR", "SINE")
DNA_orders <- c("Helitron", "Maverick", "PLE", "TIR", "YR")
Retro_superfamilies <- c('Gypsy', 'Jockey', 'L1', 'Penelope', 'tRNA')
DNA_superfamilies <- c("Academ", "CACTA", "Chapaev", "Cyrypton", "DIRS", "Dada",
                      "EnSpm", "Ginger", "Helitron", "MULE", "Maverick", "Mutator",
                      "P", "PIF-Harbinger", "PiggyBac", "Tc1-mariner", "Unknown TIR", "hAT")


### GGplot template
theme_publication <- theme_minimal() +
  theme(
    text = element_text(family = "serif", size = 12),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.position = "right"  # Changed to show legend for cluster plots
  )

### DNA vs Retro plot 
# Process the data
dvr_processed <- dna_vs_retro %>%
  mutate(Total = DNA + Retrotransposons) %>%
  mutate(DNA_Percent = DNA/Total * 100,
         Retro_Percent = Retrotransposons/Total * 100) %>%
  tidyr::pivot_longer(
    cols = c(Retro_Percent, DNA_Percent),
    names_to = "Category",
    values_to = "Percentage"
  )

# Create plot
plot <- ggplot(dvr_processed, aes(x = Species, y = Percentage, fill = Category)) +
  geom_bar(stat = "identity", position = "identity", alpha = 0.9) +
  scale_fill_manual(values = c("DNA_Percent" = "steelblue1",
                              "Retro_Percent" = "tomato"),
                    labels = c("DNA Transposons", "Retrotransposons")) +
  guides(fill = guide_legend(reverse = TRUE)) +
  labs(
    title = "DNA vs Retrotransposons Distribution by Species",
    x = "Species",
    y = "Percentage (%)"
  ) +
  theme_publication +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.margin = margin(t = 20, r = 20, b = 50, l = 20, unit = "pt"),
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_line(color = "grey95")
  ) +
  scale_y_continuous(limits = c(0, 100),
                    breaks = seq(0, 100, 25),
                    expand = c(0, 0))

# Save plot
ggsave(file.path(plot_dir, "dna_vs_retro.png"),
       plot,
       width = 12,
       height = 8,
       dpi = 300)

### End region

### Class Breakdown Plot by species
# First, correct the class vectors to match exactly
Retro_classes <- c("Retrotransposons.Autonomous", 
                  "Retrotransposons.Non.autonomous", 
                  "Retrotransposons.Unknown")

DNA_classes <- c("DNAtransposons.Subclass1", 
                "DNAtransposons.Subclass2", 
                "DNAtransposons.Unknown")

# Update color palette to match the exact column names
color_palette <- c(
  "Retrotransposons.Autonomous" = "#CB181D",
  "Retrotransposons.Non.autonomous" = "#FB6A4A",
  "Retrotransposons.Unknown" = "#FCAE91",
  "DNAtransposons.Subclass1" = "#2171B5",
  "DNAtransposons.Subclass2" = "#6BAED6",
  "DNAtransposons.Unknown" = "#BDD7E7"
)
class_processed <- Class %>%
  select(-Other, -Unknown) %>%
  mutate(Total = rowSums(select(., -Species))) %>%
  pivot_longer(
    cols = c(all_of(c(Retro_classes, DNA_classes))),
    names_to = "Category",
    values_to = "Count"
  ) %>%
  mutate(Percentage = Count/Total * 100)

# Create the plot
plot <- ggplot(class_processed, 
       aes(x = Species, y = Percentage, fill = Category)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = color_palette,
                   name = "Category") +
  labs(
    title = "Transposon Class Distribution by Species",
    x = "Species",
    y = "Percentage (%)"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    text = element_text(family = "serif"),
    plot.title = element_text(hjust = 0.5, size = 16),
    legend.position = "right"
  )

# Save the plot
ggsave(file.path(plot_dir, "class_distribution.png"), plot, width = 12, height = 7, dpi = 300)
### End region

### Class Breakdown Plot mean 
# Calculate percentages and means
# Define class vectors
Retro_classes <- c("Retrotransposons.Autonomous",
                  "Retrotransposons.Non.autonomous",
                  "Retrotransposons.Unknown")
DNA_classes <- c("DNAtransposons.Subclass1",
                "DNAtransposons.Subclass2",
                "DNAtransposons.Unknown")

# Create color palette
color_palette <- c(
  "Retrotransposons.Autonomous" = "#2171B5",
  "Retrotransposons.Non.autonomous" = "#6BAED6",
  "Retrotransposons.Unknown" = "#BDD7E7",
  "DNAtransposons.Subclass1" = "#CB181D",
  "DNAtransposons.Subclass2" = "#FB6A4A",
  "DNAtransposons.Unknown" = "#FCAE91"
)

# Calculate means and percentages
class_means <- Class %>%
  select(-Species, -Other, -Unknown) %>%
  summarise(across(everything(), mean)) %>%
  pivot_longer(
    cols = everything(),
    names_to = "Category",
    values_to = "Mean"
  ) %>%
  # Calculate relative percentages
  mutate(Percentage = (Mean / sum(Mean)) * 100) %>%
  # Sort by percentage (highest to lowest)
  arrange(desc(Percentage))

# Create the plot
class_average <- ggplot(class_means, aes(x = reorder(Category, Percentage), y = Percentage, fill = Category)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = color_palette) +
  labs(
    title = "Relative Abundance of Transposon Classes",
    x = "Transposon Class",
    y = "Relative Abundance (%)"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    text = element_text(family = "serif"),
    plot.title = element_text(hjust = 0.5, size = 16),
    legend.position = "none"
  )

# Save the plot
ggsave(file.path(plot_dir, "class_average.png"), class_average, width = 12, height = 7, dpi = 300)
### End region

### Order Breakdown Plot mean
# Calculate percentages across species
order_means <- Order %>%
  select(-Species) %>%  # Remove Species column
  mutate(RowTotal = rowSums(.)) %>%  # Calculate total in a separate column
  mutate(across(-RowTotal, ~ . / RowTotal * 100)) %>%  # Calculate percentages
  select(-RowTotal) %>%  # Remove the total column
  summarise(across(everything(), mean)) %>%  # Calculate means
  pivot_longer(
    cols = everything(),
    names_to = "Category",
    values_to = "Percentage"
  ) %>%
  # Only keep the orders we want (exclude Total)
  filter(Category != "Total") %>%  # Explicitly filter out Total
  mutate(Type = case_when(
    Category %in% Retro_orders ~ "Retrotransposon",
    Category %in% DNA_orders ~ "DNA transposon"
  )) %>%
  arrange(desc(Percentage))%>%
 filter(Percentage >= 0.05) %>%  # Filter out entries less than 1%
 arrange(desc(Percentage))

# Create the plot
OM_plot <- ggplot(order_means, aes(x = reorder(Category, -Percentage), y = Percentage, fill = Type)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c(
    "Retrotransposon" = "tomato",
    "DNA transposon" = "steelblue1"
  )) +
  scale_y_continuous(
    breaks = seq(0, max(order_means$Percentage), by = 10),  # Adjust the 'by' value to change label frequency
    expand = c(0, 0)  # Remove spacing between axis and bars
  ) +
  labs(
    title = "Mean Transposon Order Distribution Across Species",
    x = "Order",
    y = "Percentage (%)"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    text = element_text(family = "serif"),
    plot.title = element_text(hjust = 0.5, size = 16)
  )

ggsave(file.path(plot_dir, "Order_average.png"), OM_plot, width = 12, height = 7, dpi = 300)


### End region

### Superfamily Breakdown Plot together

Retro_superfamilies <- c('Gypsy', 'Jockey', 'L1', 'Penelope', 'tRNA', 'DIRS')
DNA_superfamilies <- c("Academ", "CACTA", "Chapaev", "Cyrypton", "Dada",
                      "EnSpm", "Ginger", "Helitron", "MULE", "Maverick", "Mutator",
                      "P", "PIF.Harbinger", "PiggyBac", "Tc1.mariner", "Unknown.TIR", "hAT")

# Calculate percentages across species
superfamily_means <- Superfamily %>%
 select(-Species) %>%  # Remove Species column
 mutate(RowTotal = rowSums(.)) %>%  # Calculate total in a separate column
 mutate(across(-RowTotal, ~ . / RowTotal * 100)) %>%  # Calculate percentages
 select(-RowTotal) %>%  # Remove the total column
 summarise(across(everything(), mean)) %>%  # Calculate means
 pivot_longer(
   cols = everything(),
   names_to = "Category",
   values_to = "Percentage"
 ) %>%
 mutate(Type = case_when(
   Category %in% Retro_superfamilies ~ "Retrotransposon",
   Category %in% DNA_superfamilies ~ "DNA transposon"
 )) %>%
 filter(Percentage >= 0.5) %>%  # Filter out entries less than 1%
 arrange(desc(Percentage))

# Create the plot
sm_plot <- ggplot(superfamily_means, aes(x = reorder(Category, -Percentage), y = Percentage, fill = Type)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c(
    "Retrotransposon" = "tomato",
    "DNA transposon" = "steelblue1"
  )) +
  scale_y_continuous(
    breaks = seq(0, 65, by = 5),  # Creates breaks every 5% from 0 to 65
    minor_breaks = seq(0, 65, by = 1)  # Adds minor gridlines every 1%
  ) +
  labs(
    title = "Mean Transposon Superfamily Distribution Across Species",
    x = "Superfamily",
    y = "Percentage (%)"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    text = element_text(family = "serif"),
    plot.title = element_text(hjust = 0.5, size = 16),
    panel.grid.minor = element_line(color = "gray90")  # Makes minor gridlines lighter
  )

ggsave(file.path(plot_dir, "Superfamily_average.png"), sm_plot, width = 12, height = 7, dpi = 300)
### End region

### Classified vs Unclassified Masked Bases Chart
# Process the data
# Process the data (same as before)
ku_processed <- known_vs_unknown %>%
  mutate(Total = Known + Unknown) %>%
  mutate(Known_Percent = Known/Total * 100,
         Unknown_Percent = Unknown/Total * 100) %>%
  tidyr::pivot_longer(
    cols = c(Known_Percent, Unknown_Percent),
    names_to = "Category",
    values_to = "Percentage"
  )

# Create plot
plot <- ggplot(ku_processed, aes(x = Species, y = Percentage, fill = Category)) +
  geom_bar(stat = "identity", position = "identity") +  # Changed from position = "stack"
  scale_fill_manual(values = c("Known_Percent" = "#0072B2", 
                              "Unknown_Percent" = "#E69F00"),
                   labels = c("Classified", "Unclassified")) +
  labs(
    title = "Classified vs Unclassified Percentage Distribution by Species",
    x = "Species",
    y = "Percentage (%)"
  ) +
  theme_publication +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.margin = margin(t = 20, r = 20, b = 50, l = 20, unit = "pt"),
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_line(color = "grey95")
  ) +
  scale_y_continuous(limits = c(0, 100), 
                    breaks = seq(0, 100, 25),
                    expand = c(0, 0))

# Save plot
ggsave(file.path(plot_dir, "known_vs_unknown.png"), 
       plot, 
       width = 12, 
       height = 8, 
       dpi = 300)
### End region

# Print summary of actions
cat("\nTransposable Element Landscape Visualizations\n")
cat("------------------------------------------\n")
cat("Output plots saved to:", plot_dir, "\n")
cat("Number of plots generated:", 6, "\n\n")

# List generated plots
plot_files <- list.files(plot_dir, pattern = "*.png")
cat("Generated plot files:\n")
for (f in plot_files) {
  cat("  -", f, "\n")
}
