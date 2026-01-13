#!/usr/bin/env Rscript

# Load required packages
library(ape)
library(tidyverse)
library(ggtree)
library(phytools)
library(viridis)
library(readr)
library(dplyr)
library(tidyr)

# Function to clean species names
clean_species_name <- function(name) {
  cleaned <- gsub("_.*$", "", name)
  cleaned <- gsub(" .*$", "", cleaned)
  cleaned <- gsub("^D\\.", "", cleaned)
  cleaned <- tolower(cleaned)
  return(cleaned)
}

# Load Phylogenetic Tree
tree_path <- "data/raw/Phylogeny/desmo900dated_test.tre"
tree <- read.tree(tree_path)
original_names <- tree$tip.label
tree$tip.label <- sapply(tree$tip.label, clean_species_name)

# Load Lookup Table
lookup_path <- "data/raw/lookup/lookup_table.txt"
lookup_table <- read_tsv(lookup_path, show_col_types = FALSE) %>%
  rename(SRX_ID = SRA_Accension) %>% # Rename the SRA accession column
  mutate(SRX_ID = tolower(SRX_ID), # Convert SRX_ID to lowercase
         Species = sapply(Species, clean_species_name))

# Load Diversity Data (Corrected Path and Column Selection)
diversity_path <- "results/tables/diversity/diversity_metrics.csv"
shannon_data <- read_csv(diversity_path, show_col_types = FALSE) %>%
  rename(SRX_ID = species, Shannon = shannon_entropy) %>%
  mutate(SRX_ID = tolower(SRX_ID)) %>% # Convert SRX_ID to lowercase
  select(SRX_ID, Shannon)

# Load and Process Kimura Landscape Data
kimura_path <- "data/processed/te_landscape/combined_parsed_align.tsv"
kimura_data_raw <- read_tsv(kimura_path, show_col_types = FALSE)

# Calculate percentage within each sample and find peaks
kimura_peaks <- kimura_data_raw %>%
  mutate(sample_id = tolower(sample_id)) %>% # Convert sample_id to lowercase
  group_by(sample_id) %>% # Group by sample_id first
  mutate(total_sample_bp = sum(aligned_bp, na.rm = TRUE)) %>% # Calculate total bp per sample
  filter(total_sample_bp > 0) %>% # Avoid division by zero
  mutate(perc_rep = (aligned_bp / total_sample_bp) * 100) %>% # Calculate percentage
  filter(perc_rep == max(perc_rep, na.rm = TRUE)) %>% # Find max percentage
  # In case of ties for max percentage, take the lowest Kimura value
  # Also handle potential non-numeric Kimura values if Kimura_bin was used
  mutate(Kimura_numeric = suppressWarnings(as.numeric(gsub("–.*%?", "", Kimura)))) %>% 
  filter(!is.na(Kimura_numeric)) %>%
  filter(Kimura_numeric == min(Kimura_numeric, na.rm = TRUE)) %>% 
  slice(1) %>% # Take the first row in case of remaining ties
  ungroup() %>% 
  rename(SRX_ID = sample_id) %>% # Rename sample_id to SRX_ID
  select(SRX_ID, peak_kimura = Kimura_numeric) # Select final columns

# Function to collapse duplicate tips
collapse_duplicate_tips <- function(tree, original_names) {
  name_mapping <- data.frame(
    original = original_names,
    cleaned = tree$tip.label,
    stringsAsFactors = FALSE
  )
  tip_counts <- table(tree$tip.label)
  duplicate_species <- names(tip_counts[tip_counts > 1])
  new_tree <- tree
  for(species in duplicate_species) {
    duplicate_tips <- which(tree$tip.label == species)
    if(length(duplicate_tips) > 1) {
      tips_to_remove <- duplicate_tips[-1]
      new_tree <- drop.tip(new_tree, tips_to_remove)
    }
  }
  return(new_tree)
}

# Collapse duplicates
tree_collapsed <- collapse_duplicate_tips(tree, original_names)

# Get valid species from lookup table for filtering
valid_species <- unique(lookup_table$Species)

# Filter tree
tree_filtered <- drop.tip(tree_collapsed, tree_collapsed$tip.label[!tree_collapsed$tip.label %in% valid_species])

# Save the processed tree
output_path <- "data/processed/desmognathus_processed.tre"
write.tree(tree_filtered, file = output_path)

# Print summary
cat("Original tree had", length(tree$tip.label), "tips\n")
cat("Processed tree has", length(tree_filtered$tip.label), "tips\n")
cat("Removed", length(tree$tip.label) - length(tree_filtered$tip.label), "tips\n")
cat("Processed tree saved to:", output_path, "\n")

# Merge metrics with lookup table
metrics_merged <- lookup_table %>%
  left_join(shannon_data, by = "SRX_ID") %>%
  left_join(kimura_peaks, by = "SRX_ID") %>%
  filter(!is.na(Shannon) & !is.na(peak_kimura)) # Ensure we have data for both metrics

# Aggregate metrics per species (tip label)
annotation_data <- metrics_merged %>%
  group_by(Species) %>%
  summarise(
    shannon_diversity = mean(Shannon, na.rm = TRUE),
    peak_kimura = mean(peak_kimura, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  # Ensure Species column matches tip labels
  filter(Species %in% tree_filtered$tip.label) %>%
  # Set Species as row names for ggtree
  column_to_rownames("Species")

# --- Debugging: Compare species names ---
cat("\nSpecies names in annotation data before filtering:\n")
print(sort(unique(metrics_merged$Species)))
cat("\nSpecies names in filtered tree tip labels:\n")
print(sort(unique(tree_filtered$tip.label)))
cat("----------------------------------------\n")
# ----------------------------------------

# Check if annotation data is empty
if(nrow(annotation_data) == 0) {
  stop("No matching annotation data found for the species in the filtered tree.")
}

# Function to identify monophyletic clades
identify_monophyletic_clades <- function(tree, min_size = 3) {
  internal_nodes <- unique(tree$edge[,1])
  clades <- list()
  clade_count <- 1
  for(node in internal_nodes) {
    desc <- getDescendants(tree, node)
    tips <- desc[desc <= length(tree$tip.label)]
    if(length(tips) >= min_size) {
      clades[[clade_count]] <- list(node = node, tips = tips, size = length(tips))
      clade_count <- clade_count + 1
    }
  }
  clades <- clades[order(sapply(clades, function(x) -x$size))]
  clade_data <- data.frame(
    node = sapply(clades, function(x) x$node),
    clade = paste("Clade", seq_along(clades))
  )
  return(clade_data)
}

# Identify clades
clade_data <- identify_monophyletic_clades(tree_filtered, min_size = 3)

# Setup colors
n_clades <- nrow(clade_data)
distinct_colors <- colorRampPalette(c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF"))(n_clades)
names(distinct_colors) <- clade_data$clade
plot_colors <- c("0" = "black", distinct_colors) # Use "0" for default group

# Group clades
clade_nodes <- setNames(clade_data$node, clade_data$clade)
grouped_tree <- groupClade(tree_filtered, clade_nodes)

# Create the base plot with colored branches
p_base <- ggtree(grouped_tree, aes(color = group)) +
  geom_treescale(x = 0, y = -1, width = 10, fontsize = 3, offset = 1) +
  geom_tiplab(size = 3, offset = 2, align = TRUE) +
  geom_cladelab(data = clade_data,
                mapping = aes(node = node, label = clade),
                offset = 6, barsize = 0, angle = 0,
                fontsize = 3, offset.text = 0.5) +
  scale_color_manual(values = plot_colors, name = "Clade") +
  theme_tree2() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5)) +
  ggtitle("Phylogeny with TE Metrics") +
  xlim(NA, 40) # Extend x-axis further for heatmaps

# Add Shannon Diversity Heatmap
p_shannon <- gheatmap(p_base, annotation_data[,"shannon_diversity", drop=FALSE],
                      offset = 10, width = 0.05,
                      colnames_angle = 90, colnames_offset_y = -2.5,
                      font.size = 3)

# Add Peak Kimura Heatmap
p_kimura <- gheatmap(p_shannon, annotation_data[,"peak_kimura", drop=FALSE],
                     offset = 14, width = 0.05,
                     colnames_angle = 90, colnames_offset_y = -2.5,
                     font.size = 3)

# Apply color scales to heatmaps
p_final <- p_kimura +
  scale_fill_viridis_c(option = "plasma", name = "Shannon Diversity", guide = guide_colorbar(order=1)) +
  scale_fill_viridis_c(option = "magma", name = "Peak Kimura Div.", guide = guide_colorbar(order=2))

# Save the final plot
ggsave("results/figures/phylogeny_with_metrics.png", p_final, width = 18, height = 10, dpi = 300)

# Print the time scale information
cat("\nTime scale information:\n")
print(range(branching.times(tree_filtered)))

# Print clade information
cat("\nIdentified clades:\n")
print(clade_data) 