library(ggplot2)
library(dplyr)
library(tidyr)
library(viridis)
library(tibble) # For rownames_to_column
library(scales) # For scales::percent
library(ggrepel) # For geom_text_repel

# --- Input Data ---
# Read CSVs, assuming first column is species ID (row names)
read_te_data <- function(filepath) {
  read.csv(filepath, row.names = 1, check.names = FALSE) %>% # Keep original complex names
    rownames_to_column(var = "SpeciesID") %>%
    # Pivot longer to get TE categories and percentages in columns
    pivot_longer(
      cols = -SpeciesID,
      names_to = "Category",
      values_to = "Percentage"
    ) %>%
    # Remove rows where Percentage is NA or zero (often from pivoting sparse tables)
    filter(!is.na(Percentage) & Percentage > 1e-9) # Use a very small threshold
}

# Only read superfamily data now
superfamily_df <- read_te_data("results/data/dnaPipeTE_superfamily_breakdown.csv")

# --- Output Directory ---
output_dir <- "results/figures/TE_diversity_superfamily_donuts" # Directory for single superfamily donuts
if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
}

# --- Helper Function: Prepare Data for Single Donut ---
prepare_donut_data <- function(df, species_id) {

  # Filter for the specific species
  species_df <- df %>% filter(SpeciesID == species_id)

  if (nrow(species_df) == 0) {
    warning(paste("No data found for species", species_id))
    return(NULL)
  }

  # Calculate proportions based on the sum of percentages for this species
  total_percentage <- sum(species_df$Percentage)
  if(total_percentage < 1e-9){
      message(paste("Total percentage near zero for", species_id, ". Skipping."))
      return(NULL)
  }

  df_processed <- species_df %>%
    mutate(proportion = Percentage / total_percentage) %>% # Proportion within the species
    # DO NOT filter small proportions here - plot everything
    # Sort by proportion descending for better label placement later
    arrange(desc(proportion)) %>%
    # Calculate label positions (ypos) for geom_text_repel
    mutate(
      ypos = cumsum(proportion) - 0.5 * proportion,
      Category = as.factor(Category) # Ensure category is factor
    ) %>%
    # Select relevant columns
    select(SpeciesID, Category, Percentage, proportion, ypos)

  if (nrow(df_processed) == 0) {
      message(paste("No data rows after processing for species", species_id))
      return(NULL)
  }

  return(df_processed)
}

# --- Plotting Function: Create Single Donut Plot ---
create_donut_plot <- function(plot_data, species_id, output_dir) {
  if (is.null(plot_data) || nrow(plot_data) == 0) {
      message(paste("No data to plot for:", species_id))
      return(NULL)
  }

  # Define the plot using geom_bar and coord_polar
  plot <- ggplot(plot_data, aes(x = 2, y = proportion, fill = Category)) +
    geom_bar(stat = "identity", width = 1, color = "white", linewidth = 0.2) +
    coord_polar(theta = "y", start = 0) +
    # Add labels outside the donut using ggrepel, only for proportions > 1%
    geom_text_repel(
        aes(y = ypos, label = ifelse(proportion > 0.01, scales::percent(proportion, accuracy = 0.1), NA_character_)),
        nudge_x = 0.8, # Nudge outward
        segment.color = 'grey50',
        segment.size = 0.3,
        force = 1.5, # Repulsion force
        min.segment.length = 0.2,
        box.padding = 0.15,
        point.padding = 0.1,
        size = 2.5,
        color = "black",
        show.legend = FALSE
    ) +
    # Adjust xlim to make space for labels
    xlim(0.5, 2 + 0.8 + 0.5) + # Inner hole radius (0.5) + bar width (1.5) + nudge (0.8) + buffer (0.5)
    scale_fill_viridis_d(option = "turbo", na.value = "grey80") +
    theme_void() +
    theme(
      legend.position = "right",
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      legend.key.size = unit(0.4, "cm"),
      legend.title = element_text(size=10, face="bold"),
      legend.text = element_text(size=7)
    ) +
    labs(
      title = paste(species_id, "- TE Superfamily Composition"),
      fill = "Superfamily"
    )

  filename <- file.path(output_dir, paste0(species_id, "_superfamily_donut.png"))
  ggsave(filename, plot = plot, width = 10, height = 8, dpi = 300, bg = "white")
  message(paste("Saved:", filename))
}

# --- Main Execution ---

if (!"SpeciesID" %in% colnames(superfamily_df)) {
    stop("Error: 'SpeciesID' column not found after processing superfamily_df.")
}
species_list <- unique(superfamily_df$SpeciesID)

for (species in species_list) {
  message(paste("--- Processing:", species, "---"))
  donut_data <- prepare_donut_data(superfamily_df, species)
  create_donut_plot(donut_data, species, output_dir)
}

message("Superfamily donut chart generation complete.")
