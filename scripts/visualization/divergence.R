#' Divergence Visualization Script
#'
#' Creates boxplots and phylogeny-based bubble plots for TE divergence statistics.
#'
#' @usage Rscript scripts/visualization/divergence.R

# Load required libraries (installed via conda environment)
library(ggplot2)
library(dplyr)
library(stringr)
library(treeio)       # read.tree()
library(ggtree)       # base tree plotting
library(ape)          # node.depth.edgelength
library(fs)           # path manipulation
library(scales)       # percent labels
library(yaml)         # for config loading

# --- Configuration ---
# Find project root and load config
find_project_root <- function() {
  current <- getwd()
  for (i in 1:10) {
    if (file.exists(file.path(current, "paths.yaml"))) return(current)
    current <- dirname(current)
  }
  stop("Could not find project root (paths.yaml)")
}

project_root <- find_project_root()
config <- yaml::read_yaml(file.path(project_root, "paths.yaml"))

# Define paths from config
results_data <- file.path(project_root, config$results$data)
results_figures <- file.path(project_root, config$results$figures)

input_file <- file.path(results_data, "divergence/divergence_summary_statistics_by_species.csv")

# Read the input file
data <- read.csv(input_file)

# Define output directories
base_boxplot_output_dir <- file.path(results_figures, "divergence/boxplots")
dir.create(base_boxplot_output_dir, showWarnings = FALSE, recursive = TRUE)

base_phylo_output_dir <- file.path(results_figures, "divergence_phylogeny")
dir.create(base_phylo_output_dir, showWarnings = FALSE, recursive = TRUE)

# ─────────────────────────────────────────────────────────────
# Load Tree and Get Coordinates (for Phylogeny Plots)
# ─────────────────────────────────────────────────────────────
tree_file <- "results/data/desmo900dated_test_cleaned_phylo.tre"
if (!file.exists(tree_file)) {
  stop("Tree file not found: ", tree_file)
}
tree <- read.tree(tree_file)

# Calculate max tree depth (root-to-tip distance) for axis transformation
max_depth <- max(ape::node.depth.edgelength(tree))
message(paste("Maximum tree depth (root-to-tip):", round(max_depth, 2)))

# Base tree plot for coordinates
p0_tree <- ggtree(tree)
tree_coords <- p0_tree$data %>%
  filter(isTip) %>%
  select(label, x, y)

# ─────────────────────────────────────────────────────────────
# Boxplot Generation Section
# ─────────────────────────────────────────────────────────────

# Define metrics to plot for boxplots
metrics_to_plot_boxplot <- c(
  "percent_divergence_median" = "Median Percent Divergence (%)",
  "score_median" = "Median Score",
  "percent_deletions_median" = "Median Percent Deletions (%)",
  "percent_insertions_median" = "Median Percent Insertions (%)"
)

# Filter data for threshold 0.9
data_0.9 <- data %>% filter(threshold == 0.9)

# Function to sanitize filenames
sanitize_filename <- function(filename) {
  # Replace spaces and slashes with underscores, remove other special characters
  filename <- str_replace_all(filename, "[ /]", "_")
  filename <- str_replace_all(filename, "[^a-zA-Z0-9_.-]", "")
  return(filename)
}

# Get unique group levels (class, order, superfamily)
unique_group_levels <- unique(data_0.9$group_level)

# Loop through each group level
for (level in unique_group_levels) {
  print(paste("Processing level:", level))
  
  # Create subdirectory for the level
  level_output_dir <- file.path(base_boxplot_output_dir, level)
  dir.create(level_output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Filter data for the current level
  level_data <- data_0.9 %>% filter(group_level == level)
  
  # Get unique group names (TE types) for this level
  unique_group_names <- unique(level_data$group_name)
  
  # Loop through each unique group name (TE type)
  for (te_type in unique_group_names) {
    
    # Filter data for the specific TE type
    plot_data <- level_data %>% filter(group_name == te_type)
    
    if (nrow(plot_data) > 0) {
      print(paste("  Generating plots for:", te_type))
      
      # Loop through each metric
      for (metric_col in names(metrics_to_plot_boxplot)) {
        metric_label <- metrics_to_plot_boxplot[[metric_col]]
        metric_name_for_file <- gsub("[_ ]", "", tolower(gsub("[^a-zA-Z0-9_ ]", "", metric_label))) # Sanitize metric name for filename

        # Sanitize TE type name for filename
        safe_te_name <- sanitize_filename(te_type)
        plot_filename <- file.path(level_output_dir, paste0(safe_te_name, "_", metric_name_for_file, "_0.9.png"))

        # Check if plot already exists
        if (!file.exists(plot_filename)) {
          # Create the plot
          p <- ggplot(plot_data, aes(x = Desmognathus_Species, y = .data[[metric_col]])) +
            geom_boxplot() +
            theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
            labs(title = paste(metric_label, "for", te_type, "(Threshold 0.9)"),
                 x = "Species",
                 y = metric_label)
          
          # Save the plot
          ggsave(plot_filename, plot = p, width = 8, height = 6) # Adjusted size slightly
        } else {
          message(paste("    Skipping existing boxplot:", plot_filename))
        }
      }
      
    } else {
      print(paste("  Skipping plots for:", te_type, "(No data for this TE type at threshold 0.9)"))
    }
  }
}

print("Plot generation complete. Plots saved in subdirectories under results/figures/divergence/boxplots/")


# ─────────────────────────────────────────────────────────────
# Phylogeny Bubble Plot Generation Section (Score Median)
# ─────────────────────────────────────────────────────────────

message("\nStarting Phylogeny Bubble Plot Generation...")

# Prepare divergence data for phylogeny plots
phylo_data_prep <- data_0.9 %>%
  select(Desmognathus_Species, group_level, group_name, score_median) %>%
  # Standardize species names to match tree labels (REMOVE D. prefix)
  mutate(species_label = str_replace(Desmognathus_Species, "D\\.", "")) %>%
  filter(!is.na(score_median)) # Remove rows with NA score_median

# # --- DEBUG: Print Tree Labels and Standardized Data Labels ---
# message("\n--- Tree Tip Labels ---")
# print(tree_coords$label)
# message("\n--- Standardized Divergence Data Species Labels (Unique) ---")
# print(unique(phylo_data_prep$species_label))
# message("-----------------------------------------------------------\n")
# # -----------------------------------------------------------

# Check for species name mismatch after standardization
mismatched_species_phylo <- setdiff(tree$tip.label, phylo_data_prep$species_label)
if (length(mismatched_species_phylo) > 0) {
  warning("Species name mismatch between tree and divergence data.",
          " Missing species in data (after standardization): ", paste(mismatched_species_phylo, collapse=", "))
}
mismatched_data_phylo <- setdiff(phylo_data_prep$species_label, tree$tip.label)
if (length(mismatched_data_phylo) > 0) {
  warning("Species name mismatch between divergence data and tree.",
          " Extra species in data (after standardization): ", paste(mismatched_data_phylo, collapse=", "))
}

# Get unique group levels (class, order, superfamily) again for this section
unique_group_levels_phylo <- unique(phylo_data_prep$group_level)

# Loop through each group level
for (level in unique_group_levels_phylo) {
  # # --- DEBUG: Only process superfamily ---
  # if (level != 'superfamily') next
  # # ---------------------------------------
  message(paste("Processing phylogeny plots for level:", level))

  # Create subdirectory for the level in the phylogeny output dir
  level_phylo_output_dir <- file.path(base_phylo_output_dir, level)
  dir.create(level_phylo_output_dir, showWarnings = FALSE, recursive = TRUE)

  # Filter data for the current level
  level_data_phylo <- phylo_data_prep %>%
      filter(group_level == level)

  # Get unique group names (TE types) for this level
  unique_group_names_phylo <- unique(level_data_phylo$group_name)

  # Calculate the overall max score for this level for consistent scaling
  level_max_score <- max(level_data_phylo$score_median, na.rm = TRUE)
  # Calculate the overall min POSITIVE score for this level
  level_min_pos_score <- min(level_data_phylo$score_median[level_data_phylo$score_median > 0], na.rm = TRUE)
  if (!is.finite(level_min_pos_score)) { level_min_pos_score <- 0 }

  # Calculate lower limit based on user formula (clamped)
  if (level_min_pos_score > 0) {
      intermediate_limit <- 250 - level_min_pos_score
      proposed_lower_limit <- round(intermediate_limit / 100) * 100
      actual_lower_limit <- max(0, min(proposed_lower_limit, floor(level_min_pos_score)))
  } else {
      actual_lower_limit <- 0
  }

  message(paste("  Min positive median score for level", level, ":", round(level_min_pos_score, 2)))
  message(paste("  Max median score for level", level, ":", round(level_max_score, 2)))
  message(paste("  Calculated lower scale limit for level", level, ":", actual_lower_limit))

  # Loop through each unique group name (TE type)
  for (te_type in unique_group_names_phylo) {

    # Manually exclude specific TEs
    exclude_list <- c("Chapaev", "Dada", "Other", "Retrotransposons Unknown", "DNAtransposons Unknown")
    if (te_type %in% exclude_list) {
        message(paste("  Manually skipping phylogeny plot for:", te_type))
        next
    }

    # Filter data for the specific TE type
    te_plot_data <- level_data_phylo %>%
        filter(group_name == te_type)

    # Proceed only if we have data after filtering (and passing the check above)
    if (nrow(te_plot_data) > 0) { # This check might be slightly redundant now but safe
        message(paste("  Generating phylogeny plot for:", te_type))

        # Calculate the mean score for this specific TE type
        te_mean_score <- mean(te_plot_data$score_median, na.rm = TRUE)

        # Join TE data with tree coordinates
        plot_data_joined <- tree_coords %>%
            left_join(te_plot_data, by = c("label" = "species_label")) %>%
            filter(!is.na(group_name)) # Keep only tips with data for this TE type

        # # --- DEBUG: Print joined data ---
        # message("\n--- Data for Gypsy Superfamily Plot ---")
        # print(plot_data_joined)
        # message("-------------------------------------\n")
        # # -----------------------------------

        # Create mapping for labels with scores for the current TE type
        label_map_phylo <- te_plot_data %>%
            select(species_label, score_median) %>%
            mutate(
                label2 = ifelse(
                    is.na(score_median) | score_median <= 0,
                    as.character(species_label),
                    paste0(species_label, " (", round(score_median, 2), ")") # Show score
                )
            ) %>%
            select(label = species_label, label2) # Map original label to new label2

        # Join the new labels back to the full tree data
        tree_data_updated_phylo <- p0_tree$data %>%
            left_join(label_map_phylo, by = "label") %>%
            mutate(label2 = coalesce(label2, as.character(label)))

        # Create the plot
        p_phylo <- p0_tree %+% tree_data_updated_phylo + # Update ggtree object data
            geom_tiplab(
                aes(label = label2),
                size = 2,
                align = TRUE,
                linesize = 0.1,
                offset = 0.5
            ) +
            # Plot points only for tips with data for this TE type (score > 0)
            geom_point(data = plot_data_joined,
                       aes(x = x + 0.05, y = y,
                           size = ifelse(is.na(score_median) | score_median <= 0, NA, score_median),
                           color = score_median), # Add color aesthetic mapped to score
                       alpha = 0.7) +
            scale_size_area( # Use scale_size_area for better perception
                name = "Median Score",
                max_size = 10, # Increased max bubble size slightly
                limits = c(actual_lower_limit, max(actual_lower_limit + 1, level_max_score)), # Use calculated limits
                breaks = scales::pretty_breaks(n = 5) # More breaks
                # labels = default (raw score)
            ) +
            # Add diverging color scale centered on the mean for this TE type
            scale_color_gradient2(
                name = "Score vs Mean",
                low = "blue",
                mid = "lightgrey",
                high = "red",
                midpoint = te_mean_score
            ) +
            guides(size = guide_legend(), color = guide_colorbar()) + # Ensure both legends appear
            theme_tree2() +
            scale_x_continuous(
                name = "Time (Millions of Years Ago)",
                labels = function(x) { round(max_depth - x, 1) },
                breaks = pretty_breaks(n=5)
            ) +
            coord_cartesian(clip = "off") +
            theme(
                legend.position = "right",
                legend.key.size = unit(0.5, "cm"),
                plot.margin = margin(5.5, 40, 5.5, 5.5, "pt"),
                axis.title.x = element_text(margin = margin(t = 5)),
                axis.text.x = element_text(),
                axis.line.x = element_line(),
                axis.ticks.x = element_line()
            ) +
            # Add title including the calculated mean score
            ggtitle(paste0("Median Score for ", te_type, " (", level, ", Threshold 0.9)\nMean Score: ", round(te_mean_score, 1)))

        # Sanitize TE type name for filename
        safe_te_name_phylo <- sanitize_filename(te_type)
        plot_filename_phylo <- file.path(level_phylo_output_dir, paste0(safe_te_name_phylo, "_median_score_phylo_0.9.png"))

        # Save the plot
        # # --- DEBUG: Temporarily disable save ---
        ggsave(plot_filename_phylo, plot = p_phylo, width = 10, height = 8, units = "in", dpi = 300)
        # message(paste("Skipping save for debug: ", plot_filename_phylo))
        # # ---------------------------------------

    } else {
      # This else block might not be reachable anymore due to the earlier checks, but keep for safety?
      # Let's remove it for now as the new check is more specific.
      # message(paste("  Skipping phylogeny plot for:", te_type, "(No data or score_median <= 0)"))
    }
  }
}

message("Phylogeny plot generation complete. Plots saved in subdirectories under ", base_phylo_output_dir)

