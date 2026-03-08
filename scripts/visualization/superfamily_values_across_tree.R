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

# ─────────────────────────────────────────────────────────────
# Libraries & Setup
# ─────────────────────────────────────────────────────────────
suppressPackageStartupMessages({
  library(treeio)
  library(ggtree)
  library(tidyverse)
  library(viridis)
  library(scales)
  library(fs)
})

if (!exists("is.waive", mode = "function")) {
  is.waive <- function(x) inherits(x, "waiver")
}

project_root <- find_project_root(script_dir)
config <- load_project_config(project_root)

results_data_dir <- resolve_config_path(project_root, config$results$data, "results/data")
results_figures_dir <- resolve_config_path(project_root, config$results$figures$root, "results/figures")

tree_candidates <- c(
  file.path(results_data_dir, "desmo900dated_test_cleaned_phylo.tre"),
  file.path(resolve_config_path(project_root, config$input_data$phylogeny %||% config$data$phylogeny, "input_data/phylogeny"), "desmo900dated_test.tre")
)
tree_file <- tree_candidates[file.exists(tree_candidates)][1]
if (is.na(tree_file)) {
  stop("Tree file not found for superfamily-values-across-tree workflow.", call. = FALSE)
}
tree <- read.tree(tree_file)

# Calculate max tree depth (root-to-tip distance) for axis transformation
# Use ape:: namespace to ensure function is found
max_depth <- max(ape::node.depth.edgelength(tree))
message(paste("Maximum tree depth (root-to-tip):", round(max_depth, 2)))

# Base tree plot for coordinates
p0 <- ggtree(tree)
tree_coords <- p0$data %>%
  filter(isTip) %>%
  select(label, x, y)

# ─────────────────────────────────────────────────────────────
# Function to process and plot TE data for a given category
# ─────────────────────────────────────────────────────────────
process_and_plot_te <- function(data_file, output_dir, category_name) {

  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    message("Created directory: ", output_dir)
  }

  # Read TE data
  if (!file.exists(data_file)) {
    warning("Data file not found: ", data_file, ". Skipping this category.")
    return(invisible(NULL))
  }
  te_raw <- read_csv(data_file, show_col_types = FALSE) %>%
            rename(species = 1) %>% # Rename the first column regardless of its name
            mutate(species = str_remove(species, "^D\\."))

  # Check for species name mismatch
  mismatched_species <- setdiff(tree$tip.label, te_raw$species)
  if (length(mismatched_species) > 0) {
    warning("Species name mismatch between tree and ", data_file,
            ". Missing species: ", paste(mismatched_species, collapse=", "))
  }
  mismatched_data <- setdiff(te_raw$species, tree$tip.label)
    if (length(mismatched_data) > 0) {
    warning("Species name mismatch between ", data_file, " and tree.",
            " Extra species in data: ", paste(mismatched_data, collapse=", "))
  }


  # Pivot data to long format
  te_long <- te_raw %>%
    pivot_longer(cols = -species, names_to = "category_level", values_to = "value") %>%
    mutate(value = as.numeric(value)) # Ensure value is numeric

  # Calculate max value for each category level and filter
  level_max_values <- te_long %>%
    group_by(category_level) %>%
    summarise(max_value = max(value, na.rm = TRUE), .groups = 'drop')

  levels_to_plot <- level_max_values %>%
    filter(max_value >= 0.1) %>%
    pull(category_level)

  if (length(levels_to_plot) == 0) {
    message("No categories met the 0.1% threshold in ", data_file, ". Skipping plotting.")
    return(invisible(NULL))
  }

  # Filter the long data to only include levels we will plot
  te_long_filtered <- te_long %>%
    filter(category_level %in% levels_to_plot)

  # Calculate the overall max value for this category group (for consistent scaling)
  overall_max_value <- max(te_long_filtered$value, na.rm = TRUE)
  message(paste("Overall max value for", category_name, "scaling:", round(overall_max_value, 2)))

  # Join TE data with tree coordinates, filtering out tips with no relevant data
  plot_data <- tree_coords %>%
    left_join(te_long_filtered, by = c("label" = "species")) %>%
    filter(!is.na(category_level)) # Keep only rows where join found a match

  # Get unique category levels (already filtered)
  # unique_levels <- unique(plot_data$category_level)

  # Assign to temporary variable as potential workaround for loop error
  loop_vector <- levels_to_plot

  # Loop through the filtered category levels and create a plot
  for (level in loop_vector) {
    message("Plotting: ", category_name, " - ", level)

    # Filter data for the current level
    level_data <- plot_data %>%
      filter(category_level == level)

    # Create mapping for labels with percentages for the current level
    label_map <- level_data %>%
      select(label, value) %>% # Select only label and value for this level
      mutate(
        label2 = ifelse(
          is.na(value) | value <= 0, # Check if value is NA or non-positive
          as.character(label), # Keep original label if no value or value <= 0
          paste0(label, " (", round(value, 1), "%)") # Append percentage otherwise
        )
      ) %>%
      select(label, label2) # Keep only the original label and the new label2

    # Join the new labels back to the full tree data, ensuring all tips are covered
    tree_data_updated <- p0$data %>%
      left_join(label_map, by = "label") %>% # Join by the original tip label
      # If a tip didn't have a value for this level (and thus no label2),
      # use its original label.
      mutate(label2 = coalesce(label2, as.character(label)))

    # Create the plot using the updated tree data
 
    p_level <- p0 %+% tree_data_updated + # Update ggtree object data with %+%
      geom_tiplab(
        aes(label = label2), # Use the new label2 column for tip labels
        size    = 2,
        align   = TRUE,
        linesize= 0.1,
        offset  = 0.5 # Reset offset to a smaller value first
      ) +
      # Plot points only for tips with data for this level (value > 0)
      geom_point(data = level_data %>% filter(value > 0),
                 aes(x = x + 0.05, y = y, size = value),
                 color = "steelblue",
                 alpha = 0.7) +
      scale_size_area(
        name = "Proportion (%)",
        max_size = 8,
        # Use the overall max value for consistent limits across plots in this category
        limits = c(0, max(1, overall_max_value)),
        breaks = pretty_breaks(n=4),
        labels = label_percent(scale = 1)
      ) +
      guides(size = guide_legend()) +
      # Add base theme for tree scale first and flip the tree scale
      theme_tree2() +
      # Define the reversed x-axis scale
      scale_x_continuous(
          name = "Time (Millions of Years Ago)",
          # Transform labels: max_depth - coordinate = time ago
          labels = function(x) { round(max_depth - x, 1) },
          # Specify breaks for the axis
          breaks = pretty_breaks(n=5)
      ) +
      # Allow labels to plot outside panel
      coord_cartesian(clip = "off") +
      # Apply theme overrides
      theme(
            # legend.position = "right", # Temporarily hide legend
            legend.position = "none",
            legend.key.size = unit(0.5, "cm"),
            # Use the plot.margin setting from user edit for label space
            plot.margin = margin(5.5, 40, 5.5, 5.5, "pt"),
            # Ensure the custom x-axis elements are shown
            axis.title.x = element_text(margin = margin(t = 5)), # Add space above axis title
            axis.text.x = element_text(),
            axis.line.x = element_line(),
            axis.ticks.x = element_line()
           ) +
      ggtitle(paste(category_name, ":", level))

    # Sanitize filename
    safe_level_name <- gsub("[^a-zA-Z0-9_.-]", "_", level) # Replace invalid chars with underscore
    file_name <- path(output_dir, paste0(category_name, "_", safe_level_name, "_bubble_plot.png"))

    # Save the plot
    ggsave(filename = file_name, plot = p_level,
           width = 10, height = 8, units = "in", dpi = 300)

    # Optional: print plot to viewer if running interactively
    # print(p_level)
  }
  message("Finished plotting for: ", category_name)
}

# ─────────────────────────────────────────────────────────────
# Define file paths and call the plotting function for each category
# ─────────────────────────────────────────────────────────────

# Define categories, input files, and output directories
categories <- list(
  list(name = "Superfamily",
       file = file.path(results_data_dir, "dnaPipeTE_superfamily_breakdown.csv"),
       dir = file.path(results_figures_dir, "superfamily_phylogeny")),
  list(name = "Order",
       file = file.path(results_data_dir, "dnaPipeTE_order_breakdown.csv"),
       dir = file.path(results_figures_dir, "order_phylogeny")),
  list(name = "Class",
       file = file.path(results_data_dir, "dnaPipeTE_class_breakdown.csv"),
       dir = file.path(results_figures_dir, "class_phylogeny"))
)

# Loop through each category and generate plots
for (cat in categories) {
  process_and_plot_te(data_file = cat$file,
                      output_dir = cat$dir,
                      category_name = cat$name)
}

message("All plotting complete.")
