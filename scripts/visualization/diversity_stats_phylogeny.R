library(ggplot2)
library(dplyr)
library(stringr)
library(treeio)       # read.tree()
library(ggtree)       # base tree plotting
library(ape)          # node.depth.edgelength
library(fs)           # path manipulation
library(scales)       # pretty_breaks, etc.
library(viridis)

# --- Configuration ---
BASE_DIR <- "/home/jake/Projects/Desmognathus_TE"
INPUT_DATA_DIR <- file.path(BASE_DIR, "results/data")
FIGURE_DIR <- file.path(BASE_DIR, "results/figures/diversity_phylogeny")
TREE_FILE <- file.path(INPUT_DATA_DIR, "desmo900dated_test_cleaned_phylo.tre")
ORDER_DIVERSITY_FILE <- file.path(INPUT_DATA_DIR, "long_format_diversity_order_stats.csv")
SUPERFAMILY_DIVERSITY_FILE <- file.path(INPUT_DATA_DIR, "long_format_diversity_superfamily_stats.csv")
TARGET_THRESHOLD <- 0.01

# Create output directory
dir.create(FIGURE_DIR, showWarnings = FALSE, recursive = TRUE)

# --- Check Files ---
if (!file.exists(TREE_FILE)) stop("Tree file not found: ", TREE_FILE)
if (!file.exists(ORDER_DIVERSITY_FILE)) stop("Order diversity file not found: ", ORDER_DIVERSITY_FILE)
if (!file.exists(SUPERFAMILY_DIVERSITY_FILE)) stop("Superfamily diversity file not found: ", SUPERFAMILY_DIVERSITY_FILE)

# --- Load Tree and Data ---
message("Loading tree...")
tree <- read.tree(TREE_FILE)

message("Loading diversity data...")
order_div <- read.csv(ORDER_DIVERSITY_FILE)
superfamily_div <- read.csv(SUPERFAMILY_DIVERSITY_FILE)

# Combine diversity data (optional, could process separately)
# For simplicity, we'll process separately as categories are distinct

# --- Prepare Tree Data ---
message("Preparing tree coordinates...")
max_depth <- max(ape::node.depth.edgelength(tree))
message(paste("Maximum tree depth (root-to-tip):", round(max_depth, 2)))
p0_tree <- ggtree(tree)
tree_coords <- p0_tree$data %>% 
  filter(isTip) %>% 
  select(label, x, y)

# --- Define Plotting Function ---
plot_diversity_phylo <- function(diversity_data, category_title, target_threshold, tree_coords, p0_tree, max_depth, output_dir) {
  message(paste("\n--- Processing", category_title, "Diversity ---"))
  
  # Filter for target threshold and prepare data
  div_filtered <- diversity_data %>% 
    filter(Threshold == target_threshold) %>% 
    mutate(species_label = str_replace(Species, "D\\.", "")) %>% # Standardize species names
    select(species_label, Simpson, Shannon, Pielou)
  
  # Check for species mismatches
  mismatched_tree <- setdiff(tree_coords$label, div_filtered$species_label)
  if (length(mismatched_tree) > 0) {
    warning("Species in tree but not in ", category_title, " diversity data (Threshold ", target_threshold, "): ", paste(mismatched_tree, collapse=", "))
  }
  mismatched_data <- setdiff(div_filtered$species_label, tree_coords$label)
  if (length(mismatched_data) > 0) {
    warning("Species in ", category_title, " diversity data (Threshold ", target_threshold, ") but not in tree: ", paste(mismatched_data, collapse=", "))
  }
  
  # Join diversity data with tree coordinates
  plot_data_base <- tree_coords %>% 
    left_join(div_filtered, by = c("label" = "species_label"))

  metrics_to_plot <- c("Simpson", "Shannon", "Pielou")
  
  for (metric in metrics_to_plot) {
    message(paste("  Generating plot for", metric, "index..."))
    
    metric_sym <- sym(metric) # Convert string to symbol for aes
    
    # Keep only tips with valid (non-NA) data for this metric
    plot_data_metric <- plot_data_base %>% 
        filter(!is.na(!!metric_sym))

    if (nrow(plot_data_metric) == 0) {
        message(paste("    Skipping plot for", metric, "- no valid data found for threshold", target_threshold))
        next
    }

    # Set min/max for scaling
    metric_values <- plot_data_metric %>% pull(!!metric_sym)

    if (metric == "Simpson" || metric == "Pielou") {
        min_val <- 0
        max_val <- 1
        message(paste("    Setting fixed scale [0, 1] for", metric))
    } else if (metric == "Shannon") {
        min_val <- 1
        max_val <- 1.75
        message(paste("    Setting fixed scale [1, 1.75] for", metric))
    } else { # For any other metrics
        min_val <- min(metric_values, na.rm = TRUE)
        max_val <- max(metric_values, na.rm = TRUE)
        # Ensure min_val is slightly less than max_val for scale limits if they are equal
        if (min_val >= max_val) { 
            min_val <- min_val * 0.95 # Adjust slightly
            max_val <- max_val * 1.05
        }
        # Handle case where all values are 0 or very close
        if (max_val == 0 || max_val - min_val < 1e-6) {
            max_val <- max(0.1, max_val * 1.1) # Ensure some range
        }
        message(paste("    Using dynamic data range for", metric, ":", round(min_val, 3), "-", round(max_val, 3)))
    }

    message(paste("    Data range for", metric, ":", round(min_val, 3), "-", round(max_val, 3)))

    # Create labels with metric values
    label_map <- plot_data_metric %>% 
      select(label, !!metric_sym) %>% 
      mutate(
        label2 = paste0(label, " (", round(!!metric_sym, 3), ")")
      ) %>% 
      select(label, label2)

    # Update tree data with new labels
    tree_data_updated <- p0_tree$data %>% 
      left_join(label_map, by = "label") %>% 
      mutate(label2 = coalesce(label2, as.character(label)))

    # Create plot
    p_phylo <- p0_tree %+% tree_data_updated + # Use updated data
        geom_tiplab(
            aes(label = label2),
            size = 2,
            align = TRUE,
            linesize = 0.1,
            offset = 0.5 
        ) +
        # Plot points for tips with data
        geom_point(data = plot_data_metric, 
                   aes(x = x + 0.05, y = y, 
                       size = !!metric_sym, 
                       color = !!metric_sym), 
                   alpha = 0.8) +
        scale_size_continuous( 
            name = paste(metric, "Index"),
            range = c(1, 8), # Adjust bubble size range
            limits = c(min_val, max_val),
            breaks = pretty_breaks(n = 4)
        ) + 
        scale_color_viridis( # Use viridis color scale
            name = paste(metric, "Index"), 
            option = "plasma", # Or "viridis", "magma", "inferno"
            limits = c(min_val, max_val),
            breaks = pretty_breaks(n = 4)
        ) +
        guides(size = guide_legend(title.position="top", title.hjust = 0.5),
               color = guide_colorbar(title.position="top", title.hjust = 0.5)) +
        theme_tree2() +
        scale_x_continuous(
            name = "Time (Millions of Years Ago)",
            labels = function(x) { round(max_depth - x, 1) },
            breaks = pretty_breaks(n=5)
        ) +
        coord_cartesian(clip = "off") + 
        theme(
            legend.position = "right",
            legend.box = "vertical",
            legend.key.size = unit(0.5, "cm"),
            legend.title = element_text(size=9),
            legend.text = element_text(size=8),
            plot.margin = margin(5.5, 40, 5.5, 5.5, "pt"), # Right margin for labels
            axis.title.x = element_text(margin = margin(t = 5)),
            axis.text.x = element_text(),
            axis.line.x = element_line(),
            axis.ticks.x = element_line()
        ) +
        ggtitle(paste0(category_title, "-level ", metric, " Diversity (Threshold: ", target_threshold * 100, "%)"))

    # Save plot
    plot_filename <- file.path(output_dir, paste0(tolower(category_title), "_", tolower(metric), "_phylo_thresh", target_threshold, ".png"))
    ggsave(plot_filename, plot = p_phylo, width = 10, height = 8, units = "in", dpi = 300)
    message(paste("    Plot saved to:", plot_filename))
  }
}

# --- Generate Plots ---
plot_diversity_phylo(order_div, "Order", TARGET_THRESHOLD, tree_coords, p0_tree, max_depth, FIGURE_DIR)
plot_diversity_phylo(superfamily_div, "Superfamily", TARGET_THRESHOLD, tree_coords, p0_tree, max_depth, FIGURE_DIR)

message("\n--- Script finished ---")
