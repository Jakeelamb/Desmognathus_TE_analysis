#!/usr/bin/env Rscript

library(ggplot2)
library(dplyr)
library(tidyr)
library(viridis)
library(tibble) # For rownames_to_column
library(colorspace) # For lightening colors
# Note: Labels are omitted in this hierarchical version due to complexity

# --- Classification Hierarchy ---
# Manually derived from Python CLASSIFICATION_MAP for relevant Superfamilies
hierarchy_map <- tibble::tribble(
  ~Superfamily, ~Order, ~Class,
  "Academ", "TIR", "DNAtransposons Subclass1",
  "CACTA", "TIR", "DNAtransposons Subclass1",
  "Chapaev", "TIR", "DNAtransposons Subclass1",
  "Cyrypton", "YR", "DNAtransposons Subclass1",
  "Dada", "TIR", "DNAtransposons Subclass1",
  "EnSpm", "TIR", "DNAtransposons Subclass1",
  "Ginger", "TIR", "DNAtransposons Subclass1",
  "PIF-Harbinger", "TIR", "DNAtransposons Subclass1",
  "hAT", "TIR", "DNAtransposons Subclass1",
  "P", "TIR", "DNAtransposons Subclass1",
  "PiggyBac", "TIR", "DNAtransposons Subclass1",
  "Tc1-mariner", "TIR", "DNAtransposons Subclass1",
  "MULE", "TIR", "DNAtransposons Subclass1",
  "Mutator", "TIR", "DNAtransposons Subclass1",
  # "Unknown TIR", "TIR", "DNAtransposons Subclass1", # Excluded as requested
  "Helitron", "Helitron", "DNAtransposons Subclass2",
  "Maverick", "Maverick", "DNAtransposons Subclass2",
  "L1", "LINE", "Retrotransposons Autonomous",
  "Jockey", "LINE", "Retrotransposons Autonomous",
  "LINE", "LINE", "Retrotransposons Autonomous",
  "Copia", "LTR", "Retrotransposons Autonomous",
  "Gypsy", "LTR", "Retrotransposons Autonomous",
  "LTR", "LTR", "Retrotransposons Autonomous",
  "DIRS", "DIRS", "Retrotransposons Autonomous",
  "Penelope", "PLE", "Retrotransposons Autonomous",
  "SINE", "SINE", "Retrotransposons Non-autonomous",
  "tRNA", "SINE", "Retrotransposons Non-autonomous"
) %>%
  mutate(across(everything(), as.factor))

# --- Input Data --- 
read_and_pivot_data <- function(filepath, file_level) {
  df <- read.csv(filepath, row.names = 1, check.names = FALSE) %>%
    rownames_to_column(var = "SpeciesID")

  # Filtering specific to Class level
  if (file_level == "Class") {
    cols_to_remove <- c("DNAtransposons Unknown", "Retrotransposons Unknown", "Other")
    cols_exist <- cols_to_remove[cols_to_remove %in% names(df)]
    if (length(cols_exist) > 0) {
        message(paste("Filtering Class columns:", paste(cols_exist, collapse=", ")))
        df <- df %>% select(-all_of(cols_exist))
    }
  }
  
  df_long <- df %>% 
    pivot_longer(
      cols = -SpeciesID,
      names_to = "Category",
      values_to = "Percentage"
    ) %>%
    filter(!is.na(Percentage) & Percentage > 1e-9)
  
  return(df_long)
}

# Special handling for superfamily filtering
read_filter_pivot_superfamily <- function(filepath) {
    df_wide <- read.csv(filepath, row.names = 1, check.names = FALSE) %>% 
        rownames_to_column(var = "SpeciesID")
        
    cols_with_zeros <- df_wide %>% 
        select(-SpeciesID) %>% 
        summarise(across(everything(), ~ any(. == 0, na.rm = TRUE))) %>% 
        pivot_longer(everything(), names_to = "Category", values_to = "has_zero") %>% 
        filter(has_zero) %>% 
        pull(Category)
        
    if (length(cols_with_zeros) > 0) {
        message(paste("Filtering Superfamily columns with zeros:", paste(cols_with_zeros, collapse=", ")))
        if ("Unknown TIR" %in% cols_with_zeros) {
             message("Excluding 'Unknown TIR' from Superfamily data as well.")
        }
        df_filtered_wide <- df_wide %>% select(-all_of(cols_with_zeros))
    } else {
        df_filtered_wide <- df_wide
        if ("Unknown TIR" %in% names(df_filtered_wide)) {
             message("Excluding 'Unknown TIR' from Superfamily data.")
             df_filtered_wide <- df_filtered_wide %>% select(-`Unknown TIR`)
        }
    }
    
    df_long <- df_filtered_wide %>% 
        pivot_longer(
            cols = -SpeciesID,
            names_to = "Category",
            values_to = "Percentage"
        ) %>% 
        filter(!is.na(Percentage) & Percentage > 1e-9)
        
    return(df_long)
}

# --- Load and process data ---
message("Loading and processing data...")
class_df <- read_and_pivot_data("results/data/dnaPipeTE_class_breakdown.csv", "Class")
order_df <- read_and_pivot_data("results/data/dnaPipeTE_order_breakdown.csv", "Order")
superfamily_df <- read_filter_pivot_superfamily("results/data/dnaPipeTE_superfamily_breakdown.csv")
message("Data loading complete.")

# --- Output Directory ---
output_dir <- "results/figures/TE_diversity_hierarchical_donuts"
if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
}

# --- Helper Function: Prepare Hierarchical Data ---
prepare_hierarchical_data <- function(species_id, c_df, o_df, s_df, h_map) {
    c_sub <- c_df %>% filter(SpeciesID == species_id)
    o_sub <- o_df %>% filter(SpeciesID == species_id)
    s_sub <- s_df %>% filter(SpeciesID == species_id)
    
    if (nrow(c_sub) == 0) return(NULL)
    class_data <- c_sub %>% 
        mutate(proportion = Percentage / sum(Percentage)) %>% 
        arrange(Category) %>% 
        mutate(
            Level = "Class", Class = factor(Category), Order = NA, Superfamily = NA,
            global_ymin = cumsum(lag(proportion, default = 0)), global_ymax = cumsum(proportion),
            xmin = 1.05, xmax = 1.95
        ) %>%
        mutate(fill_category = Class) %>%
        mutate(label_y = (global_ymin + global_ymax) / 2, label_x = (xmin + xmax) / 2) %>%
        select(SpeciesID, Level, Class, Order, Superfamily, fill_category, proportion, global_ymin, global_ymax, xmin, xmax, label_y, label_x)

    if (nrow(o_sub) == 0) {
        class_data_final <- class_data %>%
            rename(ymin = global_ymin, ymax = global_ymax) %>%
            select(SpeciesID, Level, Class, Order, Superfamily, fill_category, proportion, ymin, ymax, xmin, xmax, label_y, label_x) 
        return(class_data_final %>% mutate(fill_category = as.factor(fill_category)))
    }
    order_data <- o_sub %>%
        rename(Order = Category) %>%
        left_join(h_map %>% select(Order, Class) %>% distinct(), by = "Order", relationship = "many-to-many") %>%
        filter(!is.na(Class)) %>%
        left_join(class_data %>% select(Class, global_ymin, global_ymax), by = "Class") %>%
        filter(!is.na(global_ymin)) %>%
        group_by(Class, global_ymin, global_ymax) %>%
        mutate(prop_within_class = Percentage / sum(Percentage)) %>%
        arrange(Order) %>%
        mutate(
            Level = "Order", Superfamily = NA,
            ymin = global_ymin + cumsum(lag(prop_within_class, default = 0)) * (global_ymax - global_ymin),
            ymax = global_ymin + cumsum(prop_within_class) * (global_ymax - global_ymin),
            xmin = 2.05, xmax = 2.95
        ) %>%
        ungroup() %>%
        mutate(fill_category = Order) %>%
        mutate(label_y = (ymin + ymax) / 2, label_x = (xmin + xmax) / 2) %>%
        select(SpeciesID, Level, Class, Order, Superfamily, fill_category, proportion = prop_within_class, ymin, ymax, xmin, xmax, label_y, label_x)

    if (nrow(s_sub) == 0) {
         # Rename class data columns
         class_data_renamed <- class_data %>%
             rename(ymin = global_ymin, ymax = global_ymax) %>% 
             select(SpeciesID, Level, Class, Order, Superfamily, fill_category, proportion, ymin, ymax, xmin, xmax, label_y, label_x)
         # Order data already has proportion and label positions
         order_data_renamed <- order_data # No need to select again
         # Combine, ensure fill_category is factor
         combined_data <- bind_rows(class_data_renamed, order_data_renamed)
         return(combined_data %>% mutate(fill_category = as.factor(fill_category)))
    }
    superfamily_data <- s_sub %>%
        rename(Superfamily = Category) %>%
        left_join(h_map, by = "Superfamily", relationship = "many-to-many") %>%
        filter(!is.na(Order) & !is.na(Class)) %>%
        left_join(order_data %>% select(Class, Order, ymin, ymax) %>% rename(order_ymin=ymin, order_ymax=ymax),
                  by = c("Class", "Order")) %>%
        filter(!is.na(order_ymin)) %>%
        group_by(Class, Order, order_ymin, order_ymax) %>%
        mutate(prop_within_order = Percentage / sum(Percentage)) %>%
        arrange(Superfamily) %>%
        mutate(
            Level = "Superfamily",
            ymin = order_ymin + cumsum(lag(prop_within_order, default = 0)) * (order_ymax - order_ymin),
            ymax = order_ymin + cumsum(prop_within_order) * (order_ymax - order_ymin),
            xmin = 3.05, xmax = 3.95
        ) %>%
        ungroup() %>%
        mutate(fill_category = Superfamily) %>%
        mutate(label_y = (ymin + ymax) / 2, label_x = (xmin + xmax) / 2) %>%
        select(SpeciesID, Level, Class, Order, Superfamily, fill_category, proportion = prop_within_order, ymin, ymax, xmin, xmax, label_y, label_x)

    # Rename class data columns and select needed columns
    class_data_renamed <- class_data %>% 
        rename(ymin = global_ymin, ymax = global_ymax) %>% 
        select(SpeciesID, Level, Class, Order, Superfamily, fill_category, proportion, ymin, ymax, xmin, xmax, label_y, label_x)
    # Order data already has needed columns
    order_data_renamed <- order_data

    combined_data <- bind_rows(class_data_renamed, order_data_renamed, superfamily_data)
    # Ensure factor columns including the new fill_category
    combined_data <- combined_data %>% mutate(across(c(Level, Class, Order, Superfamily, fill_category), as.factor))

    return(combined_data)
}

# --- Plotting Function: Create Hierarchical Multi-Ring Donut ---
create_hierarchical_donut_plot <- function(plot_data, species_id, output_dir) {
    if (is.null(plot_data) || nrow(plot_data) == 0) {
        message(paste("No hierarchical data to plot for:", species_id))
        return(NULL)
    }

    # --- Define Manual Color Palette ---
    manual_palette <- c(
        # DNAtransposons Subclass1
        "DNAtransposons Subclass1" = "#23b8d9",
        "TIR" = "#1a90aa",
        "Academ" = "#106f83",
        "CACTA" = "#0898b6",
        "EnSpm" = "#06b7dc",
        "MULE" = "#48daf9",
        "Mutator" = "#98e4f5",
        "P" = "#53a0b2",
        "PIF-Harbinger" = "#286472",
        "PiggyBac" = "#044c5d",
        "Tc1-mariner" = "#0495b7",
        "hAT" = "#00cfff",
        "YR" = "#00aaff",
        "Cyrypton" = "#8ed9ff",
        
        # DNAtransposons Subclass2
        "DNAtransposons Subclass2" = "#23d991",
        "Helitron" = "#80efc3", # Also Superfamily
        "Maverick" = "#0d8a58", # Also Superfamily
        
        # Retrotransposons Autonomous
        "Retrotransposons Autonomous" = "#e81a1a",
        "DIRS" = "#680505", # Also Superfamily
        "LINE" = "#a93131",
        "Jockey" = "#c55151",
        "L1" = "#d13030",
        "LTR" = "#ff0000",
        "Gypsy" = "#ff0000", # Same as LTR per input
        "PLE" = "#ff9999", 
        "Penelope" = "#ff9999", # Same as PLE per input
        
        # Retrotransposons Non-autonomous
        "Retrotransposons Non-autonomous" = "#e8a31a",
        "SINE" = "#f8bd49",
        "tRNA" = "#f9ce7a",
        
        # Unknown
        "Unknown" = "#eebcf7",
        
        # Fallback for anything else unexpected
        "DEFAULT_FALLBACK" = "#808080" # Grey
    )
    
    # Get all unique categories present in the data for this species
    all_categories_in_data <- unique(as.character(plot_data$fill_category))
    
    # Select colors for the categories present, using fallback if needed
    final_colors_vec <- setNames(manual_palette[all_categories_in_data], all_categories_in_data)
    missing_from_palette <- is.na(final_colors_vec)
    if(any(missing_from_palette)){
        missing_names <- names(final_colors_vec)[missing_from_palette]
        warning(paste("Categories missing from manual palette:", paste(missing_names, collapse=", "), ". Using fallback grey."))
        final_colors_vec[missing_from_palette] <- manual_palette["DEFAULT_FALLBACK"]
    }
    # --- End Color Definition ---

    # Ensure fill_category factor levels match the final color names for reliable mapping
    # Re-level based on the names actually present and in the palette
    plot_data <- plot_data %>% 
        mutate(fill_category = factor(fill_category, levels = names(final_colors_vec)))

    max_xmax <- max(plot_data$xmax, na.rm = TRUE)
    xlim_upper <- ceiling(max_xmax) + 0.5

    plot <- ggplot(plot_data, aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax, fill = fill_category)) +
        geom_rect(color = "white", linewidth = 0.15, aes(alpha = Level)) +
        # Add text labels - filter for proportion > 0.01 (1%) to avoid clutter
        geom_text(data = . %>% filter(proportion > 0.01),
                  aes(x = label_x, y = label_y, label = sprintf("%.1f%%", proportion * 100)),
                  color = "white", size = 2.5, check_overlap = TRUE) +
        coord_polar(theta = "y", start = 0) +
        scale_fill_manual(values = final_colors_vec, 
                          guide = guide_legend(title = "Category", ncol=1),
                          na.value="#808080", # Use explicit fallback grey hex for NAs
                          drop = FALSE) + 
        scale_alpha_manual(values = c("Class" = 0.7, "Order" = 0.85, "Superfamily" = 1.0), guide="none") +
        xlim(0.5, xlim_upper) +
        theme_void() +
        theme(
            legend.position = "right",
            plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
            legend.key.size = unit(0.4, "cm"),
            legend.title = element_text(size=10, face="bold"),
            legend.text = element_text(size=7)
        ) +
        labs(
            title = paste(species_id, "- Hierarchical TE Composition")
        )

    filename <- file.path(output_dir, paste0(species_id, "_hierarchical_donut.png"))
    ggsave(filename, plot = plot, width = 12, height = 8, dpi = 300, bg = "white")
    message(paste("Saved:", filename))
}

# --- Main Execution ---

message("Starting main execution loop...")
if (!any(sapply(list(class_df, order_df, superfamily_df), function(df) "SpeciesID" %in% colnames(df)))) {
    stop("Error: 'SpeciesID' column not found after processing input data.")
}
species_list <- unique(c(class_df$SpeciesID, order_df$SpeciesID, superfamily_df$SpeciesID))
message(paste("Found", length(species_list), "species to process."))

# --- Plotting loop --- 
for (species in species_list) {
  message(paste("--- Processing:", species, "---"))
  hierarchical_data <- prepare_hierarchical_data(species, class_df, order_df, superfamily_df, hierarchy_map)
  create_hierarchical_donut_plot(hierarchical_data, species, output_dir) 
}

message("Hierarchical donut chart generation complete.") 