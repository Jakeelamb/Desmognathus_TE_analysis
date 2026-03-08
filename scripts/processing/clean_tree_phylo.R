script_args <- commandArgs(trailingOnly = FALSE)
script_path_arg <- grep("^--file=", script_args, value = TRUE)
if (length(script_path_arg) > 0) {
  script_path <- normalizePath(sub("^--file=", "", script_path_arg[1]), mustWork = FALSE)
  script_dir <- dirname(script_path)
} else {
  script_dir <- "scripts/processing"
}

source(file.path(dirname(script_dir), "R", "path_config_utils.R"))
prefer_active_conda_r_library()

# --- Libraries ---
suppressPackageStartupMessages({
  library(ape)
  library(readr)
  library(stringr)
  library(fs)
  library(dplyr)
})

# --- Configuration ---
BASE_DIR <- path_norm(find_project_root(script_dir))
config <- load_project_config(BASE_DIR)

INPUT_DATA_DIR <- path_norm(resolve_config_path(BASE_DIR, config$input_data$root, "input_data"))
OUTPUT_DATA_DIR <- path_norm(resolve_config_path(BASE_DIR, config$results$data, "results/data"))
FIGURE_DIR <- path_norm(resolve_config_path(BASE_DIR, config$results$figures$phylogeny, "results/figures/phylogeny"))

LOOKUP_PATH <- path_norm(resolve_config_path(BASE_DIR, config$input_data$lookup_table %||% config$data$lookup, "input_data/lookup_table.txt"))
INPUT_TREE_PATH <- path(INPUT_DATA_DIR, "phylogeny", "desmo900dated_test.tre")
OUTPUT_TREE_PATH <- path(OUTPUT_DATA_DIR, "desmo900dated_test_cleaned_phylo.tre")

dir_create(path_dir(OUTPUT_TREE_PATH))
dir_create(FIGURE_DIR)

message("--- Starting Tree Cleaning Script (R version) ---")

# --- Functions ---

#' Load Species Names from Lookup Table
#'
#' Reads the lookup table, extracts species names from the first column,
#' and cleans them by removing the 'D.' prefix and trimming whitespace.
#'
#' @param lookup_path Path to the lookup table file.
#' @return A character vector of unique, cleaned species names.
load_lookup_species <- function(lookup_path) {
  message(paste("Loading lookup table from:", lookup_path))
  tryCatch({
    # Read only the first column, assuming tab-separated and has header
    lookup_df <- read_tsv(lookup_path, col_names = TRUE, show_col_types = FALSE)
    
    if (ncol(lookup_df) == 0) {
      stop("Lookup table appears to be empty or incorrectly formatted.")
    }
    
    first_col_name <- names(lookup_df)[1]
    message(paste("Loaded", nrow(lookup_df), "rows from", lookup_path, 
                       ", targeting first column '", first_col_name, "'."))
                       
    # Extract the first column, convert to character, handle potential factors
    species_raw <- as.character(lookup_df[[first_col_name]])

    # Clean names: remove "D.", trim whitespace, filter out empty strings
    valid_species <- species_raw %>%
      str_trim() %>%                     # Trim whitespace
      str_remove("^D\\.\\s*") %>%        # Remove "D." prefix (and any following space)
      str_trim() %>%                     # Trim again after removal
      .[. != ""] %>%                     # Remove empty strings
      unique()                           # Get unique names

    message(paste("Loaded and cleaned", length(valid_species), "unique species names."))
    message(paste("Cleaned Lookup Species Names:", paste(sort(valid_species), collapse=", "))) # Log cleaned names
    
    return(valid_species)

  }, error = function(e) {
    # Using warning here, and then stop, similar to original logic
    warning(paste("Error loading or processing lookup table:", lookup_path, "-", e$message))
    stop("Failed to load lookup species.", call. = FALSE) # Stop script execution
  })
}


#' Find Best Matching Species Name
#'
#' Finds the longest valid species name (from the lookup list) that the
#' tree tip label starts with. Comparison is case-insensitive and ignores
#' leading/trailing whitespace.
#'
#' @param label The tree tip label (string).
#' @param valid_names A character vector of cleaned, valid species names.
#' @return The best matching valid species name (original casing from valid_names),
#'         or NULL if no match is found.
find_best_match <- function(label, valid_names) {
  if (is.null(label) || is.na(label) || str_trim(label) == "") {
    return(NULL)
  }
  
  label_processed <- str_trim(label) %>% str_to_lower()
  best_match <- NULL
  max_len <- 0
  
  for (valid_name in valid_names) {
    # Assuming valid_names are already cleaned (trimmed, no "D.")
    valid_name_processed <- str_trim(valid_name) %>% str_to_lower() # Lowercase for comparison
    
    if (str_starts(label_processed, fixed(valid_name_processed))) {
      current_match_len <- nchar(valid_name) # Use length of original valid_name
      if (current_match_len > max_len) {
        max_len <- current_match_len
        best_match <- valid_name # Return the original casing
      }
    }
  }
  
  return(best_match)
}

# Function to check if a label matches lookup table exactly
is_in_lookup <- function(label, lookup_names) {
  # Add "D." prefix to the label for comparison
  label_with_prefix <- paste0("D.", label)
  # Check if the prefixed label exists in lookup table
  label_with_prefix %in% lookup_names
}

# Function to remove suffixes (e.g., _1, _2, etc.)
remove_suffixes <- function(label) {
  # Remove any suffix that starts with underscore followed by numbers
  str_remove(label, "_\\d+$")
}

# --- Main Script Logic ---

# Ensure output directory exists (Moved up before first message)
# tryCatch({
#   dir_create(path_dir(OUTPUT_TREE_PATH))
#   message(paste("Ensured output directory exists:", path_dir(OUTPUT_TREE_PATH)))
# }, error = function(e){
#   warning(paste("Failed to create output directory:", path_dir(OUTPUT_TREE_PATH), "-", e$message))
#   stop("Directory creation failed.", call. = FALSE)
# })


# 1. Load valid species names
valid_species_names <- load_lookup_species(LOOKUP_PATH)

# 2. Load the input tree
message(paste("Loading tree from:", INPUT_TREE_PATH))
tree <- tryCatch({
  tree <- read.tree(INPUT_TREE_PATH)
  # Check if tree has branch lengths
  if (is.null(tree$edge.length)) {
    stop("Tree does not contain branch lengths!")
  }
  # Check for ultrametric tree (all tips at same distance from root)
  if (!is.ultrametric(tree)) {
    warning("Tree is not ultrametric - time scale may be inconsistent")
  }
  tree
}, error = function(e) {
  warning(paste("Failed to load or parse tree", INPUT_TREE_PATH, ":", e$message))
  stop("Tree loading failed.", call. = FALSE)
})
message(paste("Successfully loaded tree with", length(tree$tip.label), "tips."))

original_labels <- tree$tip.label
message(paste("Original Tree Tip Labels:", paste(sort(original_labels), collapse=", ")))


# 3. Clean tip labels and identify tips to keep
message("Cleaning tip labels and identifying valid species...")
original_labels <- tree$tip.label
cleaned_labels <- character(length(original_labels))
tips_to_keep <- logical(length(original_labels))  # Track which tips to keep
label_map <- list()

# First, get the raw lookup table names (with D. prefix)
lookup_df <- read_tsv(LOOKUP_PATH, col_names = TRUE, show_col_types = FALSE)
lookup_names <- as.character(lookup_df[[1]])  # Get first column with D. prefix

for (i in seq_along(original_labels)) {
  original_label <- original_labels[i]
  # Remove any suffixes
  base_label <- remove_suffixes(original_label)
  
  # Check if this label (with D. prefix) exists in lookup table
  if (is_in_lookup(base_label, lookup_names)) {
    cleaned_labels[i] <- base_label
    tips_to_keep[i] <- TRUE
    label_map[[original_label]] <- base_label
    message(paste0("Keeping valid species: ", base_label))
  } else {
    warning(paste0("Removing invalid species: ", original_label))
    tips_to_keep[i] <- FALSE
  }
}

# Apply the cleaned labels to the tree object
tree$tip.label <- cleaned_labels

# 4. Prune invalid tips
message("Pruning invalid tips...")
tips_to_remove <- which(!tips_to_keep)
if (length(tips_to_remove) > 0) {
  message(paste("Removing", length(tips_to_remove), "invalid tips"))
  tree_cleaned <- drop.tip(tree, tip = tips_to_remove)
} else {
  message("No invalid tips to remove")
  tree_cleaned <- tree
}

# 5. Visualize the cleaned tree
message("Visualizing cleaned tree...")
tryCatch({
  # Create a PNG file for the tree plot
  plot_file <- path(FIGURE_DIR, "rectangular_phylogeny.png")
  png(plot_file, width = 1200, height = 1800, res = 150)
  
  # Set up the plot with proper margins
  par(mar = c(5, 6, 4, 2))  # Increase left margin for labels
  
  # Calculate root to tip distance and time points
  root_to_tip <- max(node.depth.edgelength(tree_cleaned))
  message(paste("Maximum root to tip distance:", root_to_tip))
  
  # Plot the tree with rectangular layout
  plot(tree_cleaned, 
       type = "phylogram",
       show.tip.label = TRUE,
       cex = 0.8,           # Increase label size
       label.offset = 0.5,  # Increase label offset
       edge.width = 1.5,
       edge.color = "black",
       no.margin = FALSE,
       x.lim = c(-2, root_to_tip * 1.2))  # Add more padding on both sides
  
  # Add time scale with reversed values (0 at tips)
  time_points <- pretty(c(0, root_to_tip), n = 10)  # Get time points in original scale
  axis(1, 
       at = time_points, 
       labels = round(rev(time_points), 1),  # Reverse the labels so 0 is at tips
       las = 1)
  
  # Add scale bar at a better position, with reversed direction
  add.scale.bar(x = root_to_tip,  # Position at right side (tips)
                y = -1,           # Position below tree
                length = -5,      # Negative length to go left (backwards in time)
                cex = 0.8,        # Text size
                lwd = 1.5)        # Line width
  
  # Add title and axis label
  title(main = "Cleaned Desmognathus Phylogeny",
        xlab = "Time (Million Years Ago)")
  
  # Close the PNG device
  dev.off()
  message(paste("Tree visualization saved to:", plot_file))
}, error = function(e) {
  warning(paste("Failed to create tree visualization:", e$message))
})

# 6. Write the cleaned tree
message(paste("Writing cleaned tree to:", OUTPUT_TREE_PATH))
tryCatch({
  write.tree(tree_cleaned, file = OUTPUT_TREE_PATH)
  message("Successfully wrote cleaned tree.")
}, error = function(e) {
  warning(paste("Failed to write cleaned tree:", e$message))
  stop("Failed to write cleaned tree.", call. = FALSE)
})

message("--- Tree Cleaning Script Finished ---") 
