library(ape)
library(phytools)
library(dplyr)
library(here)
library(yaml)
library(ggtree)
library(ggplot2)
library(RColorBrewer)

# === Configuration === 
project_root <- here::here()
config_file <- file.path(project_root, "config/paths.yaml")

# Set default paths
default_tree_path <- "results/phylogeny/processed_phylogeny.nwk"
default_output_subdir <- "results/tables/phylogeny"
default_figure_subdir <- "results/figures/phylogeny"

# Initialize with defaults
tree_path <- file.path(project_root, default_tree_path)
output_dir <- file.path(project_root, default_output_subdir)
figure_dir <- file.path(project_root, default_figure_subdir)
output_csv_file <- file.path(output_dir, "clade_assignments.csv")

if (file.exists(config_file)) {
  tryCatch({
    config <- yaml::read_yaml(config_file)
    cat("Config file found. Reading paths.\n")
    
    tree_path_from_config <- file.path(config$results$phylogeny, "processed_phylogeny.nwk")
    if (!is.null(tree_path_from_config) && is.character(tree_path_from_config) && nzchar(tree_path_from_config)) {
      tree_path <- file.path(project_root, tree_path_from_config)
      cat("  Using tree path from config.\n")
    } else {
      warning("  Config missing or invalid 'results$phylogeny'. Using default tree path.")
    }
    
    output_dir_from_config <- config$results$tables$phylogeny
    if (!is.null(output_dir_from_config) && is.character(output_dir_from_config) && nzchar(output_dir_from_config)) {
      output_dir <- file.path(project_root, output_dir_from_config)
      output_csv_file <- file.path(output_dir, "clade_assignments.csv")
      cat("  Using output directory from config.\n")
    } else {
      warning("  Config missing or invalid 'results$tables$phylogeny'. Using default output dir.")
    }
    
    figure_dir_from_config <- config$results$figures$phylogeny
    if (!is.null(figure_dir_from_config) && is.character(figure_dir_from_config) && nzchar(figure_dir_from_config)) {
      figure_dir <- file.path(project_root, figure_dir_from_config)
      cat("  Using figure directory from config.\n")
    } else {
      warning("  Config missing or invalid 'results$figures$phylogeny'. Using default figure dir.")
    }
    
  }, error = function(e) {
      warning("Error reading config file: ", conditionMessage(e), ". Using defaults.")
  })
} else {
    warning("config/paths.yaml not found. Using default paths.")
}

# Ensure output directories exist
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

# Load your time-calibrated tree
cat("Loading tree from:", tree_path, "\n")
tree <- read.tree(tree_path)

# Set threshold in units matching your tree (assumes branch lengths = time in MY)
threshold_mya <- 5

# === Define Clades ===
cat("Defining clades based on", threshold_mya, "MYA threshold...\n")

# Get all internal nodes
internal_nodes <- (length(tree$tip.label) + 1):(length(tree$tip.label) + tree$Nnode)

# Calculate node heights
node_depths <- node.depth.edgelength(tree)
root_height <- max(node_depths)
node_heights <- root_height - node_depths

# Find nodes that define clades (nodes older than threshold)
old_nodes <- internal_nodes[node_heights[internal_nodes] > threshold_mya]

# Sort nodes by height (deepest first)
old_nodes <- old_nodes[order(node_heights[old_nodes], decreasing = TRUE)]

# Initialize species data frame
all_species <- data.frame(
  species = tree$tip.label,
  clade = NA_character_,
  stringsAsFactors = FALSE
)

# Function to get descendant tips
get_descendants <- function(node, tree) {
  if (node <= length(tree$tip.label)) {
    return(tree$tip.label[node])
  }
  clade <- extract.clade(tree, node)
  return(clade$tip.label)
}

# Process each old node to find clades
clade_counter <- 1
assigned_tips <- character(0)

# First, identify monophyletic groups
for (node in old_nodes) {
  # Get descendants
  tips <- get_descendants(node, tree)
  
  # Check if these tips form a new group
  unassigned_tips <- setdiff(tips, assigned_tips)
  
  # If we found unassigned tips and the node is deep enough
  if (length(unassigned_tips) > 1) {  # Only create clade if more than one species
    # Create new clade
    all_species$clade[all_species$species %in% unassigned_tips] <- paste0("Clade_", clade_counter)
    assigned_tips <- c(assigned_tips, unassigned_tips)
    clade_counter <- clade_counter + 1
  }
}

# Assign any remaining species to their own clades
unassigned <- all_species$species[is.na(all_species$clade)]
if (length(unassigned) > 0) {
  for (sp in unassigned) {
    all_species$clade[all_species$species == sp] <- paste0("Clade_", clade_counter)
    clade_counter <- clade_counter + 1
  }
}

# Create final clade assignments
clade_assignments <- all_species[, c("species", "clade")]

# Sort by species name
clade_assignments <- clade_assignments[order(clade_assignments$species), ]

# Print results
cat("\nClade assignments:\n")
print(clade_assignments)

# === Generate Phylogenetic Tree Visualization ===
cat("\nGenerating phylogenetic tree visualization...\n")

# Create a color palette for clades
n_clades <- length(unique(clade_assignments$clade))
if (n_clades <= 8) {
  clade_colors <- brewer.pal(n_clades, "Dark2")  # Using Dark2 for better visibility
} else {
  clade_colors <- colorRampPalette(brewer.pal(8, "Dark2"))(n_clades)
}

# Create named vector for colors
names(clade_colors) <- sort(unique(clade_assignments$clade))

# Create the base tree
p <- ggtree(tree, layout = "rectangular", size = 0.5)

# Get the data with tip labels
p_data <- p$data

# Add clade information to the data
p_data$clade <- NA
p_data$label_with_clade <- p_data$label  # Create new column for modified labels

# Add clade information and modify labels
for (i in seq_len(nrow(clade_assignments))) {
  idx <- which(p_data$label == clade_assignments$species[i])
  if (length(idx) > 0) {
    p_data$clade[idx] <- clade_assignments$clade[i]
    clade_num <- as.numeric(sub("Clade_", "", clade_assignments$clade[i]))
    p_data$label_with_clade[idx] <- paste0(p_data$label[idx], " (", clade_num, ")")
  }
}

# Update the tree with the new data
p$data <- p_data

# Add the basic tree elements
p <- p + 
  geom_tree(size = 0.5) +
  geom_tippoint(aes(color = clade), size = 2) +
  geom_tiplab(aes(label = label_with_clade), size = 3, hjust = -0.1) +
  scale_color_manual(values = clade_colors, guide = "none") +
  theme_tree2() +
  theme(
    axis.text.x = element_text(size = 10),
    axis.title.x = element_text(size = 12, face = "bold"),
    plot.margin = margin(20, 200, 20, 20)  # Increased right margin
  ) +
  xlim(NA, 27) +  # Allow natural left limit, extend right limit
  xlab("Time (MYA)")

# Add clade bars and labels
for (i in seq_along(old_nodes)) {
  clade_members <- extract.clade(tree, old_nodes[i])$tip.label
  if (length(clade_members) > 0) {
    clade_name <- paste0("Clade_", i)
    if (clade_name %in% names(clade_colors)) {
      p <- p + 
        geom_cladelabel(
          node = old_nodes[i],
          label = paste("Clade", i),
          color = clade_colors[clade_name],
          offset.text = 25,  # Increased offset for labels
          barsize = 2,
          fontsize = 3,
          hjust = 0,
          angle = 0,
          align = TRUE,
          offset = 0.5
        )
    }
  }
}

# Save the plot with higher resolution and expanded width
output_plot_file <- file.path(figure_dir, "clade_phylogeny.png")
ggsave(
  output_plot_file, 
  p, 
  width = 15,  # Increased width for better spacing
  height = 8, 
  dpi = 300,
  bg = "white"
)

# === Save Results ===
cat("Saving clade assignments to:", output_csv_file, "\n")
write.csv(clade_assignments, file = output_csv_file, row.names = FALSE)

cat("Script finished.\n")
