#!/usr/bin/env Rscript
#
# TE_PCA_ANALYSIS.R
#
# Description: Performs Principal Component Analysis (PCA) on transposable element data
# using the superfamily proportions dataset and visualization approach from Master_data_vis_PCA.R.
#

# Load necessary libraries
suppressPackageStartupMessages({
  library(ggrepel)
  library(ggplot2)
  library(Rtsne)
  library(umap)
  library(plotly)
  library(cluster)
  library(factoextra)
  library(gridExtra)
  library(viridis)
  library(optparse)
  library(dplyr)
})

# Parse command line arguments
option_list <- list(
  make_option(c("-i", "--input"), type="character", default="data/diversity/superfamily_proportions.csv",
              help="Input CSV file with superfamily proportions [default=%default]", metavar="character"),
  make_option(c("-o", "--output_dir"), type="character", default="data/results/pca",
              help="Output directory for PCA results [default=%default]", metavar="character"),
  make_option(c("-f", "--filter_species"), type="character", default="fuscus,folkertsi,ochrophaeus,imitator",
              help="Comma-separated list of species to remove [default=%default]", metavar="character"),
  make_option(c("-k", "--optimal_k"), type="integer", default=5,
              help="Optimal number of clusters for K-means [default=%default]", metavar="integer")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Create output directory if it doesn't exist
dir.create(opt$output_dir, recursive = TRUE, showWarnings = FALSE)

#---------------------- Helper Functions ----------------------#

# Function to clean column names
clean_colnames <- function(df) {
  cleaned_names <- gsub("_Percentage_of_sequence", "", colnames(df))
  colnames(df) <- cleaned_names
  return(df)
}

# Function to check if a column has variance
has_variance <- function(x) var(x, na.rm = TRUE) > 0

# Function to get top features
get_top_features <- function(pca_obj, pc_num, top_n = 3) {
  loadings <- pca_obj$rotation[, pc_num]
  sorted_loadings <- sort(abs(loadings), decreasing = TRUE)
  top_features <- names(sorted_loadings)[1:min(top_n, length(sorted_loadings))]
  return(paste(top_features, collapse = ", "))
}

# Publication theme
theme_publication <- theme_minimal() +
  theme(
    text = element_text(family = "serif", size = 12),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.position = "none"  # Remove legend
  )

#---------------------- Data Preparation ----------------------#

# Load data
cat("Loading TE data from:", opt$input, "\n")
df <- read.csv(opt$input)

# Filter out specified species if requested
if (!is.null(opt$filter_species) && opt$filter_species != "") {
  species_to_remove <- unlist(strsplit(opt$filter_species, ","))
  cat("Removing species:", paste(species_to_remove, collapse=", "), "\n")
  df_filtered <- df[!df$Species %in% species_to_remove, ]
} else {
  df_filtered <- df
}

# Remove "Unknown.TIR" column if it exists
if ("Unknown.TIR" %in% colnames(df_filtered)) {
  cat("Removing Unknown.TIR column\n")
  df_filtered <- df_filtered[, colnames(df_filtered) != "Unknown.TIR"]
}

# Process numeric columns
numeric_cols <- sapply(df_filtered, is.numeric)
df_numeric <- df_filtered[, numeric_cols]

# Remove constant/zero variance columns
df_with_variance <- df_numeric[, sapply(df_numeric, has_variance)]

# Handle index column if present
if ("X" %in% colnames(df_with_variance)) {
  df_with_variance <- df_with_variance[, colnames(df_with_variance) != "X"]
}

# Clean column names if needed
df_with_variance <- clean_colnames(df_with_variance)

cat("Data prepared with", ncol(df_with_variance), "TE features and", nrow(df_filtered), "samples\n")

#---------------------- PCA Analysis ----------------------#

analyze_pca_components <- function(data) {
  # Perform PCA
  cat("Performing PCA analysis...\n")
  pca_result <- prcomp(data, scale = TRUE)
  
  # Calculate variance explained
  variance_explained <- pca_result$sdev^2 / sum(pca_result$sdev^2)
  cumulative_variance <- cumsum(variance_explained)
  
  # Create dataframe for variance information
  pca_data <- data.frame(
    PC = 1:length(pca_result$sdev),
    Eigenvalue = pca_result$sdev^2,
    VarExplained = variance_explained,
    CumVarExplained = cumulative_variance
  )
  
  # Kaiser criterion
  kaiser_components <- sum(pca_result$sdev^2 > 1)
  
  # Create scree plot
  scree_plot <- ggplot() +
    geom_line(data = pca_data, aes(x = PC, y = VarExplained), color = "blue") +
    geom_point(data = pca_data, aes(x = PC, y = VarExplained), color = "blue") +
    geom_line(data = pca_data, aes(x = PC, y = CumVarExplained), color = "red", linetype = "dashed") +
    geom_hline(yintercept = c(0.8, 0.9), color = "darkgreen", linetype = "dashed") +
    scale_y_continuous(
      name = "Variance Explained",
      sec.axis = sec_axis(~., name = "Cumulative Variance Explained")
    ) +
    labs(title = "Scree Plot with Cumulative Variance", x = "Principal Component") +
    theme_publication +
    annotate("text", x = max(pca_data$PC)/2, y = c(0.82, 0.92), 
             label = c("80% threshold", "90% threshold"), color = "darkgreen")
  
  # Find important thresholds
  components_80 <- which(cumulative_variance >= 0.8)[1]
  components_90 <- which(cumulative_variance >= 0.9)[1]
  
  return(list(
    pca_result = pca_result,
    variance_explained = variance_explained,
    cumulative_variance = cumulative_variance,
    scree_plot = scree_plot,
    kaiser_components = kaiser_components,
    components_80 = components_80,
    components_90 = components_90,
    pca_data = pca_data
  ))
}

# Run the PCA analysis
pca_analysis <- analyze_pca_components(df_with_variance)
pca_result <- pca_analysis$pca_result

# Get the first 6 principal components
pc_data <- as.data.frame(pca_result$x[, 1:min(6, ncol(pca_result$x))])
pc_data$Species <- df_filtered$Species

#---------------------- Visualizations ----------------------#

# Create PCA plot by species with top feature labels in axis
pca_species_plot <- ggplot(pc_data, aes(x = PC1, y = PC2, color = Species)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_text_repel(aes(label = Species), size = 3, max.overlaps = 15) +
  labs(
    title = "PCA Plot TE Composition by Species",
    x = paste0("PC1 (", round(100 * pca_analysis$variance_explained[1], 1), "% variance)\nTop features: ", get_top_features(pca_result, 1)),
    y = paste0("PC2 (", round(100 * pca_analysis$variance_explained[2], 1), "% variance)\nTop features: ", get_top_features(pca_result, 2))
  ) +
  theme_publication

#---------------------- Clustering Analysis ----------------------#

find_optimal_clusters <- function(data, max_k = 10) {
  cat("Performing clustering analysis...\n")
  
  # WSS calculation (within-cluster sum of squares)
  wss <- sapply(1:max_k, function(k) {
    kmeans(data, centers = k, nstart = 25)$tot.withinss
  })
  
  # Silhouette calculation
  sil_scores <- sapply(2:max_k, function(k) {
    km <- kmeans(data, centers = k, nstart = 25)
    mean(silhouette(km$cluster, dist(data))[,3])
  })
  
  # Elbow plot
  elbow_plot <- ggplot(data.frame(k = 1:max_k, wss = wss), aes(x = k, y = wss)) +
    geom_line() +
    geom_point() +
    labs(title = "Elbow Plot for Optimal k",
         x = "Number of Clusters (k)",
         y = "Total Within-cluster Sum of Squares") +
    theme_publication
  
  # Silhouette plot
  sil_plot <- ggplot(data.frame(k = 2:max_k, score = sil_scores), aes(x = k, y = score)) +
    geom_line() +
    geom_point() +
    labs(title = "Silhouette Scores",
         x = "Number of Clusters (k)",
         y = "Average Silhouette Score") +
    theme_publication
  
  return(list(
    elbow_plot = elbow_plot,
    silhouette_plot = sil_plot,
    silhouette_scores = sil_scores,
    wss = wss
  ))
}

# Run clustering analysis on 6 PCs
cluster_analysis <- find_optimal_clusters(as.matrix(pc_data[, 1:6]))

# Perform final clustering with optimal k
optimal_k <- opt$optimal_k
cat("Performing K-means clustering with k =", optimal_k, "\n")
kmeans_result <- kmeans(as.matrix(pc_data[, 1:6]), centers = optimal_k, nstart = 25)
pc_data$Cluster <- as.factor(kmeans_result$cluster)

# Create final cluster plot
pca_cluster_plot <- ggplot(pc_data, aes(x = PC1, y = PC2, color = Cluster)) +
  geom_point(size = 3, alpha = 0.6) +
  stat_ellipse(level = 0.95) +
  labs(title = "PCA Plot with Clusters",
       x = paste0("PC1 (", round(100 * pca_analysis$variance_explained[1], 1), "%)"),
       y = paste0("PC2 (", round(100 * pca_analysis$variance_explained[2], 1), "%)")) +
  theme_publication +
  scale_color_viridis_d() +
  geom_text_repel(aes(label = Species), size = 3, max.overlaps = 15)

#---------------------- Save Results ----------------------#

# Save plots
cat("Saving PCA plots to:", opt$output_dir, "\n")
ggsave(file.path(opt$output_dir, "pca_scree_plot.png"), pca_analysis$scree_plot, 
       width = 10, height = 8, dpi = 300)
ggsave(file.path(opt$output_dir, "pca_by_species.png"), pca_species_plot, 
       width = 10, height = 8, dpi = 300)
ggsave(file.path(opt$output_dir, "cluster_elbow_plot.png"), cluster_analysis$elbow_plot, 
       width = 10, height = 8, dpi = 300)
ggsave(file.path(opt$output_dir, "cluster_silhouette_plot.png"), cluster_analysis$silhouette_plot, 
       width = 10, height = 8, dpi = 300)
ggsave(file.path(opt$output_dir, "pca_cluster_plot.png"), pca_cluster_plot, 
       width = 10, height = 8, dpi = 300)

# Save variance information
variance_summary <- data.frame(
  PC = 1:length(pca_analysis$variance_explained),
  Eigenvalue = pca_result$sdev^2,
  Variance_Explained = pca_analysis$variance_explained * 100,
  Cumulative_Variance = pca_analysis$cumulative_variance * 100
)
write.csv(variance_summary, file.path(opt$output_dir, "pca_variance_summary.csv"), row.names = FALSE)

# Save PC scores
pc_scores <- as.data.frame(pca_result$x)
pc_scores$Species <- df_filtered$Species
pc_scores$Cluster <- pc_data$Cluster
write.csv(pc_scores, file.path(opt$output_dir, "pca_scores.csv"), row.names = FALSE)

# Save loadings
loadings_df <- as.data.frame(pca_result$rotation)
loadings_df$Feature <- rownames(loadings_df)
write.csv(loadings_df, file.path(opt$output_dir, "pca_loadings.csv"), row.names = FALSE)

# Save cluster information
cluster_info <- data.frame(
  Species = pc_data$Species,
  Cluster = pc_data$Cluster
)
write.csv(cluster_info, file.path(opt$output_dir, "cluster_assignments.csv"), row.names = FALSE)

# Save WSS and silhouette scores
clustering_metrics <- data.frame(
  k = 1:10,
  WSS = c(cluster_analysis$wss),
  Silhouette = c(NA, cluster_analysis$silhouette_scores)
)
write.csv(clustering_metrics, file.path(opt$output_dir, "clustering_metrics.csv"), row.names = FALSE)

#---------------------- Print Summary ----------------------#

cat("\nPCA Analysis Summary:\n")
cat("---------------------\n")
cat("Total features analyzed:", ncol(df_with_variance), "\n")
cat("Total samples:", nrow(df_filtered), "\n")
cat("\nVariance Explained:\n")
for (i in 1:min(5, length(pca_analysis$variance_explained))) {
  cat(sprintf("PC%d: %.2f%% (cumulative: %.2f%%)\n", 
              i, 
              pca_analysis$variance_explained[i] * 100,
              pca_analysis$cumulative_variance[i] * 100))
}

cat("\nComponent Selection Criteria:\n")
cat("Kaiser criterion (eigenvalue > 1):", pca_analysis$kaiser_components, "components\n")
cat("Components needed for 80% variance:", pca_analysis$components_80, "\n")
cat("Components needed for 90% variance:", pca_analysis$components_90, "\n")

cat("\nTop Contributing Features:\n")
for (i in 1:min(3, ncol(pca_result$rotation))) {
  cat(sprintf("PC%d: %s\n", i, get_top_features(pca_result, i, 5)))
}

cat("\nClustering Results:\n")
cat("Optimal k used for clustering:", optimal_k, "\n")
cat("Cluster distribution:\n")
print(table(pc_data$Cluster))

cat("\nSilhouette scores for different k values:\n")
for (k in 2:10) {
  cat(sprintf("k=%d: %.3f\n", k, cluster_analysis$silhouette_scores[k-1]))
}

cat("\nAll visualizations and data files saved to:", opt$output_dir, "\n") 