options(repos = c(CRAN = "https://cloud.r-project.org"))

install.packages("ggrepel")
install.packages("ggplot2")
install.packages("Rtsne")
install.packages("umap")
install.packages("plotly")
install.packages("cluster")
install.packages("factoextra")
install.packages("gridExtra")

# Load necessary libraries
library(ggrepel)
library(ggplot2)
library(Rtsne)
library(umap)
library(plotly)
library(cluster)
library(factoextra)
library(gridExtra)

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
  top_features <- names(sorted_loadings)[1:top_n]
  return(paste(top_features, collapse = ", "))
}

# Publication theme
theme_publication <- theme_minimal() +
  theme(
    text = element_text(family = "serif", size = 12),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.position = "none"  # Removed legend
  )

#---------------------- Data Preparation ----------------------#
setwd("/home/jake/Projects/Desmognathus_RPC_Results")
# Load and filter data
df <- read.csv("/home/jake/Projects/dnaPipeTE/superfamily_proportions.csv")
species_to_remove <- c("fuscus", "folkertsi", "ochrophaeus", "imitator")
df_filtered <- df[!df$Species %in% species_to_remove, ]

# Remove specified columns
columns_to_remove <- c("Unknown.TIR")
df_filtered <- df_filtered[, !colnames(df_filtered) %in% columns_to_remove]

# Process numeric columns
numeric_cols <- sapply(df_filtered, is.numeric)
df_numeric <- df_filtered[, numeric_cols]
df_with_variance <- df_numeric[, sapply(df_numeric, has_variance)]
df_with_variance <- df_with_variance[, -1]
df_with_variance <- clean_colnames(df_with_variance)

#---------------------- PCA Analysis ----------------------#

analyze_pca_components <- function(data) {
  pca_result <- prcomp(data, scale = TRUE)
  variance_explained <- pca_result$sdev^2 / sum(pca_result$sdev^2)
  cumulative_variance <- cumsum(variance_explained)
  
  pca_data <- data.frame(
    PC = 1:length(pca_result$sdev),
    Eigenvalue = pca_result$sdev^2,
    VarExplained = variance_explained,
    CumVarExplained = cumulative_variance
  )
  
  # Kaiser criterion
  kaiser_components <- sum(pca_result$sdev^2 > 1)
  
  # Scree plot
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
  
  components_80 <- which(cumulative_variance >= 0.8)[1]
  components_90 <- which(cumulative_variance >= 0.9)[1]
  
  return(list(
    pca_result = pca_result,
    variance_explained = variance_explained,
    cumulative_variance = cumulative_variance,
    scree_plot = scree_plot,
    kaiser_components = kaiser_components,
    components_80 = components_80,
    components_90 = components_90
  ))
}

#---------------------- Run PCA ----------------------#

pca_analysis <- analyze_pca_components(df_with_variance)
pca_result <- pca_analysis$pca_result

# Get first 6 principal components
pc_data <- as.data.frame(pca_result$x[, 1:6])
pc_data$Species <- df_filtered$Species

#---------------------- Visualizations ----------------------#

# Create and save PCA species plot
pca_species_plot <- ggplot(pc_data, aes(x = PC1, y = PC2, color = Species)) +
  geom_point() +
  geom_text_repel(aes(label = Species), size = 3, nudge_x = 0.2, nudge_y = 0.2) +
  labs(
    title = "PCA Plot TE Composition by Species",
    x = paste0("PC1 (", round(100 * pca_analysis$variance_explained[1], 1), "% variance)\nTop features: ", get_top_features(pca_result, 1)),
    y = paste0("PC2 (", round(100 * pca_analysis$variance_explained[2], 1), "% variance)\nTop features: ", get_top_features(pca_result, 2))
  ) +
  theme_publication

#---------------------- Clustering Analysis ----------------------#

find_optimal_clusters <- function(data, max_k = 10) {
  # WSS calculation
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

#---------------------- Save Results ----------------------#

# Save all plots
ggsave("scree_plot_enhanced.png", pca_analysis$scree_plot, width = 10, height = 8)
ggsave("pca_plot_species_with_features.png", pca_species_plot, width = 10, height = 8)
ggsave("cluster_elbow_plot.png", cluster_analysis$elbow_plot, width = 10, height = 8)
ggsave("cluster_silhouette_plot.png", cluster_analysis$silhouette_plot, width = 10, height = 8)

# Save variance summary
variance_summary <- data.frame(
  PC = 1:length(pca_analysis$variance_explained),
  Variance_Explained = round(pca_analysis$variance_explained * 100, 2),
  Cumulative_Variance = round(pca_analysis$cumulative_variance * 100, 2)
)
write.csv(variance_summary, "pca_variance_summary.csv", row.names = FALSE)

# Print summary information
cat("\nPCA Component Analysis:\n")
cat("Kaiser criterion (eigenvalue > 1):", pca_analysis$kaiser_components, "components\n")
cat("Components needed for 80% variance:", pca_analysis$components_80, "\n")
cat("Components needed for 90% variance:", pca_analysis$components_90, "\n")

# Print silhouette scores
cat("\nSilhouette Scores for k=2 to k=10:\n")
print(round(cluster_analysis$silhouette_scores, 3))

#---------------------- Optional: Perform Final Clustering ----------------------#
# Uncomment and modify optimal_k based on your analysis of the elbow and silhouette plots

optimal_k <- 5  # Set this based on your analysis
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
  scale_color_viridis_d()

ggsave("pca_cluster_plot.png", pca_cluster_plot)