library(ggrepel)
library(tidyverse)
library(ggplot2)
library(Rtsne)
library(umap)
library(plotly)

# Load the CSV data
df <- read.csv("/home/jake/Projects/Filtered_dnaPipeTE_Percent_of_Seq_with_variance")

# Filter out specified species
species_to_remove <- c("aeneus", "fuscus", "lycos", "ocoee")
df_filtered <- df[!df$Species %in% species_to_remove, ]

# Identify numeric columns
numeric_cols <- sapply(df_filtered, is.numeric)
df_numeric <- df_filtered[, numeric_cols]

# Function to check if a column has zero variance
has_variance <- function(x) var(x, na.rm = TRUE) > 0

# Remove constant/zero variance columns
df_with_variance <- df_numeric[, sapply(df_numeric, has_variance)]
df_with_variance <- df_with_variance[, -1]

# Print information about the dataset
print(paste("Number of samples:", nrow(df_with_variance)))
print(paste("Number of features:", ncol(df_with_variance)))

# Function to perform PCA
perform_pca <- function(data) {
  pca_result <- prcomp(data, scale = TRUE)  # Center and scale data for PCA
  return(pca_result)
}

# Perform PCA
pca_result <- perform_pca(df_with_variance)

# Create scree plot
scree_data <- data.frame(PC = 1:length(pca_result$sdev), Eigenvalue = pca_result$sdev^2)
scree_plot <- ggplot(scree_data, aes(x = PC, y = Eigenvalue)) +
  geom_line() +
  geom_point() +
  labs(title = "Scree Plot", x = "Principal Component", y = "Eigenvalue")

# Save the scree plot
ggsave("scree_plot_filtered.png", scree_plot)
print("Scree plot saved as 'scree_plot_filtered.png'")

# Get the first 4 principal components
pc_data <- as.data.frame(pca_result$x[, 1:4])
pc_data$Species <- df_filtered$Species  # Use the filtered dataset for species labels

# Create PCA plot with first 4 components and no legend
pca_plot <- ggplot(pc_data, aes(x = PC1, y = PC2, color = Species)) +
  geom_point() +
  geom_text_repel(aes(label = Species), size = 3, nudge_x = 0.2, nudge_y = 0.2) +
  labs(title = "PCA Plot (First 4 Components) - Filtered Species", x = "PC1", y = "PC2") +
  theme(legend.position = "none")  # Remove the legend

# Save the PCA plot
ggsave("pca_plot_4components_filtered_no_legend.png", pca_plot)
print("PCA plot saved as 'pca_plot_4components_filtered_no_legend.png'")

# Export PCA values
pca_values <- pca_result$x[, 1:4]
write.csv(pca_values, "pca_values_filtered.csv", row.names = FALSE)
print("PCA values exported to 'pca_values_filtered.csv'")

# Calculate and print variance explained by each principal component
variance_explained <- pca_result$sdev^2 / sum(pca_result$sdev^2)
print("Variance explained by each principal component:")
print(variance_explained[1:4])

# Calculate feature importance
loadings <- pca_result$rotation[, 1:4]
feature_importance <- data.frame(
  Feature = rownames(loadings),
  Importance = rowSums(loadings^2)
)
feature_importance <- feature_importance[order(feature_importance$Importance, decreasing = TRUE), ]

# Print top 10 important features
print("Top 10 important features:")
print(head(feature_importance, 10))

# Plot feature importance
importance_plot <- ggplot(head(feature_importance, 20), aes(x = reorder(Feature, Importance), y = Importance)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(title = "Top 20 Features by Importance", x = "Feature", y = "Importance") +
  theme_minimal()

# Save the feature importance plot
ggsave("feature_importance_plot.png", importance_plot, width = 10, height = 8)
print("Feature importance plot saved as 'feature_importance_plot.png'")

# Report features associated with each PC
print("Top 5 features associated with each PC:")
for (i in 1:4) {
  pc_loadings <- loadings[, i]
  top_features <- names(sort(abs(pc_loadings), decreasing = TRUE)[1:5])
  print(paste("PC", i, ":", paste(top_features, collapse = ", ")))
}

