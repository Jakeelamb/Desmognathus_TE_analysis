#' PCA Analysis for Desmognathus TE Composition
#'
#' Performs Principal Component Analysis on TE order and superfamily
#' composition data with optional CLR transformation.
#'
#' @usage Rscript scripts/processing/pca.R

# Load shared utilities - find script directory robustly
script_dir <- tryCatch({
    args <- commandArgs(trailingOnly = FALSE)
    file_arg <- grep("^--file=", args, value = TRUE)
    if (length(file_arg) > 0) {
        dirname(normalizePath(sub("^--file=", "", file_arg[1])))
    } else {
        "scripts/processing"
    }
}, error = function(e) "scripts/processing")

source(file.path(script_dir, "pca_utils.R"))

# --- Configuration ---
config <- load_config()
data_dir <- config$results$data
output_dir <- config$results$figures
variance_threshold <- 0.80
k_range <- 2:10

# Ensure output directory exists
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# --- Load Data ---
message("Loading TE composition data...")
order_df <- load_te_data(config, "order")
superfamily_df <- load_te_data(config, "superfamily")

# --- Function: Perform PCA and Plot ---
perform_pca_and_plot <- function(df, title_prefix, id_col = names(df)[1], apply_clr = TRUE) {
    analysis_type_suffix <- if (apply_clr) "CLR" else "Raw"
    full_title_prefix <- paste0(title_prefix, "_", analysis_type_suffix)

    message(paste("\n--- Performing PCA for:", full_title_prefix, "---"))

    # Prepare data using shared utility
    prepared <- prepare_pca_data(df, id_col, apply_clr, verbose = TRUE)

    if (ncol(prepared$data) == 0 || nrow(prepared$data) == 0) {
        message(paste("No data remaining for", full_title_prefix, "PCA."))
        return(list(principal_components = NULL, n_components = 0))
    }

    # Check for zero variance columns
    variances <- apply(prepared$data, 2, var, na.rm = TRUE)
    zero_var_cols <- names(variances)[variances == 0]
    if (length(zero_var_cols) > 0) {
        warning(paste("Removing zero variance columns:", paste(zero_var_cols, collapse = ", ")))
        prepared$data <- prepared$data[, variances > 0, drop = FALSE]
    }

    if (ncol(prepared$data) == 0) {
        message(paste("No columns with variance remain for", full_title_prefix, "PCA."))
        return(list(principal_components = NULL, n_components = 0))
    }

    # Perform PCA
    pca_result <- prcomp(as.matrix(prepared$data), center = TRUE, scale. = FALSE)

    # Calculate explained variance
    explained_variance <- pca_result$sdev^2 / sum(pca_result$sdev^2)
    cumulative_variance <- cumsum(explained_variance)

    # Determine number of components for threshold
    n_components_thresh <- get_n_components(pca_result, variance_threshold)
    message(paste(full_title_prefix, ": Components for >=", variance_threshold * 100, "% variance:", n_components_thresh))

    # Display loadings
    message(paste("\n---", full_title_prefix, ": Top", min(n_components_thresh, 3), "Component Loadings ---"))
    if (n_components_thresh > 0) {
        print(round(pca_result$rotation[, 1:min(n_components_thresh, 3)], 3))
    }

    # Identify top contributors
    top_contributors_pc1 <- character(0)
    top_contributors_pc2 <- character(0)
    if (ncol(pca_result$rotation) >= 1) {
        loadings_pc1 <- abs(pca_result$rotation[, 1])
        top_contributors_pc1 <- names(sort(loadings_pc1, decreasing = TRUE)[1:min(3, length(loadings_pc1))])
        message(paste("Top 3 contributors to PC1:", paste(top_contributors_pc1, collapse = ", ")))
    }
    if (ncol(pca_result$rotation) >= 2) {
        loadings_pc2 <- abs(pca_result$rotation[, 2])
        top_contributors_pc2 <- names(sort(loadings_pc2, decreasing = TRUE)[1:min(3, length(loadings_pc2))])
        message(paste("Top 3 contributors to PC2:", paste(top_contributors_pc2, collapse = ", ")))
    }

    # Generate Scree Plot
    create_scree_plot(
        pca_result,
        title = paste(full_title_prefix, "- Scree Plot"),
        output_path = file.path(output_dir, paste0(title_prefix, "_scree_plot.png"))
    )

    # Extract principal components
    if (n_components_thresh > 0) {
        principal_components <- as.data.frame(pca_result$x[, 1:n_components_thresh, drop = FALSE])
        colnames(principal_components) <- paste0("PC", 1:n_components_thresh)
    } else {
        principal_components <- data.frame(matrix(ncol = 0, nrow = nrow(prepared$data)))
    }

    # Add identifier back
    if (nrow(principal_components) == length(prepared$ids)) {
        principal_components[[id_col]] <- prepared$ids
    }

    # Generate PCA Scatter Plot
    if (n_components_thresh >= 2 && id_col %in% names(principal_components)) {
        xlab_text <- sprintf("PC1 (%.2f%%)\nTop: %s",
                            explained_variance[1] * 100,
                            paste(top_contributors_pc1, collapse = ", "))
        ylab_text <- sprintf("PC2 (%.2f%%)\nTop: %s",
                            explained_variance[2] * 100,
                            paste(top_contributors_pc2, collapse = ", "))

        pca_scatter_plot <- ggplot(principal_components,
                                   aes(x = PC1, y = PC2, label = !!sym(id_col))) +
            geom_point(alpha = 0.8, color = "steelblue", size = 3) +
            ggrepel::geom_text_repel(size = 3, max.overlaps = 15) +
            labs(
                title = paste(full_title_prefix, "- PCA Scatter Plot"),
                x = xlab_text,
                y = ylab_text
            ) +
            theme_minimal() +
            coord_fixed()

        ggsave(file.path(output_dir, paste0(full_title_prefix, "_pca_scatter_plot.png")),
               plot = pca_scatter_plot, width = 10, height = 8, dpi = 300)
    }

    # Return results
    list(
        principal_components = principal_components,
        numeric_pcs = principal_components %>% select(starts_with("PC")),
        filtered_ids = prepared$ids,
        n_components = n_components_thresh
    )
}

# --- Main Execution ---
message("\n========== Order Diversity PCA ==========")
order_pca_raw <- perform_pca_and_plot(order_df, "Order_Diversity", apply_clr = FALSE)
order_pca_clr <- perform_pca_and_plot(order_df, "Order_Diversity", apply_clr = TRUE)

message("\n========== Superfamily Diversity PCA ==========")
superfamily_pca_raw <- perform_pca_and_plot(superfamily_df, "Superfamily_Diversity", apply_clr = FALSE)
superfamily_pca_clr <- perform_pca_and_plot(superfamily_df, "Superfamily_Diversity", apply_clr = TRUE)

message(paste("\nPCA analysis complete. Plots saved in", output_dir))
