#' Phylogenetic PCA Analysis for Desmognathus TE Composition
#'
#' Performs phylogenetically-corrected PCA on TE composition data
#' with clustering analysis and phylomorphospace visualization.
#'
#' @usage Rscript scripts/processing/phylogenetic_pca_analysis.R

# Load shared utilities (includes common libraries)
source(file.path(dirname(sys.frame(1)$ofile %||% "."), "pca_utils.R"))

# Additional libraries for phylogenetic analysis
library(ape)          # For phylogenetic tree manipulation
library(phytools)     # For phylogenetic PCA
library(RColorBrewer) # For cluster colors

# --- Configuration ---
config <- load_config()
data_dir <- config$results$data
output_dir <- config$results$figures
tree_file <- file.path(config$results$data, "desmo900dated_test_cleaned_phylo.tre")
variance_threshold <- 0.80
k_range <- 2:10

# Ensure output directory exists
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# --- Load Data ---
order_df <- load_te_data(config, "order")
superfamily_df <- load_te_data(config, "superfamily")

# --- Data Cleaning: Remove 'D.' prefix from identifiers ---
order_df <- clean_species_ids(order_df)
superfamily_df <- clean_species_ids(superfamily_df)
message("Removed 'D.' prefix from species identifiers")
# ---------------------------------------------------------

# --- Function: Perform PCA and Plot ---
perform_phylogenetic_pca_and_plot <- function(df, phy_tree, title_prefix, id_col = names(df)[1], min_value_threshold = NULL, kmeans_result = NULL, optimal_k = NULL) {
    # Adjust title prefix based on CLR, filtering, and clustering
    filter_suffix <- if (!is.null(min_value_threshold)) paste0("_MinValFiltered", min_value_threshold) else ""
    cluster_suffix <- if (!is.null(kmeans_result) && !is.null(optimal_k) && optimal_k > 0) "_Clustered" else "" # Add _Clustered if valid results passed
    # Base prefix for intermediate steps
    base_title_prefix <- paste0(title_prefix, "_CLR_pPCA", filter_suffix)
    # Final prefix for potentially clustered plot
    final_plot_title_prefix <- paste0(base_title_prefix, cluster_suffix)

    # Message depends on whether this run is for clustering overlay
    if (cluster_suffix == "") {
        message(paste("\n--- Performing Phylogenetic PCA for:", base_title_prefix, "---"))
    } else {
        message(paste("\n--- Generating Clustered Phylomorphospace Plot for:", final_plot_title_prefix, "---"))
    }

    # --- Validate id_col ---
    if (!id_col %in% names(df)) {
        stop(paste("Specified id_col '", id_col, "' not found in the dataframe."))
    }
    # Ensure identifier column is character for matching
    df[[id_col]] <- as.character(df[[id_col]])
    # Ensure unique identifiers in the specified column
    if(any(duplicated(df[[id_col]]))) {
        # Consider how to handle duplicates if they are expected or possible.
        # For now, stopping seems safest for pPCA.
        warning(paste("Duplicate identifiers found in column:", id_col, ". Consider filtering duplicates."))
        # Example: Filter duplicates, keeping the first occurrence
        # df <- df[!duplicated(df[[id_col]]), ]
        # stop(paste("Duplicate identifiers found in column:", id_col, ". Cannot proceed with pPCA."))
    }
    # --- Original ID vector (cleaned) ---
    original_ids <- df[[id_col]]
    # -----------------------------------

    # Select numeric column names
    numeric_cols <- names(df)[sapply(df, is.numeric)]

    # --- Match Tree Tips and Data Identifiers --- #
    message("Matching tree tips with data identifiers...")
    tree_tips <- phy_tree$tip.label
    # Use the actual identifier column from the dataframe for matching
    data_ids_from_col <- df[[id_col]]

    # <<< Debugging: Print first few IDs and tips >>>
    # message("First 5 cleaned data identifiers (from column):")
    # print(head(data_ids_from_col, 5))
    # message("First 5 tree tip labels:")
    # print(head(tree_tips, 5))
    # <<< End Debugging >>>

    common_taxa <- intersect(tree_tips, data_ids_from_col)
    if (length(common_taxa) == 0) {
        stop("No common taxa found between the phylogenetic tree and the dataset.")
    }
    message(paste("Found", length(common_taxa), "common taxa."))

    # --- Subset Data and Tree --- #
    # Subset the tree first
    pruned_tree <- ape::keep.tip(phy_tree, common_taxa)

    # Subset the original dataframe to keep only rows with common taxa
    # Match based on the identifier column
    df_subset <- df[df[[id_col]] %in% common_taxa, , drop = FALSE]

    # Select only numeric columns from the subsetted dataframe
    numeric_df <- df_subset[, numeric_cols, drop = FALSE]

    # Extract the identifier vector from the *subsetted* data
    id_vector_subset <- df_subset[[id_col]]

    # Reorder the numeric data frame rows to match the tree tip order
    # This is crucial for phytools::phyl.pca
    match_order <- match(pruned_tree$tip.label, id_vector_subset)
    numeric_df_ordered <- numeric_df[match_order, , drop = FALSE]
    id_vector_ordered <- id_vector_subset[match_order]

    # Assign final ordered IDs as row names to the numeric matrix later (right before PCA)
    final_ids_for_pca <- id_vector_ordered

    message(paste("Subsetting and ordering data to", nrow(numeric_df_ordered), "rows matching the pruned tree."))
    # ---------------------------- #

    # --- Filtering Steps on the Ordered Numeric Data --- #
    numeric_df <- numeric_df_ordered # Use the ordered data for filtering
    original_colnames <- colnames(numeric_df)
    message("Starting column filtering...")

    # --- Step 1: Filter columns with any NA values ---
    cols_before_na <- colnames(numeric_df)
    na_mask <- colSums(is.na(numeric_df)) == 0
    numeric_df <- numeric_df[, na_mask, drop = FALSE]
    cols_after_na <- colnames(numeric_df)
    removed_na_cols <- setdiff(cols_before_na, cols_after_na)
    if(length(removed_na_cols) > 0) {
        message(paste("Removed columns due to NA values:", paste(removed_na_cols, collapse=", ")))
    } else {
        message("No columns removed due to NA values.")
    }
    # --------------------------------------------------

    # --- Step 1.5: Filter columns based on minimum value threshold ---
    if (!is.null(min_value_threshold) && ncol(numeric_df) > 0) {
        message(paste("Applying minimum value threshold filter:", min_value_threshold))
        cols_before_min_val <- colnames(numeric_df)
        # Calculate minimum for each column, ignoring NA (though NAs should be gone)
        col_mins <- apply(numeric_df, 2, min, na.rm = TRUE)
        # Keep columns where the minimum is greater than or equal to the threshold
        min_val_mask <- col_mins >= min_value_threshold
        numeric_df <- numeric_df[, min_val_mask, drop = FALSE]
        cols_after_min_val <- colnames(numeric_df)
        removed_min_val_cols <- setdiff(cols_before_min_val, cols_after_min_val)
        if(length(removed_min_val_cols) > 0) {
            message(paste("Removed columns due to minimum value <", min_value_threshold, ":", paste(removed_min_val_cols, collapse=", ")))
        } else {
            message("No columns removed based on minimum value threshold.")
        }
    } else if (!is.null(min_value_threshold)) {
         message("Skipping minimum value threshold check as no columns remain after NA filtering.")
    }
    # -----------------------------------------------------------

    # --- Step 2: Filter columns with any Zero values --- (Reinstated)
    if (ncol(numeric_df) > 0) {
        cols_before_zero <- colnames(numeric_df)
        zero_mask <- sapply(numeric_df, function(col) !any(col == 0, na.rm = TRUE))
        numeric_df <- numeric_df[, zero_mask, drop = FALSE]
        cols_after_zero <- colnames(numeric_df)
        removed_zero_cols <- setdiff(cols_before_zero, cols_after_zero)
        if(length(removed_zero_cols) > 0) {
            message(paste("Removed columns due to zero values:", paste(removed_zero_cols, collapse=", ")))
        } else {
            message("No columns removed due to zero values.")
        }
    } else {
        message("Skipping zero column check as no columns remain after NA filtering.")
    }
    # --------------------------------------------------

    # Check if any columns remain after filtering
    if (ncol(numeric_df) == 0) {
        message("No numeric columns remain after filtering NA and zero values.")
        return(list(phylogenetic_principal_components = NULL, numeric_pcs = NULL, filtered_ids = final_ids_for_pca, n_components = 0))
    }
    message(paste("Columns remaining after filtering:", paste(colnames(numeric_df), collapse=", ")))

    # --- Filtering Step: Remove rows with NA in remaining numeric columns ---
    # Note: This requires re-aligning the tree and IDs if rows are dropped here.
    original_row_count <- nrow(numeric_df)
    if (ncol(numeric_df) > 0) {
        rows_before_na <- nrow(numeric_df)
        complete_rows_mask <- complete.cases(numeric_df)

        if (any(!complete_rows_mask)) {
             message(paste("Found rows with NA values in remaining numeric columns after initial filtering."))
             # Filter the numeric data
             numeric_df <- numeric_df[complete_rows_mask, , drop = FALSE]
             # Filter the corresponding IDs
             final_ids_for_pca <- final_ids_for_pca[complete_rows_mask]
             rows_after_na <- nrow(numeric_df)
             message(paste("Removed", rows_before_na - rows_after_na, "rows containing NA values."))

             # Re-prune the tree to match the filtered data
             pruned_tree <- ape::keep.tip(pruned_tree, final_ids_for_pca)
             message("Re-pruned tree to match rows after NA removal.")
        } else {
             message("No rows removed due to NA values in final numeric columns.")
        }
    } else {
        message("Skipping row NA check as no numeric columns remain.")
         # If no numeric columns, pPCA cannot proceed
         message(paste("No numeric columns remain for", base_title_prefix, "pPCA."))
         return(list(phylogenetic_principal_components = NULL, n_components = 0, filtered_ids = character(0)))
    }
    # ----------------------------------------------------------------------- #

    # Final check for sufficient data
    if (nrow(numeric_df) < 2 || ncol(numeric_df) < 1) {
        stop(paste("Insufficient data (rows < 2 or cols < 1) remaining after filtering for", base_title_prefix))
    }

    # Check for zero variance columns
    numeric_mat_pre_pca <- as.matrix(numeric_df)
    variances <- apply(numeric_mat_pre_pca, 2, var, na.rm = TRUE)
    zero_var_cols <- names(variances[variances == 0])
    if (length(zero_var_cols) > 0) {
        warning(paste("Columns with zero variance found and removed for", base_title_prefix, ":", paste(zero_var_cols, collapse=", ")))
        numeric_df <- numeric_df[, !(names(numeric_df) %in% zero_var_cols), drop = FALSE]
         # Re-check if columns remain
        if (ncol(numeric_df) == 0) {
            message(paste("No numeric columns with variance remain for", base_title_prefix, "pPCA."))
            return(list(phylogenetic_principal_components = NULL, n_components = 0, filtered_ids = final_ids_for_pca))
        }
        numeric_mat_pre_pca <- as.matrix(numeric_df) # Update matrix if cols removed
    }

    # Prepare final matrix for PCA
    numeric_mat <- numeric_mat_pre_pca # Use the potentially filtered matrix

    # --- Apply CLR transformation --- # Always apply
    message(paste("--- Applying CLR Transformation before pPCA for", base_title_prefix, "---"))
    # Check for non-positive values before CLR transformation
    if (any(numeric_mat <= 0, na.rm = TRUE)) {
       warning(paste("Non-positive values found in", base_title_prefix, "data. Applying a small pseudo-count for CLR."))
       # Use the robust pseudo-count method from before
       positive_finite_values <- numeric_mat[numeric_mat > 0 & is.finite(numeric_mat)]
       min_positive <- if (length(positive_finite_values) > 0) min(positive_finite_values, na.rm = TRUE) else NA
       pseudo_count <- if (is.finite(min_positive) && min_positive > 0) min_positive / 2 else 1e-9
       numeric_mat[numeric_mat <= 0 | !is.finite(numeric_mat)] <- pseudo_count
       message(paste("Applied pseudo-count for CLR:", format(pseudo_count, scientific = TRUE)))
    }
    # Apply CLR transformation
    message("Applying CLR transformation...")
    comp_data <- compositions::acomp(numeric_mat)
    pca_input_data <- compositions::clr(comp_data) # Returns a standard matrix
    message("Proceeding to pPCA on CLR-transformed data...")
    # -------------------------------- #

    # --- Set row names on the final matrix RIGHT BEFORE PCA --- #
    rownames(pca_input_data) <- final_ids_for_pca
    # ---------------------------------------------------------- #

    # <<< Debugging: Examine matrix before pPCA >>>
    # message("Dimensions of matrix input to PCA/pPCA:")
    # print(dim(pca_input_data))
    # message("Variance of columns input to PCA/pPCA:")
    # print(apply(pca_input_data, 2, var))
    # message("Correlation matrix of CLR data (rounded):")
    # print(round(cor(pca_input_data), 2))
    # <<< End Debugging >>>

    # --- Try standard PCA first --- #
    message("\nAttempting standard PCA on CLR data...")
    std_pca_result <- NULL
    tryCatch({
        std_pca_result <- prcomp(pca_input_data, center = TRUE, scale. = FALSE) # Centering is typical
        message("Standard PCA successful.")
    }, error = function(e_std) {
        message(paste("Standard PCA failed:", e_std$message))
        warning(paste("Standard PCA failed for", base_title_prefix, ". Cannot proceed to pPCA."))
        std_pca_result <<- NULL # Ensure it's NULL on failure
    })
    # ---------------------------- #

    # --- Proceed to pPCA only if standard PCA worked --- #
    if (is.null(std_pca_result)) {
        message("Skipping phylogenetic PCA because standard PCA failed.")
        ppca_result <- NULL
    } else {
        message("\nStandard PCA succeeded, proceeding to Phylogenetic PCA...")
        # Perform Phylogenetic PCA using phytools
        message("Performing phylogenetic PCA...")
        # Use method="BM" as a common choice, could also be "Pagel"
        # Ensure the input data matrix has row names matching tree tips
        ppca_result <- NULL # Initialize result to NULL
        tryCatch({
            # First attempt: mode="cov"
            ppca_result <- phytools::phyl.pca(pruned_tree, pca_input_data, method="BM", mode="cov") # Revert to BM & COV first
            message("Phylogenetic PCA completed (mode='cov', method='BM').")
            # <<< Debugging: Print structure of result >>>
            # message("Structure of ppca_result object:")
            # print(str(ppca_result))
            # <<< End Debugging >>>
        }, error = function(e_cov) {
            message(paste("phyl.pca (mode='cov', method='BM') failed:", e_cov$message))
            warning("phyl.pca (method='BM') failed with mode='cov', trying mode='corr'.")
            tryCatch({
                # Second attempt: mode="corr"
                ppca_result <<- phytools::phyl.pca(pruned_tree, pca_input_data, method="BM", mode="corr") # Revert to BM & Fallback CORR
                message("Phylogenetic PCA completed (mode='corr', method='BM').")
                 # <<< Debugging: Print structure of result >>>
                 # message("Structure of ppca_result object (fallback mode='corr'):")
                 # print(str(ppca_result))
                 # <<< End Debugging >>>
            }, error = function(e_corr) {
                 # Both modes failed
                 message(paste("phyl.pca (fallback mode='corr', method='BM') also failed:", e_corr$message))
                 warning(paste("Phylogenetic PCA (method='BM') failed for", base_title_prefix, ". Returning NULL."))
                 # Set ppca_result to NULL explicitly (already initialized, but good practice)
                 ppca_result <<- NULL
                 # DO NOT stop() here, allow script to continue
            })
        })
    } # End if std_pca_result is not NULL
    # ----------------------------------------------------- #

    # --- Check if pPCA succeeded --- #
    if (is.null(ppca_result)) {
        message(paste("Skipping further processing for", base_title_prefix, "due to pPCA failure."))
        # Return a list indicating failure, consistent with the success structure
        return(list(phylogenetic_principal_components = NULL,
                    numeric_pcs = NULL,
                    filtered_ids = final_ids_for_pca, # Return IDs that went into the attempt
                    n_components = 0))
    }
    # ----------------------------- #

    # Calculate explained variance from eigenvalues (StDev^2)
    # Standard deviations are in ppca_result$StDev
    # Eigenvalues might be directly in ppca_result$Eval

    eigenvalues <- NULL # Initialize

    if (!is.null(ppca_result$StDev)) {
        message("Using ppca_result$StDev to calculate eigenvalues.")
        # <<< Debugging: Print StDev >>>
        # message("Raw Standard Deviations (ppca_result$StDev):")
        # print(ppca_result$StDev)
        # <<< End Debugging >>>
        eigenvalues <- ppca_result$StDev^2
    } else if (!is.null(ppca_result$Eval) && is.numeric(ppca_result$Eval)) {
        message("ppca_result$StDev is NULL. Attempting to use eigenvalues from ppca_result$Eval.")
        # <<< Debugging: Print Eval object >>>
        # message("Contents of ppca_result$Eval:")
        # print(ppca_result$Eval)
        # <<< End Debugging >>>
        # Eval might be a matrix (diagonal are eigenvalues) or a vector
        if (is.matrix(ppca_result$Eval)) {
            eigenvalues <- diag(ppca_result$Eval)
        } else {
            eigenvalues <- ppca_result$Eval # Assume it's already the vector of eigenvalues
        }
        # Ensure eigenvalues are non-negative (can happen with numerical issues)
        eigenvalues[eigenvalues < 0] <- 0
        message("Extracted eigenvalues from Eval.")
    } else {
        warning(paste("Could not find valid standard deviations ($StDev) or eigenvalues ($Eval) in ppca_result for", base_title_prefix))
        # Proceed to the check below, which will likely trigger the return
    }

    # <<< Debugging: Print Eigenvalues >>>
    # message("Calculated/Extracted Eigenvalues:")
    # print(eigenvalues)
    # <<< End Debugging >>>

    # Check if eigenvalues sum is valid for variance calculation
    eigenvalue_sum <- sum(eigenvalues, na.rm = TRUE) # Added na.rm just in case
    if (is.null(eigenvalues) || length(eigenvalues) == 0 || is.na(eigenvalue_sum) || eigenvalue_sum <= 0) {
        warning(paste("Sum of eigenvalues is zero, negative, non-finite, or eigenvalues could not be determined for", base_title_prefix, ". Cannot calculate explained variance. Skipping downstream analysis."))
        return(list(phylogenetic_principal_components = NULL, numeric_pcs = NULL, filtered_ids = final_ids_for_pca, n_components = 0))
    }

    explained_variance <- eigenvalues / eigenvalue_sum
    cumulative_variance <- cumsum(explained_variance)

    # Determine number of components for the threshold
    n_components_thresh <- if(length(cumulative_variance) > 0) which.max(cumulative_variance >= variance_threshold) else 0
    if (n_components_thresh == 0 && length(cumulative_variance) > 0 && cumulative_variance[1] >= variance_threshold) {
        n_components_thresh <- 1 # Handle case where first component meets threshold
    }
    message(paste(base_title_prefix, ": Number of components explaining >=", variance_threshold * 100, "% variance:", n_components_thresh))

    # Display Principal Component Loadings (Rotation Matrix)
    message(paste("\n---", base_title_prefix, ": Phylogenetic Principal Component Loadings (Top", n_components_thresh, "Components) ---"))
    # Loadings are in ppca_result$L (n_features x n_components)
    if (n_components_thresh > 0 && n_components_thresh <= ncol(ppca_result$L)) {
        print(round(ppca_result$L[, 1:n_components_thresh, drop = FALSE], 3))
    } else if (n_components_thresh == 0) {
        message("No components selected based on variance threshold.")
    } else {
         message("Cannot display loadings (invalid n_components_thresh).")
    }
    message("---------------------------------------------------")

    # Identify and print top contributors to PC1 and PC2
    top_contributors_pc1 <- character(0)
    top_contributors_pc2 <- character(0)
    if (ncol(ppca_result$L) >= 1 && n_components_thresh >=1) {
        loadings_pc1 <- abs(ppca_result$L[, 1])
        top_contributors_pc1 <- names(sort(loadings_pc1, decreasing = TRUE)[1:min(3, length(loadings_pc1))])
        message(paste(base_title_prefix, "Top 3 contributors to pPC1 (based on absolute loading):"))
        print(top_contributors_pc1)
    }
    if (ncol(ppca_result$L) >= 2 && n_components_thresh >=2) {
        loadings_pc2 <- abs(ppca_result$L[, 2])
        top_contributors_pc2 <- names(sort(loadings_pc2, decreasing = TRUE)[1:min(3, length(loadings_pc2))])
        message(paste(base_title_prefix,"Top 3 contributors to pPC2 (based on absolute loading):"))
        print(top_contributors_pc2)
    }

    # Scree Plot Data
    scree_data <- data.frame(
        PC = 1:length(explained_variance),
        Variance = explained_variance,
        CumulativeVariance = cumulative_variance
    )

    # Generate Scree Plot (remains largely the same logic)
    scree_plot <- ggplot(scree_data, aes(x = PC)) +
        geom_line(aes(y = Variance, linetype = "Individual"), color = "blue") +
        geom_point(aes(y = Variance), color = "blue") +
        geom_line(aes(y = CumulativeVariance, linetype = "Cumulative"), color = "red") +
        geom_point(aes(y = CumulativeVariance), color = "red") +
        geom_hline(yintercept = variance_threshold, linetype = "dashed", color = "darkgreen") +
        geom_vline(xintercept = n_components_thresh, linetype = "dashed", color = "purple") +
        scale_linetype_manual(name = "Variance Type", values = c("Individual" = "dashed", "Cumulative" = "solid")) +
        labs(
            title = paste(base_title_prefix, "- Scree Plot"),
            x = "Number of Phylogenetic Principal Components",
            y = "Proportion of Variance Explained"
        ) +
        annotate("text", x = n_components_thresh, y = 0, label = paste(n_components_thresh, "Components"), vjust = -0.5, hjust = -0.1, color = "purple") +
        annotate("text", x = max(scree_data$PC), y = variance_threshold, label = paste0(variance_threshold * 100, "% Threshold"), hjust = 1, vjust = -0.5, color = "darkgreen") +
        theme_minimal() +
        theme(legend.position = "bottom")

    ggsave(file.path(output_dir, paste0(base_title_prefix, "_scree_plot.png")), plot = scree_plot, width = 10, height = 6)

    # Extract phylogenetic principal components (scores) up to the threshold
    # Scores are in ppca_result$S (n_samples x n_components)
    if (n_components_thresh > 0 && n_components_thresh <= ncol(ppca_result$S)) {
        phylogenetic_pcs <- as.data.frame(ppca_result$S[, 1:n_components_thresh, drop = FALSE])
        colnames(phylogenetic_pcs) <- paste0("pPC", 1:n_components_thresh)
    } else {
        phylogenetic_pcs <- data.frame(matrix(ncol = 0, nrow = nrow(pca_input_data)))
    }

    # Add identifier back to principal components
    phylogenetic_pcs[[id_col]] <- rownames(phylogenetic_pcs) # Use rownames which are the matched IDs

    # Generate PCA Scatter Plot (pPC1 vs pPC2) if possible
    # Use the final title prefix which might include _Clustered
    plot_filename <- file.path(output_dir, paste0(final_plot_title_prefix, "_phylomorphospace_plot.png"))
    if (n_components_thresh >= 2) {
        message(paste("Generating phylomorphospace plot (pPC1 vs pPC2) for", final_plot_title_prefix))
        # Construct axis labels with variance and contributors
        xlab_text <- sprintf("pPC1 (%.2f%% Variance)", explained_variance[1] * 100)
        if (length(top_contributors_pc1) > 0) {
            xlab_text <- paste0(xlab_text, "\nTop Contributors: ", paste(top_contributors_pc1, collapse=", "))
        }
        ylab_text <- sprintf("pPC2 (%.2f%% Variance)", explained_variance[2] * 100)
        if (length(top_contributors_pc2) > 0) {
            ylab_text <- paste0(ylab_text, "\nTop Contributors: ", paste(top_contributors_pc2, collapse=", "))
        }

        # Use phytools::phylomorphospace
        scores_matrix <- as.matrix(phylogenetic_pcs[, c("pPC1", "pPC2")])
        rownames(scores_matrix) <- final_ids_for_pca

        # --- Clustering visualization settings (if applicable) ---
        cluster_colors <- NULL
        # Start with the base title, add cluster info if applicable
        plot_main_title <- paste(base_title_prefix, "- Phylomorphospace (pPC1 vs pPC2)")
        point_colors <- "black" # Default point color

        if (!is.null(kmeans_result) && !is.null(optimal_k) && optimal_k > 0) {
            # Ensure cluster assignments length matches scores
            if (length(kmeans_result$cluster) == nrow(scores_matrix)) {
                # Generate colors using RColorBrewer
                palette <- brewer.pal(n = max(3, optimal_k), name = "Set1")
                if(optimal_k > length(palette)) palette <- colorRampPalette(palette)(optimal_k)

                point_colors <- palette[kmeans_result$cluster] # Use cluster colors for points
                # Update title to reflect clustering
                plot_main_title <- paste(final_plot_title_prefix, "- Phylomorphospace (pPC1 vs pPC2) with", optimal_k, "Clusters")
            } else {
                warning("Cluster assignments length mismatch for 2D plot. Using default point colors.")
                # Keep base title if colors aren't applied
            }
        } else {
             # If not clustering, ensure the non-clustered title is used
             plot_main_title <- paste(base_title_prefix, "- Phylomorphospace (pPC1 vs pPC2)")
        }
        # -----------------------------------

        # Open PNG device with increased width to accommodate labels
        png(filename = plot_filename, width = 12, height = 8, units = "in", res = 300)

        # Save current graphical parameters and set new margins
        original_par <- par(no.readonly = TRUE)
        par(mar = c(5, 8, 4, 2) + 0.1) # Increase left margin

        # Plot the phylogeny and morphospace structure (including points colored by cluster if applicable)
        # Use the 'col' argument for phytools points if possible, otherwise overlay points
        # phylomorphospace doesn't directly support per-point color vectors easily.
        # Plot structure first, then overlay colored points and labels.

        phytools::phylomorphospace(pruned_tree, scores_matrix,
                                   label = "off",        # Turn off default labels
                                   node.size = c(0, 0.5), # Keep small node visibility
                                   ftype = "off",        # Turn off default points/text at tips
                                   xlab = xlab_text,
                                   ylab = "", # Empty ylab to add manually
                                   main = plot_main_title)

        # Add colored points (use simple points, pch=19 or similar)
        points(scores_matrix[,1], scores_matrix[,2], pch = 19, col = point_colors, cex = 1.2)

        # Add labels using text with an offset (always add labels)
        text(scores_matrix[,1], scores_matrix[,2],
             labels = rownames(scores_matrix),
             pos = 4,
             cex = 0.6,
             offset = 0.5)

        # Add y-axis label with adjusted position
        mtext(ylab_text, side = 2, line = 5.5, cex = 0.9) # Increased line

        # Restore original graphical parameters
        par(original_par)

        message(paste("Saved phylomorphospace plot to:", plot_filename))

    } else if (n_components_thresh == 1) {
        message(paste("Generating phylomorphospace plot (pPC1 only) for", final_plot_title_prefix))
        # Construct axis label
        xlab_text <- sprintf("pPC1 (%.2f%% Variance)", explained_variance[1] * 100)
        if (length(top_contributors_pc1) > 0) {
            xlab_text <- paste0(xlab_text, "\nTop Contributors: ", paste(top_contributors_pc1, collapse=", "))
        }

        # Use phytools::phylomorphospace1D
        scores_vector <- phylogenetic_pcs[,"pPC1"] # Extract PC1
        names(scores_vector) <- final_ids_for_pca # Ensure names match

        # --- Clustering visualization (1D) ---
        cluster_colors <- NULL
        # Start with base title
        plot_main_title <- paste(base_title_prefix, "- Phylomorphospace (pPC1)")

        if (!is.null(kmeans_result) && !is.null(optimal_k) && optimal_k > 0) {
            # Ensure cluster assignments length matches scores
            if (length(kmeans_result$cluster) == length(scores_vector)) {
                palette <- brewer.pal(n = max(3, optimal_k), name = "Set1") # Use Set1 palette
                if(optimal_k > length(palette)) palette <- colorRampPalette(palette)(optimal_k)
                cluster_colors <- palette[kmeans_result$cluster]
                # Add cluster info to title
                plot_main_title <- paste(final_plot_title_prefix, "- Phylomorphospace (pPC1) with", optimal_k, "Clusters")
            } else {
                warning("Cluster assignments length mismatch for 1D plot. Skipping coloring.")
                # Keep base title
            }
        } else {
             # Ensure non-clustered title
             plot_main_title <- paste(base_title_prefix, "- Phylomorphospace (pPC1)")
        }
        # -----------------------------------

        # Open PNG device with increased width
        png(filename = plot_filename, width = 12, height = 6, units = "in", res = 300)

        # Save current graphical parameters and set new margins
        original_par <- par(no.readonly = TRUE)
        par(mar = c(5, 8, 4, 2) + 0.1) # Increase left margin from 6 to 8

        # Plot 1D phylomorphospace (basic structure)
        # Plotting points is handled separately if clustering
        phytools::phylomorphospace1d(pruned_tree, scores_vector,
                                     label = "off",
                                     node.size = c(0, 0), # Set node size to 0
                                     ftype="off", # Turn off default points/labels on tips
                                     xlab = xlab_text,
                                     ylab = "", # Empty ylab
                                     main = plot_main_title)

        # Add colored points if clustering was successful
        if (!is.null(cluster_colors)) {
             # Add points along the axis. The y-coords from phylomorphospace1d are not easily accessible.
             # We can approximate by plotting points at a fixed y offset or try using plot.phylo with points.
             # For simplicity, we will overlay points directly on the x-axis (scores_vector) at y=1 (arbitrary low value)
             # This won't align with the phylogeny visually but shows the cluster colors.
             points(scores_vector, rep(par("usr")[3] * 0.95, length(scores_vector)), # Plot near bottom edge
                    pch = 19, col = cluster_colors, cex = 1.2)
             message("Colored points for 1D plot added near bottom axis due to plotting limitations.")
        } else {
             # Add default points if not clustering
             points(scores_vector, rep(par("usr")[3] * 0.95, length(scores_vector)),
                     pch = 19, col = "black", cex = 1.2)
        }

        # Add y-axis label (e.g., "Phylogeny")
        mtext("Phylogeny", side = 2, line = 5.5, cex = 0.9) # Increased line from 4.5 to 5.5

        # Add sample labels (always, regardless of clustering)
        # We need the coordinates from the plot object... difficult.
        # Let's skip labels on the 1D plot for now to avoid complexity.
        # text(scores_vector, ???, labels=names(scores_vector), pos=4, cex=0.6)

        # Restore original graphical parameters
        par(original_par)

        # Close PNG device
        dev.off()

        message(paste("Saved 1D phylomorphospace plot to:", plot_filename))
        message(paste(base_title_prefix, ": Only 1 phylogenetic principal component explains >=", variance_threshold * 100, "% variance. Plot shows pPC1 vs phylogeny."))
    } else {
        message(paste(base_title_prefix, ": Not enough phylogenetic principal components (", n_components_thresh, ") for a phylomorphospace plot."))
    }

    # --- Prepare return value --- #
    # Separate numeric pPCs from the ID column for clarity downstream
    numeric_pcs_df <- phylogenetic_pcs %>% select(starts_with("pPC"))

    return(list(phylogenetic_principal_components = phylogenetic_pcs, # Includes ID col
                numeric_pcs = numeric_pcs_df,                      # Numeric pPCs only
                filtered_ids = final_ids_for_pca,                  # IDs used in the final pPCA
                n_components = n_components_thresh))
}

# --- Function: Perform Clustering Analysis ---
perform_clustering_analysis <- function(principal_components, title_prefix) {
    message(paste("
--- Performing Clustering Analysis for:", title_prefix, "---"))

    # Check if there are enough data points and components for clustering
    if (nrow(principal_components) < 2) {
        message("Not enough data points for clustering.")
        return(NULL)
    }
     if (ncol(principal_components) == 0) {
        message("No principal components available for clustering.")
        return(NULL)
    }

    # Determine a reasonable k_max based on data size
    max_possible_k <- nrow(principal_components) - 1
    current_k_max <- max(k_range)
    effective_k_max <- min(current_k_max, max_possible_k)
    if (effective_k_max < min(k_range)) {
        message(paste(title_prefix, ": Not enough data points (", nrow(principal_components), ") to test the minimum k (", min(k_range), "). Skipping clustering."))
        return(list(optimal_k = NA, kmeans_result = NULL))
    }
    effective_k_range <- min(k_range):effective_k_max
    message(paste(title_prefix, ": Testing k values from", min(effective_k_range), "to", max(effective_k_range)))


    # Calculate WCSS (Within-Cluster Sum of Squares) using Elbow Method
    # Use factoextra's fviz_nbclust for convenience, or calculate manually
    message(paste(title_prefix, ": Generating Elbow plot..."))
    elbow_plot_obj <- fviz_nbclust(principal_components, kmeans, method = "wss", k.max = effective_k_max) +
                        ggtitle(paste(title_prefix, "- Elbow Method For Optimal k")) +
                        labs(subtitle = "Within-Cluster Sum of Squares")
    ggsave(file.path(output_dir, paste0(title_prefix, "_elbow_plot.png")), plot = elbow_plot_obj, width = 10, height = 6)

    # Calculate Silhouette Scores
    # Use factoextra's fviz_nbclust or calculate manually
    message(paste(title_prefix, ": Generating Silhouette plot..."))
    # Ensure k.max > 1 for silhouette
    if (effective_k_max > 1) {
        silhouette_plot_obj <- fviz_nbclust(principal_components, kmeans, method = "silhouette", k.max = effective_k_max) +
                                 ggtitle(paste(title_prefix, "- Silhouette Score For Optimal k"))
        ggsave(file.path(output_dir, paste0(title_prefix, "_silhouette_plot.png")), plot = silhouette_plot_obj, width = 10, height = 6)
    } else {
        message(paste(title_prefix, ": Skipping Silhouette plot (k.max <= 1)."))
    }

    # Extract silhouette scores if needed for determining optimal k programmatically
    silhouette_scores <- sapply(effective_k_range, function(k) {
        # Check added: Ensure k > 1 for silhouette calculation and enough data points
        if (k > 1 && nrow(principal_components) >= k) {
             # Added error handling for kmeans/silhouette
             tryCatch({
                 km <- kmeans(principal_components, centers = k, nstart = 25) # nstart improves stability
                 # Check if silhouette can be calculated (requires >1 unique cluster distance)
                 sil_result <- try(silhouette(km$cluster, dist(principal_components)), silent = TRUE)
                 if (inherits(sil_result, "try-error")) {
                     warning(paste(title_prefix, ": Silhouette calculation failed for k=", k, ". Returning NA."))
                     NA
                 } else {
                     mean(sil_result[, 3]) # Return mean silhouette width
                 }
             }, error = function(e) {
                 warning(paste(title_prefix, ": K-means failed for k=", k, ". Error:", e$message, ". Returning NA."))
                 NA
             })
        } else {
            NA # Not enough samples for this k or k=1
        }
    })

    silhouette_data <- data.frame(k = effective_k_range, score = silhouette_scores) %>% drop_na()

    if(nrow(silhouette_data) > 0) {
        optimal_k_silhouette <- silhouette_data$k[which.max(silhouette_data$score)]
        message(paste(title_prefix, ": Optimal number of clusters based on Silhouette Score:", optimal_k_silhouette))
    } else {
         message(paste(title_prefix, ": Could not determine optimal k via Silhouette Score (possibly insufficient data points relative to k or calculation errors). Setting optimal k to NA."))
         optimal_k_silhouette <- NA # Explicitly set to NA
    }

    # --- Perform final clustering with optimal k and return result ---
    final_kmeans_result <- NULL
    if (!is.na(optimal_k_silhouette) && optimal_k_silhouette > 0 && nrow(principal_components) >= optimal_k_silhouette) {
        message(paste(title_prefix, ": Performing final k-means clustering with k=", optimal_k_silhouette))
        set.seed(123) # for reproducibility
        tryCatch({
             final_kmeans_result <- kmeans(principal_components, centers = optimal_k_silhouette, nstart = 25)
             message(paste(title_prefix, ": Final k-means clustering successful."))
        }, error = function(e) {
             warning(paste(title_prefix, ": Final k-means clustering failed for k=", optimal_k_silhouette, ". Error:", e$message))
             final_kmeans_result <<- NULL # Ensure it's NULL on error
        })
    } else {
        message(paste(title_prefix, ": Cannot perform final clustering (optimal k not determined, is NA/0, or insufficient data). Optimal k:", optimal_k_silhouette))
        optimal_k_silhouette <- NA # Ensure k is NA if clustering isn't performed
    }
    # ------------------------------------------------------------------

    # Note: Elbow method often requires visual inspection or algorithms like Kneedle.
    # fviz_nbclust provides a visual guide.

    return(list(optimal_k = optimal_k_silhouette, kmeans_result = final_kmeans_result)) # Return list
}

# --- Function: Plot K-means Clusters --- 
plot_kmeans_clusters <- function(numeric_pcs, kmeans_result, optimal_k, filtered_ids, title_prefix) {
    message(paste("\n--- Generating K-means Cluster Plot for:", title_prefix, "(Optimal K=", optimal_k, ") ---"))

    # Ensure we have necessary components
    if (is.null(numeric_pcs) || nrow(numeric_pcs) < 2 || ncol(numeric_pcs) < 2) {
        message(paste(title_prefix, ": Insufficient PCA data (need >= 2 components and >= 2 samples) for cluster plot."))
        return()
    }
    if (is.null(kmeans_result) || is.null(optimal_k) || optimal_k <= 0) {
        message(paste(title_prefix, ": K-means clustering results not available or optimal K is invalid. Skipping cluster plot."))
        return()
    }
    if (length(filtered_ids) != nrow(numeric_pcs)) {
         message(paste(title_prefix, ": ID vector length mismatch. Labels might be incorrect on cluster plot."))
         # Proceeding, but labels might be wrong.
    }

    # Prepare data for plotting (use only PC1 and PC2 for fviz_cluster)
    plot_data <- numeric_pcs[, 1:2, drop = FALSE]
    rownames(plot_data) <- filtered_ids # Assign IDs as rownames for fviz_cluster labeling

    # Use Convex Hull
    plot_filename <- file.path(output_dir, paste0(title_prefix, "_kmeans_cluster_plot.png"))
    cluster_plot <- fviz_cluster(kmeans_result, # Use the passed optimal result
                                 data = plot_data, # Use PC1 & PC2
                                 geom = c("point", "text"), # Show points and labels
                                 repel = TRUE,           # Use repel for labels
                                 ellipse.type = "convex",  # Use convex hull
                                 ggtheme = theme_minimal(),
                                 main = paste(title_prefix, "- K-means Clusters (pPC1 vs pPC2, Optimal K=", optimal_k, ")"))

    ggsave(plot_filename, plot = cluster_plot, width = 10, height = 8)
    message(paste("Saved K-means cluster plot to:", plot_filename))
}

# --- Function: Plot Phylogeny with Colored Tips based on Clusters ---
plot_clustered_phylogeny <- function(full_phy_tree, kmeans_result, optimal_k, filtered_ids, title_prefix) {
    plot_title <- paste(title_prefix, "- Phylogeny with Optimal Clusters (K=", optimal_k, ")")
    message(paste("\n--- Generating Phylogeny Plot with Cluster-Colored Tips for:", plot_title, "---"))

    # --- Input Validation ---
    if (is.null(full_phy_tree)) {
        message("Input tree is NULL. Skipping phylogeny plot.")
        return()
    }
    if (is.null(kmeans_result) || is.null(optimal_k) || optimal_k <= 0) {
        message(paste(title_prefix, ": K-means clustering results not available or optimal K is invalid. Skipping phylogeny plot."))
        return()
    }
    if (is.null(filtered_ids) || length(filtered_ids) == 0) {
        message(paste(title_prefix, ": Filtered ID vector is empty. Skipping phylogeny plot."))
        return()
    }
    if (length(kmeans_result$cluster) != length(filtered_ids)) {
        message(paste(title_prefix, ": K-means cluster vector length does not match filtered ID vector length. Skipping phylogeny plot."))
        return()
    }
    # ----------------------

    # Prune tree to match the tips included in the clustering analysis
    common_tips <- intersect(full_phy_tree$tip.label, filtered_ids)
    if (length(common_tips) < 2) {
        message(paste(title_prefix, ": Fewer than 2 common tips between tree and clustered data. Skipping phylogeny plot."))
        return()
    }
    pruned_tree <- ape::keep.tip(full_phy_tree, common_tips)
    message(paste("Pruned tree to", length(pruned_tree$tip.label), "tips for plotting."))

    # Create a mapping from filtered_ids to clusters
    id_to_cluster_map <- setNames(kmeans_result$cluster, filtered_ids)

    # Generate color palette
    palette <- brewer.pal(n = max(3, optimal_k), name = "Set1")
    if(optimal_k > length(palette)) palette <- colorRampPalette(palette)(optimal_k)

    # Create a vector of tip colors in the order of the pruned tree's tip labels
    tip_colors <- sapply(pruned_tree$tip.label, function(tip) {
        cluster_num <- id_to_cluster_map[[tip]]
        if (!is.null(cluster_num) && cluster_num > 0 && cluster_num <= length(palette)) {
            palette[cluster_num]
        } else {
            "black" # Default color if mapping fails (shouldn't happen)
        }
    })

    # Set up plot filename
    plot_filename <- file.path(output_dir, paste0(title_prefix, "_phylogeny_optimal_clusters.png"))

    # Plotting parameters
    max_age <- max(nodeHeights(pruned_tree)[,2]) # Use pruned tree height
    age_breaks <- seq(0, ceiling(max_age), by = 5)

    # Generate plot
    png(filename = plot_filename, width = 12, height = 10, units = "in", res = 300)
    par(mar = c(5, 4, 4, 2) + 0.1) # Adjust margins slightly
    plot(pruned_tree, type = "phylogram", direction = "right",
         show.tip.label = TRUE,
         tip.color = tip_colors, # Apply colors here
         cex = 0.7,              # Adjust label size if needed
         label.offset = 0.5,
         edge.width = 1.5,
         main = plot_title)
    axis(1, at = max_age - age_breaks, labels = age_breaks, las = 1) # Adjust axis for rightwards plot
    mtext("Time (Million Years Ago)", side = 1, line = 3)
    abline(v = max_age - age_breaks, lty = 2, col = "gray90") # Adjust grid lines
    dev.off()

    message(paste("Saved phylogeny plot with cluster-colored tips to:", plot_filename))
}

# --- Main Execution ---

# Load the phylogenetic tree
message(paste("Loading phylogenetic tree from:", tree_file))
if (!file.exists(tree_file)) {
    stop("Phylogenetic tree file not found at the specified path: ", tree_file)
}
phy_tree <- ape::read.tree(tree_file)
message("Phylogenetic tree loaded successfully.")

# Perform Phylogenetic PCA for order_df
# Raw pPCA # Removed
# order_pca_result_raw <- perform_phylogenetic_pca_and_plot(order_df, phy_tree, "Order_Diversity", apply_clr = FALSE)
# CLR pPCA
order_pca_result_clr <- perform_phylogenetic_pca_and_plot(order_df, phy_tree, "Order_Diversity")

# CLR pPCA with Min Value Filter (0.01)
order_pca_result_clr_filtered <- perform_phylogenetic_pca_and_plot(
    order_df,
    phy_tree,
    "Order_Diversity", # Base title prefix
    min_value_threshold = 0.01
)

# Perform Phylogenetic PCA for superfamily_df
# Raw pPCA # Removed
# superfamily_pca_result_raw <- perform_phylogenetic_pca_and_plot(superfamily_df, phy_tree, "Superfamily_Diversity", apply_clr = FALSE)
# CLR pPCA
superfamily_pca_result_clr <- perform_phylogenetic_pca_and_plot(superfamily_df, phy_tree, "Superfamily_Diversity")

# CLR pPCA with Min Value Filter (0.01)
superfamily_pca_result_clr_filtered <- perform_phylogenetic_pca_and_plot(
    superfamily_df,
    phy_tree,
    "Superfamily_Diversity", # Base title prefix
    min_value_threshold = 0.01
)

# --- Clustering Analysis --- 
# Run clustering on the numeric principal components from each PCA result.

message("\n--- Performing Clustering Analysis on CLR-Transformed pPCA Results ---")

# Initialize clustering result variables
clustering_result_order_clr <- NULL
clustering_result_order_clr_filtered <- NULL
clustering_result_superfamily_clr <- NULL
clustering_result_superfamily_clr_filtered <- NULL

# Clustering for Order Diversity (CLR)
if (!is.null(order_pca_result_clr) && !is.null(order_pca_result_clr$numeric_pcs) && nrow(order_pca_result_clr$numeric_pcs) > 1 && ncol(order_pca_result_clr$numeric_pcs) > 0) {
    clustering_result_order_clr <- perform_clustering_analysis(
        order_pca_result_clr$numeric_pcs,
        "Order_Diversity_CLR_pPCA" # Match PCA title prefix logic
    )
} else {
    message("Skipping clustering for Order Diversity (CLR pPCA) due to insufficient PCA results.")
}

# Clustering for Order Diversity (CLR, Filtered)
if (!is.null(order_pca_result_clr_filtered) && !is.null(order_pca_result_clr_filtered$numeric_pcs) && nrow(order_pca_result_clr_filtered$numeric_pcs) > 1 && ncol(order_pca_result_clr_filtered$numeric_pcs) > 0) {
    clustering_result_order_clr_filtered <- perform_clustering_analysis(
        order_pca_result_clr_filtered$numeric_pcs,
        "Order_Diversity_CLR_pPCA_MinValFiltered0.01" # Match PCA title prefix logic
    )
} else {
    message("Skipping clustering for Order Diversity (CLR pPCA, Filtered) due to insufficient PCA results.")
}

# Clustering for Superfamily Diversity (CLR)
if (!is.null(superfamily_pca_result_clr) && !is.null(superfamily_pca_result_clr$numeric_pcs) && nrow(superfamily_pca_result_clr$numeric_pcs) > 1 && ncol(superfamily_pca_result_clr$numeric_pcs) > 0) {
    clustering_result_superfamily_clr <- perform_clustering_analysis(
        superfamily_pca_result_clr$numeric_pcs,
        "Superfamily_Diversity_CLR_pPCA" # Match PCA title prefix logic
    )
} else {
    message("Skipping clustering for Superfamily Diversity (CLR pPCA) due to insufficient PCA results.")
}

# Clustering for Superfamily Diversity (CLR, Filtered)
if (!is.null(superfamily_pca_result_clr_filtered) && !is.null(superfamily_pca_result_clr_filtered$numeric_pcs) && nrow(superfamily_pca_result_clr_filtered$numeric_pcs) > 1 && ncol(superfamily_pca_result_clr_filtered$numeric_pcs) > 0) {
    clustering_result_superfamily_clr_filtered <- perform_clustering_analysis(
        superfamily_pca_result_clr_filtered$numeric_pcs,
        "Superfamily_Diversity_CLR_pPCA_MinValFiltered0.01" # Match PCA title prefix logic
    )
} else {
    message("Skipping clustering for Superfamily Diversity (CLR pPCA, Filtered) due to insufficient PCA results.")
}

# --- Generate K-means Cluster Plots (Separate from Phylogeny) ---

message("\n--- Generating K-means Cluster Plots (PC1 vs PC2) ---")

# Plot Order Diversity (CLR) Clusters
if (!is.null(clustering_result_order_clr) && !is.null(order_pca_result_clr)) {
    plot_kmeans_clusters(
        numeric_pcs = order_pca_result_clr$numeric_pcs,
        kmeans_result = clustering_result_order_clr$kmeans_result,
        optimal_k = clustering_result_order_clr$optimal_k,
        filtered_ids = order_pca_result_clr$filtered_ids,
        title_prefix = "Order_Diversity_CLR_pPCA"
    )
} else {
    message("Skipping cluster plot for Order Diversity (CLR pPCA) - Clustering or PCA failed.")
}

# Plot Order Diversity (CLR, Filtered) Clusters
if (!is.null(clustering_result_order_clr_filtered) && !is.null(order_pca_result_clr_filtered)) {
    plot_kmeans_clusters(
        numeric_pcs = order_pca_result_clr_filtered$numeric_pcs,
        kmeans_result = clustering_result_order_clr_filtered$kmeans_result,
        optimal_k = clustering_result_order_clr_filtered$optimal_k,
        filtered_ids = order_pca_result_clr_filtered$filtered_ids,
        title_prefix = "Order_Diversity_CLR_pPCA_MinValFiltered0.01"
    )
} else {
    message("Skipping cluster plot for Order Diversity (CLR pPCA, Filtered) - Clustering or PCA failed.")
}

# Plot Superfamily Diversity (CLR) Clusters
if (!is.null(clustering_result_superfamily_clr) && !is.null(superfamily_pca_result_clr)) {
    plot_kmeans_clusters(
        numeric_pcs = superfamily_pca_result_clr$numeric_pcs,
        kmeans_result = clustering_result_superfamily_clr$kmeans_result,
        optimal_k = clustering_result_superfamily_clr$optimal_k,
        filtered_ids = superfamily_pca_result_clr$filtered_ids,
        title_prefix = "Superfamily_Diversity_CLR_pPCA"
    )
} else {
    message("Skipping cluster plot for Superfamily Diversity (CLR pPCA) - Clustering or PCA failed.")
}

# Plot Superfamily Diversity (CLR, Filtered) Clusters
if (!is.null(clustering_result_superfamily_clr_filtered) && !is.null(superfamily_pca_result_clr_filtered)) {
    plot_kmeans_clusters(
        numeric_pcs = superfamily_pca_result_clr_filtered$numeric_pcs,
        kmeans_result = clustering_result_superfamily_clr_filtered$kmeans_result,
        optimal_k = clustering_result_superfamily_clr_filtered$optimal_k,
        filtered_ids = superfamily_pca_result_clr_filtered$filtered_ids,
        title_prefix = "Superfamily_Diversity_CLR_pPCA_MinValFiltered0.01"
    )
} else {
    message("Skipping cluster plot for Superfamily Diversity (CLR pPCA, Filtered) - Clustering or PCA failed.")
}

# --- Generate Phylomorphospace Plots with Optimal Cluster Coloring ---

message("\n--- Generating Phylomorphospace Plots with Optimal Cluster Coloring ---")

# Re-run PCA plot generation, passing optimal clustering results for coloring

# Order Diversity (CLR) - Clustered Phylomorphospace
if (!is.null(clustering_result_order_clr) && !is.null(order_pca_result_clr)) {
    perform_phylogenetic_pca_and_plot(
        order_df,
        phy_tree,
        "Order_Diversity",
        min_value_threshold = NULL, # Match original PCA settings
        kmeans_result = clustering_result_order_clr$kmeans_result,
        optimal_k = clustering_result_order_clr$optimal_k
    )
} else {
    message("Skipping clustered phylomorphospace plot for Order Diversity (CLR pPCA) - Clustering or PCA failed.")
}

# Order Diversity (CLR, Filtered) - Clustered Phylomorphospace
if (!is.null(clustering_result_order_clr_filtered) && !is.null(order_pca_result_clr_filtered)) {
    perform_phylogenetic_pca_and_plot(
        order_df,
        phy_tree,
        "Order_Diversity",
        min_value_threshold = 0.01, # Match filtered PCA settings
        kmeans_result = clustering_result_order_clr_filtered$kmeans_result,
        optimal_k = clustering_result_order_clr_filtered$optimal_k
    )
} else {
    message("Skipping clustered phylomorphospace plot for Order Diversity (CLR pPCA, Filtered) - Clustering or PCA failed.")
}

# Superfamily Diversity (CLR) - Clustered Phylomorphospace
if (!is.null(clustering_result_superfamily_clr) && !is.null(superfamily_pca_result_clr)) {
    perform_phylogenetic_pca_and_plot(
        superfamily_df,
        phy_tree,
        "Superfamily_Diversity",
        min_value_threshold = NULL, # Match original PCA settings
        kmeans_result = clustering_result_superfamily_clr$kmeans_result,
        optimal_k = clustering_result_superfamily_clr$optimal_k
    )
} else {
    message("Skipping clustered phylomorphospace plot for Superfamily Diversity (CLR pPCA) - Clustering or PCA failed.")
}

# Superfamily Diversity (CLR, Filtered) - Clustered Phylomorphospace
if (!is.null(clustering_result_superfamily_clr_filtered) && !is.null(superfamily_pca_result_clr_filtered)) {
    perform_phylogenetic_pca_and_plot(
        superfamily_df,
        phy_tree,
        "Superfamily_Diversity",
        min_value_threshold = 0.01, # Match filtered PCA settings
        kmeans_result = clustering_result_superfamily_clr_filtered$kmeans_result,
        optimal_k = clustering_result_superfamily_clr_filtered$optimal_k
    )
} else {
    message("Skipping clustered phylomorphospace plot for Superfamily Diversity (CLR pPCA, Filtered) - Clustering or PCA failed.")
}

# --- Generate Phylogeny Plots with Optimal Cluster Tip Coloring ---

message("\n--- Generating Phylogeny Plots with Optimal Cluster Tip Coloring ---")

# Plot Order Diversity (CLR) - Phylogeny with Colored Tips
if (!is.null(clustering_result_order_clr) && !is.null(order_pca_result_clr)) {
    plot_clustered_phylogeny(
        full_phy_tree = phy_tree, # Use the original full tree
        kmeans_result = clustering_result_order_clr$kmeans_result,
        optimal_k = clustering_result_order_clr$optimal_k,
        filtered_ids = order_pca_result_clr$filtered_ids,
        title_prefix = "Order_Diversity_CLR_pPCA"
    )
} else {
    message("Skipping phylogeny tip color plot for Order Diversity (CLR pPCA) - Clustering or PCA failed.")
}

# Plot Order Diversity (CLR, Filtered) - Phylogeny with Colored Tips
if (!is.null(clustering_result_order_clr_filtered) && !is.null(order_pca_result_clr_filtered)) {
    plot_clustered_phylogeny(
        full_phy_tree = phy_tree,
        kmeans_result = clustering_result_order_clr_filtered$kmeans_result,
        optimal_k = clustering_result_order_clr_filtered$optimal_k,
        filtered_ids = order_pca_result_clr_filtered$filtered_ids,
        title_prefix = "Order_Diversity_CLR_pPCA_MinValFiltered0.01"
    )
} else {
    message("Skipping phylogeny tip color plot for Order Diversity (CLR pPCA, Filtered) - Clustering or PCA failed.")
}

# Plot Superfamily Diversity (CLR) - Phylogeny with Colored Tips
if (!is.null(clustering_result_superfamily_clr) && !is.null(superfamily_pca_result_clr)) {
    plot_clustered_phylogeny(
        full_phy_tree = phy_tree,
        kmeans_result = clustering_result_superfamily_clr$kmeans_result,
        optimal_k = clustering_result_superfamily_clr$optimal_k,
        filtered_ids = superfamily_pca_result_clr$filtered_ids,
        title_prefix = "Superfamily_Diversity_CLR_pPCA"
    )
} else {
    message("Skipping phylogeny tip color plot for Superfamily Diversity (CLR pPCA) - Clustering or PCA failed.")
}

# Plot Superfamily Diversity (CLR, Filtered) - Phylogeny with Colored Tips
if (!is.null(clustering_result_superfamily_clr_filtered) && !is.null(superfamily_pca_result_clr_filtered)) {
    plot_clustered_phylogeny(
        full_phy_tree = phy_tree,
        kmeans_result = clustering_result_superfamily_clr_filtered$kmeans_result,
        optimal_k = clustering_result_superfamily_clr_filtered$optimal_k,
        filtered_ids = superfamily_pca_result_clr_filtered$filtered_ids,
        title_prefix = "Superfamily_Diversity_CLR_pPCA_MinValFiltered0.01"
    )
} else {
    message("Skipping phylogeny tip color plot for Superfamily Diversity (CLR pPCA, Filtered) - Clustering or PCA failed.")
}

message(paste("\nPhylogenetic PCA and Clustering analysis complete. Plots saved in", output_dir))

