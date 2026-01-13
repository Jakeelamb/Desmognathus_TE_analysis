#!/usr/bin/env Rscript
# Install missing R packages for Desmognathus TE Analysis

cat("Checking and installing required R packages...\n\n")

# List of required packages
required_packages <- c(
    "stats",       # Base R - for prcomp, kmeans
    "dplyr",       # Data manipulation
    "ggplot2",     # Plotting
    "cluster",     # For silhouette
    "factoextra",  # PCA visualization
    "readr",       # Reading CSV files
    "tidyr",       # Data tidying
    "compositions",# CLR transformation
    "ggrepel",     # Text labels in plots
    "ape",         # Phylogenetic analyses
    "phytools",    # Phylogenetic tools
    "geiger",      # Phylogenetic comparative methods
    "viridis",     # Color palettes
    "gridExtra",   # Grid layouts for plots
    "tidyverse",   # Collection of data science packages
    "ggtree",      # Phylogenetic tree visualization
    "phangorn"     # Phylogenetic analysis
)

# Check which packages are already installed
installed <- installed.packages()[, "Package"]
to_install <- required_packages[!required_packages %in% installed]

if (length(to_install) == 0) {
    cat("✓ All required R packages are already installed!\n")
} else {
    cat("Installing missing packages:\n")
    cat(paste("  -", to_install, collapse = "\n"), "\n\n")
    
    # Set CRAN mirror
    options(repos = c(CRAN = "https://cloud.r-project.org"))
    
    # Install packages from CRAN
    cran_packages <- setdiff(to_install, c("ggtree"))
    if (length(cran_packages) > 0) {
        install.packages(cran_packages, dependencies = TRUE, quiet = FALSE)
    }
    
    # Install ggtree from Bioconductor if needed
    if ("ggtree" %in% to_install) {
        cat("\nInstalling ggtree from Bioconductor...\n")
        if (!requireNamespace("BiocManager", quietly = TRUE)) {
            install.packages("BiocManager")
        }
        BiocManager::install("ggtree")
    }
    
    cat("\n✓ Installation complete!\n")
}

# Verify all packages can be loaded
cat("\nVerifying package installation...\n")
failed <- c()
for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        failed <- c(failed, pkg)
        cat("  ✗", pkg, "- FAILED\n")
    } else {
        cat("  ✓", pkg, "\n")
    }
}

if (length(failed) > 0) {
    cat("\n⚠ The following packages failed to install:\n")
    cat(paste("  -", failed, collapse = "\n"), "\n")
    cat("\nYou may need to install system dependencies or try manual installation.\n")
    quit(status = 1)
} else {
    cat("\n✓ All packages verified successfully!\n")
    quit(status = 0)
}

