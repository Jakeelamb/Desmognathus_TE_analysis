#!/usr/bin/env Rscript
# Check R package availability for Desmognathus TE Analysis.
#
# Dependencies should be installed through the versioned Dusky environment, not
# by mutating a user's R library during analysis runs.

cat("Checking required R packages...\n\n")

required_packages <- c(
  "stats",
  "dplyr",
  "ggplot2",
  "cluster",
  "factoextra",
  "readr",
  "tidyr",
  "compositions",
  "ggrepel",
  "ape",
  "phytools",
  "geiger",
  "viridis",
  "gridExtra",
  "tidyverse",
  "ggtree",
  "phangorn"
)

failed <- c()
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    failed <- c(failed, pkg)
    cat("  FAIL", pkg, "\n")
  } else {
    cat("  OK  ", pkg, "\n")
  }
}

if (length(failed) > 0) {
  cat("\nMissing R packages:\n")
  cat(paste("  -", failed, collapse = "\n"), "\n")
  cat("\nUpdate the Conda environment with `conda env update -f Dusky.yml`.\n")
  quit(status = 1)
}

cat("\nAll required R packages are available.\n")
quit(status = 0)
