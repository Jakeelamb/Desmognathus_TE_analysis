#!/usr/bin/env Rscript
# Quick PGLS analysis for key hypotheses

suppressPackageStartupMessages({
  library(ape); library(caper); library(dplyr); library(readr); library(ggplot2); library(yaml)
})

project_root <- '/home/jake/Projects/Desmognathus_TE'
config <- yaml::read_yaml(file.path(project_root, 'paths.yaml'))
data_dir <- file.path(project_root, config$results$data)
figures_dir <- file.path(project_root, config$results$figures)
phylo_dir <- file.path(project_root, config$input_data$phylogeny)
output_dir <- file.path(figures_dir, 'pgls')
output_data_dir <- file.path(data_dir, 'pgls')
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(output_data_dir, showWarnings = FALSE, recursive = TRUE)

# Load data
tree <- read.tree(file.path(phylo_dir, 'desmo900dated_test.tre'))
order_df <- read_csv(file.path(data_dir, 'dnaPipeTE_order_breakdown.csv'), show_col_types = FALSE)

# Prepare data
order_df[[1]] <- gsub('^D\\.', '', order_df[[1]])
order_df <- order_df %>% rename(species = 1) %>% filter(species %in% tree$tip.label)
pruned_tree <- drop.tip(tree, setdiff(tree$tip.label, order_df$species))
comp_data <- comparative.data(pruned_tree, as.data.frame(order_df), names.col = species, vcv = TRUE)

results <- list()

# Test multiple models
models <- list(
  c("LINE", "LTR"),
  c("TIR", "LTR"),
  c("TIR", "LINE"),
  c("DIRS", "LTR"),
  c("Helitron", "TIR")
)

for (m in models) {
  response <- m[1]
  predictor <- m[2]
  formula_str <- paste(response, "~", predictor)

  tryCatch({
    pgls_result <- pgls(as.formula(formula_str), data = comp_data, lambda = "ML")
    s <- summary(pgls_result)

    results[[formula_str]] <- data.frame(
      model = formula_str,
      lambda = pgls_result$param["lambda"],
      r_squared = s$r.squared,
      coefficient = s$coefficients[predictor, "Estimate"],
      std_error = s$coefficients[predictor, "Std. Error"],
      t_value = s$coefficients[predictor, "t value"],
      p_value = s$coefficients[predictor, "Pr(>|t|)"]
    )
    cat("Completed:", formula_str, "\n")
  }, error = function(e) {
    cat("Failed:", formula_str, "-", e$message, "\n")
  })
}

result_df <- bind_rows(results)
result_df$significant <- result_df$p_value < 0.05
write_csv(result_df, file.path(output_data_dir, 'pgls_key_results.csv'))

# Create plot
if (nrow(result_df) > 0) {
  p <- ggplot(result_df, aes(x = reorder(model, -p_value), y = -log10(p_value), fill = significant)) +
    geom_bar(stat = "identity") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
    coord_flip() +
    scale_fill_manual(values = c("FALSE" = "gray70", "TRUE" = "steelblue")) +
    labs(title = "PGLS Results: TE Order Correlations",
         x = "Model", y = "-log10(p-value)") +
    theme_minimal()
  ggsave(file.path(output_dir, "pgls_results.png"), p, width = 8, height = 5, dpi = 300)
}

cat("\n=== PGLS Results ===\n")
print(result_df)
cat("\nSaved to:", file.path(output_data_dir, 'pgls_key_results.csv'), "\n")
