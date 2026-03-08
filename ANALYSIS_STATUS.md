# Desmognathus TE Analysis Status

**Last Updated:** 2026-03-08
**Project:** PhD Research - Transposable Element Evolution in Desmognathus Salamanders

For the current paper-facing frozen snapshot, start with `PAPER_FREEZE_MANIFEST.md`.

---

## Completed Analyses

### Data Processing Pipeline

| Script | Purpose | Status |
|--------|---------|--------|
| `scripts/processing/dnaPipe.py` | Process dnaPipeTE output, classify TEs into Class/Order/Superfamily hierarchy | Complete |
| `scripts/processing/repeatmask.py` | Parse RepeatMasker alignments, extract divergence metrics | Complete |
| `scripts/processing/ec.py` | Ectopic recombination analysis via LTR depth ratios | Complete |
| `scripts/processing/divergence.py` | Calculate divergence statistics grouped by TE category | Complete |
| `scripts/processing/diversity_stats.py` | Shannon, Simpson, Pielou's evenness indices | Complete |

### Phylogenetic Analysis

| Script | Purpose | Status |
|--------|---------|--------|
| `scripts/processing/clean_tree_phylo.R` | Validate and clean phylogenetic tree | Complete |
| `scripts/processing/pca.R` | Canonical compositional PCA on frozen TE breakdown tables | Rebuilt |
| `scripts/processing/pca_utils.R` | Shared compositional PCA utilities and matrix builders | Rebuilt |
| `scripts/processing/phylogenetic_pca_analysis.R` | Supplementary phylogenetic PCA from canonical CLR matrices | Rebuilt |
| `scripts/processing/analyze_phylogenetic_signal.R` | Pagel's Lambda, Blomberg's K | Complete |
| `scripts/processing/analyze_phylogenetic_correlogram.R` | Moran's I correlograms from binned phylogenetic distances | Complete |

### Visualization

| Script | Purpose | Status |
|--------|---------|--------|
| `scripts/visualization/divergence.R` | Divergence boxplots and phylogeny bubble plots | Complete |
| `scripts/visualization/hierarchical_donut_TE_diversity.R` | Nested donut charts per species | Complete |
| `scripts/visualization/plot_te_landscape_analysis.R` | TE landscape summary plots from frozen TE tables | Complete |
| `scripts/visualization/plot_phylogenetic_signal.R` | Phylogenetic signal visualization | Complete |
| `scripts/visualization/plot_ectopic_recombination.R` | Canonical ectopic recombination plot + species summary + nonparametric tests | Complete |
| `scripts/visualization/diversity_stats_phylogeny.R` | Phylogeny colored by diversity | Complete |
| `scripts/visualization/donut_TE_diversity.R` | Simple donut charts | Complete |
| `scripts/visualization/plot_simple_phylogeny.R` | Basic phylogeny plots | Complete |
| `scripts/R/visualization/plot_te_landscape.R` | Single species TE landscape | Complete |
| `scripts/visualization/plot_all_te_landscapes.R` | All species TE landscapes | Complete |
| `scripts/visualization/plot_phylogeny_with_te_landscape.R` | Combined phylo + landscape | Repaired; output richness depends on populated `results/landscapes/` |
| `scripts/visualization/superfamily_values_across_tree.R` | Superfamily traits on tree | Complete |

---

## Analyses Added (2026-01-31)

### 1. Trait Evolution Modeling
**Script:** `scripts/processing/trait_evolution.R`
**Purpose:** Compare Brownian Motion vs Ornstein-Uhlenbeck models for TE trait evolution
**Methods:**
- Model fitting with `geiger::fitContinuous()` for BM, OU, EB, and white noise models
- AIC/AICc model comparison with automatic best-model selection
- Ancestral state reconstruction using `phytools::fastAnc()` with confidence intervals
- ContMap visualization for continuous trait evolution on phylogeny
- Rate parameter estimation (sigma^2 for BM, alpha for OU)
**Outputs:**
- `results/data/trait_evolution/evolutionary_model_comparison.csv`
- `results/figures/trait_evolution/model_comparison_summary.png`
- `results/figures/trait_evolution/delta_aicc_bm_ou.png`
- `results/figures/trait_evolution/ancestral_*.png` (per-trait ancestral reconstructions)

### 2. LTR Insertion Age Estimation
**Script:** `scripts/processing/ltr_age_estimation.py`
**Purpose:** Estimate LTR retrotransposon insertion times using intra-element divergence
**Methods:**
- Parse RepeatMasker .align files for LTR element annotations
- Identify paired 5' and 3' LTRs based on proximity and length similarity
- Calculate intra-element divergence as proxy for insertion age
- Convert divergence to age using substitution rate (default: 1.3e-8 subst/site/year)
- Generate age distribution histograms and boxplots by species/class
**Outputs:**
- `results/data/ltr_age/ltr_insertion_ages.csv`
- `results/data/ltr_age/ltr_age_summary_by_species.csv`
- `results/data/ltr_age/ltr_age_summary_by_class.csv`
- `results/figures/ltr_age/ltr_age_distribution_all.png`
- `results/figures/ltr_age/ltr_age_by_class.png`
- `results/figures/ltr_age/ltr_age_by_species.png`
- `results/figures/ltr_age/ltr_divergence_landscapes.png`

### 3. PGLS Regression Framework
**Script:** `scripts/processing/pgls_analysis.R`
**Purpose:** Phylogenetic Generalized Least Squares regression for phylogenetically-corrected correlations
**Methods:**
- `caper::pgls()` with maximum likelihood lambda estimation
- Pairwise correlation tests between all TE order/superfamily pairs
- BH-corrected p-values for multiple testing
- Specific hypothesis tests (LINE~LTR, TIR~Retrotransposons)
- Integration with diversity metrics if available
**Outputs:**
- `results/data/pgls/pgls_order_pairwise.csv`
- `results/data/pgls/pgls_superfamily_pairwise.csv`
- `results/figures/pgls/pgls_lambda_distribution.png`
- `results/figures/pgls/pgls_volcano_plot.png`
- `results/figures/pgls/pgls_rsq_vs_lambda.png`

### 4. Statistical Group Comparisons (PERMANOVA)
**Script:** `scripts/processing/permanova_analysis.R`
**Purpose:** Formal statistical tests for compositional differences between phylogenetic clades
**Methods:**
- Automatic clade definition by tree-cutting (k=5 groups)
- `vegan::adonis2()` PERMANOVA on Bray-Curtis and Euclidean (CLR) distances
- Beta dispersion tests with `betadisper()` to check PERMANOVA assumptions
- Built-in pairwise `adonis2` comparisons with BH-adjusted p-values
- PCoA ordination colored by clade with 95% ellipses
- Distance matrix heatmaps ordered by clade
**Outputs:**
- `results/data/permanova/permanova_summary.csv`
- `results/data/permanova/permanova_order_pairwise.csv`
- `results/data/permanova/permanova_superfamily_pairwise.csv`
- `results/data/permanova/species_clade_assignments.csv`
- `results/data/permanova/distance_matrix_order_bray.csv`
- `results/data/permanova/distance_matrix_superfamily_bray.csv`
- `results/figures/permanova/distance_heatmap_order_bray.png`
- `results/figures/permanova/distance_heatmap_superfamily_bray.png`
- `results/figures/permanova/pcoa_order_bray.png`
- `results/figures/permanova/pcoa_superfamily_bray.png`
- `results/figures/permanova/betadispersion_order.png`
- `results/figures/permanova/betadispersion_superfamily.png`

---

## Planned Future Additions

### Pending Data
- [ ] Genome size integration (data in progress)
- [ ] Morphological trait correlations (data collection needed)

### Potential Analyses
- [ ] Horizontal transfer detection (dS/dN analysis)
- [ ] TE burst timing comparisons across species
- [ ] Network visualization of TE similarity

---

## Data Inventory

| Directory | Contents | Size |
|-----------|----------|------|
| `input_data/dnaPipeTE/` | 35 species TE classifications | ~1.1 GB |
| `input_data/repeatmasker/` | 38 RepeatMasker alignments | ~31 GB |
| `input_data/phylogeny/` | Time-calibrated tree | 1.3 KB |
| `input_data/ectopic_recombination/` | LTR depth files | ~40 GB |
| `input_data/morphological/` | Empty - awaiting data | 0 |

---

## Environment

- **Python:** 3.8+ (pandas, numpy, scipy, biopython, dask, matplotlib, seaborn, tqdm, pyyaml)
- **R:** 4.x (tidyverse, ggplot2, ape, phytools, geiger, caper, vegan, compositions, nlme, pheatmap, factoextra, ggrepel, broom, cluster)
- **Conda:** `Dusky.yml` (updated 2026-01-31 with new dependencies)

### Running the New Analyses

```bash
# Update conda environment (if needed)
conda env update -f Dusky.yml

# Activate environment
conda activate Dusky

# 1. Trait Evolution Modeling (BM vs OU comparison)
Rscript scripts/processing/trait_evolution.R

# 2. LTR Insertion Age Estimation
python scripts/processing/ltr_age_estimation.py

# 3. PGLS Regression Analysis
Rscript scripts/processing/pgls_analysis.R

# 4. PERMANOVA Analysis
Rscript scripts/processing/permanova_analysis.R
```

### Notes on New Scripts
- Active scripts accept either repo-root `paths.yaml` or `config/paths.yaml`
- Output directories are created automatically
- Scripts can be run independently in any order
- R entrypoints now prepend the active conda R library when `CONDA_PREFIX` is set
- R packages will attempt to install from CRAN if missing
