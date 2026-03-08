# Data Dictionary

This document maps the conceptual nodes in the path-analysis plan to the current observed variables.

## Core Variables

| Concept | Current column | Current source | Status | Notes |
|---|---|---|---|---|
| Genome size | `genome_size_pg` | `cellprofiler_test/output/qc_report_blockbalanced/final_species_results.csv` | Provisional | Current primary estimate is area-derived. Keep for TE-mechanism work and planning, but avoid final `genome -> nucleus` claims until the independent final estimates are ready. |
| Genome-size uncertainty | `genome_size_se_pg` | same as above | Available | Useful for screening species and sensitivity analyses. |
| Genome QC status | `genome_result_status`, `genome_flag_summary` | same as above | Available | Use for stable-only or caution-exclusion sensitivity runs. |
| Nucleus size | `morph_nucleus_area_um2` | `cellprofiler_test/output/publication_analysis/species_morphology_summary.csv` | Available | Median nucleus area from linked cell+nucleus morphology summary. |
| Cell size | `morph_cell_area_um2` | same as above | Available | Median cell area from linked morphology summary. |
| Cytoplasm area | `morph_cytoplasm_area_um2` | same as above | Available | Optional downstream morphology proxy. |
| N:C ratio | `morph_nc_ratio` | same as above | Available | Derived quantity. Do not include it in the same DAG as both cell area and nucleus area. |
| TE composition, order level | `order_*` columns | `results/data/dnaPipeTE_order_breakdown.csv` | Available | Raw percentages, currently the cleanest TE-composition block. |
| TE composition, superfamily level | `superfamily_*` columns | `results/data/dnaPipeTE_superfamily_breakdown.csv` | Available | Higher-dimensional block. Better for PCA or secondary analyses than for first-pass path models. |
| TE diversity | `order_shannon`, `order_simpson`, `order_pielou` | `results/data/diversity_order_stats.csv` | Available | `order_pielou` is the current preferred evenness summary. |
| Ectopic recombination proxy | `ectopic_mean_ratio`, `ectopic_median_ratio` | aggregated from `results/data/ectopic_recombination_filtered_3000bp_5+domains_no_unknown_species.csv` | Available | Species-level summary of the terminal:internal depth ratio. Keep mechanistic interpretation cautious. |
| Ectopic support | `ectopic_n_elements`, `ectopic_complete_fraction` | same as above | Available | Helpful for filtering weak species. |

## Deferred Variables

| Concept | Current state | Why deferred |
|---|---|---|
| True independent genome size | In progress in `cellprofiler_test` | Needed before formal `genome -> nucleus -> cell` interpretation. |
| DNA loss rate | Not yet consolidated into a stable species-level path-analysis proxy | The divergence/deletion outputs exist, but the exact species-level summary still needs a deliberate definition. |
| TE divergence rate | Not yet consolidated into a stable species-level path-analysis proxy | Same issue as DNA loss. |
| Body size | Not yet staged in this repo | Mentioned in slides, but no curated comparative table lives in `Desmognathus_TE` yet. |
| Life-history strategy | Not yet staged in this repo | Same. |
| Developmental rate | Not yet staged in this repo | Same. |
| Lifestyle / morphology class | Not yet staged in this repo | Can be added later as discrete or ordinal covariates once curated. |

## Preferred First-Pass Derived Variables

These are the planned observed variables used by the current R scaffold:

| Derived variable | Meaning | Construction |
|---|---|---|
| `ltr_balance` | Interpretable TE-composition summary | `log((order_LTR + p) / (order_LINE + p))`, where `p` is a small pseudocount |
| `te_evenness` | TE distribution evenness | Current order-level Pielou evenness |
| `ectopic_index` | Species-level ectopic recombination proxy | `log10(ectopic_mean_ratio)` after positivity checks |
| `gs` | Scaled genome-size variable | Z-score of `log10(genome_size_pg)` |
| `ns` | Scaled nucleus-size variable | Z-score of `log10(morph_nucleus_area_um2)` |
| `cs` | Scaled cell-size variable | Z-score of `log10(morph_cell_area_um2)` |

## Datasets Produced By The Builder

The builder writes the following path-ready tables to `path_analysis/data/derived/`:

| File | Intended use |
|---|---|
| `master_species_table.csv` | Union of all currently available species with availability flags and merged variables |
| `dataset_te_genome.csv` | Core TE-mechanism models |
| `dataset_te_genome_ectopic.csv` | TE + ectopic mechanism models |
| `dataset_genome_morphology.csv` | Nucleotypic planning models |
| `dataset_te_genome_morphology.csv` | Integrated TE -> genome -> morphology planning models |
| `dataset_te_genome_ectopic_morphology.csv` | Most restrictive integrated subset |
| `dataset_overlap_summary.csv` | Species counts and species lists for each staged dataset |
