# TE Compositional PCA Workflow

This document defines the canonical TE PCA workflow for the repository.

The main rule is simple: TE PCA is an exploratory or supplementary ordination
layer, not the primary manuscript predictor block. Primary comparative models
should continue to use interpretable TE variables such as `ltr_line_logratio`,
`order_pielou`, divergence, deletion, and ectopic summaries.

## Canonical Scripts

- `scripts/processing/pca.R`
  Standard compositional PCA entrypoint.
- `scripts/processing/phylogenetic_pca_analysis.R`
  Supplementary phylogenetic PCA entrypoint.
- `scripts/processing/pca_utils.R`
  Shared data-loading, filtering, zero-replacement, CLR, and output helpers.

## Input Data

The PCA workflow uses the frozen repo-local TE breakdown tables:

1. `results/data/dnaPipeTE_order_breakdown.csv`
2. `results/data/dnaPipeTE_superfamily_breakdown.csv`

These are treated as the canonical ordination inputs.

## Analysis Sets

The standard PCA script writes four predefined analyses:

1. `superfamily_primary`
   Superfamily CLR PCA keeping features present in at least 3 species.
2. `superfamily_presence5`
   Sensitivity version keeping features present in at least 5 species.
3. `order_major`
   Sensitivity CLR PCA restricted to `LTR`, `LINE`, `TIR`, `DIRS`, and `SINE`.
4. `order_all`
   Sensitivity CLR PCA on all order-level features.

Interpretation rule:

- `superfamily_primary` is the main supplementary ordination.
- `order_major` and `order_all` are sensitivity ordinations only.
- No raw-proportion PCA is canonical anymore.

## Zero Handling

The workflow does not drop every feature that contains a zero somewhere.

Instead:

1. retain features according to the predefined presence threshold
2. convert rows to closed compositions
3. if zeros remain, replace each zero with `0.5 * min_positive_value`
   measured across the retained composition matrix
4. re-close the composition
5. apply CLR

This keeps the transformation explicit and reproducible.

## PCA Geometry

The canonical ordination geometry is:

- compositional closure
- deterministic zero replacement
- CLR transformation
- `prcomp(center = TRUE, scale. = FALSE)`

Post-CLR variance scaling is intentionally not used, because that would
upweight low-variance rare features.

## Output Products

The standard PCA workflow writes:

- `results/tables/pca/te_pca_analysis_manifest.csv`
- `results/tables/pca/te_pca_pc_metric_correlations.csv`
- `results/tables/pca/te_pca_stability_summary.csv`
- `results/tables/pca/te_pca_stability_detail.csv`
- per-analysis tables for:
  - zero-replaced compositions
  - CLR matrices
  - PCA scores
  - PCA loadings
  - explained variance

The phylogenetic PCA workflow reads those saved CLR matrices directly and writes:

- `results/tables/pca/te_ppca_analysis_manifest.csv`
- per-analysis phylogenetic PCA scores
- per-analysis phylogenetic PCA loadings
- per-analysis phylogenetic explained variance

## Relationship To The Main Paper

The canonical interpretation is:

- use direct TE summary variables for main comparative inference
- use standard compositional PCA as a supplementary ordination summary
- use phylogenetic PCA as a supplementary phylogeny-aware ordination only

If a principal component is discussed in the manuscript, it should be justified
with the saved loading tables, stability summaries, and correlations to direct
TE summary variables.
