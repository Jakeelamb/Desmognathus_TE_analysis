# Input Preparation Plan

This document is the working plan for preparing path-analysis inputs while the final independent genome-size estimates are still in progress. It focuses on two things:

1. trait acquisition from literature and online data sources
2. reduction of the transposable-element data into interpretable, model-usable variables

The goal is not to assemble the largest possible trait table. The goal is to build a small, sourceable, phylogenetically defensible input set that can support a manuscript-quality path analysis.

## Current Position

The current consolidated workspace already has these blocks:

- TE composition and diversity
  - `results/data/dnaPipeTE_class_breakdown.csv`
  - `results/data/dnaPipeTE_order_breakdown.csv`
  - `results/data/dnaPipeTE_superfamily_breakdown.csv`
  - `results/data/diversity_order_stats.csv`
  - `results/data/diversity_superfamily_stats.csv`
- TE turnover / deletion metrics
  - `results/data/divergence/divergence_summary_statistics_by_species.csv`
- Ectopic recombination proxy
  - `results/data/ectopic_recombination_filtered_3000bp_5+domains_no_unknown_species.csv`
- Genome and morphology
  - `~/Projects/cellprofiler_test/output/qc_report_blockbalanced/final_species_results.csv`
  - `~/Projects/cellprofiler_test/output/publication_analysis/species_morphology_summary.csv`

What is still missing for the chapter-level path analysis is the organismal trait layer: body size, external morphology, development, life history, lifestyle, and basic range/ecology covariates.

## Design Rules

- Use one observed proxy per conceptual block in the main DAGs.
- Prefer DOI-backed datasets, revision appendices, and peer-reviewed supplements over tertiary species accounts.
- Treat recent taxonomy as a first-class problem. Many current Desmognathus names are absent from older trait databases because of 2022-2023 revisions.
- Keep raw extraction separate from curated species summaries.
- Do not mix measurement standards without flags. Total length, SVL, and body mass are not interchangeable.
- Do not let a variable into the main path model unless it has adequate coverage across the target species subset.

## Immediate Deliverables

These should exist before we do any serious trait mining:

1. `data/templates/species_taxonomy_crosswalk.csv`
2. `data/templates/trait_registry.csv`
3. `data/templates/literature_trait_extraction.csv`
4. `data/templates/source_manifest.csv`

Those templates are added alongside this plan.

## Traceability Rules

This project needs publication-grade provenance.

That means:

- every source gets a stable `source_id`
- every extracted trait row must carry that `source_id`
- every curated species summary row must retain the source ids that contributed to it
- every taxonomic name remapping must cite the revision source that justified it
- every web-derived source should keep the exact URL and the access date
- if a trait is inferred rather than directly reported, the inference must be labeled in notes

No trait should enter the final path-analysis table without a traceable provenance chain.

### Required provenance fields

At minimum, every extracted record should preserve:

- `source_id`
- source citation
- DOI or URL
- source type
- original taxon name
- page, table, appendix, or dataset location when available
- access date for web sources
- curation notes

### Source manifest workflow

1. Add the source once to `source_manifest.csv`.
2. Reuse the same `source_id` across all extraction rows.
3. When collapsing raw rows into species summaries, keep a semicolon-delimited list of contributing `source_id` values.
4. If a value came from a revision supplement and was validated against a species account, retain both ids.

This avoids the common failure mode where a summary table exists but no one can reconstruct which source produced which row.

## Trait Blocks To Build

### Tier 1: highest priority for the main manuscript model

These are the organismal blocks most likely to be both biologically useful and realistically fillable across the current overlap set.

#### 1. Body size

Preferred variable:

- `adult_svl_mm`

Preferred hierarchy:

1. adult female mean SVL
2. adult mean SVL, sex pooled
3. maximum adult SVL
4. total length only if SVL is unavailable, and only with a permanent flag that it is not directly comparable

Why this matters:

- body size is a likely upstream covariate for cell size and life history
- it is widely reported across Desmognathus and should be more complete than fecundity or age at maturity

#### 2. Lifestyle / habitat association

Preferred variables:

- `aquaticity_index`
- `microhabitat_class`

Recommended coding:

- `aquaticity_index`
  - `0` = terrestrial direct developer
  - `1` = streamside / seep-associated, mostly terrestrial adult
  - `2` = semi-aquatic stream adult
  - `3` = swamp / wetland / strongly aquatic adult
- `microhabitat_class`
  - `mountain_woodland`
  - `streamside`
  - `stream_aquatic`
  - `swamp_coastal_plain`
  - `cave_or_subterranean` if ever relevant

Why this matters:

- this is likely a better first-pass ecological covariate than many sparse life-history variables
- it links naturally to development mode, body size, and potential selection on cell physiology

#### 3. Developmental mode

Preferred variables:

- `development_mode`
- `larval_period_months` if available

Recommended coding:

- `development_mode`
  - `direct_development`
  - `aquatic_larva`

Why this matters:

- direct development is already central in your slides and in plethodontid life-history thinking
- it is easier to score consistently than many quantitative reproductive traits

#### 4. Elevation

Preferred variables:

- `elevation_min_m`
- `elevation_max_m`
- `elevation_mid_m`

Why this matters:

- elevation is often available even when detailed demography is not
- it gives a simple environmental axis that may align with microhabitat and developmental strategy

### Tier 2: strong candidates if coverage is adequate

#### 5. Reproductive output

Preferred variables:

- `clutch_size_mean`
- `egg_diameter_mm`

These are biologically attractive but will likely be patchy. They should be collected now, but only promoted into the main DAG if coverage is high enough.

#### 6. Developmental timing

Preferred variables:

- `age_at_maturity_months`
- `incubation_or_hatching_time_days`

Same rule as above: collect now, include later only if the coverage is real.

### Tier 3: useful, but probably better as supplement or sensitivity analysis

#### 7. External morphology beyond size

Preferred approach:

- derive a size-corrected species-level shape axis rather than putting many raw linear measurements into the DAG

Recommended source type:

- revision datasets with individual-level morphometrics

Preferred derived variable:

- `shape_pc1_resid_svl`

Use this only if enough species can be assembled with comparable raw measurements. Otherwise, keep morphology as categorical lifestyle classes plus body size.

## Source Hierarchy

### Primary source reservoirs

These are the first places to mine because they are explicit, citable, and likely to be reusable.

1. AmphiBIO
   - broad amphibian ecological and reproductive coverage
   - DOI-backed and downloadable
   - especially useful for coarse habitat, breeding strategy, development, and body-size coverage
   - caveat: for caudates, body size is often total length rather than SVL, and missing values mean "not available", not "absent"

2. Recent Desmognathus revision papers and their supplements
   - best source for modern names and for morphometric data on recently split taxa
   - especially important because older databases may still use pre-split names
   - likely sources for SVL, shape measurements, and ecology notes

3. Revision-linked data repositories
   - Zenodo / Dryad appendices for morphometrics and specimen metadata
   - these may be more usable than the papers themselves for constructing species means

4. Species-specific life-history literature
   - Bruce and related Desmognathus natural-history papers are likely the best source for clutch size, age structure, and reproductive timing

### Secondary source reservoirs

Use these to fill gaps or validate coarse categories:

- AmphibiaWeb
- IUCN Red List species accounts
- USGS range products
- GBIF / VertNet occurrence data for elevation, after cleaning

Secondary sources are appropriate for habitat category, elevation range, and taxonomic cross-checks, but not as the only source for a quantitative manuscript variable when a primary source exists.

## Taxonomy Problem: handle this first

This is not optional. Many current Desmognathus species in the path-analysis overlap set were described or elevated recently. A trait value reported under an older name may belong to:

- the exact current species
- the broader pre-split complex
- a mixed lineage concept that is no longer taxonomically clean

Therefore, the first table to populate is `species_taxonomy_crosswalk.csv`.

Each row must capture:

- current species name used in this repo
- primary search name
- legacy or broader search names
- species complex
- revision source used to justify the mapping
- whether the older source is acceptable as species-level, complex-level, or unusable

If a trait comes from a pre-split source and cannot be localized confidently, it should either:

- be excluded from the main analysis, or
- be carried with a `complex_level_only` flag for sensitivity analyses

## Trait Curation Rules

- Adults only unless the trait is explicitly developmental.
- Wild individuals preferred over captive.
- Record sex, life stage, and summary statistic when available.
- Record the original taxon name exactly as printed in the source.
- Keep units in the raw extraction table exactly as reported.
- Standardize only in the curated species summary table.
- If multiple populations are reported, retain locality-level raw rows first, then calculate a species summary with provenance.
- Never silently convert total length to SVL.
- For categorical traits, keep the source wording in notes before converting to a scored state.
- If a value is backfilled from a secondary source, mark it explicitly.
- If a value is inferred from text, mark it explicitly and keep the original wording.

## Coverage Targets

For a trait to enter the main path model:

- target at least `20` species in the `te_genome` overlap set of `27`
- ideally at least `14` species in the integrated `te_genome_morphology` set of `18`
- no obvious taxonomic ambiguity for more than a small minority of those rows

Traits that do not hit these targets can still be valuable for:

- exploratory PGLS
- reduced-subset sensitivity analyses
- supplementary descriptive figures

## TE Data: how to make it model-usable

### What the current audit says

From the current local data:

- order means are dominated by `LTR` (`~61.9%`), `LINE` (`~17.3%`), and `TIR` (`~17.3%`)
- naive CLR PCA on all nine orders is dominated by rare components like `Maverick` and `Helitron`
- even a filtered major-order PCA remains a relatively abstract balance axis
- weighted order-level divergence and weighted deletion metrics are highly correlated (`r ~ 0.92`)

This means:

- a default "PCA of everything" is not a good main-path predictor
- turnover and DNA-loss proxies should not both enter the same small DAG as independent predictors

### TE concepts to represent

The TE block should be split into conceptual variables, with one observed proxy per concept:

1. composition
2. diversity / evenness
3. turnover / age structure
4. DNA loss
5. ectopic recombination

### Recommended TE variables

#### Main composition variable

Primary:

- `ltr_line_logratio`

Why:

- simple
- interpretable
- driven by abundant orders rather than rare ones
- already aligns with current scaffold work

Sensitivity alternatives:

- `retro_dna_logratio`
- filtered CLR PC1 using only major orders such as `LTR`, `LINE`, `TIR`, `DIRS`, `SINE`

Do not use the full nine-order CLR PC1 as the main manuscript predictor.

#### Diversity variable

Primary:

- `order_pielou`

Why:

- already computed
- interpretable as compositional evenness
- lower dimensional than any raw order block

#### Turnover / age variable

Primary candidate:

- weighted species mean of order-level `percent_divergence_median` at threshold `0.9`

Sensitivity alternatives:

- `LTR percent_divergence_median`
- `LINE percent_divergence_median`

Keep the threshold fixed across all species and models.

#### DNA loss variable

Primary candidate:

- weighted species mean of order-level `percent_deletions_median` at threshold `0.9`

Important:

- because weighted deletion and weighted divergence are so strongly correlated in the current data, use one or the other in a given candidate family, not both

#### Ectopic recombination variable

Primary:

- `log10(ectopic_mean_ratio)`

QC filters:

- minimum element support threshold
- complete-element fraction threshold
- explicit note that the source file is tab-delimited despite the `.csv` extension

### Presentation strategy for TE data

The TE block should have a separate figure strategy from the model-input strategy.

For the manuscript main text:

- order-level stacked composition bars or heatmap
- one simple divergence / age-structure figure
- one ectopic recombination figure if that mechanism remains competitive

For supplement:

- superfamily composition plots
- phylogenetic PCA / phylomorphospace figures
- full TE landscape plots

The existing `scripts/processing/phylogenetic_pca_analysis.R` and the figures in `results/figures/*pPCA*` are still useful, but they should stay exploratory or supplementary unless we deliberately redefine the TE PCA to avoid rare-order domination.

## Recommended V1 Path-Analysis Input Set

This is the strongest first-pass target once the final genome estimates are ready:

- `genome_true`
- `nucleus_area`
- `cell_area`
- `adult_svl_mm`
- `aquaticity_index`
- `development_mode`
- `elevation_mid_m`
- `ltr_line_logratio`
- `order_pielou`
- one of:
  - `weighted_te_divergence`
  - `weighted_te_deletions`
- `ectopic_index`

That is already enough for a serious staged comparative analysis.

## Workflow

1. Populate the taxonomy crosswalk.
2. Populate the trait registry only with variables we truly intend to curate.
3. Enter literature values into the long-form extraction table.
4. Collapse raw extractions into a curated species summary with flags.
5. Derive TE features in a dedicated script rather than inside the model-fitting script.
6. Join curated traits to the existing `master_species_table.csv`.
7. Freeze a versioned `path_analysis_input_master.csv`.
8. Build path-analysis subsets from that frozen table.

## Next Scripts To Add

When we move from planning to implementation, these should live in `path_analysis/scripts/`:

- `prepare_literature_traits.py`
- `prepare_te_features.py`
- `build_path_input_master.py`

The current `path_model_scaffold.R` should remain focused on transformed analysis inputs and model comparison, not on raw data assembly.
