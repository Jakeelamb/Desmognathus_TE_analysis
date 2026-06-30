# TE Data Audit

This document records the current audit status of the transposable element data
used by the `path_analysis/` workflow.

It is intentionally a local-data audit, not an HPC recomputation log. The goal
is to determine whether the repo-local TE outputs already present in this repo are
internally consistent, traceable, and suitable for research use and for the
current phylogenetic path-analysis layer.

## Scope

Audited on 2026-03-08 against the repo-local TE inputs already on disk:

- `results/data/dnaPipeTE_order_breakdown.csv`
- `results/data/dnaPipeTE_superfamily_breakdown.csv`
- `results/data/diversity_order_stats.csv`
- `results/data/diversity_superfamily_stats.csv`
- `results/data/divergence/divergence_summary_statistics_by_species.csv`
- `results/data/ectopic_recombination_filtered_3000bp_5+domains_no_unknown_species.csv`

Downstream products audited:

- `path_analysis/data/derived/te_path_features.csv`
- `path_analysis/data/derived/te_model_feature_panel.csv`
- `path_analysis/data/derived/path_input_master.csv`
- `path_analysis/data/derived/analysis_panel_summary.csv`
- `path_analysis/data/derived/panels/*.csv`

Reproducible audit artifacts:

- `path_analysis/scripts/audit_te_data.py`
- `path_analysis/data/derived/te_audit_summary.csv`
- `path_analysis/data/derived/te_audit_edge_cases.csv`

## Bottom Line

The current TE data are in good enough shape to use for the current comparative
and path-analysis pipelines, with a few explicit caveats.

Current audit outcome:

- `19` pass
- `1` warn
- `0` fail

The main conclusion is that the stored TE feature tables can be reconstructed
exactly from the repo-local source tables without rerunning the original
HPC analyses.

The ectopic support cleanup performed during this audit did not change current
linked-genome panel membership. The TE + linked genome and TE + linked genome
+ ectopic panel sizes are `18` and `16`, respectively.

## What Is Validated

- Source hashes match the values stored in the derived TE tables for the order,
  diversity, divergence, and ectopic inputs.
- The order-breakdown table is numerically sane: `34` species, no negative
  values, and every species row sums to `100` within floating-point noise.
- `Shannon_Diversity`, `Pielou_Evenness`, and `Simpson_Diversity` in
  `results/data/diversity_order_stats.csv` match the order table exactly.
- The thresholded divergence table used by the path-analysis features is
  complete for the current scope: `34` species, `10` TE order bins, `340` rows, no
  missing turnover/deletion/insertion medians, and no species-order duplicates.
- `path_analysis/data/derived/te_path_features.csv` reconstructs exactly from
  the stored source tables for all compared derived TE columns.
- The ectopic support summaries now distinguish:
  - total ectopic rows seen per species: `ectopic_n_rows_total`
  - usable ratio-bearing rows: `ectopic_n_elements`
  - rows with known complete/not-complete state: `ectopic_n_complete_known`
- `path_analysis/data/derived/te_model_feature_panel.csv` matches the expected
  subset and readiness flags from `te_path_features.csv`.
- `path_analysis/data/derived/path_input_master.csv` preserves the TE feature
  values exactly after merge. For overlapping provenance fields, the canonical
  TE-feature versions live in suffixed `*_te` columns.
- `path_analysis/data/derived/analysis_panel_summary.csv` matches the current
  panel files on disk.
- All TE source ids are present in
  `path_analysis/data/templates/source_manifest.csv`, and TE-specific
  traceability gaps are absent from
  `path_analysis/data/derived/source_traceability_gaps.csv`.
- The overlap datasets and the current medium-plus panel files still
  agree on species membership for the core TE families.

## Current Caveats

### 1. One ectopic row has a missing ratio

`aureatus` has one row in the filtered ectopic source table with
`ratio_terminal_internal = NA`.

Current implication:

- it is retained in `ectopic_n_rows_total`
- it is excluded from `ectopic_n_elements`
- it does not affect `ectopic_mean_ratio`

This does not currently break any path-analysis panel, but it is worth knowing
if `ectopic_n_elements` is reported as a support metric.

### 2. One ectopic row has `Complete = unknown`

`anicetus` has one row in the filtered ectopic source table where `Complete` is
`unknown`.

Current implication:

- the current aggregation uses only known `yes` and `no` rows in the
  denominator
- the row is therefore excluded from `ectopic_complete_fraction`
- `ectopic_n_complete_known` makes the denominator explicit

This is minor and now handled explicitly, but it is still worth documenting if
`ectopic_complete_fraction` is discussed in the results.

## Missingness That Is Explicit Rather Than Suspicious

The extended TE panel is missing ectopic support for exactly four species:

- `catahoula`
- `kanawha`
- `lycos`
- `valtos`

These species are absent from the stored ectopic source table itself. This is
clean missingness, not a downstream merge failure.

Current panel consequences:

- TE core panel species: `34`
- TE extended panel species: `30`
- TE + linked genome overlap: `18`
- TE + linked genome + ectopic overlap: `16`

## Canonical Files To Trust

For current comparative and path-analysis use, treat these as the canonical TE layer:

1. `path_analysis/data/derived/te_path_features.csv`
2. `path_analysis/data/derived/te_model_feature_panel.csv`
3. `path_analysis/data/derived/path_input_master.csv`
4. `path_analysis/data/derived/panels/*.csv`

Important detail:

- for overlapping provenance fields inside `path_input_master.csv`, prefer the
  suffixed `*_te` columns when you specifically want the TE-feature-table
  provenance rather than the master-dataset provenance

## What This Audit Does Not Prove

This audit does not revalidate the original computationally expensive upstream
TE analyses from raw assemblies or HPC jobs.

What it does prove is narrower and still useful:

- the repo-local TE source tables are present
- the derived path-analysis TE tables are reproducible from those repo-local inputs
- the joins, panel counts, and provenance references currently line up

So the remaining trust boundary is upstream of these stored source tables, not
within the current path-analysis staging layer.

## Recommended Next Steps

- For methods writing, use the wording in `TE_METHODS_LANGUAGE.md` rather
  than paraphrasing these caveats from memory.
- Treat `te_path_features.csv` and `te_model_feature_panel.csv` as versioned
  TE inputs unless a deliberate upstream TE rerun is planned.
- If `ectopic_n_elements` or `ectopic_complete_fraction` will be emphasized in
  the results, report their current semantics explicitly:
  usable-ratio counts exclude missing ratios, and complete fractions exclude
  `Complete = unknown`.
- Keep the current distinction clear between validated local post-processing and
  not-rerun upstream HPC generation.
- If you want to go one layer deeper, the next audit target is the upstream
  provenance chain from assemblies/HPC jobs into the repo-local summary
  tables.
