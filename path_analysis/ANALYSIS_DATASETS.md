# Analysis Datasets

This note defines the intended panel files for the comparative path-analysis stage.

## Panel philosophy

- `all` panels preserve the full overlap set for each model family.
- `primary` panels require the core organismal block to be present.
- `primary_mediumplus` panels additionally require medium-or-better confidence for body size, development, and lifestyle.
- `primary_strict_body` panels additionally remove total-length and mixed-stage body-size proxies.

These panels are intended to support a reproducible workflow where:

- the `all` set is the broad sensitivity analysis
- the `primary_mediumplus` set is the most defensible organismal extension set

## Phylogenetic imputation layer

The path-input master table may also contain `phylo_` columns from
`data/derived/organismal_traits_phylo_inference.csv`.

Those columns are not part of the default panel logic.

They exist only to support:

- sensitivity analyses for missing or low-confidence organismal traits
- explicit audits showing where the phylogeny supports or contradicts weak observations

Rules:

- never replace a source-backed observed value silently with a `phylo_` value
- treat `phylo_inferred_for_missing` as inferred, not reported
- treat `low_confidence_observed_supported_by_phylogeny` as corroboration, not a source upgrade by itself

## Generated files

All panel files are written to `data/derived/panels/`.
The panel summary is written to `data/derived/analysis_panel_summary.csv`.
The per-species inclusion audit is written to `data/derived/analysis_species_readiness.csv`.

- `te_genome_all.csv`
  Tree + TE + genome overlap.
- `te_genome_primary.csv`
  `te_genome_all` filtered to complete core organismal traits.
- `te_genome_primary_mediumplus.csv`
  `te_genome_primary` filtered to medium-plus confidence for body size, development, and lifestyle.
- `te_genome_organismal_primary_mediumplus.csv`
  Organismal-augmented TE + genome family on the medium-plus organismal subset.
- `te_genome_primary_phylofill.csv`
  TE + genome subset using explicit phylogenetic sensitivity fills for otherwise missing or low-confidence organismal traits.
- `te_genome_organismal_primary_phylofill.csv`
  Organismal-augmented TE + genome family with explicit phylogenetic sensitivity fills.
- `te_genome_ectopic_all.csv`
  Tree + TE + genome + ectopic overlap.
- `te_genome_ectopic_primary_mediumplus.csv`
  Ectopic family subset with medium-plus organismal confidence.
- `te_genome_ectopic_organismal_primary_mediumplus.csv`
  Ectopic-plus-organismal TE + genome family on the medium-plus organismal subset.
- `te_genome_ectopic_primary_phylofill.csv`
  Ectopic family subset with explicit phylogenetic sensitivity fills.
- `te_genome_ectopic_organismal_primary_phylofill.csv`
  Ectopic-plus-organismal family with explicit phylogenetic sensitivity fills.
- `te_genome_ectopic_primary_strict_body.csv`
  Ectopic family subset restricted to adult-oriented body-size proxies.
- `te_genome_ectopic_organismal_primary_strict_body.csv`
  Ectopic-plus-organismal family restricted to adult-oriented body-size proxies.
- `te_genome_morphology_all.csv`
  Tree + TE + genome + morphology overlap.
- `te_genome_morphology_primary_mediumplus.csv`
  Morphology family subset with medium-plus organismal confidence.
- `te_genome_morphology_primary_phylofill.csv`
  Morphology family subset with explicit phylogenetic sensitivity fills.
- `te_genome_morphology_primary_strict_body.csv`
  Morphology family subset restricted to adult-oriented body-size proxies.
- `te_genome_primary_strict_body.csv`
  TE + genome subset restricted to adult-oriented body-size proxies.
- `te_genome_organismal_primary_strict_body.csv`
  Organismal-augmented TE + genome family restricted to adult-oriented body-size proxies.

The current comparison file is:

- `phylofill_panel_comparison.csv`
  Panel-level comparison of observed-only versus `primary_phylofill` species sets and phylo-assisted species counts.

## Embedded flags

Each panel includes:

- confidence columns for body size, development, and lifestyle
- `uses_total_length_body_proxy`
- `uses_mixed_stage_body_proxy`
- `is_low_confidence_body_proxy`
- `uses_total_length_body_proxy`
- source ids for organismal, TE, genome, and morphology blocks

These flags are meant to drive prespecified sensitivity analyses rather than post hoc species dropping.

## Readiness audit

`analysis_species_readiness.csv` stores the panel logic in a reproducibility-auditable form. For each species, it records:

- whether the species is eligible for each panel
- the exclusion reasons when it is not eligible
- the specific body-size proxy class and confidence tier that drove the decision

Current exclusion reason codes include:

- `missing_tree_tip`
- `missing_te`
- `missing_genome`
- `missing_ectopic`
- `missing_morphology`
- `missing_body_size`
- `missing_development`
- `missing_aquaticity`
- `missing_microhabitat`
- `low_body_confidence`
- `low_development_confidence`
- `low_lifestyle_confidence`
- `total_length_body_proxy`
- `mixed_stage_body_proxy`
