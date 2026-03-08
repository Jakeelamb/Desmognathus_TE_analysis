# Analysis Datasets

This note defines the intended panel files for the comparative path-analysis stage.

## Panel philosophy

- `all` panels preserve the full overlap set for each model family.
- `primary` panels require the core organismal block to be present.
- `primary_mediumplus` panels additionally require medium-or-better confidence for body size, development, and lifestyle.
- `primary_strict_body` panels additionally remove total-length and mixed-stage body-size proxies.

These panels are intended to support a manuscript workflow where:

- the `all` set is the broad sensitivity analysis
- the `primary_mediumplus` set is the most defensible organismal extension set

## Generated files

All panel files are written to `data/derived/panels/`.

- `te_genome_all.csv`
  Tree + TE + genome overlap.
- `te_genome_primary.csv`
  `te_genome_all` filtered to complete core organismal traits.
- `te_genome_primary_mediumplus.csv`
  `te_genome_primary` filtered to medium-plus confidence for body size, development, and lifestyle.
- `te_genome_ectopic_all.csv`
  Tree + TE + genome + ectopic overlap.
- `te_genome_ectopic_primary_mediumplus.csv`
  Ectopic family subset with medium-plus organismal confidence.
- `te_genome_ectopic_primary_strict_body.csv`
  Ectopic family subset restricted to adult-oriented body-size proxies.
- `te_genome_morphology_all.csv`
  Tree + TE + genome + morphology overlap.
- `te_genome_morphology_primary_mediumplus.csv`
  Morphology family subset with medium-plus organismal confidence.
- `te_genome_morphology_primary_strict_body.csv`
  Morphology family subset restricted to adult-oriented body-size proxies.
- `te_genome_primary_strict_body.csv`
  TE + genome subset restricted to adult-oriented body-size proxies.

## Embedded flags

Each panel includes:

- confidence columns for body size, development, and lifestyle
- `uses_total_length_body_proxy`
- `uses_mixed_stage_body_proxy`
- `is_low_confidence_body_proxy`
- `uses_total_length_body_proxy`
- source ids for organismal, TE, genome, and morphology blocks

These flags are meant to drive prespecified sensitivity analyses rather than post hoc species dropping.
