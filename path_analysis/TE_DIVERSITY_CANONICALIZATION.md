# TE Diversity Canonicalization

This document resolves the generation path for the TE diversity summary tables
used by the analysis-facing path-analysis layer and records the promotion of that
path into the upstream processing script.

It is intentionally non-destructive. No files in `results/data/` were replaced
as part of this step.

## Bottom Line

The canonical diversity summary tables are now provenance-resolved locally.

As of this cleanup, `scripts/processing/diversity_stats.py` has been promoted as
the canonical upstream writer for:

- `results/data/diversity_order_stats.csv`
- `results/data/diversity_superfamily_stats.csv`

The existing snapshots were not overwritten during this audit step. The
non-destructive candidate builder remains the verification path.

Observed construction path:

1. `scripts/processing/diversity_stats.py` reads:
   - `results/data/dnaPipeTE_order_breakdown.csv`
   - `results/data/dnaPipeTE_superfamily_breakdown.csv`
2. It writes the threshold grids:
   - `results/data/comparison_diversity_order_stats_granular_0_5pct.csv`
   - `results/data/comparison_diversity_superfamily_stats_granular_0_5pct.csv`
3. The canonical summary tables are the corresponding breakdown tables with the
   threshold `0.0` metrics appended as:
   - `Simpson_Diversity <- Simpson_0.0`
   - `Shannon_Diversity <- Shannon_0.0`
   - `Pielou_Evenness <- Pielou_0.0`

That reproduces the current `results/data/diversity_order_stats.csv` and
`results/data/diversity_superfamily_stats.csv` within tight numeric tolerance
without rerunning any HPC work.

## Non-Destructive Proof

Scratch candidates:

- `path_analysis/data/derived/diversity_order_stats_candidate.csv`
- `path_analysis/data/derived/diversity_superfamily_stats_candidate.csv`

Audit summary:

- `path_analysis/data/derived/diversity_canonicalization_audit.csv`

Builder:

- `path_analysis/scripts/build_canonical_diversity_tables.py`

Current audit interpretation:

- row counts match for both datasets
- column layouts match for both datasets
- species order matches for both datasets
- no candidate rows are missing any diversity metric
- candidate tables match the current canonical snapshots within `1e-12`
- the `0.0` threshold metrics also match direct recomputation from the
  breakdown tables within `1e-12`

The candidate CSV text is not byte-identical to the stored snapshots in every
case because of negligible floating-point string formatting differences in some
cells. That is a representation issue, not a scientific-data mismatch.

## Why The Provenance Previously Looked Ambiguous

The ambiguity came from a missing writer step in the visible pipeline docs.

Historically, `scripts/processing/diversity_stats.py` wrote the long-format and
comparison-threshold files, but not the canonical summary snapshots:

- `results/data/diversity_order_stats.csv`
- `results/data/diversity_superfamily_stats.csv`

The removed ad hoc helper path is not the right writer for the current
canonical files because it wrote an untracked `superfamily_proportions.csv`
side path rather than the audited order/superfamily diversity snapshots.

## Current Recommendation

Treat the diversity provenance issue as resolved for analysis and
path-analysis use.

Current safe operating rule:

1. keep the current snapshots as-is unless a deliberate TE rerun is intended
2. use `scripts/processing/diversity_stats.py` for future canonical diversity
   generation
3. use `path_analysis/scripts/build_canonical_diversity_tables.py` when you
   want a non-destructive reconstruction audit against the repo-local snapshots
4. treat the old standalone diversity-helper path as retired; the repo now
   keeps only the canonical writer and the non-destructive audit builder
