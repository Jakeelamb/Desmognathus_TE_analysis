# Supplemental dataset release

`DATASET_MANIFEST.csv` is the authoritative index. It records the analysis,
title, description, release status, source path, source and output SHA-256
hashes, dimensions, and byte size of every exported CSV.

## Interpretation classes

- `main_candidate`: suitable for a main descriptive or comparative result after
  final manuscript review.
- `supplementary`: valid supporting analysis that is not a primary claim.
- `sensitivity_only`: useful for robustness or hypothesis generation but not a
  stand-alone validated measurement/causal claim.
- `audit_only`: provenance, exclusion, quality, or calibration evidence.

## Non-CSV files

Phylogenetic trees remain in standard Newick/NEXUS formats under `trees/`.
Reviewer-readable node and edge CSVs are provided as S38 and S39. Images, masks,
TE libraries, sequence alignments, and very large hit-level tables should be
deposited in the final DOI-bearing repository rather than forced into lossy or
impractically large spreadsheet files.

The 1.8-GB final-18 and 3.4-GB all-resource hit-level RepeatMasker tables are not
duplicated here. The compact S12 landscape and S13 inventory conserve the
reported species/order summaries; the hit-level tables require external
repository deposition with their existing manifests and hashes.
