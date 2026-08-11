# Supplemental dataset release

`DATASET_MANIFEST.csv` is the authoritative index. It records the analysis,
title, description, release status, primary source path, source and output
SHA-256 hashes, dimensions, and byte size of every exported CSV. For joined
products, `additional_source_paths` and `additional_source_sha256` record every
other input in the same order.

## Interpretation classes

- `main_candidate`: suitable for a main descriptive or comparative result after
  final manuscript review.
- `supplementary`: valid supporting analysis that is not a primary claim.
- `sensitivity_only`: useful for robustness or hypothesis generation but not a
  stand-alone validated measurement/causal claim.
- `audit_only`: provenance, exclusion, quality, or calibration evidence.

S40 is the exact 21-species TE34 ∩ Path24 sensitivity table. It joins
classified-superfamily Shannon entropy to relative nuclear IOD and its
conditional interval. It contains no pg/Gb conversion and does not represent
an independently validated absolute genome-size measurement.

S08 contains four complete 34-species diversity strata: order and superfamily
summaries under classified-conditional and unresolved-mass-aware composition.
Its manifest row binds both compact composition matrices and the mass-accounting
table. Natural-log Shannon entropy and the explicitly named Gini-Simpson index
are the only paper-facing diversity indices. The classified-conditional
order-level stratum is primary; superfamily and unresolved-mass-aware strata are
sensitivity/audit views. Observed richness is an audit/support count rather than
a diversity endpoint. A positive `Unresolved` bin is an accounting category,
not a biological TE taxon. Rare features excluded from the shared-feature PCA
remain in S08 diversity calculations.

S14-S17 are bound to
`analyses/01_transposable_elements/provenance/ltr_terminal_internal_release.json`.
S14 retains all 1,086 usable rows from the historical pipe-split selection for
audit. It corrects `domain_count` by counting whitespace-delimited TEsorter
`DOMAIN|MODEL` annotations and explicitly marks 408 LTR/Gypsy elements with at
least five domains as eligible; the other 678 usable rows are audit-only. The
within-species, inclusive two-sided Tukey IQR rule retains 381 eligible
elements. The primary S15 branch intersects them with the independent gate
requiring positive depth across at least 80% of left-LTR, right-LTR, and
internal positions, yielding 380 elements across all 30 species. S16 separates
historical pipe-split counts from corrected selection and support counts. S17's
two corrupt historical rows each have four real domains and are not exclusions
from the current cohort. These tables describe a mapping/deletion-footprint
sensitivity, not a direct ectopic-recombination, solo-LTR, or DNA-loss rate. A
later compound one-element screen is preserved only as a retired reconstruction
and is not the original IQR method.

S04, S26, and S27-S37 express the Path24 image phenotype as relative nuclear IOD.
The preserved model source files used a conditional *D. fuscus*-scaled pg
alias, but log transformation followed by standardization makes those saved
fits numerically invariant to that one positive scale factor. Publication
exports therefore use `relative_iod`/`iod` labels and retain the exact saved
coefficients, rankings, and sensitivity results without implying an absolute
genome-size measurement.

S19 includes animal support for every morphology estimate:
`n_size_specimens`, `largest_specimen_n`, and `largest_specimen_fraction`.
Seven species are represented by one animal, so its paired-object bootstrap
intervals are conditional on the observed animals and images rather than
population-level among-animal intervals.

S23 is named `Supplementary_Data_S23_nuclear_iod_by_image.csv` and contains
unnormalized image-level nuclear-IOD summaries. Its `specimen_id` field denotes
one animal; it is not a ratio table. S22 is the species-level relative index.
S24 exposes the declared image-quality balance gate and its four failures;
single-image support and failed-balance species remain visible in Figure 6.

## Non-CSV files

Phylogenetic trees remain in standard Newick/NEXUS formats under `trees/`.
Reviewer-readable node and edge CSVs are provided as S38 and S39. Images, masks,
TE libraries, sequence alignments, and very large hit-level tables should be
deposited in the final DOI-bearing repository rather than forced into lossy or
impractically large spreadsheet files.

The 1.8-GB historical final-18 and 3.4-GB all-resource hit-level RepeatMasker
tables are not duplicated here. The compact S12 landscape and S13 inventory
conserve the reported species/order summaries; the hit-level tables require
external repository deposition with their existing manifests and hashes.
