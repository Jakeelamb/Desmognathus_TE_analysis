# Supplemental dataset release

`DATASET_MANIFEST.csv` is the authoritative index. It records the analysis,
title, description, release status, primary source path, source and output
SHA-256 hashes, dimensions, and byte size of every exported CSV. For joined
products, `additional_source_paths` and `additional_source_sha256` record every
other input in the same order.

`../Desmognathus_figure_data.xlsx` is a generated, indexed cross-sheet view of
the exact non-tree tables directly consumed by the current ggplot figure
builder. Tree inputs stay here in standard Newick and reviewer-readable CSV
formats because those representations do not work well as spreadsheet tabs.
The workbook does not replace this manifest, the S01-S43 CSV release, or the
standard tree files.

## Interpretation classes

- `main_candidate`: suitable for a main descriptive or comparative result after
  final manuscript review.
- `supplementary`: valid supporting analysis that is not a primary claim.
- `sensitivity_only`: useful for robustness or hypothesis generation but not a
  stand-alone validated measurement/causal claim.
- `audit_only`: provenance, exclusion, quality, or calibration evidence.

S40 is the exact 21-species TE34 ∩ Path24 sensitivity table. It joins
canonical top-ten-superfamily Shannon entropy to conditional *D. fuscus*-anchored
genome size and its interval. The table supports a descriptive association on
the selected pg/1C scale; it does not represent an independently measured or
validated absolute genome size.

S41 is a 118-row descriptive cross-study plotting table: 84 exact CC0 rows
from the Wang et al. 2025 Dryad Figure 2 workbook and 34 current
*Desmognathus* rows recomputed over each species' ten most abundant positive
TE superfamilies. The published RepeatModeler/RepeatMasker data and current
dnaPipeTE data come from different genomic inputs and annotation pipelines.
They are displayed for context only; S41 supports no pooled regression,
hypothesis test, or claim that the estimators are equivalent. Current genome
sizes in the 21-species scatter subset are conditional *D. fuscus*-anchored
IOD-derived sensitivities, not direct C-values.

S42 records the 34 Bonsai tip labels and the post-hoc reproductive-strategy
annotation used in Figure 2. S43 is the compact three-row validation surface
for that display: observed exact-side localization, fixed-tree permutation
tests, 200 bounded feature-resampling refits, and the global Fitch statistic.
The strategy labels were supplied after review of the initial PCA, so these
tests are explicitly post hoc.

S08 is the canonical 34-row, four-column diversity table. For each species, the
ten most abundant strictly positive classified superfamilies are selected
independently, with abundance ties resolved deterministically by feature name,
then reclosed to one. The table contains `species`,
`n_superfamilies_included`, natural-log `shannon_entropy`, and the explicitly
named `gini_simpson` index; the included count is ten for every species. Wang et
al. (2025) specify the top-ten-superfamily scope but not whether the retained
abundances were reclosed, so this study states that normalization explicitly.
The Bonsai display and retained PCA audit follow a separate feature rule and use
the 24 superfamilies positive in all 34 species. Unresolved aligned-base mass
enters neither calculation and remains only in the S06-S07 dnaPipeTE
quality-control tables.

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

S04, S26-S37, and S40 use conditional *D. fuscus* assembly-anchored genome size
as the default analysis trait; S38-S39 contain the phylogeny edge and node
tables. The pg/1C values are a positive rescaling of the
underlying nuclear-IOD signal. Log transformation followed by standardization
makes the saved Path24 coefficients, rankings, and sensitivity results
numerically invariant to that scale factor. The genome-size labels therefore
do not imply an independent assay or validated absolute C-value, and the
representative DAG arrows do not establish a unique causal direction.

S22 records the complete derivation of the analysis-facing picogram fields from
the selected published pooled
*D. fuscus* assembly estimate: `16.1 Gbp / 0.978 Gbp per pg = 16.462167689
pg/1C`. Process_413/specimen 32469 and Process_414/specimen 32470 are
author-confirmed standards stained in the same experimental runs as the
unknowns, but their exact standard-to-target run map and the approximately 41%
difference between their median IOD values remain unresolved. The assembly-
derived scale is therefore conditional. S21-S25 retain the raw/object,
image-level, relative-index, image-quality, and calibration evidence; they do
not define a second relative-IOD analysis branch.

S19 includes animal support for every morphology estimate:
`n_size_specimens`, `largest_specimen_n`, and `largest_specimen_fraction`.
Seven species are represented by one animal, so its paired-object bootstrap
intervals are conditional on the observed animals and images rather than
population-level among-animal intervals.

S23 is named `Supplementary_Data_S23_nuclear_iod_by_image.csv` and contains
unnormalized image-level nuclear-IOD summaries. Its `specimen_id` field denotes
one animal; it is not a ratio or genome-size table. S22 is the species-level
calibration bridge. S24 exposes the declared image-quality balance gate and its
four failures; single-image support and failed-balance species remain visible
in the genome-size Figure 6 because they limit the underlying IOD evidence.

## Non-CSV files

Two phylogenetic trees remain in standard Newick/NEXUS formats under `trees/`,
with reviewer-readable Path24 node and edge CSVs provided as S38 and S39. The
third Newick file is the rooted display serialization of the TE34 Bonsai
composition representation used in Figure 2. Its branch lengths are fitted
diffusion times, its root is only a display root, and it is not a species
phylogeny. Images, masks, TE libraries, sequence alignments, and very large
hit-level tables should be deposited in the final DOI-bearing repository rather
than forced into lossy or impractically large spreadsheet files.

The 1.8-GB historical final-18 and 3.4-GB all-resource hit-level RepeatMasker
tables are not duplicated here. The compact S12 landscape and S13 inventory
conserve the reported species/order summaries; the hit-level tables require
external repository deposition with their existing manifests and hashes.
