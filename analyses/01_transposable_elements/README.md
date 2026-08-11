# 01 — Transposable elements

This directory is the active TE34 analysis surface. `data/` contains the compact
species-level and element-level products used by the paper; `explore.ipynb`
provides a small editable view over them.

## PCA: start here

```bash
# Read-only end-to-end verification; allow about three minutes.
make te-pca

# Open scores, loadings, variance, diagnostics, and formal test tables.
make te-pca-view
```

`make te-pca` is the sole public PCA check. It first derives feature prevalence,
the complete order- and superfamily-level CLR table, and the final superfamily
ordination from the two compact composition inputs. It then replays the formal
clustering, trait-association, and phylogenetic-signal audit. It compares every
result with the canonical files in `data/` and does not write scratch output.

| Stage | Canonical files under `data/` |
| --- | --- |
| Compact inputs | `te_order_composition.csv`, `te_superfamily_composition.csv` |
| Feature and transform audit | `te_feature_prevalence.csv`, `te_pca_clr_matrix.csv` |
| Final ordination | `te_pca_scores.csv`, `te_pca_variance.csv`, `te_pca_loadings.csv` |
| Formal result audit | `te_pca_clustering_diagnostics.csv`, `te_pca_clustering_stability.csv`, `te_pca_candidate_cluster_assignments.csv`, `te_pca_trait_association_tests.csv`, `te_pca_phylogenetic_signal_tests.csv` |
| Reproducibility record | `te_pca_clustering_analysis_manifest.json` |

The candidate-assignment table is rejected audit evidence, not a set of final
biological cluster labels. There is no second workbook or output directory that
competes with these CSVs.

## Canonical products

- `te_order_composition.csv` and `te_superfamily_composition.csv`: closed,
  classified-conditional composition matrices after equal-weight averaging of
  the two *D. orestes* technical runs.
- `te_diversity.csv`: natural-log Shannon entropy and explicitly named
  Gini-Simpson index at order and superfamily levels under both the
  classified-conditional and mass-aware unresolved-bin views. Observed
  richness is retained only as an audit/support count.
- `te_feature_prevalence.csv`, `te_pca_clr_matrix.csv`, `te_pca_scores.csv`,
  `te_pca_variance.csv`, and `te_pca_loadings.csv`: the shared-feature
  compositional ordination audit trail.
- `te_pca_clustering_diagnostics.csv`, `te_pca_clustering_stability.csv`, and
  `te_pca_candidate_cluster_assignments.csv`: formal audit of discrete
  structure in the final PCA. Candidate assignments are explicitly rejected;
  no publication cluster labels were selected.
- `te_pca_trait_association_tests.csv` and
  `te_pca_phylogenetic_signal_tests.csv`: post-hoc organismal-trait tests and
  provisional tree-based signal tests on the full compositional geometry.
- `te_pca_clustering_analysis_manifest.json`: input hashes, seeds, software,
  acceptance gates, and interpretation limits for the structure audit.
- `repeatmasker_divergence_landscape.csv` and
  `repeatmasker_species_inventory.csv`: compact summaries of the archived
  hit-level RepeatMasker evidence.
- `ltr_*`: element metrics, species robustness, source coverage, and explicit
  exclusions for the LTR terminal:internal proxy. Their compact active
  provenance is `provenance/ltr_terminal_internal_release.json`.
- `assembly_quality.csv`, `dnapipete_input_quality.csv`, and
  `dnapipete_mass_accounting.csv`: quality and mass-accounting evidence.

## What is editable here

Open `explore.ipynb` to filter species, inspect features, change descriptive
summaries, or prototype plots. The canonical publication tables remain
unchanged until deliberately replaced and `make publication` is run.

The PCA uses closure, features positive in every included TE34 species, CLR
transformation, column centering, and SVD without post-CLR variance scaling.
Axes are oriented deterministically so the largest absolute loading is
positive. The classified-conditional order-level stratum is the primary
paper-facing diversity analysis. Superfamily-level and mass-aware
unresolved-bin estimates are sensitivity/audit strata. Rare features remain in
all diversity summaries even when they are not eligible for the shared-feature
PCA.

The two compact composition matrices reproduce the full released prevalence,
CLR, scores, variance, and loadings trail without the archived hit-level/HPC
data. The complete upstream genomic workflow remains frozen evidence, as
described below.

## Final PCA structure audit

The primary clustering representation uses all 23 nonzero, unscaled released
PC axes. This exactly preserves Euclidean distances in the 24-superfamily CLR
matrix, so it retains the full Aitchison geometry rather than selecting groups
from the visually cleaner PC1–PC2 projection.

The audit fits k-means for `k=1…8` with 1,000 starts. The primary gap rule is
the Tibshirani one-standard-error selector with squared Euclidean dispersion,
1,000 scaled-PCA null datasets, and seed `20260810`. Conservative release gates
also require adequate silhouette, recurring cluster sizes, leave-one-feature-
out stability, 1,000 fixed 28-of-34 species subsamples, and agreement with the
first-six-PC sensitivity. These audit gates were adopted during the rebuild;
they were not preregistered.

No discrete solution passes. The gap statistic selects `k=1`. Unconstrained
silhouette maximization selects `k=7`, but its mean silhouette is only 0.228
and two clusters are single species. The sole singleton-free candidate is
`k=2` (9/25 species), but its mean silhouette is 0.179 and its median
species-subsample adjusted Rand index is -0.042. Both candidates remain in the
audit table so the rejected result is visible rather than silently deleted.

Continuous composition nevertheless contains structure. Multivariate
phylogenetic signal is `K_mult=0.6581` with a whole-vector tip-permutation
`p=0.00001`; because `K_mult < 1`, the signal is weaker than Brownian-motion
expectation even though it exceeds randomized tip assignments. The expert
three-level reproductive strategy explains 26.2% of Aitchison variation by
species-label PERMANOVA (`p=0.00001`) without a detected dispersion difference.
That trait hypothesis is post hoc and the ordinary permutations are not
phylogenetically corrected. All tree results remain provisional until the
46-tip collaborator tree's citation and calibration provenance are supplied.

The formal audit is included in `make te-pca`. Its manifest binds the compact
inputs, PCA outputs, both implementation scripts, Python lock, R environment,
organismal traits, and focal tree by SHA-256. Canonical result replacement is a
deliberate maintenance operation, not part of normal exploration.

## LTR terminal:internal release

Start with [`provenance/ltr_terminal_internal_release.json`](provenance/ltr_terminal_internal_release.json).
It hash-binds the four canonical LTR tables, the archived corrected-release
builder and manifest, the recovered mapping facts, unresolved provenance, the
release counts and gates, the two historical fail-closed exclusions, the
TEsorter domain-count correction, and the original IQR method evidence. Rebuild
the local release from S14 or verify it read-only with:

```bash
uv run python analyses/01_transposable_elements/recompute_ltr_release.py
uv run pytest -q tests/test_ltr_release.py
```

The archived parser selected 1,088 rows because it split the TEsorter `Domains`
field on `|`, which is the delimiter *inside* each `DOMAIN|MODEL` annotation.
Consequently, four real annotations were stored as five and five were stored as
six. Two of those historical rows have corrupt partial depth files; S14 retains
the other 1,086 usable rows as a non-destructive audit substrate. Correctly
counting whitespace-delimited annotations yields 409 five-domain rows. One is
classified as Copia, leaving 408 LTR/Gypsy candidates at least 3,000 bp long
with at least five domains across 30 species. The remaining 678 usable rows stay
visible as audit-only evidence rather than being deleted.

The original analysis used an inclusive, within-species, two-sided Tukey
filter: retain ratios from `Q1 - 1.5 × IQR` through `Q3 + 1.5 × IQR`. Reapplied
within the corrected 408-element cohort, it retains 381 elements and flags 27
(26 upper and one lower outlier). Independently, 407 eligible elements pass the
0.8 regional positive-depth gate. Their intersection is the 380-element primary
Figure 4 branch and retains all 30 species.

The 0.8 gate requires the fraction of positions with positive depth to be at
least 0.8 separately in the left LTR, right LTR, and internal region. It is not
a gate on the fraction of positions reported in a depth file; all 1,086 usable
files report every expected position. The IQR fences are calculated first over
the 408 correctly selected ratios; the 80% gate did not exist in the original
notebook and is therefore recorded as a separate intersection.

The recoverable mapping workflow used paired, trimmed,
mitochondrial-filtered reads against per-species LTR contigs; Bowtie2 used
`-L 20 --very-sensitive-local`, followed by SAMtools BAM conversion, sorting,
indexing, and depth. The read-preprocessing commands and versions, Bowtie2 and
SAMtools versions, MAPQ threshold, multimapper policy, secondary/supplementary
alignment policy, and duplicate policy remain unresolved. These details must
not be inferred from the preserved outputs.

S14 records the corrected and historical domain counts, eligibility reason,
quartiles, fences, retention decision, exclusion direction, and primary-branch
membership for every usable historical row. Thus neither the 678 ineligible
audit rows nor the 27 IQR outliers are destructively erased. S15's declared
branches all start from the corrected 408-element cohort, while
`primary_iqr_filtered_coverage_ge_80pct` is the paper's primary LTR branch. A
later July 2026 compound one-element screen is retained only as a provenance
warning: it was not the original IQR method and must not be substituted for it.

## Frozen and archive-only upstream work

dnaPipeTE, RepeatModeler/custom-library construction, RepeatMasker, DANTE/LTR
annotation, TEsorter, Bowtie2, SAMtools, and depth generation are not rerun by
this paper repository. Their HPC products—including all `.align` files and
depth tables—are checksum-preserved in the pre-cleanup archive. The active LTR
release record binds the archived corrected-release builder and manifest by
exact path and SHA-256 without claiming that missing upstream commands have
been recovered. The historical custom `dedupe_telib.fasta` is absent, so a
full upstream regeneration cannot be claimed.

The multi-gigabyte hit-level tables are also archive-only. The compact active
tables conserve the reported summaries but do not pretend to replace primary
hit-level evidence for a future independent reanalysis.

## Interpretation boundary

The terminal:internal LTR statistic is a sensitivity-only
mapping/deletion-footprint proxy. Do not call it a direct ectopic-recombination
rate, solo-LTR rate, or DNA-loss rate. The planiceps-to-fuscus correction is
restricted to `SRX20497025 / GCA_032353935.1` and does not define general
synonymy; see the [shared identity-resolution contract](../../docs/IDENTITY_RESOLUTION.md).

Run `make te-pca-view` for interactive inspection and `make validate` for the
repository-wide hash and panel checks.
