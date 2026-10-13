# Desmognathus Genome Biology and Evolution publication package

This is the controlled release surface for the planned *Genome Biology and Evolution*
article. It is built from the compact canonical data under `data/`
and `analyses/`, keeping paper outputs separate from frozen upstream evidence
and archive-only computation.

## Layout

- `datasets/` — S01–S43 CSV supplements, two phylogenetic tree files, one
  non-phylogenetic Bonsai composition-tree representation, column inventory,
  SHA-256 manifest, and release summary.
- `figures/` — seven numbered ggplot masters plus five unnumbered
  review references: two Wang-style candidate layouts, the Path24 tree aligned
  to development mode and phenotype estimates, the Path24 equivalence-class
  DAG gallery, and a descriptive Wang et al. 2025 cross-study TE-diversity reference. The numbered and
  auxiliary records are separated in `FIGURE_MANIFEST.csv` and
  `AUXILIARY_FIGURE_MANIFEST.csv`.
- `methods/` — internal factual scaffold and evidence files plus one clean,
  author-edited researcher-facing Methods DOCX and synchronized Markdown mirror.
- `journal_requirements/` — GBE/OUP requirements, with
  `GBE_REQUIREMENTS.md` as the official requirements authority; recheck the
  live GBE instructions immediately before submission.
- `Desmognathus_study_data.Rmd` and `.html` — shareable data, filter, equation,
  frozen-statistics, and ggplot inspection surface.
- `Desmognathus_figure_review.Rmd` and `.html` — manifest-driven visual gallery
  with all seven numbered and all five auxiliary figures, legends, alt text, final
  dimensions, and critique checklist.
- `Desmognathus_figure_data.xlsx` — indexed multi-sheet workbook containing the
  16 exact non-tree CSV tables consumed by the numbered and auxiliary ggplot figures. Tree
  inputs remain under `datasets/` in Newick and reviewer-readable CSV formats;
  the auxiliary development track additionally reads the canonical
  `analyses/03_phylogenetic_path/data/organismal_traits.csv`. The workbook is a
  generated view, not a new authority.

## Rebuild and validate

From the repository root:

```bash
# Canonical analysis tables -> S01-S43, tree copies, manifests, and catalog
make publication

# Seven numbered figures plus five auxiliary review figures, each in
# PDF/PNG/TIFF
make figures

# Self-contained researcher-facing HTML report
make report

# Self-contained browser gallery for figure-by-figure visual review
make figure-review

# Single labeled workbook of all direct figure inputs
make figure-data

# Verify the author-edited Materials and Methods Word draft
make methods-word

# Dataset identity, dimensions, hashes, panels, taxonomy decision, and tree tips
make test
```

Both figure production and the Path24 R analysis use the single `desmognathus`
environment defined in `environment.yml`. Python uses `pyproject.toml` and
`uv.lock`. There is no separate publication-figure environment.

## Figure style and review

Use ggplot2 wherever possible with the shared theme and Arial throughout; new
builds must fail if Arial is unavailable rather than silently use another font.
The 88-mm single-column and 180-mm double-column widths are project house-style
choices, not GBE mandates. Base text is 8 pt with a 7-pt minimum target at final
size; GBE's hard minimum is 6 pt. Preserve the shared palette and canonical
species ordering documented in `figures/GBE_FIGURE_PLAN.md`.

New exports comprise vector PDF masters with embedded Arial, review PNGs, and
300-ppi CMYK TIFF revision derivatives for color figures. Existing exports are
pending one-at-a-time re-rendering and review under this contract. Work through
each figure with the author, agree its reference layout and purpose, then
review palette, species order, axes and units, typography, and export
properties at final dimensions before proceeding to the next figure.

## Release rule

`main_candidate`, `supplementary`, `sensitivity_only`, and `audit_only` are not
interchangeable. `datasets/DATASET_MANIFEST.csv` records the class for every
supplement. Conditional *D. fuscus* assembly-anchored genome size is the default
analysis trait in S04, S26-S37, and S40; S38-S39 are phylogeny tables. The trait
is IOD-derived, not an independently
measured genome-size assay or validated absolute C-value; S21-S25 preserve its
IOD measurement and calibration provenance. The LTR terminal:internal statistic
is not a measured ectopic-recombination rate, and the Path24 graph comparisons
do not establish a unique causal direction.

S08 is the canonical TE-diversity table. Each species contributes natural-log
Shannon and Gini-Simpson estimates calculated after its ten most abundant
positive superfamilies are selected and explicitly reclosed to one. This
species-specific top-ten rule is distinct from the shared-feature rule used by
the Bonsai representation and retained PCA audit, which retain the 24
superfamilies detected in all 34 species. Figure 2 uses S42 tip annotations,
S43 validation statistics, and the manifest-bound Bonsai Newick file; S09–S11
remain supplementary PCA geometry outputs.

The DAG gallery labels its first node **Genome size**, meaning the conditional
*D. fuscus* assembly-anchored pg/1C sensitivity derived from IOD. Because the
Path24 fits log-transform and standardize the trait, the saved fit is invariant
to this positive scale conversion. The label does not turn the estimate into an
independent genome-size assay and does not identify a unique causal arrow. IOD
and genome size are not separate DAG nodes; the former is the measurement basis
for the latter.

`methods/METHODS_BULLET_SCAFFOLD.md` is the internal factual Methods authority.
`methods/MATERIALS_AND_METHODS_RESEARCHER_DRAFT.docx` is the single shareable,
author-edited narrative. The neighboring Markdown file mirrors its visible prose,
and `make methods-word` verifies synchronization without replacing the Word file.
It contains no internal audit labels or action-item placeholders.
Complete the unresolved rows in `methods/METHODS_REQUIRED_AUTHOR_INPUT.csv`
before treating the draft as submission-final prose.
