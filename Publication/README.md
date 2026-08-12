# Desmognathus GBE publication package

This is the controlled release surface for the planned *Genome Biology and
Evolution* article. It is built from the compact canonical data under `data/`
and `analyses/`, keeping paper outputs separate from frozen upstream evidence
and archive-only computation.

## Layout

- `datasets/` — S01–S40 CSV supplements, two phylogenetic tree files, column
  inventory, SHA-256 manifest, and release summary.
- `figures/` — current GBE-styled ggplot masters, legends, alt text, and figure
  manifest, plus one unnumbered DAG reference asset
  (`Path_DAG_equivalence_class_reference.{pdf,png,tif}`) kept outside
  `FIGURE_MANIFEST.csv` for review only.
- `methods/` — factual Methods scaffold, evidence matrix, missing author-input
  checklist, and two synchronized Word working exports for review.
- `journal_requirements/` — the current project copy of the GBE/OUP submission
  requirements.
- `Desmognathus_study_data.Rmd` and `.html` — shareable data, filter, equation,
  frozen-statistics, and ggplot inspection surface.
- `Desmognathus_figure_review.Rmd` and `.html` — manifest-driven visual gallery
  with every figure, legend, alt text, final dimensions, and critique checklist.
- `Desmognathus_figure_data.xlsx` — indexed multi-sheet workbook containing the
  19 exact non-tree CSV tables directly consumed by the current ggplot figure
  builder. Tree inputs remain under `datasets/` in Newick and reviewer-readable
  CSV formats; the workbook is a generated view, not a new authority.

## Rebuild and validate

From the repository root:

```bash
# Canonical analysis tables -> S01-S40, tree copies, manifests, and catalog
make publication

# Eight GBE ggplot figures and PDF/PNG/TIFF derivatives
make figures

# Self-contained researcher-facing HTML report
make report

# Self-contained browser gallery for figure-by-figure visual review
make figure-review

# Single labeled workbook of all direct figure inputs
make figure-data

# Dataset identity, dimensions, hashes, panels, taxonomy decision, and tree tips
make test
```

Both figure production and the Path24 R analysis use the single `desmognathus`
environment defined in `environment.yml`. Python uses `pyproject.toml` and
`uv.lock`. There is no separate publication-figure environment.

## Release rule

`main_candidate`, `supplementary`, `sensitivity_only`, and `audit_only` are not
interchangeable. `datasets/DATASET_MANIFEST.csv` records the class for every
supplement. Relative nuclear IOD is not an independently validated absolute
genome size, the LTR terminal:internal statistic is not a measured
ectopic-recombination rate, and the Path24 graph comparisons do not establish a
unique causal direction.

`methods/METHODS_BULLET_SCAFFOLD.md` is the sole active factual Methods
authority. The Word-ready scaffold and narrative DOCX are reviewed derivatives,
not independent authorities or submission-final prose. Complete every
unresolved row in `methods/METHODS_REQUIRED_AUTHOR_INPUT.csv` before finalizing
the manuscript. Dataset rebuilds intentionally do not regenerate Word prose.
