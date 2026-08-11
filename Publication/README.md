# Desmognathus GBE publication package

This is the controlled release surface for the planned *Genome Biology and
Evolution* article. It is built from the compact canonical data under `data/`
and `analyses/`, keeping paper outputs separate from frozen upstream evidence
and archive-only computation.

## Layout

- `datasets/` — S01–S40 CSV supplements, two phylogenetic tree files, column
  inventory, SHA-256 manifest, and release summary.
- `figures/` — current GBE-styled ggplot masters, legends, alt text, and figure
  manifest.
- `methods/` — factual Methods scaffold, evidence matrix, and missing
  author-input checklist.
- `journal_requirements/` — the current project copy of the GBE/OUP submission
  requirements.

## Rebuild and validate

From the repository root:

```bash
# Canonical analysis tables -> S01-S40, tree copies, manifests, and catalog
make publication

# Eight GBE ggplot figures and PDF/PNG/TIFF derivatives
make figures

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
authority. Complete every unresolved row in
`methods/METHODS_REQUIRED_AUTHOR_INPUT.csv` before manuscript prose and DOCX
production. Dataset rebuilds intentionally do not generate manuscript prose.
