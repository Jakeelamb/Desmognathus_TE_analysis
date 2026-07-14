# Desmognathus GBE publication package

This directory is the controlled release surface for the planned *Genome
Biology and Evolution* Article. It is generated from frozen, accession-explicit
artifacts and keeps manuscript-ready evidence separate from sensitivity and
audit-only products.

## Layout

- `journal_requirements/` — current first-party GBE/OUP requirements and the
  project figure house specification.
- `methods/` — Methods source matrix, required author-supplied wet-lab details,
  Markdown source, and generated DOCX.
- `datasets/` — plain-CSV supplemental tables where structurally possible,
  standard phylogenetic files where CSV would lose scientific structure, and a
  SHA-256 release manifest.
- `figures/` — GBE-styled ggplot figure masters, legends, and alt text.

## Release rule

`main_candidate`, `supplementary`, `sensitivity_only`, and `audit_only` are not
interchangeable. The dataset manifest carries one of these labels for every
file. In particular, the current image assay supports a relative nuclear-IOD
sensitivity index but does not yet establish an absolute genome size in pg.
The terminal:internal LTR statistic is an exploratory deletion-footprint proxy,
not a measured ectopic-recombination rate.

## Rebuild

```bash
# Frozen supplemental CSVs and tree views
uv run --with biopython scripts/publication/build_publication_datasets.py

# One-time R/ggplot environment creation
conda env create -f Publication/figure_environment.yml

# Seven GBE-styled ggplot figures and all derivatives
scripts/publication/run_gbe_figure_build.sh

# Double-spaced, line-numbered DOCX with embedded figures
uv run --with python-docx python scripts/publication/build_methods_docx.py

# Release validation
uv run python -m unittest \
  scripts.python.tests.test_publication_package \
  scripts.python.tests.test_publication_outputs
```

The DOCX is a working Methods draft, not a submission-ready assertion that all
methods are known. Complete every row in
`methods/METHODS_REQUIRED_AUTHOR_INPUT.csv` before manuscript freeze.
