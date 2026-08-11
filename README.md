# Desmognathus transposable elements, morphology, and phylogenetic path analysis

This is the compact working repository for finishing the paper. It contains the
analysis-ready datasets, the decisions needed to interpret them, three editable
notebooks, the current R analyses, and the controlled publication package. The
large HPC products and the complete pre-cleanup working tree are preserved in a
separate archival checkout rather than mixed into the paper workspace.

## Start here

```bash
# Install the Python environment and the single R environment
make setup

# Check the scientific panel contract and every publication-table hash
make test

# Inspect every source alias, exact tree representative, and accession exception
make identity

# Reconstruct and formally verify the final TE34 PCA
make te-pca

# See the compact data inventory
make status

# Open one notebook per analysis component
make notebooks
```

The active scientific contract is:

- **TE34** — 34 species with vetted accession-linked genomic resources.
- **Path24** — 24 species with finalized cell, nucleus, and relative nuclear-IOD
  measurements.
- **overlap21** — the exact intersection used for the descriptive TE–IOD
  comparison. The three Path24-only species are *D. brimleyorum*, *D.
  folkertsi*, and *D. ochrophaeus*; their TE values are not imputed.

## Repository map

```text
analyses/
  01_transposable_elements/   TE34 summaries, LTR proxy, and TE notebook
  02_morphology/              reviewed objects, relative IOD, and provenance
  03_phylogenetic_path/       Path24 traits, trees, sensitivities, and R model
data/
  identity/                   shared panels, taxonomy, accessions, and sources
  DATA_CATALOG.csv            hash and shape inventory of every canonical file
Publication/                  supplemental tables, figures, and manuscript aids
scripts/
  project.py                  catalog, identity, publication, status, validation
  build_notebooks.py          regenerates the three small notebooks
  publication/                GBE ggplot figure code
tests/                        compact repository contract tests
```

Each analysis directory has its own README and one `explore.ipynb`. The
notebooks intentionally read only the small canonical tables beside them, so
they remain useful after the HPC-era directory tree is removed.

## What can be rerun here

```bash
# Rebuild Publication/datasets from the canonical analysis tables
make publication

# Rebuild the catalog after a deliberate data edit, then validate the release
make catalog
make validate

# Show the identity decisions and verify them directly against raw sources
make identity

# Refit the Path24 analysis quickly, writing disposable output under analyses/.../output
make path-quick

# Reconstruct the complete CLR/PCA audit trail and replay the formal structure audit
make te-pca

# Open the focused interactive PCA/data notebook
make te-pca-view

# Rebuild the current GBE ggplot figures
make figures
```

Python is managed by `uv` from `pyproject.toml` and `uv.lock`. R and the system
figure converters are defined in the single `desmognathus` Conda environment
from `environment.yml`. Historical `Dusky` and figure-only environments are no
longer repository dependencies. Because `phylopath` is not packaged by the
configured Conda channels, `scripts/setup_r.sh` installs its checksum-pinned
CRAN 1.3.1 source after Conda resolves the rest of the environment.

## Evidence boundaries

- Relative nuclear IOD is an image-derived relative phenotype. It is not an
  independently validated absolute genome size or C-value.
- The LTR terminal:internal statistic is a mapping/deletion-footprint proxy. It
  is not a direct rate of ectopic recombination, solo-LTR formation, or DNA loss.
- The Path24 models are exploratory phylogenetic association and sensitivity
  analyses. They do not uniquely establish causal direction.
- The implemented Path24 fit currently contains relative nuclear IOD, nucleus
  area, and cell area only. A potential final model integrating TE composition
  and SVL has not yet been designed or fitted.
- The two *D. orestes* dnaPipeTE runs are technical retries and are averaged
  with equal weight before species-level summaries.
- The three-level `reproductive_strategy` field in the curated SVL/organismal
  table is expert-vetted from Alex Pyron's interpretation of the initial PCA.
  It is a post hoc descriptive trait, not an a priori test of PCA clustering,
  and it is not interchangeable with adult aquaticity.
- The final 24-superfamily CLR composition has significant continuous
  phylogenetic and post-hoc reproductive-strategy structure, but no k-means
  solution passes the gap, silhouette, cluster-size, and resampling gates.
  Candidate labels are retained only as rejected audit evidence; there are no
  accepted biological PCA clusters.
- Only the public resource `SRX20497025 / GCA_032353935.1`, labelled *D.
  planiceps* by its source, is treated as *D. fuscus*. This accession-specific
  decision, the independent SVL aliases, and exact tree-tip selection are
  defined in [docs/IDENTITY_RESOLUTION.md](docs/IDENTITY_RESOLUTION.md).
- *D. aeneus*, *D. orestes*, and *D. wrighti* remain in Path24 after manual
  review despite falling below an automatic target in at least one measurement
  stream.

See [docs/REPRODUCIBILITY.md](docs/REPRODUCIBILITY.md) for the line between
locally rerunnable, frozen, and archive-only work, and
[docs/ARCHIVE.md](docs/ARCHIVE.md) for recovery instructions. The exact
replacement and deletion audit is in
[docs/CLEANUP_LEDGER.md](docs/CLEANUP_LEDGER.md).

## Publication surface

`Publication/datasets/DATASET_MANIFEST.csv` is the authoritative S01–S40 index.
`make publication` copies one canonical source into each supplement and updates
its SHA-256, dimensions, column inventory, tree inventory, and release summary.
`Publication/methods/METHODS_BULLET_SCAFFOLD.md` is the sole active factual
Methods authority. Submission prose will be created from that scaffold after
the remaining author inputs are resolved; the active repository does not carry
a stale parallel manuscript draft.

The focal 46-tip collaborator tree is preserved, but its final source and
calibration citation is still unresolved. That is a publication blocker, not a
reason to invent provenance.

## Pre-cleanup archive

The exact pre-cleanup branch, Git history, dirty working state, ignored HPC
outputs, intermediates, and legacy analyses are preserved at:

```text
/home/jake/Projects/Desmognathus_TE_archive
```

That checkout is the recovery and provenance layer. Do not use it as the active
paper workspace. Its `PRE_CLEANUP_ARCHIVE.md` and `ARCHIVE_MANIFEST.csv` record
the snapshot and per-file checksums.
