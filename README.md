# Desmognathus transposable elements, morphology, and phylogenetic path analysis

This is the compact working repository for finishing the paper. It contains the
analysis-ready datasets, the decisions needed to interpret them, one shareable
R Markdown data report, the current R analyses, and the controlled publication package. The
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

# Verify the publication-facing TE34 Bonsai representation and validation
make te-bonsai

# Reconstruct the retained supplementary PCA geometry audit
make te-pca

# See the compact data inventory
make status

# Render the shareable data, code, methods, and ggplot report
make report

# Build the single indexed workbook of exact ggplot figure inputs
make figure-data

# Verify the author-edited Methods Word draft matches its Markdown companion
make methods-word
```

The active scientific contract is:

- **TE34** — 34 species with vetted accession-linked genomic resources.
- **Path24** — 24 species with finalized cell, nucleus, and conditional
  genome-size estimates derived from nuclear IOD.
- **overlap21** — the exact intersection used for the descriptive TE–genome-size
  comparison. The three Path24-only species are *D. brimleyorum*, *D.
  folkertsi*, and *D. ochrophaeus*; their TE values are not imputed.

## Repository map

```text
analyses/
  01_transposable_elements/   TE34 summaries, Bonsai representation, LTR proxy, and PCA audit
  02_morphology/              reviewed objects, IOD calibration, genome size, and provenance
  03_phylogenetic_path/       Path24 traits, trees, sensitivities, and R model
data/
  identity/                   shared panels, taxonomy, accessions, and sources
  DATA_CATALOG.csv            hash and shape inventory of every canonical file
Publication/                  data report, figure-input workbook, supplemental tables, figures, and manuscript aids
scripts/
  project.py                  catalog, identity, publication, status, validation
  publication/                GBE ggplot figure code
tests/                        compact repository contract tests
```

`Publication/Desmognathus_study_data.Rmd` is the researcher-facing data and
methods inspection surface. It reads the compact controlled release, exposes
the filters and equations, links every supplemental CSV, and embeds the ggplot
figures without becoming a second data authority. The companion
`Publication/Desmognathus_figure_review.Rmd` is a manifest-driven gallery for
critiquing every current figure, legend, alt text, and final-size layout.
`Publication/Desmognathus_figure_data.xlsx` is the single indexed workbook of
the 20 non-tree CSV tables directly consumed by those figures. Tree inputs stay
under `Publication/datasets/` in Newick and reviewer-readable CSV formats rather
than being forced into spreadsheet tabs. The workbook is a generated
convenience view; the manifest-bound release files remain authoritative.

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

# Verify the publication-facing Bonsai tree and strategy-localization artifacts
make te-bonsai

# Reconstruct the complete retained CLR/PCA audit trail
make te-pca

# Render the shareable study-data report
make report

# Render the visual critique gallery with every current figure
make figure-review

# Build and verify the labeled multi-sheet figure-data workbook
make figure-data

# Verify the researcher-facing Methods DOCX and Markdown remain synchronized
make methods-word

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

- Conditional *D. fuscus* assembly-anchored genome size in pg/1C is the default
  Path24 analysis trait. It is derived from nuclear IOD using the selected
  16.1-Gbp assembly anchor, not measured by an independent genome-size assay or
  validated as an absolute C-value. S21–S25 retain the underlying IOD
  measurement, image-quality, and calibration provenance.
- The LTR terminal:internal statistic is a mapping/deletion-footprint proxy. It
  is not a direct rate of ectopic recombination, solo-LTR formation, or DNA loss.
- The Path24 models are exploratory phylogenetic association and sensitivity
  analyses. They do not uniquely establish causal direction.
- The implemented Path24 fit currently contains conditional genome size,
  nucleus area, and cell area only. Positive IOD-to-pg rescaling leaves its
  standardized fit unchanged. A potential final model integrating TE
  composition and SVL has not yet been designed or fitted.
- The two *D. orestes* dnaPipeTE runs are technical retries and are averaged
  with equal weight before species-level summaries.
- Canonical TE diversity uses each species' ten most abundant positive
  superfamilies, reclosed to one before natural-log Shannon and Gini-Simpson
  calculation. Bonsai and the retained PCA audit instead use the same 24
  superfamilies detected in all 34 species and the same unscaled CLR geometry;
  these feature rules are intentionally separate from diversity.
- The three-level `reproductive_strategy` field in the curated SVL/organismal
  table is expert-vetted from Alex Pyron's interpretation of the initial PCA.
  It is a post hoc descriptive trait, not an a priori test of PCA clustering,
  and it is not interchangeable with adult aquaticity.
- The publication-facing Bonsai tree is a representation of TE-composition
  similarity, not a species phylogeny. Its fixed fitted tree strongly localizes
  post-hoc reproductive-strategy labels, but the exact aquatic-egg and
  direct-development edges recur in only 27.5% and 35.0% of 200 bounded
  feature-resampling refits. PCA remains a supplementary geometry and rejected-
  clustering audit; neither method establishes causal or phylogenetically
  independent trait association.
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

`Publication/datasets/DATASET_MANIFEST.csv` is the authoritative S01–S43 index.
`make publication` copies one canonical source into each supplement and updates
its SHA-256, dimensions, column inventory, tree inventory, and release summary.
`Publication/methods/METHODS_BULLET_SCAFFOLD.md` is the internal factual Methods
authority. The author-edited `MATERIALS_AND_METHODS_RESEARCHER_DRAFT.docx` and
its synchronized Markdown companion provide the collaborator-facing narrative
without audit tags, internal paths, or action-item placeholders. The verification
target does not overwrite Word edits. Remaining author inputs are
tracked separately and must be resolved before submission.

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
