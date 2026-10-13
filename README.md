# Desmognathus transposable elements, morphology, and phylogenetic path analysis

This repository contains the analysis-ready datasets, local analysis code,
Methods, figures, and publication tables for the Desmognathus study. Large HPC
products and the earlier working tree are preserved in a separate archive.

## Start here

Install `uv` and Conda, then run these commands from a fresh checkout. Figure
rendering also requires Arial; the environment check verifies that R can find
the font used by the retained figures.

```bash
make setup          # Install locked Python dependencies and create a missing R environment
make check-env      # Check the installed environments without changing packages
make test           # Validate the release and run the tests
make identity       # Inspect taxonomy, accession exceptions, and tree representatives
make te-bonsai      # Check retained TE34 Bonsai artifacts
make methods-word   # Check the Methods Word draft against its Markdown companion
make status         # Inspect the data inventory
```

Python is pinned to 3.13.11 in `.python-version`; `uv.lock` pins its packages.
The R environment defaults to `desmognathus`. On Linux x86-64, setup uses
`environment-linux-64.explicit.txt`, which pins the Conda package builds and
their SHA-256 hashes. Other platforms resolve `environment.yml`. Setup installs
checksum-verified CRAN `phylopath` 1.3.1 only when creating an environment.
An existing R environment is validated and left unchanged.

To preserve an older environment, choose separate destinations before setup:

```bash
export UV_PROJECT_ENVIRONMENT=/absolute/path/to/new/python-environment
export DESMOGNATHUS_R_ENV=/absolute/path/to/new/r-environment
make setup
make check-env
```

`DESMOGNATHUS_R_ENV` also accepts a Conda environment name. Routine Make targets
use the installed Python environment without syncing packages. If a check
fails, read its diagnostic before running setup or changing environments.

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

## Local checks and regeneration

`make validate`, `make identity`, `make status`, `make methods-word`, and
`make te-bonsai` inspect retained artifacts. `make lint` checks Python source.
`make test` runs validation and the test suite; test-generated artifacts use
temporary directories.

The following commands write analysis or publication artifacts. Run them in a
separate checkout when comparing a reconstruction with the retained release.

| Command | Artifacts it creates or updates |
| --- | --- |
| `make publication` | Supplementary exports, their manifest, and the data catalog |
| `make catalog` | The canonical data catalog |
| `make te-pca` | The reconstructed CLR/PCA audit |
| `make path-quick` / `make path` | Path24 release inputs and model output |
| `make figures` | GBE figure files |
| `make report` | Figures and the study-data HTML report |
| `make figure-review` | Figures and the figure-review HTML gallery |
| `make figure-data` | The figure-input workbook |

The retained analyses start from frozen tables and upstream outputs. They do
not rerun dnaPipeTE, RepeatModeler, RepeatMasker, segmentation, or image review.
See [docs/REPRODUCIBILITY.md](docs/REPRODUCIBILITY.md) for the available inputs
and remaining gaps. Historical `Dusky` and figure-only environments are not
required by the local entrypoints.

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
- The original Path24 path suite contains conditional genome size, nucleus
  area, and cell area. Later analyses include maximum SVL, an eight-class
  body-role comparison, and TE-adjusted regressions in overlap21. These remain
  exploratory comparisons, without a selected all-data causal network.
  Positive IOD-to-pg rescaling leaves standardized cellular fits unchanged.
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
