# Desmognathus transposable elements, morphology, and phylogenetic path analysis

This repository contains the analysis-ready datasets, local analysis code,
Methods, figures, and publication tables for the Desmognathus study. Large HPC
products and the earlier working tree are preserved in a separate archive.

## Start here

Install Pixi 0.77.1 or newer, then run these commands from the repository root.
The locked environments support Linux x86-64 with glibc 2.39 or newer. Figure
rendering requires Arial; the environment check verifies that R can find it.

```bash
pixi install --locked --all
pixi run --clean-env --locked --no-install setup
pixi run --clean-env --locked --no-install test
pixi run --clean-env --locked --no-install identity
pixi run --clean-env --locked --no-install te-bonsai
pixi run --clean-env --locked --no-install methods-word
pixi run --clean-env --locked --no-install status
```

`pixi.toml` defines the local tasks and dependencies. `pixi.lock` pins all
environments. The default environment contains Python 3.13.11, R 4.4.3, the
analysis packages, Node, GNU utilities, Pandoc, and figure converters. The
Bonsai environment uses Python 3.14 and the pinned upstream requirements, which
differ from the main analysis packages. The `image-qc` environment adds pinned
image-reading and quality-audit packages to the default Python dependencies.
All three environments share the retained R package builds.

`setup` installs checksum-verified CRAN `phylopath` 1.3.1 into the private
`.pixi/r-library` directory. It leaves existing UV and Conda environments
unchanged. Repeated setup validates the completed installation. A damaged or
incomplete private library fails validation without being overwritten.

```bash
make setup       # Install all locked Pixi environments and prepare phylopath
make check-env   # Check the installed default environment without installing packages
make test        # Run the same Pixi test task
```

Make targets delegate to the tasks in `pixi.toml`; they do not maintain separate
command definitions. Routine Make targets use `--locked --no-install` and
`--clean-env` and refuse missing environments. Include `--clean-env` in direct
Pixi runs to exclude inherited Python and R settings. A plain
`pixi run --locked <task>` may install locked packages before running the task.
Use `make check-env` or
`pixi run --clean-env --locked --no-install --offline check-env` for a check that
cannot install or repair packages.

The older `pyproject.toml`, `uv.lock`, `environment.yml`, and Linux explicit
snapshot remain unchanged as provenance. Some are hashed inputs to retained
scientific manifests. They are not the active setup. The local workflow rules
are recorded in [AGENTS.md](AGENTS.md).

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
the non-tree CSV tables directly consumed by those figures. Tree inputs stay
under `Publication/datasets/` in Newick and reviewer-readable CSV formats rather
than being forced into spreadsheet tabs. The workbook is a generated
convenience view; the manifest-bound release files remain authoritative.

## Local checks and regeneration

`make validate`, `make identity`, `make status`, `make methods-word`, and
`make te-bonsai` inspect retained artifacts. `make lint` checks Python source.
`make test` runs validation and the test suite; test-generated artifacts use
temporary directories. `make te-pca` recomputes and checks the retained PCA
audit without promoting data; its permutation and resampling checks are
substantial and are separate from the quick checks above.

The following commands write analysis or publication artifacts. Run them in a
separate checkout when comparing a reconstruction with the retained release.

| Command | Artifacts it creates or updates |
| --- | --- |
| `make publication` | Supplementary exports, their manifest, and the data catalog |
| `make catalog` | The canonical data catalog |
| `make path-quick` / `make path` | Path24 model output after checking release inputs |
| `make figures` | GBE figure files |
| `make report` | Figures and the study-data HTML report |
| `make figure-review` | Figures and the figure-review HTML gallery |
| `make figure-data` | The figure-input workbook |

For new Bonsai fits, supply a new run directory. These tasks use the separate
locked Bonsai environment and fetch the exact upstream Git revision into a
private source cache. The frozen historical shell runners remain as evidence.

```bash
pixi run --clean-env --locked --no-install -e bonsai bonsai-rerun --run-dir /absolute/path/to/new-fit
pixi run --clean-env --locked --no-install -e bonsai bonsai-bootstrap-rerun --run-dir /absolute/path/to/new-bootstrap
```

The bootstrap defaults to the retained 200-replicate campaign. A bounded run
with `--replicates 2 --jobs 1` tests execution only and does not produce the
final campaign analysis. Resume requires `--resume` and a matching run manifest.
After a forced termination, a retained `.bonsai-run.lock` prevents concurrent
resume; verify that the run has stopped before removing that scratch lock.

The retained analyses start from frozen tables and upstream outputs. Local
setup does not rerun dnaPipeTE, RepeatModeler, RepeatMasker, segmentation, or
image review. See [docs/REPRODUCIBILITY.md](docs/REPRODUCIBILITY.md) for available
inputs and gaps. Historical `Dusky` and figure-only environments are not used.

The [HPC source audit](docs/HPC_SCRIPT_AUDIT.md) recovered 12 historical job
files from local Git history, plus command notes for preprocessing and TE
classification. Seven job files add to the five previously recovered sources.
Exact copies and checksums are retained as evidence. The custom TE library,
container images, and complete job-to-output records remain unresolved; these
historical files are not active Pixi tasks.

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

Local Methods audits are registered as `specimen-records` in the default
environment and `standard-qc` in `image-qc`. See the
[specimen audit](docs/submission_work_20261013/specimen_records/README.md) and
[standard-review workflow](docs/STANDARD_QC.md) for source requirements and
fresh output directories. The image audit requires the preserved microscopy
archive; the small review gallery is included in the repository.

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
