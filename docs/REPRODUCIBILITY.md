# Reproducibility boundary

Local analyses start from retained upstream outputs. The sections below
distinguish runnable code, frozen evidence, and material held in the archive.

## 1. Locally rerunnable

- data catalog, hash validation, and S01–S43 publication exports;
- the R Markdown study-data report and manifest-driven figure-review gallery,
  each with a self-contained HTML rendering;
- the single labeled figure-data workbook, rebuilt from the manifest-bound CSV
  inputs other than the tree files with exact sheet/source parity checks; tree
  inputs remain in their standard Newick and reviewer-readable CSV formats;
- the single researcher-facing Materials and Methods DOCX, deterministically
  author-edited in Word and checked against its Markdown mirror with
  `make methods-word`, which does not overwrite the DOCX;
- the complete TE34 order/superfamily CLR-PCA reconstruction and formal final
  PCA audit from the preserved compact composition tables, including an exact
  replay of the 34-row, four-column canonical diversity table calculated after
  species-specific selection and reclosure of the ten most abundant positive
  superfamilies;
- compact morphology and genome-size release validation from the frozen object
  tables, including the underlying IOD estimator, species estimates, image
  aggregation, pixel scale, animal support, and the declared quality-balance
  gate;
- the arithmetic conditional *D. fuscus* reference conversion from the published
  pooled assembly (`16.1 Gbp / 0.978 Gbp per pg = 16.462167689 pg/1C`), which
  defines the default Path24 genome-size trait while retaining IOD as its
  measurement and calibration provenance;
- the corrected LTR/Gypsy five-domain gate, primary IQR/coverage branch, and
  species/resource summaries from the frozen 1,086-row zero-aware audit table;
- the Path24 R analysis from the preserved compact traits and trees;
- the current GBE ggplot figure bundle.
- the S41 cross-study diversity table from an exact hash-bound Wang et al. 2025
  Dryad snapshot plus current top-ten-superfamily recomputation; this is a
  descriptive non-harmonized comparison with no pooled inference.

### Pixi setup and checks

`pixi.toml` and `pixi.lock` are the active local environment and task definitions.
Install Pixi 0.77.1 or newer. The verified platform is Linux x86-64 with glibc
2.39 or newer. Other platforms are not declared. Arial is an external font
prerequisite checked by the project font resolver.

```bash
pixi install --locked --all
pixi run --clean-env --locked --no-install setup
make check-env
make test
```

The default environment contains Python 3.13.11, R 4.4.3, the paper-analysis
packages, Node, GNU utilities, Pandoc, and figure converters. The separate
`bonsai` environment contains Python 3.14 and the upstream requirements at the
pinned Git revision.
Its NumPy and pandas versions differ from the main environment, so the two
Python dependency sets are solved separately in the same lockfile. They share
the retained R package builds.

The `image-qc` environment uses the default Python, NumPy and pandas pins plus
Pillow, tifffile, roifile, scikit-image, SciPy and Matplotlib. It supports the
existing-mask standards audit without segmentation or changes to released
measurements. Run its installed-package check with
`pixi run -e image-qc --clean-env --locked --no-install check-env`.

The [Methods recovery records](submission_work_20261013/PLAN.md) distinguish
recovered specimen and model evidence from unresolved calibration,
camera and upstream-runtime records. Local audit tasks require fresh output
directories and preserve the released tables.

The [Methods audit verification](submission_work_20261013/VERIFICATION.md)
records the current three-environment checks, actual local task replays and
scientific preservation. The default suite passed 365 tests with 11 native
archive-dependent skips; all skipped cases passed in the separate 23-test
image-QC run. Existing default and Bonsai lock records were unchanged.

`setup` checksum-verifies CRAN `phylopath` 1.3.1 and installs it into a private
`.pixi/r-library` directory. Setup checks a staged library before publishing it.
Existing completed libraries are checked without reinstalling. Incomplete or
damaged libraries fail without being overwritten. Neither setup nor the R
wrapper selects a legacy UV or Conda environment.

Activation isolates Python imports and R startup files from user settings.
Use `pixi run --clean-env` when invoking tasks directly. A live rejection
control showed that inherited `PYTHONHOME` can survive Pixi's activation delta
and prevent Python startup; the clean environment removes it before startup.
Make targets include this flag.
Both `R_LIBS_USER` and `R_LIBS_SITE` point to the private R library. The site
setting preserves access to `phylopath` when the frozen path script resets its
library search path. R and Python executables must belong to the selected
project Pixi prefix.

`make check-env` guards the installed prefix and invokes Pixi with
`--locked --no-install --offline`. It compares actual installed Conda metadata
and Python distributions with the lock inventory, then checks R package pins,
Cairo, Arial, the private library, and environment-local converters. It cannot
create or repair an environment. A plain `pixi run --locked check-env` may
install packages before the check. Other Make targets also use
`--locked --no-install`; `make setup` is the explicit provisioning command.

The older `pyproject.toml`, `uv.lock`, `environment.yml`, and
`environment-linux-64.explicit.txt` are retained snapshots. Scientific manifests
hash several of them, so they remain byte-for-byte unchanged. The old local
`.venv` and Conda environments are also preserved. Local command policy is
recorded in [AGENTS.md](../AGENTS.md).

The [README](../README.md) distinguishes checks from artifact writers. Use a
complete scratch checkout for publication exports, figures, reports, workbooks,
path output, and new Bonsai fits. `make te-pca` checks the reconstructed PCA
audit without promoting its results. Its resampling settings remain those in
the frozen script.

Bonsai refits use `scripts/run_bonsai.py` in the `bonsai` environment. It fetches
and checks the exact upstream revision, or accepts an explicit clean pinned
checkout. Run directories must be new, outside the scientific directories, or
below `.pixi/runs`. Bootstrap resume requires a matching run manifest and an
exclusive run lock. The default is the retained 200-replicate campaign. Bounded
runs are execution smoke checks and omit final campaign analysis.

### Migration verification

The migration was exercised on 13 October 2026 in a fresh, complete export of
the tracked repository, with separate Pixi prefixes and a new private R
library. Both environments passed installed-package checks. All 261 historical
R-environment package URLs and SHA-256 hashes matched the default Pixi lock.
Repeated setup left all 38 files in the completed `phylopath` library unchanged.

| Workflow | Observed result |
| --- | --- |
| Environment and test gates | Both locked environments passed installed-package checks; the final scratch suite passed all 337 tests with no skips, and Ruff passed |
| Publication exports and catalog | S01–S43 validated; TE34, Path24, and overlap21 membership preserved |
| Figure-input workbook | 16 input tables and 19 sheets; exact retained workbook hash reproduced |
| Figure generation | Seven numbered figures and five auxiliary figures generated in scratch |
| Study report and figure-review gallery | Both fresh HTML renders passed source-fingerprint, manifest-order, narrative-boundary, and twelve-embedded-PNG checks |
| Read-only release checks | Identity, status, Methods Word parity, and retained Bonsai reconstruction checks passed |
| TE PCA | Complete reconstruction and formal permutation/resampling audit passed |
| Path24 | Default full analysis completed: 250 measurement bootstraps, 200 published bootstrap trees, 24 leave-one-species-out fits, and 800 simulations; zero failed fits |
| Bonsai | Default full refit and all 200 feature-bootstrap refits, plus matched calibration, completed with their analysis and plots; retained strategy-support fractions reproduced |

All artifact-writing commands ran in scratch directories. Retained datasets,
Methods, figures, HTML reports, workbooks, and analysis outputs were checked
against the pre-migration SHA-256 inventory rather than replaced by these runs.
Independent review confirmed 12,280 protected entries unchanged, with no missing
or added entries; the only two source exceptions were the report templates.
The two maintained R Markdown templates were corrected to accept both the
retained auxiliary-figure manifest and the current figure builder's manifest.
They now read the selected panel's narrative from that manifest and record
explicit source fingerprints. Historical rendered reports and their checksum
records remain unchanged.

## 2. Frozen analysis evidence

The active repo retains compact outputs that are needed to inspect, filter, and
report the analysis but whose expensive upstream computation is not honestly
reproducible here:

- TE34 composition, diversity, PCA, RepeatMasker landscape, and resource QC;
- element-level and robustness tables for the LTR terminal:internal proxy;
- finalized reviewed cell/nucleus objects and nuclear-IOD objects used to derive
  the conditional genome-size estimates;
- microscopy selection decisions and production-lineage snapshot;
- the collaborator source tree and published-tree sensitivity inputs.

Every canonical file is listed with a SHA-256 in `data/DATA_CATALOG.csv`.

## 3. Archive-only evidence and history

The adjacent pre-cleanup archive preserves raw/large HPC outputs, remote genome
caches, hit-level substrates, generated review galleries, intermediates,
retired analysis18 branches, old environments, and every pre-cleanup script and
document. A per-file checksum manifest distinguishes irreplaceable HPC evidence
from regenerable intermediates and historical worktree content.

The [13 October 2026 HPC source audit](HPC_SCRIPT_AUDIT.md) recovered 12
historical job files from local Git history. Five were already retained during
submission preparation; seven additional exact copies and a command notebook
are now preserved with [source identities and checksums](hpc_script_audit/manifest.json).
The notebook includes fastp, classifier, Dfam and BAM-to-FASTQ examples, with
transcription errors and settings that require reconciliation. Historical
sources are evidence of configured commands, not proof that each retained
output came from those commands, and are not maintained Pixi tasks.

The full genomic pipeline cannot currently be claimed reproducible. The
historical custom TE library `dedupe_telib.fasta` and container images were not
located in the audited local trees; complete HPC runtime and job-to-output
records remain unresolved. The archive preserves observed outputs. Recovery
does not regenerate or alter those outputs.

The canonical diversity table itself can be promoted from the frozen
superfamily composition matrix with
`pixi run --clean-env --locked --no-install --executable python analyses/01_transposable_elements/recompute_pca.py --write`.
`make te-pca` must then pass as the read-only reconstruction gate. Promotion
with `--write` belongs in a separate checkout when comparing reconstructions.

This upstream boundary does not prevent local PCA verification: `make te-pca`
starts at the frozen, species-level composition matrices and reconstructs the
complete transform and ordination trail before replaying the formal audit.

The microscopy object releases can be reanalyzed locally. Their raw imagery,
segmentation models, and production runs remain in the sibling
`cellprofiler_test` project or external storage. The copied production-lineage
JSON/Markdown identifies that upstream boundary. Absolute upstream filesystem
paths in the frozen object tables are provenance strings and may not remain
resolvable on another machine.

The author confirmed that Process_413/specimen 32469 and
Process_414/specimen 32470 were the *D. fuscus* DNA-content standards and were
stained with all retained unknowns in one Feulgen run. The retained standard
images were 40× brightfield without oil, using the green channel and the same
fixation protocol as the other slides. The author-confirmed run map is
retained in `outputs/calibration_search_20261013/`. Standard-nucleus acceptance
and the approximately 41% difference between the two standard-slide median IOD
values still require technical QC. The selected 16.462167689-pg assembly-derived
scale remains a conditional conversion, not an independently validated absolute
C-value; its denominator remains the reviewed target-species IOD estimate.

## Release checks

`make test` verifies:

- ordered S01–S43 files equal their canonical sources byte-for-byte;
- manifest hashes and dimensions are current;
- natural-log Shannon entropy and the Gini-Simpson index reproduce for all 34
  TE34 species after deterministic selection and explicit reclosure of each
  species' ten most abundant positive superfamilies; the separate PCA uses the
  24 superfamilies positive in every species, while unresolved aligned-base mass
  is excluded from both analyses and retained only as S06-S07 dnaPipeTE
  quality-control evidence;
- LTR release hashes, formulas, whitespace-delimited TEsorter domain count,
  408-element eligible cohort, original within-species two-sided IQR filter,
  regional positive-depth gate, 380-element primary intersection, and
  historical source-corrupt rows reproduce;
- morphology, IOD, and conditional genome-size species/image summaries, animal
  support, physical scale, and quality-balance decisions reproduce;
- the selected 16.1-Gbp pooled-assembly anchor converts to 16.462167689 pg/1C
  with the declared 0.978-Gbp-per-pg relation; the Path24 analysis uses that
  conditional genome-size scale, while S21-S25 preserve the unscaled IOD and
  calibration audit;
- TE34, Path24, and overlap21 membership;
- the exact SVL aliases and registered max-of-maxima aggregation;
- the Path24 tree is the exact representative-tip prune of the source tree;
- the accession-specific fuscus identity decision across every active panel;
- the canonical data catalog's paths, byte sizes, and SHA-256 hashes.

Run `make identity` for the corresponding human-readable, read-only report.
The complete contract and current exceptions are documented in
[IDENTITY_RESOLUTION.md](IDENTITY_RESOLUTION.md).

The unresolved focal-tree citation, DNA-standard QC and absolute calibration,
and remaining wet-lab author inputs remain visible publication blockers even
when the computational checks pass.
