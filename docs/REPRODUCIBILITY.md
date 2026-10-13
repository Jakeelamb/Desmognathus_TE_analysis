# Reproducibility boundary

The project now has three explicit reproducibility lanes. This distinction is
scientific provenance, not merely storage policy.

## 1. Locally rerunnable

- data catalog, hash validation, and S01–S41 publication exports;
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

Python uses one UV lock. R and figure conversion use one Conda environment.
Run `make setup`, then `make test`. The R setup checksum-verifies the exact
`phylopath` 1.3.1 CRAN source because that package is not available from the
configured Conda channels.

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

The full genomic pipeline cannot currently be claimed reproducible because the
historical custom TE library `dedupe_telib.fasta` and complete HPC software/run
provenance are absent. The archive preserves observed outputs rather than
fabricating a rerun path.

The canonical diversity table itself can be promoted from the frozen
superfamily composition matrix with
`uv run --frozen --no-sync python analyses/01_transposable_elements/recompute_pca.py --write`;
`make te-pca` must then pass as the read-only end-to-end gate.

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

- ordered S01–S41 files equal their canonical sources byte-for-byte;
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
