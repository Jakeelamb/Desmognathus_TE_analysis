# Reproducibility boundary

The project now has three explicit reproducibility lanes. This distinction is
scientific provenance, not merely storage policy.

## 1. Locally rerunnable

- data catalog, hash validation, and S01–S40 publication exports;
- the R Markdown study-data report and manifest-driven figure-review gallery,
  each with a self-contained HTML rendering;
- the single labeled figure-data workbook, rebuilt from the manifest-bound CSV
  inputs other than the tree files with exact sheet/source parity checks; tree
  inputs remain in their standard Newick and reviewer-readable CSV formats;
- the complete TE34 order/superfamily CLR-PCA reconstruction and formal final
  PCA audit from the preserved compact composition tables, including an exact
  replay of the 68-row, five-column classified-only diversity table (34 order
  and 34 superfamily rows) and its observed-richness audit counts;
- compact morphology/IOD release validation from the frozen object tables,
  including species estimates, image aggregation, pixel scale, animal support,
  and the declared quality-balance gate;
- the arithmetic conditional *D. fuscus* reference conversion from the published
  pooled assembly (`16.1 Gbp / 0.978 Gbp per pg = 16.462167689 pg/1C`), while
  relative IOD remains the primary phenotype;
- the corrected LTR/Gypsy five-domain gate, primary IQR/coverage branch, and
  species/resource summaries from the frozen 1,086-row zero-aware audit table;
- the Path24 R analysis from the preserved compact traits and trees;
- the current GBE ggplot figure bundle.

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
- finalized reviewed cell/nucleus objects and nuclear-IOD objects;
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

This upstream boundary does not prevent local PCA verification: `make te-pca`
starts at the frozen, species-level composition matrices and reconstructs the
complete transform and ordination trail before replaying the formal audit.

The microscopy object releases can be reanalyzed locally. Their raw imagery,
segmentation models, and production runs remain in the sibling
`cellprofiler_test` project or external storage. The copied production-lineage
JSON/Markdown identifies that upstream boundary. Absolute upstream filesystem
paths in the frozen object tables are provenance strings and may not remain
resolvable on another machine.

Process_413/specimen 32469 and Process_414/specimen 32470 are author-confirmed
*D. fuscus* DNA-content standards that were stained in the same experimental
runs as the unknowns. The exact standard-to-target run map has not yet been
recovered, and the approximately 41% difference between the two standard-slide
median IOD values still requires technical QC. Consequently, the selected
16.462167689-pg assembly-derived scale is reproducible as a conditional
conversion but is not an independently validated absolute C-value.

## Release checks

`make test` verifies:

- ordered S01–S40 files equal their canonical sources byte-for-byte;
- manifest hashes and dimensions are current;
- Shannon entropy, the Gini-Simpson index, observed-richness audit counts, and
  both 34-species TE34 strata reproduce after classified categories are reclosed
  to one; order-level results are primary and superfamily-level results are
  sensitivity/audit, while unresolved aligned-base mass is excluded from
  diversity and retained only as S06-S07 dnaPipeTE quality-control evidence;
- LTR release hashes, formulas, whitespace-delimited TEsorter domain count,
  408-element eligible cohort, original within-species two-sided IQR filter,
  regional positive-depth gate, 380-element primary intersection, and
  historical source-corrupt rows reproduce;
- morphology and IOD species/image summaries, animal support, physical scale,
  and quality-balance decisions reproduce;
- the selected 16.1-Gbp pooled-assembly anchor converts to 16.462167689 pg/1C
  with the declared 0.978-Gbp-per-pg relation, without changing the primary
  relative-IOD values;
- TE34, Path24, and overlap21 membership;
- the exact SVL aliases and registered max-of-maxima aggregation;
- the Path24 tree is the exact representative-tip prune of the source tree;
- the accession-specific fuscus identity decision across every active panel;
- the canonical data catalog's paths, byte sizes, and SHA-256 hashes.

Run `make identity` for the corresponding human-readable, read-only report.
The complete contract and current exceptions are documented in
[IDENTITY_RESOLUTION.md](IDENTITY_RESOLUTION.md).

The unresolved focal-tree citation, exact DNA-standard run map/QC, and remaining
wet-lab author inputs remain visible publication blockers even when the
computational checks pass.
