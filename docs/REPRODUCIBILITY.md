# Reproducibility boundary

The project now has three explicit reproducibility lanes. This distinction is
scientific provenance, not merely storage policy.

## 1. Locally rerunnable

- data catalog, hash validation, and S01–S40 publication exports;
- the three compact exploration notebooks;
- the complete TE34 order/superfamily CLR-PCA reconstruction and formal final
  PCA audit from the preserved compact composition tables, including an exact
  replay of all 136 Shannon/Gini-Simpson index rows and their observed-
  richness audit counts;
- compact morphology/IOD release validation from the frozen object tables,
  including species estimates, image aggregation, pixel scale, animal support,
  and the declared quality-balance gate;
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

## Release checks

`make test` verifies:

- ordered S01–S40 files equal their canonical sources byte-for-byte;
- manifest hashes and dimensions are current;
- Shannon entropy, the Gini-Simpson index, observed-richness audit counts, and
  all four TE34 strata reproduce; classified-conditional order-level results
  are primary, while superfamily and mass-aware results are sensitivity/audit;
- LTR release hashes, formulas, whitespace-delimited TEsorter domain count,
  408-element eligible cohort, original within-species two-sided IQR filter,
  regional positive-depth gate, 380-element primary intersection, and
  historical source-corrupt rows reproduce;
- morphology and IOD species/image summaries, animal support, physical scale,
  and quality-balance decisions reproduce;
- TE34, Path24, and overlap21 membership;
- the exact SVL aliases and registered max-of-maxima aggregation;
- the Path24 tree is the exact representative-tip prune of the source tree;
- the accession-specific fuscus identity decision across every active panel;
- the canonical data catalog's paths, byte sizes, and SHA-256 hashes.

Run `make identity` for the corresponding human-readable, read-only report.
The complete contract and current exceptions are documented in
[IDENTITY_RESOLUTION.md](IDENTITY_RESOLUTION.md).

The unresolved focal-tree citation and any missing wet-lab author inputs remain
visible publication blockers even when the computational checks pass.
