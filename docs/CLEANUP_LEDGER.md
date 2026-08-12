# Repository cleanup ledger

Cleanup date: 2026-08-10
Pre-cleanup branch: `codex/largest-100-cells-exploration`
Pre-cleanup HEAD: `7c235b808efc08117ffda246890e425e464fdb74`

## Preservation gate

Before any large deletion, the complete dirty checkout was copied to
`/home/jake/Projects/Desmognathus_TE_archive`. Its manifest records 35,731
non-cache working files totaling 64,377,618,470 bytes. Every manifest entry was
independently rehashed and matched. The manifest SHA-256 is:

```text
c07f866f561f6e48e9b8f07e31b592de46790888c3a625866abf78bfe2e90b6f
```

The archive also retains its full `.git` directory and the exact dirty working
state. Transient caches and `.git` internals are intentionally outside the
per-file working-data manifest.

## Replacement map

| Retired active hierarchy | Current authority |
| --- | --- |
| `input_data/lookup_table.txt` and scattered panels | `data/identity/` |
| `results/data/corrected/te34*` and TE review products | `analyses/01_transposable_elements/data/` |
| microscopy releases spread across results/path workspaces | `analyses/02_morphology/data/` and `analyses/02_morphology/provenance/` |
| mixed `path_analysis/` inputs, outputs, and viewers | `analyses/03_phylogenetic_path/` |
| eight overlapping review notebooks | one data/methods report plus one manifest-driven figure-review gallery under `Publication/` |
| several publication builders | `make publication` |
| `Dusky.yml` plus a figure-only Conda spec | `environment.yml` plus `scripts/setup_r.sh` |

## Removed from the active checkout

The following remain in the checksum archive but no longer exist in the paper
workspace:

- approximately 42 GB of genomic/HPC input data, including dnaPipeTE,
  RepeatMasker `.align`, LTR/depth evidence, genome caches, and old source trees;
- approximately 7.8 GB of dnaPipeTE/divergence intermediates;
- approximately 11 GB of mixed results, hit-level substrates, audit figures,
  and retired analysis18 products;
- approximately 376 MB of mixed path-analysis downloads, generated galleries,
  alternate panels, and model branches;
- legacy notebooks, plans, root-level figures and decision-table duplicates;
- old processing/visualization scripts whose inputs and contracts no longer
  define the paper release;
- local tool caches and obsolete editor/project settings.

The active working-file surface is approximately 25.5 MiB. `.git` history and
the local `.venv` are regenerable/storage layers and are not part of that
analysis-data figure.

## Deletion audit

- The old publication dataset builder is replaced by `scripts/project.py`; it
  is removed from active code.
- The separate figure environment and historical mixed Python/R/HPC environment
  are replaced by one R environment plus one UV Python lock; both old installed
  project environments were removed after the replacement ran figures and the
  Path24 quick suite successfully.
- Multi-layer frozen registries are replaced for day-to-day navigation by
  `data/DATA_CATALOG.csv`; the detailed historical registries remain in the
  archive, while the compact microscopy decision lineage remains active.
- Retired Cell21/final18/integrated-path branches are absent from the active
  workflow. Their results and code remain archive-recoverable.
- Generated path output and test/tool caches are ignored and were removed after
  validation.
- The five-row `organismal_trait_coverage.csv` snapshot was removed when the
  expert reproductive-strategy field was added. Its saved 18-species overlap
  counts predated the TE34/Path24/overlap21 contract, no active code consumed
  it, and current coverage is directly computable from the canonical
  one-row-per-species `organismal_traits.csv`; the historical file remains in
  the checksum archive.
- The wide `path_input_master.csv` and `analysis_species_readiness.csv`
  snapshots were removed after the *D. orestes* SVL alias correction made
  their duplicated *D. orestes* fields stale. Neither had an active builder or
  consumer; their source values remain in the canonical component tables and
  both historical snapshots remain recoverable from the checksum archive.
- The one-off `outputs/te_pca_structure_20260810/TE_PCA_STRUCTURE_AUDIT.xlsx`
  workbook and its empty parent directories were removed. It was untracked,
  uncataloged, unreferenced, had no generator or unique result, duplicated the
  canonical PCA audit CSVs, and named a nonexistent CLR source file. The
  canonical CSVs and the shareable R Markdown report preserve every result it exposed.
- The new `Publication/Desmognathus_figure_data.xlsx` does not restore that
  one-off PCA workbook or create a parallel data authority. One reproducible
  builder collects the exact non-tree direct ggplot inputs into labeled sheets.
  The unusable Newick, tree-edge, and tree-node spreadsheet tabs were removed;
  their manifest-bound CSV and Newick sources remain authoritative under
  `Publication/datasets/`.
- The PCA verifier's disposable `--write` path was removed. `make te-pca` now
  performs one read-only reconstruction and formal replay against the canonical
  files, while `make report` renders the data/methods workbench and
  `make figure-review` renders the visual critique gallery.
- The old `te-structure` Make target was replaced by `te-pca`; no compatibility
  alias or parallel workflow was retained. The verifier now reconstructs the
  complete order and superfamily prevalence/CLR trail rather than checking only
  the final superfamily scores.
- The three generated `explore.ipynb` files and their Python notebook builder
  were replaced by one shareable R Markdown/HTML report. The `notebooks` and
  `te-pca-view` Make targets, JupyterLab, `nbformat`, and their lockfile-only
  dependency tree were deleted; canonical CSVs and formal replay commands
  remain the analysis authorities.
- The redundant `zero_replaced_proportion` field was deleted from
  `te_pca_clr_matrix.csv`. It was byte-for-byte duplicate stage metadata even
  though the retained shared-positive matrix required no zero replacement;
  `closed_proportion` remains the single compositional input authority, and all
  PCA scores, loadings, variance, and formal statistical results are unchanged.
- Source-label logic is now localized behind one read-only `make identity`
  interface. Exact aliases moved out of taxonomy/prose into
  `source_species_aliases.csv`, maximum-SVL aggregation moved into the structured
  trait registry, and the old inline tree/fuscus checks moved out of
  `scripts/project.py` into `scripts/identity.py`.
- The inaccurate crosswalk claim that `orestes_B` was collapsed was deleted.
  The active tree keeps exact `orestes` and `fuscus` representatives and prunes
  their alternates; no generic suffix mapper, duplicate R implementation,
  generated identity report, compatibility command, or output directory was
  added.
- The misleading `relative_nuclear_iod_by_image.csv` name and matching S23
  publication filename were removed. The one canonical image-level table is now
  `nuclear_iod_by_image.csv`, its animal field is `specimen_id`, and its
  metadata explicitly states that values are unnormalized rather than ratios.
- Diversity replay was added to the existing PCA verifier instead of creating a
  second diversity builder or output directory. The two compact composition
  inputs now reproduce one 68-row, five-column classified-only table: 34 order
  and 34 superfamily rows after classified categories are reclosed to one. The
  order stratum is primary, the superfamily stratum is sensitivity/audit, and
  observed richness remains only an audit/support count. The unresolved-bin
  diversity mode, its redundant `composition_mode` column, and the resulting 68
  rows were deleted; unresolved aligned-base mass remains only as S06-S07
  dnaPipeTE quality-control evidence. Redundant Hill-number and Pielou fields
  and claims were removed rather than carried forward.
- S19 was extended in place with animal count and concentration fields derived
  from S18. No parallel morphology summary was retained.
- One compact LTR release record now binds S14-S17 and the sibling archive
  evidence by hash. S14 keeps all 1,086 historically selected usable elements
  for audit and corrects the archived pipe-delimiter domain-count error in
  place. The 408 eligible LTR/Gypsy five-domain elements, original
  within-species two-sided 1.5-IQR decisions, and independent regional
  positive-depth gate are explicit row fields. Their 380-element intersection
  is the primary S15 branch; no second element table or destructive filtered
  copy was added. S16 now separates historical and corrected counts, S17 keeps
  the two four-domain corrupt historical rows as audit evidence, and the later
  compound one-element reconstruction remains explicitly retired.
- The blank `GBE_methods_reference.docx` shell and the obsolete
  `Materials_and_Methods_GBE_draft.docx`/`.md` pair were removed from the
  active Methods directory. The old prose encoded superseded panel sizes and
  object counts and competed with the factual scaffold; exact copies remain in
  the checksum archive and every tracked version remains recoverable from Git
  history. `METHODS_BULLET_SCAFFOLD.md` remains the sole active factual Methods
  authority. The current Word-ready scaffold and narrative DOCX are synchronized
  review exports and do not create a second factual authority.
- The unsupported 16.36-pg *D. fuscus* convention was replaced in place by one
  documented conditional reference: the Myers et al. pooled assembly estimate
  (`16.1 Gbp / 0.978 Gbp per pg = 16.462167689 pg/1C`). No parallel absolute-
  genome-size dataset was added. Process_413/specimen 32469 and
  Process_414/specimen 32470 are confirmed same-run standards, but the exact
  standard-to-target run map and standard-slide QC remain unresolved; relative
  IOD therefore remains the primary phenotype.

## Known blockers preserved rather than hidden

- the historical custom `dedupe_telib.fasta` is missing, so full genomic
  upstream regeneration cannot be claimed;
- focal 46-tip tree citation/calibration provenance remains unresolved;
- wet-lab, specimen, and image-acquisition author inputs remain incomplete;
- relative nuclear IOD is not an independently validated absolute C-value; the
  selected 16.462167689-pg assembly-derived scale remains conditional pending
  exact standard-to-target run mapping and standard QC;
- archived HPC and microscopy outputs are evidence, not a promise of complete
  rerunnability.
