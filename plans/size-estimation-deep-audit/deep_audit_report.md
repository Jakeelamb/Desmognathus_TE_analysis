# Size Estimation Deep Audit

Date: 2026-07-08

Scope: read-only `/improve` audit of the current cell-size, nucleus-size, genome-size, top-50 review, import, notebook, model, and provenance surfaces. No source files were changed during this audit. Audit-only plots and tables were generated under this `plans/` directory.

## Bottom Line

The current top-50 dataset is real and internally coherent enough for exploratory work: it has 1,050 selected rows, 21 species, and 50 linked cell/nucleus pairs per species. All 21 current genome/morphology species have matching tree tips.

The analysis is not yet publication-grade. The main problem is not one bad plot. The main problem is that the evidence chain can still produce confident-looking numbers from weak or ambiguous conditions:

- the Desmognathus import can default to a broader non-top50 verified dataset or silently fall back to legacy reconstruction;
- the genome QC shift metric mixes image filtering with reference recalibration;
- one-image/one-specimen species can still receive object-bootstrap CIs that look precise;
- final top-50 artifacts can contain uncertain review states and stale decision columns;
- primary comparative panels are mostly sensitivity-limited genome estimates;
- model warnings and stale stability summaries are not surfaced strongly enough.

The fix is a short sequence of fail-closed provenance checks, estimator semantics fixes, uncertainty relabeling, and stronger notebook/model diagnostics.

## Audit Artifacts

- [current_size_qc_audit_table.csv](/home/jake/Projects/Desmognathus_TE/plans/size-estimation-deep-audit/current_size_qc_audit_table.csv)
- [phylogeny_size_qc_audit.png](/home/jake/Projects/Desmognathus_TE/plans/size-estimation-deep-audit/figures/phylogeny_size_qc_audit.png)
- [size_qc_variable_relationships.png](/home/jake/Projects/Desmognathus_TE/plans/size-estimation-deep-audit/figures/size_qc_variable_relationships.png)
- [genome_qc_sensitivity_ranked.png](/home/jake/Projects/Desmognathus_TE/plans/size-estimation-deep-audit/figures/genome_qc_sensitivity_ranked.png)
- [make_audit_plots.py](/home/jake/Projects/Desmognathus_TE/plans/size-estimation-deep-audit/make_audit_plots.py)

## Current Numeric Snapshot

- Top-50 source: `selection_mode=exact_selected_pairs`, 1,050 selected rows, target 50 linked pairs per species in [summary.json](/home/jake/Projects/cellprofiler_test/output/runs/mixed_cellpose_yolo_full_dataset_v1_bgclean/verified_species_dataset_top50_latest/summary.json:4), [summary.json](/home/jake/Projects/cellprofiler_test/output/runs/mixed_cellpose_yolo_full_dataset_v1_bgclean/verified_species_dataset_top50_latest/summary.json:16), and [summary.json](/home/jake/Projects/cellprofiler_test/output/runs/mixed_cellpose_yolo_full_dataset_v1_bgclean/verified_species_dataset_top50_latest/summary.json:24).
- Non-top50 verified latest source: 2,065 selected rows and target 100 linked pairs per species in [summary.json](/home/jake/Projects/cellprofiler_test/output/runs/mixed_cellpose_yolo_full_dataset_v1_bgclean/verified_species_dataset_latest/summary.json:14) and [summary.json](/home/jake/Projects/cellprofiler_test/output/runs/mixed_cellpose_yolo_full_dataset_v1_bgclean/verified_species_dataset_latest/summary.json:22).
- Current top-50 row source mix: 440 manual selections and 610 auto selections.
- Genome result status: 15 sensitivity-limited, 6 minor-caution, 0 stable.
- Genome source: 19 image-QC-pass estimates, 2 all-selected fallbacks.
- Genome primary support: 13 limited, 4 medium, 2 low, 2 missing primary.
- Image-QC-pass genome shift vs all-selected: 19 species, minimum -22.35%, median -12.40%, maximum 0%.
- Largest image-QC shifts: D. auriculatus -22.35%, D. perlapsus -21.73%, D. welteri -19.59%, D. kanawha -19.01%, D. amphileucus -18.29%.
- Raw Spearman relationships from the current notebook sidecar: genome-cell 0.317, genome-nucleus 0.491, genome-N:C 0.153, cell-nucleus 0.618, cell-N:C -0.343, nucleus-N:C 0.426.

Interpretation boundary: those correlations are raw species-level relationships. They should be shown beside phylogenetic-corrected diagnostics before being used as biological evidence.

## P0 Findings

### P0.1 The import is not fail-closed to the curated top-50 dataset

Impact: A normal documented refresh can silently import a different data universe than the reviewed top-50 dataset.

Evidence:

- The bridge default is `verified_species_dataset_latest`, not `verified_species_dataset_top50_latest`: [bridge.py](/home/jake/Projects/Desmognathus_TE/path_analysis/scripts/cellprofiler_bridge/bridge.py:26).
- The path-analysis README still documents the non-top50 `verified_species_dataset_latest`: [README.md](/home/jake/Projects/Desmognathus_TE/path_analysis/README.md:85).
- The provenance audit says the current snapshot is the top50 source: [CELLPROFILER_PROVENANCE_AUDIT.md](/home/jake/Projects/Desmognathus_TE/path_analysis/CELLPROFILER_PROVENANCE_AUDIT.md:22).
- Missing verified files return `present: False` instead of failing: [bridge.py](/home/jake/Projects/Desmognathus_TE/path_analysis/scripts/cellprofiler_bridge/bridge.py:579).
- The rebuild path then falls back to legacy genome and morphology builders: [bridge.py](/home/jake/Projects/Desmognathus_TE/path_analysis/scripts/cellprofiler_bridge/bridge.py:1780) and [bridge.py](/home/jake/Projects/Desmognathus_TE/path_analysis/scripts/cellprofiler_bridge/bridge.py:1787).
- The top50 artifact is 1,050 rows; the default latest verified artifact is 2,065 rows.

Required fix:

- Make top50 exact-selected input the default, or require `--verified-species-dir` explicitly.
- Fail if `selection_mode != exact_selected_pairs`, `target_linked_per_species != 50`, row count is not 1,050, or the selected-pairs hash does not match the imported provenance.
- Move legacy reconstruction behind an explicit `--allow-legacy-cellprofiler-fallback` flag.

### P0.2 Genome QC sensitivity currently mixes filtering effects with reference recalibration

Impact: The reported percent shift cannot be interpreted as the effect of dropping low-quality images alone.

Evidence:

- `all_selected`, `image_qc_pass`, and `high_iod_qc` are each re-analyzed through `safe_analyze_subset`: [build_verified_species_dataset.py](/home/jake/Projects/cellprofiler_test/scripts/build_verified_species_dataset.py:405), [build_verified_species_dataset.py](/home/jake/Projects/cellprofiler_test/scripts/build_verified_species_dataset.py:413), and [build_verified_species_dataset.py](/home/jake/Projects/cellprofiler_test/scripts/build_verified_species_dataset.py:421).
- `pct_shift_vs_all_selected` compares already recalibrated outputs: [build_verified_species_dataset.py](/home/jake/Projects/cellprofiler_test/scripts/build_verified_species_dataset.py:242).
- The D. fuscus scale changes across states: all-selected 0.017128 pg/IOD, image-QC-pass 0.015004 pg/IOD, high-IOD-QC 0.015568 pg/IOD in [summary.json](/home/jake/Projects/cellprofiler_test/output/runs/mixed_cellpose_yolo_full_dataset_v1_bgclean/verified_species_dataset_top50_latest/summary.json:39), [summary.json](/home/jake/Projects/cellprofiler_test/output/runs/mixed_cellpose_yolo_full_dataset_v1_bgclean/verified_species_dataset_top50_latest/summary.json:48), and [summary.json](/home/jake/Projects/cellprofiler_test/output/runs/mixed_cellpose_yolo_full_dataset_v1_bgclean/verified_species_dataset_top50_latest/summary.json:57).
- D. anicetus has 50 pairs, 1 image, 1 specimen in both all-selected and image-QC-pass, yet reports -12.400463% shift: [species_genome_sensitivity.csv](/home/jake/Projects/cellprofiler_test/output/runs/mixed_cellpose_yolo_full_dataset_v1_bgclean/verified_species_dataset_top50_latest/species_genome_sensitivity.csv:5) and [species_genome_sensitivity.csv](/home/jake/Projects/cellprofiler_test/output/runs/mixed_cellpose_yolo_full_dataset_v1_bgclean/verified_species_dataset_top50_latest/species_genome_sensitivity.csv:7).
- The bridge has an `alternate_genome_pg` channel but leaves `cross_state_pct_diff` blank: [bridge.py](/home/jake/Projects/Desmognathus_TE/path_analysis/scripts/cellprofiler_bridge/bridge.py:522).

Required fix:

- Emit two separate sensitivity columns:
  - fixed-reference filtering effect;
  - fully re-estimated subset effect.
- Populate `cross_state_pct_diff` with the chosen primary definition.
- Make the notebook plot both definitions and label them plainly.

### P0.3 One-image and one-specimen estimates can look too precise

Impact: Some species have weak independent support but still show tight-looking CIs and large object-level effective n.

Evidence:

- Effective n is row-weight based: [build_balanced_species_estimates.py](/home/jake/Projects/cellprofiler_test/scripts/build_balanced_species_estimates.py:402).
- Bootstrap resampling falls back when there is only one independent specimen/image level: [build_balanced_species_estimates.py](/home/jake/Projects/cellprofiler_test/scripts/build_balanced_species_estimates.py:430).
- Percentile CIs are written for all species: [build_balanced_species_estimates.py](/home/jake/Projects/cellprofiler_test/scripts/build_balanced_species_estimates.py:592).
- D. anicetus has 1 image and 1 specimen, but `genome_primary_effective_n=49.800929` and a genome CI of 12.580091-13.763685 pg: [species_estimates_verified.csv](/home/jake/Projects/cellprofiler_test/output/runs/mixed_cellpose_yolo_full_dataset_v1_bgclean/verified_species_dataset_top50_latest/species_estimates_verified.csv:3).

Required fix:

- Report independent support units separately: `n_pairs`, `n_images`, `n_specimens`, `object_effective_n`, and `independent_support_n`.
- For one-image or one-specimen species, label intervals as object-resampling-only or suppress them from publication-facing summaries.
- In plots, encode support level visually instead of making all error bars look equivalent.

### P0.4 Final top-50 review states are not cleanly final

Impact: The final verified top-50 can contain uncertain states and stale decision columns that disagree with the effective decision.

Evidence:

- The freeze validation catches duplicate review keys but does not re-enforce the full linked/non-edge/non-overlap invariant set: [freeze_species_grid_top50.py](/home/jake/Projects/cellprofiler_test/scripts/freeze_species_grid_top50.py:113).
- `maybe` is mapped into `auto_grid_reviewed_maybe`, not excluded: [build_verified_species_dataset.py](/home/jake/Projects/cellprofiler_test/scripts/build_verified_species_dataset.py:109).
- Current final curated top50 has 5 `maybe` rows.
- Current final curated top50 has 207 rows where `decision` and `grid_effective_decision` disagree.
- Example uncertain row: final CSV line 301 has `grid_effective_decision=maybe`: [final_curated_top50_linked_pairs.csv](/home/jake/Projects/cellprofiler_test/output/runs/mixed_cellpose_yolo_full_dataset_v1_bgclean/top50_linked_pair_review/final_curated_top50_latest/final_curated_top50_linked_pairs.csv:301).
- Example stale-decision row: final CSV line 2 has `decision=unlabeled` and `grid_effective_decision=keep`: [final_curated_top50_linked_pairs.csv](/home/jake/Projects/cellprofiler_test/output/runs/mixed_cellpose_yolo_full_dataset_v1_bgclean/top50_linked_pair_review/final_curated_top50_latest/final_curated_top50_linked_pairs.csv:2).

Required fix:

- Treat `maybe` as non-final. It should require explicit keep/discard before freeze.
- Canonicalize one final decision column in frozen outputs; move legacy decision columns to provenance-only names or drop them.
- Re-run invariant validation at freeze and import: linked cell+nucleus, one-to-one nucleus overlap, non-edge cell, non-cut tile, positive area, mask path present, and no multi-cell overlap.

### P0.5 All-selected fallbacks are underclassified

Impact: Two species with missing image-QC-pass primary genomes are still labeled only `minor_caution`.

Evidence:

- Final species estimates fall back to all-selected when the primary is missing: [build_verified_species_dataset.py](/home/jake/Projects/cellprofiler_test/scripts/build_verified_species_dataset.py:336).
- D. folkertsi is `all_selected_fallback`, `genome_primary_missing=True`, but `result_status=minor_caution`: [cellprofiler_final_species_results.csv](/home/jake/Projects/Desmognathus_TE/path_analysis/data/external/derived/cellprofiler_final_species_results.csv:9).
- D. ochrophaeus is also `all_selected_fallback`, `genome_primary_missing=True`, `result_status=minor_caution`, and has a literal `nan; all_selected_genome_fallback` flag: [cellprofiler_final_species_results.csv](/home/jake/Projects/Desmognathus_TE/path_analysis/data/external/derived/cellprofiler_final_species_results.csv:17).
- The bridge appends `all_selected_genome_fallback` when `genome_primary_missing` is set, then still allows medium support to map to `minor_caution`: [bridge.py](/home/jake/Projects/Desmognathus_TE/path_analysis/scripts/cellprofiler_bridge/bridge.py:216) and [bridge.py](/home/jake/Projects/Desmognathus_TE/path_analysis/scripts/cellprofiler_bridge/bridge.py:221).

Required fix:

- Any `genome_primary_missing` or `all_selected_fallback` row should be at least `sensitivity_limited`.
- Build warning strings with explicit non-null checks so missing values do not become literal `nan`.

### P0.6 Primary genome panels are mostly sensitivity-limited

Impact: The current "primary" comparative panels can sound stronger than they are.

Evidence:

- Panel masks require `has_genome` but do not gate on `genome_result_status`: [build_analysis_panels.py](/home/jake/Projects/Desmognathus_TE/path_analysis/scripts/build_analysis_panels.py:307).
- Genome status counts are summarized after panel creation: [build_analysis_panels.py](/home/jake/Projects/Desmognathus_TE/path_analysis/scripts/build_analysis_panels.py:509).
- `te_genome_primary` currently has 18 species, 0 stable genomes, 4 minor-caution genomes, and 14 sensitivity-limited genomes: [analysis_panel_summary.csv](/home/jake/Projects/Desmognathus_TE/path_analysis/data/derived/analysis_panel_summary.csv:3).

Required fix:

- Rename broad panels as sensitivity-inclusive, or add stable/minor-only panels and make those the conservative primary analysis.
- In every figure/model table, show `n_stable`, `n_minor`, and `n_sensitivity_limited`.

## P1 Findings

### P1.1 Top-50 source sidecars are not fully in the general source registry

Impact: The source traceability audit can pass while omitting key top50 evidence.

Evidence:

- The general registry script registers only `cellprofiler_final_species_results` and `cellprofiler_species_morphology_summary`: [audit_source_traceability.py](/home/jake/Projects/Desmognathus_TE/path_analysis/scripts/audit_source_traceability.py:78) and [audit_source_traceability.py](/home/jake/Projects/Desmognathus_TE/path_analysis/scripts/audit_source_traceability.py:82).
- The current registry mirrors only those two CellProfiler records: [source_file_registry.csv](/home/jake/Projects/Desmognathus_TE/path_analysis/data/derived/source_file_registry.csv:55) and [source_file_registry.csv](/home/jake/Projects/Desmognathus_TE/path_analysis/data/derived/source_file_registry.csv:56).
- Discovery JSON does include genome sensitivity and image IOD quality sidecars: [cellprofiler_source_discovery.json](/home/jake/Projects/Desmognathus_TE/path_analysis/data/external/derived/cellprofiler_source_discovery.json:66) and [cellprofiler_source_discovery.json](/home/jake/Projects/Desmognathus_TE/path_analysis/data/external/derived/cellprofiler_source_discovery.json:72).

Required fix:

- Promote top50 upstream files, selected-pair source, genome sensitivity, image IOD quality, and reconstruction audit hashes into `source_file_registry.csv`.
- Make the traceability audit fail if any required top50 sidecar is absent.

### P1.2 Pair-level reviewed provenance is aggregated away

Impact: Desmognathus artifacts cannot independently reconstruct which exact reviewed cell/nucleus pairs made each species estimate.

Evidence:

- The bridge groups selected pairs by species and filename: [bridge.py](/home/jake/Projects/Desmognathus_TE/path_analysis/scripts/cellprofiler_bridge/bridge.py:287).
- It stores selected-pairs source path/hash after aggregation: [bridge.py](/home/jake/Projects/Desmognathus_TE/path_analysis/scripts/cellprofiler_bridge/bridge.py:294).

Required fix:

- Add an imported pair-level trace sidecar with `review_key`, species, source image/tile, cell object id, nucleus label, area/IOD metrics, effective decision, QC flags, and upstream source hash.

### P1.3 Model warnings are hidden from the main evidence surfaces

Impact: Model rankings look cleaner than the phylogenetic fitting diagnostics.

Evidence:

- `phylo_path` warnings are written only to sidecar warning files: [path_model_scaffold.R](/home/jake/Projects/Desmognathus_TE/path_analysis/scripts/path_model_scaffold.R:510) and [path_model_scaffold.R](/home/jake/Projects/Desmognathus_TE/path_analysis/scripts/path_model_scaffold.R:513).
- Current morphology model warnings include lambda-bound warnings: [genome_morphology_model_warnings.txt](/home/jake/Projects/Desmognathus_TE/path_analysis/results/genome_morphology_model_warnings.txt:1).
- Current TE-morphology model warnings include the same class: [te_genome_morphology_model_warnings.txt](/home/jake/Projects/Desmognathus_TE/path_analysis/results/te_genome_morphology_model_warnings.txt:1).
- The notebook reads rankings and PDFs but not model-warning sidecars.

Required fix:

- Parse warnings into every ranking summary: `warning_count`, `warning_class`, `has_lambda_boundary_warning`.
- Add a notebook warning table and mark affected model runs in figures.

### P1.4 Stale model-stability artifacts can contradict current rankings

Impact: A summary can state a different winning model than the current result files.

Evidence:

- Current `te_genome_organismal_primary_mediumplus_model_ranking.csv` winner is `te_baseline`, while the older stability audit reports `body_size_additive`: [te_genome_organismal_primary_mediumplus_model_ranking.csv](/home/jake/Projects/Desmognathus_TE/path_analysis/results/te_genome_organismal_primary_mediumplus_model_ranking.csv:2) and [MODEL_STABILITY_AUDIT.md](/home/jake/Projects/Desmognathus_TE/path_analysis/MODEL_STABILITY_AUDIT.md:36).
- Current `te_genome_ectopic_organismal_primary_mediumplus_model_ranking.csv` winner is `ectopic_baseline`, while the older stability audit reports `te_body_size_baseline`: [te_genome_ectopic_organismal_primary_mediumplus_model_ranking.csv](/home/jake/Projects/Desmognathus_TE/path_analysis/results/te_genome_ectopic_organismal_primary_mediumplus_model_ranking.csv:2) and [MODEL_STABILITY_AUDIT.md](/home/jake/Projects/Desmognathus_TE/path_analysis/MODEL_STABILITY_AUDIT.md:50).
- Current `te_genome_ltr_history_primary_mediumplus_model_ranking.csv` winner is `te_baseline`, while the older stability audit reports `history_additive`: [te_genome_ltr_history_primary_mediumplus_model_ranking.csv](/home/jake/Projects/Desmognathus_TE/path_analysis/results/te_genome_ltr_history_primary_mediumplus_model_ranking.csv:2) and [MODEL_STABILITY_AUDIT.md](/home/jake/Projects/Desmognathus_TE/path_analysis/MODEL_STABILITY_AUDIT.md:64).

Required fix:

- Regenerate stability artifacts from current panels.
- Stamp stability summaries with input hashes and generation timestamps.
- Make stale stability files fail the verification gate.

### P1.5 The notebook needs matched raw and phylogenetic relationship views

Impact: Readers can conflate raw cross-species correlations with phylogenetic evidence.

Evidence:

- The notebook produces raw scatterplots and raw Spearman summaries.
- The audit plot [phylogeny_size_qc_audit.png](/home/jake/Projects/Desmognathus_TE/plans/size-estimation-deep-audit/figures/phylogeny_size_qc_audit.png) shows current traits on the tree, but this is not yet integrated into the official notebook.

Required fix:

- Add a tree-pruned trait heatmap/bar plot for genome, cell, nucleus, N:C, QC status, and uncertainty.
- Add phylogenetic residual or PGLS diagnostic plots beside raw scatterplots.
- Label every panel as raw, phylogenetic residual, PGLS, or path-model evidence.

### P1.6 There is no one-command full evidence-chain gate

Impact: Tests can pass while the human-facing analysis notebook, path models, imports, and CellProfiler source contract are stale.

Evidence:

- `DATA_MANIFEST.yml` lists workflow steps but the CellProfiler import command does not include the top50 verified directory: [DATA_MANIFEST.yml](/home/jake/Projects/Desmognathus_TE/DATA_MANIFEST.yml:126).
- `scripts/run_full_pipeline.sh` completes the TE landscape pipeline, not the full size/genome/path/notebook evidence chain: [run_full_pipeline.sh](/home/jake/Projects/Desmognathus_TE/scripts/run_full_pipeline.sh:25).

Required fix:

- Add a separate `verify_size_evidence_chain` command that checks, in order:
  - top50 artifact row/species/hash contract;
  - CellProfiler import source and sidecar hashes;
  - master dataset tree match;
  - panel membership and genome status counts;
  - model input species after transforms;
  - model warning parsing;
  - notebook execution and required figure/table outputs.

## P2 Findings

### P2.1 Extend QC sensitivity to cell and nucleus morphology

Genome sensitivity is visible, but the same QC subsets should be shown for cell area, nucleus area, and N:C ratio. This will answer whether image quality changes only genome estimates or also the morphology conclusions.

### P2.2 Make manual overrides stricter and more auditable

Manual keep/repair rows should preserve their manual status but still display whether they fail strict auto QC, edge-touching, shape, trim, or overlap gates. Any override should carry a reason.

### P2.3 Add a duplicate/dropped tree-tip audit

Current genome/morphology species are tree-matched, but the pipeline should still emit retained/dropped/duplicate normalized tips for Python and R tree handling. This protects against later taxonomy or tree updates.

### P2.4 Make IOD background references path-aware

The top50 run reports 43 images and 15 images without cache background references in [summary.json](/home/jake/Projects/cellprofiler_test/output/runs/mixed_cellpose_yolo_full_dataset_v1_bgclean/verified_species_dataset_top50_latest/summary.json:33) and [summary.json](/home/jake/Projects/cellprofiler_test/output/runs/mixed_cellpose_yolo_full_dataset_v1_bgclean/verified_species_dataset_top50_latest/summary.json:36). Background matching should be keyed by stable source identity, not filename alone, and missing cache rows should be explicit sensitivity conditions.

## Scientific Interpretation Risks

### Genome size is not independent of the nuclear measurement

The genome estimate is calibrated from linked nucleus IOD using D. fuscus as reference. That can be a useful comparative index, but it is not the same evidentiary status as independent flow cytometry or assembly-based genome size. Any model claiming `genome -> nucleus -> cell` must be labeled as sensitivity/model-structure exploration, not independent causal proof.

### The best near-term authenticity upgrade is not more plotting

The highest-value upgrade is making every number traceable to exact selected pairs, exact source hashes, exact QC state, and exact tree/model inputs. After that, the notebook can be visually rich without being misleading.

## Recommended Execution Order

1. Lock the import contract: top50 exact-selected only, no silent fallback, source hash checks.
2. Clean final review semantics: no `maybe`, one canonical final decision, full invariant revalidation.
3. Fix genome sensitivity math: fixed-reference effect vs recalibrated-subset effect.
4. Fix support semantics: independent support units, one-image/one-specimen CI labeling, fallback severity.
5. Rebuild panels with explicit stable/minor/sensitivity-inclusive names.
6. Add pair-level provenance and register all top50 sidecars.
7. Surface model warnings and regenerate/stamp model stability outputs.
8. Upgrade the notebook with tree-trait plots, phylogenetic-corrected relationship diagnostics, and warning/support overlays.
9. Add the full evidence-chain verification command.

## Verification Performed During Audit

- Ran Desmognathus unit tests: `scripts/run_tests.sh` passed 7 tests.
- Ran CellProfiler top50/review tests: `uv run pytest tests/test_build_top50_linked_pair_worklist.py tests/test_freeze_species_grid_top50.py tests/test_build_verified_species_dataset.py tests/test_species_grid_review.py -q` passed 15 tests.
- Generated audit plots with `scripts/run_in_dusky.sh python plans/size-estimation-deep-audit/make_audit_plots.py`.
- Verified the three audit PNGs are nonblank by image dimensions and sampled pixel diversity.

## Plan Candidates To Write Next

Recommended first implementation plans:

1. `P0-top50-import-contract.md`
2. `P0-genome-sensitivity-semantics.md`
3. `P0-review-freeze-finality-and-invariants.md`
4. `P0-support-uncertainty-and-fallback-severity.md`
5. `P1-notebook-phylo-warning-evidence-surface.md`

These should be small, agent-executable plans with explicit acceptance tests and no broad refactors.
