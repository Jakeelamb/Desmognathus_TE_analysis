# Standard measurement flow review — 2026-10-13

**ISSUES: existing evidence does not show standards 413/414 passing the same target nuclear-IOD QC.** Measurement inputs and masks are available for a bounded remeasurement. No producer was rerun. Sources, Methods, figures, frozen decisions and existing outputs were not changed.

## Revision and source pins

- Paper checkout HEAD: `15dfb3c18755cc8659b626fa213087c62c6db02f`; extensive dirty work retained.
- CellProfiler sibling HEAD: `d7c37fdaf2869245964b6710ca76696306c2daba`.
- Paper archive HEAD: `7c235b808efc08117ffda246890e425e464fdb74`.
- Current opened source SHA256: `cellprofiler_test/scripts/run_yolo_tile_measurements.py` = `b92e52219b3953d139dca4183e808a44f4b467dc5449f4b26d292d777e7273d3`; `scripts/audit_iod_measurement_quality.py` = `41f205ed52ae353a3e836b8b0a1951b256444823ff99561b6189d8ba82cdc696`; archive `path_analysis/scripts/build_image_quality_matched_genome_iod_panel.py` = `4f5f18e6da535db83360e55bc9e4fc716bc7b188603c12f1b6e9d0735efc8df4`.
- Standard staging receipt mtime: 2026-07-13 14:46:36 -0700; YOLO nucleus receipt mtime: 2026-07-13 15:08:58 -0700. These are filesystem times, not independently recorded execution timestamps. Current measurement source mtime is 2026-07-14; source hash at historical execution is not recorded in that receipt.

## Existing standard evidence

Root: `/home/jake/Projects/cellprofiler_test/output/runs/life_history_pmount_reference_v1/`.

| Evidence | Process 413 / specimen 32469 | Process 414 / specimen 32470 |
|---|---:|---:|
| Kept staging tiles | 11 | 13 |
| YOLO measured nuclei | 31 | 16 |
| Median nuclear IOD, all YOLO objects | 594.078400 | 837.701417 |
| Median mean OD | 0.209082 | 0.257512 |
| Image background cache | 255 | 255 |

Read directly from `staging/source_inventory.csv`, `staging/tile_manifest.csv`, `staging/background_cache.csv`, `nucleus_measurements/all_measurements.csv`, and their `summary.json` files. Standard 414 median / 413 median is approximately 1.4101; this is an unreviewed, unmatched measurement comparison.

The standard nucleus receipt uses the same model path and shape thresholds as the target bgclean nucleus receipt: minimum area 500 pixels, circularity .55, solidity .92, aspect ratio 2.8. The currently available model `runs/segment/output/yolo_nucleus_training_final_area500_shape/weights/best.pt` hashes to `36d20b16f4bef51c6e3657f655471e7dedb6bba149a996ded1f2161c676f7f19`. Both receipts require an image background cache; standards have 24/24 tiles with `image_cache`, targets 324/324. Target receipt image type is `brightfield`; standards are `pmount`.

All 24 standard tile inputs exist. All 47 nucleus rows have existing `source_image_path` and `mask_path`. The 212 cell measurements and their 24/24 prediction receipt are preserved. Standard YOLO linkage/triage/quality-matched/manual-freeze outputs were not present under this run root; its `imagej_linkage` is an alternate branch, not the target YOLO branch.

The canonical 805-row `analyses/02_morphology/data/nuclear_iod_objects.csv` contains **zero** Process 413/414 rows. `provenance/iod_source_manifest.json` explicitly says these standards are not the released denominator and still require target-QC review/remeasurement. Final denominator remains D. fuscus target Processes 314,315,320,321.

Raw source inventory records:

- 413: `/run/media/jake/easystore/20250521/Slide glass_1321/Process_413.vsi`, SHA256 `2fe307f8cb2fafdb58b9cfa8dcfb76489f3f6969db30222dc910870cdb950179`; companion `_Process_413_/stack1/frame_t_0.ets`.
- 414: `/run/media/jake/easystore/20250521/Slide glass_1322/Process_414.vsi`, SHA256 `480f0182a9f41c065f1303c5ceaabd4d4ebc1d3865180c985a782eeeb0d99047`; companion `_Process_414_/stack1/frame_t_0.ets`.

Both VSI paths and both ETS paths are unavailable at their recorded locations. This review did not search other mounts; retained tiles permit mask-level remeasurement but not complete source restaging.

## Exact measurement and eligibility path

1. `cellprofiler_test/scripts/stage_vsi_green_tiles.py`: `bfconvert_full_command` extracts Bio-Formats series 0, channel index 1; `extract_tiles_from_ome` keeps tiles under the recorded `auto` filter. `microscopy_common.classify_tile` rejects mean >230 or <10 or SD <5; labels SD <15 or mean >200 yellow. `compute_image_background_from_ome` samples a grid, sums uint8 histograms, and selects the 95th percentile; this is the standard image cache. The full-image target cache exists at `mixed_cellpose_yolo_full_dataset_v1_bgclean/background_cache/brightfield_image_ibg.csv`, but its original producing command was not located in opened source.
2. `scripts/run_yolo_tile_measurements.py:grayscale_uint8` takes channel 0 for RGB arrays or planes, shifts uint16 right eight bits, otherwise clips/casts to uint8. `tile_background` requires the image cache; optional tile-p95 fallback changes the background regime and must not be enabled for parity.
3. `measure_objects` measures all label pixels from the segmentation mask: `pixel_vals = clip(gray[mask],1,None)`; `OD = log10(max(I_bg,1)/pixel_vals)`; `IOD = sum(OD)`; `mean_od = mean(OD)`. Negative pixel OD is retained, not clipped to zero. Area is pixels × 0.0144 μm². Centroids are stored in global image coordinates. Output IOD precision is six decimal places.
4. `apply_pair_review_repairs.py:measure_repaired_iod` uses the same logarithmic formula with stored `nuc_i_bg` on repaired masks. A standard mask alteration would require a retained new repair/review record rather than changing old masks.
5. `build_cell_nucleus_linkage_report.py:annotate_pair_qc` establishes physical one-to-one pairs, then within-image/species trims cell area, nucleus area, IOD, NC ratio to 5–95% (no trimming when fewer than ten eligible pairs). `keep_strict_core` requires all four robust log-z scores ≤2.5 and combined distance ≤4.5. Robust z = .6745×(value−median)/MAD; fewer than three or zero MAD gives zero scores. Therefore **eligibility already uses IOD**, although subsequent technical matching excludes it.
6. `score_new_linked_pairs_for_top50.py` enriches shape/focus with `run_linked_pair_review.compute_shape_qc_for_row`, fits the historical keep classifier, and writes `all_scored_masks.csv`; `shape_is_suspect` records jagged/ellipse/concavity reasons. A classifier recommendation is not manual acceptance.
7. Archive `build_balanced_genome_iod_sensitivity.py:prepare_candidate_frame`, called by `build_image_quality_matched_genome_iod_panel.py:build`, requires nonempty species/filename, `has_cell_match`, `physical_pair_ok`, `one_to_one_cell`, `keep_mask_pair`, `keep_strict_core`, no `cell_edge_touch`, positive cell/nucleus area and IOD, overlap label count ≤1 and fraction ≥.98, no discard/nucleus_only/maybe decision, and no suspect shape if that column exists.
8. `audit_iod_measurement_quality.py:compute_row_iod_metrics`: erode mask one pixel for interior; outer ring is dilation six minus dilation one; use first raw channel scaled by its crop 99.9th percentile; ring median minus inner median is local contrast. Gaussian gradient magnitude σ=1 over outer boundary gives gradient; edge sharpness = boundary gradient / max(abs(contrast),1e−6). Ring noise is Q3−Q1. `raw_gray01` casts to float32 before testing uint8, so its uint8 /255 branch is currently unreachable; preserve the actual code behavior for parity.
9. Archive `add_matching_features`: `log1p(edge_sharpness)` and `log1p(ring_IQR / max(ring_median,1e−6))`, finite values required. `robust_standardize`: historical 1st/99th percentile clipping, median and IQR scaling. `match_panel`: medoid species D. apalachicolae, image-balanced 40-object template, Hungarian assignment. Common support requires maximum absolute standardized mean difference ≤.10 and KS distance ≤.25. No IOD, darkness, contrast magnitude, area, or cell size is used as a matching feature.
10. `load_review_registry` excludes `problem` keys from matching; final freeze retains only explicit `keep` keys in common support. `reviewed_panel_freeze_ready` initially requires ≥30 keeps/species and the expected common-support species set. Later `review_finalization.json` adds reviewed Path24 extensions; release is 805 nuclei/24 species and deliberately includes documented lower support/balance failures. This is not a universal rule that all final species passed common support.

## Coordinate contract requiring runtime check

Observed source/data distinction: standard `measure_objects` rows have global centroids and nonzero tile origins. `build_cell_nucleus_linkage_report.link_by_tile` subtracts origin for sampling (lines 413–415) but retains global output centroid columns (477–478). The standard run's only paired CSV is `imagej_linkage/linked_nucleus_pairs.csv.gz`, e.g. Process 413 x=18061, origin x=16384; this is not the target YOLO pairing branch.

Target repaired rows can use mask-local coordinates: `apply_pair_review_repairs.measure_binary_metrics` adds only local mask slice offsets, and `apply_repairs_to_tile`/`main` overwrite repaired centroids without adding original tile origins. The opened target source first row, Process 312, has x=1834.4, y=811.13 while tile offsets are 4096,8192. Parent independently observed canonical Process 393 with local centroid and zero stored origins. No explicit centroid/origin normalization was located in `prepare_candidate_frame`, matching, or linkage output.

`audit_iod_measurement_quality.row_local_masks` sets raw origin to zero for a tile artifact but uses the row's origins for the mask. If the initial mask crop contains no label it falls back to the full mask and raw tile. This fallback can recover a repaired/local target row even with nonzero stored origins. It does not establish correct sampling for global-centroid standard rows where the initial mask crop does contain the label.

For a pixel-QC **copy** of a global standard linked row, subtract original x/y tile origins from nuclear centroids and set tile origins to zero, keeping original global coordinates separately. That aligns both arrays with the existing helper contract; validate one known label at nonzero origin against direct full-tile extraction. This is an adapter/schema issue, not evidence that canonical target QC is defective. Do not silently compare incorrectly located crops or replace historical target metrics with a corrected implementation without separate evidence.

## Archived-to-canonical mapping

Historical `path_analysis/data/external/derived/image_quality_matched_genome_iod/` now exists under `/home/jake/Projects/Desmognathus_TE_archive/` for matched caches, templates, sensitivity and finalization tables. Do not recreate old active-paper paths.

- Historical `image_quality_matched_nuclei_frozen_reviewed.csv.gz` → active `analyses/02_morphology/provenance/iod_pre_extension_frozen_objects.csv.gz`: verified SHA256 `f22d014f9f395c5aadec8fe2538f6f651d1328884ef618e0baf0a97a1e00dcbe` matches source manifest.
- Historical root `image_quality_matched_all_species_decisions.csv` → `provenance/iod_review_decisions.csv`: verified `9f79b25be26df38bb74c5f28b6a643754d853c2ddbc9bdef04c9b49d6b4eb3bb`.
- Historical replacement decisions → `provenance/iod_replacement_decisions.csv`: verified `5ad3dae5a848e93c11f1a9534ece30e18c6587f3c7dde35d7f68724f5d6b51f5`.
- Historical finalized reviewed nuclear table → canonical `data/nuclear_iod_objects.csv` (current compact release, not a byte-identical compressed-file mapping).
- Target candidate source remains live in sibling `pair_review_postrepair_strict/triage/all_scored_masks.csv`; verified manifest hash `dc3eb40d556999a34637c7d19dc0f94985af74532594fb04ca4681aafca62cf8`.

## Bounded commands and missing seam

Run from `/home/jake/Projects/Desmognathus_TE`; always use a new output directory because the measurement CLIs clear their output contents. These commands are reviewed interfaces, **not executed**, and do not by themselves satisfy full target QC.

```bash
export UV_PROJECT_ENVIRONMENT=/tmp/desmo-release-audit-20261013/env
CELL=/home/jake/Projects/cellprofiler_test
STD="$CELL/output/runs/life_history_pmount_reference_v1"
OUT=/home/jake/Projects/Desmognathus_TE/outputs/standard_calibration_review_20261013/remeasurement_v1

# Existing-mask first path: reuse preserved cells/nuclei; establish YOLO linkage.
uv run --frozen --no-sync --python 3.13 "$CELL/scripts/build_cell_nucleus_linkage_report.py" \
  --cell-csv "$STD/cell_measurements/all_measurements.csv" \
  --nucleus-csv "$STD/nucleus_measurements/all_measurements.csv" \
  --image-type pmount --output-dir "$OUT/linkage"
uv run --frozen --no-sync --python 3.13 "$CELL/scripts/score_new_linked_pairs_for_top50.py" \
  --pairs-csv "$OUT/linkage/linked_nucleus_pairs.csv.gz" \
  --training-scored-csv "$CELL/output/runs/mixed_cellpose_yolo_full_dataset_v1_bgclean/pair_review_postrepair_strict/triage/all_scored_masks.csv" \
  --selection mask_pair --output-dir "$OUT/scored"
```

The actual bounded nucleus-producing CLI, if fresh segmentation is needed, is:

```bash
uv run --frozen --no-sync --python 3.13 "$CELL/scripts/run_yolo_tile_measurements.py" \
  --manifest "$STD/staging/tile_manifest.csv" \
  --model "$CELL/runs/segment/output/yolo_nucleus_training_final_area500_shape/weights/best.pt" \
  --output-dir "$OUT/nucleus_measurements" --object-kind nucleus --image-type pmount \
  --background-cache "$STD/staging/background_cache.csv" \
  --patch-size 1024 --stride 1024 --imgsz 1024 --batch 8 --conf .10 --iou .50 --max-det 512 \
  --min-mask-area 500 --min-circularity .55 --min-solidity .92 --max-aspect-ratio 2.8 --device 0
```

It consumes exactly 24 preserved tiles and regenerates YOLO masks and `all_measurements.csv`, `tile_manifest.csv`, `summary.json`. Prefer a mask-only adapter calling `grayscale_uint8` and `measure_objects` once per existing tile/mask, or `measure_repaired_iod` per accepted binary mask, before any segmentation producer rerun. There is **no existing mask-only nucleus CLI** in the inspected scripts.

There is also **no complete existing CLI that adds standards to the exact frozen target matching template**. The archived matching CLI takes `--candidates-csv`, `--output-dir`, `--target-per-species`, `--force-pixel-qc`, repeatable `--decisions-csv`, but `match_panel` re-estimates scaling and template from its supplied pool. Passing only standards, or relabeling them as one/two artificial species, creates a different QC comparison. Passing them as D. fuscus among targets risks selecting target D. fuscus rather than designated standards.

The missing small adapter must: preserve process/specimen identity; recompute IOD from existing masks; call the exact strict eligibility and pixel functions; project onto the **frozen** source-manifest clipping/scaling and D. apalachicolae template; report each standard separately and any shortfall; retain acceptance review keys/decisions; aggregate within-image medians and equal-image/run-matched summaries only after acceptance. Do not derive a new standard calibration from unreviewed 47-row measurements.

Environment preflight via the original `/tmp/desmo-release-audit-20261013/env` uv environment found numpy/pandas but no scipy, scikit-image, tifffile, Pillow, torch, ultralytics, joblib or cv2. Thus the image/QC commands are currently blocked in that environment until a reproducible image-capable environment is provided. No packages were installed by this investigator. Parent reports a separate `/tmp/desmo-standard-qc-20261013/env` is being prepared; its readiness was not checked here. `uv --python 3.13` is the actual CLI syntax.

## Acceptance gaps

- **PASS (availability only):** retained 24 raw tiles, 47 nucleus rows and masks; production model exists and matches pinned production hash; candidate/frozen/decision hashes match declared source provenance.
- **ISSUES (same QC):** no YOLO standard strict-eligibility, frozen-template balance, or accepted-object review evidence; only 31 and 16 total nuclei before further rejection; cannot assume 40 or ≥30 accepted per standard.
- **ISSUES (measurement comparability):** pmount vs brightfield image type; staging/acquisition/background equivalence and linearity not established; global-standard versus local-target coordinate contract needs a genuine runtime fixture.
- **ISSUES (calibration):** exact target→staining run→standard mapping absent; ~41% unreviewed median difference unresolved; pooled-assembly-derived 16.462167689 pg/1C anchor is conditional and its uncertainty is not included.

First pickup: provision the approved image-capable uv environment, validate local-coordinate adapter handling on one preserved nonzero-origin standard tile against its mask, then measure/review the two standard sets into a fresh output directory with the fixed target template. This review provides no new absolute-C-value or publication-readiness acceptance.
