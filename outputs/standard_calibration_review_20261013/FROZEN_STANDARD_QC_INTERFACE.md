# Frozen standard QC interface — 2026-10-13

**A genuine frozen-template seam exists.** Reuse archive `extend_finalized_microscopy_reviewers.py:apply_frozen_scaling` (101–112) and `match_species_to_template` (134–180), with the original 40 D. apalachicolae rows from `image_quality_matched_nuclei_all_species_sensitivity.csv.gz`. Do not invoke its whole CLI: `NEW_SPECIES` is hard-coded to aeneus/orestes/wrighti (58), `build_iod_extension` filters to them (342), and it composes historical review pages. A small standards-only wrapper is required; no new estimator, template or classifier is necessary.

Scope: source reads and pandas verification only; this report is the sole new artifact from this investigator. No segmentation, matching execution, source mutation, or manual acceptance occurred. The requested `cellprofiler_test/archive` directory does not exist; the actual archive is `/home/jake/Projects/Desmognathus_TE_archive`. Read MEASUREMENT_FLOW, AUDIT_DESIGN and AUTHOR_IMAGE_CONFIRMATION once. Parent confirms same staining run; author confirmation records 40x brightfield without oil, green channel and shared fixation protocol. Those are not outstanding generic acquisition-type questions. Photometric linearity and absolute C-value remain unresolved.

## Preserved inputs, tested live

Archive HEAD `7c235b808efc08117ffda246890e425e464fdb74`; paper/sibling HEAD pins remain those in MEASUREMENT_FLOW. Current source hashes govern reuse even if their working files differ from HEAD.

- Canonical `analyses/02_morphology/provenance/iod_source_manifest.json`, lines 105–117, SHA256 `bb692d5cea7c84e50a97164c5a05476853546cad1717ba2edf9c1d0838067c35`: exact scalar parameters, not a fitted sklearn/pickle scaler.
- Archive `path_analysis/data/external/derived/image_quality_matched_genome_iod/image_quality_matched_nuclei_all_species_sensitivity.csv.gz`, SHA256 `d07c05e989ea821efe7ec34fd3344c70076ebbdfc415adc917340487906e5a53`, exactly matches canonical manifest `matched_panel_sha256`. Select `species == 'D. apalachicolae'`, sort `quality_template_id`. Observed 40 unique IDs 1–40; 20 rows Process_356_raw_green.ome.tiff, 20 Process_357_raw_green.ome.tiff. Preserved columns include both `match_log_*`, both `z_match_log_*`, template IDs, review keys and identity. This is the template array; no separate `.npy` or binary scaler is needed.
- Archive `manifest.json` has exactly equal `feature_scaling` to canonical manifest. Reapplying canonical scalars to preserved template features reproduces stored z columns with maximum absolute differences `1.1102230246251565e-16` and `9.540979117872439e-17`. These are pandas round-trip differences, not a newly fit scaler.
- Live standard measurement identities: exactly 47 nuclei/212 cells, all D. fuscus, `image_type=pmount`; filenames `Process_413_green.ome.tiff`/specimen32469 and `Process_414_green.ome.tiff`/specimen32470. Preserve this historical image_type to read the data despite author-confirmed brightfield imaging.
- `extend_finalized_microscopy_reviewers.py` SHA256 `0563b69ea8be8e911cfe49af85821c7d9d42f3995f1d37c1485205f1e724dd3c`; mtime `2026-07-14 11:47:47.913560275 -0700` (filesystem date, not execution provenance).
- Archive `build_image_quality_matched_genome_iod_panel.py` SHA256 `4f5f18e6da535db83360e55bc9e4fc716bc7b188603c12f1b6e9d0735efc8df4`; `build_balanced_genome_iod_sensitivity.py` SHA256 `a36b7044ee5358b2c8d1a4a22cabd08dbd47c2668d416cd6825daddc4399584a`.
- Sibling `build_cell_nucleus_linkage_report.py` SHA256 `0756768118793bc607ef3b91175fcf8bf351d01e1a5d53e212256941b74fa5be`; `run_linked_pair_review.py` SHA256 `fab5563394846c10e6d31409844ffe10fdf6eb51db218a8efa5b484b545568ab`.

Read-only probe executed from paper root with `UV_PROJECT_ENVIRONMENT=/tmp/desmo-standard-qc-20261013/env uv run --frozen --no-sync --python 3.13 python -`. Template/schema/scaler verification passed. Importing the genuine extension module failed with `ModuleNotFoundError: No module named 'scipy'`; therefore native matching/rejection self-controls were not executed. No package changes were made. Its top-level imports also load archive viewer/morphology helpers; provision their genuine dependencies rather than extracting functions with AST or silently rewriting them.

Coordinator subsequently reports successful pixel-audit v1 reconstruction for all47 standards/four target controls and target eight-metric parity below1e-12. `pixel_audit_v1/objects.csv` exists (path verified here; receipt and rows not re-audited by this investigator). Its technical metrics do not establish linked strict eligibility or manual acceptance. Coordinator's image review found partial masks, edge cases and debris/overlap, so retain diagnostics as diagnostics. Image dependencies are available through an ephemeral `uv run --frozen --python 3.13 --with Pillow --with tifffile --with scipy --with scikit-image --with matplotlib --with roifile` environment, rather than the base `--no-sync` environment probed above. No further execution is authorized for this slice.

## Minimal linked-pair contract

Existing-mask producer: sibling `build_cell_nucleus_linkage_report.py` consumes existing standard `cell_measurements/all_measurements.csv` and `nucleus_measurements/all_measurements.csv`, with `--image-type pmount --window-radius 1`. Its output `linked_nucleus_pairs.csv.gz` carries these fields. Keep original linked/global rows for physical linkage and shape/focus computation; supply the separately validated pixel-local view for nuclear pixel QC.

| Stage | Required fields and invariants |
|---|---|
| Identity allowlist | Exact (`filename`,`specimen_id`) pairs above; `species='D. fuscus'`; immutable `standard_process_id` 413/414 and original `nucleus_object_id`. Never infer the standard group from species alone. Reject any other process/specimen, duplicated nucleus object ID, or rows not present in the original47 identity inventory. |
| Native eligibility (`prepare_candidate_frame`,241–297) | Nonempty `species`,`filename`,`review_key`; `nuc_iod`,`nuc_area_um2`,`cell_area_um2`; `nucleus_mask_overlap_label_count`,`nucleus_mask_best_overlap_fraction`; all six flags `has_cell_match`,`physical_pair_ok`,`one_to_one_cell`,`keep_mask_pair`,`keep_strict_core`,`cell_edge_touch`. Require valid numeric/bool input rather than relying on permissive bool coercion. |
| Shape/focus | `shape_is_suspect` must exist and be measured for this adapter, even though native eligibility makes it optional. Native `run_linked_pair_review.compute_shape_qc_for_row` supplies cell/nucleus smoothness, ellipse IoU, solidity, shape reason flags and raw focus fields. Retain `mask_label_id`, both label paths, cell/nucleus centroids and both tile origins needed by that helper. Validate actual crops/mask labels; absent data can produce NaNs without automatically making `shape_is_suspect` true. |
| Nuclear pixel QC | `nucleus_source_image_path`,`nucleus_mask_path`,`nucleus_label`,`nuc_area_px`,`nuc_centroid_x`,`nuc_centroid_y`,`tile_x0`,`tile_y0`; original `nuc_i_bg`, IOD and full-label pixel audit linkage. `audit_iod_measurement_quality.row_local_masks`137–167 / `compute_row_iod_metrics`170 supply metrics. Required matching fields are `iod_edge_sharpness`,`iod_ring_noise_iqr`,`iod_ring_intensity_median`. Preserve all pixel metrics/crop diagnostics from the chosen mask audit. |
| Review decisions | Native registry schema `species`,`review_key`,`decision`; valid values `keep`,`problem`,`unsure`, repeated files use last key. Native key helper `review_key_for_row`284 returns `nucleus_object_id` when available, otherwise filename/tile/nucleus/cell composite. Standard IDs remain explicit additional fields. Exclude `problem` before matching; no `keep` means unreviewed, not accepted. |

`prepare_candidate_frame` fills absent numeric auxiliary columns with NaN and calls `add_qc_scores`; these clarity/shape percentile scores are unused by the technical matcher. Do not call `apply_quality_floor` or `select_balanced_by_image`, which implement a different selection path. Native shape computation can be called directly; fitting `score_new_linked_pairs_for_top50`'s keep classifier is unnecessary for the frozen matching interface. If that existing scoring CLI is reused for operational convenience, its predictions remain triage, not manual decisions.

## Exact gates and frozen projection

- Native linkage `physical_pair_ok`494–501: matched cell, positive areas and strictly cell area > nucleus area. `annotate_pair_qc`545–604: one nucleus per cell; within (`filename`,`species`) eligibility computes 5th–95th percentile gates on cell area, nucleus area, IOD and NC area ratio. `trim_mask`234 skips trim below10; robust z241 is `.6745*(log value - median)/MAD`, zero below3 or zero MAD. Strict requires every |z|≤2.5 and Euclidean combined≤4.5. These historical eligibility calculations use IOD; frozen technical matching does not. No preserved per-image trim/MAD cutoff artifact is supplied by the canonical scaler. Reuse the same eligibility rule on standard linked rows; do not claim frozen numerical eligibility cutoffs from the target pool.
- Native strict candidate gate241–294: required positive IOD/areas, all five positive linkage/keep flags, no cell-edge-touch, overlap-label count≤1, best-overlap fraction≥.98, no `discard`,`nucleus_only`,`maybe`, and no suspect shape. Empty strict pool raises `ValueError`; adapter should retain a zero-candidate status/count and skip pixel/matcher calls for that group.
- Shape threshold source91–101: cell smoothness≥.90, ellipseIoU≥.78, solidity≥.93; nucleus≥.92/.82/.95. Native reason generation840–849 only flags finite values below thresholds: adapter must reject missing/nonfinite required shape/focus evidence rather than treating it as clean.
- `quality_match.add_matching_features`74–93: `log1p(max(edge_sharpness,0))`; `log1p(max(ring_IQR,0)/max(ring_median,1e-6))`; drops nonfinite features. Record counts/reasons of dropped rows.
- `extension.apply_frozen_scaling`101: `z=(clip(feature,clip_low,clip_high)-median)/scale_iqr`. For edge use `(0.0918853819511443,0.14351591028167543,0.1162730078971146,0.01408501476837741)`; noise `(0.05683348912867326,0.2093502121396634,0.1094338192883331,0.03819577458050219)`. Never call `robust_standardize`, `quality_medoid_species`, `build_common_template` or `match_panel` for the standards.
- `extension.match_species_to_template`134–180: Euclidean cdist in those two preserved z features; rectangular Hungarian assignment; preserves original template IDs; at most40 unique candidates. No hard pair-distance cutoff or candidate image quota. Empty template/candidates returns zero with `insufficient_candidates`. Otherwise status `common_support` iff maximum |mean z difference|≤.10 AND maximum two-sample KS≤.25.

**Important shortfall semantics:** this extension helper compares balance against only the assigned template IDs (157–166), not all40 when candidates<40. It can label an arbitrarily small subset `common_support`; this says nothing about40-object coverage or review readiness. Retain that native result as `assigned_subset_balance`, and separately calculate/report both feature mean/KS differences against all40 preserved rows. Never shrink/re-estimate the template to16, relabel target D. fuscus, or silently pool the standards. Original `optimal_quality_match`163–184 in the main matching module assumes every template slot has an assignment and would fail on fewer40 candidates; use the extension helper.

The starting maxima31/16 mean standard413 has at least9 missing40 slots, standard414 at least24, before any QC. Original reviewed freeze guard278–292 requires≥30 keeps per common-support species; apply/report this cardinality separately per designated standard group if evaluating equivalent review sufficiency. Standard414 cannot meet30 from its current16. A pooled group with47 objects is an optional explicitly separate diagnostic; it cannot certify either standard independently or establish a run-matched denominator. No relaxation or replacement with target Processes314/315/320/321 is permitted.

## Smallest handoff

First produce the standards-only linkage in a fresh directory (existing native CLI; not executed here):

```bash
CELL=/home/jake/Projects/cellprofiler_test
STD="$CELL/output/runs/life_history_pmount_reference_v1"
OUT=/home/jake/Projects/Desmognathus_TE/outputs/standard_calibration_review_20261013/frozen_standard_qc_v1
UV_PROJECT_ENVIRONMENT=/tmp/desmo-standard-qc-20261013/env \
uv run --frozen --no-sync --python 3.13 "$CELL/scripts/build_cell_nucleus_linkage_report.py" \
  --cell-csv "$STD/cell_measurements/all_measurements.csv" \
  --nucleus-csv "$STD/nucleus_measurements/all_measurements.csv" \
  --image-type pmount --window-radius 1 --output-dir "$OUT/linkage"
```

Proposed wrapper (not created; these arguments are an implementable contract, not an existing CLI):

```bash
UV_PROJECT_ENVIRONMENT=/tmp/desmo-standard-qc-20261013/env \
uv run --frozen --no-sync --python 3.13 apply_frozen_standard_qc.py \
  --linked-pairs-csv "$OUT/linkage/linked_nucleus_pairs.csv.gz" \
  --standard-inventory-csv "$STD/nucleus_measurements/all_measurements.csv" \
  --pixel-audit-results /home/jake/Projects/Desmognathus_TE/outputs/standard_calibration_review_20261013/pixel_audit_v1/objects.csv \
  --manifest /home/jake/Projects/Desmognathus_TE/analyses/02_morphology/provenance/iod_source_manifest.json \
  --matched-panel /home/jake/Projects/Desmognathus_TE_archive/path_analysis/data/external/derived/image_quality_matched_genome_iod/image_quality_matched_nuclei_all_species_sensitivity.csv.gz \
  --standard 413:32469 --standard 414:32470 --target-count 40 --minimum-reviewed 30 \
  --output-dir "$OUT/projection"
# Optional repeatable --decisions-csv containing standards' actual manual keep/problem/unsure registry.
```

`pixel-audit-results` must carry unique original nucleus identity, validated IOD/full-label area, pixel metrics, local-coordinate view and source pins, as defined by the parent's audit implementation; check actual field names and successful receipt before integration. Preserve originals, use a copy for local pixel input, and enforce a one-to-one identity merge. Compute native shape/focus on verified paired crops; merge those metrics into the linked rows, call `prepare_candidate_frame(...,max_overlap_label_count=1,min_overlap_fraction=.98)`, exclude problem decisions, attach pixel metrics, `add_matching_features`, frozen scaling, then loop explicit413/414 groups through `match_species_to_template` against the same40 rows. The helper ignores species for assignment, so keep the real D. fuscus label and group outside it. Replace the displayed base-environment `uv` invocations with the coordinator's ephemeral `--with` dependency invocation above when executing; no environment installation or execution was performed here.

Required outputs: per-object gate/rejection history and review status; unchanged40-row template copy and scalar manifest with pins; selected rows with template IDs/distances; each standard's input/linked/strict/finite/matched/reviewed counts, missing template IDs,40- and30-shortfalls, native assigned-subset balance and full40-template balance; receipt retaining exact argv/interpreter/dependency versions/source hashes. Review acceptance must remain distinct from balance and cardinality. Aggregate accepted IOD only after manual review, retain each standard's median separately, and do not replace the released denominator.

Boundary/rejection controls required before using the new wrapper: unexpected target Process314/315/320/321 row; mismatched specimen/process; duplicate nucleus ID; wrong template hash or missing/duplicate template ID; missing/scaler nonfinite or nonpositive scale; invalid/nonfinite required metric; missing mask label/path/crop audit; zero-candidate group;16-candidate shortfall; unreviewed row falsely counted as keep. Preserve source checksums before/after and reject existing destination directories. Native helper tests remain unqualified until scipy/import dependencies are available.
