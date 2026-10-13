# Retained-mask pixel audit execution — 2026-10-13

**Reconstruction and historical metric parity PASS. Calibration remains unqualified.** Both bounded runs exited 0. All 47 inventoried standard rows and the first canonical row for each target Process 314/315/320/321 reproduced whole-label IOD at six-decimal precision and mask area exactly. The four targets' 32 stored pixel metrics reproduced with maximum absolute delta **9.020562075079397e-17**, below the fixed `1e-12` tolerance. No old measurement, mask, raw tile, Methods, figure or frozen review record changed.

## Artifacts and exact checks

- Implementation: [audit_mask_pixels.py](audit_mask_pixels.py), SHA256 `18839f8f40a0059cd785bfad627f0cf6dea788deeb2f67572f32c253b7f19ded`.
- First run: [receipt](pixel_audit_v1/receipt.json), [full command log](pixel_audit_v1/command.log), [outer clocks/exit](pixel_audit_v1/command_receipt.json), [exact command](pixel_audit_v1/command.txt), [objects CSV](pixel_audit_v1/objects.csv), [full original rows, mappings and deltas](pixel_audit_v1/objects.json), [falsifications](pixel_audit_v1/controls.json), [all-object gallery](pixel_audit_v1/gallery.html).
- Deterministic repeat: [receipt](pixel_audit_v2/receipt.json), [full log](pixel_audit_v2/command.log), [outer clocks/exit](pixel_audit_v2/command_receipt.json), [exact command](pixel_audit_v2/command.txt), [gallery](pixel_audit_v2/gallery.html).
- Final artifact checks: [24-check evidence](pixel_audit_v2/validation.json), [full log](pixel_audit_v2/validation.log), [exact replayable command](pixel_audit_v2/validation_command.txt), [clocks/exit](pixel_audit_v2/validation_command_receipt.json). Exit 0; 24/24 checks passed.
- Existing destination refusal: [full log](pixel_audit_v2/existing_output_refusal.log), [command, clocks, exit, complete before/after file pins](pixel_audit_v2/existing_output_refusal_receipt.json). Exit **2**; all **161** first-run files retained identical hashes, bytes and mtimes.

Both runs use exactly:

```bash
UV_PROJECT_ENVIRONMENT=/tmp/desmo-standard-qc-20261013/env \
uv run --frozen --python 3.13 \
  --with Pillow --with tifffile --with scipy --with scikit-image \
  --with matplotlib --with roifile \
  python outputs/standard_calibration_review_20261013/audit_mask_pixels.py \
  --out outputs/standard_calibration_review_20261013/pixel_audit_v1
```

The repeat substitutes only `pixel_audit_v2`. The refusal repeats the existing v1 destination. First-run command clock: `2026-10-13T17:55:58.871419583Z`–`2026-10-13T17:56:17.132915858Z`, exit 0. Second-run clock: `2026-10-13T17:57:03.192865925Z`–`2026-10-13T17:57:21.484528260Z`, exit 0. Clocks are execution receipts rather than performance claims.

The final checks independently verified the exact inventoried standard records, unchanged first canonical target records (CSV lines 313/325/333/338), all 51 reconstructions, all 32 historical metrics, and complete label coverage in both metric and gallery crops. Target areas are **3117/3430/2960/2901** pixels. All **153** saved RAW/MASK/OVERLAY images match their actual retained source pixels and recorded labels. All **157** numerical/control/gallery artifacts have identical hashes across the two runs. Source/helper/script pins and resolved package versions also agree. Receipt clocks, destinations and interpreter temporary paths are allowed execution differences.

The six controls all rejected the copied alteration: absent label, unconverted standard global coordinates, changed expected source hash, changed stored background, one removed positive-OD foreground pixel, and changed stored IOD. The removed-pixel control retained area delta **-1** and IOD delta **-0.1512680076535844**; the changed background (255→254) retained IOD delta **-5.09208835335005**. Symlink ancestors, input/output ancestry overlap in both directions and a missing destination parent also rejected. All controls leave original pixels and rows untouched.

## Runtime and source pins

Paper HEAD: `15dfb3c18755cc8659b626fa213087c62c6db02f`. Sibling HEAD: `d7c37fdaf2869245964b6710ca76696306c2daba`. Dirty work was preserved. Each run pins **62** consumed input/source files before and after computation, including **seven** native helper modules; no pin changed. Native `grayscale_uint8`, `measure_repaired_iod` and `compute_row_iod_metrics` were called directly. No Torch was loaded; no segmentation/model/HPC producer ran.

Resolved runtime: CPython **3.13.11**; numpy **2.5.2**; pandas **3.0.5**; Pillow **12.3.0**; tifffile **2026.9.20**; scipy **1.18.1**; scikit-image **0.26.0**; matplotlib **3.11.2**; roifile **2026.9.22**. Exact interpreter paths and all file hashes are in each receipt. Native core helper SHA256 pins: `apply_pair_review_repairs.py` = `ae6e82575b63b50b237c5c85c99e375c14f1636171d905beecf8a76d0eae3327`; `audit_iod_measurement_quality.py` = `41f205ed52ae353a3e836b8b0a1951b256444823ff99561b6189d8ba82cdc696`.

Standards retain their original global coordinates and tile origins alongside their explicit local pixel-QC copy; historical target rows were passed unchanged and all four use the existing full-tile fallback. The gallery uses whole-label bounds plus 64 pixels independently of the metric helper's crop. RAW images use their original uint8 values without contrast changes. MASK images show only the recorded label; OVERLAY adds a red inner boundary. All original coordinates, bounds, clipping, complete-label coverage and file mappings remain reviewable in `objects.json`.

## Review flags and limits

There are 31 Process 413 and 16 Process 414 candidates. All appear in original CSV order, without IOD selection. Three standard gallery crops clip at tile bounds:

| Gallery ordinal | Recorded key | Observation |
|---:|---|---|
| 17 | `Process_413_green.ome.tiff::tile_y008192_x045056.tiff::4` | **Mask itself touches the lower tile edge**; full recorded label is covered, but the image does not establish a complete biological nucleus. |
| 37 | `Process_414_green.ome.tiff::tile_y000000_x028672.tiff::1` | Padding clipped; recorded mask does not touch tile edge. |
| 41 | `Process_414_green.ome.tiff::tile_y020480_x040960.tiff::1` | Padding clipped; 593-pixel recorded mask does not touch tile edge. |

Systematic visual review and independent final acceptance remain with the coordinator/reviewer. This audit proves faithful reconstruction of retained candidates, not acceptance under the frozen target QC criteria. It does not resolve YOLO/ImageJ selection disagreement, acquisition response/linearity, prior museum storage history, specimen-specific C-value or unavailable original full VSI/ETS recovery. Author confirmation covers 40× brightfield without oil, green channel and shared fixation protocol; it supplies none of those additional guarantees. No new denominator, genome-size estimate or calibration PASS is authorized by these numerical results.

## Retained checker failure

The first post-run artifact-checker invocation **failed, exit 1**, because its independent reader attempted to open a target PNG with `tifffile.imread`. The audit correctly used the existing native extension-aware loader. The failed checker [log](pixel_audit_v2/validation_failed_v1.log), [exact command](pixel_audit_v2/validation_failed_v1_command.txt) and [clocks/exit](pixel_audit_v2/validation_failed_v1_command_receipt.json) are preserved. Only the checker reader was corrected; retained audit artifacts were validated again without another scientific execution. Final artifact validation passed, exit 0.

Owned changes are confined to this report, `audit_mask_pixels.py`, and the fresh `pixel_audit_v1/` and `pixel_audit_v2/` evidence directories. No commit was made.
