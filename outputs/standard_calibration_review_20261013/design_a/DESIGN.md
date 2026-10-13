# Candidate A: direct retained-mask audit

## Caller usage

Proposed one-file CLI, **design only**; no implementation is supplied:

```bash
export UV_PROJECT_ENVIRONMENT=/tmp/desmo-standard-qc-20261013/env
uv run --frozen --python 3.13 --with Pillow --with tifffile --with scipy --with scikit-image \
  audit_existing_standard_masks.py \
  --inventory outputs/standard_calibration_review_20261013/source_inventory.json \
  --targets analyses/02_morphology/data/nuclear_iod_objects.csv \
  --output-dir outputs/standard_calibration_review_20261013/mask_audit_01
```

Repeat with `mask_audit_02` to compare deterministic measurements, gallery membership, and falsification outcomes. Reusing `mask_audit_01` must fail before writing. Every run audits all 47 inventoried YOLO standard rows and the first canonical CSV row for each of Processes 314/315/320/321; no discretionary subsets or background overrides.

## Shape and signatures

One CLI owns selection, adaptation, evidence, and gallery publication. Private immutable shapes:

```text
RowPin = (table_path, table_sha256, csv_line, object_key, original_record)
StandardRow = (pin, raw_path, mask_path, label, GlobalPoint, TileOrigin, stored_iod, i_bg)
TargetRow = (pin, raw_path, mask_path, label, NativeTargetPoint, recorded_origin, stored_iod, i_bg)
AuditView = (pin, TileLocalPoint, zero_origin, selected_label, original_coordinates)
Result = (pin, area, reconstructed_iod, iod_delta, eight_pixel_metrics,
          metric_deltas, crop_bounds, coordinate_frame, outcome)

audit(inventory: Path, targets: Path, fresh_output: Path) -> AuditReceipt
load_pinned_helpers(source_pins: SourcePins) -> ExistingHelpers
standard_view(row: StandardRow) -> AuditView
measure(row: StandardRow | TargetRow, helpers: ExistingHelpers) -> Result
publish(results: tuple[Result, ...], fresh_output: ExclusiveDirectory) -> AuditReceipt
# Bodies: not implemented.
```

Pin standards by CSV order and full `filename::tile_name::label`, never label alone. The inventory and live standard table must agree exactly. Target controls are canonical lines 313/325/333/338, labels 20/2/41/34, respectively; selection is first-row order, not IOD. Preserve every original record and source/global field, including absent values.

Boundary checks require existing raw/mask files, aligned 2D shapes after existing grayscale conversion, integer labels, nonempty selected labels, stored area agreement, and finite stored background. Reconstruct IOD on complete TIFF label pixels using the existing `grayscale_uint8` and `measure_objects`; compare decimal values quantized to six places. Measure each unique standard tile/mask once, then join by label. No background estimation or replacement masks.

Standard pixel views subtract original tile offsets from global centroids and zero the view origins. Target baseline rows go unchanged into `compute_row_iod_metrics`, preserving its historical crop/fallback behavior and float32 percentile normalization. Validate standard crop label coverage against the complete mask; retain clipping/padding diagnostics. Metric parity tolerance is absolute `1e-12`, fixed before execution.

## Rationale and alternatives

The smallest viable seam reuses exact helper bodies. A fail-closed AST loader selects only named definitions from the SHA256-pinned YOLO and pixel-audit scripts, with explicit dependencies/constants; it compiles unchanged bodies, never their imports or `main`. Import the existing lightweight `microscopy_common` and `run_linked_pair_review` helpers normally. Retain whole-file/function hashes, filenames, and environment versions. Source drift or an unresolved dependency is an error, not a copied formula. This costs a small loader but avoids actual import-only OpenCV and matplotlib dependencies; the specified image environment lacks both. No Torch or segmentation backend is imported.

Full linkage/scoring followed by the legacy audit CLI hides more QC, but requires cell pairing/classifier state, creates unrelated eligibility decisions, and exposes destructive output directories and calibration summaries. It cannot preserve a frozen target template merely by passing standards. It loses for this bounded audit. Normal wholesale helper imports are simpler but fail in the authorized environment. Synthesis remains the coordinator's decision; this is an independent native same-family candidate.

## Evidence, outputs, and falsification

Read-only probes reconstructed all four target IOD values exactly and their eight retained metrics within `9.03e-17`; their native rows invoke whole-tile fallback. Zero-origin copies also matched these four controls. First standard rows reconstructed 607.722738/140.478295 exactly after coordinate adaptation; Torch was absent from loaded modules. All 47 rows remain unverified until implementation runs.

Publish `rows.csv`, `controls.json`, `provenance.json`, and a static gallery containing all 47 standards in pinned order, each with raw crop, binary mask, contour overlay, key, bounds, source links/hashes, stored/reconstructed IOD, and diagnostics. Crops use identical bounds; display normalization never feeds measurement.

Always run repeatable in-memory falsifications: missing/incorrect label, shifted mask, changed background, unconverted nonzero-origin centroid, and source hash mismatch. Also prove existing-output refusal. Record expected refusals/discrepancies, exit status, command, and receipt time. A failed target control retains deltas and byte/coordinate/dependency diagnostics; it blocks parity success without changing targets or tolerances.

One process exclusively creates a fresh nonsymlink direct child of this review root; resolved output/source overlap is forbidden. Hash sources before/after, write only that child, publish the completion marker last, and retain incomplete attempts. Concurrent invocations with the same name cannot both acquire ownership. Outcomes mean reconstruction/metric audit only: manual acceptance, full frozen-template QC, acquisition comparability, and calibration remain unqualified.

Source anchors: retained `SOURCE_INVENTORY.md`, `MEASUREMENT_FLOW.md`, `source_inventory.json`; sibling `run_yolo_tile_measurements.py:144,176`, `audit_iod_measurement_quality.py:97,137,170`, `run_linked_pair_review.py:385,456,585`; canonical table SHA256 `fffc91fa5f53a374711b7b71baf2292733e6185fa7d86bb02dc189a9cfc6d4bf`. Next implementation step: build the source-bound loader and reproduce the four pinned controls before producing the complete audit.
