# Candidate B: direct native mask-helper audit

## Usage first

Proposed single new file: `outputs/standard_calibration_review_20261013/audit_mask_pixels.py`. Its only caller choice is a fresh output directory; its bounded input contract is the inventoried 47 YOLO nuclei and the first canonical row for each Process 314/315/320/321.

```bash
UV_PROJECT_ENVIRONMENT=/tmp/desmo-standard-qc-20261013/env \
uv run --frozen --python 3.13 \
  --with Pillow --with tifffile --with scipy --with scikit-image \
  --with matplotlib --with roifile \
  python outputs/standard_calibration_review_20261013/audit_mask_pixels.py \
  --out outputs/standard_calibration_review_20261013/pixel_audit_v1
```

Second call: identical command with `pixel_audit_v2`, comparing deterministic artifact hashes. Third call: repeat `pixel_audit_v1`; it must refuse before any write. Every invocation includes falsification controls. Exit 0 means reconstruction/parity checks passed; exit 1 retains issues; exit 2 means preflight or output-path refusal. These outcomes confer no manual acceptance, frozen-template QC, calibration, or publication approval.

## Data shape and signatures

All types and functions below belong to that one file; bodies are **not implemented**.

```python
StandardGlobal = (global_centroid_xy, tile_origin_xy)
TargetHistorical = (original_centroid_xy, original_origin_xy)
CoordinateEncoding = StandardGlobal | TargetHistorical
PinnedRow = (role, csv_sha256, csv_line, row_key, immutable_source_record,
             coordinate_encoding, raw_path_sha256, mask_path_sha256)
AuditResult = (pinned_row, mask_area, recomputed_iod, iod_delta,
               eight_pixel_metrics, stored_metric_deltas, issues)

def audit(out: FreshOutputDirectory) -> AuditReceipt:
    raise NotImplementedError
def inspect(row: PinnedRow, helpers: ExistingNativeHelpers) -> AuditResult:
    raise NotImplementedError
def verify_controls(rows: tuple[PinnedRow, ...], helpers) -> ControlResults:
    raise NotImplementedError
```

Private parsing validates one inventory-backed row set, then groups by raw/mask paths. Derive gallery order, row keys, counts and summaries from that set; retain table order, CSV line, Process/specimen identity and complete original records. Standards use `filename::tile_name::label`. Targets retain their existing `review_key` and original field values.

## Measurement ownership

Import existing native `apply_pair_review_repairs.grayscale_uint8` and `measure_repaired_iod`; the latter calls the former. Read the actual TIFF labels, select the recorded label over the full tile, and call `measure_repaired_iod` with a full-tile slice and stored background. Require exact formatted six-decimal IOD equality, exact area, matching dimensions and positive existing label; retain the unrounded delta. Never substitute backgrounds or omit negative OD pixels.

Call unchanged `audit_iod_measurement_quality.compute_row_iod_metrics`. For a standard-only copied row, project column names, subtract its original tile origin from global centroids and set helper origins to zero. Preserve both coordinate systems in output and validate the projected centroid against actual mask pixels. For targets, pass the unchanged historical row: their repaired local centroids with nonzero origins intentionally trigger the existing full-tile fallback. Converting targets to the standard crop regime changes metric normalization and cannot establish historical parity.

## Evidence and alternatives

Grounding: `SOURCE_INVENTORY.md`, `MEASUREMENT_FLOW.md`, `source_inventory.json`; paper HEAD `15dfb3c18755cc8659b626fa213087c62c6db02f`, sibling HEAD `d7c37fdaf2869245964b6710ca76696306c2daba`.

Read-only native-helper probe selected canonical lines 313/325/333/338, with key suffixes `314::tile_y024576_x024576.tiff::20`, `315::tile_y016384_x024576.tiff::2`, `320::tile_y004096_x081920.tiff::41`, `321::tile_y004096_x028672.tiff::34` (full filename prefixes retained). All four reconstructed IOD at six decimals; areas were 3117/3430/2960/2901; eight-metric maximum absolute delta was below `9.03e-17`; Torch was not loaded.

Accept two additional runtime dependencies: native helper imports transitively require matplotlib and roifile. The supplied four-extra-package environment alone fails import on matplotlib. Reject AST extraction: it adds dependency reconstruction and a second execution mechanism. Reject full linkage/legacy CLIs: they expose selection/model/matching stages, may clear output directories, and exceed this mask audit. Native direct helpers provide the deeper, smaller interface. Coordinator synthesis remains pending; this is an independent native same-model-family candidate.

## Protection and acceptance checks

Resolve and reject symlinked output ancestors, input/output ancestry overlap and existing destinations; create exclusively. Never call a legacy `main` or `clear_dir`. Hash every consumed source, helper and input before/after; pin interpreter/package versions, command, execution status and every row-to-file mapping in `receipt.json`.

Produce `objects.csv`, `controls.json`, `receipt.json`, and a static gallery showing all 47 standard raw crops beside mask/overlay crops, plus four target controls. Preserve uint8 display scale, annotate crop bounds, clipping, coordinates and pinned keys; no IOD-based filtering.

Deterministically reject an absent label, an inconsistent coordinate projection, and an altered input hash. Remove the first positive-OD foreground pixel from an in-memory mask copy: area and reconstruction must fail. Alter a copied stored IOD by `0.000010`: precision comparison must fail. Sources remain untouched. If any target fails, retain exact metric deltas, paths/hashes, encoding and fallback evidence; report reconstruction versus metric parity separately and never relax tolerances. First implementation step is the fresh-output guard and pinned-row parser.
