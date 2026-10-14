# Microscopy recovery workflow

This pilot helps review nuclear boundaries and test measurement sensitivity.
It preserves the released data and does not calculate or promote a new genome
size. The [measurement research](research/feulgen_recovery_measurement_evidence.md)
and [software comparison](research/feulgen_recovery_software_options.md) explain
the optical assumptions, alternative tools, and remaining acquisition evidence.

The default panel contains all 47 designated standard nuclei and the 35 released
*D. fuscus* target nuclei from Processes 314, 315, 320 and 321. All 82 start
unreviewed in this pilot. Historical target acceptance remains provenance, not
a new decision. The 40-object matching template is not a biological gold standard.

## Install and check

Use the repository's locked Pixi environments. The `recovery` environment adds
Rust 1.99.0 and OpenJDK 21.0.10 to the existing analysis and image-QC dependencies.
It leaves the other environments' package records unchanged.

```bash
pixi install --locked -e recovery
pixi run -e recovery --clean-env --locked --no-install recovery-setup
pixi run -e recovery --clean-env --locked --no-install check-env
pixi run -e recovery --clean-env --locked --no-install microscopy-recovery --help
```

Setup downloads checksum-pinned ImageJ 1.54p, compiles the Java ROI bridge and
Rust numerical reducer, and runs genuine geometry and photometry controls.
Tools and receipts live under `.pixi/recovery-tools/`. A failed setup is retained;
it is not repaired in place. A completed toolchain is checked before use. Changes
to source, lock, runtime, or qualification evidence invalidate it.

ImageJ is the contour editor and canonical ROI rasterizer. Bio-Formats is an
optional, separately installed source reader described below. This environment
does not install Fiji's plugin collection, operate the camera or run learned
segmentation.

## Read original microscope files

Install the optional checksum-pinned Bio-Formats 8.5.0 JAR through Pixi:

```bash
pixi run -e recovery --clean-env --locked --no-install source-reader-setup
```

The reader uses the locked Java runtime and keeps its JAR under
`.pixi/source-reader-tools/`. A failed cache is retained and refused rather than
repaired in place. Setup also accepts `--jar /path/to/bioformats_package.jar` as
a local seed; its published SHA-256 must match before it can be used.

Inspect native metadata without decoding image pixels:

```bash
pixi run -e recovery --clean-env --locked --no-install microscopy-source-reader metadata \
  --source '/run/media/jake/easystore/20250521/Slide glass_1321/Process_413.vsi' \
  --out .pixi/runs/microscopy-recovery/native-metadata-v1
```

Recover a small green-channel region spanning a retained tile boundary:

```bash
pixi run -e recovery --clean-env --locked --no-install microscopy-source-reader crop \
  --source '/run/media/jake/easystore/20250521/Slide glass_1322/Process_414.vsi' \
  --x 41975 --y 20416 --width 237 --height 179 \
  --out .pixi/runs/microscopy-recovery/native-context-v1
```

Coordinates are slide-global pixels with a top-left origin. The crop selects
flattened series 0, green sample 1, Z=0 and T=0. Its area is limited to one million
pixels. The reader rejects outputs with extra TIFF planes, the wrong dimensions,
RGB samples or a dtype other than uint8. It performs no autoscaling, casting or
pixel normalization. `-noflat` is deliberately absent from conversion: the
upstream converter would otherwise attempt additional pyramid resolutions.

Output directories must be fresh and below `.pixi/runs/`. The command retains
native stdout, stderr, exit status, tool/source pins and a receipt even if
conversion or validation fails. The reader opens source files for reading and
disables the upstream update check; it does not request reader caches beside the
source. A writable source mount is not described as read-only. VSI files are
hashed before and after reading; large ETS companions are not fully hashed by
this command.

Use recovered context as a separate diagnostic input. It does not replace packet
pixels, enlarge an existing ROI or transfer a prior decision. The
[Bio-Formats conversion documentation](https://bio-formats.readthedocs.io/en/stable/users/comlinetools/conversion.html)
describes the native region selection. The
[connected-source recovery record](recovery_work_20261014/NATIVE_SOURCE_RECOVERY.md)
documents the actual crops and their byte comparisons.

Recovery packets bind to their toolchain identity, including `pixi.toml`.
After registry changes, run explicit `recovery-setup` and create a fresh inventory
and packet in a new run directory. Preserve the old packet and its historical
qualification; do not rewrite it or transfer its decisions into the fresh packet.

## Build a new review packet

The original imagery and label masks must be available in the sibling
`cellprofiler_test` checkout. Absolute paths in frozen tables are checked against
those sources. A portable clone without the archive cannot fabricate this input.

Every output path must be fresh. These commands use a new run directory:

```bash
mkdir -p .pixi/runs/microscopy-recovery/pilot-v1
pixi run -e recovery --clean-env --locked --no-install microscopy-recovery inventory \
  --paper-root . --cell-root ../cellprofiler_test \
  --out .pixi/runs/microscopy-recovery/pilot-v1/inventory
pixi run -e recovery --clean-env --locked --no-install microscopy-recovery pack \
  --inventory .pixi/runs/microscopy-recovery/pilot-v1/inventory/inventory.json \
  --out .pixi/runs/microscopy-recovery/pilot-v1/packet
pixi run -e recovery --clean-env --locked --no-install microscopy-recovery validate \
  --packet .pixi/runs/microscopy-recovery/pilot-v1/packet
```

Crops are native uint8 green-channel pixels, with fixed-margin context and
original binary masks. There is no quantitative contrast change, interpolation,
denoising or sharpening. Validation checks actual crop pixels and mask membership
against their sources, as well as hashes, areas and recorded-background IOD.
Standard coordinates are slide-global; target coordinates are tile-local.

## Review together

Serve **only** the packet's `reviewer` folder. The packet root contains a private
source crosswalk and previous measurements and must not be served.

```bash
pixi run -e recovery --clean-env --locked --no-install --executable python \
  -m http.server 8765 --bind 127.0.0.1 \
  --directory .pixi/runs/microscopy-recovery/pilot-v1/packet/reviewer
```

Open `http://127.0.0.1:8765/`. Case labels are randomized pseudonyms. Toggle the
mask and inspect the raw image at native scale and with display zoom. Judge
whether the complete nucleus is visible, its boundary is supported, and nearby
cells, debris or a source edge make it uncertain. Metadata blinding cannot
conceal all visual differences between images.

If the available crop cannot show the complete boundary, mark it unsure. A new
packet can use a larger declared `pack --margin`; retain the old packet. A larger
crop cannot recover pixels beyond the original source edge.

Decisions default to blank. A named reviewer must explicitly choose `keep`,
`problem` or `unsure` and give reasons. Exports contain only authored events;
local browser storage is a draft cache. Independent reviewers use separate
exports. Conflicting decisions remain unresolved in the report.

Import an export into a fresh event directory:

```bash
pixi run -e recovery --clean-env --locked --no-install microscopy-recovery decision \
  --packet .pixi/runs/microscopy-recovery/pilot-v1/packet \
  --event /path/to/reviewer-export.json \
  --out .pixi/runs/microscopy-recovery/pilot-v1/reviewer-one-events
```

## Propose a contour correction

Prepare an editable copy for one case. Replace `CASE_ID` with its pseudonym:

```bash
pixi run -e recovery --clean-env --locked --no-install microscopy-recovery prepare-editor \
  --packet .pixi/runs/microscopy-recovery/pilot-v1/packet --case CASE_ID \
  --out .pixi/runs/microscopy-recovery/pilot-v1/editor-CASE_ID
```

The workspace contains copied raw pixels, an original-mask reference, a launch
argument list in `editor.json`, and a runnable command in `README.txt`. The launch
uses private preferences and disables GTK file dialogs because the tested locked
Java runtime hung while loading the host GTK dialog. Draw an integer rectangle,
polygon or freehand area on the raw copy, add the selection to ROI Manager and
rename it `nucleus`. Composite area ROIs are also supported. Oval, rotated,
rounded, spline and stack/channel-positioned ROIs are rejected. Optional
blank-field selections are named `background_001`, `background_002`, and so on.
Save a single selection as `nucleus.roi`, or save multiple selections as a ZIP.
The importer accepts the single native file by wrapping its exact bytes into a
canonical archive. Do not add a background selection merely to obtain a ZIP.
Do not save display contrast changes over the quantitative raw copy.

Complete the blank author, time and rationale fields in
`authorship.template.json`. Time must include a timezone. Its fixed bindings
identify the packet, case, crop and original state; do not change them.

```bash
pixi run -e recovery --clean-env --locked --no-install microscopy-recovery import-roi \
  --packet .pixi/runs/microscopy-recovery/pilot-v1/packet \
  --roi /path/to/rois.zip --sidecar /path/to/authorship.json \
  --out .pixi/runs/microscopy-recovery/pilot-v1/proposal-one
```

The importer uses ImageJ's own rasterization, including composite holes. It
rejects unsupported geometry, stale bindings, off-crop shapes and overlapping
nuclear/background selections. Imported contours are proposals. A proposal
does not create a human acceptance decision. A decision for a proposal must bind
to its new state hash; an original-mask decision does not transfer automatically.

To review a proposal in the browser, create a policy JSON:

```json
{
  "name": "reviewer-one-contour-pilot",
  "state_selection": {"CASE_ID": "PROPOSAL_STATE_SHA256"}
}
```

Build a fresh queue with that policy and the explicit proposal directories:

```bash
pixi run -e recovery --clean-env --locked --no-install microscopy-recovery review-queue \
  --packet .pixi/runs/microscopy-recovery/pilot-v1/packet \
  --proposals .pixi/runs/microscopy-recovery/pilot-v1/proposal-one \
  --policy /path/to/review-policy.json \
  --out .pixi/runs/microscopy-recovery/pilot-v1/proposal-review
```

Serve only the new queue's `reviewer` folder. Unselected cases still show their
original masks. The selected case shows the proposed mask, and exported decisions
bind to that proposal state. Pass its proposal directory to `decision` with
`--proposals` when importing the export. An existing original-state review remains
independent.

## Compare measurements

Create a baseline report with no new decisions or proposals:

```bash
pixi run -e recovery --clean-env --locked --no-install microscopy-recovery report \
  --packet .pixi/runs/microscopy-recovery/pilot-v1/packet \
  --out .pixi/runs/microscopy-recovery/pilot-v1/baseline-report
```

Later reports take explicit `--events` and `--proposals` directories. A policy
JSON with `name` and `state_selection` can select proposal states for named cases.
Keep each report and input set separate. Missing, unsure and conflicting reviews
stay visible; no last-writer-wins merge is used.

The Rust reducer measures unchanged source pixels as
`sum(log10(background / max(pixel, 1)))`. Negative optical-density terms remain.
The historical replay uses 255 for the standards and target-image backgrounds
255, 218 or 216 as recorded. A shared-255 calculation is labeled sensitivity.
Candidate-background means and medians are alternatives to inspect, not validated
incident-light references. Geometry comparisons include original/proposed area
and added/removed pixel counts. IOD changes are reported separately from those
counts. Dice, Jaccard and boundary-distance metrics in the research plan remain
future comparisons.

Reports contain no pg/1C conversion or pooled calibration denominator. Software
cannot restore clipped pixels, establish historical camera linearity, or replace
an independent DNA-content reference. Mask review is the first experiment;
paired reimaging, camera/illumination checks and reference validation follow if
the evidence warrants them.

## Verification and preservation

Portable backend and queue tests are part of `make test`. Image controls explicitly
skip in the default environment, which has neither Pillow nor tifffile. They run
in `image-qc` or `recovery`; native Java/Rust controls require `recovery`.
Run the complete recovery checks explicitly:

```bash
pixi run -e recovery --clean-env --locked --no-install --executable python \
  -m pytest tests/test_recovery_backend.py tests/test_microscopy_recovery.py \
  tests/test_recovery_review.py -q
```

Native-dependent skips in another environment are not native acceptance. GUI
checks use an isolated desktop and synthetic annotation fixtures. Infrastructure
tests do not author decisions or contour repairs in the real 82-case packet.

Raw crops, masks, ROI archives, downloaded JARs and private crosswalks remain in
ignored local runs. Commit small code, research and verification records. Preserve
all retained scientific manifests, Methods, figures, measurements and dependency
snapshots. See [REPRODUCIBILITY.md](REPRODUCIBILITY.md) for the broader boundary.

The [executed recovery verification](recovery_work_20261013/VERIFICATION.md)
records native producer handoffs, GUI tests, explicit skips and preservation.

The [AI pilot review](recovery_work_20261014/AI_VISUAL_REVIEW.md) records the
82-object first pass, wider-context findings and two unaccepted contour trials.
It keeps AI assessments separate from human decisions and the retained analysis.
