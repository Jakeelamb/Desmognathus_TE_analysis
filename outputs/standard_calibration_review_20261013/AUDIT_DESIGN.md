# Chosen retained-mask audit

Use [candidate B](design_b/DESIGN.md), the direct native helper design. The [independent design judgment](DESIGN_JUDGMENT.md) scored B 12/12 and A 11/12 under the same six criteria. The coordinator independently read both candidates and agrees. Native imports add two runtime dependencies but reuse the actual helper implementations without reconstructing their dependency graph through an AST loader.

Graft candidate A's fixed absolute target-pixel-metric tolerance of `1e-12` and complete-label crop coverage diagnostics. Retain both original standard coordinates and their explicitly normalized tile-local pixel-QC view. Pass historical target rows unchanged to the existing metric function and retain its actual crop/fallback behavior.

The executable input contract is exactly 47 inventoried standard nucleus rows and four target controls, chosen as the first canonical row for Processes 314, 315, 320 and 321. Sources are read-only. All results derive from full recorded mask labels and stored backgrounds. This is a reconstruction and pixel audit, not a manual nucleus acceptance process or a completed calibration.

Required output is one script, `audit_mask_pixels.py`, and exclusively created run directories containing object results, failure controls, source and runtime pins, and a systematic image gallery. Reject existing output directories before any write. Repeated successful runs must reproduce numerical results and gallery membership. Differences in clocks and fresh destination paths are receipt differences rather than scientific nondeterminism.

Required falsification covers an absent mask label, invalid unconverted standard coordinates, an altered source hash, a changed background, a removed positive-OD mask pixel, and an altered stored IOD. Original pixels and masks remain unchanged. IOD must match at recorded six-decimal precision, mask area must match exactly, and each target's eight stored pixel metrics must match within `1e-12`.

The new author [image confirmation](AUTHOR_IMAGE_CONFIRMATION.json) identifies standards 413/414 as 40× brightfield without immersion oil, using the green channel and the same fixation protocol as other slides. It does not supply camera-response tests, previous museum storage history, manual standard acceptance or a specimen-specific C-value.

No native alternative model family is available, so this is independent same-family design review. Full runtime and final evidence acceptance remain pending.
