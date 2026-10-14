# AI review of the 82-object calibration pilot

Two models inspected every raw crop and original mask in the calibration pilot.
The coordinator then inspected every disagreement and every case either model
flagged as a visible defect. This is a diagnostic review of the pilot, not a
review of the entire microscopy dataset. It produced no changes to measurements,
source masks, Methods, figures or the retained release.

## Review and evidence

The reviewers were configured as `gpt-6.1-sol` (`pstack_deep`) and `gpt-6-luna`
(`pstack_helper`). They worked separately, without the other review, the private
source crosswalk or previous IOD values. Each inspected 21 contact sheets covering
all 82 objects. The sheets used unchanged uint8 crops beside the same crops with
a cyan mask overlay, enlarged threefold with nearest-neighbor display sampling.
No contrast adjustment or quantitative pixel change was applied.

The criterion was whether the mask enclosed the visible dark nuclear region.
Angular contours alone were insufficient for a defect call. Faint cell rims and
ambiguous halos were not automatically included as nuclear material. Folding,
overlap, weak contrast and missing context remained reasons for uncertainty.

| Joint assessment | Objects |
|---|---:|
| Both apparently supported the original mask | 28 |
| Both were uncertain | 34 |
| Both flagged a visible defect | 5 |
| Disagreed | 15 |
| Total | 82 |

Twelve objects were flagged as defective by at least one model. Agreement is
descriptive: it does not establish a biological boundary or camera response.
Qualitative confidence labels are not calibrated probabilities.

The coordinator's targeted inspection covered 20 objects: all 15 disagreements
plus the five jointly flagged defects. It identified five mask-repair candidates,
two objects requiring missing tile context, eight uncertain cases and five
apparently supported masks. A repair candidate still requires a complete visible
boundary before a replacement contour is drafted.

| Case | Coordinator finding | Next action |
|---|---|---|
| `N11f41a0edb7dac8d23bb` | Dark pointed lower tip extends below the mask | Draft a separate contour proposal |
| `Nf816b73e83bc3d3e7f61` | Small mask encloses only part of a larger dark region | Inspect wider context before drafting |
| `N823287fa03ce457343f0` | Flat lower edge crosses an extended dark region | Inspect wider context before drafting |
| `Ncf47a0f18f718c508252` | Separate mask island above-left of the main region | Check the island and main contour separately |
| `N7fce100ee144e5a90c39` | Contiguous upper-left dark lobe is outside the mask | Inspect wider context before drafting |
| `N815f697d43b9efed92c4` | Nuclear material is clipped at the top of a retained tile | Recover the north neighbor of the Process 414 tile |
| `N3291d5e0a6bc3fe2af23` | Nucleus reaches the bottom of a retained tile | Recover the south neighbor of the Process 413 tile |

The [coordinator findings](coordinator_findings.json) preserve reasons and next
actions for every targeted object. The original
[Sol](sol_assessment.json) and [Luna](luna_assessment.json) records preserve all
82 assessments each; disagreements were not replaced by a majority vote.

## Wider context and contour trials

The coordinator subsequently inspected wider source-pixel context for all twelve
flagged objects. Each display added 64 pixels on every available side; its
original packet region matched source raw and mask pixels exactly. The
[context receipt](wider_context_receipt.json) pins these displays and sources.
The [follow-up findings](context_findings.json) preserve the initial assessments
and record what the wider context changed: seven cases remained ambiguous or
overlapped, one retained an isolated mask-component problem, two required missing
neighboring tiles and two supported a contour trial.

Two integer-polygon ROIs were authored with genuine ImageJ 1.54p, imported through
the maintained recovery CLI and measured with its Rust reducer. Contours were
chosen from visible morphology before their IOD was calculated. They remain
AI-authored proposals in separate ignored storage, with zero acceptance events.

| Case | Original → proposed area, pixels | Original → proposed signed IOD | Follow-up |
|---|---:|---:|---|
| `N11f41a0edb7dac8d23bb` | 1,682 → 1,810 | 316.054 → 334.617 | Both model reviews found better enclosure of the lower tip; still unaccepted |
| `N7fce100ee144e5a90c39` | 842 → 1,136 | 221.321 → 289.120 | Model reviews disagreed about joining the upper-left body; held as ambiguous |

IOD uses the recorded background of 255 in both cases. These changes measure
contour sensitivity; an increase in IOD does not establish a more accurate
DNA-content estimate. The [trial record](contour_trials_summary.json) retains
integer vertices, source/state hashes, native producer receipts and measurement
results. Separate [Sol](sol_proposal_visual_review.json) and
[Luna](luna_proposal_visual_review.json) assessments used only the comparison
image, without IOD values. No original mask was overwritten or trial promoted
into the paper analysis.

## Recovered acquisition records and missing context

A read-only audit recovered Olympus metadata from retained Process 320/321 TIFFs:
DP23 #1 camera, nominal 40× objective, air, approximately 0.6 NA, 0.12 µm pixel
scale, and exposures of 45.467 and 44.315 ms. The actual TIFFs are single-channel
uint8, while their copied XML describes RGB and pyramid IFDs. These fields record
acquisition lineage; they do not prove that the exports retain native RGB samples
or a linear camera response. Process 314/315 rewritten green-plane XML retains
pixel scale without camera or exposure settings. Standard Process 413/414 staging
records preserve original VSI/ETS pointers, but no standard exposure settings
were recovered. The original source drive is currently absent.

Both clipped objects lie at internal tile seams within the recorded full
acquisitions. Their context may be recoverable from the original microscope
files; it is not available in the retained neighboring tiles checked locally.

| Case | Retained tile | Required neighbor |
|---|---|---|
| `N815f697d43b9efed92c4` | Process 414 `tile_y020480_x040960.tiff` | North: `tile_y016384_x040960.tiff` |
| `N3291d5e0a6bc3fe2af23` | Process 413 `tile_y008192_x045056.tiff` | South: `tile_y012288_x045056.tiff` |

Neither neighbor occurs in the retained staging manifest, and no other listed
tile overlaps the required continuation strip. Temporary full green OME exports
are absent. Bounded checks of related local image folders and expected project
mirror paths on the mounted backup recovered no substitute. This was not an
exhaustive filesystem search. No missing pixels were synthesized.

The audit also traced how retained backgrounds enter signed IOD measurement.
The current standard staging helper pools an 8×8 grid of image boxes and selects
the 95th-percentile intensity. Target cache values are preserved, but their
original generating receipt was not recovered. No new background value or
absolute genome-size calibration was selected.

## Verification and retained artifacts

[Verification](verification.json) checked 82 unique case/state pairs in each
review, all 21 display pages and all 165 source asset pins. The maintained recovery
validator passed pixel equality, mask equality and numerical replay for the
original packet. All 12,282 protected entries (434,825,677 bytes of regular files)
matched the pre-tooling baseline, with zero changes.

Packet SHA-256:
`fd9afdd0fc827d304553d31a003dc9cef6c46a0caaeed0680bbed106650baa94`.
The [display manifest](display_manifest.json) pins the original assets and sheets.
Raw crops, contact sheets and detailed private acquisition receipts remain in
ignored local run storage. Small review records are retained here.

The [independent verification](independent_initial_verification.json) passed 25
checks and detected five deliberate in-memory negative controls. It pins the
initial report before this wider-context and proposal follow-up; its coverage
does not extend to the later trial records.

These AI assessments are separate from authored human decisions. The initial
review created no replacement ROI. The subsequent trials created two proposals
and zero acceptance events, and did not promote the conditional genome-size
calibration. See
[MICROSCOPY_RECOVERY.md](../MICROSCOPY_RECOVERY.md) for proposal and measurement
comparison commands.
