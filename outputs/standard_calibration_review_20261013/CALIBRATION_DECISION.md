# Calibration decision

**Do not replace the released denominator with the two unreviewed standard medians. Absolute genome-size calibration remains unresolved.** The designated slides have confirmed single-run membership and imaging/preparation context, and their stored measurements reproduce from retained pixels. Their candidate masks have visible problems and have not completed target-equivalent review. The assembly-derived scale remains conditional.

## What is now established

The author identified Processes 413 and 414, specimens 32469 and 32470, as 40× brightfield images without immersion oil, using the green channel and the same fixation protocol as the other slides. The earlier author answer placed both standards and all retained unknowns in one Feulgen staining run. These are author-attested experimental facts; exact staining dates, original camera response and lifetime museum preservation are not established by those answers.

The retained-mask audit reconstructed all 47 standard candidates and four historical target controls. IOD agrees with stored values at six-decimal precision and mask area agrees exactly. The target controls' eight historical pixel metrics agree within `1e-12`. Full raw/mask/outline crops provide a reproducible object census. These checks establish arithmetic and source consistency, not accepted nuclei or DNA stoichiometry.

| Current quantity | Value | Interpretation |
|---|---:|---|
| Process 413 candidates | 31; median IOD 594.078400 | Unreviewed mask measurements |
| Process 414 candidates | 16; median IOD 837.701417 | Unreviewed mask measurements |
| Directional median difference | 41.01%, using Process 413 as denominator | Reproducible discrepancy, cause unresolved |
| Released *D. fuscus* denominator | 943.889140125 | Equal-image mean of reviewed target-image medians from Processes 314/315/320/321 |
| Conditional reference factor | 16.1 Gbp / 0.978 Gbp per pg | Assembly-derived scale, not a validated specimen C-value |

## Why a standard mean is not yet usable

The whole-label computation can faithfully reproduce an incomplete or unsuitable mask. The visual review found a partial nuclear mask in Process 414, a nearly rectangular partial outline in another candidate, a debris-rich selected region in Process 413, and one mask touching the source tile edge. Padding clipping in a second gallery crop is distinct from a source-mask edge defect. No exclusions or mask repairs have been applied to the released data.

A genuine frozen quality template is available: 40 preserved *D. apalachicolae* objects and their canonical clipping/scaling parameters. Native extension helpers can reuse it without fitting a new template. They support fewer than 40 candidates, but their balance calculation then compares only assigned template slots. A wrapper must report full-template balance, unmatched slots, linked-pair eligibility and manual acceptance separately. Calling a small subset `common_support` would not establish a completed standard review.

The historical 30-kept-object guard belongs to the species-level review workflow. It is not an independently justified requirement of 30 nuclei per standard image. The current 31/16 starting pools nevertheless cannot each cover the entire 40-object template, and further eligibility review can only reduce their counts. The four existing target reference images themselves contribute 12, 8, 5 and 10 reviewed nuclei. Any accepted standard summary therefore needs an explicit aggregation and uncertainty policy, rather than an invented per-image minimum.

## Choice of next calibration state

| Option | Decision | Evidence needed or reason |
|---|---|---|
| Retain released IOD-derived, assembly-anchored values with the limitation stated | Retain as conditional research outputs | Current calculation is traceable; it does not certify absolute C-values. |
| Average the two existing standard candidate medians and recalibrate | Reject at this stage | Unreviewed masks, incomplete eligibility/balance evidence and unresolved between-slide difference. |
| Switch to a historical published *D. fuscus* C-value | Do not select one yet | Literature values do not establish identity and preparation suitability for these specimens. |
| Use accepted designated-standard measurements and a justified DNA-content reference | Required for a validated absolute scale | Review/repair complete nuclei using preserved source evidence, apply fixed technical-QC projection with honest shortfalls, and establish the reference value and calibration uncertainty. |

The primary-source anchor review confirms that the chosen 16.1-Gbp estimate is a pooled assembly from different animals and that repeat resolution limits its interpretation as true DNA content. A conversion between bp and pg changes units; it does not establish assembly completeness or a specimen-specific calibrant. Two standards of the same species also do not provide a multi-DNA-content photometric linearity test.

## Publication implication

This review resolves staining-run membership, retained standard imaging mode and shared preparation protocol. It does not resolve absolute calibration, supply accepted standard-object decisions, or make the manuscript submission-ready. Figures and released scientific tables were preserved. Proposed Methods corrections are provided separately so that the manuscript can describe the actual denominator and newly confirmed facts accurately.

The next measurement unit is a standards-only existing-mask linkage and boundary-review queue, followed by the frozen-template projection with separate subset/full-template balance and manual decisions. No segmentation rerun, new reference number, scientific model refit or figure regeneration is justified solely by the present arithmetic checks. The remaining manuscript-wide author and submission requirements remain in the prior submission-preparation verdict.

Evidence: [source inventory](SOURCE_INVENTORY.md), [pixel execution](PIXEL_AUDIT_EXECUTION.md), [visual review](VISUAL_REVIEW.md), [frozen interface](FROZEN_STANDARD_QC_INTERFACE.md), [anchor evidence](ANCHOR_EVIDENCE.md), [current updates](CURRENT_EVIDENCE_ADDENDUM.md), [proposed Methods wording](METHODS_UPDATE_PROPOSED.md), and [independent final verification](VERIFICATION.md).
