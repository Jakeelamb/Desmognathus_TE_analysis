# Standard-mask visual review

The coordinator inspected all 47 standard candidates and four historical target controls using the raw and outline images in `pixel_audit_v1/gallery.html`. Five labeled contact sheets in `visual_review/` provided the complete census; raw and overlay crops for ordinals 3, 17, 32 and 45 were also inspected separately at their retained resolution. Contact sheets resize images for viewing but do not change intensity or contrast. The source-matched gallery crops are the measurement evidence.

This review flags visible problems. It does not supply accepted nucleus decisions, repair masks, or replace the target pipeline's linked-pair eligibility and frozen image-quality matching.

| Ordinal | Preserved review key | Observation | Consequence |
|---|---|---|---|
| 3 | Process 413, `tile_y000000_x032768`, label 2 | Selected dark region lies in a crowded, irregular, debris-rich area. An isolated intact nucleus is not established by this crop. | Requires exclusion or explicit boundary review before calibration. |
| 17 | Process 413, `tile_y008192_x045056`, label 4 | The object and outline reach the lower source-image edge. | Needs the same edge exclusion used for targets. |
| 32 | Process 414, `tile_y000000_x024576`, label 1 | A small polygon captures only part of the larger visible dark nuclear region. Stored area is 696 pixels and IOD is 140.478295. | Reproducing this mask's IOD does not show that a complete nucleus was measured. |
| 41 | Process 414, `tile_y020480_x040960`, label 1 | Gallery padding is clipped at the upper edge; the full-label audit retains the complete mask. | This display observation does not establish that the mask touches the source tile edge. Requires boundary review, without an automatic edge rejection based on the contact sheet. |
| 45 | Process 414, `tile_y028672_x065536`, label 1 | A narrow, nearly rectangular outline occupies part of the dark region, with surrounding overlap. Stored area is 1,692 pixels and IOD is 484.269010. | Requires boundary/overlap review; do not interpret as a complete accepted nucleus. |

The four historical target controls show comparatively complete outlines around visibly identifiable nuclei in these specific crops. Four examples do not establish that every target nucleus is comparable with the standard images.

The 47 candidate medians are reproducible. Their between-image difference has not been causally assigned to segmentation, fixation, imaging, staining, specimen identity or DNA content. The visible defects demonstrate why the unreviewed medians are unsuitable as a final standard denominator; they do not quantify how much of the discrepancy each defect contributes.

All flags refer to existing source pixels and masks. No biological observations were accepted or removed from the released datasets in this review.
