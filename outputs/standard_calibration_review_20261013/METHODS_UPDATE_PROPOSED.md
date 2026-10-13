# Proposed Methods correction from the calibration review

This text records the newly confirmed preparation facts and describes the current released calculation. It does not implement a new standard denominator. The frozen submission-preparation package and online manuscript are unchanged.

For the specimen preparation and image acquisition paragraph:

> The two designated *D. fuscus* standard slides, from specimens 32469 and 32470 (Processes 413 and 414), and the retained unknown slides were stained together in one Feulgen run. The author identified the retained standard images as 40× brightfield images acquired without immersion oil, using the green channel, and confirmed that the standards underwent the same fixation protocol as the other slides. All images retained for morphometry and IOD were acquired by brightfield microscopy on an Olympus APX100 with a 40× objective. Images from the subsequent 40× oil-immersion acquisition were excluded because their quality was poor.

For the conditional genome-size calibration paragraph:

> The released conditional genome-size estimate was calculated as (target equal-image IOD / *D. fuscus* equal-image IOD) × 16.462167689 pg/1C. The computational denominator was the finalized *D. fuscus* species estimate from target Processes 314, 315, 320 and 321, rather than the two designated standard images. The scale factor was derived from the selected pooled *D. fuscus* assembly estimate of 16.1 Gbp divided by 0.978 Gbp per pg (Myers et al. 2025). Although the author confirmed a shared staining run, the designated standard nuclei have not yet completed the same eligibility, image-quality matching and acceptance process as the target nuclei. The existing candidate measurements from the two standard images therefore do not establish a validated calibration denominator. The released values remain IOD-derived, assembly-anchored species summaries rather than independently validated C-values. Reported intervals include sampled IOD variation but exclude uncertainty in the assembly anchor and calibration.

The shared fixation answer does not establish that every specimen, including previously supplied dried blood smears, had an identical lifetime preservation history. Exact reagent concentrations and execution times are not supplied by that answer. A generic staining protocol is not a record that each listed setting was used in this experimental run.

Evidence: [author image identification](AUTHOR_IMAGE_CONFIRMATION.json), [single-run confirmation](../calibration_search_20261013/AUTHOR_STAINING_CONFIRMATION.json), [53-row run map](../calibration_search_20261013/STAINING_RUN_MAP_AUTHOR_CONFIRMED.csv), [current denominator and candidate summaries](SOURCE_INVENTORY.md), and [primary-source anchor review](ANCHOR_EVIDENCE.md). Current Methods wording is preserved in `Publication/methods/MATERIALS_AND_METHODS_RESEARCHER_DRAFT.md`.
