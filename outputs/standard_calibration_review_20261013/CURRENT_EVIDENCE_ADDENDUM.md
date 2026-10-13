# Updates to the earlier investigator snapshots

`MEASUREMENT_FLOW.md` is a retained snapshot of its investigator's bounded review. The following statements supersede its earlier gap list without rewriting that historical investigation.

- Target-to-run membership is now author-confirmed: all retained unknowns and both designated standards belonged to one staining run. See `../calibration_search_20261013/AUTHOR_STAINING_CONFIRMATION.json` and the 53-row `STAINING_RUN_MAP_AUTHOR_CONFIRMED.csv`. The calendar date was not recovered.
- The author identified Process 413 and 414 as 40× brightfield, without oil, using the green channel, and the same fixation protocol as the other slides. Metadata label `pmount` is not an unresolved objective/illumination-mode discrepancy. See `AUTHOR_IMAGE_CONFIRMATION.json`.
- A reproducible image-capable uv invocation has executed the retained-mask pixel audit. Missing image packages in the earlier investigator environment are no longer an execution blocker for that bounded audit. See `PIXEL_AUDIT_EXECUTION.md` when finalized.
- Standard coordinate normalization was tested on real retained pixels and masks, while historical target rows were passed unchanged. The reconstruction and target metric-parity results are in `pixel_audit_v1/` and `pixel_audit_v2/`.
- The archived extension script supplies native frozen-scaling and shortfall-aware matching functions, and its 40-row template is preserved. The complete standards-only eligibility/review wrapper remains unimplemented. See `FROZEN_STANDARD_QC_INTERFACE.md`; shortfall balance against assigned slots must not be mistaken for balance against the full template.

Photometric linearity, exact fixation execution details, previous specimen preservation, accepted standard nuclei, frozen-template quality matching and an independently validated standard C-value remain separate questions. Availability and arithmetic reconstruction do not resolve them.
