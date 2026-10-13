# Independent mapping review, 2026-10-13

Reviewed `LOCAL_SEARCH.md`, `SOURCE_LINKS.md`, and the selected original-workbook cell extraction in `local_candidates.json`. The mapping CSVs and IOD source manifest were opened during the preceding source investigation. Original XLSX/ODS packages were opened by the separate local investigator, not independently reopened in this review.

## Conclusion

No reviewed field establishes a staining-run assignment for the retained unknown slides. The record joins image → physical glass/specimen, but the required join physical slide → staining event → reference slide remains missing. The documented author confirmation that standards were stained in the same experimental runs is relevant but does not identify each run's membership.

Recovered dates have distinct meanings:

- `scan_date` 20250521 and 20250523 explicitly identify acquisition dates for standards' image/glass records; they do not identify staining events.
- `Date` 5/18 on target image-ID records is undefined in the CSV. It might refer to imaging or another event; it cannot be assigned a staining meaning without clarification.
- `Date` 1/12/2023 in the historical inventory is also undefined. Its shared value does not prove simultaneous staining.
- `Verbatim Date` 7-Apr-37 and September 2023 field dates occur with specimen collection/locality records, not explicit Feulgen-run membership.
- File modification and Git commit dates describe digital records, not experimental dates.

`Fixed`, `Formalin`, `Fresh`, `Dried Blood`, mounting flags, image types, relative microscope folder names, and sequence-accession joins do not encode staining batches. The old taxonomic/fixation notes need separate reconciliation; this review does not overturn the author's designated-standard identities or infer a full preservation history.

## Minimum missing mapping fields

1. An unambiguous physical slide/glass identifier for every retained target image, including repeat images of the same slide or specimen.
2. A staining-event identifier or explicitly documented staining date/batch for that physical slide.
3. The reference slide/glass identifier(s) stained in that event, tied to specimen 32469 or 32470 and the relevant acquisition record. Distinguish multiple standard slides from repeat imaging of one slide.
4. Evidence for this assignment: worksheet/notebook row, contemporaneous run list, or a clearly labeled author reconstruction with its uncertainty.

Standards QC, acquisition comparability, the approximately 41% IOD difference, and the conditional assembly anchor remain independent calibration issues even if run membership is recovered. A complete batch map alone does not validate the absolute pg scale.

## Highest-value single author question after correct-account search

“Were all the retained brightfield slides stained together in one Feulgen run with both standards (32469 and 32470), or were there several runs? If there were several, which slide IDs were in each run and which standard was included?”

This first establishes whether a run-by-run map is needed at all and, if so, asks for the actual membership. A notebook or worksheet link can supply the answer without requiring a new recollection-based table. Do not lead with scan dates or presume that image numbering encodes batches.
