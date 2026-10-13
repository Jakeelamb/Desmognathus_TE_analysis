# Specimen record associations

The supplied answer contains 97 MVZ records and 47 RAP records. The CSVs retain original species names, suggested names, comments and all source fields. `raw_fields_json` preserves the original field counts, including absent trailing empty fields.

Released specimen IDs have no MVZ/RAP namespace. A unique numeric correspondence is a candidate association, pending an explicit namespace record or author confirmation. Species agreement is reported separately and never changes the analysis taxonomy. The two reference standards have an author-designated MVZ association.

| Assay | Confirmed namespace | Unique numeric candidate | Unmatched | Ambiguous |
| --- | ---: | ---: | ---: | ---: |
| morphometry | 0 | 41 | 1 | 0 |
| nuclear_iod | 0 | 49 | 1 | 0 |
| reference_standard | 2 | 0 | 0 | 0 |

Unmatched analysis specimen IDs: 3463. Do not substitute RAP3462 for analysis specimen 3463.

MVZ `Date` is an ambiguous spreadsheet/lab date. Collection information comes from `Verbatim Date`; incomplete dates remain verbatim. RAP dates retain separate day, month and year fields. Collection records and historical `Dried Blood`/`Fixed` labels do not establish lifetime fixation and preservation.

Process_413/MVZ32469 and Process_414/MVZ32470 are the intended standards. Both source records report 7 Apr 1937, 6 mi E Calderwood, Tennessee. Released target D. fuscus measurements use other animals, including numeric IDs 3372/3373 associated provisionally with RAP3372/RAP3373. The standards are not members of the released 42/50-animal target panels.

Candidate RAP permit counts for the 50-animal IOD inventory:

| Permit verbatim | Candidate animals |
| --- | ---: |
| AGFC 053020233 | 2 |
| FWC LSSC-17-00031  | 1 |
| GADNR 124431310 | 9 |
| KDFWR SC2311207 | 5 |
| NCWRC 23-SC00374 | 16 |
| NCWRC SWC-24-000073 | 4 |
| WVDNR 2023.001 | 5 |

Exact species-name disagreements:

| Assay | Analysis ID | Analysis species | Source species | Suggested species |
| --- | --- | --- | --- | --- |
| morphometry | 109787 | D. orestes | D. ochrophaeus | D. orestes |
| morphometry | 91572 | D. orestes | D. ochrophaeus | D. orestes |
| nuclear_iod | 109787 | D. orestes | D. ochrophaeus | D. orestes |
| nuclear_iod | 169445 | D. orestes | D. ochrophaeus | D. orestes |
| nuclear_iod | 91572 | D. orestes | D. ochrophaeus | D. orestes |

RAP3518 retains the supplied `Contaminated???` comment. This is a source flag, not an automatic exclusion. MVZ suggested reidentifications remain in the source CSV and association table.

The receipt pins the source answer, metadata snapshot, released inventories and script. The script refuses existing output directories and writes no retained analysis artifact.
