# Local calibration-record search, 2026-10-13

Local spreadsheet originals were accessible and read without modification. Several relevant specimen workbooks were located, but none of the inspected originals supplies an unknown-slide-to-staining-run-to-standard mapping. The recovered dates must not be substituted for staining dates.

## Strongest original records

| Source and exact cells | Observed record | What it establishes |
| --- | --- | --- |
| [Desmognathus_slides_pyron_mueller.ods](/home/jake/PhD_Research/Desmognathus_slides_pyron_mueller.ods), Sheet1 I179:I180 | 32469; 32470 | The designated standards' museum specimen IDs occur in the old lab inventory. |
| Same workbook, Sheet1 D179:D180, A179:A180, C179:C180, E179:E180, G179:G180 | `Formalin`; `conanti`; `???`; `Not convinced it isn't monticola`; `D. fuscus` | Archival preservation and identity notes conflict with a simple uniformly ethanol-fixed D. fuscus interpretation. These are unresolved old notes, not grounds to override the current author's standard assignment. |
| Same workbook, Sheet1 F179:F180 and M179:M180 | `1/12/2023`; `7-Apr-37` | F is headed `Date` (F96), M is headed `Verbatim Date` (M96). Neither is headed staining date. The latter is in the specimen locality/collection metadata block. |
| Same workbook, Sheet2 A84:C85, E84:E85, G84:G85 | `conanti`; `???`; `Not convinced it isn't monticola`; `D. fuscus`; 32469/32470 | A second sheet preserves the same identity concern. A1=`Photo Notes`, B1=`RAP Species`, C1=`RAP Notes`, E1=`Species`, G1=`MVZ:ID`. |
| [Berkely_Salamander_Lab_Work_MVZ_RAP_2.xlsx](/home/jake/PhD_Research/Berkely_Salamander_Lab_Work_MVZ_RAP_2.xlsx), sheet `Berkely_Salamander_Lab_Work_MVZ`, A84:G85 | A=`conanti`, B=`???`, C=`Not convinced it isn't monticola`, E=`D. fuscus`, F=`D. fuscus`, G=32469/32470 | The workbook preserves source labels, proposed labels and the unresolved note together. D84:D85 stores Excel date serial 44938 (2023-01-12); K84:K85 stores 13612 (1937-04-07), consistent with the ODS displayed dates. No staining-run column. |
| [Lamb Slide Specimens.xlsx](/home/jake/Projects/misc_desmog_scripts/Lamb%20Slide%20Specimens.xlsx), Sheet1 A4:D5 and K4:P5 | RAP3372 and RAP3373, genus Desmognathus, species fuscus, clade fuscus B; collection 24 September 2023; collectors RAP, MAS, JMH; permit WVDNR 2023.001; `Whole + Tissue` | Supports the specimen identities of the ordinary D. fuscus denominator animals. Field collection dates and specimen type do not establish staining or fixation history. |
| ODS above, Sheet1 A47:D49 and remaining RAP rows | Header `ID`, `Species`, `Preperation`; RAP3442/3437 and other Pyron specimens labelled `Fresh` | Provides a preparation distinction from the old museum block labelled `Formalin`. `Fresh` does not establish a full fixation timeline or exclude later fixation. |

The PhD_Research ODS modification time is 2024-12-11 18:40:50 UTC; the Berkeley XLSX is 2024-12-11 17:32:19 UTC; the misc_desmog_scripts Lamb workbook is 2025-02-19 20:59:10 UTC. These are filesystem modification times, not experimentally recorded dates or trustworthy authorship timestamps. No author identity was established from the workbooks. SHA-256 hashes and compact selected cell records are retained in `local_candidates.json`.

## Previously known image mappings confirmed

| Source | Exact record | Limits |
| --- | --- | --- |
| [master_image_metadata.csv](/home/jake/Projects/cellprofiler_test/master_image_metadata.csv:123), lines 123–124 | Process_413/414 → 32469/32470 → D. fuscus; preservation `Fixed`, mounting `dry`, image_type `pmount`, 0.12 µm/pixel | Image/specimen joins only. `Fixed` does not specify a fixative. |
| Same file, lines 178–180 | Processes 468/469 → 32469; Process 470 → 32470; mounting `pmount` | Later acquisition records for the same specimens. |
| [slide_glass_mapping.csv](/home/jake/Projects/cellprofiler_test/slide_glass_mapping.csv:83), lines 83–84 | Glass 1321/1322 → Processes 413/414, `scan_date`=20250521, original VSI names and relative drive paths | Explicit scan dates, not staining dates. The `drive_path` values are relative directory labels, not online Google Drive links. |
| [microscope_to_imageid_desmognathus.csv](/home/jake/Projects/cellprofiler_test/microscope_to_imageid_desmognathus.csv:127), lines 127–128 and 182–184 | 413/414 unmounted FALSE; 468–470 mounted TRUE, matching specimen IDs | Supports mounting/acquisition reconciliation. No run or standard pairing. |

The master metadata additionally associates both 32469/32470 rows with SRX20497025 / GCA_032353935.1, the accessions shown for D. planiceps in the supplied manuscript table. This is another archival metadata inconsistency; the accession column should not be used as evidence of these standards' taxonomic identity without reconciliation.

## Other candidates inspected

- [conserved_slides.xlsx](/home/jake/Projects/misc_desmog_scripts/Ectopic_recombination/Desmo_trait_evo/conserved_slides.xlsx), Sheet1, contains field series, species, clade, locality, collection date, collectors, permit and specimen fields. No staining-run fields in the inspected cells.
- [Desmognathus_slides_pyron_mueller.ods](/home/jake/Projects/misc_desmog_scripts/Desmognathus_slides_pyron_mueller.ods), Sheet1 rows 133–134 and Sheet2 rows 84–85, preserves the same 32469/32470 `Formalin` and identity concerns. Its layout differs from the PhD_Research copy. These copies are not byte-identical.
- [Lamb Slide Specimens.xlsx](/home/jake/Projects/Desmognathus_Results/Lamb%20Slide%20Specimens.xlsx) is another nonidentical copy. Sheet names and workbook XML text were inspected for staining, Feulgen, batch and the two standard IDs, with no matching text. This was a bounded text check, not a complete cell-by-cell semantic audit.
- `/home/jake/Documents/Attachments/Documents/untitled_0.ods`, Sheet1, is an assembly/accession table rather than a slide record.

No Google Drive/Sheets hyperlinks were found in the five fully decoded candidate workbook packages. The later Results XLSX copy was not claimed as a fully decoded hyperlink audit.

## Search boundary

Filename discovery used `rg --files --hidden` for `.xlsx`, `.xls`, `.ods`, `.csv`, `.tsv`, `.gsheet` across `/home/jake`, excluding caches, node_modules, Git internals, local tool state, virtual environments, generated outputs and the main financial repositories. That first broad pass reached 2,497 paths and included some remaining installed-library fixtures and a mounted financial dataset; these were filtered out of candidate review. No financial data contents were read. A targeted second inventory covered Documents, Downloads, Desktop, cellprofiler_test and Desmognathus_TE (398 paths), followed by PhD_Research, misc_desmog_scripts and Desmognathus_Results. The dedicated genome-size .xls files and plotted figure-data workbooks were not opened because their names indicate trait/analysis tables rather than staining records.

`findmnt` showed no additional mounted personal research disk. `/run/user/1000/gvfs` yielded no spreadsheet paths. The sabertooth SSH financial mount and encrypted Bernd knowledge volume were excluded from content inspection. Targeted text searches of the identified research folders found the existing image mapping files and manuscript notes, but no named staining-run or batch record. This bounded absence does not establish that the missing mapping is absent from online storage or other files.

Next useful online search terms from these originals: `Desmognathus_slides_pyron_mueller`, `Berkely_Salamander_Lab_Work_MVZ_RAP_2`, `Lamb Slide Specimens`, `32469`, `32470`, `1321`, `1322`, plus staining/Feulgen/batch. The original identity/preservation notes warrant an author query after the correct Drive account is connected.
