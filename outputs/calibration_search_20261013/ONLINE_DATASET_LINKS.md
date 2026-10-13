# Online dataset link extraction

Source: browser-exported workbook `Desmognathus dataset-52841d9f-57d0-46ab-8166-49d734618bbf.xlsx` (SHA-256 `79793579f5175f51db8a1b57dc31a7f839c39472a937aab30d1d627b57ae750c`), exported from [the online spreadsheet](https://docs.google.com/spreadsheets/d/1CaMYntCP1eDSa1kX297Y36b6hhb9Jcps_1SQqnE2M6E/edit).

The XLSX contains these tabs: `Sheet1`, `List of 38 Desmognathus species`, `dnaPipeTE results`, and `Links`. The second title is the 31-character XLSX form of the live Google Sheets tab `List of 38 Desmognathus species and Axolotl and TE` reported in the browser, consistent with Excel's sheet-name limit.

The only URL-bearing cells found are plain text on `Links`:

| Title cell | URL cell | Title | URL |
|---|---|---|---|
| `Links!A1` | `Links!B1` | phylogeny | https://pubmed.ncbi.nlm.nih.gov/35222955/ |
| `Links!A2` | `Links!B2` | SRA | https://www.ncbi.nlm.nih.gov/sra |
| `Links!A3` | `Links!B3` | NCBI | https://www.ncbi.nlm.nih.gov/datasets/genome/?taxon=43045 |
| `Links!A4` | `Links!B4` | Animal Genome DB | https://www.genomesize.com/results.php?page=1 |

These are general phylogeny, sequence, genome-assembly, and genome-size resources. None is identified in the workbook as a microscopy, blood, slide, staining, Feulgen, batch, pixel, or picogram source. The workbook has no Excel hyperlink objects or hyperlink relationships; its four URLs are literal cell text. Because this is an XLSX export, absent native Google Sheets link annotations cannot be ruled out from this file alone.

No cells matched microscopy, blood, slide, stain/staining, Feulgen, batch, metadata, pixel/pixels, picogram/picograms, specimen, `32469`, `32470`, process, or `413414`. Thus the export contains no record-level staining event linkage or target specimen/process identifiers.

The structural headers found in the candidate tabs are `Sheet1!A1` (`Link`), `Sheet1!C1` (`SRA Acession`), `Sheet1!D1:E1` (`Genome Accession`), `Sheet1!F1` (`dnaPipeTE`), `Sheet1!G1` (`Genome Sequence size from NCBI`), `Sheet1!H1` (`Contig N51`), `Sheet1!I1` (`C-value`), `Sheet1!J1` (`TE_curate_alpine`), `Sheet1!K1` (`TE_curate`), and `Sheet1!N1` (`RepeatModeler`). `List of 38 Desmognathus species` has `A1:G1` fields `Species`, `family`, `GenBank assembly accession`, `Contig N50 (kb)`, `Scaffold N50 (kb)`, `Scaffolds count`, and `GC percent`. No microscopy or staining field appears among these headers.

No external links were fetched. The compact cell-level extraction is in `ONLINE_DATASET_LINKS.json`.
