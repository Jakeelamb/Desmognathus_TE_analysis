# Microscopy source-link tracing, 2026-10-13

Read-only source investigation. No source data, Methods, image, mapping, or cloud file was edited. This report is the only output owned by this investigation.

## Finding

No Google Drive or Google Sheets URL/file identifier was found in the inspected current project provenance, analysis notebooks, documentation, or source scripts. Four local CSV names are confirmed production metadata inputs and useful search terms for the correct Drive account:

- `master_image_metadata.csv`
- `microscope_to_imageid_desmognathus.csv`
- `slide_species_mapping.csv`
- `slide_glass_mapping.csv`

The source code reads local CSVs; inspected code does not establish which cloud workbook generated those CSVs. Names and local joins are observations. A Google Sheets export origin remains an inference until an original workbook/export record is located.

## Exact sources and data flow

- `/home/jake/Projects/cellprofiler_test/docs/ACTIVE_WORKFLOW_INVENTORY.md:12-15` lists all four CSVs as active data inputs. Its canonical-production notice is dated 2026-07-14.
- `/home/jake/Projects/cellprofiler_test/data/README.md:12-23` identifies metadata aliases and explicitly says these are symlinks, not independent copies. Observed `data/metadata/` symlinks resolve to the four project-root CSVs.
- `/home/jake/Projects/cellprofiler_test/src/cellprofiler_tools/pipeline_runs.py:70`, `load_metadata_lookup()`, reads `data/metadata/master_image_metadata.csv` using `csv.DictReader`, keyed by `(image_type, filename)`. `load_manifest_rows()` joins that lookup to manifest rows, with nonempty manifest values overriding metadata values. This is an input-consumption path, not an import from an online spreadsheet.
- `/home/jake/Projects/Desmognathus_TE/analyses/02_morphology/provenance/iod_source_manifest.json:73-103` records designated standards 413/414, specimens 32469/32470, and a separate current denominator comprising target Processes 314/315/320/321. The anchor decision is dated 2026-08-12. Exact staining-run mapping remains listed as unresolved.

## Mapping facts retained in the local CSVs

| Source | Exact rows | Observation |
|---|---|---|
| `master_image_metadata.csv` | 4-5, 10-11 | Target Processes 314/315/320/321 map to specimens 3372/3373/3401/3402; species D. fuscus; image type brightfield; mounting dry; preservation Dried Blood. |
| `microscope_to_imageid_desmognathus.csv` | 5-6, 11-12 | The same four target pairs are labeled date 5/18, Pmount FALSE, Dried Blood. The column is simply Date; it must not be relabeled as a staining date. |
| `master_image_metadata.csv` | 123-124 | Processes 413/414 map to standards 32469/32470; image_type pmount but mounting dry; preservation Fixed. These fields have different meanings and should not be conflated. |
| `microscope_to_imageid_desmognathus.csv` | 127-128 | Processes 413/414 map to 32469/32470; Pmount FALSE; Fixed; Date blank. |
| `slide_glass_mapping.csv` | 83-84 | Glass 1321/1322, Processes 413/414, scan_date 20250521, original files Process_413.vsi/Process_414.vsi. |
| `slide_glass_mapping.csv` | 138-140 | Glass 1417/1417/1418, Processes 468/469/470, scan_date 20250523. |
| `slide_species_mapping.csv` | 174-176 | Processes 468/469 are specimen 32469, and 470 is specimen 32470. These are repeat imaging records of the two standards, not three additional animals. |
| `microscope_to_imageid_desmognathus.csv` | 182-184 | Processes 468/469/470 are Pmount TRUE and Fixed, with Date blank. |

`slide_glass_mapping.csv` has a column named `drive_path`; observed values are relative microscope folders such as `20250521/Slide glass_1321`. These are neither Google URLs nor verified cloud folder IDs. Scan dates and relative folder paths do not establish staining batches.

## Source revisions and pins

Cellprofiler checkout HEAD observed: `d7c37fdaf2869245964b6710ca76696306c2daba`, commit date 2026-03-13T12:33:48-07:00. Checkout has existing modifications and untracked files. The four mappings' latest shared recorded commit was `8b8bc5a6d94e9ace5ea0eee9d1d019ce4ca3def1`, 2026-03-09T20:08:29-07:00; this is a repository commit date, not a spreadsheet creation or staining date.

Observed SHA256 values:

```
bb164b21dfd5226a5e35879c4095ceb78474837b1510b05354fd165de4593629  master_image_metadata.csv
b987be5eafe79ef1492c8353b1ecbc605c50883e081c27ddb0660bb43a089595  microscope_to_imageid_desmognathus.csv
3d2c1ecefbd8b54d65c218fb0910b43c44a3a3c6660ac89e4e5b6eedf066d5dc  slide_species_mapping.csv
aef95aa3cb3a964cb19ff86a6bc5112f1057c747466ee220397478b238617d99  slide_glass_mapping.csv
```

## Search scope and limitations

Targeted `rg` searches covered current Desmognathus_TE docs, analyses, scripts, and notebooks, and cellprofiler_test docs, scripts, source, and data metadata references. Searches included ignored source files but excluded environments, package caches, generated figure trees, and credential/browser-profile data. A targeted `git grep` of historical Desmognathus commit `643b5c09420e67adaa54712231fdf14351ecc11b` (2026-08-11) covered the historical research-review notebooks and morphology analysis; no Google/Sheets link was found there either. These negative results do not prove that the cloud workbook does not exist.

The archive location cited in the cellprofiler documentation, `/home/jake/Projects/_archive/Desmognathus_cellprofiler_nonproduction_20260714`, was unavailable at search time (No such file or directory). That source slice was not searched.

Root owns cloud discovery and the separate local spreadsheet inventory. No cloud read was performed by this investigator. The parent reports the connected Drive profile is a different person's account, so empty searches from that account do not establish absence from Jake's Drive. Correct-account search and opening the newly discovered original local workbooks are the next steps.
