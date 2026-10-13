# Calibration-record discovery, 2026-10-13

The bounded local and online record searches are complete. Jake's existing signed-in Chrome provided access to the research Drive account. Relevant inventories, the GBE manuscript, and staining/imaging protocol documents were recovered. Jake's direct answer, **“One staining run with both standards”**, resolves membership of the retained brightfield slides and designated standards by author attestation.

The new [author-confirmed run map](STAINING_RUN_MAP_AUTHOR_CONFIRMED.csv) contains 51 retained image rows and two designated standard rows. Its run identifier is an internal label; the actual staining date and original experimental run identifier remain unknown. This record does not validate the standards or the absolute genome-size calibration.

## Cloud access

The connected Google Drive profile returned `gsuhr02@gmail.com` (Garrett Suhr), rather than the research account. No related file contents were fetched from that account. Searches for `Desmognathus`, `salamander`, `Feulgen`, and `staining` returned no matches there; those results do not establish absence from Jake's Drive.

Jake identified the research account as `jakelamb74@gmail.com` and offered browser sign-in. Google rejected sign-in in the separate agent desktop `desmo-drive-records-3d76c8d7fe74`, displaying `This browser or app may not be secure` and advising a different browser. The error was independently read from the visible page. The failed desktop was stopped; no bypass was attempted.

Jake then explicitly authorized use of his existing signed-in browser. Native cua-driver window discovery returned no windows; the browser extension exposed the existing Chrome session. The visible Drive account was verified as `jakelamb74@gmail.com`. The separate Drive connector remains linked to the other account; no account change was needed for browser access.

Read-only browser exports of five research files were retained, with source URLs and SHA256 pins in [ONLINE_SOURCE_PINS.json](ONLINE_SOURCE_PINS.json). Searches in the verified research account used `Desmognathus`, `slide`, and `microscope`. Their findings establish these recovered sources, not an exhaustive absence of other records. No cloud files were edited.

## Recovered local leads

- `/home/jake/PhD_Research/Desmognathus_slides_pyron_mueller.ods`
- `/home/jake/PhD_Research/Berkely_Salamander_Lab_Work_MVZ_RAP_2.xlsx`
- `/home/jake/Projects/misc_desmog_scripts/Lamb Slide Specimens.xlsx`
- Four image/slide mapping CSVs in `/home/jake/Projects/cellprofiler_test`; exact rows and hashes are retained in the source-link report.

The original inventories contain unresolved historical species and preservation notes for designated standards 32469 and 32470. These notes have been preserved without changing the current standard assignment. Collection dates, scan dates, and file modification times are not staining-batch evidence.

## Online findings and author confirmation

The [online Lab Work review](ONLINE_LAB_WORK.md) records standard specimens 32469/32470 at `Sheet1!A84:K85` and target D. fuscus collection records at `Sheet1 (1)!A4:P5` and `A10:P11`. These are specimen inventories, not experimental staining-event fields. The [dataset link review](ONLINE_DATASET_LINKS.md) recovered genomic and phylogenetic resource links, without a microscopy batch mapping.

[AUTHOR_STAINING_CONFIRMATION.json](AUTHOR_STAINING_CONFIRMATION.json) preserves Jake's exact question and answer. The resulting run map keeps the original identifiers and roles. Existing local glass records join 22 of its 53 rows; missing glass identifiers stay blank. Scan dates remain separate from staining dates.

The manuscript and protocol locations are listed in [ONLINE_RECORD_RECOVERY.md](ONLINE_RECORD_RECOVERY.md). Protocol descriptions do not establish that every listed step or setting was used in the retained measurements.

## Remaining calibration work

The current computational denominator uses D. fuscus target Processes 314, 315, 320, and 321, rather than designated standards 413 and 414. Standard-slide QC, acquisition comparability, historical identity/preservation notes, and the absolute DNA-content anchor remain unresolved. The single-run answer closes the batch-membership question; publication-ready absolute genome-size estimates still require these checks.

Original spreadsheets, mappings, Methods, analyses, the frozen submission preparation, and selected figure artifacts were unchanged during this search. New evidence is retained separately in this directory. Independent verification is recorded in [ONLINE_RECOVERY_VERIFICATION.md](ONLINE_RECOVERY_VERIFICATION.md).
