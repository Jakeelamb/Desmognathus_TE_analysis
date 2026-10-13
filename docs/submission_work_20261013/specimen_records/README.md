# Methods specimen records

The supplied Google Sheets answer is preserved verbatim in `AUTHOR_SUPPLIED_RECORDS.txt`; `SOURCE.json` records its origin and checksum. Direct Sheets metadata access returned HTTP 403. This audit uses the supplied paste.

`MASTER_IMAGE_METADATA.csv` is an exact snapshot from the sibling microscopy project. `METADATA_SOURCE.json` records its path and checksum. The audit uses its image/specimen identities and historical preservation labels. Its genomic-accession fields are retained but do not enter the joins.

The verified outputs are in `audit_final/`:

- `mvz_source_records.csv` contains all 97 MVZ records.
- `rap_source_records.csv` contains all 47 RAP records.
- `assay_specimen_associations.csv` retains each assay/animal association, source row, image IDs, original and suggested taxa, dates, permits and comments.
- `REPORT.md` states the match coverage and unresolved identities.
- `receipt.json` pins every input and the script.

Run the registered Pixi task after integration, using a new destination:

```bash
pixi run --clean-env --locked --no-install specimen-records --output-dir .git/specimen-records-review
pixi run --clean-env --locked --no-install python -m pytest tests/test_specimen_records.py -q
```

The CLI also accepts `--repo-root` and `--source-dir`. It refuses an existing output directory, changed source checksums and writes inside retained scientific artifact directories. It needs only Python's standard library.

The target-panel specimen IDs have no MVZ/RAP prefixes in the released objects or local metadata. The audit therefore reports 41 of 42 morphometry animals and 49 of 50 IOD animals as unique numeric candidate associations, pending namespace confirmation. The author-designated standards have confirmed MVZ associations. The supplied MVZ `Date` field remains separate from collection information in `Verbatim Date`.

Analysis specimen 3463 is absent from the supplied tables. It contributes 50 morphometry objects and 31 IOD nuclei from Process_331. RAP3462 contributes five IOD nuclei from Process_324 and remains a separate animal. RAP3518 retains the source comment `Contaminated???` and contributes seven IOD nuclei from Process_355. The audit preserves the released measurements and makes no exclusion decision.

The author reports that 3463 likely arose from a transcription error during analysis entry. The intended source record has not been confirmed, so this clarification is recorded separately from the machine-generated associations. No identifier, animal count or measurement has been changed.

The source CSVs preserve MVZ suggested reidentifications without changing analysis taxonomy. In the used panel, the supplied original names for three D. orestes specimens are D. ochrophaeus, with D. orestes recorded separately as the suggested species.
