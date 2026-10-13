# Specimen record checks

Checks ran in the isolated specimens worktree through the coordinator's locked default Pixi manifest. No new environment or dependency was installed.

```bash
pixi run --manifest-path /home/jake/Projects/Desmognathus_TE/pixi.toml --locked --no-install ruff check scripts/publication/reconcile_specimen_records.py tests/test_specimen_records.py
pixi run --manifest-path /home/jake/Projects/Desmognathus_TE/pixi.toml --locked --no-install python -m pytest tests/test_specimen_records.py -q
pixi run --manifest-path /home/jake/Projects/Desmognathus_TE/pixi.toml --locked --no-install python scripts/publication/reconcile_specimen_records.py --output-dir docs/submission_work_20261013/specimen_records/audit_final
```

Final results: lint passed; 10 tests passed; real audit exited 0. The tests run the CLI on both released inventories, compare two independent audit directories byte for byte, and verify refusal of an existing output directory, altered metadata, altered author source and destinations under protected scientific directories. Namespace collisions, duplicate rows, the missing 3463 identity, original/suggested taxa, ambiguous MVZ dates and the RAP contamination comment have explicit checks.

Initial lint failed because the new Python file had an unused executable shebang and two test subprocess calls omitted `check=False`. Removing the shebang and declaring the intended subprocess error handling fixed those failures. Final source was formatted and checked again.

The receipt records 1,152 morphometry objects from 42 animals and 805 IOD nuclei from 50 animals. All 94 assay/animal groups, including the two standards, have an exact filename/specimen match in the retained microscopy metadata. There are 41/42 and 49/50 unique numeric source candidates; 3463 is unmatched in both panels. Two standards have author-supported MVZ namespaces. Final outputs do not validate target namespaces, specimen preservation histories or absolute genome-size calibration.
