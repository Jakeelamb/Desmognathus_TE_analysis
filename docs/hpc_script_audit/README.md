# Historical HPC source evidence

See [the audit](../HPC_SCRIPT_AUDIT.md) for findings and remaining gaps.
[manifest.json](manifest.json) records the exact source identities and SHA-256
checksums of 12 historical jobs and one command notebook. Five jobs reuse
existing copies under `outputs/submission_preparation_20261013/recovered_commands`;
the other seven and the notebook are copied here without changing their bytes.
[scope.json](scope.json) records the bounded local search.
[verification.json](verification.json) records the independent source and
preservation checks and the project validation result.

These files are historical evidence. They have not been executed or validated
as an upstream reproduction pipeline and are not registered as Pixi tasks.
The newly recovered job copies have non-executable permissions. They depend on
original HPC paths, libraries, environments and containers that have not all
been recovered.

Preserve the source bytes when correcting or reconstructing a workflow later.
The notebook contains transcription errors and parameter-scan examples.
`combined_pipeline.sh` references Necturus, and `Calculate_depth.sh` contains
calculation defects described in the audit. Neither should replace the retained
Desmognathus analysis or be run over its data.
