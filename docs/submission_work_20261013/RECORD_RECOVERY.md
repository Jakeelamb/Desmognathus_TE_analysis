# Methods record recovery — 13 October 2026

This memo records a bounded, read-only search for laboratory and analysis records. Its findings support the Methods audit. Data, images, masks, figures and historical outputs were preserved.

## Original acquisition sources

The original Olympus sources are addressed under `/run/media/jake/easystore`. That directory was absent during this check. No VSI or ETS files were found within the bounded sibling-project search. An unavailable external source is not evidence that it has been lost.

The sibling `cellprofiler_test/slide_glass_mapping.csv` already supplies scan dates and relative drive paths. More specific retained source identities are available in `cellprofiler_test/output/runs/life_history_pmount_reference_v1/staging/slides/Process_413_green.ome/summary.json` and the corresponding Process_414 record:

| Standard | Original VSI path beneath easystore | Recorded VSI bytes / SHA256 | Recorded companion ETS bytes | Height × width, pixels |
|---|---|---|---|---|
| Process_413 | `20250521/Slide glass_1321/Process_413.vsi` | 581,841 / `2fe307f8cb2fafdb58b9cfa8dcfb76489f3f6969db30222dc910870cdb950179` | 8,484,460,228 | 48,786 × 50,334 |
| Process_414 | `20250521/Slide glass_1322/Process_414.vsi` | 417,611 / `480f0182a9f41c065f1303c5ceaabd4d4ebc1d3865180c985a782eeeb0d99047` | 20,591,410,196 | 78,680 × 78,126 |

Each companion path is `_Process_ID_/stack1/frame_t_0.ets` within its slide-glass directory. Both summaries record background 255 and temporary full green OME extraction. Their filesystem modification times are 13 July 2026, 21:36:34 and 21:46:36 UTC; these are summary-file timestamps, not acquisition or staining dates. Source sizes and hashes above were read from retained summaries, not reverified against inaccessible originals.

The 22 three-species target staging summaries in `cellprofiler_test/output/runs/life_history_target_vsi_v1/staging/slides/*/summary.json` likewise retain original VSI paths, hashes and sizes, companion paths and sizes, dimensions, background and tile-selection records. They do not record camera model, exposure, gain, gamma, white balance, sensor bit depth or detector linearity. Their recorded source sizes total 11,439,478 VSI bytes and 220,168,910,464 ETS bytes; this is historical metadata, not a current storage inventory.

`cellprofiler_test/scripts/stage_vsi_green_tiles.py:187` calls Bio-Formats `showinf -nopix`, then extracts dimensions. Its summary writer at line 534 preserves the fields listed above. The filename search did not recover saved showinf reports, acquisition metadata files or OME XML sidecars from the sibling output tree. The raw acquisition settings may still be recoverable from the original VSI headers when external storage is available. Do not rerun image conversion or read whole mosaics merely to inspect headers. Local green OME-TIFF derivatives exist under `cellprofiler_test/data/brightfield`; header-only inspection may supply further clues, but copied OME metadata is explicitly not trusted by the staging reader at line 287.

## Segmentation lineage

The current sibling production contract is `cellprofiler_test/docs/production_cell_nucleus_lineage_v1.json`, lineage `desmognathus_cell_nucleus_production_v1_20260714`. It is byte-identical to [the active repository copy](../../analyses/02_morphology/provenance/production_cell_nucleus_lineage.json), SHA256 `754be4859f416750f383c589e4e6c80f34575db98d51e1012b8cc78607f121c6`. The sibling file's modification timestamp is 14 July 2026, 15:00:52 UTC. The sibling Git HEAD is `d7c37fdaf2869245964b6710ca76696306c2daba` (author timestamp 13 March 2026); the newer filesystem records must not be attributed to that older revision without Git tracking evidence.

All three referenced model files exist. Their actual bytes were hashed read-only through the project's locked Pixi Python runtime and matched the contract:

| Artifact beneath cellprofiler_test | Bytes | Verified SHA256 |
|---|---|---|
| `output/tile_training_round_v1/train/models/desmognathus_tile_round1` | 609,388,583 | `e04ee76afc0bd1e1526a3bb455aac51e5389e87442129eb45ac87653014a7c14` |
| `runs/segment/output/yolo_nucleus_training_final_area500_shape/weights/best.pt` | 6,617,757 | `36d20b16f4bef51c6e3657f655471e7dedb6bba149a996ded1f2161c676f7f19` |
| `output/runs/mixed_cellpose_yolo_full_dataset_v1_bgclean/top50_linked_pair_review/pair_quality_model/pair_quality_model.joblib` | 13,434 | `e47b06dfc2657899ed363fc67aecd9af894790eae345985be0630a02d460cda8` |

The contract includes training manifests, training counts, inference thresholds and a 22-slide, 230-tile three-species application. These are useful reproducibility records. Its explicit provenance boundary, also stated in `cellprofiler_test/docs/PRODUCTION_CELL_NUCLEUS_LINEAGE.md:27`, says that it does **not** establish which custom Cellpose checkpoint generated every older 21-species cell mask. The recommended candidate counts in that contract are not the final manually accepted morphometry release. Preserve the distinction; do not substitute the contract's candidate counts or checkpoint attribution into released Methods.

## Specimens and fixation

The local microscope crosswalk and master metadata distinguish `Dried Blood` and `Fixed`, but do not supply concentrations, fixation times or collection protocols. The sibling project `Dried_blood_all_samples` contains 327 JPEGs with RAP specimen labels in filenames; no VSI/ETS sources or permit/voucher documents were recovered there in this bounded search. Filenames do not establish lifetime preservation or handling.

New collection records supplied by the author are reconciled in the specimen association audit. They supersede the earlier collection-record search gap for matched specimens. Original source values remain preserved.

## Focal tree and HPC boundary

The existing [focal-tree input](../../analyses/03_phylogenetic_path/trees/source_time_tree_46.tre) is the recoverable computational input; direct provision by Alex Pyron is recorded in `data/identity/source_manifest.csv`. The opened archival `plans/publication-readiness-deep-audit/phylogeny_release_audit_analysis18_v1.md` repeats the unresolved publication/archive identifier, original filename/hash linkage, units, estimation method, calibration identities/priors and uncertainty sample. Its older analysis18 rounding and root-age statements are not current Path24 provenance. The published Stewart–Wiens sensitivity sources do not identify the focal Pyron tree. Ask the provider for the original tree and associated Methods/calibration record, then compare exact bytes before changing any description.

The [completed HPC audit](../HPC_SCRIPT_AUDIT.md) remains the source registry. Its scopes already include the archive, backups, research notes and results roots; this investigation did not repeat that wholesale. No additional exact custom library, container or submitted-job receipt was recovered from the newly examined microscopy metadata and project-note locations. The next upstream recovery target remains original HPC home/scratch or an exact backup. Source recovery is distinct from executing or rebuilding a producer.

## Bounds and next actions

Sources examined: active Methods/provenance records; sibling `cellprofiler_test` metadata, manifests, production-lineage documents, staging summaries and staging source code; bounded original-source/model filenames under the sibling project and `Dried_blood_all_samples`; Desmognathus Documents notebooks/workflow notes; and the focal-tree archival audit. No whole-home text search, unrelated conversations, image-pixel parsing, segmentation inference, HPC workflow execution or source editing was performed.

1. Make the original easystore available read-only, locate both standard VSI headers plus representative retained targets and their ETS companions, and extract header metadata into a separate audit record. Match hashes to retained summaries before interpreting settings.
2. Trace final-release object IDs to the actual inference run manifests for legacy masks. Existing trained artifacts are available and hash-verified; their availability does not close checkpoint-to-mask attribution.
3. Join the author-supplied specimen records to retained assays and retain any unresolved fixation, handling or identity details.
4. Obtain the focal tree's original file and source/calibration record from its provider. Preserve the present tree and all outputs.

The active repository starting revision for this audit was `d64f5a13981e1ff0edc2e30d66940110ac04be06` (13 October 2026, 14:08:21 −07:00). Only this new recovery memo is written by this investigation.
