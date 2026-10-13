# Retained standard-source inventory (2026-10-13)

Verdict: **comparison remains unqualified**. Existing standard tiles, masks and candidate measurements are available for bounded QC. They are not the finalized target QC release. No scientific producer was rerun, original VSI/ETS bytes read, cloud accessed, or Methods/figure edited.

Repository HEAD `15dfb3c18755cc8659b626fa213087c62c6db02f`; upstream cellprofiler_test HEAD `d7c37fdaf2869245964b6710ca76696306c2daba`. Exact source hashes, mtimes, CSV row numbers, full source records and representative paths are in `source_inventory.json`.

## Reconstructing the approximately 41% claim

Source: `/home/jake/Projects/cellprofiler_test/output/runs/life_history_pmount_reference_v1/nucleus_measurements/all_measurements.csv`, 47 existing YOLO candidate nuclei, CSV lines 2–48. Process 413/specimen 32469 has **31 nuclei, median IOD 594.0784**; Process 414/specimen 32470 has **16 nuclei, median IOD 837.701417**. Difference = 243.623017000; `100*(837.701417/594.0784 - 1)` = **41.008563348%**, with 413's median as denominator. Symmetric percent difference is 34.030793577%; difference relative to the larger median is 29.082321225%. Thus “41%” is directional and requires its denominator.

IOD is integrated optical density summed over nuclear-mask pixels (arbitrary OD-pixel units), not pg or mean OD. The recorded pixel scale is 0.12 µm/pixel (area 0.0144 µm²/pixel). Both candidate tables store background `i_bg=255`; that is a processing value, not proof of original acquisition bit depth or linearity.

A distinct `/home/jake/Projects/cellprofiler_test/output/runs/life_history_pmount_reference_v1/imagej_nucleus_measurements/all_measurements.csv` contains 310/191 candidates and medians **220.6166/202.7581**. Its `imagej_linkage/image_summary.csv` strict-pair medians are **328.2552/144.97835** (10/8 strict pairs), and both `analysis_ready_image=False`. These are separate segmentations/selections; pooling them or substituting them for the 47-nucleus table is unsupported.

## Current denominator

`analyses/02_morphology/data/nuclear_iod_by_image.csv` target D. fuscus rows:

| Process | Specimen | Retained n | Image median IOD |
|---|---|---:|---:|
| Process_314_raw_green.ome.tiff | 3372 | 12 | 1079.68817100 |
| Process_315_raw_green.ome.tiff | 3373 | 8 | 1061.30965400 |
| Process_320_raw_green.ome.tiff | 3401 | 5 | 833.27677700 |
| Process_321_raw_green.ome.tiff | 3402 | 10 | 801.28195850 |

Equal-image mean of those four medians = **943.889140125 IOD**, from 35 reviewed target nuclei. Descriptive ratios of candidate standard medians to that denominator are 0.629394253 (413) and 0.887499794 (414); their equal-image candidate median mean 715.889908500 is 0.758447023 times the target denominator. These are arithmetic comparisons, not validated recalibration or new genome-size estimates.

## Acquisition, preservation and alternate images

`master_image_metadata.csv:123–124` records 413/414 as `image_type=pmount`, `mounting=dry`, `preservation=Fixed`, pixel size 0.12. Targets 314/315/320/321 (lines 4/5/10/11) are `brightfield`, `dry`, `Dried Blood`. The parent relays the author correction: all retained target measurements were 40× brightfield; oil 40× was excluded; the online 60× statement is wrong. **pmount does not establish oil**, and no original objective/exposure/gain/linearity fields are present in these rows.

Metadata lines 178–180 link alternate 468/469 to specimen 32469 and 470 to 32470; `mounting=pmount`, preservation Fixed. `slide_glass_mapping.csv:138–140` links 468/469 to glass 1417 and 470 to 1418 on directory date 20250523, versus designated standards glass1321/1322 on 20250521 (lines83/84). Those dates are scan directory labels, not staining dates. Both legacy pmount image indexes record 468/469/470 as `skipped_converged`; within the targeted all_measurements/image_summary/image_index/source_inventory search no alternate-image object table was found. This search is not an exhaustive archive search.

## Retained files and unavailable original sources

`/home/jake/Projects/cellprofiler_test/output/runs/life_history_pmount_reference_v1/staging/source_inventory.csv:2–3` retains exact original VSI hashes, companion ETS paths/byte counts, full image dimensions and extraction mode `temporary_full_green_ome`. Original paths under `/run/media/jake/easystore` are unavailable now (mount root absent). Legacy full OME-TIFF paths under `/home/jake/Projects/cellprofiler_test/data/pmount` are absent. Recorded hashes were preserved; original bytes were not rehashed. Existing reference staging has 24 retained tiles; all referenced candidate source tiles/masks exist. Headers on deterministic first-object tiles show 4096×4096, **8 bits/sample**. This only establishes retained tile depth; original acquisition depth remains unknown. Full raw image recovery is a separate source-availability step.

The first candidate per image is pinned with its exact tile/mask/object key in JSON, chosen by CSV order rather than IOD. Table fields are `mask_path` and `source_image_path`, not `nucleus_mask_file` or `image_file`; labels are tile-local.

## Author confirmation and strongest next action

`outputs/calibration_search_20261013/AUTHOR_STAINING_CONFIRMATION.json` (2026-10-13) records “One staining run with both standards”; staining date and original run identifier remain null. `ONLINE_SOURCE_PINS.json` pins read-only snapshots inspected 2026-10-13; snapshots do not recover native revision history. This inventory preserves their hashes and does not reinterpret the erroneous online magnification statement as acquisition evidence.

Next: inspect retained standard tiles and masks systematically, then determine whether the existing candidates can pass the **same frozen review criteria as targets**. Verify objective/acquisition comparability and recover original source availability before an approved remeasurement. Resolve segmentation/population disagreement between YOLO and ImageJ, acquisition response, and preservation/identity history before choosing a standard denominator. Single-run membership resolves batch linkage only.
