# Morphology: cells, nuclei, and IOD-derived genome size

This directory is the active, compact morphology component of the paper repo.
It starts from finalized manual-review decisions and does **not** rerun image
segmentation, model training, or the large review galleries.

## Canonical data

The study panel is `data/identity/path24_panel.csv`: 24 species with finalized
cell/nucleus morphology and conditional genome-size estimates. Nuclear IOD is
the measurement and calibration substrate for those estimates and remains fully
available for audit.

| File | Unit and final size | Role |
| --- | --- | --- |
| `data/cell_nucleus_objects.csv` | 1,152 reviewed cell/nucleus objects; 24 species | Canonical object-level morphology table |
| `data/cell_nucleus_species_estimates.csv` | 24 species | Canonical species-level cell/nucleus estimates and biological support |
| `data/nuclear_iod_objects.csv` | 805 reviewed nuclei; 24 species | Canonical object-level nuclear-IOD table |
| `data/relative_nuclear_iod_species.csv` | 24 species | Canonical species-level IOD summaries and conditional *D. fuscus*-anchored genome-size estimates |
| `data/nuclear_iod_by_image.csv` | 51 images; 24 species | Unnormalized image-level nuclear-IOD summary |
| `data/cell_mask_review_audit.csv` | 129 review-audit rows | Audit of the manually extended morphology release |
| `data/nuclear_iod_quality_balance.csv` | 24 species | Per-species IOD support and balance |
| `data/nuclear_iod_quality_diagnostics.csv` | 2 diagnostics | Release-level IOD quality checks |

Publication supplements S18-S25 are copies of these canonical tables. Do not
edit the publication copies as a second source of truth.

## What can be changed locally

Use `Publication/Desmognathus_study_data.Rmd` to change display filters,
summaries, or ggplot views. It reads only the compact tables above, so none of
this requires Cellpose, YOLO, raw microscopy images, or an HPC cluster.

From the repository root:

```bash
make setup-python
uv run python analyses/02_morphology/validate_release.py
make validate
make report
```

`validate_release.py` replays the compact release without raw imagery. It checks
the S19 species summaries against S18, the S23 image summaries and both the S22
IOD audit fields and conditional genome-size fields against S21, the fixed
physical-area conversion, the S24 balance rule, and the declared animal-support
totals.

If a reviewed canonical table is intentionally revised, update its matching
decision/provenance record, then refresh the publication copies and catalog:

```bash
make refresh
```

Treat the object-level review decisions as frozen evidence, not tuning knobs.
A new object inclusion or exclusion is a new reviewed release and must retain
the old decision, source hash, and reason.

## Frozen provenance and upstream boundary

`provenance/` must remain with the paper repo. It contains the pre-extension
frozen object tables, manual review decisions, review-finalization manifest,
IOD source manifest, legacy hash registry, and the production model-lineage
documents needed to audit the active CSVs. These files document how the release
was made; the report should not rewrite them.

`provenance/review_finalization.json` is the authority for the final reviewed
S18/S21 release. `provenance/production_cell_nucleus_lineage.json` is an earlier
candidate-generation and model-lineage snapshot; its `ready_for_final_freeze`
status describes that pre-final application stage and does not supersede the
later review-finalization manifest.

The following remain upstream in the sibling `../cellprofiler_test` repository
or in managed archival storage and are not active paper-repo inputs:

- raw OME-TIFF/VSI microscopy files and external-drive source inventories;
- tiles, masks, candidate objects, HTML/PNG review galleries, and intermediate
  CellProfiler/Cellpose/YOLO runs;
- Cellpose and YOLO weights, pair-quality models, training labels, and GPU
  environments;
- superseded 21-species/final18 tables, alternate sensitivity galleries, and
  non-production rebuild attempts.

The full pre-cleanup repository is preserved in the sibling
`../Desmognathus_TE_archive` directory. Do not copy historical tables back into
this directory merely because an archived script expects their old paths. Only
the compact finalized tables and provenance above define the paper release.

## Interpretation boundary

Conditional *D. fuscus*-anchored genome size in pg/1C is the default downstream
Path24 analysis trait. It is calculated from the target-to-*D. fuscus*
equal-image IOD ratio using the published pooled-assembly estimate of `16.1
Gbp`, converted with `1 pg = 0.978 Gbp` to `16.462167689 pg/1C`. This is not an
independent genome-size assay or an independently validated absolute C-value.

The author confirmed that Process_413/specimen 32469 and
Process_414/specimen 32470 were the *D. fuscus* standards and were stained with
all retained unknowns in one Feulgen run. The retained standard images were
40× brightfield without oil, using the green channel and the same fixation
protocol as the other slides. The released conditional estimate still uses
the finalized *D. fuscus* target-species IOD estimate as its denominator.
Standard-nucleus acceptance, the approximately 41% difference between the two
standard-slide median IOD values, and absolute calibration remain unresolved.
The [calibration decision](../../outputs/standard_calibration_review_20261013/CALIBRATION_DECISION.md)
records these limits. S21-S25 retain the IOD measurement, quality, and
calibration audit.

The stored physical-area conversion is exactly `0.0144 µm²/pixel`. Per-nucleus
integrated optical density (IOD) is the sum of optical density over nuclear-mask
pixels; it is distinct from mean optical density. The primary species estimand
is the equal-image mean of within-image median nuclear IOD. Dividing those 24
species estimates by their Path24 median gives the retained audit field
`relative_iod_index`, whose Path24 median is one. Dividing each estimate by the
*D. fuscus* estimate and multiplying by `16.462167689 pg/1C` gives the default
conditional genome-size trait and interval. S23 contains the unnormalized image
medians and means used in that aggregation; it is not itself a relative-IOD,
genome-size, or ratio table.

Genomic and microscopy records are joined at species level. They are not
measurements from the same specimen, and object-level image paths do not prove
genomic-sample identity.

## Review support and naming hazards

The review targets were 50 morphology objects per newly reviewed species and
40 IOD nuclei per newly reviewed species. Final retained counts below those
targets are intentional and must remain visible:

| Species | Morphology objects | IOD nuclei | Note |
| --- | ---: | ---: | --- |
| *D. aeneus* | 44 | 21 | Manual Path24 inclusion; below both targets |
| *D. orestes* | 50 | 22 | Manual Path24 inclusion; IOD below target |
| *D. wrighti* | 8 | 9 | Manual Path24 inclusion; below both targets |
| *D. ochrophaeus* | legacy reviewed panel | 32 | Existing Path24 species; IOD below target |

`specimen_id` denotes one animal. The morphology stream contains 42 animals,
the IOD stream contains 50, and their union contains 51. Morphology support is
one to three animals per species; *D. anicetus*, *D. apalachicolae*,
*D. auriculatus*, *D. bairdi*, *D. gvnigeusgwotli*, *D. intermedius*, and
*D. perlapsus* each have one animal. S19 therefore reports `n_size_specimens`,
`largest_specimen_n`, and `largest_specimen_fraction` alongside object and image
counts. Its bootstrap intervals remain conditional paired-object intervals, not
population-level among-animal uncertainty.

For nuclear IOD, *D. anicetus*, *D. bairdi*, and *D. gvnigeusgwotli* each have
one image from one animal. The S24 balance rule is
`max_abs_standardized_mean_difference <= 0.10` and `max_ks_distance <= 0.25`;
*D. aeneus*, *D. ochrophaeus*, *D. orestes*, and *D. wrighti* fail that rule.
These limitations are part of the released sensitivity evidence, not missing
values to fill or reasons to alter frozen review decisions.

Keep these identity rules explicit:

- Active morphology tables use labels such as `D. orestes`, while the identity
  panel uses bare names such as `orestes`. Normalize through
  `data/identity/species_taxonomy_crosswalk.csv`; do not join by ad hoc string
  replacement.
- The collaborator-SVL aliases and their max-of-maxima rule are defined once in
  the [shared identity-resolution contract](../../docs/IDENTITY_RESOLUTION.md).
- The public genomic resource SRX20497025 / GCA_032353935.1 is treated as
  *D. fuscus* by an accession-specific expert decision despite its public
  *D. planiceps* label. This is **not** a general synonymy and does not rename
  any morphology specimen or morphology row labeled *D. planiceps*.
- Historical manifests may contain old absolute paths or spelling mistakes.
  Preserve them as provenance, but use the current crosswalk and canonical
  component filenames for analysis.
