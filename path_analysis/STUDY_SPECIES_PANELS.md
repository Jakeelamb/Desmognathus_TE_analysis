# Study species panels

The project has three deliberately separate species panels. A species label is
not permission to substitute specimens, read accessions, assemblies, or image
measurements across sources.

| Panel | N | Purpose | Evidence boundary |
|---|---:|---|---|
| `te_resource_panel34_v1` | 34 | TE composition, RepeatMasker landscapes, diversity/PCA, and TE-only descriptive work | The active vetted WGS/SRA lookup only. It is not filtered by microscopy availability. |
| `cell_linked_panel21_v1` | 21 | Linked cell/nucleus morphology and relative nuclear-IOD descriptive work | Current frozen linked-cell and IOD summaries. It is not filtered by availability of a trusted WGS resource. |
| `integrated_path_panel18_v1` | 18 | Comparative/path models joining TE and cell evidence | Exact intersection of the prior two panels. No trait is imputed to enlarge it. |

## Versioned inputs and outputs

- Source declaration: `data/templates/analysis_species_panels.csv`
- Derived panel tables: `data/derived/panels/study_*_v1.csv`
- Immutable panel manifest: `data/derived/panels/study_species_panels_v1.manifest.json`
- TE34 descriptive cache: `../results/data/corrected/te34/`
- Cell21 descriptive cache: `../results/data/corrected/cell21/`

The builders are downstream only:

```bash
scripts/run_in_dusky.sh python scripts/processing/build_study_species_panels.py
scripts/run_in_dusky.sh python scripts/processing/build_te34_descriptive_bundle.py
scripts/run_in_dusky.sh python scripts/processing/build_cell21_descriptive_bundle.py
```

They do not run RepeatModeler, RepeatMasker, dnaPipeTE, mapping, segmentation,
or model training.

## Taxonomic/resource decisions

- `planiceps` is excluded: it is a mislabelled *D. fuscus* resource, not an
  independent tip.
- The active *D. fuscus* accession pair is retained as the expert-reidentified
  resource in the lookup table.
- Genomic resources labelled `folkertsi`, `imitator`, and `ochrophaeus` are not
  used as active TE resources. The cell data for *folkertsi* and *ochrophaeus*
  remain valid Cell21 observations; they are simply not in the integrated
  intersection.
- `aeneus`, `organi`, and `wrighti` are active TE34 resources. Their legacy
  nucleus-only image artifacts are not treated as current linked microscopy,
  and therefore do not enlarge Cell21 or path18.
- dnaPipeTE retry labels (for example `SRX19952890R2`) are not additional
  biological samples. Where multiple runs exist for one species, each run is
  first normalized independently and then equal-weight averaged into one
  species estimate, preventing double counting or depth-weighting a retry.

## Interpretation boundaries

TE34 figures are descriptive genomic-resource analyses. Cell21 figures are
descriptive linked-cell analyses. Only path18 figures may combine both layers,
and their relative nuclear-IOD variable remains an uncalibrated image phenotype
rather than absolute genome size or a C-value.
