# Canonical Cell and Nucleus Workflow

The active production lineage is now limited to three trained artifacts:

1. Cell masks: custom Cellpose checkpoint `output/tile_training_round_v1/train/models/desmognathus_tile_round1`.
2. Nucleus masks and optical-density measurements: final YOLO checkpoint `runs/segment/output/yolo_nucleus_training_final_area500_shape/weights/best.pt`.
3. Linked-pair quality ranking: July 8 logistic-regression model `output/runs/mixed_cellpose_yolo_full_dataset_v1_bgclean/top50_linked_pair_review/pair_quality_model/pair_quality_model.joblib`.

The machine-readable contract, hashes, training-set sizes, inference settings, and filtering thresholds are in `docs/production_cell_nucleus_lineage_v1.json`. The only supported entrypoint for applying the stack to new VSI data is `scripts/run_production_cell_nucleus_pipeline.py`.

## Current Three-Species Application

The exact stack has been applied to all 22 staged slides and 230 tiles in `output/runs/life_history_target_vsi_v1`. Its YOLO-only top-candidate output is `output/runs/life_history_target_vsi_v1/top50_yolo_exact_v1`.

| Species | Eligible | Selected | Slides | Specimens | Median model probability | Shortfall |
|---|---:|---:|---:|---:|---:|---:|
| *D. aeneus* | 70 | 50 | 3 | 2 | 0.520 | 0 |
| *D. orestes* | 414 | 50 | 2 | 2 | 0.973 | 0 |
| *D. wrighti* | 22 | 22 | 4 | 3 | 0.129 | 28 |

The wrighti shortfall is an evidence limitation. The pipeline must not fill it with ImageJ-derived candidates or lower-quality nonphysical pairings under the label of the canonical method. More valid wrighti imagery, or an explicitly documented sensitivity analysis using a changed inclusion rule, is required to reach 50.

These are model recommendations, not a final freeze. Aeneus has 46/50 selected cells from one specimen and orestes has 48/50 from one specimen. That concentration is compatible with deliberately selecting the largest cells, but specimen-balanced and leave-one-specimen-out summaries are required before treating the species estimates as insensitive to specimen composition. Wrighti's lower classifier probability makes visual review especially important.

## Provenance Boundary

The July pair-quality classifier was trained on reviewed linked-pair features from the 21-species workflow. It is valid to reuse that classifier without a species feature on new linked pairs. We do **not** claim that the custom Cellpose checkpoint generated every older archived 21-species cell mask; the production claim here applies to the new three-species inference run and to the July pair-quality selection model.

ImageJ segmentation outputs and all earlier masks, images, crops, annotations, measurements, and review decisions remain preserved as historical data. They are excluded from the production candidate path, not deleted.

## Reproduce or Verify

Print the exact command plan and verify all three checkpoint hashes:

```bash
uv run python scripts/run_production_cell_nucleus_pipeline.py
```

Execute the entire canonical path:

```bash
uv run python scripts/run_production_cell_nucleus_pipeline.py --execute
```

Resume only the linkage/selection portion from existing segmentation outputs:

```bash
uv run python scripts/run_production_cell_nucleus_pipeline.py \
  --from-step linkage \
  --execute
```

The wrapper deliberately contains no ImageJ segmentation step.
