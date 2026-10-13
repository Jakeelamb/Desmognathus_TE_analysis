# Figure 4 source joins and validation

All inputs are snapshotted under this isolated run. No model was refitted.
S14 primary_ltr_element supplies all 380 raw ratios. S15 primary_iqr_filtered_coverage_ge_80pct supplies the 30 released medians, counts and bootstrap intervals.
S04 tree_tip joins to LTR species. Display values are pg/1C; corresponding Gbp values remain in the exported table for exact comparison to the stored regression inputs. Missing values are displayed as n/a without imputation.
Row order is the shared canonical source-tree tip-y order, matching the existing Figure 4 distribution.
Stored inputs/ltr_pgls species_data.csv is joined by species against S15 and S04; counts, medians, both interval endpoints and genome sizes must match within 1e-12.
Stored model_summary.csv and fitted_relationship.csv are copied byte-for-byte to plotted_data. They supply all PGLS inference and predictions.
Plot-only genome coordinates, interval endpoints, and label positions in Gbp are divided by 0.978 to display pg/1C. pgls_predictions_display_pg.csv and pgls_points_and_label_layout.csv add explicit pg columns. Ratio predictions, ratio intervals, model coefficients, and original CSV copies remain unchanged. A constant predictor-unit conversion changes an intercept but leaves the log-log slope and its inference unchanged; plotting the transformed x coordinates preserves the stored fitted line exactly.
Ratio density is estimated in log2 space for display, matching the prior figure. The PGLS axes are log10, matching the stored model. Both axes label the original ratios directly.
The collage is 180 x 280 mm; isolated panels are 180 x 180 and 180 x 120 mm. All text is Arial and at least 7.5 pt except the shared panel-tag settings, which are larger.
