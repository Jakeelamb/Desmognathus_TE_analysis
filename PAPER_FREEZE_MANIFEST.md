# Paper Freeze Manifest

This file records the current paper-facing freeze state for the Desmognathus TE project.

## Frozen Computational State

- Source code commit used for the frozen analysis state: `8c98c12f83d072561ed212e7f06de3f6563baf28`
- Manifest generated at: `2026-03-08 21:58:51 UTC`
- Tracked inventory files: `paper_freeze/key_file_manifest.csv`, `paper_freeze/results_inventory.csv`, and `paper_freeze/results_summary.csv`
- Ignored generated outputs inventoried locally: `555` files totaling `7.0 GB`

## Canonical Entry Points

- `scripts/processing/diversity_stats.py`
- `scripts/processing/pca.R`
- `scripts/processing/phylogenetic_pca_analysis.R`
- `scripts/processing/parse_repeatmasker_landscape.py`
- `scripts/processing/pgls_analysis.R`
- `scripts/processing/permanova_analysis.R`
- `scripts/processing/trait_evolution.R`
- `scripts/processing/analyze_phylogenetic_signal.R`
- `scripts/processing/analyze_phylogenetic_correlogram.R`
- `scripts/visualization/plot_all_te_landscapes.R`
- `scripts/R/visualization/te_phylo_landscape.R`
- `path_analysis/scripts/build_analysis_panels.py`
- `path_analysis/scripts/path_model_scaffold.R`

## Analysis Tiers

### Primary Paper Results

These are the most defensible current comparative outputs and should anchor the main text.

- `te_genome_primary_mediumplus`: `27` species. TE plus genome observed-only medium-plus.
- `te_genome_primary_strict_body`: `16` species. TE plus genome adult-oriented strict-body.
- `te_genome_ectopic_primary_mediumplus`: `24` species. TE plus genome plus ectopic observed-only medium-plus.
- `te_genome_ectopic_primary_strict_body`: `15` species. TE plus genome plus ectopic adult-oriented strict-body.
- `te_genome_organismal_primary_mediumplus`: `27` species. TE plus genome plus organismal observed-only medium-plus.
- `te_genome_organismal_primary_strict_body`: `16` species. TE plus genome plus organismal adult-oriented strict-body.
- `te_genome_ectopic_organismal_primary_mediumplus`: `24` species. TE plus genome plus ectopic plus organismal observed-only medium-plus.
- `te_genome_ectopic_organismal_primary_strict_body`: `15` species. TE plus genome plus ectopic plus organismal adult-oriented strict-body.

Primary winner snapshots:
- `path_analysis/results/te_genome_primary_mediumplus_model_ranking.csv` -> `mediated_evenness` (weight `0.999999617311239`)
- `path_analysis/results/te_genome_primary_strict_body_model_ranking.csv` -> `mediated_evenness` (weight `0.9713727673463166`)
- `path_analysis/results/te_genome_ectopic_primary_mediumplus_model_ranking.csv` -> `ectopic_only` (weight `0.8883015230056893`)
- `path_analysis/results/te_genome_ectopic_primary_strict_body_model_ranking.csv` -> `ectopic_only` (weight `0.8724612835313356`)
- `path_analysis/results/te_genome_organismal_primary_mediumplus_model_ranking.csv` -> `body_size_additive` (weight `0.44112883954290916`)
- `path_analysis/results/te_genome_organismal_primary_strict_body_model_ranking.csv` -> `body_size_additive` (weight `0.813699168971215`)
- `path_analysis/results/te_genome_ectopic_organismal_primary_mediumplus_model_ranking.csv` -> `te_body_size_baseline` (weight `0.9306266761692223`)
- `path_analysis/results/te_genome_ectopic_organismal_primary_strict_body_model_ranking.csv` -> `te_body_size_baseline` (weight `0.7348469778841338`)

### Supplementary Results

- Compositional PCA and phylogenetic PCA under results/tables/pca and results/figures/pca.
- Phylogenetic signal, correlogram, and TE-landscape figures under results/tables/phylogenetic_signal and results/figures/phylo_signal.
- Landscape CSVs and landscape/phylo-landscape figures under results/landscapes, results/figures/landscape, and results/figures/phylo_landscape.

### Provisional or Sensitivity-Only Results

- Morphology-linked path families remain planning or sensitivity analyses only.
- Phylogenetic trait-imputation columns remain sensitivity-only and should not replace observed values.
- Order-level all-feature ordinations are supplementary only and should not carry the main biological interpretation.

## Key Frozen Files

These are the most important frozen inputs and downstream tables to reference in methods, supplements, and reproduction notes.

- `input_data/lookup_table.txt`  size=`1377` bytes  sha256=`f55c5c0561d374dc1dedb53c85e286a9f5a159c6b15bf29d0c617dde230a8485`
- `input_data/phylogeny/desmo900dated_test.tre`  size=`1322` bytes  sha256=`1cf561526405c794b70d5767182eefe84dfa1552f9f9c0be6b36c465b19301aa`
- `results/data/dnaPipeTE_merged_classifications.csv`  size=`1813859288` bytes  sha256=`c75a62fefe6eab9952f32ddf6071bda44c457bf4ba99587ff694affe0ba88169`
- `results/data/dnaPipeTE_order_breakdown.csv`  size=`6306` bytes  sha256=`47bc6e616a8d69e035ea05e7b07bf09759ea36ef419f0f22b6db3570f3c06a57`
- `results/data/dnaPipeTE_superfamily_breakdown.csv`  size=`13862` bytes  sha256=`6e31027d69e45ec9426644d3a8c245b6c490fdaee6406946a9f0485f8cdf58f4`
- `results/data/diversity_order_stats.csv`  size=`8158` bytes  sha256=`bbd53d9292263c7edd17504051d26de01b9e67f134d24648cced209dd6a711c9`
- `results/data/diversity_superfamily_stats.csv`  size=`15320` bytes  sha256=`b1e63672e87f5041b0386a6097a8ecd0e50abee8e78578cde16ebb163513fa3d`
- `results/data/ectopic_recombination_filtered_3000bp_5+domains_no_unknown_species.csv`  size=`206324` bytes  sha256=`5ab614026e6f63dac908301151186de0b1cf3def81b5a4513dcb64e6947ed22b`
- `results/tables/pca/te_pca_analysis_manifest.csv`  size=`1359` bytes  sha256=`64c551d9925b1dda8f8038622cd01181909e30030c08de6d9e7726613f4e6636`
- `results/tables/pca/te_ppca_analysis_manifest.csv`  size=`794` bytes  sha256=`936a131a112a5ba1381027c34ffaf1dbaf8241c73946c2619cf1195c439fe5c0`
- `path_analysis/data/derived/te_path_features.csv`  size=`48141` bytes  sha256=`f29ab1a988332e0d5d5ed0066c403d77dd700be6e76d56ba3539ec1f1864d08d`
- `path_analysis/data/derived/te_model_feature_panel.csv`  size=`31967` bytes  sha256=`19b6f0d56e4ef44b4d2cf34c46a146660bb55ba223d1c42ada2142a15beab071`
- `path_analysis/data/derived/path_input_master.csv`  size=`175597` bytes  sha256=`7945cd8e88745f9dbe52b86f41afbc116fca46f07b7457fc1e3676bcca98a7b2`
- `path_analysis/data/derived/analysis_panel_summary.csv`  size=`4951` bytes  sha256=`f16c47922686a86b497402469631d49193762825723ef1d1761bc9b1d70eccab`
- `path_analysis/data/derived/analysis_species_readiness.csv`  size=`35473` bytes  sha256=`41f79e7eb50f403d04039d2f5d80cf950e977e81799d59d21b37193109c41ff4`
- `path_analysis/data/derived/phylofill_panel_comparison.csv`  size=`700` bytes  sha256=`ba9d3096ab0dcb82e386602932f2b46ef8fc6e19b023051785dda9b536983770`

## Ignored Results Inventory Summary

- `results/analysis_report.html`: `1` files, `19.0 KB`
- `results/data`: `35` files, `6.9 GB`
- `results/data/divergence`: `1` files, `3.1 MB`
- `results/data/permanova`: `6` files, `44.6 KB`
- `results/data/pgls`: `3` files, `12.1 KB`
- `results/data/trait_evolution`: `1` files, `6.0 KB`
- `results/figures`: `38` files, `8.2 MB`
- `results/figures/TE_diversity_hierarchical_donuts`: `34` files, `10.9 MB`
- `results/figures/TE_diversity_multi_donuts`: `34` files, `11.2 MB`
- `results/figures/TE_diversity_superfamily_donuts`: `34` files, `6.9 MB`
- `results/figures/class_phylogeny`: `7` files, `1.6 MB`
- `results/figures/divergence`: `156` files, `36.8 MB`
- `results/figures/divergence_phylogeny`: `34` files, `10.6 MB`
- `results/figures/diversity`: `1` files, `0 B`
- `results/figures/diversity_phylogeny`: `6` files, `1.8 MB`
- `results/figures/diversity_threshold_plots`: `6` files, `2.7 MB`
- `results/figures/landscape`: `4` files, `1000.1 KB`
- `results/figures/order_phylogeny`: `8` files, `1.8 MB`
- `results/figures/pca`: `17` files, `2.4 MB`
- `results/figures/permanova`: `8` files, `1.2 MB`
- `results/figures/pgls`: `4` files, `340.4 KB`
- `results/figures/phylo_landscape`: `3` files, `754.7 KB`
- `results/figures/phylo_signal`: `3` files, `1.1 MB`
- `results/figures/phylogeny`: `3` files, `662.6 KB`
- `results/figures/superfamily_phylogeny`: `13` files, `2.8 MB`
- `results/figures/te_landscape`: `6` files, `1.3 MB`
- `results/figures/trait_evolution`: `4` files, `397.7 KB`
- `results/landscapes`: `35` files, `2.9 MB`
- `results/phylogeny`: `1` files, `1.1 KB`
- `results/tables/diversity`: `1` files, `0 B`
- `results/tables/pca`: `45` files, `237.9 KB`
- `results/tables/phylogenetic_signal`: `2` files, `14.1 KB`
- `results/tables/phylogeny`: `1` files, `845 B`

## Scope Notes

- `results/` remains git-ignored. The freeze state is therefore represented by the tracked checksum inventories rather than by versioning the generated files themselves.
- Observed trait values remain primary. `phylo_` columns remain explicit sensitivity-only inference.
- Morphology-linked path families remain excluded from the primary paper claims until independent final genome-size estimates are available.

