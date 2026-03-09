# Assembly Bias Audit

This note asks whether the core TE and paired-LTR history variables are
dominated by assembly quality / fragmentation metrics rather than biology.

## Source Basis

- `results/data/ltr_age/genome_assembly_manifest.csv`
- `path_analysis/data/derived/te_path_features.csv`
- `path_analysis/data/derived/ltr_history_features.csv`

## Coverage

- Assemblies with local verified FASTA files audited: `34`
- `ltr_line_logratio` non-missing species: `34`
- `order_pielou` non-missing species: `34`
- `weighted_te_divergence_p90` non-missing species: `34`
- `weighted_te_deletions_p90` non-missing species: `34`
- `ectopic_log10_mean_ratio` non-missing species: `30`
- `ltr_history_median_k2p_distance` non-missing species: `30`
- `ltr_history_n_pairs_estimated` non-missing species: `30`

## Bottom Line

- Most core TE and paired-LTR history variables do not show a strong assembly-bias signal in this screen.
- The main exception is `ltr_line_logratio`, which is negatively associated with `assembly_top10_bp_fraction` (rho `-0.590`, BH p `0.01531`).
- Treat the LTR:LINE balance feature as technically sensitive enough to require extra caution and cross-method validation, but do not downgrade the full TE feature set wholesale.

## Strongest Screened Associations

- `ltr_line_logratio` vs `assembly_top10_bp_fraction`: rho `-0.590`, raw p `0.0002431`, BH p `0.01531`, n `34`
- `ltr_line_logratio` vs `assembly_n_contigs_ge_100kb`: rho `-0.451`, raw p `0.007492`, BH p `0.1702`, n `34`
- `ltr_line_logratio` vs `assembly_longest_contig_bp`: rho `-0.447`, raw p `0.008104`, BH p `0.1702`, n `34`
- `order_pielou` vs `assembly_top10_bp_fraction`: rho `0.410`, raw p `0.01614`, BH p `0.2542`, n `34`
- `ltr_line_logratio` vs `assembly_total_bp`: rho `-0.359`, raw p `0.03694`, BH p `0.3564`, n `34`
- `ltr_line_logratio` vs `assembly_local_bytes_gz`: rho `-0.359`, raw p `0.03694`, BH p `0.3564`, n `34`
- `ltr_line_logratio` vs `assembly_n_contigs`: rho `-0.355`, raw p `0.0396`, BH p `0.3564`, n `34`
- `ltr_line_logratio` vs `assembly_n_base_fraction`: rho `-0.330`, raw p `0.05627`, BH p `0.411`, n `34`
- `ltr_line_logratio` vs `assembly_n90_bp`: rho `-0.327`, raw p `0.05872`, BH p `0.411`, n `34`
- `order_pielou` vs `assembly_longest_contig_bp`: rho `0.285`, raw p `0.1024`, BH p `0.61`, n `34`

## Files

- `results/data/assembly_bias/assembly_quality_metrics.csv`
- `results/data/assembly_bias/assembly_bias_correlations.csv`
- `results/data/assembly_bias/assembly_bias_models.csv`
