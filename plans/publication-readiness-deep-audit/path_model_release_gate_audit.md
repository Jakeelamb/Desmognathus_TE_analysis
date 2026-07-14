# Phylogenetic path-model release-gate audit

This audit is non-destructive: it evaluates stored ranking tables and does not rewrite historical model, edge, or figure outputs.

The gate follows the phylopath interpretation: a significant Fisher-C p-value rejects a causal model; CICc and weights must be finite; and supported models within delta CICc <= 2 form a competitive set rather than a unique winner.

Primary method source: [van der Bijl 2018](https://pmc.ncbi.nlm.nih.gov/articles/PMC5923215/).

## Result

- Families audited: 19.
- Families passing the narrow unique-top ranking gate: 10.
- Competitive supported sets requiring model-set or averaging language: 0.
- Blocked families: 9.
- Publication winner claims currently allowed: 0; every stored ranking uses historical/pre-audit inputs.
- `supported_unique_top_model`: 10.
- `blocked_nonfinite_ranking`: 5.
- `blocked_top_model_rejected`: 2.
- `blocked_no_globally_supported_model`: 2.

## Family-level decisions

| Family | n | Stored top | Fisher-C p | Status | Unique-top ranking gate |
|---|---:|---|---:|---|---|
| `genome_morphology` | 21 | `mediated_cell_size` | 0.860189 | `supported_unique_top_model` | yes |
| `te_genome_ectopic_organismal_primary_mediumplus` | 16 | `ectopic_baseline` | 0.0495809 | `blocked_top_model_rejected` | no |
| `te_genome_ectopic_organismal_primary_phylofill` | 16 | `ectopic_baseline` | 0.0495809 | `blocked_top_model_rejected` | no |
| `te_genome_ectopic_organismal_primary_strict_body` | 8 | `ectopic_baseline` | 0.2453 | `blocked_nonfinite_ranking` | no |
| `te_genome_ectopic_primary_mediumplus` | 16 | `ectopic_only` | 0.0022365 | `blocked_no_globally_supported_model` | no |
| `te_genome_ectopic_primary_phylofill` | 16 | `ectopic_only` | 0.0022365 | `blocked_no_globally_supported_model` | no |
| `te_genome_ectopic_primary_strict_body` | 8 | `ectopic_only` | 0.360235 | `blocked_nonfinite_ranking` | no |
| `te_genome_ltr_history_primary_mediumplus` | 15 | `te_baseline` | 0.425512 | `supported_unique_top_model` | yes |
| `te_genome_ltr_history_primary_strict_body` | 8 | `history_direct` | 0.155326 | `blocked_nonfinite_ranking` | no |
| `te_genome_morphology` | 18 | `te_evenness_path` | 0.695889 | `supported_unique_top_model` | yes |
| `te_genome_morphology_primary_mediumplus` | 18 | `te_evenness_path` | 0.700044 | `supported_unique_top_model` | yes |
| `te_genome_morphology_primary_phylofill` | 18 | `te_evenness_path` | 0.700044 | `supported_unique_top_model` | yes |
| `te_genome_morphology_primary_strict_body` | 9 | `te_direct_to_genome_nucleus` | 0.673949 | `blocked_nonfinite_ranking` | no |
| `te_genome_organismal_primary_mediumplus` | 18 | `te_baseline` | 0.633226 | `supported_unique_top_model` | yes |
| `te_genome_organismal_primary_phylofill` | 18 | `te_baseline` | 0.633226 | `supported_unique_top_model` | yes |
| `te_genome_organismal_primary_strict_body` | 9 | `te_baseline` | 0.286185 | `blocked_nonfinite_ranking` | no |
| `te_genome_primary_mediumplus` | 18 | `mediated_evenness` | 0.0896064 | `supported_unique_top_model` | yes |
| `te_genome_primary_phylofill` | 18 | `mediated_evenness` | 0.0896064 | `supported_unique_top_model` | yes |
| `te_genome_primary_strict_body` | 9 | `genome_null` | 0.167396 | `supported_unique_top_model` | yes |

The final column means only that the stored table passes the Fisher-C/finite-CICc/delta-CICc ranking gate. It does not authorize a publication claim from historical inputs. A lower-ranked supported model is reported only as a diagnostic when the stored top model is rejected; this audit does not automatically promote it. Candidate-family definitions and corrected inputs must be frozen before rerunning phylopath.
