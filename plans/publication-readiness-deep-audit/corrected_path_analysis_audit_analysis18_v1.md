# Corrected phylogenetic path-analysis audit — final 18 species (v1)

## Verdict

The corrected path implementation is technically reproducible and substantially stronger than the historical analysis: **958 fits completed with 0 failures** across 18 measurement specifications, the published main time tree, 200 published time-tree replicates, and leave-one-species-out analyses. Every candidate set now has an explicit null, every basis set is exported, rejected global models are fail-closed, and all outputs identify the image node as a **relative nuclear-IOD proxy**, never absolute genome size.

The path analysis is nevertheless **not approved for a publication claim about genome size or causal mechanisms**. The decisive reasons are upstream measurement validity, not software execution: relative IOD is uncalibrated and algebraically includes nuclear area; segmentation has no held-out final-panel species validation; the largest-50 trait is an upper-tail estimand; within-species measurement error is not propagated; and the candidate DAGs were formalized during this audit rather than genuinely preregistered before viewing the data. The terminal:internal family is additionally rejected by its global-fit test in every full-panel specification.

## Anchor result

Anchor specification: equal-image frozen selected-50 morphology, image-QC-passing relative IOD, and the published main Stewart–Wiens time tree.

| family | model | global_p | weight | gate |
| --- | --- | --- | --- | --- |
| Integrated TE–IOD–morphology | LTR → evenness → IOD → nucleus → cell | 0.789 | 0.998 | supported unique top model |
| Relative IOD → morphology | IOD → nucleus → cell | 0.772 | 0.974 | supported unique top model |
| TE composition → relative IOD | LTR → evenness → IOD | 0.118 | 0.998 | supported unique top model |
| Terminal:internal proxy → relative IOD | Null | 0.007 | 0.603 | blocked no globally supported model |

The selected integrated graph has these standardized direct estimates (approximate Wald intervals from `phylopath`):

| edge | estimate | interval_excludes_zero |
| --- | --- | --- |
| even → iod | -0.431 [-0.873, 0.011] | False |
| iod → ns | 0.552 [0.243, 0.860] | True |
| ltr → even | -0.654 [-0.877, -0.432] | True |
| ns → cs | 0.625 [0.242, 1.007] | True |

The `even → IOD` interval includes zero even though the mediated graph has overwhelming CICc weight. Model weight compares whole conditional-independence structures; it is not a substitute for uncertainty on a particular arrow. The positive `IOD → nucleus area` path is not an independent biological validation because the stored IOD is exactly nuclear pixel area multiplied by mean optical density.

## Robustness coverage

- Data specification: 84 family/specification/tree fits. The integrated TE-evenness chain ranked first in 36/36; the IOD–nucleus–cell chain ranked first in 36/36.
- Tree uncertainty: each family was refit on the published main tree plus 200 published time trees. The same top model was selected on every tree for all four families; the terminal:internal winner still failed global fit on every tree.
- Species influence: the same TE-IOD, morphology, and integrated winners remained top after omitting each species. Terminal:internal models failed global fit in all 16 analyzable omissions.
- Measurement uncertainty remains incompletely propagated. Switching estimands and IOD QC subsets is a sensitivity analysis, not an error-aware hierarchical SEM.

| phase | family_label | n_fits | modal_top_model_label | modal_top_model_rate | global_fit_pass_rate | minimum_top_weight | maximum_top_weight |
| --- | --- | --- | --- | --- | --- | --- | --- |
| data | Integrated TE–IOD–morphology | 36 | LTR → evenness → IOD → nucleus → cell | 1.000 | 1.000 | 0.997 | 1.000 |
| data | Relative IOD → morphology | 36 | IOD → nucleus → cell | 1.000 | 1.000 | 0.561 | 0.996 |
| data | TE composition → relative IOD | 6 | LTR → evenness → IOD | 1.000 | 1.000 | 0.998 | 1.000 |
| data | Terminal:internal proxy → relative IOD | 6 | Null | 1.000 | 0.000 | 0.603 | 0.634 |
| trees | Integrated TE–IOD–morphology | 201 | LTR → evenness → IOD → nucleus → cell | 1.000 | 1.000 | 0.998 | 0.999 |
| trees | Relative IOD → morphology | 201 | IOD → nucleus → cell | 1.000 | 1.000 | 0.973 | 0.974 |
| trees | TE composition → relative IOD | 201 | LTR → evenness → IOD | 1.000 | 1.000 | 0.998 | 0.999 |
| trees | Terminal:internal proxy → relative IOD | 201 | Null | 1.000 | 0.000 | 0.603 | 0.603 |
| loo | Integrated TE–IOD–morphology | 18 | LTR → evenness → IOD → nucleus → cell | 1.000 | 1.000 | 0.695 | 0.999 |
| loo | Relative IOD → morphology | 18 | IOD → nucleus → cell | 1.000 | 1.000 | 0.922 | 0.989 |
| loo | TE composition → relative IOD | 18 | LTR → evenness → IOD | 1.000 | 1.000 | 0.987 | 0.999 |
| loo | Terminal:internal proxy → relative IOD | 16 | Null | 0.938 | 0.000 | 0.443 | 0.833 |

## Actual-tree simulation calibration

Independent Brownian residual traits and the observed chain were simulated on the published final-18 tree. Null scenarios use 200 replicates per family; signal scenarios use 100 replicates per family at 0.5×, 1×, and 1.5× the observed standardized coefficients. “Recovery” requires the generating model to rank first and pass the unique-support gate. “False +” uses the matching strict definition: a non-null model ranks first and receives the unique-support gate under independent null traits. Wilson 95% intervals quantify Monte Carlo/binomial uncertainty. Coefficient bias and approximate interval coverage are exported for every generating edge.

| family | scenario | n_successful | modal_top_model | rate_95_ci | mean_top_weight |
| --- | --- | --- | --- | --- | --- |
| Integrated TE–IOD–morphology | independent null | 200 | integrated_null | 0.095 [0.062, 0.144] | 0.892 |
| Integrated TE–IOD–morphology | observed chain ×0.5 | 100 | integrated_null | 0.100 [0.055, 0.174] | 0.786 |
| Integrated TE–IOD–morphology | observed chain ×1 | 100 | te_evenness_chain | 0.500 [0.404, 0.596] | 0.807 |
| Integrated TE–IOD–morphology | observed chain ×1.5 | 100 | te_evenness_chain | 0.810 [0.722, 0.875] | 1.000 |
| Relative IOD → morphology | independent null | 200 | morphology_null | 0.200 [0.150, 0.261] | 0.575 |
| Relative IOD → morphology | observed chain ×0.5 | 100 | iod_nucleus_cell_chain | 0.190 [0.125, 0.278] | 0.625 |
| Relative IOD → morphology | observed chain ×1 | 100 | iod_nucleus_cell_chain | 0.810 [0.722, 0.875] | 0.894 |
| Relative IOD → morphology | observed chain ×1.5 | 100 | iod_nucleus_cell_chain | 0.850 [0.767, 0.907] | 1.000 |
| TE composition → relative IOD | independent null | 200 | proxy_null | 0.190 [0.142, 0.250] | 0.529 |
| TE composition → relative IOD | observed chain ×0.5 | 100 | proxy_null | 0.190 [0.125, 0.278] | 0.563 |
| TE composition → relative IOD | observed chain ×1 | 100 | mediated_evenness | 0.790 [0.700, 0.858] | 0.849 |
| TE composition → relative IOD | observed chain ×1.5 | 100 | mediated_evenness | 0.930 [0.863, 0.966] | 1.000 |

Simulation is a diagnostic under the stated Brownian data-generating process, not proof that the real DAG is causal. Poor recovery would block model discrimination; good recovery cannot rescue invalid trait meaning or unmodeled measurement error.

The strict unique false-selection rates were 19.0% for TE–IOD, 20.0% for IOD–morphology, and 9.5% for the integrated family. A non-null model ranked first under the null even more often (49.5%, 45.5%, and 12.0%, respectively), usually with model-set competition. At the observed coefficient scale, strict generating-model recovery was 79.0%, 81.0%, and only 50.0%. Thus the integrated pattern is empirically stable in the observed dataset but only moderately identifiable at `n=18`; the false-selection surface is too large for confirmatory language.

Across generating edges and signal scales, approximate 95% interval coverage ranged from 0.860 to 1.000; mean coefficient bias ranged from -0.040 to 0.055. See the calibration CSV for each edge and effect scale.

## Publication release matrix

| component | status | permitted_language |
| --- | --- | --- |
| LTR deletion / ectopic recombination | BLOCKED_AS_RATE | Exploratory terminal:internal deletion-footprint proxy only. |
| TE diversity indices | APPROVED_DESCRIPTIVE | Descriptive Shannon, Gini-Simpson, Hill, richness, and evenness patterns. |
| TE compositional PCA | APPROVED_DESCRIPTIVE | Descriptive ordination; no unrestricted PERMANOVA group claim. |
| Time-tree sensitivity | APPROVED_FOR_SENSITIVITY | Tree-set sensitivity, with source limitation stated. |
| Cell/nucleus linkage | APPROVED_PROVENANCE | Linked upper-tail morphometry provenance. |
| Segmentation generalization | BLOCKED | No final-panel segmentation-accuracy claim. |
| Morphology trait | SENSITIVITY_ONLY | Upper-tail cell/nucleus area sensitivity, not typical erythrocyte size. |
| Absolute genome size | BLOCKED | Relative nuclear-IOD proxy only; never pg or C-value. |
| Phylogenetic path analysis | EXPLORATORY_ONLY | Exploratory phylogenetically corrected association-model sensitivity; no causal claim. |

## Comparison with leading methods

- `phylopath` requires a biologically justified model set, basis-set/global-fit evidence, CICc comparison, and explicit model uncertainty (van der Bijl 2018, [doi:10.7717/peerj.4718](https://doi.org/10.7717/peerj.4718)). Those computational reporting surfaces are now present.
- The complete model comparison is propagated over all 200 time trees supplied for comparative analyses by Stewart and Wiens (2025, [doi:10.1016/j.ympev.2024.108272](https://doi.org/10.1016/j.ympev.2024.108272)), matching current tree-uncertainty expectations.
- Actual-tree Monte Carlo addresses the small-panel model-selection problem highlighted by Boettiger, Coop, and Ralph (2012, [doi:10.1111/j.1558-5646.2011.01574.x](https://doi.org/10.1111/j.1558-5646.2011.01574.x)), but a hierarchical error-aware model remains the stronger next analysis.
- The image assay does not meet Feulgen densitometry standards in Hardie, Gregory, and Hebert (2002, [doi:10.1177/002215540205000601](https://doi.org/10.1177/002215540205000601)) or the same-batch reference practice used in salamanders by Mueller et al. (2008, [doi:10.1016/j.zool.2007.07.010](https://doi.org/10.1016/j.zool.2007.07.010).)
- Recent salamander TE-diversity studies report clearly defined Shannon/Gini-Simpson metrics and phylogenetically controlled trait tests; matching their quality means transparent definitions and uncertainty, not forcing the same direction of result (Decena-Segarra and Rovito 2024, [doi:10.1093/molbev/msae225](https://doi.org/10.1093/molbev/msae225)).

## Figures and exact data

- `results/figures/corrected/path_analysis/corrected_path_anchor_model_comparison_analysis18_v1.png` — all anchor candidate models, weights, global fit.
- `results/figures/corrected/path_analysis/corrected_path_anchor_dag_analysis18_v1.png` — standardized anchor direct paths and approximate intervals.
- `results/figures/corrected/path_analysis/corrected_path_measurement_model_weights_analysis18_v1.png` — morphology/IOD specification sensitivity.
- `results/figures/corrected/path_analysis/corrected_path_tree_uncertainty_analysis18_v1.png` — direct paths across published time trees.
- `results/figures/corrected/path_analysis/corrected_path_leave_one_out_influence_analysis18_v1.png` — species influence on the integrated path.
- `results/figures/corrected/path_analysis/corrected_path_simulation_calibration_analysis18_v1.png` — recovery, false selection, and interval coverage.
- `results/figures/corrected/path_analysis/publication_release_gate_matrix_analysis18_v1.png` — permitted-language release matrix.
- `results/data/corrected/path_analysis/corrected_path_anchor_model_comparison_analysis18_v1.csv`, `results/data/corrected/path_analysis/corrected_path_anchor_edges_analysis18_v1.csv`, `results/data/corrected/path_analysis/corrected_path_model_stability_analysis18_v1.csv`, `results/data/corrected/path_analysis/corrected_path_tree_edge_uncertainty_analysis18_v1.csv` — review-sized exact tables.
- All raw model rankings, basis sets, fitted best edges, simulation replicates, and manifests remain in `results/data/corrected/path_analysis/`.

## Remaining work required for confirmatory promotion

1. Establish an absolute genome-size trait with a documented DNA-stoichiometric assay, same-batch reference standard, OD equation, slide/batch/specimen replication, and external validation; otherwise reframe the paper around relative IOD and morphology without genome-size claims.
2. Produce specimen/slide-held-out and final-panel-species segmentation truth with object-level split/merge/miss metrics, then rerun the morphometry from the validated production model.
3. Make specimen-balanced central tendency primary and keep largest-50 as a named upper-tail sensitivity; acquire biological replicates where only one image/specimen exists.
4. Propagate specimen-level uncertainty and TE uncertainty through an error-aware phylogenetic SEM; repeat tree and species influence analyses.
5. Treat the current DAG family as exploratory. Freeze a biologically justified confirmatory set before collecting/reprocessing the validation data.

No historical output was deleted or overwritten by this audit.
