# Publication-readiness remediation status

**Updated:** 2026-07-09
**Preservation rule:** no historical analysis, figure, or derived table has been deleted. Corrected products use new paths until they pass scientific and output-parity review.

| Gate | Current status | Evidence | Remaining before promotion |
|---|---|---|---|
| *D. fuscus* resource identity | Decision encoded | `source_manifest.csv`, `species_taxonomy_crosswalk.csv`, and provenance tests distinguish the current reidentified resource, excluded untrusted draft, and 2025 validation assembly | Propagate `te_*accession` fields during the next full derived-table rebuild |
| Genomic versus microscopy samples | Guarded | Data dictionary and tests prohibit genomic accession fields in share-facing CellProfiler tables | Keep this namespace rule in the release notebook |
| Current versus 2025 *fuscus* comparison | Diagnostic benchmark complete | `fuscus_resource_benchmark.csv` and `.md` preserve methods and denominators | Obtain missing EarlGrey/LTR_FINDER assets from Alex for annotation-level validation |
| RepeatMasker hit classification | Corrected 18-species table complete; historical outputs preserved | Shared classifier, 9.53-million-hit audit, and `repeatmasker_detailed_classification_hit_level_analysis18_v1.csv` with 5,132,397 conserved hits and independent readback verification | Rebuild divergence summaries into a separate output branch |
| RepeatMasker divergence/gap summaries | Corrected 18-species sensitivity complete; historical output preserved | Corrected summary conserves all 5,132,397 hits in an all-hit branch and 4,715,646 context-eligible hits in the legacy-threshold branch; feature sensitivity is recorded in `repeatmasker_corrected_feature_sensitivity.md` | Keep gap metrics out of direct DNA-loss language; decide whether the legacy dnaPipeTE threshold belongs in the confirmatory model |
| RepeatMasker divergence landscape | Corrected compact cache and figure complete | One-percent bins conserve 5,132,397 hits and 1,155,346,511 inclusive query-coordinate aligned bp; PNG/PDF small multiples use a within-species aligned-hit-bp denominator | Do not label the y-axis genome or assembly percentage; overlapping hit intervals are not deduplicated |
| Current path outputs | Corrected exploratory branch complete | 958 fits across data specifications, 201-tree sensitivity, and species omissions completed without fit failures; every basis set/ranking/edge is exported and fail-closed | Do not promote to a genome-size or causal claim until measurement and prospective-model gates pass |
| dnaPipeTE mass and denominator | Relative mass complete; configured absolute-load sensitivity recovered | Git-history parameters recover a 1.5-Gb quantification denominator; the 18-species table conserves TE/Other/Unknown/unresolved mass and gives 60.619%–68.196% repeat-aligned load | Keep sensitivity-only until runtime logs, container/library hashes, and sampling uncertainty are recovered or the final panel is rerun reproducibly |
| Terminal:internal LTR depth | Corrected robustness/QC branch complete; biological validation blocked | 567 exact elements across 16 panel species, zero-aware means, four robust branches, bootstrap intervals, influence diagnostics, and four validation figures; 565/567 pass regional coverage | Recover mapping provenance and validate against structural solo:intact-LTR calls; do not use the arithmetic mean or call the ratio an ectopic-recombination rate |
| TE diversity indices | Descriptive final-panel branch approved | Historical tables reproduce to <5e-16 as Shannon, Gini-Simpson, and Pielou; corrected output adds richness, dominance, Hill numbers, mass sensitivity, and exact resource identity | Use classified-conditional values with unresolved-bin sensitivity; add phylogenetic regression only after the tree gate |
| TE order CLR PCA | Descriptive ordination approved | 18 species x 10 nonzero order features; PC1/PC2 56.7%/28.9%; minimum LOO score correlations 0.990/0.953 | Keep descriptive unless a phylogenetic ordination/direct multivariate sensitivity passes |
| TE superfamily CLR PCA | Supplementary/sensitivity only | 27 features, 32 zeros; zero-replacement changes variance and PC2 LOO score stability falls to 0.498 | Do not promote superfamily axes to the confirmatory path set |
| Clade PERMANOVA | Not approved | Historical unrestricted species-label permutations violate phylogenetic exchangeability | Replace with a phylogenetic multivariate model or retain as explicitly noninferential sensitivity |
| Final-panel focal time tree | Structurally approved; focal source provenance open | Exact 18-tip prune, no substitutions, positive/bifurcating branches, two-year maximum rounding correction, hashed tree/crosswalk/patristic outputs | Obtain collaborator tree citation, original identifier, calibration record, and tree type |
| Published tree uncertainty | Approved sensitivity set | Stewart-Wiens 2025 main tree plus 200 time-calibrated bootstrap trees; all 18 exact taxa; crown-age median 13.859 Myr (95% 12.166-16.403); one rooted bipartition differs from focal | Propagate this set through final path models; state that bootstrap trees do not model reticulation |
| TE order phylogenetic PCA | Approved descriptive across tree-time uncertainty | Fixed-tree pPCA and all 200 published time trees retain order geometry; minimum ordinary-vs-pPCA score correlations 0.988/0.910 and loading correlations 0.955/0.723 | Keep axes descriptive; use declared log-ratios/direct multivariate phylogenetic models for mechanistic inference |
| Morphology estimand and review state | Upper-tail sensitivity complete; validation blocked | Six estimators, rank stability, hierarchical intervals, review coverage, selection-depth figures, and 900 linked pairs are frozen; 574/900 are model-ranked/unlabeled | Finish blinded review, add specimen-balanced central tendency, and validate the production masks on held-out final-panel species |
| Image-IOD genome size | Direct C-value interpretation blocked; conditional calibrated estimates complete | The frozen quality-matched panel contains 721 reviewed nuclei from 20 species and reports equal-image IOD ratios rescaled to *D. fuscus* = 16.36 pg; the anchor provenance remains unresolved and no same-batch DNA standard was imaged | Present these only as conditional *fuscus*-anchored genome-size estimates; obtain a DNA-stoichiometric same-batch reference and external validation before calling them direct C-values |
| Phylogenetic path analysis | Exploratory corrected audit complete | Explicit-null candidate sets, basis sets, global fit, CICc/weights, standardized edges, 200-tree uncertainty, leave-one-out influence, and actual-tree simulation are frozen | Propagate specimen/TE measurement error and freeze a prospective DAG set after obtaining a valid independent genome-size trait |
| Eight-notebook research-review bundle | Executed and fail-closed | Eight independent notebooks, optional microscopy HTML viewers, portable PNG/PDF figures, and SHA-256/execution/error gates are recorded in `notebooks/research_review/research_review_manifest.json` | Use the per-notebook claim boundaries during collaborator review; do not promote blocked measurement or causal claims |

## Verified commands

```bash
scripts/run_tests.sh
scripts/run_in_dusky.sh python scripts/processing/build_fuscus_resource_benchmark.py
scripts/run_in_dusky.sh python scripts/processing/audit_repeatmasker_hit_classification.py
scripts/run_in_dusky.sh python scripts/processing/rebuild_repeatmasker_hit_classification.py \
  --output results/data/corrected/repeatmasker_detailed_classification_hit_level_analysis18_v1.csv
scripts/run_in_dusky.sh python scripts/processing/build_corrected_divergence_summary.py
scripts/run_in_dusky.sh python scripts/processing/audit_dnapipete_mass_accounting.py
scripts/run_in_dusky.sh python scripts/processing/build_dnapipete_absolute_load_sensitivity.py
scripts/run_in_dusky.sh python scripts/processing/build_corrected_ectopic_recombination.py
scripts/run_in_dusky.sh python scripts/processing/build_corrected_te_diversity_pca.py
scripts/run_in_dusky.sh python scripts/processing/audit_phylogeny_release.py
scripts/run_in_dusky.sh Rscript scripts/processing/audit_published_tree_uncertainty.R
scripts/run_in_dusky.sh Rscript scripts/processing/audit_order_phylogenetic_pca.R
scripts/run_in_dusky.sh Rscript scripts/processing/audit_order_pca_tree_uncertainty.R
scripts/run_in_dusky.sh python scripts/processing/audit_microscopy_release.py
scripts/run_in_dusky.sh python scripts/processing/build_corrected_path_inputs.py
scripts/run_in_dusky.sh Rscript scripts/processing/audit_corrected_path_models.R --phase all
scripts/run_in_dusky.sh Rscript scripts/processing/simulate_corrected_path_calibration.R
scripts/run_in_dusky.sh python scripts/processing/summarize_corrected_path_audit.py
scripts/run_in_dusky.sh python scripts/processing/build_publication_audit_notebook.py
scripts/run_in_dusky.sh jupyter nbconvert --to notebook --execute notebooks/Desmognathus_publication_audit_analysis18_v1.ipynb --inplace --ExecutePreprocessor.timeout=600
scripts/run_in_dusky.sh python scripts/processing/audit_publication_notebook.py
scripts/run_in_dusky.sh python scripts/processing/build_corrected_repeat_landscape.py
scripts/run_in_dusky.sh python scripts/processing/build_research_review_notebooks.py
scripts/run_in_dusky.sh Rscript scripts/processing/build_audited_historical_style_figures.R
scripts/run_in_dusky.sh Rscript path_analysis/scripts/run_cell_nucleus_genome_phylogenetic_path_analysis.R
scripts/run_in_dusky.sh python path_analysis/scripts/build_cell_nucleus_genome_path_presentation.py
scripts/run_in_dusky.sh python scripts/processing/audit_research_review_notebooks.py
```

The corrected RepeatMasker writer refuses to overwrite
`results/data/repeatmasker_detailed_classification_combined.csv`. It first
writes a `.partial` file and promotes it only after all rows complete.

The dnaPipeTE mass-accounting branch is restricted to the declared 18-species
panel. It excludes *D. orestes* and every other out-of-panel staged resource,
preserves aligned-repeat totals and unresolved mass at each taxonomic level,
and labels its percentages as relative dnaPipeTE aligned-repeat composition.

## Resource request for the 2025 assembly

Ask Alex for the EarlGrey `*_summaryFiles/` directory, filtered-repeat GFF,
high-level count table, divergence summary, combined repeat library, raw
LTR_FINDER GFF, EarlGrey 4.4.2 command/log, and checksums. The public Zenodo
record currently exposes only the compressed FASTA and AGP.
