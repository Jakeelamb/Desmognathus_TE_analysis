# Body size and cellular traits: isolated exploratory analysis

All new files are confined to this manuscript run. Existing figures, data and notes remain separate.

## Data and transformations

The 24 species come unchanged from the run's frozen Supplementary Data S04. Exact tree_tip matches join the canonical organismal_traits.csv species key; no taxon alias is newly inferred. All 24 have positive complete genome size, cell area, nucleus area and maximum SVL. Sources and SHA256 hashes are copied into inputs/body_size before fitting. Canonical and frozen traits and trees are checked for equality.
SVL is collaborator-provided species maximum snout-vent length (mm), not the body length of the animals supplying microscopy. Genome size is the existing IOD-derived, D. fuscus assembly-anchored estimate. Cell and nucleus size are the existing species medians among retained upper-tail paired erythrocyte objects (areas in square micrometres), not unbiased species-wide means. No specimen-level body-size relationship is inferred.
Every trait is log10-transformed and standardized using the full24 mean and sample SD. The subset22 retains exactly that scaling, so displayed coefficients use a common reference. scaling.csv preserves the constants. The comparison excludes exactly D. aeneus and D. wrighti, the two smallest sampled species and the two direct developers. This is a joint influence/life-history sensitivity, not an independent small-body group, not a test that distinguishes body size from development, and not an outlier cleanup. No threshold or species search was performed.

## Regressions

b = body size; g = genome size; c = cell area; n = nucleus area. Regressions in each panel are g~b, c~b, n~b; c~g, n~c, n~g; c~g+b, n~c+b, n~g+b; n~g+c; and n~g+c+b. Pagel-lambda regressions use phylolm with lambda constrained to [1e-7,1], maximum likelihood, no measurement-error model. Lambda bounds and warnings are retained. Confidence intervals are beta +/- qt(.975,n-p)*SE, with p the number of regression coefficients including the intercept, matching phylolm's t-test degrees of freedom. These intervals are conditional approximations and do not integrate uncertainty in estimated lambda.
Pairwise curves use the fixed-effect regression mean on the standardized log scale, with variance X V X'. Curves and intervals are back-transformed without bias correction, representing the fitted conditional median on the original trait scale, not individual prediction intervals. Their range is restricted to the SVL range observed in each panel. coefficient_covariance.csv retains the full coefficient covariance matrices.

## Four prespecified hypotheses

All four hypotheses impose b->g, b->c and b->n. body_only adds no cellular edge; cell_bridge adds g->c and c->n; nucleus_bridge adds g->n and n->c; genome_bridge adds g->c and g->n. These are a finite exploratory comparison conditional on treating body size as a common predictor. The set is neither exhaustive over four-node DAGs nor a test proving body size exogenous. Each graph represents a Markov-equivalence class; the directions are hypotheses, not uniquely identified causal effects. A saturated graph is not included because it has no testable d-separation claim.
phylopath::phylo_path and its summary implement d-separation, Fisher C, q and CICc. The package-selected d-separation regression orientation and estimated lambda are exported. k is the number of tested independence claims; q is the package's parameter count. Weights normalize within each candidate set and species panel, not across them. Absolute CICc values must not be compared between full24 and subset22 or between the four- and three-trait sets. Good global fit does not validate causal direction.
Component path coefficients are fit with phylolm from each candidate's structural equations to preserve common full24 standardization. phylopath::est_DAG is intentionally not used for coefficients because it re-standardizes numeric variables in each subset. All candidate edges are exported, including weak or uncertain edges. No model-averaged causal effect is reported.
The original ten testable three-trait equivalence classes are also rerun in each panel to distinguish changed species composition from body adjustment. The full24 ranking is numerically checked against frozen S30 before verification passes.

## Limits and verification

This is an exploratory, species-level analysis with n=24 or22. Genome-size calibration uncertainty, morphology measurement uncertainty, SVL uncertainty, tree uncertainty, and within-species sampling error are not propagated here. Direct-developer removal simultaneously changes body-size range, development composition and taxon sampling. It cannot attribute any change uniquely to one cause. Predictions and t intervals, exact membership, positive complete traits, common scaling, source preservation, tree matching, Fisher C, weights and full24 reproduction are checked by assertions in the analysis script. Warnings are printed and saved without user-code suppression.

Reproduce from the run directory:

    PATH=/home/jake/miniconda3/envs/desmognathus/bin:$PATH Rscript --no-environ scripts/body_size_analysis.R

CSV schemas: regression_coefficients.csv includes intercepts and uses standardized_coefficient, standard_error, ci_low, ci_high, p_value, fitted_phylogenetic_parameter, n_species and residual_df. path_edge_estimates.csv adds candidate_set, model, parent and child. path_model_comparison.csv identifies both candidate_set and panel. pairwise_predictions.csv contains max_svl_mm, fitted, ci_low and ci_high on raw scales.
