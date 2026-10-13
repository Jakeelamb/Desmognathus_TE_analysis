# Exploratory LTR ratio / estimated genome-size association

Nineteen species have both released genome-size estimates and primary LTR ratios (239 retained elements); 11 LTR species without genome sizes are excluded without imputation. Each hollow point is the released median of untransformed element ratios for one species. Gray bars are estimated genome-size intervals and 95% within-species element-bootstrap median intervals. Pagel-lambda PGLS fits log10(median ratio) as the response to log10(genome Gbp), by maximum likelihood with lambda constrained to [0,1], on the original branch lengths of the source tree pruned to 19 species. Each species contributes one observation; elements are not biological replicates and do not act as regression weights. The blue line is the back-transformed fitted geometric mean; its band is a 95% t-based coefficient interval, not an individual-species prediction interval. The dotted ratio-one line denotes equal terminal and internal depth. Slope 0.024196 (95% CI -0.560865 to 0.609257), P = 0.931489; lambda = 6.2e-09 at a parameter boundary. No detectable association in these 19 species; the interval permits effects in either direction.

Leave-one-species-out diagnostics are retained as analysis tables; the author omitted their plot from the review gallery. Each diagnostic fit omits one species, prunes the tree, and re-estimates lambda with the remaining 18 species. Leave-one-out slopes range from -0.110282 to 0.156522; every P value is at least 0.542409. These deletions assess influence, not independent replications or 19 additional primary hypothesis tests.

Exploratory companion, not a manuscript figure assignment. Genome sizes are D. fuscus assembly-anchored IOD-derived estimates converted at 0.978 Gbp per pg, not independent species C-values. The LTR ratio is a mapping/deletion-footprint proxy, not a deletion rate. The supplied dated tree's publication and calibration provenance remain unresolved. Coefficient intervals and P values are conditional on fitted lambda and the supplied tree; the fit does not propagate measurement uncertainty in either variable.

## Model specification

Primary model: log10(ratio_median) ~ log10(genome_size_gbp); phylolm, model=lambda, ML, bounds [0,1].
Slope is a log-log regression coefficient, not a Pearson correlation coefficient. A 1% increase in genome size corresponds approximately to a slope-percent change in fitted ratio.
The medians are calculated on the original ratio scale before taking logarithms. All data are positive; no pseudocount is added.
Standard errors use phylolm's coefficient covariance. Two-sided t tests and 95% t intervals use n-2 residual degrees of freedom, conditional on estimated lambda. The fitted lambda reaches its lower boundary; this is reported rather than suppressed as a successful interior estimate.
No element-count weighting, measurement-error model, causal interpretation, or uncertainty over tree topology/calibration is included.
Model sensitivities (predefined before plotting results): Brownian covariance (lambda=1), ordinary regression, and an untransformed PGLS. They are recorded in model_summary.csv and are not selected by significance.
The 19 leave-one-species-out fits each re-estimate lambda, retain source branch lengths, and use 18 species. All outputs are downstream derivatives; canonical publication tables remain unchanged.

## Files

species_data.csv: exact included species, original medians, counts, converted genome values, both uncertainty intervals, and transformed variables.
excluded_species.csv: 11 LTR species without a genome estimate, with explicit exclusion reason.
model_summary.csv: primary and sensitivity slopes, t intervals, P values, lambda boundaries, and captured model warnings.
leave_one_species_out.csv: all 19 influence fits; fitted_relationship.csv: plotted geometric fit and confidence ribbon.
pgls_19_species.tre: actual pruned tree with source branch lengths; input_checksums.csv and session_info.txt: provenance.

Method reference: [phylolm package documentation](https://search.r-project.org/CRAN/refmans/phylolm/html/phylolm.html); [coefficient inference and prediction implementation](https://github.com/lamho86/phylolm/blob/master/R/phylolm.R).
