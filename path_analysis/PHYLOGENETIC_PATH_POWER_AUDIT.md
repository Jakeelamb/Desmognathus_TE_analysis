# Statistical power audit for the final-18 phylogenetic path analysis

## Bottom line

Eighteen species is a small comparative sample. It can support a compact,
three-variable phylogenetic path analysis, but it is much better suited to
detecting large, stable relationships than to distinguishing several similar
causal structures. There is no defensible universal minimum number of species
for phylogenetic path analysis: power depends on effect sizes, the tested
independence claims, phylogenetic topology and branch lengths, phylogenetic
signal, measurement error, and the exact candidate set.

The **number of candidate DAGs does not mechanically reduce the power of an
individual PGLS path test**. Holding the basis-test orientation and model
settings fixed, adding another candidate does not change the data or the
component regression for an existing conditional-independence claim. It does,
however, increase the number of ways to obtain a competitive or spuriously
top-ranked model, can dilute CICc weights, and makes selection less decisive.
The **number of parameters inside each DAG** has a direct small-sample cost.

Accordingly, this analysis should report effect sizes and uncertainty, global
fit, CICc rankings and weights, competitive models, and simulation-calibrated
selection rates. A single top-model label or a collection of nominally
significant paths would overstate what 18 species can establish.

## What phylogenetic correction costs and protects

Species are not independent replicates. PGLS models their expected covariance
instead of treating all tips as independent. This generally protects the Type I
error rate; it does not create additional observations. The information loss is
not a single fixed discount such as “18 species equals 9 independent species.”
It is parameter- and tree-dependent. In a formal analysis of hierarchically
correlated comparative data, Ane showed that some location and lineage
parameters can retain bounded information as more related tips are added,
whereas random-covariate regression effects can remain consistently and
efficiently estimable. Thus, a generic effective sample size is not an adequate
power calculation for every path coefficient ([Ane 2008](https://doi.org/10.1214/08-AOAS173)).

The original phylogenetic confirmatory path-analysis simulations are strong
evidence for using PGLS, but **not evidence that an 18-species analysis has high
power**. Those simulations used 1,000 datasets per condition on trees of 100
species. PPA Type I error was 0.047--0.072 across phylogenetic-signal conditions,
whereas path analysis that ignored phylogeny reached a Type I error of 0.916
under Brownian signal. PPA power remained 0.951--0.973 in that particular
100-species design ([von Hardenberg and Gonzalez-Voyer 2013](https://doi.org/10.1111/j.1558-5646.2012.01790.x);
[open manuscript](https://digital.csic.es/handle/10261/65583)). These power
values cannot be transferred to the final-18 tree.

More broadly, simulations of phylogenetic model choice show that power depends
strongly on taxon count, tree shape, and parameter values, and that AIC/AICc can
have high error rates even when corrected for small samples. The authors
explicitly recommend analysis-specific Monte Carlo calibration rather than a
minimum-taxon rule ([Boettiger, Coop, and Ralph 2012](https://doi.org/10.1111/j.1558-5646.2011.01574.x);
[full text](https://pmc.ncbi.nlm.nih.gov/articles/PMC3448289/)).

## Candidate-model count versus within-model complexity

`phylopath` first fits the PGLS regressions needed for the unique basis-set
claims and then assembles each candidate's Fisher C statistic. It calculates

`CICc = C + 2q[n / (n - 1 - q)]`,

where `n` is species count and `q` is the number of nodes plus directed edges.
Relative likelihoods are `exp(-0.5 * delta_CICc)`, and weights are normalized
over every model in the supplied set. Therefore:

- adding a candidate does not enter the CICc formula for an existing DAG or
  consume residual degrees of freedom in that DAG's component regression;
- it can change every CICc **weight**, because weights have a candidate-set
  denominator;
- searching a larger, post hoc model set gives noise more opportunities to
  produce a plausible winner and introduces model-selection bias; and
- adding an edge within a DAG increases `q`, incurs a larger CICc penalty, and
  estimates another relationship from the same 18 observations.

These implementation details are explicit in the official
[`phylopath` source](https://github.com/Ax3man/phylopath/blob/fe062d28166d10fec8162903a1c67e9e83e752d8/R/phylopath.R#L114-L205)
and its
[CICc/weight functions](https://github.com/Ax3man/phylopath/blob/fe062d28166d10fec8162903a1c67e9e83e752d8/R/internal.R#L155-L163).
The package paper describes the method as confirmatory, calls candidate-set
definition the crucial step, and recommends averaging supported competitive
models rather than forcing a winner ([van der Bijl 2018](https://doi.org/10.7717/peerj.4718);
[full text](https://pmc.ncbi.nlm.nih.gov/articles/PMC5923215/)). Large
data-driven model searches are known to produce spurious selected effects and
biased post-selection inference when observations are not abundant relative to
the search space ([Lukacs, Burnham, and Anderson 2010](https://doi.org/10.1007/s10463-009-0234-4)).

For three variables at `n = 18`, the CICc penalty component alone is:

| DAG structure | `q` | Number of basis claims | CICc penalty |
|---|---:|---:|---:|
| no edges | 3 | 3 | 7.71 |
| one edge | 4 | 2 | 11.08 |
| two edges | 5 | 1 | 15.00 |
| saturated three-edge DAG | 6 | 0 | not testable by d-separation |

Thus, adding the first and second edges costs about 3.36 and 3.92 CICc units,
respectively, before considering improvement in Fisher's C. Conversely, a
denser DAG makes fewer conditional-independence predictions and is less
falsifiable. The saturated class makes none and cannot be ranked by the
d-separation method. Shipley's derivation establishes information-criterion
comparison for d-separation models but also states that path analysis assumes a
sufficient sample; it does not supply a universal minimum
([Shipley 2013](https://doi.org/10.1890/12-0976.1)).

## Direct evidence from this project at `n = 18`

The most relevant power evidence is the repository's analysis-specific
calibration, not the 100-species benchmark above. The script
[`simulate_corrected_path_calibration.R`](../scripts/processing/simulate_corrected_path_calibration.R)
simulates traits on the published final-18 time tree, uses Brownian residual
processes, and refits the actual candidate sets. For the analogous three-trait
`iod_morphology` family (four candidates), the frozen results in
[`corrected_path_simulation_calibration_analysis18_v1.csv`](../results/data/corrected/path_analysis/corrected_path_simulation_calibration_analysis18_v1.csv)
show:

| Generating scenario | Replicates | Expected model ranked first | Expected model uniquely supported |
|---|---:|---:|---:|
| independent null | 200 | 54.0% | 18.0% |
| observed chain, 0.5x effects | 100 | 33.0% | 19.0% |
| observed chain, 1.0x effects | 100 | 93.0% | 81.0% |
| observed chain, 1.5x effects | 100 | 100.0% | 85.0% |

Under the independent null, a non-null model ranked first in 45.5% of
simulations and was a uniquely supported false non-null in 20.0%. These values
are properties of this candidate set, calibration rule, tree, and simulation
model--not universal false-positive rates. They demonstrate the key point:
with 18 species, model recovery is weak for half-sized effects and reasonably
strong for observed-sized effects, but even observed-sized effects do not make
selection certain.

## Reporting and analysis recommendations

1. Keep the primary candidate set small, biologically justified, and specified
   before inspecting its ranking. Treat any exhaustive set as a sensitivity or
   identifiability audit, not as ten independent confirmatory tests.
2. Do not Bonferroni-correct merely because several DAGs are ranked; CICc model
   comparison is not ten independent null-hypothesis tests. Do control the
   scientific search space and avoid post hoc model proliferation.
3. Report every globally supported model within `delta_CICc <= 2`, its weight,
   and the stability of rank/weight across measurement, tree, and leave-one-out
   analyses. Do not interpret failure to reject a DAG as proof that it is true.
4. When models are competitive, use supported-model averaging where the target
   path has a coherent interpretation, or report the competing structures
   without selecting one. Exclude globally rejected models from averaging, as
   recommended by `phylopath`.
5. Report standardized coefficients with bootstrap or simulation intervals.
   Separate uncertainty in the path estimates from uncertainty in model
   selection.
6. Retain the actual-tree parametric calibration and expand it around
   scientifically meaningful effect sizes. Report correct-model top rate,
   competitive-set inclusion, unique-support rate, and null false-selection
   rate with Monte Carlo intervals.
7. Phrase the final result as exploratory model support. With 18 species and
   observational cross-species data, the analysis can distinguish broad
   covariance structures when effects are strong; it cannot guarantee causal
   direction, especially for Markov-equivalent DAGs.
