# Reference-method audit for the cell–nucleus–genome path analysis

## Target paper

Yu, Liu, Mai, and Liao (2020), “Genome size variation is associated with
life-history traits in birds,” *Journal of Zoology* 310:255–260,
[doi:10.1111/jzo.12755](https://doi.org/10.1111/jzo.12755).

The paper analyzed 240 bird species with phylogenetic comparative regressions
and confirmatory phylogenetic path analysis. Its supporting-information record
describes candidate causal diagrams, a CICc-ranked model table (Table S5), and a
visual candidate-model set (Figure S2). The article cites `phylopath` as the
analysis implementation. The publisher exposes one Word supporting-information
file, but no article-specific analysis-code archive was listed on the article
page, in the supplement description, or in targeted GitHub searches performed
on 2026-07-10.

Primary source: [publisher article and supporting-information record](https://zslpublications.onlinelibrary.wiley.com/doi/10.1111/jzo.12755).

## Exact method being mirrored

The current `phylopath` implementation follows phylogenetic confirmatory path
analysis based on d-separation:

1. State a finite set of biologically motivated directed acyclic graphs before
   examining their ranking.
2. Convert every graph into its minimal basis set of conditional-independence
   claims.
3. Test each claim with phylogenetic generalized least squares.
4. Combine component probabilities with Fisher's C statistic. A significant
   global p-value rejects the graph's implied independence structure.
5. Rank non-rejected graphs with the small-sample C-statistic information
   criterion (CICc), report delta CICc and CICc weights, and retain model
   uncertainty when multiple graphs are competitive.
6. Estimate standardized direct paths only after the model-comparison step.

Primary sources:

- van der Bijl (2018), [`phylopath`: Easy phylogenetic path analysis in R](https://doi.org/10.7717/peerj.4718).
- von Hardenberg and Gonzalez-Voyer (2013), [Disentangling evolutionary cause-effect relationships with phylogenetic confirmatory path analysis](https://doi.org/10.1111/j.1558-5646.2012.01790.x).
- [`phylopath` official function reference](https://ax3man.github.io/phylopath/reference/index.html) and [source repository](https://github.com/Ax3man/phylopath).

The locally installed analysis version is `phylopath` 1.3.1. Its implemented
criterion is `CICc = C + 2q[n / (n - 1 - q)]`, where `q` is the node-plus-edge
parameter count. The package refuses fully connected graphs because they imply
no conditional independencies and therefore cannot be tested by d-separation.

## Consequences for the Desmognathus analysis

The frozen primary overlap contains 20 species, not 240. Small-sample correction,
global-fit checks, leave-one-species-out analysis, measurement bootstrapping,
evolutionary-model sensitivity, and actual-tree simulation are therefore primary
requirements rather than optional embellishments.

There are 25 labeled acyclic DAGs on the three nodes genome size (GS), nucleus
size (NS), and cell size (CS), but only 11 Markov-equivalence classes. Ten classes
contain testable conditional-independence claims; the six fully connected DAGs
form one untestable saturated class.

Most importantly, the three mechanisms proposed for this study are
observationally equivalent:

- GS → NS → CS
- CS → NS → GS
- NS → GS and NS → CS

All three have the same GS–NS–CS skeleton, no collider, and the same testable
claim: `GS independent of CS conditional on NS`. A cross-sectional covariance
matrix, including a phylogenetically corrected one, cannot orient those arrows.
The analysis can test whether nucleus size is the best-supported *bridge*, but it
cannot decide which of those three arrow directions is causal without additional
temporal, experimental, or externally anchored biological information.

## Figure design decision

The path figures will mirror the reference's visual grammar: white rectangular
nodes, black outlines, strong black arrows, standardized coefficients printed
beside paths, and a compact variable key. Unlike the reference slide, each figure
will also state the equivalence-class limitation, global-fit result, CICc weight,
and coefficient uncertainty so visual polish cannot be mistaken for causal
identification.
