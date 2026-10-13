# Yu et al. (2020): path-analysis methods audit

Checked 24 September 2026. This is a read-only literature audit; no figures, observations, or models were changed.

## What is verified for this exact paper

The intended DOI is [10.1111/jzo.12755](https://doi.org/10.1111/jzo.12755): Yu, Liu, Mai and Liao, *Genome size variation is associated with life-history traits in birds*, *Journal of Zoology* 310:255–260 (2020), first published online 25 December 2019. The publisher describes a comparative study of 240 bird species, combining phylogenetic comparative analysis with path analysis. [Publisher record](https://zslpublications.onlinelibrary.wiley.com/doi/10.1111/jzo.12755)

The publisher's supplementary-material description verifies **PGLS** analyses and **phylogenetic path analysis of genome size, body mass and hatching time**. Supplementary Table S5 ranks candidate models by **CICc**; Figure S2 depicts the candidates; Figure S1 contains the 240-species phylogeny. The bibliography includes the original phylogenetic confirmatory path-analysis paper and the `phylopath` software paper. [Publisher, Supporting Information and References](https://zslpublications.onlinelibrary.wiley.com/doi/10.1111/jzo.12755)

## What remains unverified

The accessible publisher record does not expose the main Methods or the path-figure legend. Therefore this audit cannot certify:

- The exact number or identities of candidate DAGs and whether reverse causal directions were included.
- Whether all 240 birds had complete data for the path analysis.
- The software/function/version actually run; a software citation alone does not prove its use.
- Log transformation, standardization, chosen evolutionary covariance model, or phylogeny construction/summarization.
- Whether the illustrated graph is one selected model or a model average.
- The coefficients' exact estimation method, uncertainty method, bootstrap count, or numerical model weights. Numerical arrow labels are now visible in the user-supplied screenshot described below.

Do not fill these gaps by attributing ordinary `phylopath` defaults to Yu et al.

The public [supporting document link](https://zslpublications.onlinelibrary.wiley.com/doi/suppl/10.1111/jzo.12755/supinfo/jzo12755-sup-0001-Supinfo.doc) returned a Cloudflare 403 to direct retrieval; web-tool retrieval also failed. Publisher full-text/PDF routes were inaccessible through the web tool. [OpenAlex's DOI record](https://api.openalex.org/works/doi:10.1111/jzo.12755) listed no open-access location; a lawful author/repository full-text copy was not found in the bounded search. These are access findings, not evidence that the material does not exist.

## The method family explains the attractive diagrams

The primary software paper describes phylogenetic confirmatory path analysis using DAG-implied conditional-independence tests and phylogenetic regressions. Its summary separates goodness of fit (Fisher's C and its P value), CICc ranking, and model weights. It then estimates **standardized regression coefficients for individual paths**. Bootstrap confidence intervals and averaging across candidate models are available but require choices by the analyst. [van der Bijl 2018](https://pmc.ncbi.nlm.nih.gov/articles/PMC5923215/)

These are two distinct quantities:

| Quantity | Meaning | Suitable display |
|---|---|---|
| Model weight | Relative support for an entire candidate graph within the tested set | Beside the graph or in a ranking panel |
| Standardized path coefficient | Size and sign of an association conditional on that graph's specified predictors | Number beside an arrow, with arrow width encoding magnitude |

The software author's tutorial explicitly demonstrates arrow labels and widths based on standardized coefficients, and discusses selecting one model or averaging competing models. This verifies how to make that display; it does **not** identify Yu et al.'s particular settings. [Official phylopath tutorial](https://ax3man.github.io/phylopath/articles/intro_to_phylopath.html)

## Implication for the manuscript figure

An arrow diagram with quantitative path estimates can make the current cellular hypothesis easier to understand. A prominent diagram can show a biologically motivated direction **as an explicitly assumed causal model**, accompanied by standardized coefficients and intervals. Competing, observationally equivalent orientations should remain visible nearby or in the supplement, with their shared model-class support clearly distinguished from direction-specific evidence. Neither a large model weight nor a large path coefficient should be relabelled a probability that one arrow is causally correct.

This last paragraph is an analytical recommendation, not a finding reported by Yu et al. Exact imitation of their statistical choices requires their Methods, main path-figure legend, and Table S5/Figure S2.

## Additional visual evidence supplied by the user

The user subsequently supplied `/tmp/codex-clipboard-0be9ede5-87c6-4700-80d1-839668087e9a.png`, showing the article heading and a seven-variable network. Labels define age at sexual maturity (ASM), body mass (BM), maximal longevity (ML), hatching time (HT), clutch size (CS), fledging time (FT), and genome size (GS). The diagram contains signed numerical labels, black positive paths, red negative paths, and visibly different arrow widths. For example, BM -> GS is labelled 0.80, HT -> GS 0.07, and ML -> GS -0.45. It retains multiple direct and indirect paths in one connected diagram.

This screenshot supports those visual observations, but does not supply the estimation methods or full figure legend. It therefore does not resolve standardization, model averaging, uncertainty settings, or candidate-model selection. The separate local triangular candidate borrows the display idea while using our own verified S28 regressions; it does not claim to replicate the paper's exact fitting procedure.
