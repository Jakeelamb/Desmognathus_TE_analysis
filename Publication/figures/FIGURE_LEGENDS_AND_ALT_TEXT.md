# Figure legends and alt text

## Figure_1

The declared panels contain 34 genomic TE species, 24 finalized phenotype/path species, and an exact 21-species overlap for integrated TE-phenotype summaries. Assembly contiguity and dnaPipeTE classification accounting are shown for TE34; the dashed line is the median unresolved order fraction.

**Alt text:** Three panels show the 34-species genomic panel, 24-species finalized phenotype and path panel, and their 21-species overlap; assembly span versus contig N50; and repeat-aligned read fraction versus unresolved TE-order mass.

## Figure_2

Order-level Shannon entropy and the Gini-Simpson index (1 − Σᵢpᵢ²) are higher when unresolved order mass is retained as an explicit category than when diversity is conditioned on classified order mass. Panel A compares these two denominators across the same 34 species; boxes show medians and interquartile ranges, and points are species. Observed richness is not plotted as a diversity index. The retry-averaged 24-feature CLR-PCA summarizes species differences among superfamilies detected in all 34 species; only the ten most distant PCA scores are labeled.

**Alt text:** Panel A shows order-level Shannon entropy and Gini-Simpson index distributions for 34 species under classified-conditional and mass-aware-unresolved denominator choices. Panel B shows species scores on the first two superfamily CLR principal components, and panel C shows variance explained by the first ten components.

## Figure_3

RepeatMasker divergence profiles vary among reported repeat categories and species. Curves are equal-weight means across all 34 species after zero-filling absent category/bin combinations; `Unclassified`, when shown, is an annotation category rather than a biological TE order. The species panel sums repeat-aligned bases in bins below 5% divergence.

**Alt text:** The left panel shows mean repeat abundance across divergence bins for the six most abundant reported categories and all other categories combined. The right panel ranks species by the fraction of aligned repeats below five percent divergence.

## Figure_4

The zero-aware terminal:internal LTR depth proxy varies among species but has uneven element support. The analysis starts from 408 LTR/Gypsy candidates at least 3,000 bp long with at least five correctly parsed TEsorter domain annotations. Points are species medians from the 380-element primary branch: ratios retained inclusively within each species' two-sided 1.5-IQR fences and elements whose left LTR, right LTR, and internal region each have at least 80% positive-depth coverage. Intervals are 2,000-replicate element-bootstrap intervals; the dashed line marks equal terminal and internal depth.

**Alt text:** A forest plot ranks species by median log2 terminal-to-internal LTR read depth and a matched bar plot shows the exact number of elements contributing to each displayed estimate.

## Figure_5

Reviewed upper-tail erythrocyte cell and corresponding nucleus areas differ among 24 finalized species. Points are species medians of all accepted reviewed pairs (up to 50 per species); triangles identify the seven species represented by one animal. Intervals are conditional paired-object bootstrap intervals and do not represent population-level among-individual uncertainty.

**Alt text:** Two aligned forest plots rank 24 Desmognathus species by cell area and show corresponding nucleus-area estimates with bootstrap intervals; point shape distinguishes one-animal from multiple-animal support.

## Figure_6

Relative nuclear IOD, nucleus area, and cell area are positively associated in Pagel-lambda PGLS sensitivity fits. Circles identify balanced multiple-image IOD support, triangles identify balanced estimates from one image and animal, and diamonds identify four species that fail the declared quality-balance gate. Bars are conditional measurement-bootstrap intervals; lines and ribbons are PGLS fits and model confidence intervals. The IOD index is not an independently calibrated absolute genome size.

**Alt text:** Three scatterplots show relative nuclear IOD versus nucleus area, relative nuclear IOD versus cell area, and nucleus area versus cell area, each with measurement intervals and a phylogenetically corrected fitted line; point shape exposes single-image and quality-balance limitations in the IOD panels.

## Figure_7

The finalized 24-species relative-IOD path-model ranking favors one Markov-equivalence class, while simulation recovery varies by generating class and residual regime. Delta CICc is relative to the best Path24 model; color and shape record the global d-separation fit gate. The IOD node is standardized log10 relative nuclear IOD, not absolute genome size, and these results do not orient causal arrows within an equivalence class.

**Alt text:** The left panel compares candidate relative-IOD equivalence classes in the finalized 24-species analysis by delta CICc and global-fit status. The right panel shows true-class recovery rates under independent and Brownian residual simulations.

## Figure_S1

Across the exact 21-species TE34-by-Path24 overlap, classified-superfamily Shannon entropy is shown against the relative nuclear-IOD index. Horizontal intervals are conditional 95% IOD intervals; the dashed line marks the Path24 median anchor of one. Triangles identify the three manually reviewed below-target Path24 inclusions. Relative IOD is an image-derived phenotype, not an independently validated absolute genome size.

**Alt text:** Scatterplot of 21 Desmognathus species comparing relative nuclear IOD with TE superfamily Shannon entropy. Three manually retained species use orange triangles; the remaining species use blue circles.
