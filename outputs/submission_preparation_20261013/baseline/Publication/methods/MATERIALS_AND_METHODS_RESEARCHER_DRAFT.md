# Materials and Methods

## Study design and data sources

Analyses were conducted on data-complete species panels defined for each data type. The transposable-element analyses comprised 34 species, the microscopy and phylogenetic-comparative analyses comprised 24 species, and 21 species occurred in both panels. *Desmognathus brimleyorum*, *D. folkertsi*, and *D. ochrophaeus* were represented in the microscopy panel but lacked transposable-element data, and no values were imputed for these species. Genomic and microscopy measurements were joined at the species level and were not assumed to have been obtained from the same animals.

Genome assemblies and raw sequencing resources were obtained from the National Center for Biotechnology Information Assembly and Sequence Read Archive databases. Accession-specific identity records were used to reconcile source names with the analysis taxonomy. Two dnaPipeTE outputs were available for *D. orestes*. These represented technical runs rather than biological replicates, so each run was normalized separately and the two compositions were averaged with equal weight.

Maximum snout–vent length, body-size, and life-history data were provided by Alex Pyron. Maximum SVL was used rather than total length; for D. orestes, the larger of two reported values (52 and 50mm) was retained. Reproductive strategy was categorized as direct development, aquatic eggs, or terrestrial eggs with aquatic larvae.

## Genomic read processing and de novo TE library construction

Mitochondrial reads were removed by indexing a mitochondrial reference (genbank [AY728227.1](https://www.ncbi.nlm.nih.gov/nuccore/AY728227.1/), aligning paired reads with BWA-MEM using eight threads, converting and sorting the alignments with SAMtools, and applying samtools view -f 12 -F 256 to retain pairs for which both mates were unmapped to the mitochondrial reference. These nonmapping pairs formed the nuclear-read input for subsequent TE analyses.

A *Desmognathus*-specific repeat library was generated with RepeatModeler. The workflow used the TEtools 1.87.sif container under Singularity 3.7.4, built a database for each assembly, and ran RepeatModeler with 32 threads, -LTRStruct, and -genomeSampleSizeMax 10000000000. Redundant consensus sequences were clustered with CD-HIT-EST using -d 0 -aS 0.8 -c 0.8 -G 0 -g 1 -b 500 -T 64 -M 500000. Subsequent dnaPipeTE analyses used the custom library as its TE database.

## TE abundance estimation and classification

Transposable-element abundance was estimated from the nuclear-filtered short reads with dnaPipeTE. The workflow used Singularity 4.1.2, 16 CPUs, mate 1 of the nuclear-filtered reads, a configured 15-Gb genome-size denominator, 0.1× target coverage, RM_t = 0.15, two internal samples, and custom TE library. The 15-Gb value was used only to configure read sampling and was not interpreted as a measured genome size or absolute TE load.

The dnaPipeTE \_reads_per_component_and_annotation output was parsed to assign repeat components to class, order, and superfamily. Within each taxonomic level, aligned bases assigned to each classified TE category were divided by the total classified aligned bases for that species. Total aligned-base mass, retained classified mass, and unresolved mass were preserved as quality-control quantities. Unresolved mass did not enter the diversity calculations. RepeatMasker .align files used for divergence profiles described dnaPipeTE Trinity repeat-contig assemblies rather than the accession-linked genome assemblies.

## TE diversity and compositional ordination

TE diversity was calculated at the superfamily level following the top-ten scope used by Wang et al. (2025). For each species, the ten most abundant strictly positive TE superfamilies were selected independently; abundance ties were resolved by superfamily name. The selected abundances were reclosed to sum to one before natural-log Shannon entropy, H' = −Σ p<sub>i</sub> ln(p<sub>i</sub>), and the Gini-Simpson index, 1 − Σ p<sub>i</sub><sup>2</sup>, were calculated, where p<sub>i</sub> is the reclosed proportion of classified aligned bases assigned to superfamily i. Wang et al. (2025) do not state whether their retained top-ten abundances were reclosed; reclosure is therefore specified explicitly here.

PCA was conducted on the classified superfamily compositions of the 34 TE species. After technical-run averaging, superfamilies with positive abundance in every species were retained, leaving 24 shared superfamilies. Ginger, Merlin, Dada, Chapaev, and CR1 were excluded from the ordination because they were not detected in every species. Each retained composition was reclosed to sum to one and transformed by the centered log ratio, CLR(*p*<sub>i</sub>) = ln(*p*<sub>i</sub>) - mean\[ln(*p*)\]. The retained matrix contained no zeros, so no zero replacement was applied. CLR columns were centered without variance scaling, and principal components were obtained by singular-value decomposition. Axis signs were oriented deterministically so that the loading with the largest absolute value on each axis was positive.

## Repeat divergence profiles

Repeat divergence profiles were reconstructed from the RepeatMasker .align files associated with the dnaPipeTE Trinity repeat-contig assemblies. Inclusive aligned hit length was calculated as \|query end - query start\| + 1. Rows with invalid divergence or length values were removed, and reported divergence was assigned to one-percentage-point bins with a terminal 50% bin. Within each species, aligned bases in each TE category and divergence bin were divided by the total inclusive aligned RepeatMasker hit bases for that species. Missing species-by-category-by-bin combinations were filled with zero before equal-weight averaging across the 34 species.

These fractions describe the composition of aligned repeat annotations rather than genome occupancy, insertion counts, insertion ages, or evolutionary rates. Unclassified repeats were retained as an annotation category rather than treated as a biological TE order.

## LTR terminal-to-internal depth proxy

DANTE, LTRharvest, and LTRdigest were used to identify candidate full-length long terminal repeat (LTR) elements and annotate their 5' LTR, internal region, and 3' LTR. Candidate sequences were joined to TEsorter annotations by exact sequence identifier and coordinates. Elements were eligible when they were at least 3,000 bp long, were classified by TEsorter as Order = LTR and Superfamily = Gypsy, and contained at least five whitespace-delimited DOMAIN\|MODEL annotations. Correct parsing produced 408 eligible LTR/Gypsy elements across 30 species.

Trimmed, mitochondrial-filtered paired reads from each species were mapped to the corresponding LTR contigs with Bowtie2 using -L 20 --very-sensitive-local. SAMtools was used to convert, sort, and index the alignments and to calculate per-position depth. Positions absent from the depth output were assigned a depth of zero. Mean terminal depth was calculated as the length-weighted mean across the 5' and 3' LTRs, and mean internal depth was calculated over the annotated internal region. The terminal-to-internal ratio was obtained by dividing mean terminal depth by mean internal depth. The archived records do not preserve the Bowtie2 and SAMtools versions or the policies for mapping quality, multimappers, secondary or supplementary alignments, and duplicate reads.

Before species-level summaries were calculated, an inclusive two-sided Tukey filter was applied independently within each species to the untransformed element ratios. Elements were retained when Q1 − 1.5 × IQR ≤ ratio ≤ Q3 + 1.5 × IQR. The quartiles and fences were calculated from the 408 eligible elements before any coverage filtering, and 381 elements passed the IQR rule. An independent coverage gate required positive depth at 80% or more of the expected positions separately in the 5' LTR, 3' LTR, and internal region; 407 elements passed this gate. The intersection of the two filters contained 380 elements across all 30 species. Species were summarized by the median ratio with 2,000 element-bootstrap replicates.

The terminal-to-internal statistic was used as an indirect mapping and deletion-footprint proxy motivated by the hypothesis that recombination between terminal repeats can remove internal sequences.

## Specimen preparation, staining, and image acquisition

Each microscopy specimen identifier represented one individual animal; multiple images and measured objects from the same specimen were repeated observations rather than independent biological replicates. All animals were ethanol-fixed before erythrocyte collection. Erythrocytes were collected from an incision through the heart with a pipette and deposited at the center of a microscope slide. Slides were stained with Feulgen-Schiff reagent following Hardie et al. (2002) and imaged on an Olympus APX100 with a 40× objective.

## Cell and nucleus segmentation and morphometry

Cell boundaries were segmented with Cellpose and nuclei were segmented with a custom trained YOLO model. Masks retained source-image, tile, centroid, overlap, and object-level identifiers. Objects were eligible when they formed a valid one-to-one cell-nucleus pair, did not touch an image edge, were not overlapping other cells or debris, and passed manual review.

Eligible cells were ranked within species by descending cell area. The target was the 50 largest accepted cells per species, and rejected objects were replaced by the next-largest eligible accepted cell. This procedure estimates the upper tail of detected eligible erythrocyte size rather than the population-average cell. The final morphometry dataset contained 1,152 reviewed cell-nucleus pairs across 24 species and 42 animals. Twenty-two species contributed 50 cells, *D. aeneus* contributed 44, and *D. wrighti* contributed 8.

Cell and nucleus areas were calculated from mask pixels using a linear calibration of 0.12 micrometers per pixel, equivalent to 0.0144 square micrometers per pixel. Species-level cell and nucleus estimates were medians of the accepted objects. Uncertainty intervals were obtained by paired-object bootstrap resampling conditional on the observed images and animals. These intervals do not estimate population-level variation among animals.

## Nuclear IOD measurement and conditional genome-size calibration

Pixel optical density was calculated as log10(background intensity / pixel intensity), and nuclear integrated optical density (IOD) was calculated by summing pixel optical density over the nuclear mask. Images were matched using log edge sharpness and log relative ring noise; IOD, darkness, nucleus area, and cell area were excluded from the matching variables. Balance was assessed with a maximum absolute standardized mean difference threshold of 0.10 and a maximum Kolmogorov-Smirnov distance threshold of 0.25. Species that failed either threshold remained in the release and were flagged rather than removed.

The final IOD dataset contained 805 reviewed nuclei from 51 images, 50 animals, and 24 species. For each image, the median per-nucleus IOD was calculated, and the species estimate was the arithmetic mean of the image medians so that images rather than nuclei received equal weight. Uncertainty was estimated with 2,000 hierarchical-bootstrap replicates that resampled images and then nuclei within images. Dividing each species estimate and interval by the median estimate across the 24-species panel produced a relative-IOD index used to document the measurement and calibration path.

Two *D. fuscus* reference slides, prepared from specimens 32469 and 32470, served as DNA-content standards and were stained in the same experimental runs as the unknowns. Conditional genome size was calculated as target genome size (pg/1C) = (target equal-image IOD / *D. fuscus* equal-image IOD) × 16.462167689 pg/1C. The reference value was obtained from the selected pooled *D. fuscus* assembly estimate as 16.1 Gbp / 0.978 Gbp per pg = 16.462167689 pg/1C (Myers et al. 2025).

## Phylogenetic comparative analyses

The focal time-calibrated *Desmognathus* phylogeny was supplied by Alex Pyron. Taxon names were reconciled with the study species, exact unsuffixed representative tips were retained, and alternate tips and species without complete data were pruned with ape::drop.tip. Single-descendant internal nodes produced by pruning were collapsed, and duplicate species tips were not permitted. The primary morphology and conditional-genome-size analyses used 24 tips. Tree sensitivity was evaluated by comparing the focal tree pruned to the same 18-species taxon set with a published 18-species topology, thereby separating taxon-set effects from differences between trees.

Positive trait values were log10-transformed and sample-standardized before comparative fitting. Pairwise associations among conditional genome size, nucleus area, and cell area were estimated with Pagel-lambda phylogenetic generalized least squares. Lambda was estimated by maximum likelihood over 10<sup>-7</sup> to 1. Coefficient tests used residual degrees of freedom and the *t* distribution, and 95% confidence intervals and generalized *R*<sup>2</sup> were reported. Forward, reverse, and two-predictor models were retained as sensitivity analyses rather than interpreted as causal directions.

Phylogenetic path analysis was conducted on conditional genome size, nucleus area, and cell area. The candidate set contained 25 directed acyclic graphs representing 11 Markov-equivalence classes, ten of which contained testable conditional-independence claims. Pagel-lambda models were primary, alpha was 0.05, and models within 2 CICc units of the best-supported class were considered competitive. Markov-equivalent graphs were not used to assign a unique causal direction.

Model sensitivity was examined with 250 measurement-bootstrap datasets that resampled paired cell-nucleus objects and the image and nucleus hierarchy underlying genome-size estimates. Additional analyses comprised 24 leave-one-species-out fits, 200 published bootstrap trees, and four evolutionary models: Brownian motion, Pagel lambda, Ornstein-Uhlenbeck with a fixed root, and Ornstein-Uhlenbeck with a random root. Simulation calibration used 800 datasets generated from four candidate causal classes under independent and Brownian residual regimes. Multiplication of IOD by a positive calibration constant cancels after log transformation and standardization, so fitted coefficients are invariant to the IOD-to-pg conversion. This mathematical invariance does not validate an absolute genome-size scale.

The fitted path model contains only conditional genome size, nucleus area, and cell area. TE composition and maximum snout-vent length have not been incorporated into this model, and the current three-trait analysis was not treated as a completed all-data causal model.

## Software and data availability

Downstream analyses were conducted from preserved RepeatModeler, dnaPipeTE, RepeatMasker, LTR mapping, and reviewed microscopy outputs. Local processing included species-identity joins, technical-run averaging, filtering, diversity calculations, CLR-PCA, phylogenetic analyses, and figure generation. The phylogenetic analyses used R 4.4.3, ape 5.8.1, phylolm 2.6.5, and phylopath 1.3.1. The released LTR and reviewed-morphology tables retain row-level eligibility, exclusion, and quality-control fields. Other raw and high-performance-computing outputs remain in frozen archival storage. Released files are accompanied by source paths, dimensions, column definitions, and SHA-256 checksums.
