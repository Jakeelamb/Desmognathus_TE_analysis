# Materials and Methods bullet scaffold

**Evidence checked:** 12 August 2026

- Author-facing factual outline; convert these bullets into the authors' own prose.
- `[VERIFIED]`: supported by the current repository.
- `[HISTORICAL COMMAND]`: recovered from version history, but runtime provenance may remain incomplete.
- `[AUTHOR CONFIRMED]`: supplied directly by the author during Methods reconciliation.
- `[PLANNED — NOT PERFORMED]`: intended future work that must not be written as a completed method or result.
- `[AUTHOR INPUT]`: must come from notebooks, collection records, or original command logs.
- `[BOUNDARY]`: interpretation that must remain explicit.
- Data Availability remains a separate manuscript section.

## 1. Study design and data sources

- `[VERIFIED]` Declare analysis-specific denominators rather than one universal species panel.
- `[VERIFIED]` TE34: 34 species with a vetted SRA resource, assembly accession, matching tree tip, and TE output; source file `data/identity/te34_panel.csv`.
- `[VERIFIED]` Path24: 24 species with finalized cell, nucleus, nuclear-IOD, and phylogenetic-path evidence; source file `data/identity/path24_panel.csv`.
- `[VERIFIED]` Overlap21: exact `TE34 × Path24` intersection used only when both evidence streams are required; source file `analyses/03_phylogenetic_path/data/te_relative_iod_overlap.csv`.
- `[VERIFIED]` Overlap21 species: *D. aeneus*, *D. amphileucus*, *D. anicetus*, *D. apalachicolae*, *D. auriculatus*, *D. bairdi*, *D. campi*, *D. fuscus*, *D. gvnigeusgwotli*, *D. intermedius*, *D. kanawha*, *D. marmoratus*, *D. mavrokoilius*, *D. monticola*, *D. ocoee*, *D. orestes*, *D. perlapsus*, *D. tilleyi*, *D. valtos*, *D. welteri*, and *D. wrighti*.
- `[VERIFIED]` *D. brimleyorum*, *D. folkertsi*, and *D. ochrophaeus* remain in Path24 but lack an active TE34 resource; no TE value was imputed.
- `[VERIFIED]` *D. aeneus*, *D. orestes*, and *D. wrighti* were manually reviewed and retained below an automatic sample target; reduced support remains explicit.
- `[BOUNDARY]` A species-level genomic–microscopy join does not establish that both measurements came from the same biological individual.
- `[VERIFIED]` SRA and assembly resources are joined through the accession-explicit active lookup; the expert-reidentified *D. fuscus* resource is `SRX20497025 / GCA_032353935.1`.
- `[VERIFIED]` dnaPipeTE retry labels are technical runs, not biological replicates; the two *D. orestes* runs were normalized separately and equal-weight averaged.
- `[VERIFIED]` Collaborator maximum SVL comes from `analyses/03_phylogenetic_path/data/collaborator_max_svl.xlsx`; total length and literature body-size proxies are not silently substituted.
- `[VERIFIED]` Developmental mode, aquaticity, microhabitat, and elevation were curated from manual literature extraction, NC Biodiversity accounts, and AmphiBIO with source-priority flags.
- `[VERIFIED]` Exact collaborator-SVL source aliases `orestes-ac` (52 mm) and `orestes-b` (50 mm) map to *D. orestes*; the registered maximum-of-maxima rule gives 52 mm. All Overlap21 species consequently have TE, morphology/IOD, and maximum-SVL data; the final integrated model remains unselected and unfitted.
- `[AUTHOR INPUT]` State whether panel rules were specified before outcome inspection.
- `[AUTHOR INPUT]` Supply individual counts, specimen matching, locality/date, collector, voucher/tissue identifier, sex/stage, permit, ethics approval, and Nagoya/access-compliance details.

## 2. Genomic processing and transposable-element characterization

- `[HISTORICAL COMMAND]` Mitochondrial filtering: BWA reference indexing; `bwa mem -t 8`; SAMtools conversion/sorting; `samtools view -f 12 -F 256` to retain read pairs unmapped to the mitochondrial reference.
- `[AUTHOR INPUT]` Confirm read-download procedure/date, library layout, adapter and quality trimming, contamination screening, mate handling, BAM-to-FASTQ conversion, and software versions.
- `[HISTORICAL COMMAND]` RepeatModeler: TEtools `1.87.sif`, Singularity 3.7.4, `BuildDatabase`, and `RepeatModeler -threads 32 -LTRStruct -genomeSampleSizeMax 10000000000`.
- `[HISTORICAL COMMAND]` CD-HIT-EST: `-aS 0.8 -c 0.8 -G 0 -g 1 -b 500 -T 64 -M 500000`.
- `[AUTHOR INPUT]` Confirm the checksum lineage from the CD-HIT product to `/nfs/home/jlamb/TE_libs/dedupe_telib.fasta`; the repository does not prove that identity.
- `[HISTORICAL COMMAND]` dnaPipeTE: Singularity 4.1.2; 16 CPUs; nuclear-filtered mate 1; configured 15-Gb genome-size denominator; 0.1× target coverage; `RM_t=0.15`; two internal samples; custom `dedupe_telib.fasta`.
- `[BOUNDARY]` The configured 15-Gb dnaPipeTE denominator is a sensitivity setting, not a measured genome size or validated absolute TE load.
- `[VERIFIED]` Stored dnaPipeTE components were parsed to class, order, and superfamily; aligned bases define total, retained, and unresolved repeat mass for S06-S07 dnaPipeTE quality control. Unresolved aligned-base mass was excluded from diversity calculations.
- `[VERIFIED]` RepeatMasker `.align` inputs describe dnaPipeTE `Trinity.fasta` repeat-contig assemblies, not the accession-linked GCA assemblies.
- `[AUTHOR INPUT]` Supply dnaPipeTE version/commit, container digest, custom-library checksum, random-seed/iteration details, Trinity settings, RepeatMasker version/search engine/options/library release, and complete command logs.

## 3. TE composition, diversity, repeat landscapes, and LTR proxy

- `[VERIFIED]` Diversity was calculated separately at order and superfamily level after each species' classified categories were reclosed to one. S08 is one 68-row, five-column table with 34 order-level and 34 superfamily-level rows. The order-level stratum is the primary paper-facing diversity analysis, and the superfamily-level stratum is a sensitivity/audit view.
- `[VERIFIED]` Only exact-positive bins entered the formulas. Natural-log Shannon entropy was `H′ = −Σpᵢ ln(pᵢ)` and the explicitly named Gini-Simpson index was `1 − Σpᵢ²`. These are the only paper-facing diversity indices. Observed richness is retained as an audit/support count, `S = count(pᵢ > 0)`, not as a diversity endpoint.
- `[BOUNDARY]` Unresolved aligned-base mass remains only in S06-S07 dnaPipeTE quality-control tables and does not enter S08 diversity. Rare classified features excluded from the shared-feature PCA remain in both order- and superfamily-level diversity calculations.
- `[VERIFIED]` TE34 superfamily PCA: retain features with positive abundance in all 34 species after same-species retry averaging; close rows; apply CLR; center columns; fit without post-CLR variance scaling; orient each axis so its largest absolute loading is positive. The retained matrix contains no zeros, so no zero replacement is required.
- `[VERIFIED]` S09 scores, S10 variance, S11 loadings, and Figure 2 now come from the same retry-averaged 24-feature CLR-PCA fit.
- `[VERIFIED]` Formal clustering audit used all 23 nonzero, unscaled PCA axes, which preserve the complete 24-feature CLR/Aitchison distances; PC1–PC2 and PC1–PC6 fits were sensitivity analyses rather than the selection surface.
- `[VERIFIED]` K-means diagnostics covered `k = 1…8` with 1,000 starts. The primary gap statistic used 1,000 scaled-PCA null datasets, squared Euclidean dispersion, the Tibshirani one-standard-error selector, and seed 20260810.
- `[VERIFIED]` No discrete cluster solution passed the rebuild gates. Gap selected `k=1`; maximum silhouette selected a weak `k=7` solution with two singleton clusters; the only singleton-free candidate was `k=2` (9/25 species, mean silhouette 0.179) and was unstable under 28-of-34 species subsampling (median adjusted Rand index -0.042).
- `[BOUNDARY]` The clustering gates were adopted during this post-hoc rebuild, not preregistered. Candidate `k=2` and `k=7` assignments are retained as rejected audit outputs and must not be described as supported biological clusters.
- `[VERIFIED]` Full 24-feature composition showed multivariate phylogenetic signal under Adams' `K_mult` with 99,999 complete-vector tip permutations: `K_mult = 0.6581`, `p = 0.00001`. Because `K_mult < 1`, signal was weaker than Brownian-motion expectation even though it exceeded the randomized-tip null.
- `[VERIFIED]` Post-hoc reproductive-strategy PERMANOVA used the full Aitchison distance matrix for 33 coded TE34 species and 99,999 label permutations: `R² = 0.2624`, `p = 0.00001`; the spatial-median, bias-adjusted PERMDISP diagnostic was not significant (`p = 0.6763`).
- `[VERIFIED]` Holm-adjusted post-hoc contrasts supported direct development versus aquatic larvae (`R² = 0.1217`, adjusted `p = 0.00045`) and aquatic versus terrestrial eggs among aquatic-larval species (`R² = 0.1691`, adjusted `p = 0.00002`).
- `[VERIFIED]` In a marginal size sensitivity on 33 complete species, reproductive strategy remained associated after log10 maximum SVL (`R² = 0.1991`, `p = 0.00001`), whereas SVL added no marginal association after strategy (`R² = 0.0205`, `p = 0.5717`).
- `[BOUNDARY]` Reproductive strategy was supplied after expert review of the initial PCA, so these are post-hoc species-level associations. Ordinary label permutations do not correct for shared phylogeny, and the source tree's citation/calibration provenance remains unresolved.
- `[BOUNDARY]` PCA is an exploratory ordination, not a causal predictor. Continuous phylogenetic or trait structure does not imply stable discrete clusters.
- `[VERIFIED]` Corrected RepeatMasker substrate: 9,533,364 hits and 2,135,413,681 inclusive aligned bp across 34 species.
- `[VERIFIED]` Hit length: `abs(query_end − query_start) + 1`; invalid divergence/length rows excluded; reported divergence floored into one-percentage-point bins with a terminal 50% bin.
- `[VERIFIED]` Each species was normalized by its total inclusive aligned RepeatMasker hit bp.
- `[VERIFIED]` Across-species landscape means zero-fill every species/category/bin combination before equal-weight averaging across TE34.
- `[BOUNDARY]` Landscape fractions describe aligned repeat annotations—not genome occupancy, insertion counts, direct insertion ages, interval-deduplicated coverage, or evolutionary rates.
- `[BOUNDARY]` `Unclassified` is an accounting/annotation category, not a biological TE order.
- `[VERIFIED]` LTR candidates: at least 3,000 bp; exact sequence-coordinate TEsorter joins; TEsorter `Order = LTR`, `Superfamily = Gypsy`, and at least five whitespace-delimited `DOMAIN|MODEL` annotations.
- `[HISTORICAL COMMAND]` Paired, trimmed, mitochondrial-filtered reads were mapped to the corresponding per-species LTR contigs with Bowtie2 `-L 20 --very-sensitive-local`; SAMtools was then used for BAM conversion, sorting, indexing, and per-position depth.
- `[AUTHOR INPUT]` Bowtie2 and SAMtools versions, MAPQ threshold, multimapper handling, secondary/supplementary-alignment handling, and duplicate policy remain unresolved; do not infer them from the stored depth products.
- `[VERIFIED]` Regional depth retains unreported and explicit zero-depth positions as zero.
- `[VERIFIED]` Historical LTR audit: an archived delimiter bug treated 1,088 rows as five/six-domain candidates; two truncated depth files were excluded, leaving 1,086 usable rows preserved in S14. Correct parsing yields 408 eligible LTR/Gypsy elements with at least five domains across 30 species.
- `[VERIFIED]` Primary LTR branch: 381 of 408 eligible elements are retained by the required within-species two-sided 1.5-IQR rule; 407 independently pass the 80% regional positive-depth gate; their intersection contains 380 primary elements across all 30 species.
- `[VERIFIED]` The 80% branch requires positive read depth at `≥ 0.8` of expected positions separately in the left LTR, right LTR, and internal region. It is not a reporting-completeness filter; every usable depth file reports all expected positions.
- `[VERIFIED]` Species summaries: median, geometric mean, trimmed/winsorized estimates, 2,000 bootstrap replicates, and leave-one-element-out influence.
- `[VERIFIED]` The original analysis notebook used an inclusive, within-species, two-sided Tukey filter (`Q1 − 1.5 × IQR ≤ ratio ≤ Q3 + 1.5 × IQR`). In the corrected zero-aware release, fences are recalculated within species over the 408 eligible elements before intersection with the independently declared 80% gate. S14 retains all 1,086 usable historical rows and marks the 678 ineligible rows explicitly; filtering changes primary branch membership rather than deleting audit evidence.
- `[BOUNDARY]` A later July 2026 compound one-element screen is not the original IQR method and is retained only as a retired provenance reconstruction.
- `[BOUNDARY]` Terminal:internal depth is a mapping/deletion-footprint proxy—not an ectopic-recombination, solo-LTR, or DNA-loss rate.
- `[VERIFIED]` Integrated overlap21 sensitivity: full classified-superfamily Shannon entropy joined to relative nuclear IOD; historical salamander points are excluded from the publication-facing table.

## 4. Specimen preparation, microscopy, and image-based phenotyping

- `[AUTHOR CONFIRMED]` Each `specimen_id` represents one individual animal; images and measured objects sharing a `specimen_id` are repeated observations from that animal, not independent biological replicates.
- `[VERIFIED]` The active releases contain 42 animals in the cell–nucleus morphology stream and 50 in the nuclear-IOD stream: 41 occur in both streams, for 51 distinct animals across their union.
- `[AUTHOR CONFIRMED]` All animals used for erythrocyte collection were fixed in ethanol before collection; there was no live-specimen collection group.
- `[AUTHOR CONFIRMED]` Erythrocytes were collected from a heart incision using a pipette.
- `[AUTHOR CONFIRMED]` The collected erythrocytes were pipetted onto the center of the slide.
- `[AUTHOR INPUT]` Ethanol concentration and fixation/storage duration; anticoagulant/concentration; volume; time to smear; storage temperature.
- `[AUTHOR CONFIRMED]` Slides were stained with Feulgen–Schiff reagent following the protocol in Hardie, Gregory, and Hebert (2002), *From Pixels to Picograms: A Beginners' Guide to Genome Quantification by Feulgen Image Analysis Densitometry* (doi: `10.1177/002215540205000601`).
- `[AUTHOR CONFIRMED]` Slides were imaged with an Olympus APX100 using a 40× objective.
- `[AUTHOR CONFIRMED]` Cellpose was used for cell segmentation and YOLO was used for nucleus segmentation.
- `[VERIFIED COMPUTATIONAL]` The preserved masks retain source-image, tile, mask, centroid, overlap, and object-level traceability.
- `[AUTHOR INPUT]` Confirm Cellpose/YOLO versions, exact model weights/checksums, training lineage, inference thresholds, and which checkpoint generated each archived mask set.
- `[VERIFIED]` Eligibility required a valid physical one-to-one cell–nucleus pair, accepted/core flags, no edge touch, positive area/IOD fields, high overlap support, and no suspect/discard flags.
- `[VERIFIED]` Eligible cells were ranked within species by descending cell area.
- `[AUTHOR CONFIRMED]` Cell and nucleus masks were manually reviewed, and rejected candidates were replaced by the next-largest eligible accepted cell.
- `[AUTHOR CONFIRMED]` The sampling target was the 50 largest accepted cells per species.
- `[BOUNDARY]` Cell selection estimates the upper tail of detected eligible erythrocytes, not the population-average cell.
- `[AUTHOR CONFIRMED]` Accepted masks were used to calculate cell area, nucleus area, and nuclear IOD.
- `[VERIFIED]` Final morphology panel: 1,152 reviewed cell–nucleus pairs across Path24; 21 species contributed 50, *D. aeneus* 44, *D. orestes* 50, and *D. wrighti* 8.
- `[VERIFIED]` Cell and nucleus areas derive from mask pixels using the stored `0.12 µm/pixel` linear scale (`0.0144 µm²` per pixel); species summaries use medians and conditional paired-object bootstrap intervals.
- `[VERIFIED]` Morphology contains 42 animals. Seven species are represented by one animal (*D. anicetus*, *D. apalachicolae*, *D. auriculatus*, *D. bairdi*, *D. gvnigeusgwotli*, *D. intermedius*, and *D. perlapsus*); S19 records animals per species and the largest single-animal contribution.
- `[BOUNDARY]` Morphology intervals resample accepted objects conditional on the observed images and animals. They do not estimate population-level among-animal uncertainty, and high contribution from one animal remains visible rather than being treated as independent replication.
- `[VERIFIED]` Pixel optical density: `log10(I_background / intensity)`; nuclear IOD: sum across nuclear-mask pixels.
- `[VERIFIED]` IOD matching used log edge sharpness and log relative ring noise; IOD, darkness, nucleus area, and cell area were excluded from matching.
- `[VERIFIED]` Final IOD panel: 805 reviewed nuclei across 24 species, 51 images, and 50 specimens.
- `[VERIFIED]` Species IOD estimate: arithmetic mean of image-specific median per-nucleus IOD values; hierarchical bootstrap resampled images and then nuclei within images; 2,000 replicates; seed 20260710. Dividing all 24 estimates and interval bounds by the Path24 median anchor (`950.176112815`) gives the primary relative-IOD index.
- `[VERIFIED]` S23 (`nuclear_iod_by_image.csv`) is the unnormalized 51-image aggregation audit. Its image medians and means are not relative indices or ratios.
- `[VERIFIED]` The declared IOD quality-balance gate is maximum absolute standardized mean difference `≤ 0.10` and maximum Kolmogorov–Smirnov distance `≤ 0.25`. *D. aeneus*, *D. ochrophaeus*, *D. orestes*, and *D. wrighti* fail this gate and remain visibly flagged rather than silently removed.
- `[VERIFIED]` *D. anicetus*, *D. bairdi*, and *D. gvnigeusgwotli* each have one IOD image from one animal. Point shape exposes this support limitation in Figure 6.
- `[AUTHOR CONFIRMED]` Process_413/specimen 32469 and Process_414/specimen 32470 were the *D. fuscus* DNA-content standards, and the standards were stained in the same experimental runs as the unknowns. The exact standard-to-target run map and the approximately 41% difference between the two standard-slide median IOD values remain unresolved calibration/QC issues.
- `[BOUNDARY]` Primary variable is relative image-derived nuclear IOD. The separate descriptive scale uses the published pooled *D. fuscus* assembly estimate of `16.1 Gbp`, converted with `1 pg = 0.978 Gbp` to `16.462167689 pg/1C`. This assembly-derived anchor remains conditional and is not an independently validated C-value.

## 5. Phylogenetic comparative analyses

- `[AUTHOR CONFIRMED]` The focal time-calibrated *Desmognathus* phylogeny was provided directly by Alex Pyron.
- `[VERIFIED]` Primary phenotype/path analysis uses all 24 finalized species and an exact 24-tip trim of the focal dated tree.
- `[VERIFIED]` Tree-tip labels were normalized, tips without complete analysis data were removed with `ape::drop.tip`, and resulting single-descendant internal nodes were collapsed; duplicate species tips were not permitted.
- `[VERIFIED]` The published-tree sensitivity uses 18 species because six Path24 taxa are absent; it does not redefine the primary panel.
- `[AUTHOR INPUT]` Supply the focal-tree publication/archive identifier, calibration method, and branch-length provenance; the collaborator source is now confirmed but these citation details remain unresolved.
- `[BOUNDARY]` Do not restore or use the archived `results/phylogeny/processed_phylogeny.nwk` for inference; its exact 0.8 branch scaling is undocumented. The active exact tree is `analyses/03_phylogenetic_path/trees/path24_time_tree.nwk`.
- `[VERIFIED]` Positive traits were log10-transformed and sample-standardized before comparative fitting.
- `[VERIFIED BOUNDARY]` Express the path node as relative nuclear IOD. Its standardized log10 values are numerically identical to those from the conditional *D. fuscus*-scaled pg column because multiplication by one positive constant cancels during standardization; this does not validate an absolute genome size.
- `[VERIFIED]` Pagel λ was fitted by maximum likelihood over `10^-7…1`; coefficient uncertainty used residual degrees of freedom, t tests, 95% intervals, and generalized R².
- `[VERIFIED]` Path design: 25 labeled DAGs; 11 Markov-equivalence classes; ten testable classes; λ primary model; α=0.05; competitive threshold `ΔCICc ≤ 2`.
- `[VERIFIED]` Sensitivity suite already preserved: 250 measurement bootstraps; 24 leave-one-species-out fits; 200 published bootstrap trees; four evolutionary models; 800 simulations.
- `[VERIFIED]` Saved path runtime: phylopath 1.3.1; phylolm 2.6.5; ape 5.8.1; R 4.4.3.
- `[BOUNDARY]` Path results are exploratory association-model sensitivities; Markov-equivalent arrows cannot be causally oriented from these data.
- `[BOUNDARY]` Analyses using both TE and phenotype evidence must use overlap21; Path24-only morphology/IOD models remain 24 species.
- `[VERIFIED]` The currently implemented Path24 model contains only relative nuclear IOD, nucleus area, and cell area; TE composition and SVL are not fitted predictors.
- `[PLANNED — NOT PERFORMED]` The final analysis step may integrate TE composition and collaborator maximum SVL with the morphology/IOD traits on the exact complete-case overlap and a correspondingly pruned tree. Its variables, DAGs, analysis panel, and statistical specification have not yet been selected or fitted.
- `[BOUNDARY]` Do not describe the current three-trait Path24 results as the completed all-data phylogenetic path analysis.

## 6. Software and reproducibility

- `[VERIFIED]` HPC-derived dnaPipeTE, RepeatModeler, RepeatMasker, LTR discovery/mapping/depth, segmentation, and model-training products were treated as immutable inputs.
- `[VERIFIED]` Local work was limited to panel joins, retry averaging, summaries, CLR-PCA, statistics, packaging, checksums, and figure rendering.
- `[VERIFIED]` Release artifacts carry SHA-256 hashes, row/column counts, source paths, and manuscript-use classifications.
- `[VERIFIED]` Publication figure system: Nimbus Sans; at least 7-pt target text; accessible colors plus shape/line-type redundancy; embedded-font CMYK PDF; review PNG; 300-ppi CMYK TIFF.
- `[AUTHOR INPUT]` Supply unresolved upstream software versions, original command logs, focal-tree citation, laboratory protocols, and final repository DOI.
- `[AUTHOR ACTION]` Archive code and underlying data in a persistent DOI-issuing repository; provide accessions/DOIs and dataset citations in the separate Data Availability section.
- `[AUTHOR ACTION]` GBE currently requires disclosure of any AI-assisted coding, data processing, proofreading, or language work in the cover letter and Methods and/or Acknowledgements; authors must write and approve the final disclosure.

## Analyses not to rerun

- SRA preprocessing or mitochondrial filtering.
- RepeatModeler, TE-library construction, or CD-HIT.
- dnaPipeTE or RepeatMasker.
- DANTE, LTRharvest, LTRdigest, or TEsorter.
- Bowtie2/SAMtools LTR mapping or depth generation.
- Cell/nucleus model training or segmentation inference.
- Existing full path-analysis bootstrap, tree-sensitivity, leave-one-out, or simulation suite unless a separate rerun is explicitly authorized.
