# D. fuscus calibration anchor: primary-source review

Reviewed 2026-10-13 for checkout `15dfb3c18755cc8659b626fa213087c62c6db02f`. This is an evidence review, not a new measurement. Only this file was authored for this task.

**Recommendation: retain 16.1/0.978 only as an explicitly conditional assembly anchor; it does not establish the absolute 1C DNA content of standards 32469 or 32470.** High confidence in that distinction. The current evidence does not justify replacing it with a historical species value and calling the standards independently calibrated. No qualified specimen-specific C-value was located in this bounded search; that is not proof that none exists.

Session facts supplied by the user/coordinator: standards are specimens 32469/32470, Processes 413/414; all retained unknowns and both standards were stained in one Feulgen run; retained unknown images are 40× brightfield. Historical fixation and specimen identity remain uncertain. These are session facts, not facts established by the papers below. Raw-image behavior and the implemented denominator are reviewed by other investigators.

## What the assembly establishes

The final [Myers et al. 2025 article](https://academic.oup.com/g3journal/article/15/9/jkaf157/8188250) was published 7 July and corrected/typeset 22 July 2025. Its genomic animals were AMNH A-194068–9, collected 20 May 2023 at Bell Branch, Wilkes County, North Carolina, assigned to fuscus B3 (Sampling and sequencing; Fig. 1). Hifiasm used an expected 14-Gb size input. Pooling used four-haplotype settings; the initial 24-Gb assembly contained duplicated haplotigs. Two purge_dups rounds produced the 16.1-Gb primary assembly (Assembly and scaffolding; Results and discussion). The authors allow a true size near 14–15 Gb because repeats remain unresolved. Their Methods supply no independent per-nucleus DNA-content measurement. Extracted DNA yield is not such a measurement.

Precise machine-readable passages were read through the official [NCBI BioC copy](https://www.ncbi.nlm.nih.gov/research/bionlp/RESTful/pmcoa.cgi/BioC_xml/PMC12405886/unicode): offsets 6192 (sampling), 12528/13166 (assembly assumptions/pooling), and 19166/20239 (purging and true-size caveat). Publisher PDF indexing places Assembly and scaffolding on p. 3 and the pertinent Results and discussion on p. 5. No independent cytometric C-value table was identified.

The current [NCBI assembly record](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_050004315.1/) and its [official dataset report](https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/GCA_050004315.1/dataset_report), retrieved 2026-10-13, identify:

| Field | Retained evidence |
|---|---|
| Accession/version | GCA_050004315.1, current; aDesFus1-2.1; released 2025-04-29 |
| Assembly type | Haploid |
| Source | Pooled AMNH A-194068–9; SAMN40589320 derived from SAMN38201181 and SAMN38201182 |
| Total sequence length | 16,117,580,045 bp |
| Ungapped length | 16,117,153,145 bp |

Inference: “haploid” describes the intended assembly representation; it does not certify a DNA-content standard. The accession supports a traceable assembled-length estimate from other animals, with the limitations documented above. It does not transfer a measured C-value to 32469/32470.

## Conversion and the remaining scientific distinction

The original [Doležel et al. 2003 letter](https://onlinelibrary.wiley.com/doi/pdf/10.1002/cyto.a.10013), p. 127, right column, gives `DNA content (pg) = genome size (bp)/(0.978 × 10^9)`. Table 1 and the preceding derivation concern nucleotide-pair mass, assuming equal AT/GC contributions and neglecting modified nucleotides. This is a units conversion, not a test of assembly completeness or a cytometric measurement of D. fuscus.

Calculated here: `16.1/0.978 = 16.46216768916155 pg`. The stored 16.462167689 value follows the paper’s rounded span; converting the exact accession span gives approximately 16.48014 pg. Extra arithmetic digits do not increase the precision of the rounded input or establish calibrant accuracy. Keep these two provenance choices distinct.

A conditional target estimate remains `(target IOD/reference IOD) × assumed reference 1C value`. For equal 2C erythrocyte populations, the common factor of two cancels. Successful comparator measurement and a trustworthy assumed mass are separate requirements; improvement of one does not establish the other. The literature interval is a set of alternative assumptions, not a confidence interval for these two specimens.

## Hardie 2002: what same-run staining resolves

Original [Hardie et al. 2002 PDF](https://journals.sagepub.com/doi/pdf/10.1177/002215540205000601), pp. 743–746 and Appendix 2 (p. 747):

| Locator | Relevant guidance |
|---|---|
| p. 743 | Separate standard slides in each staining run are acceptable; match cells and use reliable C-values. |
| pp. 743–744, Fig. 3; equations 4–5 | IOD ratios use reference C-values; standard curves check linearity. |
| pp. 744–745, Figs. 4–5 | Slide age/fixation affect stain and hydrolysis. Prolonged MFA fixation addresses some age differences. |
| pp. 745–746, Fig. 6 | Compaction/cell type affect uptake; hydrolysis differences can remain within cell types. |
| p. 741 | 40× is discussed; fewer pixels can reduce accuracy. 100× is generally favored except for large nuclei. |
| p. 746 | Independent confirmation can include complete sequencing. |

Case inference: single-run confirmation resolves run membership; fixation history and standard mass remain open. Two same-species slides do not establish a DNA-range curve. Within-preparation CV guidance supplies no universal acceptable difference between reference-slide medians. No quantitative gate is proposed.

## Independently measured historical values: provenance limits

Historical D. fuscus measurements are indexed, but **no original record inspected here establishes a suitable measured C-value for 32469/32470**. The following are finding-aid records, not newly verified original tables:

| Indexed 1C value | Indexed method/standard | Original source and access outcome |
|---|---|---|
| 15.13 pg | Feulgen erythrocytes; Rana esculenta 6.00 pg ([record 556](https://genomesize.com/result_species.php?id=556)) | [Olmo 1974](https://doi.org/10.1080/11250007409430082), pp. 29–33. Publisher [issue listing](https://www.tandfonline.com/toc/tizo19/41/1) identifies the article and later online posting (14 September 2009); original PDF retrieval returned 403. No table/voucher/locality/fixation verification. |
| 18.00 pg | Feulgen erythrocytes; Rana pipiens 6.70 pg/various ([record 553](https://genomesize.com/result_species.php?id=553)) | [Bachmann 1970](https://doi.org/10.1007/BF00277456), pp. 289–293; [PubMed metadata](https://pubmed.ncbi.nlm.nih.gov/5433173/) confirmed. Full text unavailable through attempted routes. |
| 17.70 pg | Flow cytometry erythrocytes; Ambystoma laterale 29.20 pg ([record 555](https://genomesize.com/result_species.php?id=555)) | [Licht and Lowcock 1991](https://doi.org/10.1016/0305-0491(91)90089-V), pp. 83–92. Publisher full text returned 403; individual provenance unverified. |
| 15.00 pg, two records | Feulgen erythrocytes; specified frog/human standards in [record 554](https://genomesize.com/result_species.php?id=554), unspecified standard in [record 557](https://genomesize.com/result_species.php?id=557) | Indexed to Morescalchi and Olmo 1974/1982. Original per-specimen evidence not recovered; do not count the two records as demonstrated independent biological calibrations. |

Confidence is moderate that historical measured values exist as described by the index, but low in their specimen-level applicability here. In particular, 15.13 pg is a method-matched candidate, not a demonstrated uniquely superior value: several records use Feulgen erythrocytes. Choosing among them requires original standard, preparation, population, and taxonomic evidence.

[Hally et al. 1986](https://link.springer.com/article/10.1007/BF00494802), pp. 185–192, was checked at its original publisher Summary: sampled D. quadramaculatus, D. monticola, D. ochrophaeus, and D. wrighti, using blood-smear Feulgen cytophotometry with Xenopus nuclei. **D. fuscus was not among those sampled species**, so this paper cannot supply an independently measured fuscus reference value. A broader search also found Sharbel et al. 1995 ([author-posted original scan](https://www.researchgate.net/profile/David-Green-40/publication/270465955_Partial_Genetic_Compatibility_and_Unidirectional_Hybridization_in_Syntopic_Populations_of_the_Salamanders_Desmognathus_fuscus_and_D_ochrophaeus/links/5617be1108ae3eadf96af8c1/Partial-Genetic-Compatibility-and-Unidirectional-Hybridization-in-Syntopic-Populations-of-the-Salamanders-Desmognathus-fuscus-and-D-ochrophaeus.pdf); [DOI](https://doi.org/10.2307/1446911)); the scanned text was not readable through available retrieval, and no C-value is adopted from it.

The apparent numeric conflict is therefore unresolved evidence from different measurement methods, reference scales, populations, historical species concepts, and assembly processing—not a controlled comparison showing which value belongs to these slides. The scientifically supportable output remains conditional scaling plus explicit uncertainty. High confidence that neither an assembly accession nor shared staining alone closes the absolute-calibration question; no claim of a validated absolute C-value is warranted by the evidence recovered here.
