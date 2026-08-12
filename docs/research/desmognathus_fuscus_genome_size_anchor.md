# *Desmognathus fuscus* genome-size calibration anchor

**Research conclusion (checked 2026-08-12):** no source located supports **16.36 pg** as a published *D. fuscus* C-value. The study now uses the modern pooled assembly estimate of **16.1 Gbp** as its explicitly assembly-derived reference. With **1 pg = 0.978 Gbp**, this gives **16.1 / 0.978 = 16.462167689 pg per 1C genome**. This is a conditional calibration anchor, not an independently measured C-value. The method-matched published Feulgen value of 15.13 pg remains a literature comparator rather than the selected anchor or specimen-specific evidence for standards 32469 and 32470.

## Published C-values

The [Animal Genome Size Database FAQ](https://genomesize.com/faq.php) defines genome size/C-value as haploid nuclear DNA content (1C). It contains five *D. fuscus* records, none equal to 16.36 pg:

| 1C value (pg) | Method and cells | Primary source |
|---:|---|---|
| 15.00 | Feulgen densitometry; red blood cells | Morescalchi and Olmo (1974), as indexed in [database record 554](https://genomesize.com/result_species.php?id=554) |
| 15.00 | Feulgen densitometry; red blood cells | Morescalchi and Olmo (1982), as indexed in [database record 557](https://genomesize.com/result_species.php?id=557) |
| **15.13** | **Feulgen densitometry; red blood cells; *Rana esculenta* standard** | [Olmo (1974), DOI 10.1080/11250007409430082](https://doi.org/10.1080/11250007409430082); [database record 556](https://genomesize.com/result_species.php?id=556) |
| 17.70 | Flow cytometry; red blood cells | [Licht and Lowcock (1991), DOI 10.1016/0305-0491(91)90089-V](https://doi.org/10.1016/0305-0491(91)90089-V); [database record 555](https://genomesize.com/result_species.php?id=555) |
| 18.00 | Feulgen densitometry; red blood cells | [Bachmann (1970), DOI 10.1007/BF00277456](https://doi.org/10.1007/BF00277456); [database record 553](https://genomesize.com/result_species.php?id=553) |

The database itself cautions that discrepant records can reflect experimental error and recommends judging values by method, standard, and matching cell type. On those criteria, 15.13 pg is the closest methodological match to this study's Feulgen–Schiff erythrocyte workflow. It is nevertheless an old species-level measurement made before modern recognition of the *D. fuscus* complex, so it does not establish the C-value of specimens 32469 or 32470.

## NCBI assemblies and Ensembl

The live [NCBI *D. fuscus* genome dataset](https://www.ncbi.nlm.nih.gov/datasets/genome/?taxon=52100) reports the following haploid assembly spans:

| Accession | Assembly | Span (bp) | Converted span (pg at 0.978 Gbp/pg) |
|---|---|---:|---:|
| [GCA_050004315.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_050004315.1/) | aDesFus1-2.1; pooled AMNH A-194068/69 | 16,117,580,045 | 16.480 |
| [GCA_049999565.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_049999565.1/) | Dfus2.0; AMNH A-194069 | 15,792,375,770 | 16.148 |
| [GCA_046863835.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_046863835.1/) | GW_Dfus_1.0; AMNH A-194068 | 15,877,946,558 | 16.235 |
| [GCA_030265095.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_030265095.1/) | ASM3026509v1; highly fragmented | 1,079,203,562 | 1.103 |

The conversion uses the database's cited relation **1 pg = 0.978 Gbp**. Assembly span is not a cytometric C-value: it can be altered by unresolved repeats, collapsed sequence, duplicated haplotypes, contamination, and incompleteness. The modern assembly paper reports a combined assembly of about 16.1 Gb and explicitly notes that it may be slightly larger than the true genome because repeats remain unresolved ([Myers et al. 2025](https://doi.org/10.1093/g3journal/jkaf157)). These assemblies also derive from different animals than standards 32469 and 32470.

A live search of the official [Ensembl species inventory](https://rest.ensembl.org/info/species?content-type=application/json) and [Ensembl Rapid Release species list](https://rapid.ensembl.org/species.html) found no *Desmognathus fuscus* entry. Ensembl therefore provides no independent C-value or assembly-based support for 16.36 pg.

## Selected anchor and interpretation of the old 16.36-pg value

Using the accepted conversion, `16.000 Gbp / 0.978 Gbp per pg = 16.360 pg`. The exact match makes it likely that 16.36 was obtained by converting a rounded 16.0-Gb genome/assembly estimate. This is an inference, not documented provenance, and should not be described as a published *D. fuscus* C-value.

The selected reference instead follows the published rounded pooled-assembly estimate:

`16.1 Gbp / 0.978 Gbp per pg = 16.462167689 pg per 1C genome`

The exact NCBI assembly span is retained above for provenance, but the calibration decision uses the paper's stated 16.1-Gbp value rather than silently substituting the exact assembly-span conversion.

## Defensible calibration and remaining validation

For erythrocytes, the measured nucleus is diploid (2C), whereas reported genome size is 1C. When target and reference are the same cell type and ploidy, the factor of two cancels:

`target genome size (1C pg) = (target IOD / reference IOD) × reference genome size (1C pg)`

The author confirmed that Process_413/specimen 32469 and Process_414/specimen 32470 were the *D. fuscus* DNA-content standards and that the standards were stained in the same experimental runs as the unknowns. Use **16.462167689 pg 1C**, derived from the selected 16.1-Gbp pooled assembly, only as the named conditional anchor. The historical 15.00–18.00-pg C-values remain useful sensitivity comparators. Do not present the resulting estimates as independently validated absolute C-values until all of the following are resolved:

- recover the exact standard-to-target staining and imaging run map for Process_413/specimen 32469 and Process_414/specimen 32470;
- calibrate each target only against the standard processed in the same run, rather than averaging standards across runs;
- investigate the approximately 41% difference between the two standard-slide median IOD values and exclude saturation, focus, segmentation, background, or preservation artifacts;
- confirm voucher identity and locality for 32469/32470 so their modern *D. fuscus* lineage can be assessed;
- report the chosen 1C anchor, its source, the calibration equation, nuclei retained/excluded, and a sensitivity interval in the final data release.
