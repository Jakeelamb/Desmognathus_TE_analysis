# Feulgen/IOD genome-size validation and QC gate

## Decision

Use two deliberately separate measurement sets.

1. **Cell and nucleus morphology:** the visually vetted 50 largest linked cell--nucleus pairs per species estimate an explicitly upper-tail morphology endpoint. They are not an estimate of the typical cell, typical nucleus, or genome-size distribution.
2. **Genome size:** estimate each species from a larger, predeclared population of QC-passing, isolated, in-focus, non-clipped **2C/G1 nuclei**, with reference cells co-processed in the same staining run. Do not choose the largest cells, largest nuclei, or highest-IOD nuclei for this analysis.

The closest taxon-specific precedent is Feulgen cytophotometry of *Desmognathus* blood smears calibrated with *Xenopus laevis* erythrocytes as an internal reference ([Hally et al. 1986](https://doi.org/10.1007/BF00494802)). It supports the use of an internal biological reference, but it does not validate a current calibration without current co-processing and QC.

## What turns IOD into an absolute genome-size estimate

For a nucleus, IOD is the sum of per-pixel optical density relative to a clear-field incident-light measurement. It becomes pg DNA only through a validated, DNA-stoichiometric Feulgen assay and reference calibration; intensity from an arbitrary image is not genome size ([Hardie, Gregory & Hebert 2002](https://doi.org/10.1177/002215540205000601)).

For each staining run, include a reference series with independently trusted DNA contents, ideally the same cell type and biological state as the unknowns, spanning below and above the expected *Desmognathus* range. The standards and unknowns must share the staining run; a standard curve is both the conversion and the run-level check. Hardie et al. recommend a broad multi-standard series, note that cell-type mismatch can cause error, and show why a standard curve must test staining across the range rather than merely rescale one mean ([Hardie et al. 2002](https://doi.org/10.1177/002215540205000601)).

The reported conversion must be retained per run:

\[
\widehat{G}_{u,r}=f_r(\mathrm{IOD}_{u,r}),
\]

where \(f_r\) is the archived standard-curve fit for run \(r\), not a single global scale reused without evidence. Record each standard's assigned C-value source, number of accepted nuclei, modal/mean IOD, curve coefficients, residuals, and unknown-to-standard bracketing status.

## Required QC evidence

| Level | Acceptance criterion | Why it is required |
|---|---|---|
| Biological state | Predefine the diploid target population (normally mature erythrocyte 2C/G1 nuclei); inspect an IOD histogram per specimen/slide and exclude S-phase, 4C/G2/M, aggregates, and unidentified mixtures from the 2C estimate. When a 4C population is present, use it as a stoichiometry check, not as observations to pool with 2C. | A cross-laboratory Feulgen image-cytometry validation proposed peak CV <6% and 4C/2C = 1.9--2.1 as slide QC ([Vilhar et al. 2001](https://doi.org/10.1006/anbo.2001.1394)). The 4C/2C test is applicable only where both populations are biologically expected and resolvable; it must be validated for this sample type, not assumed. |
| Reference calibration | Co-fix/stain unknowns and standards together; use a bracketed multi-standard curve for every run. A single external historical value is not a substitute. | Feulgen staining depends on fixation and hydrolysis conditions; co-fixation of sample and standard is specifically recommended ([Greilhuber & Temsch 2001](https://hrcak.srce.hr/160855)). Internal-standard image analysis has shown close agreement with an orthogonal method when linearity was demonstrated ([Voglmayr & Greilhuber 1998](https://doi.org/10.1006/fgbi.1998.1097)). |
| Camera and illumination | Before biological measurement, verify optical-density linearity with neutral-density filters/density wedges, spatial uniformity (shading), and drift. Lock exposure/gain/white balance; retain calibration images and test results for each session. | A calibration protocol supplies practical targets for OD linearity (slope 0.9--1.1, \(R^2>0.990\), SE <5%), field shading (slope CV <2%), and repeated-nucleus field-position IOD (CV <3%) ([Vilhar & Dermastia 2002](https://hrcak.srce.hr/file/5821)). A long-running clinical QC program likewise checks calibration, spatial uniformity, and camera linearity and documents focus, CCD, lamp, dust, and gain/offset failures ([Chiu et al. 2004](https://doi.org/10.1155/2004/794021)). |
| Image/field | Accept only a homogeneous, blank incident-light area; no clipped/saturated pixels in a measured nucleus; stable focus; and comparable field occupancy. Flag every excluded image with its failure reason. | Focus, glare, diffraction, sampling and computation all create Feulgen-OD error ([Smeulders & ten Kate 1987](https://doi.org/10.1364/AO.26.003249)). In vertebrate Feulgen image analysis, Hardie et al. use a non-clipped green-channel acquisition guideline (~190 maximum on an 8-bit scale) and found field-density changes could move mean IOD by >10% ([Hardie et al. 2002](https://doi.org/10.1177/002215540205000601)). Treat 190 as a method-specific starting point, not a universal threshold: the present camera must pass the density-wedge test. |
| Mask/object | Include only one intact nucleus per mask: entirely inside the field, non-overlapping, non-touching/merged, free of debris and obvious segmentation error. Preserve object-level thumbnails, mask, image ID, reviewer, and exclusion reason. | Manual exclusion of damaged, overlapping, edge-touching, anomalous, and non-nuclear objects is part of the animal Feulgen image-analysis method ([Hardie et al. 2002](https://doi.org/10.1177/002215540205000601)). Glare and nuclear morphology can generate spurious DNA differences if not controlled ([Bedi & Goldstein 1976](https://doi.org/10.1083/jcb.71.1.68)). |
| Run/batch | Require standards to pass before accepting an unknown run. Keep run, slide, field/image, specimen, and operator as explicit variables; never silently pool raw IOD across runs. | A 13-laboratory intercomparison found large variation in nominally diploid IOD despite calibration, motivating internal controls and standardization ([Thunnissen et al. 1997](https://pubmed.ncbi.nlm.nih.gov/9000581/)). |

## Minimum validation experiment before reporting pg/C-values

1. **Instrument qualification.** Archive density-wedge linearity, flat-field/shading, drift, focus, and clipping checks before acquisition and at a prespecified cadence. Fail a session instead of correcting an unknown sample after the fact.
2. **Stain-run qualification.** Co-process standards with every unknown slide/run. Show the standard curve, residuals, and 2C peak CV for every accepted run. Repeat/withhold any run outside its prespecified acceptance band.
3. **Biological-state qualification.** Demonstrate the target 2C population with per-slide IOD histograms. If a 4C population exists, show its 4C/2C ratio; if it does not, state why that control is unavailable and add a suitable same-run control.
4. **Independent accuracy check.** For a prespecified set of species spanning the observed IOD range, compare Feulgen estimates with an independent method (preferably internal-standard flow cytometry) or independently trusted co-processed reference material. Assess agreement with paired bias, limits of agreement, and proportional-bias/residual plots—not correlation alone. The cross-lab validation by Vilhar et al. found image cytometry comparable to photometric cytometry and flow cytometry across a 100-fold range ([Vilhar et al. 2001](https://doi.org/10.1006/anbo.2001.1394)); that is the standard of evidence to reproduce locally.
5. **Replicate the biological unit.** Estimate uncertainty by resampling at the highest available independent level (specimen, then slide/run, then image), retaining nucleus-level variation as a lower-level component. A tight bootstrap over 50 nuclei from one image/specimen is not an independent biological confidence interval.

## Reporting contract

For every species and run, publish a machine-readable QC ledger with: specimen, slide, staining run, reference identity/value, image/field, raw and accepted nuclear counts, 2C-gate rule, IOD histogram/peak CV, 4C/2C result or `not_applicable`, calibration fit/residuals, clipping/focus/field-uniformity status, mask-review decision, exclusion reason, and the exact aggregation/uncertainty model.

Report at least four uncertainty components separately: between-specimen, between-run/slide, between-image/field, and between-nucleus. A species with one specimen or one accepted image can have a descriptive calibrated value, but it must be labeled as having no estimable biological replication; it should not receive the same inferential status as a replicated species.

Until the co-processed-standard, linearity, biological-state, and independent-accuracy gates pass, label the trait **relative or reference-scaled nuclear IOD**, not an absolute pg DNA/C-value genome-size estimate.

## Sources

- Hally MK, Rasch EM, Mainwaring HR, Bruce RC. 1986. Cytophotometric evidence of variation in genome size of desmognathine salamanders. *Histochemistry* 85:185--192. [https://doi.org/10.1007/BF00494802](https://doi.org/10.1007/BF00494802)
- Vilhar B, Greilhuber J, Dolenc Koce J, Temsch EM, Dermastia M. 2001. Plant genome size measurement with DNA image cytometry. *Annals of Botany* 87:719--728. [https://doi.org/10.1006/anbo.2001.1394](https://doi.org/10.1006/anbo.2001.1394)
- Hardie DC, Gregory TR, Hebert PDN. 2002. From pixels to picograms: a beginners' guide to genome quantification by Feulgen image analysis densitometry. *Journal of Histochemistry & Cytochemistry* 50:735--749. [https://doi.org/10.1177/002215540205000601](https://doi.org/10.1177/002215540205000601)
- Vilhar B, Dermastia M. 2002. Standardisation in DNA image cytometry. *Acta Stereologica* 21:1--15. [full text](https://hrcak.srce.hr/file/5821)
- Greilhuber J, Temsch EM. 2001. Feulgen densitometry: some observations relevant to best practice in genome size research. *Acta Botanica Croatica* 60:285--298. [journal record](https://hrcak.srce.hr/160855)
- Voglmayr H, Greilhuber J. 1998. Genome size determination in the agarics *Agaricus bisporus* and *A. blazei* by flow cytometry and Feulgen image analysis. *Fungal Genetics and Biology* 25:108--115. [https://doi.org/10.1006/fgbi.1998.1097](https://doi.org/10.1006/fgbi.1998.1097)
- Smeulders AWM, ten Kate TK. 1987. A study on the information loss in DNA image cytometry. *Applied Optics* 26:3249--3254. [https://doi.org/10.1364/AO.26.003249](https://doi.org/10.1364/AO.26.003249)
- Bedi KS, Goldstein DJ. 1976. Apparent measurement of DNA content by Feulgen microdensitometry: a case of changing a constant. *Journal of Cell Biology* 71:68--77. [https://doi.org/10.1083/jcb.71.1.68](https://doi.org/10.1083/jcb.71.1.68)
- Chiu KY, Bulten W, et al. 2004. Quality assurance in DNA image cytometry. *Analytical Cellular Pathology* 26:55--70. [https://doi.org/10.1155/2004/794021](https://doi.org/10.1155/2004/794021)
- Thunnissen FBJM, den Bakker MA, et al. 1997. A European interlaboratory comparison of DNA image cytometry. *Cytometry* 28:123--132. [PubMed](https://pubmed.ncbi.nlm.nih.gov/9000581/)
