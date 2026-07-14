# dnaPipeTE configured-denominator absolute-load sensitivity

Git history recovers the deleted upstream run script at commit `cc64ffb134968250671ff0e3456255063a1aba97` (blob `738c7cece113ce03b34f7523bf2139431646fac1`). It configured a 15,000,000,000-bp genome size and 0.1× coverage, implying an independent 1,500,000,000-bp quantification sample for every run. dnaPipeTE documents `aligned_bases` as bases from that independent sample mapped back to annotated repeat contigs.

Dividing conserved repeat-aligned bases by the configured quantification-sample size gives:

- 18 species, exactly the declared final TE/genome panel
- repeat-aligned fraction range: 60.619%–68.196%
- median repeat-aligned fraction: 65.744%
- *D. fuscus* (`SRX20497025` / `GCA_032353935.1`): 66.010% repeat-aligned and 62.536% assigned to a DNA- or retrotransposon class

This is a **sensitivity estimate**, not yet a confirmatory predictor. The command was recovered from Git history rather than runtime logs; the container version/digest and custom-library checksum are absent; all species used the same configured 15-Gb genome size; only the nuclear-filtered first read mate was supplied; and replicate quantification samples or sampling uncertainty were not retained. The image-IOD genome-size measurements are from independent microscopy specimens and are not substituted into this denominator.

Primary method sources: [Goubert et al. dnaPipeTE paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC4419797/) and [author tutorial](https://tehub.org/tutorials/docs/dnaPipeTE).
