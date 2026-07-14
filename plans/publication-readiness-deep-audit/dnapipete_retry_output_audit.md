# dnaPipeTE computational retry-output audit

## Provenance decision

Suffixes appended after an 11-character SRX accession (for example, `2`, `3`,
or `R2`) identify computational dnaPipeTE retries after failed or incomplete
runs. They do not identify new specimens, new SRA accessions, read mates, or
biological replicates.

The repository currently stages one suffixed output:

- `SRX19952890_reads_per_component_and_annotation`
- `SRX19952890R2_reads_per_component_and_annotation`

Both map to *D. orestes* / `SRX19952890`.

## File comparison

| Metric | Base output | `R2` retry |
|---|---:|---:|
| Component rows | 370,066 | 369,904 |
| Summed `#reads` | 7,404,782 | 7,404,964 |
| Summed aligned bases | 971,774,483 | 971,590,301 |
| Summed RepeatMasker hit length | 106,739,459 | 106,857,392 |

Additional comparison results:

- aligned-base profiles by native RepeatMasker classification have Pearson
  correlation `0.9999768`
- the largest classification-composition difference is `0.25685` percentage
  points
- the files share only `126,509` exact contig names (contig-name Jaccard
  `0.2062`) and `359` exact full rows

These are two full-scale stochastic computational realizations of the same
biological input, not non-overlapping shards that should be concatenated.

## Analysis consequence

The historical all-resource merge includes both outputs and therefore
double-weights `SRX19952890`. The final TE/genome analysis panel excludes
*D. orestes*, so neither output enters the corrected analysis-18 mass ledger.
The corrected writer also requires exactly one source label for every included
species and fails if a retry suffix would be summed.

## Release rule

For any future rerun, retain suffixes as `computational_run_id` provenance and
use a small selection manifest to nominate exactly one completed output per
SRA. Selection should be based on completion status, parameters, logs, and
checksums. Retry outputs must never be counted as independent biological
replicates or silently summed.
