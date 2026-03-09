# LTR Age Audit

This note records the current status of the local LTR insertion-age branch.

## Bottom Line

- The previous local approach was not suitable for paper use because it tried to infer insertion age from RepeatMasker divergence-to-consensus hits rather than divergence between true paired 5' and 3' LTR sequences.
- The canonical local paired-element source is the ectopic master table produced by `scripts/processing/ec.py`, not the RepeatMasker adjacency heuristic.
- Local genome FASTA assemblies are now available, so this branch computes true 5' and 3' LTR sequence divergence directly from assembly coordinates.
- The default local output is sequence divergence, not absolute age in years, because no substitution rate is imposed automatically.

## Current Local Readiness

- Valid paired LTR elements in canonical ectopic master table: `1371`
- Species with at least one paired LTR element: `31`
- Pairs lacking a lookup-table genome accession: `80` across `1` species bucket(s)
- Complete (`Complete = yes`) LTR elements: `247`
- LTR pairs with both LTRs >= `100` bp: `1371`
- LTR pairs with five or more annotated domains: `1158`
- Recommended high-confidence LTR pairs (`Complete = yes` and `5+` domains): `247`
- Dominant superfamily among paired LTR elements: `Gypsy` (`1363` elements)
- Species currently ready for local sequence extraction with genome FASTA present: `30`

## Sequence-Based Divergence

- LTR pairs attempted for sequence estimation: `1371`
- LTR pairs with successful sequence-based divergence estimates: `1291`
- Species with at least one successful sequence estimate: `30`
- High-confidence pairs with successful sequence estimates: `238`
- Median ungapped comparable sites per successful alignment: `363`
- Median p-distance across successful alignments: `0.06701`
- Median K2P distance across successful alignments: `0.070564`

## Highest-Coverage Species

- `D.tilleyi`: `142` successful sequence estimates, median K2P `0.0674`, median comparable sites `340`
- `D.abditus`: `98` successful sequence estimates, median K2P `0.0702`, median comparable sites `372`
- `D.welteri`: `94` successful sequence estimates, median K2P `0.0721`, median comparable sites `400`
- `D.valentinei`: `82` successful sequence estimates, median K2P `0.0731`, median comparable sites `346`
- `D.perlapsus`: `53` successful sequence estimates, median K2P `0.0690`, median comparable sites `439`
- `D.adatsihi`: `50` successful sequence estimates, median K2P `0.0755`, median comparable sites `392`
- `D.santeetlah`: `45` successful sequence estimates, median K2P `0.0713`, median comparable sites `328`
- `D.campi`: `44` successful sequence estimates, median K2P `0.0768`, median comparable sites `352`
- `D.carolinensis`: `43` successful sequence estimates, median K2P `0.0695`, median comparable sites `297`
- `D.wrighti`: `43` successful sequence estimates, median K2P `0.0671`, median comparable sites `358`

## Interpretation

- The paired-LTR substrate exists locally and is substantial, so the biological question is tractable without falling back to divergence-to-consensus heuristics.
- The current local branch now supports direct 5'/3' LTR divergence estimation from assembly coordinates, which is the defensible basis for LTR insertion-age inference.
- Absolute age claims should remain conditional on an externally justified substitution rate; the default repo output stays at the divergence level.
- The primary-literature calibration layer now lives in `LTR_SUBSTITUTION_RATE_CALIBRATION.md`, with generated sensitivity tables under `results/data/ltr_age/`.
