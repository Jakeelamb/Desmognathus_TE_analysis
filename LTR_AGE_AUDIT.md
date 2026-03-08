# LTR Age Audit

This note records the audit status of the local LTR insertion-age branch.

## Bottom Line

- The previous local approach was not suitable for paper use because it tried to infer insertion age from RepeatMasker divergence-to-consensus hits rather than divergence between true paired 5' and 3' LTR sequences.
- The canonical local paired-element source is the ectopic master table produced by `scripts/processing/ec.py`, not the RepeatMasker adjacency heuristic.
- Local sequence-based age estimation is currently blocked because matching genome FASTA assemblies are not present in the repository workspace.

## Current Local Readiness

- Valid paired LTR elements in canonical ectopic master table: `1371`
- Species with at least one paired LTR element: `31`
- Complete (`Complete = yes`) LTR elements: `247`
- LTR pairs with both LTRs >= `100` bp: `1371`
- LTR pairs with five or more annotated domains: `1158`
- Dominant superfamily among paired LTR elements: `Gypsy` (`1363` elements)
- Species currently ready for true sequence-based age estimation with local genome FASTA present: `0`

## Highest-Coverage Species

- `D.tilleyi`: `142` paired LTR elements, `36` complete, median internal length `4754.5` bp
- `D.abditus`: `98` paired LTR elements, `19` complete, median internal length `4756.0` bp
- `D.welteri`: `94` paired LTR elements, `14` complete, median internal length `4689.5` bp
- `D.valentinei`: `82` paired LTR elements, `12` complete, median internal length `4833.5` bp
- `Unknown`: `80` paired LTR elements, `9` complete, median internal length `4863.0` bp
- `D.perlapsus`: `53` paired LTR elements, `7` complete, median internal length `4945.0` bp
- `D.adatsihi`: `50` paired LTR elements, `4` complete, median internal length `4897.5` bp
- `D.santeetlah`: `45` paired LTR elements, `11` complete, median internal length `4836.0` bp
- `D.campi`: `44` paired LTR elements, `11` complete, median internal length `4916.5` bp
- `D.carolinensis`: `43` paired LTR elements, `4` complete, median internal length `4985.0` bp

## Interpretation

- The paired-LTR substrate exists locally and is substantial, so the biological question is still tractable.
- What is missing is not paired-element annotation but the sequence-access layer required to compare the 5' and 3' LTRs directly.
- Until genome FASTA assemblies are available locally and wired into a sequence-extraction workflow, this branch should be treated as `blocked for paper-ready age inference`.

