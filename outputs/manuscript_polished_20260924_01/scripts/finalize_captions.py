"""Apply consistent finding-first captions and define displayed abbreviations."""
import csv
import re
from pathlib import Path

RUN = Path(__file__).resolve().parents[1]
FINDINGS = {
    1: 'Classified repeat alignments share an LTR-dominated architecture across the 34 sampled species.',
    2: 'Within-LTR low-divergence representation differs among reproductive groups in an exploratory phylogenetic comparison.',
    3: 'Dominant-superfamily evenness tends to decline with genome size across the 21 overlapping species.',
    4: 'The selected LTR terminal/internal depth proxy does not resolve a relationship with genome size.',
    5: 'Genome size covaries with erythrocyte area, while nuclear and cellular dimensions are not uniformly proportional.',
    6: 'Compatible models preserve a genome–cell–nucleus connection while leaving causal direction and body-size placement unresolved.',
    7: 'Genome–cell associations persist across body-size adjustment, animal weighting and the specified species sensitivity.',
    8: 'The measured LTR divergence signature does not account for the genome–cell association in the fitted models.',
}
DEFINITIONS = {
    'TE': 'TE, transposable element',
    'LTR': 'LTR, long terminal repeat',
    'LINE': 'LINE, long interspersed nuclear element',
    'TIR': 'TIR, terminal inverted repeat',
    'SVL': 'SVL, snout–vent length',
    'PGLS': 'PGLS, phylogenetic generalized least squares',
    'CI': 'CI, confidence interval',
    'CICc': 'CICc, C-statistic information criterion with small-sample correction',
}

def finalize():
    for n in range(1, 9):
        out = RUN / f'figures/figure_{n}'
        for manifest in (out / 'FIGURE_MANIFEST.csv', out / 'panels/PANEL_MANIFEST.csv'):
            with manifest.open(newline='') as handle:
                reader = csv.DictReader(handle)
                fields = reader.fieldnames
                rows = list(reader)
            for row in rows:
                legend = row['legend']
                if n == 6:
                    legend = legend.removeprefix('Figure 6. Cellular associations and alternative causal hypotheses. ')
                if not legend.startswith(FINDINGS[n]):
                    legend = FINDINGS[n] + ' ' + legend
                if 'Abbreviations:' not in legend:
                    glossary = [definition for term, definition in DEFINITIONS.items() if re.search(r'\b' + term + r'\b', legend)]
                    if glossary:
                        legend += ' Abbreviations: ' + '; '.join(glossary) + '.'
                row['legend'] = legend
            with manifest.open('w', newline='') as handle:
                writer = csv.DictWriter(handle, fieldnames=fields)
                writer.writeheader(); writer.writerows(rows)
            if manifest.name == 'FIGURE_MANIFEST.csv':
                row = rows[0]
                (out / 'FIGURE_LEGENDS_AND_ALT_TEXT.md').write_text(f'# Fig. {n}. {row["title"]}\n\n{row["legend"]}\n\nAlt text: {row["alt_text"]}\n')
            else:
                (out / 'panels/PANEL_LEGENDS_AND_ALT_TEXT.md').write_text('\n\n'.join(
                    f'# {r["panel_id"]}. {r["title"]}\n\n{r["legend"]}\n\nAlt text: {r["alt_text"]}' for r in rows) + '\n')

if __name__ == '__main__':
    finalize()
