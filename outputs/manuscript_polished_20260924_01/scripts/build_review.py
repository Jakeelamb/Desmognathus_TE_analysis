"""Build eight polished collages and isolated panels with an independent note store."""
import csv
import html
import importlib.util
import json
import sys
from pathlib import Path
from finalize_captions import finalize

RUN = Path(__file__).resolve().parents[1]
SOURCE = RUN / 'inputs/source_run'
path = SOURCE / 'inputs/publication_scripts/build_candidate_review.py'
spec = importlib.util.spec_from_file_location('polished_gallery', path)
module = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = module
spec.loader.exec_module(module)
target = RUN / 'index.html'
finalize()
module.build_candidate_review([RUN / f'figures/figure_{n}' for n in range(1, 9)], target)
doc = target.read_text()
doc = doc.replace('Figure review — collages and isolated panels', 'Desmognathus — eight polished figures')
doc = doc.replace('<h1>Figure review: collages and isolated panels</h1>', '<h1>Genome size, cellular architecture and repeat evolution</h1><p class="eyebrow">Eight-figure publication refinement · separate review version</p>')
doc = doc.replace('draft figures ·', 'main figures ·')
doc = doc.replace('Review numbers are provisional; final selection and refinement follow this draft pass.', 'Arial throughout · fixed species order and biological colors · frozen numerical results.')
doc = doc.replace('to this project and remain available after refreshing or rebuilding the gallery.', 'to this eight-figure version. Earlier notes remain in their original galleries.')
doc = doc.replace('Loading project notes…', 'Loading this version’s notes…')
doc = doc.replace('desmognathus.figure-review.notes.v1', 'desmognathus.polished-eight.20260924.notes.v1')
intro = '''<aside class="story"><p><strong>Central result:</strong> Genome-size estimates and erythrocyte dimensions covary beyond maximum body size. Nuclear and cellular dimensions are not uniformly proportional. Repeat-divergence profiles differ among reproductive groups, while the measured repeat summary does not account for the cellular association.</p>
<p><strong>Reading the arrows:</strong> Figure 6 displays coefficients within a proposed causal orientation and the other directions compatible with its independence claims. Relative model support is not a probability that an arrow is causal.</p>
<p>Every observation, fit and exclusion rule is preserved. The collages use selected species labels for legibility; named isolated panels and frozen source tables retain full inspection detail. No heatmap or outlier filtering has been added.</p>
<p class="links"><a href="http://127.0.0.1:8767/">Earlier 13-figure gallery and notes</a><a href="http://127.0.0.1:8766/">Original gallery</a><a href="FINDINGS_AND_CAPTIONS.md">Findings and captions</a><a href="NOTE_DECISIONS.md">How prior notes were applied</a><a href="README.md">Sources and rebuild</a></p>
<p>Eight figures use GBE’s full main-text display allowance; the existing tables are retained as supporting material. <a href="http://127.0.0.1:8767/#table-1">Review all five tables</a>.</p></aside>'''
doc = doc.replace('<nav aria-label="Figures and panels">', intro + '<nav aria-label="Figures and panels">')
if (RUN / 'Manuscript_figures_1-8.zip').exists():
    doc = doc.replace('<a href="FINDINGS_AND_CAPTIONS.md">', '<a href="Manuscript_figures_1-8.zip">Download all eight collages</a><a href="FINDINGS_AND_CAPTIONS.md">')
messages = {
  1: ('A shared repeat architecture', 'LTRs dominate classified repeat alignments across all 34 species. These percentages describe classified aligned bases, not whole-genome occupancy.'),
  2: ('Reproductive groups differ in LTR profiles', 'Within-LTR normalization separates the low-divergence signature from total LTR representation. Strategy association remains exploratory; the adjusted direct-development contrast is uncertain.'),
  3: ('A tendency toward lower dominant-family evenness', 'Both indices describe the same ten-superfamily composition. Their agreement is not independent replication; the Shannon interval crosses zero.'),
  4: ('The terminal/internal proxy does not resolve a genome-size relationship', 'All selected elements are retained, including extremes. The internal-region definition and mapping limitations prevent interpreting this as a calibrated deletion rate.'),
  5: ('Cellular covariation without fixed proportionality', 'Genome–cell association is the clearest pairwise pattern. Nuclear area increases less steeply in the displayed conditional regression; paired nuclear fractions vary among species.'),
  6: ('Causal models constrain explanations, while direction remains unresolved', 'The leading expanded classes preserve the cellular connection but differ in body-size placement. Weighted arrows express a specified hypothesis; equivalent orientations remain possible.'),
  7: ('Body context and animal sampling', 'The two smallest direct developers influence organismal associations. Positive genome–cell and cell–nucleus coefficients persist in the stated body-adjusted sensitivities.'),
  8: ('The measured TE signature does not explain the cellular relationship', 'The exact 21-species comparison retains the genome–cell coefficient after TE adjustment. Wide TE intervals leave historical insertion and deletion mechanisms unresolved.'),
}
captions = ['# Eight-figure findings and captions', '', 'Refinement of frozen results; no models refitted. Figures are descriptive or exploratory as specified below.', '']
for n, (finding, interpretation) in messages.items():
    marker = f'<section id="figure-{n}"'
    banner = f'<aside class="finding"><strong>{html.escape(finding)}</strong><p>{html.escape(interpretation)}</p></aside>'
    doc = doc.replace(marker, banner + marker, 1)
    row = next(csv.DictReader((RUN / f'figures/figure_{n}/FIGURE_MANIFEST.csv').open()))
    captions += [f'## Fig. {n}. {row["title"]}', '', '**Finding:** ' + interpretation, '', row['legend'], '', 'Alt text: ' + row['alt_text'], '']
doc = doc.replace('</style>', '.eyebrow{color:#586674}.story,.finding{padding:16px 22px;border-left:3px solid #0072b2;background:#fff;margin:20px 0}.links{display:flex;gap:12px 22px;flex-wrap:wrap}.finding p{margin-bottom:0}.finding{max-width:1056px;margin:20px auto}.caption{font-size:15px}img{max-width:1000px}nav{border-top:1px solid #dce3e7;padding-top:12px}</style>')
target.write_text(doc)
(RUN / 'FINDINGS_AND_CAPTIONS.md').write_text('\n'.join(captions))
notes = RUN / 'review_notes.json'
if not notes.exists():
    notes.write_text(json.dumps({'schema_version': 1, 'notes': {}}, indent=2) + '\n')
print(f'Built 8 figures with {doc.count("data-note-id=")} independent note fields')
