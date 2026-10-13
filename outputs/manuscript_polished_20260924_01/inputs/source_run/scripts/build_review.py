"""Build this run's review gallery from its own manifests and tables."""
import html
import csv
import hashlib
import importlib.util
import json
import re
import sys
from pathlib import Path

RUN = Path(__file__).resolve().parents[1]
module_path = RUN / "inputs/publication_scripts/build_candidate_review.py"
spec = importlib.util.spec_from_file_location("isolated_gallery_builder", module_path)
module = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = module
spec.loader.exec_module(module)
destination = RUN / "index.html"
module.build_candidate_review([RUN / f"figures/figure_{n}" for n in range(1,14)], destination)
document = destination.read_text()
document = document.replace("Figure review — collages and isolated panels", "Manuscript candidates — isolated run")
document = document.replace("<h1>Figure review: collages and isolated panels</h1>", "<h1>Genome size, cellular architecture, and repeat evolution</h1><p class=run-label>Manuscript candidates · separate run · 24 September 2026</p>")
document = document.replace("Review numbers are provisional; final selection and refinement follow this draft pass.", "Five consolidated manuscript tables follow the figures. Numbers remain provisional for review; the earlier gallery retains its original numbering.")
document = document.replace("to this project and remain available after refreshing or rebuilding the gallery.", "only within this separate run and remain available after refreshing or rebuilding this page.")
document = document.replace("Loading project notes…", "Loading this run’s notes…")
document = document.replace("</h1><p class=run-label>", "</h1><p class=run-label>")
intro = '''<div class="run-intro"><p><strong>Central question:</strong> Which hypotheses can explain the relationships among genome size, nuclear size, cell size, and organismal body size?</p>
<p>Figures 9–13 add animal-level sampling sensitivity, nuclear–cell scaling, an exact 21-species TE–cellular comparison, reproductive-strategy differences in within-LTR divergence, and alternative positions of body size in the path analysis. Figures 1–8 remain available for comparison.</p>
<p class="run-links"><a href="#figure-9">New: animal-level robustness</a><a href="#figure-10">New: nuclear–cell scaling</a><a href="#figure-11">New: TE–cellular integration</a><a href="#figure-12">New: divergence and development</a><a href="#figure-13">New: alternative body-size paths</a><a href="#support-composition">Composition supporting figures</a><a href="#table-1">Five manuscript tables</a><a href="http://127.0.0.1:8766/" target="_blank" rel="noopener">Earlier figure gallery ↗</a><a href="tables/integrated/Manuscript_tables_integrated.xlsx">Download consolidated workbook</a></p>
<p class="separation">The earlier figures, data, and notes are preserved. Notes on this page belong to this run.</p></div>'''
document = document.replace('<nav aria-label="Figures and panels">', intro + '<nav aria-label="Figures and panels">')
document = document.replace('</nav></header>', '<div class="nav-group"><strong>Supporting figures</strong><a href="#support-composition">S1 · Repeat composition</a></div><div class="nav-group"><strong>Tables</strong><a href="#table-1">1 · Data coverage</a><a href="#table-2">2 · Focal associations</a><a href="#table-3">3 · Path comparisons</a><a href="#table-4">4 · Species measurements</a><a href="#table-5">5 · Sampling sensitivity</a></div></nav></header>')
document = document.replace('<div class="nav-group"><strong>Tables</strong>', '<div class="nav-group"><strong>Alternative path display</strong><a href="#figure-6-network">Direct + indirect collage</a><a href="#figure-6-networka">Network</a><a href="#figure-6-networkb">Intervals</a></div><div class="nav-group"><strong>Tables</strong>')
calibration = 'figures/figure_6/validation/Simulation_calibration'
calibration_version = hashlib.sha256((RUN / (calibration + '.png')).read_bytes()).hexdigest()[:12]
contrast = 'figures/figure_6/validation/Marginal_conditional_comparison'
contrast_version = hashlib.sha256((RUN / (contrast + '.png')).read_bytes()).hexdigest()[:12]
weighted_revision = '''<aside class="path-support"><strong>Updated Figure 6: weighted path diagrams</strong>
<p>Arrow labels show standardized coefficients; the shared class weight shows model support. Panel c gives the coefficient intervals and conditional comparisons.</p>
<p><a href="revisions/figure_6_before_weighted_arrows/index.html" target="_blank" rel="noopener">View the previous Figure 6 and saved notes ↗</a></p></aside>'''
document = document.replace('<section id="figure-6"', weighted_revision + '<section id="figure-6"', 1)
path_support = f'''<aside class="path-support"><details><summary>Path-analysis supporting results: simulation recovery and full numerical tables</summary>
<p>Recovery rates measure how often each generating model class ranks first in 100 simulations per residual setting. These are a check of discrimination, separate from the empirical model weights.</p>
<p><a href="{calibration}.pdf">Simulation PDF</a> · <a href="{calibration}.png?v={calibration_version}">Full-size PNG</a> · <a href="{calibration}.tif">CMYK TIFF</a> · <a href="inputs/datasets/Supplementary_Data_S30_path_model_comparison.csv">All model scores (CSV)</a> · <a href="inputs/datasets/Supplementary_Data_S28_pairwise_phylogenetic_regressions.csv">Regression coefficients (CSV)</a> · <a href="inputs/datasets/Supplementary_Data_S36_path_simulation_calibration.csv">Simulation results (CSV)</a></p>
<a href="{calibration}.png?v={calibration_version}"><img loading="lazy" src="{calibration}.png?v={calibration_version}" alt="Recovery of four generating path-model classes under Brownian and independent residuals, with 100 simulations per setting." /></a>
<h3>Original marginal and conditional comparison</h3>
<p><a href="{contrast}.pdf">PDF</a> · <a href="{contrast}.png?v={contrast_version}">Full-size PNG</a> · <a href="{contrast}.tif">CMYK TIFF</a> · <a href="figures/figure_6/analysis/dag_edges.csv">All displayed arrow estimates (CSV)</a></p>
<a href="{contrast}.png?v={contrast_version}"><img loading="lazy" src="{contrast}.png?v={contrast_version}" alt="Marginal and conditional genome–nucleus and genome–cell regression coefficients and intervals." /></a>
</details></aside>'''
document = document.replace('<section id="figure-6-a"', path_support + '<section id="figure-6-a"', 1)
payload = json.loads((RUN / "tables/integrated/workbook_payload.json").read_text())

def esc(value):
    return html.escape(str(value))

def format_num(value, places=4):
    if value is None:
        return "—"
    return f"{value:.{places}f}"

def notes(anchor,title):
    return f'''<div class="figure-notes"><label for="note-{anchor}">Notes for {esc(title)}</label>
<textarea id="note-{anchor}" class="note-input" data-note-id="{anchor}" rows="4" maxlength="50000" aria-describedby="note-status-{anchor}" disabled placeholder="What should change?"></textarea>
<p class="note-status" id="note-status-{anchor}" role="status">Loading saved notes…</p>
<div class="note-conflict" hidden><p>A newer note was saved elsewhere. The saved version is:</p><pre class="remote-note"></pre><div class="note-actions"><button type="button" data-note-action="use-saved">Use saved version</button><button type="button" data-note-action="keep-draft">Save my draft instead</button></div></div></div>'''

network_root = 'figures/figure_6/alternatives/direct_indirect'
network_legend = (RUN / network_root / 'LEGEND_AND_ALT_TEXT.md').read_text().split('\n\n')[1]
network_alt = 'Genome-first triangle: genome to cell 0.72, cell to nucleus 0.57, and genome to nucleus 0.008. The direct genome–nucleus path is dashed and has a wide interval spanning zero.'
network_cards = []
for anchor, title, base in [
    ('figure-6-network', 'Alternative Figure 6 — Direct and indirect paths', 'Direct_indirect_candidate'),
    ('figure-6-networka', 'Alternative Figure 6a — Reference-style network', 'Direct_indirect_network'),
    ('figure-6-networkb', 'Alternative Figure 6b — Intervals for the network paths', 'Direct_indirect_intervals'),
]:
    prefix = f'{network_root}/{base}'
    version = hashlib.sha256((RUN / (prefix + '.png')).read_bytes()).hexdigest()[:12]
    network_cards.append(f'''<section id="{anchor}" aria-labelledby="{anchor}-title"><h2 id="{anchor}-title">{esc(title)}</h2>
<p>Separate candidate inspired by the supplied paper figure. The assumed ordering is genome → cell → nucleus, with the direct genome → nucleus path retained. This fully connected model is not assigned the reduced chain's model weight.</p>
<p class="downloads"><a href="{prefix}.pdf?v={version}">PDF</a> · <a href="{prefix}.png?v={version}">Full-size PNG</a> · <a href="{prefix}.tif?v={version}">CMYK TIFF</a> · <a href="{network_root}/displayed_path_estimates.csv">Coefficients CSV</a><button type="button" data-note-action="focus">Add / edit note</button></p>
<a href="{prefix}.png?v={version}"><img src="{prefix}.png?v={version}" alt="{esc(network_alt)}" loading="lazy" /></a>
{notes(anchor, title)}<p class="caption">{esc(network_legend)}</p></section>''')
document = document.replace('<section id="figure-6-b"', ''.join(network_cards) + '<section id="figure-6-b"', 1)
document = re.sub(r'(\d+) draft figures · (\d+) isolated panels\.',
    lambda m: f'{m[1]} main figure candidates + 2 supporting collages · {int(m[2]) + 5} isolated panels.', document)
body_root = 'figures/figure_8/analysis'
body_support = f'''<aside class="path-support"><strong>New body-size analysis: inputs, models and numerical results</strong>
<p>Maximum SVL is a species-level body-size proxy. The full panel contains 24 species; the 22-species sensitivity omits only D. aeneus and D. wrighti. It does not fit a separate two-species path model. All four candidate hypotheses treat body size as a shared predictor; weights compare the remaining genome–cell–nucleus structure within each panel.</p>
<p><a href="{body_root}/traits_joined.csv">Joined traits CSV</a> · <a href="{body_root}/regression_coefficients.csv">Regression coefficients CSV</a> · <a href="{body_root}/path_model_comparison.csv">Candidate comparison CSV</a> · <a href="{body_root}/path_edge_estimates.csv">Path coefficients CSV</a> · <a href="{body_root}/README.md">Methods and interpretation</a></p>
</aside>'''
document = document.replace('<section id="figure-7"', body_support + '<section id="figure-7"', 1)

new_support = {
    9: ('Animal-level sensitivity', 'point_coefficients.csv', 'README.md'),
    10: ('Paired nuclear–cell scaling', 'allometry_model_summary.csv', 'README.md'),
    11: ('Exact 21-species integration', 'regression_coefficients.csv', 'README.md'),
    12: ('Within-LTR divergence and reproductive strategy', 'strategy_omnibus_tests.csv', 'README.md'),
    13: ('Alternative roles for body size', 'model_comparison.csv', 'METHODS.md'),
}
for number, (title, csv_file, methods_file) in new_support.items():
    root = f'figures/figure_{number}/analysis'
    support = f'<aside class="path-support"><strong>{esc(title)}</strong><p><a href="{root}/{csv_file}">Numerical results CSV</a> · <a href="{root}/{methods_file}">Methods and interpretation</a></p></aside>'
    document = document.replace(f'<section id="figure-{number}"', support + f'<section id="figure-{number}"', 1)

# Reuse the preserved composition candidate, with its own note identifiers.
support_root = RUN / 'figures/support_composition'
with (support_root / 'FIGURE_MANIFEST.csv').open() as handle:
    composition = next(csv.DictReader(handle))
with (support_root / 'panels/PANEL_MANIFEST.csv').open() as handle:
    composition_panels = list(csv.DictReader(handle))
composition_cards = []
for anchor, title, folder, row in [('support-composition', 'Supporting S1 — Repeat composition and reproductive strategy', support_root, composition)] + [
    (f'support-composition{r["panel_id"]}', f'Supporting S1{r["panel_id"]} — {r["title"]}', support_root / 'panels', r) for r in composition_panels
]:
    prefix = str(folder.relative_to(RUN))
    version = hashlib.sha256((folder / row['png']).read_bytes()).hexdigest()[:12]
    links = ' · '.join(f'<a href="{prefix}/{row[key]}?v={version}">{label}</a>' for key, label in [('pdf','PDF'),('png','Full-size PNG'),('tiff','CMYK TIFF')])
    composition_cards.append(f'<section id="{anchor}"><h2>{esc(title)}</h2><p>Preserved composition candidate. The Bonsai topology summarizes repeat composition; it is not the species phylogeny. Group annotations were supplied after the initial ordination.</p><p class="downloads">{links}<button type="button" data-note-action="focus">Add / edit note</button></p><a href="{prefix}/{row["png"]}?v={version}"><img loading="lazy" src="{prefix}/{row["png"]}?v={version}" alt="{esc(row["alt_text"])}" /></a>{notes(anchor,title)}<p class="caption">{esc(row["legend"])}</p></section>')

def table_card(number,title,headers,body,caption,filename):
    anchor=f"table-{number}"
    cells = ''.join('<tr>'+''.join(f'<td>{c}</td>' for c in row)+'</tr>' for row in body)
    return f'''<section id="{anchor}" aria-labelledby="{anchor}-title"><h2 id="{anchor}-title">Table {number} — {esc(title)}</h2>
<p class="downloads"><a href="tables/integrated/{filename}">CSV</a> · <a href="tables/integrated/Manuscript_tables_integrated.xlsx">Editable workbook</a><button type="button" data-note-action="focus">Add / edit note</button></p>
<div class="table-scroll"><table><thead><tr>{''.join('<th>'+esc(h)+'</th>' for h in headers)}</tr></thead><tbody>{cells}</tbody></table></div>
<p class="caption">{esc(caption)}</p>{notes(anchor,'Table '+str(number))}</section>'''

def display_cell(value, number_format=None):
    if value is None:
        return '—'
    if isinstance(value, (int, float)):
        if number_format and '%' in number_format:
            return f'{value:.1%}'
        if isinstance(value, int):
            return f'{value:,}'
        if number_format and 'E' in number_format:
            return f'{value:.3g}'
        if value and abs(value) < .001:
            return f'{value:.3g}'
        if number_format and '0.0000' in number_format:
            return f'{value:.4f}'
        return f'{value:.3f}'
    return esc(value)

table_cards=[]
for number, table in enumerate(payload['sheets'], 1):
    body=[[display_cell(v, table['formats'].get(str(c))) for c,v in enumerate(row)] for row in table['rows']]
    table_cards.append(table_card(number, table['title'].split('. ',1)[1], table['headers'], body, table['caption'], table['csv']))
document=document.replace('</main><script>', ''.join(composition_cards)+''.join(table_cards)+'</main><script>')
style='''
.run-label{color:#53606b;margin-top:4px}.run-intro{border-left:3px solid #0072b2;padding:2px 18px;margin:22px 0}.run-links{display:flex;gap:10px 24px;flex-wrap:wrap}.separation{font-size:14px;color:#53606b}
.path-support{padding:20px 28px;background:white;margin:20px 0;border:1px solid #dce3e7}.path-support summary{cursor:pointer;font-weight:bold}.path-support img{max-width:100%;height:auto}
.table-scroll{overflow-x:auto;margin:22px 0}table{width:100%;border-collapse:collapse;font-size:14px}th{text-align:left;padding:10px;background:#edf1f3;border-bottom:1px solid #7e8e99}td{padding:12px 10px;border-bottom:1px solid #e2e7eb;vertical-align:top}td:nth-child(2){text-align:right}tbody tr:last-child td{border-bottom:1px solid #7e8e99}
'''
document=document.replace('</style>',style+'</style>')
# A distinct browser draft namespace supplements the separate server-side store.
document=document.replace('desmognathus.figure-review.notes.v1', 'desmognathus.manuscript-story.20260924-01.notes.v1')
# Stable supplementary anchors also satisfy the unchanged note server contract.
document=document.replace('support-composition', 'figure-1-composition')
destination.write_text(document)
notes_file=RUN/'review_notes.json'
if not notes_file.exists():
    notes_file.write_text('{"schema_version": 1, "notes": {}}\n')
assert document.count('data-note-id=') == 64
assert 'Table 2 — Focal comparative associations' in document
print('Built review gallery: 13 main collages, alternative path, supporting composition, isolated panels, five tables, 64 note fields')
