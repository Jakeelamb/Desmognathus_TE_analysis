"""Build this run's review gallery from its own manifests and tables."""
import html
import hashlib
import importlib.util
import json
import sys
from pathlib import Path

RUN = Path(__file__).resolve().parents[1]
module_path = RUN / "inputs/publication_scripts/build_candidate_review.py"
spec = importlib.util.spec_from_file_location("isolated_gallery_builder", module_path)
module = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = module
spec.loader.exec_module(module)
destination = RUN / "index.html"
module.build_candidate_review([RUN / f"figures/figure_{n}" for n in range(1,7)], destination)
document = destination.read_text()
document = document.replace("Figure review — collages and isolated panels", "Manuscript candidates — isolated run")
document = document.replace("<h1>Figure review: collages and isolated panels</h1>", "<h1>Genome size, cellular architecture, and repeat evolution</h1><p class=run-label>Manuscript candidates · separate run · 24 September 2026</p>")
document = document.replace("Review numbers are provisional; final selection and refinement follow this draft pass.", "Two manuscript tables follow the figures. These are new manuscript numbers; the earlier gallery retains its original numbering.")
document = document.replace("to this project and remain available after refreshing or rebuilding the gallery.", "only within this separate run and remain available after refreshing or rebuilding this page.")
document = document.replace("Loading project notes…", "Loading this run’s notes…")
document = document.replace("</h1><p class=run-label>", "</h1><p class=run-label>")
intro = '''<div class="run-intro"><p><strong>Central question:</strong> Which causal hypotheses can explain the relationships among estimated genome size, nuclear size, and cell size?</p>
<p>The 24-species phylogenetic path comparison favors the class in which cell size connects genome-size and nuclear-size variation. Figure 6 shows the competing models, the conditional associations behind that result, and the causal directions that remain compatible under these models. Figures 1–4 establish the repeat-evolution context; Figure 5 shows the pairwise cellular relationships.</p>
<p class="run-links"><a href="#figure-6">Start with path analysis</a><a href="#figure-2">Repeat divergence</a><a href="#table-1">Manuscript tables</a><a href="http://127.0.0.1:8766/" target="_blank" rel="noopener">Earlier figure gallery ↗</a><a href="tables/Manuscript_tables.xlsx">Download editable tables</a></p>
<p class="separation">The earlier figures, data, and notes are preserved. Notes on this page belong to this run.</p></div>'''
document = document.replace('<nav aria-label="Figures and panels">', intro + '<nav aria-label="Figures and panels">')
document = document.replace('</nav></header>', '<div class="nav-group"><strong>Tables</strong><a href="#table-1">1 · Data coverage</a><a href="#table-2">2 · Comparative associations</a></div></nav></header>')
calibration = 'figures/figure_6/validation/Simulation_calibration'
calibration_version = hashlib.sha256((RUN / (calibration + '.png')).read_bytes()).hexdigest()[:12]
path_support = f'''<aside class="path-support"><details><summary>Path-analysis supporting results: simulation recovery and full numerical tables</summary>
<p>Recovery rates measure how often each generating model class ranks first in 100 simulations per residual setting. These are a check of discrimination, separate from the empirical model weights.</p>
<p><a href="{calibration}.pdf">Simulation PDF</a> · <a href="{calibration}.png?v={calibration_version}">Full-size PNG</a> · <a href="{calibration}.tif">CMYK TIFF</a> · <a href="inputs/datasets/Supplementary_Data_S30_path_model_comparison.csv">All model scores (CSV)</a> · <a href="inputs/datasets/Supplementary_Data_S28_pairwise_phylogenetic_regressions.csv">Regression coefficients (CSV)</a> · <a href="inputs/datasets/Supplementary_Data_S36_path_simulation_calibration.csv">Simulation results (CSV)</a></p>
<a href="{calibration}.png?v={calibration_version}"><img loading="lazy" src="{calibration}.png?v={calibration_version}" alt="Recovery of four generating path-model classes under Brownian and independent residuals, with 100 simulations per setting." /></a>
</details></aside>'''
document = document.replace('<section id="figure-6-a"', path_support + '<section id="figure-6-a"', 1)
payload = json.loads((RUN / "tables/tables.json").read_text())

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

def table_card(number,title,headers,body,caption,filename):
    anchor=f"table-{number}"
    cells = ''.join('<tr>'+''.join(f'<td>{c}</td>' for c in row)+'</tr>' for row in body)
    return f'''<section id="{anchor}" aria-labelledby="{anchor}-title"><h2 id="{anchor}-title">Table {number} — {esc(title)}</h2>
<p class="downloads"><a href="tables/{filename}">CSV</a> · <a href="tables/Manuscript_tables.xlsx">Editable workbook</a><button type="button" data-note-action="focus">Add / edit note</button></p>
<div class="table-scroll"><table><thead><tr>{''.join('<th>'+esc(h)+'</th>' for h in headers)}</tr></thead><tbody>{cells}</tbody></table></div>
<p class="caption">{esc(caption)}</p>{notes(anchor,'Table '+str(number))}</section>'''

coverage = [[esc(v if v is not None else '—') for v in row] for row in payload['coverage']]
table1=table_card(1,'Data coverage and measurement units',['Component','Species','Objects','Animals','Measurement','Source'],coverage,payload['coverage_notes'],'Table_1_data_coverage.csv')
associations=[]
for row in payload['associations']:
    p = row['p_value']
    p_text=f"{p:.2g}" if p < .001 else f"{p:.3f}"
    ci = '—' if row['ci_low'] is None else f"{row['ci_low']:.4f} to {row['ci_high']:.4f}"
    associations.append([esc(row['comparison']),str(row['n_species']),esc(row['effect_type']),format_num(row['estimate']),ci,p_text])
table2=table_card(2,'Focal comparative associations',['Comparison','Species','Statistic','Estimate','95% interval','P'],associations,payload['association_notes'],'Table_2_comparative_associations.csv')
document=document.replace('</main><script>',table1+table2+'</main><script>')
style='''
.run-label{color:#53606b;margin-top:4px}.run-intro{border-left:3px solid #0072b2;padding:2px 18px;margin:22px 0}.run-links{display:flex;gap:10px 24px;flex-wrap:wrap}.separation{font-size:14px;color:#53606b}
.path-support{padding:20px 28px;background:white;margin:20px 0;border:1px solid #dce3e7}.path-support summary{cursor:pointer;font-weight:bold}.path-support img{max-width:100%;height:auto}
.table-scroll{overflow-x:auto;margin:22px 0}table{width:100%;border-collapse:collapse;font-size:14px}th{text-align:left;padding:10px;background:#edf1f3;border-bottom:1px solid #7e8e99}td{padding:12px 10px;border-bottom:1px solid #e2e7eb;vertical-align:top}td:nth-child(2){text-align:right}tbody tr:last-child td{border-bottom:1px solid #7e8e99}
'''
document=document.replace('</style>',style+'</style>')
# A distinct browser draft namespace supplements the separate server-side store.
document=document.replace('desmognathus.figure-review.notes.v1', 'desmognathus.manuscript-story.20260924-01.notes.v1')
destination.write_text(document)
notes_file=RUN/'review_notes.json'
if not notes_file.exists():
    notes_file.write_text('{"schema_version": 1, "notes": {}}\n')
assert document.count('data-note-id=') >= 24
assert 'Table 2 — Focal comparative associations' in document
print('Built independent review gallery with six collages, individual panels, two tables, and separate notes')
