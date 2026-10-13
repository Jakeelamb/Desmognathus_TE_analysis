"""Verify figure exports, source preservation and local gallery targets."""
import argparse
import csv
import hashlib
import json
import re
import subprocess
from html.parser import HTMLParser
from pathlib import Path
from urllib.parse import unquote, urlsplit

import pdfplumber
from PIL import Image

RUN = Path(__file__).resolve().parents[1]
parser = argparse.ArgumentParser()
parser.add_argument('--publication-tiffs', action='store_true')
args = parser.parse_args()
errors = []
results = []
previous = json.loads((RUN / 'qa/source_baseline.json').read_text())
note_changes = []
for record in previous:
    p = Path(record['source'])
    same = p.exists() and hashlib.sha256(p.read_bytes()).hexdigest() == record['sha256']
    if not same:
        if p.name == 'review_notes.json':
            note_changes.append(str(p))
        else:
            errors.append(f'Source changed: {p}')
    if record['copied']:
        frozen = RUN / 'inputs/source_run' / record['relative']
        if hashlib.sha256(frozen.read_bytes()).hexdigest() != record['sha256']:
            errors.append(f'Frozen input changed: {frozen}')

for n in range(1, 9):
    out = RUN / f'figures/figure_{n}'
    main = list(csv.DictReader((out / 'FIGURE_MANIFEST.csv').open()))
    panels = list(csv.DictReader((out / 'panels/PANEL_MANIFEST.csv').open()))
    assert len(main) == 1 and panels
    assert [r['panel_id'] for r in panels] == list('abcdefghijklmnopqrstuvwxyz'[:len(panels)])
    assert float(main[0]['width_mm']) == 180 and float(main[0]['height_mm']) <= 245
    for folder, row in [(out, main[0])] + [(out / 'panels', r) for r in panels]:
        pdf, png, tif = [folder / row[k] for k in ('pdf','png','tiff')]
        if args.publication_tiffs:
            subprocess.run(['gs','-dSAFER','-dBATCH','-dNOPAUSE','-sDEVICE=tiff32nc',
                '-sCompression=lzw','-r600',f'-sOutputFile={tif}',str(pdf)],
                check=True, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
        with Image.open(png) as image:
            image.load()
            png_size = image.size
        with Image.open(tif) as image:
            image.load()
            assert image.mode == 'CMYK', (tif, image.mode)
            dpi = tuple(float(value) for value in image.info.get('dpi', (0, 0)))
            assert min(dpi) >= 599, (tif, dpi)
            tif_size = image.size
        with pdfplumber.open(pdf) as doc:
            assert len(doc.pages) == 1
            page = doc.pages[0]
            text = page.extract_text() or ''
            assert not re.search(r'(conditional|estimated)\s+genome\s+size', text, re.I), pdf
            upright = [c for c in page.chars if c.get('upright') and c['text'].strip()]
            min_font = min(c['size'] for c in upright)
            if min_font < 7.35:
                errors.append(f'Text below 7.5pt target: {pdf} ({min_font:.3f})')
            clipped = [c for c in page.chars if c['text'].strip() and
                (c['x0'] < -.5 or c['x1'] > page.width + .5 or c['top'] < -.5 or c['bottom'] > page.height + .5)]
            if clipped:
                errors.append(f'Text outside page: {pdf} ({len(clipped)} characters)')
            # Cairo rounds its media box to whole PDF points (less than 0.36 mm).
            assert abs(page.width - float(row['width_mm']) / 25.4 * 72) < 1
            assert abs(page.height - float(row['height_mm']) / 25.4 * 72) < 1
        fonts = subprocess.check_output(['pdffonts',str(pdf)],text=True).splitlines()[2:]
        assert fonts
        for line in fonts:
            if 'Arial' not in line or not re.search(r'\byes\s+yes\b', line):
                errors.append(f'Nonembedded or non-Arial font in {pdf}: {line}')
        expected = (round(page.width / 72 * 300), round(page.height / 72 * 300))
        assert max(abs(a-b) for a,b in zip(png_size, expected)) <= 1, (png, png_size, expected)
        expected_tiff = (round(page.width / 72 * dpi[0]), round(page.height / 72 * dpi[1]))
        assert max(abs(a-b) for a,b in zip(tif_size, expected_tiff)) <= 1, (tif, tif_size, expected_tiff)
        results.append({'pdf':str(pdf.relative_to(RUN)), 'width_mm':float(row['width_mm']),
            'height_mm':float(row['height_mm']), 'actual_pdf_width_mm':page.width / 72 * 25.4,
            'actual_pdf_height_mm':page.height / 72 * 25.4, 'minimum_upright_font_pt':min_font,
            'font_families':'Arial embedded', 'png_dimensions':png_size, 'tiff_dimensions':tif_size,
            'tiff_color':'CMYK','tiff_dpi':dpi,'text_within_page':not clipped})

class Links(HTMLParser):
    def __init__(self):
        super().__init__(); self.links=[]; self.ids=[]; self.notes=[]
    def handle_starttag(self, tag, attrs):
        a=dict(attrs)
        for k in ('href','src'):
            if k in a:self.links.append(a[k])
        if 'id' in a:self.ids.append(a['id'])
        if 'data-note-id' in a:self.notes.append(a['data-note-id'])

links = Links(); links.feed((RUN / 'index.html').read_text())
assert len(links.ids) == len(set(links.ids))
assert len(links.notes) == len(results) and len(set(links.notes)) == len(links.notes)
checked=0
for link in links.links:
    u=urlsplit(link)
    if u.scheme or u.netloc:continue
    if u.path:
        if not (RUN / unquote(u.path)).exists():errors.append(f'Broken local link: {link}')
        checked+=1
    elif u.fragment and u.fragment not in links.ids:errors.append(f'Broken anchor: {link}')
for script in (RUN / 'scripts').glob('*.R'):
    executable='\n'.join(x for x in script.read_text().splitlines() if not x.lstrip().startswith('#'))
    if re.search(r'\b(?:phylolm|phylo_path|glm|lm)\s*\(',executable):errors.append(f'Unexpected model fit in {script}')
report={'status':'PASS' if not errors else 'FAIL','main_figures':8,'figure_and_panel_exports':len(results),
    'preserved_source_files':len(previous),'frozen_input_files':sum(r['copied'] for r in previous),
    'user_note_changes':note_changes,'note_fields':len(links.notes),'local_asset_links':checked,
    'errors':errors,'artifacts':results}
(RUN / 'qa/verification.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps({k:v for k,v in report.items() if k!='artifacts'},indent=2))
if errors:raise SystemExit(1)
