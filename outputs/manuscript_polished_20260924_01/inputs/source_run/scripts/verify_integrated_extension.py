"""Independent integration checks: preservation, exports, tables and model scores."""
import csv
import hashlib
import json
import math
import re
import zipfile
from collections import defaultdict
from pathlib import Path
import xml.etree.ElementTree as ET

RUN = Path(__file__).resolve().parents[1]
def digest(path):
    with path.open('rb') as handle:
        return hashlib.file_digest(handle, 'sha256').hexdigest()
def read_csv(path):
    with (RUN/path).open() as handle:
        return list(csv.DictReader(handle))

baseline=json.loads((RUN/'qa/integration_before.json').read_text())
changed=[name for name,value in baseline.items() if not (RUN/name).is_file() or digest(RUN/name)!=value]
assert not changed, changed
document=(RUN/'index.html').read_text()
ids=re.findall(r'data-note-id="([^"]+)"',document)
assert len(ids)==len(set(ids))==64
assert all(re.fullmatch(r'(?:figure-[1-9][0-9]*(?:-[a-z][a-z0-9]*)?|table-[1-9][0-9]*)',i) for i in ids)

scores=read_csv('figures/figure_13/analysis/model_comparison.csv')
tests=read_csv('figures/figure_13/analysis/dsep_tests.csv')
claims=defaultdict(list)
for row in tests:
    claims[row['panel'],row['test_order'],row['model']].append(float(row['p_value']))
groups=defaultdict(list)
for row in scores:
    n,k,q=map(int,[row['n_species'],row['k'],row['q']])
    ps=claims[row['panel'],row['test_order'],row['model']]
    assert len(ps)==k
    C=-2*sum(math.log(p) for p in ps)
    CICc=C+2*q*n/(n-q-1)
    # Even-df chi-square survival, avoiding reuse of R's implementation.
    global_p=math.exp(-C/2)*sum((C/2)**j/math.factorial(j) for j in range(k))
    assert abs(C-float(row['C']))<1e-9
    assert abs(CICc-float(row['CICc']))<1e-9
    assert abs(global_p-float(row['p']))<1e-10
    groups[row['panel'],row['test_order']].append(row)
for rows in groups.values():
    best=min(float(r['CICc']) for r in rows)
    raw=[math.exp(-.5*(float(r['CICc'])-best)) for r in rows]
    assert all(abs(float(r['w'])-w/sum(raw))<1e-12 for r,w in zip(rows,raw))

# Check the exported workbook itself, not only its build payload or preview.
payload=json.loads((RUN/'tables/integrated/workbook_payload.json').read_text())
ns={'m':'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
typed_cells=0
def column(index):
    s=''
    while index:
        index,r=divmod(index-1,26)
        s=chr(65+r)+s
    return s
with zipfile.ZipFile(RUN/'tables/integrated/Manuscript_tables_integrated.xlsx') as z:
    strings=[]
    if 'xl/sharedStrings.xml' in z.namelist():
        strings=[''.join(si.itertext()) for si in ET.fromstring(z.read('xl/sharedStrings.xml'))]
    workbook=ET.fromstring(z.read('xl/workbook.xml'))
    names=[s.attrib['name'] for s in workbook.find('m:sheets',ns)]
    assert names==[s['name'] for s in payload['sheets']]
    used_styles=set()
    for index,spec in enumerate(payload['sheets'],1):
        root=ET.fromstring(z.read(f'xl/worksheets/sheet{index}.xml'))
        cells={c.attrib['r']:c for c in root.findall('.//m:c',ns)}
        used_styles.update(int(c.attrib.get('s',0)) for c in cells.values()
                           if c.find('m:v',ns) is not None or c.find('m:is',ns) is not None)
        for rownum,row in enumerate(spec['rows'],6):
            for cnum,expected in enumerate(row,1):
                c=cells[f'{column(cnum)}{rownum}']
                assert c.find('m:f',ns) is None
                v=c.find('m:v',ns)
                if c.attrib.get('t')=='s': actual=strings[int(v.text)]
                elif c.attrib.get('t')=='inlineStr': actual=''.join(c.find('m:is',ns).itertext())
                else: actual=v.text
                if isinstance(expected,(int,float)):
                    assert c.attrib.get('t') not in ('s','inlineStr')
                    assert math.isclose(float(actual),expected,rel_tol=1e-12,abs_tol=1e-13),(spec['name'],rownum,cnum)
                    typed_cells+=1
                else:
                    assert actual.strip()==('n.a.' if expected is None else expected),(spec['name'],rownum,cnum,actual,expected)
        assert not any(c.attrib.get('t')=='e' for c in cells.values())
        if len(spec['rows'])>15:
            assert root.find('.//m:pane',ns) is not None
    style=ET.fromstring(z.read('xl/styles.xml'))
    fonts=list(style.find('m:fonts',ns)); formats=list(style.find('m:cellXfs',ns))
    assert all(fonts[int(formats[i].attrib.get('fontId',0))].find('m:name',ns).attrib['val']=='Arial' for i in used_styles)

result={'passed':True,'prior_files_preserved':len(baseline),'changed_prior_files':changed,
        'note_fields':len(ids),'expanded_path_score_rows_verified':len(scores),
        'workbook_sheets':len(payload['sheets']),'workbook_numeric_cells_verified':typed_cells,
        'workbook_values_and_missingness_match_payload':True,'workbook_populated_cells_all_Arial':True}
(RUN/'qa/integrated_extension_verification.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
