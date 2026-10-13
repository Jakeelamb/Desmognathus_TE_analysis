"""Independent stdlib audit: run with uv run --frozen --no-sync python <this file>."""
from pathlib import Path
from collections import defaultdict
import csv, hashlib, json, math, statistics, subprocess
out = Path(__file__).resolve().parent
fig = out.parent
run = fig.parent.parent

def rows(p):
    with Path(p).open(newline='') as f:
        return list(csv.DictReader(f))

def near(a,b,tol=1e-8):
    assert abs(float(a)-float(b)) <= tol, (a,b)

snap = run/'inputs/animal_robustness'
for r in rows(snap/'MANIFEST.csv'):
    for key in ('source_path','snapshot_path'):
        assert hashlib.sha256(Path(r[key]).read_bytes()).hexdigest() == r['sha256']

support=rows(out/'species_support.csv')
for assay, filename in [('Paired morphology','morphology_animal_summaries.csv'),('Nuclear IOD','iod_animal_summaries.csv')]:
    animals=rows(out/filename)
    for s in (r for r in support if r['assay']==assay):
        a=[r for r in animals if r['tree_tip']==s['tree_tip']]
        counts=[int(r['n_objects']) for r in a]
        assert len(counts)==int(s['n_animals'])
        assert sum(counts)==int(s['n_objects'])
        near(max(counts)/sum(counts),s['max_animal_fraction'])
        if len(counts)>1: assert float(s['max_animal_fraction'])<1
assert next(float(r['max_animal_fraction']) for r in support if r['tree_tip']=='orestes' and r['assay']=='Paired morphology')==.98

sets=[('point_model_comparison.csv','point_dsep.csv',4),('bootstrap_model_comparison.csv','bootstrap_dsep.csv',250),('leave_one_species_out_model_comparison.csv','leave_one_species_out_dsep.csv',24)]
for rankfile,dsfile,nfits in sets:
    rankings=rows(out/rankfile); claims=rows(out/dsfile)
    byid=defaultdict(list)
    for r in rankings: byid[r['id']].append(r)
    assert len(byid)==nfits
    for identifier,group in byid.items():
        assert len(group)==4
        assert sorted(int(r['rank']) for r in group)==[1,2,3,4]
        near(sum(float(r['w']) for r in group),1)
        best=min(float(r['CICc']) for r in group)
        likelihoods=[math.exp(-(float(r['CICc'])-best)/2) for r in group]
        for r,likelihood in zip(group,likelihoods):
            ds=[c for c in claims if c['id']==identifier and c['model']==r['model']]
            k=int(r['k']);q=int(r['q']);n=int(r['n_species'])
            assert len(ds)==k
            C=-2*sum(math.log(float(c['p_value'])) for c in ds)
            near(C,r['C'])
            p=math.exp(-C/2)*sum((C/2)**j/math.factorial(j) for j in range(k))
            near(p,r['p'])
            near(C+2*q+2*q*(q+1)/(n-q-1),r['CICc'])
            near(likelihood/sum(likelihoods),r['w'])

boots=rows(out/'bootstrap_species_traits.csv'); grouped=defaultdict(list)
for r in boots: grouped[r['replicate']].append(r)
assert len(grouped)==250
for group in grouped.values():
    assert len(group)==24 and len({r['tree_tip'] for r in group})==24
    assert all(float(r[k])>0 for r in group for k in ('b','g','c','n'))
    near(next(r['g'] for r in group if r['tree_tip']=='fuscus'),16.1/.978)
assert len(rows(out/'failures.csv'))==0

pdfs=list(fig.glob('*.pdf'))+list((fig/'panels').glob('*.pdf'))
assert len(pdfs)==5
for pdf in pdfs:
    fonts=subprocess.check_output(['pdffonts',str(pdf)],text=True)
    records=fonts.splitlines()[2:]
    assert records and all('Arial' in r for r in records)
    text=subprocess.check_output(['pdftotext',str(pdf),'-'],text=True)
    assert 'Estimated' not in text and 'Conditional genome' not in text
    for ext in ('.png','.tif'):
        assert pdf.with_suffix(ext).stat().st_size>1000
    desc=subprocess.check_output(['identify','-format','%[colorspace],%x,%y',str(pdf.with_suffix('.tif'))],text=True)
    assert desc=='CMYK,300,300',desc
report=dict(passed=True,source_hashes=True,support_shares_recomputed_from_animal_counts=True,
    dsep_C_P_CICc_weights_replayed_fits=278,bootstrap_species_panels=250,
    pdfs=5,embedded_Arial=True,CMYK_300dpi_tiffs=5,rendered_label_check=True)
(out/'independent_verification.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report))
