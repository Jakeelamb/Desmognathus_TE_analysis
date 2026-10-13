# Independent online record recovery verification

Verdict: **PASS for retained record recovery and author-attested one-run membership; HOLD for genome calibration and standards QC.**

Verified repository HEAD: `15dfb3c18755cc8659b626fa213087c62c6db02f`. The workspace is dirty and the evidence directory is untracked, so this verdict identifies the actual artifacts by the SHA256 values below rather than treating HEAD as their complete revision. Final scientific-record checks ran at 2026-10-13T17:12:27.173260+00:00 to 2026-10-13T17:12:27.449764+00:00; current coordinator summaries were read and checked at 2026-10-13T17:13:39.671811+00:00.

Only this evidence report was written by this verifier. No browser/cloud calls, subagents, record builder, scientific producer, source/Methods/analysis/figure mutation, or acceptance-test change was performed. The populated CSV was inspected independently rather than regenerated.

## Requirement evidence

| Requirement | Observed result | Exact evidence |
|---|---|---|
| All retained targets and both standards | 53 unique filenames: 47 target + 4 current species denominator + 2 designated standard rows; exact 51-image coverage of nuclear_iod_by_image.csv, 50 retained specimen IDs, 805 retained nuclei | STAINING_RUN_MAP_AUTHOR_CONFIRMED.csv, pinned template, current nuclear_iod_by_image.csv |
| Original identifiers and roles preserved | Every species/filename/specimen_id/role tuple matches the original template in order; no original identifiers or roles changed | STAINING_RUN_MAP_TO_COMPLETE.csv SHA256 in check output |
| Four current denominator images | Processes 314/315/320/321, specimens 3372/3373/3401/3402, remain current species denominator | Six exact denominator/standard rows in summary check output below |
| Two designated standards | Process_413_green.ome.tiff → specimen 32469 → glass 1321; Process_414_green.ome.tiff → specimen 32470 → glass 1322 | master_image_metadata.csv lines 123/124; slide_glass_mapping.csv lines 83/84 |
| One-run assignment and provenance | All 53 rows use author_confirmed_single_run, author attestation, and AUTHOR_STAINING_CONFIRMATION.json; both standards paired for every row | Recorded exact answer: “One staining run with both standards”; supplied question covers all retained brightfield slides and both specimens |
| No invented experimental date/name | All 53 staining_date values blank; actual_staining_date and original_experimental_run_identifier null; run_id_is_internal_label true | CSV and AUTHOR_STAINING_CONFIRMATION.json |
| Exact image/acquisition joins | One unique exact metadata filename match per CSV row, matching specimen/species/process; every line locator resolves to that row; 22 existing glass records match, 31 unavailable glass IDs stay blank | CSV acquisition_record_location; master_image_metadata.csv; slide_glass_mapping.csv |
| Scan date separate from staining | Available scan dates are 20250516, 20250521, 20250522; standard records show 20250521. None promoted to staining_date | Separate scan_date/staining_date columns and negative control |
| Corrected supporting extraction | 51 exact filename metadata joins, 12 separately identified alternate acquisitions, matching retained IOD values, and exact author quote | online_lab_work_extraction.json at final pinned revision |
| Workbook evidence accurately bounded | 772 selected inventory cells across 51 ranges match value/date conversion, number format, style and fill; counts/headers/colors/hyperlinks match; no experimental-event terms or formulas found in populated cells | Salamander_slides_Lab_Work_source_snapshot.xlsx and extraction |
| Inventory gaps retained | Specimen 3463 / Process_331_raw_green.ome.tiff has no online inventory range; physical glass identities remain missing for 31 map rows, including the four denominator images | Explicit missing inventory record and blank CSV fields |
| Dataset link extraction | All 4 literal URLs and their titles/cell locations match; inspected headers match; no corresponding microscopy/staining/target terms or Excel hyperlink objects found | Desmognathus_dataset_source_snapshot.xlsx; ONLINE_DATASET_LINKS.json |
| Source custody | All 5 pinned snapshots hash correctly and each is byte-identical to its still-available original browser export | ONLINE_SOURCE_PINS.json; original export checks below |
| Inputs and pinned preparation preserved | All 5 combined input hashes match before/after checks; 136 preparation/evidence pin pairs and 614 protected figure hashes match retained pins | ONLINE_SOURCE_PINS.json; extraction input_pins; submission preparation FINAL_PINS.json and baseline.json |
| Calibration claims remain bounded | calibration_validated false; standard role says QC unresolved; row notes retain unresolved QC/calibration; source manifest retains conditional_assembly_derived_not_measured_c_value and unchanged target denominator | CSV, ONLINE_SOURCE_PINS.json, iod_source_manifest.json |
| Historical notes and current summaries | LOCAL_SEARCH.md retains conanti/???/“Not convinced it isn't monticola”/Formalin concerns; current ONLINE_LAB_WORK.md, SEARCH_STATUS.md and ONLINE_RECORD_RECOVERY.md say these remain unresolved and do not override the author assignment | Manually read current summaries; hashed source manifest and preserved preparation evidence |

The record resolves membership through author attestation. It does not establish a contemporaneous experimental event from workbook cells. Blank workbook staining_event/paired_standard fields remain appropriate; the separate author-attested CSV supplies the assignment.

## Failure sensitivity and corrected revision

Six in-memory rejection controls were detected: changed specimen, scan date used as staining date, invented physical glass, workbook-derived provenance, missing second standard pairing, and altered snapshot bytes. No on-disk input or scientific output was changed for these controls.

The initial strict read-only check exited 1 with “AssertionError: extracted image metadata join.” At the initial extraction revision (SHA256 a3f27b99a8d6b5653f5ee037005020f9831a805c6a99adeaee27e2d8a994f92f), 12 retained raw-image entries included both raw/brightfield and mounted metadata under local_image_metadata without filenames. The copied author answer also appended a period. The first continuation retained the failures rather than declaring a pass.

The original extraction owner corrected these evidence outputs. This verifier did not edit their implementation or data. The final strict rerun preserves the exact-filename equality check, additionally checks the newly explicit filenames and separately named alternate-acquisition records, and requires the exact author answer. It passes against extraction SHA256 f65d499ba9823b9f98b9c86cd76f3753123be0c3ad9cf154568077c631c6656f. The CSV hash stayed 085cee6d6a35ee4b7a1db0bd4058fde8f6664ea741199f5621b0805c9aa2f34d across this correction.

## What remains unproven

- Standards segmentation/remeasurement and QC, the reported approximately 41% standard-slide IOD disagreement, acquisition linearity and comparability, selection of a valid standard denominator, anchor uncertainty, and an independently validated absolute DNA-content scale were not tested. Scientific producers were deliberately not rerun. The manifest still names the target denominator and conditional assembly-derived anchor.
- The scientific truth of the recalled one-run event is author-attested, not independently observable from these exports. No actual staining date or original experimental run identifier was recovered.
- Online inventory specimen joins do not supply missing physical glass identities, complete preservation history, or a final resolution of historical taxonomic/accession discrepancies. Original archival workbooks were not reopened for a fresh semantic audit.
- Hashes and local byte identity verify retained custody and preservation; they do not independently prove the live browser account, cloud read-only history, or native document revision. Browser/account/read-only assertions are coordinator-reported. This verifier made no browser or cloud calls.
- The exports do not establish complete native Google Sheets comments/notes/formatting/revision history. Lab-work XLSX lacks docProps/core.xml; source-authored/revision timestamps remain unknown. Bounded cell/term searches do not prove there are no other records or indirect meanings.
- Protocol snapshots are protocol descriptions, not proof of actual slide settings. Current manuscript protocol wording still needs author reconciliation; its checksum is verified, not scientific correctness.
- This is record-recovery acceptance, not manuscript release, submission approval, calibrated genome-size acceptance, or proof that the run-record builder will behave correctly on arbitrary future inputs. Its current implementation was read but not executed.

No tool was broken during this verification. The initial failing command detected a real supporting-evidence ambiguity that was subsequently corrected by its owner. No missing required execution evidence remains for this bounded read-only check. The missing scientific evidence listed above keeps calibration/QC on HOLD.

## Current artifact and input pins

The full exact hashes, observed counts, negative-control messages, clock and interpreter are in the raw final output below. The two final summaries have these SHA256 values:

| Artifact | SHA256 |
|---|---|
| SEARCH_STATUS.md | 6a80a9059f8225c1698e239d667a3cb89914835ccfe2f42f15eff47a2aa03aeb |
| ONLINE_RECORD_RECOVERY.md | a40cd63a8910c47570d2f6c6a5201f0210f106f91c9d23d24694f7a0b6cf907c |

Summary local links resolved except links to this not-yet-written evidence file at the instant of the summary check. This report closes those two links; an explicit post-write check is retained at the end.

## Exact final verification command and observed result

Working directory: /home/jake/Projects/Desmognathus_TE. Project UV environment: /tmp/desmo-release-audit-20261013/env. No dependency synchronization or scientific producer execution.

```bash
UV_PROJECT_ENVIRONMENT=/tmp/desmo-release-audit-20261013/env PYTHONDONTWRITEBYTECODE=1 uv run --frozen --no-sync --python 3.13 python - <<'PY'
import copy,csv,datetime,hashlib,json,re,subprocess,sys,zipfile
from collections import Counter,defaultdict
from pathlib import Path
from openpyxl import load_workbook
start=datetime.datetime.now(datetime.timezone.utc).isoformat()
root=Path('outputs/calibration_search_20261013').resolve()
project=root.parent.parent
sha=lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
read=lambda p:json.loads(Path(p).read_text())
def csvrows(p):
    with Path(p).open(newline='') as f:return list(csv.DictReader(f))
def require(ok,label):
    if not ok:raise AssertionError(label)
rows=csvrows(root/'STAINING_RUN_MAP_AUTHOR_CONFIRMED.csv')
template=csvrows(project/'outputs/submission_preparation_20261013/STAINING_RUN_MAP_TO_COMPLETE.csv')
retained=csvrows(project/'analyses/02_morphology/data/nuclear_iod_by_image.csv')
attestation=read(root/'AUTHOR_STAINING_CONFIRMATION.json')
pins=read(root/'ONLINE_SOURCE_PINS.json')
extraction=read(root/'online_lab_work_extraction.json')
dataset=read(root/'ONLINE_DATASET_LINKS.json')
metadata_path=Path('/home/jake/Projects/cellprofiler_test/master_image_metadata.csv')
glass_path=metadata_path.with_name('slide_glass_mapping.csv')
metadata=csvrows(metadata_path);glasses=csvrows(glass_path)
metadata_index=defaultdict(list);glass_index=defaultdict(list)
for line,row in enumerate(metadata,2):metadata_index[row['filename']].append((line,row))
for line,row in enumerate(glasses,2):glass_index[row['process_id']].append((line,row))
protected=read(project/'outputs/submission_preparation_20261013/FINAL_PINS.json')
baseline=read(project/'outputs/submission_preparation_20261013/baseline.json')
artifact_names=['ONLINE_LAB_WORK.md','online_lab_work_extraction.json','ONLINE_DATASET_LINKS.md','ONLINE_DATASET_LINKS.json','AUTHOR_STAINING_CONFIRMATION.json','STAINING_RUN_MAP_AUTHOR_CONFIRMED.csv','ONLINE_SOURCE_PINS.json','build_author_run_record.py']
artifact_hashes={name:sha(root/name) for name in artifact_names}
input_pins=dict(extraction['input_pins'])
for path,digest in pins['local_inputs_unchanged'].items():
    require(path not in input_pins or input_pins[path]==digest,'overlapping input pins disagree')
    input_pins[path]=digest
require(all(sha(path)==digest for path,digest in input_pins.items()),'input hash drift')
require(sha(root/'STAINING_RUN_MAP_AUTHOR_CONFIRMED.csv')==pins['mapping_sha256'],'mapping pin mismatch')
require(sha(root/'AUTHOR_STAINING_CONFIRMATION.json')==pins['author_confirmation_sha256'],'attestation pin mismatch')
for snap in pins['read_only_source_snapshots']:
    require(sha(root/snap['path'])==snap['sha256'],'snapshot pin mismatch: '+snap['path'])
require(sha(extraction['snapshot'])==extraction['sha256'],'extraction snapshot pin mismatch')
require(sha(root/'Desmognathus_dataset_source_snapshot.xlsx')==dataset['source']['sha256'],'dataset snapshot pin mismatch')
pin_pairs=list(protected['final_files_sha256'].items())+list(protected['evidence_sources_sha256'].items())
protected_failures=[path for path,digest in pin_pairs if not (project/path).is_file() or sha(project/path)!=digest]
require(not protected_failures,'submission preparation pin drift: '+repr(protected_failures))
figure_failures=[record['path'] for record in baseline['protected_figure_run'] if sha(project/record['path'])!=record['sha256']]
require(not figure_failures,'protected figure pin drift')
require(attestation['answer_verbatim']=='One staining run with both standards','author answer altered')
require(attestation['actual_staining_date'] is None and attestation['original_experimental_run_identifier'] is None,'experimental event details invented')
require(attestation['run_id_is_internal_label'] is True,'run label not declared internal')
require(attestation['scope']=='All retained brightfield slides and both designated reference slides','author scope altered')
require(attestation['standard_specimen_ids']==['32469','32470'],'author standard scope changed')
def verify_map(candidate):
    require(len(candidate)==53,'row count')
    require([tuple(r[k] for k in ('species','filename','specimen_id','role')) for r in candidate]==[tuple(r[k] for k in ('species','filename','specimen_id','role')) for r in template],'original identifiers or roles changed')
    require(len({r['filename'] for r in candidate})==53,'duplicate images')
    target=[r for r in candidate if r['role'] in {'target','current species denominator'}]
    standards=[r for r in candidate if r['role']=='designated standard; QC unresolved']
    require(len(target)==51 and len(standards)==2,'target/standard counts')
    require(Counter(r['filename'] for r in target)==Counter(r['filename'] for r in retained),'retained coverage')
    require({re.fullmatch(r'Process_(\d+).*',r['filename']).group(1) for r in candidate if r['role']=='current species denominator'}=={'314','315','320','321'},'denominator membership')
    require({(re.fullmatch(r'Process_(\d+).*',r['filename']).group(1),r['specimen_id']) for r in standards}=={('413','32469'),('414','32470')},'designated standards changed')
    for row in candidate:
        require(row['staining_run_id']==attestation['staining_run_id'],'run membership')
        require(row['staining_date']=='','scan/collection date promoted to staining date')
        require(row['standard_process_for_run']=='Process_413;Process_414','paired standards')
        require(row['record_location']=='AUTHOR_STAINING_CONFIRMATION.json' and row['confirmation_basis']=='author attestation','provenance assigned to workbook')
        require(row['confirmation_date']=='2026-10-13' and row['run_id_is_internal_label']=='true','internal run or confirmation date altered')
        require('internal label' in row['notes'] and 'Staining date unknown' in row['notes'] and 'QC and calibration remain unresolved' in row['notes'],'record claims validation')
        require(len(metadata_index[row['filename']])==1,'ambiguous image metadata join')
        line,record=metadata_index[row['filename']][0]
        require(record['specimen_id']==row['specimen_id'] and record['species']==row['species'],'image identity join')
        process=re.fullmatch(r'Process_(\d+).*',row['filename']).group(1)
        require(record['slide_id']==process,'process identity join')
        require(len(glass_index[process])<=1,'ambiguous glass join')
        glass_line,glass=glass_index[process][0] if glass_index[process] else (None,{})
        for dest,source in [('physical_slide_glass_id','slide_glass_id'),('scan_date','scan_date'),('original_vsi_filename','vsi_filename')]:
            require(row[dest]==glass.get(source,''),'wrong glass/acquisition field: '+dest)
        expected=str(metadata_path)+':'+str(line)+((';'+str(glass_path)+':'+str(glass_line)) if glass_line else '')
        require(row['acquisition_record_location']==expected,'wrong source line locator')
verify_map(rows)
negative=[]
for label,key,value in [
    ('wrong specimen','specimen_id','999999'),
    ('scan date as staining date','staining_date',rows[0]['scan_date']),
    ('invented physical glass','physical_slide_glass_id','9999'),
    ('workbook-derived provenance','confirmation_basis','workbook-derived'),
    ('missing standard pairing','standard_process_for_run','Process_413'),
]:
    changed=copy.deepcopy(rows);changed[0][key]=value
    try:verify_map(changed)
    except AssertionError as err:negative.append({'case':label,'detected_by':str(err)})
    else:raise AssertionError('negative not detected: '+label)
require(hashlib.sha256((root/pins['read_only_source_snapshots'][0]['path']).read_bytes()+b'changed').hexdigest()!=pins['read_only_source_snapshots'][0]['sha256'],'hash negative not detected')
negative.append({'case':'snapshot bytes changed in memory','detected_by':'SHA256 differs from snapshot pin'})
wb=load_workbook(root/'Salamander_slides_Lab_Work_source_snapshot.xlsx',data_only=False)
def normal(value):
    return value.isoformat() if isinstance(value,(datetime.date,datetime.datetime)) else value
selected_count=0
for loc,cells in extraction['selected_inventory_rows'].items():
    sheet=loc.split('!',1)[0]
    for saved in cells:
        cell=wb[sheet][saved['cell']]
        for key,observed in [('value',normal(cell.value)),('number_format',cell.number_format),('style_id',cell.style_id),('fill_type',cell.fill.fill_type),('fill_rgb',cell.fill.fgColor.rgb if cell.fill.fgColor.type=='rgb' else None)]:
            require(saved[key]==observed,'extracted cell mismatch: '+sheet+'!'+saved['cell']+' '+key)
        selected_count+=1
nonempty={};fills={};links={}
for saved in extraction['sheets']:
    ws=wb[saved['sheet']]
    cells=[c for line in ws for c in line if c.value is not None]
    nonempty[ws.title]=len(cells)
    require(len(cells)==saved['nonempty_cells'],'nonempty cell count')
    require(ws.calculate_dimension()==saved['dimension'],'worksheet dimension')
    require(all(normal(ws[coord].value)==v for coord,v in saved['headers'].items()),'exported header mismatch')
    observed_fills=[{'cell':c.coordinate,'value':normal(c.value),'fill':c.fill.fgColor.rgb} for c in cells if c.fill.fill_type=='solid' and c.fill.fgColor.type=='rgb' and c.fill.fgColor.rgb not in {'FFFFFFFF','00000000','00FFFFFF'}]
    require(observed_fills==saved['nonwhite_fills'],'colored-marker extraction mismatch')
    fills[ws.title]=len(observed_fills)
    observed_links=[{'cell':c.coordinate,'url':c.hyperlink.target} for c in cells if c.hyperlink]
    require(observed_links==saved['hyperlinks'],'link extraction mismatch')
    links[ws.title]=len(observed_links)
    require(not saved['comments'] and not any(c.comment for c in cells),'unexpected workbook comments')
    require(len(ws.conditional_formatting)==saved['conditional_formatting_count']==0,'conditional formatting')
    require(not any(c.data_type=='f' for c in cells),'unexpected workbook formulas')
    terms=r'feulgen|stain|batch|standard|hydrolysis|fixation|formalin|glass|process'
    require(not any(re.search(terms,str(c.value),re.I) for c in cells),'experimental event field found')
with zipfile.ZipFile(root/'Salamander_slides_Lab_Work_source_snapshot.xlsx') as z:
    require('docProps/core.xml' not in z.namelist(),'source dates might exist in package')
    require(not any(re.search(r'(?:threadedComments|comments\d*\.xml)',name,re.I) for name in z.namelist()),'comment members found')
extracted_by_filename={r['filename']:r for r in extraction['retained_images']}
require(len(extracted_by_filename)==len(retained)==51,'extracted retained coverage')
missing=[]
for row in retained:
    record=extracted_by_filename[row['filename']]
    require(all(record[k]==v for k,v in row.items()),'retained IOD fields changed in extraction')
    require(record['staining_event'] is None and record['paired_standard'] is None,'workbook assigned experimental membership')
    expected_meta=[{k:m[k] for k in ['filename','slide_id','specimen_id','species','preservation','mounting','image_type']} for _,m in metadata_index[row['filename']]]
    process=re.fullmatch(r'Process_(\d+).*',row['filename']).group(1)
    require(record['local_image_metadata']==expected_meta,'extracted image metadata join')
    expected_alternates=[{k:m[k] for k in ['filename','slide_id','specimen_id','species','preservation','mounting','image_type']} for m in metadata if m['slide_id']==process and m['filename']!=row['filename']]
    require(record['local_image_metadata_join']=='exact filename','extraction join label')
    require(record['alternate_acquisitions_same_process']==expected_alternates,'alternate acquisition provenance')
    require(record['local_glass_mapping']==[g for _,g in glass_index[process]],'extracted glass join')
    loc=record['inventory_range']
    if loc:
        sheet,extent=loc.split('!',1);line=re.search(r'\d+',extent).group()
        id_cell=('D' if sheet=='Sheet1' else 'A')+line
        identity=str(wb[sheet][id_cell].value)
        if identity.startswith('RAP'):identity=identity[3:]
        if identity.endswith('.0'):identity=identity[:-2]
        require(identity==row['specimen_id'],'inventory identity join')
    else:missing.append((row['filename'],row['specimen_id']))
require(missing==[('Process_331_raw_green.ome.tiff','3463')],'unexpected missing specimen')
for standard in extraction['designated_standards']:
    require((standard['process'],standard['specimen']) in [('413','32469'),('414','32470')],'standard inventory linkage')
    sheet,extent=standard['range'].split('!',1);line=re.search(r'\d+',extent).group()
    require(str(int(wb[sheet]['D'+line].value))==standard['specimen'],'standard specimen inventory join')
dw=load_workbook(root/'Desmognathus_dataset_source_snapshot.xlsx',data_only=False)
require(dw.sheetnames==['Sheet1','List of 38 Desmognathus species','dnaPipeTE results','Links'],'dataset tabs')
urls=[];hyperlinks=[];matches=[]
terms=r'microscopy|blood|slide|stain|feulgen|batch|metadata|pixels?|picograms?|specimen|32469|32470|process|413414'
for ws in dw:
    for line in ws:
        for c in line:
            if c.hyperlink:hyperlinks.append((ws.title,c.coordinate))
            if c.value is not None:
                if re.search(r'https?://',str(c.value)):urls.append((ws.title,c.coordinate,str(c.value)))
                if re.search(terms,str(c.value),re.I):matches.append((ws.title,c.coordinate,str(c.value)))
require(urls==[('Links',r['url_cell'].split('!',1)[1],r['url']) for r in dataset['links']],'literal dataset URL coverage')
require(not matches and not hyperlinks,'dataset link/record absence disagrees')
for saved in dataset['links']:
    require(dw['Links'][saved['title_cell'].split('!',1)[1]].value==saved['title'],'dataset title cell')
for saved in dataset['inspected_headers']:
    sheet,coord=saved['cell'].split('!',1)
    require(dw[sheet][coord].value==saved['value'],'dataset header')
require(dataset['relevant_record_search']['matched_cells']==[] and dataset['links_fetched'] is False,'dataset falsely reports fetched records')
require(extraction['author_attestation_received_after_extraction']['exact_answer']==attestation['answer_verbatim'],'extraction author quote changed')
require(pins['calibration_validated'] is False and pins['cloud_writes_performed'] is False,'false calibration or cloud write assertion')
require(len(rows)==pins['target_image_rows']+pins['designated_standard_rows'],'receipt row count')
require(sum(int(r['n_nuclei']) for r in retained)==pins['retained_nuclei']==805,'receipt nuclei')
require(sum(bool(r['physical_slide_glass_id']) for r in rows)==pins['images_with_glass_records']==22,'receipt glass count')
manifest=read(project/'analyses/02_morphology/provenance/iod_source_manifest.json')
anchor=manifest['conditional_genome_size_anchor']
require(anchor['status']=='conditional_assembly_derived_not_measured_c_value' and anchor['anchor_uncertainty_included'] is False,'conditional calibration claim changed')
require('standards 413 and 414 are not yet the denominator' in anchor['current_release_denominator'],'denominator silently changed')
require(all(sha(path)==digest for path,digest in input_pins.items()),'inputs changed during verification')
require(artifact_hashes=={name:sha(root/name) for name in artifact_names},'artifacts changed during verification')
builder=(root/'build_author_run_record.py').read_text()
require('workbook-derived' not in builder,'builder claims workbook-derived provenance')
result={
 'started_utc':start,'ended_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
 'python':sys.version.split()[0],'git_head':subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip(),
 'artifact_sha256':artifact_hashes,'checked_input_sha256':input_pins,
 'snapshot_pins_verified':len(pins['read_only_source_snapshots']),
 'preparation_pin_pairs_verified':len(pin_pairs),'protected_figure_pins_verified':len(baseline['protected_figure_run']),
 'map_rows':len(rows),'roles':dict(Counter(r['role'] for r in rows)),
 'retained_images':len(retained),'retained_specimens':len({r['specimen_id'] for r in retained}),'retained_nuclei':805,
 'glass_rows':sum(bool(r['physical_slide_glass_id']) for r in rows),'missing_glass_rows':sum(not r['physical_slide_glass_id'] for r in rows),
 'scan_dates':sorted({r['scan_date'] for r in rows if r['scan_date']}),'staining_dates':sorted({r['staining_date'] for r in rows}),
 'selected_inventory_cells_verified':selected_count,'selected_inventory_ranges_verified':len(extraction['selected_inventory_rows']),
 'lab_nonempty_cells':nonempty,'colored_marker_counts':fills,'lab_hyperlink_counts':links,
 'dataset_literal_urls_verified':len(urls),'inventory_missing':missing,'negative_controls':negative,
 'author_quote_matches_extraction':extraction['author_attestation_received_after_extraction']['exact_answer']==attestation['answer_verbatim'],
 'extraction_exact_filename_joins_verified':len(extraction['retained_images']),'alternate_acquisitions_verified':sum(len(r['alternate_acquisitions_same_process']) for r in extraction['retained_images']),'calibration_validated':False,'all_assertions_passed':True
}
print(json.dumps(result,indent=2))
PY
```

Observed exit status: 0; tool wall time: 0.408090481 seconds.

```json
{
  "started_utc": "2026-10-13T17:12:27.173260+00:00",
  "ended_utc": "2026-10-13T17:12:27.449764+00:00",
  "python": "3.13.11",
  "git_head": "15dfb3c18755cc8659b626fa213087c62c6db02f",
  "artifact_sha256": {
    "ONLINE_LAB_WORK.md": "5f24b181e3b9c17e3d882a9e7e177abbc5f4897bf1d3b3b54d61d57a7527d63d",
    "online_lab_work_extraction.json": "f65d499ba9823b9f98b9c86cd76f3753123be0c3ad9cf154568077c631c6656f",
    "ONLINE_DATASET_LINKS.md": "8e31e8dd2d19b0fe94c5ceb5b434879f0b52774dce40181a32cd5adc7ae5751e",
    "ONLINE_DATASET_LINKS.json": "5bb85d34d7f41a981e29771ea2d5a1a2728e9a67615eccb821f5f169cc4f8e6e",
    "AUTHOR_STAINING_CONFIRMATION.json": "14c82b6d07e272d7ba120bf8cca62cdc2422b32f017e0a75f82878dfc192efd5",
    "STAINING_RUN_MAP_AUTHOR_CONFIRMED.csv": "085cee6d6a35ee4b7a1db0bd4058fde8f6664ea741199f5621b0805c9aa2f34d",
    "ONLINE_SOURCE_PINS.json": "14571e1cee158f30ac348426a640f90563c7c51356f33fa0940e2a5e5520f773",
    "build_author_run_record.py": "4ae6269296da19773ed408ce30ee93c588497405fcb1ff68658cda301ecb2f87"
  },
  "checked_input_sha256": {
    "/home/jake/Projects/Desmognathus_TE/analyses/02_morphology/data/nuclear_iod_by_image.csv": "577cc0d074715a91a2304320b134d474c75bdd2f4e04b5d93ce321b0c3e32033",
    "/home/jake/Projects/cellprofiler_test/master_image_metadata.csv": "bb164b21dfd5226a5e35879c4095ceb78474837b1510b05354fd165de4593629",
    "/home/jake/Projects/cellprofiler_test/slide_glass_mapping.csv": "aef95aa3cb3a964cb19ff86a6bc5112f1057c747466ee220397478b238617d99",
    "/home/jake/Projects/Desmognathus_TE/analyses/02_morphology/provenance/iod_source_manifest.json": "bb692d5cea7c84e50a97164c5a05476853546cad1717ba2edf9c1d0838067c35",
    "/home/jake/Projects/Desmognathus_TE/outputs/submission_preparation_20261013/STAINING_RUN_MAP_TO_COMPLETE.csv": "1fbdafd20d983e9bd201b9a0c57472b8ffeb0e407caaf45740139c00ae6dcf53"
  },
  "snapshot_pins_verified": 5,
  "preparation_pin_pairs_verified": 136,
  "protected_figure_pins_verified": 614,
  "map_rows": 53,
  "roles": {
    "target": 47,
    "current species denominator": 4,
    "designated standard; QC unresolved": 2
  },
  "retained_images": 51,
  "retained_specimens": 50,
  "retained_nuclei": 805,
  "glass_rows": 22,
  "missing_glass_rows": 31,
  "scan_dates": [
    "20250516",
    "20250521",
    "20250522"
  ],
  "staining_dates": [
    ""
  ],
  "selected_inventory_cells_verified": 772,
  "selected_inventory_ranges_verified": 51,
  "lab_nonempty_cells": {
    "Sheet1": 1071,
    "Sheet2": 64,
    "Sheet1 (1)": 769
  },
  "colored_marker_counts": {
    "Sheet1": 0,
    "Sheet2": 1,
    "Sheet1 (1)": 29
  },
  "lab_hyperlink_counts": {
    "Sheet1": 0,
    "Sheet2": 11,
    "Sheet1 (1)": 0
  },
  "dataset_literal_urls_verified": 4,
  "inventory_missing": [
    [
      "Process_331_raw_green.ome.tiff",
      "3463"
    ]
  ],
  "negative_controls": [
    {
      "case": "wrong specimen",
      "detected_by": "original identifiers or roles changed"
    },
    {
      "case": "scan date as staining date",
      "detected_by": "scan/collection date promoted to staining date"
    },
    {
      "case": "invented physical glass",
      "detected_by": "wrong glass/acquisition field: physical_slide_glass_id"
    },
    {
      "case": "workbook-derived provenance",
      "detected_by": "provenance assigned to workbook"
    },
    {
      "case": "missing standard pairing",
      "detected_by": "paired standards"
    },
    {
      "case": "snapshot bytes changed in memory",
      "detected_by": "SHA256 differs from snapshot pin"
    }
  ],
  "author_quote_matches_extraction": true,
  "extraction_exact_filename_joins_verified": 51,
  "alternate_acquisitions_verified": 12,
  "calibration_validated": false,
  "all_assertions_passed": true
}

```

## Exact summary and export identity check

```bash
UV_PROJECT_ENVIRONMENT=/tmp/desmo-release-audit-20261013/env PYTHONDONTWRITEBYTECODE=1 uv run --frozen --no-sync --python 3.13 python - <<'PY'
import ast,csv,datetime,hashlib,json,re
from pathlib import Path
root=Path('outputs/calibration_search_20261013')
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
pins=json.loads((root/'ONLINE_SOURCE_PINS.json').read_text())
tree=ast.parse((root/'build_author_run_record.py').read_text())
exports=[]
for node in tree.body:
    if isinstance(node,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='exports' for t in node.targets):exports=ast.literal_eval(node.value)
exports.append(('/tmp/browser-use/exports/Salamander slides Lab Work-fa268342-88bc-4753-be8a-3d71dc12961a.xlsx','Salamander_slides_Lab_Work_source_snapshot.xlsx','https://docs.google.com/spreadsheets/d/1qmsaXn9pffPq68fphp9Wz5Tcvup43oas2AN8DKQh8QY/edit'))
byte_identity=[]
for path,name,url in exports:
    present=Path(path).is_file()
    match=present and sha(path)==sha(root/name)
    if present:assert match,(path,'export snapshot differs')
    byte_identity.append({'retained':name,'original_export_exists':present,'bytes_match':match})
summary_hashes={}
summary_links=[]
for name in ['SEARCH_STATUS.md','ONLINE_RECORD_RECOVERY.md']:
    text=(root/name).read_text();summary_hashes[name]=sha(root/name)
    assert 'One staining run with both standards' in text
    assert 'internal label' in text
    assert 'remain' in text and ('unknown' in text or 'unresolved' in text)
    assert 'Scan dates' in text
    assert '314' in text and '315' in text and '320' in text and '321' in text
    assert '413' in text and '414' in text
    assert 'identity/preservation' in text
    assert 'calibration' in text
    for target in re.findall(r'\]\(([^)]+)\)',text):
        if target.startswith('https://'):
            assert target in {s['source_url'] for s in pins['read_only_source_snapshots']},target
        else:
            exists=(root/target).is_file()
            assert exists or target=='ONLINE_RECOVERY_VERIFICATION.md',target
            summary_links.append({'summary':name,'target':target,'exists_at_check':exists})
with (root/'STAINING_RUN_MAP_AUTHOR_CONFIRMED.csv').open(newline='') as f:rows=list(csv.DictReader(f))
selected=[{k:r[k] for k in ['filename','specimen_id','role','staining_run_id','staining_date','physical_slide_glass_id','scan_date','acquisition_record_location']} for r in rows if r['role']!='target']
print(json.dumps({'checked_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'summary_sha256':summary_hashes,'summary_local_links':summary_links,'export_snapshot_byte_identity':byte_identity,'denominator_and_standard_rows':selected},indent=2))
PY
```

Observed exit status: 0.

```json
{
  "checked_utc": "2026-10-13T17:13:39.671811+00:00",
  "summary_sha256": {
    "SEARCH_STATUS.md": "6a80a9059f8225c1698e239d667a3cb89914835ccfe2f42f15eff47a2aa03aeb",
    "ONLINE_RECORD_RECOVERY.md": "a40cd63a8910c47570d2f6c6a5201f0210f106f91c9d23d24694f7a0b6cf907c"
  },
  "summary_local_links": [
    {
      "summary": "SEARCH_STATUS.md",
      "target": "STAINING_RUN_MAP_AUTHOR_CONFIRMED.csv",
      "exists_at_check": true
    },
    {
      "summary": "SEARCH_STATUS.md",
      "target": "ONLINE_SOURCE_PINS.json",
      "exists_at_check": true
    },
    {
      "summary": "SEARCH_STATUS.md",
      "target": "ONLINE_LAB_WORK.md",
      "exists_at_check": true
    },
    {
      "summary": "SEARCH_STATUS.md",
      "target": "ONLINE_DATASET_LINKS.md",
      "exists_at_check": true
    },
    {
      "summary": "SEARCH_STATUS.md",
      "target": "AUTHOR_STAINING_CONFIRMATION.json",
      "exists_at_check": true
    },
    {
      "summary": "SEARCH_STATUS.md",
      "target": "ONLINE_RECORD_RECOVERY.md",
      "exists_at_check": true
    },
    {
      "summary": "SEARCH_STATUS.md",
      "target": "ONLINE_RECOVERY_VERIFICATION.md",
      "exists_at_check": false
    },
    {
      "summary": "ONLINE_RECORD_RECOVERY.md",
      "target": "ONLINE_LAB_WORK.md",
      "exists_at_check": true
    },
    {
      "summary": "ONLINE_RECORD_RECOVERY.md",
      "target": "ONLINE_DATASET_LINKS.md",
      "exists_at_check": true
    },
    {
      "summary": "ONLINE_RECORD_RECOVERY.md",
      "target": "Desmognathus_GBE_manuscript_source_snapshot.md",
      "exists_at_check": true
    },
    {
      "summary": "ONLINE_RECORD_RECOVERY.md",
      "target": "Cell_Stain_and_image_methodology_source_snapshot.md",
      "exists_at_check": true
    },
    {
      "summary": "ONLINE_RECORD_RECOVERY.md",
      "target": "Staining_Procedure_source_snapshot.md",
      "exists_at_check": true
    },
    {
      "summary": "ONLINE_RECORD_RECOVERY.md",
      "target": "ONLINE_SOURCE_PINS.json",
      "exists_at_check": true
    },
    {
      "summary": "ONLINE_RECORD_RECOVERY.md",
      "target": "AUTHOR_STAINING_CONFIRMATION.json",
      "exists_at_check": true
    },
    {
      "summary": "ONLINE_RECORD_RECOVERY.md",
      "target": "STAINING_RUN_MAP_AUTHOR_CONFIRMED.csv",
      "exists_at_check": true
    },
    {
      "summary": "ONLINE_RECORD_RECOVERY.md",
      "target": "ONLINE_RECOVERY_VERIFICATION.md",
      "exists_at_check": false
    }
  ],
  "export_snapshot_byte_identity": [
    {
      "retained": "Desmognathus_dataset_source_snapshot.xlsx",
      "original_export_exists": true,
      "bytes_match": true
    },
    {
      "retained": "Desmognathus_GBE_manuscript_source_snapshot.md",
      "original_export_exists": true,
      "bytes_match": true
    },
    {
      "retained": "Cell_Stain_and_image_methodology_source_snapshot.md",
      "original_export_exists": true,
      "bytes_match": true
    },
    {
      "retained": "Staining_Procedure_source_snapshot.md",
      "original_export_exists": true,
      "bytes_match": true
    },
    {
      "retained": "Salamander_slides_Lab_Work_source_snapshot.xlsx",
      "original_export_exists": true,
      "bytes_match": true
    }
  ],
  "denominator_and_standard_rows": [
    {
      "filename": "Process_314_raw_green.ome.tiff",
      "specimen_id": "3372",
      "role": "current species denominator",
      "staining_run_id": "author_confirmed_single_run",
      "staining_date": "",
      "physical_slide_glass_id": "",
      "scan_date": "",
      "acquisition_record_location": "/home/jake/Projects/cellprofiler_test/master_image_metadata.csv:4"
    },
    {
      "filename": "Process_315_raw_green.ome.tiff",
      "specimen_id": "3373",
      "role": "current species denominator",
      "staining_run_id": "author_confirmed_single_run",
      "staining_date": "",
      "physical_slide_glass_id": "",
      "scan_date": "",
      "acquisition_record_location": "/home/jake/Projects/cellprofiler_test/master_image_metadata.csv:5"
    },
    {
      "filename": "Process_320_raw_green.ome.tiff",
      "specimen_id": "3401",
      "role": "current species denominator",
      "staining_run_id": "author_confirmed_single_run",
      "staining_date": "",
      "physical_slide_glass_id": "",
      "scan_date": "",
      "acquisition_record_location": "/home/jake/Projects/cellprofiler_test/master_image_metadata.csv:10"
    },
    {
      "filename": "Process_321_raw_green.ome.tiff",
      "specimen_id": "3402",
      "role": "current species denominator",
      "staining_run_id": "author_confirmed_single_run",
      "staining_date": "",
      "physical_slide_glass_id": "",
      "scan_date": "",
      "acquisition_record_location": "/home/jake/Projects/cellprofiler_test/master_image_metadata.csv:11"
    },
    {
      "filename": "Process_413_green.ome.tiff",
      "specimen_id": "32469",
      "role": "designated standard; QC unresolved",
      "staining_run_id": "author_confirmed_single_run",
      "staining_date": "",
      "physical_slide_glass_id": "1321",
      "scan_date": "20250521",
      "acquisition_record_location": "/home/jake/Projects/cellprofiler_test/master_image_metadata.csv:123;/home/jake/Projects/cellprofiler_test/slide_glass_mapping.csv:83"
    },
    {
      "filename": "Process_414_green.ome.tiff",
      "specimen_id": "32470",
      "role": "designated standard; QC unresolved",
      "staining_run_id": "author_confirmed_single_run",
      "staining_date": "",
      "physical_slide_glass_id": "1322",
      "scan_date": "20250521",
      "acquisition_record_location": "/home/jake/Projects/cellprofiler_test/master_image_metadata.csv:124;/home/jake/Projects/cellprofiler_test/slide_glass_mapping.csv:84"
    }
  ]
}

```

## Retained failure observation

Initial strict check: exit 1, tool wall time 0.315714891 seconds.

```text
Traceback (most recent call last):
  File "<stdin>", line 137, in <module>
  File "<stdin>", line 13, in require
AssertionError: extracted image metadata join

```

The continuation kept this exact-image mismatch as a failure and completed the remaining checks; exit 0, all_assertions_passed=false. It found 12 affected raw-image entries: Processes 354–359 and 362–367, and author_quote_matches_extraction=false. Its extraction hash was a3f27b99a8d6b5653f5ee037005020f9831a805c6a99adeaee27e2d8a994f92f. This superseded failure is not an acceptance result. The full final strict command above tests the corrected revision without relaxing exact-image matching.

## Final revision and evidence-link check

```bash
UV_PROJECT_ENVIRONMENT=/tmp/desmo-release-audit-20261013/env PYTHONDONTWRITEBYTECODE=1 uv run --frozen --no-sync --python 3.13 python - <<'PY'
import datetime,hashlib,json
from pathlib import Path
root=Path('outputs/calibration_search_20261013')
expected={"ONLINE_LAB_WORK.md":"5f24b181e3b9c17e3d882a9e7e177abbc5f4897bf1d3b3b54d61d57a7527d63d","online_lab_work_extraction.json":"f65d499ba9823b9f98b9c86cd76f3753123be0c3ad9cf154568077c631c6656f","ONLINE_DATASET_LINKS.md":"8e31e8dd2d19b0fe94c5ceb5b434879f0b52774dce40181a32cd5adc7ae5751e","ONLINE_DATASET_LINKS.json":"5bb85d34d7f41a981e29771ea2d5a1a2728e9a67615eccb821f5f169cc4f8e6e","AUTHOR_STAINING_CONFIRMATION.json":"14c82b6d07e272d7ba120bf8cca62cdc2422b32f017e0a75f82878dfc192efd5","STAINING_RUN_MAP_AUTHOR_CONFIRMED.csv":"085cee6d6a35ee4b7a1db0bd4058fde8f6664ea741199f5621b0805c9aa2f34d","ONLINE_SOURCE_PINS.json":"14571e1cee158f30ac348426a640f90563c7c51356f33fa0940e2a5e5520f773","build_author_run_record.py":"4ae6269296da19773ed408ce30ee93c588497405fcb1ff68658cda301ecb2f87","SEARCH_STATUS.md":"6a80a9059f8225c1698e239d667a3cb89914835ccfe2f42f15eff47a2aa03aeb","ONLINE_RECORD_RECOVERY.md":"a40cd63a8910c47570d2f6c6a5201f0210f106f91c9d23d24694f7a0b6cf907c"}
assert all(hashlib.sha256((root/name).read_bytes()).hexdigest()==digest for name,digest in expected.items()),'reviewed artifact revision changed'
report=root/'ONLINE_RECOVERY_VERIFICATION.md'
assert report.is_file() and 'Verdict: **PASS' in report.read_text()
print(json.dumps({'checked_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'reviewed_artifact_pins_match':len(expected),'evidence_file_exists':report.is_file(),'summary_links_to_verification_resolved':True,'report_bytes':report.stat().st_size},indent=2))
PY
```

Observed exit status: 0.

```json
{
  "checked_utc": "2026-10-13T17:16:17.554844+00:00",
  "reviewed_artifact_pins_match": 10,
  "evidence_file_exists": true,
  "summary_links_to_verification_resolved": true,
  "report_bytes": 40507
}

```


