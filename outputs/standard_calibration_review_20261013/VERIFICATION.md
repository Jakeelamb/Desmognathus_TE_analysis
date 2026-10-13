# Independent final verification — 2026-10-13

**PASS for retained-mask reconstruction, target metric parity, gallery/source consistency, repeat determinism and protected-artifact integrity. Absolute calibration remains unqualified.** The final decision correctly retains the released conditional scale and rejects recalibration from the two unreviewed candidate medians.

I inspected the implementation and actual final files, recomputed from retained source pixels, reran the native metrics and refusal controls, and checked the manifests against current bytes. No implementation, acceptance check, source image, mask, scientific table, figure, Methods document or frozen preparation file was edited. This report is my only assigned output. The reviewer is a separate native agent from the same model family; this is not cross-family review.

## Identified revision

| Checkout | Observed HEAD |
|---|---|
| Paper | `15dfb3c18755cc8659b626fa213087c62c6db02f` |
| CellProfiler sibling | `d7c37fdaf2869245964b6710ca76696306c2daba` |
| Paper archive | `7c235b808efc08117ffda246890e425e464fdb74` |

The working tree contains extensive preexisting changes. These HEADs alone do not identify the reviewed working artifacts; the hashes below and retained manifests do.

## Observed results

The independent full-pixel command below ran from `/home/jake/Projects/Desmognathus_TE`, from `2026-10-13T18:03:33.851279+00:00` to `2026-10-13T18:03:51.303862+00:00`, exit **0**.

| Check against actual final artifact | Observed result |
|---|---|
| Original standard CSV and inventory identity | 47 exact source records: Process 413=31, Process 414=16; original order preserved |
| First canonical target rows | Unchanged rows at CSV lines 313/325/333/338 for Processes 314/315/320/321 |
| Independently calculated whole-label IOD and area | All 51 areas exact; all 51 IOD values agree at six decimals; maximum unrounded difference from rounded CSV=4.895912297797622e-07 |
| Native metrics rerun | All 408 current metric values agree exactly with objects.json; 32 historical target values agree within 1e-12, maximum difference=9.020562075079397e-17 |
| Metric and display label coverage | All 51 complete recorded labels retained; all four historical target rows exercise the native full-tile fallback |
| Actual RAW/MASK/OVERLAY arrays | 153/153 saved PNG arrays equal independently extracted source crops, label masks and red inner boundaries; no contrast/intensity adjustment |
| Actual gallery membership/order | 51 articles, 153 unique images, original object order, no extra crop files |
| Consumed inputs/helpers | Both runs' 62 before/after pins agree and match current source SHA256; includes seven native helper modules |
| Actual repeat consistency | 157/157 numerical, control and gallery files are byte-identical; 51 object records numerically identical; receipt clocks/destinations are distinct |
| Six falsifications rerun | Absent label, unconverted global coordinates, wrong expected source hash, altered background, removed positive-OD pixel and changed stored IOD all rejected |
| Output guards rerun | Existing destination, symlink ancestor, ancestry overlap in both directions and absent destination parent all refused |
| Actual native CLI refusal | Exit 2; all 161 existing v1 files retained identical SHA256, size and mtime |
| Selected artwork/source baseline | All 614 `baseline.json.protected_figure_run` entries match current bytes |
| Frozen submission-preparation pins | All 58 `FINAL_PINS.json.final_files_sha256` and 78 `evidence_sources_sha256` values match current bytes |
| Selected submission inventory | All 34 rows in the actual category/item/path/sha256/disposition schema match current files |
| Final worker evidence manifest | All 333 list entries with path/bytes/sha256 match; root integration reports outside that manifest are separately pinned below |

The source arrays are uint8 grayscale with uint32 labels: 50 rows use 4096×4096 arrays and one target row uses 4096×1182 arrays. This establishes the retained arrays' depth; it does not recover original camera bit depth.

The falsifications were actual calls on copied rows/masks, not only assertions about retained PASS flags. On the first genuine standard object, background 255→254 changes IOD by -5.09208835335005; removing the positive-OD pixel at [2804,1760] changes area by -1 and IOD by -0.1512680076535844. A stored-IOD change of 0.000010 fails six-decimal equality. The hash control changes the expected digest, leaving source bytes intact.

The separate scope/manifest command below ran from `2026-10-13T18:06:29.059384+00:00` to `2026-10-13T18:06:29.119514+00:00`, exit **0**. The maintained read-only morphology validator also exited **0**:

```bash
UV_PROJECT_ENVIRONMENT=/tmp/desmo-standard-qc-20261013/env uv run --frozen --python 3.13 --with Pillow --with tifffile --with scipy --with scikit-image --with matplotlib --with roifile python analyses/02_morphology/validate_release.py
```

Its observed output was `PASS: morphology release replay agrees with frozen canonical tables; conditional 16.1-Gbp assembly anchor validated`. Here “validated” means agreement of the implemented conditional arithmetic with frozen tables; it does not mean independently validated DNA content.

## Author facts, denominator and frozen seam

The two author records retain these exact answers:

> yes these are the 40x brightfield images without oil. they are also the green channel images. they were fixed with the same protocol as the others

> One staining run with both standards

The live run-map schema distinguishes **47 target**, **4 current species denominator**, and **2 designated standard; QC unresolved** rows. Together the first two roles equal exactly the 51 released target filenames; the two standards are specimens 32469/32470. All 53 use `author_confirmed_single_run`, explicitly an internal label, with a blank staining date. This is author attestation, not a recovered contemporaneous date/log.

`AUTHOR_IMAGE_CONFIRMATION.json`, `METHODS_UPDATE_PROPOSED.md` and `CALIBRATION_DECISION.md` correctly limit the confirmation to 40× brightfield/no oil/green channel/shared fixation protocol and one staining run. They do not infer identical lifetime museum preservation, exact fixation execution, camera response or photometric linearity. `CURRENT_EVIDENCE_ADDENDUM.md` explicitly supersedes the older measurement-flow snapshot's missing run-map/acquisition-mode/environment/runtime-fixture statements; the earlier snapshot remains retained and must be read with this addendum.

The current 805-row nuclear table contains no Process 413/414. Its D. fuscus reference rows are:

| Target process | Reviewed nuclei | Actual image median IOD |
|---|---:|---:|
| 314 | 12 | 1079.688171 |
| 315 | 8 | 1061.3096540000001 |
| 320 | 5 | 833.276777 |
| 321 | 10 | 801.2819585 |

Their equal-image mean is **943.889140125**, matching the released species summary. Actual standard medians are **594.078400** and **837.701417**; the latter is **41.00856334786789%** above the former. No accepted standard denominator or causal explanation of that difference follows from these measurements.

I read the real `extend_finalized_microscopy_reviewers.py:apply_frozen_scaling/match_species_to_template` and `build_image_quality_matched_genome_iod_panel.py:reviewed_panel_freeze_ready`. The retained matched-panel SHA256 matches the canonical source manifest. It contains exactly 40 D. apalachicolae template IDs 1–40, 20 each from Processes 356/357. Reapplying the preserved scalar parameters reproduces the two stored z features within 1.1102230246251565e-16 and 9.540979117872439e-17. The native extension computes balance against assigned template slots when fewer than 40 candidates are supplied. Consequently subset common support does not prove full-template balance. I did not execute standard linkage or matching.

The historical ≥30 guard groups by **species**, not standard image. The frozen-interface report's proposed per-standard cardinality policy is a future equivalent-workflow diagnostic, not an established DNA-cytometry minimum; the final decision correctly states this distinction. Manual keep/problem/unsure decisions, linked strict eligibility, fixed-template projection and honest shortfall reporting remain missing for the standards.

The canonical anchor metadata remains `conditional_assembly_derived_not_measured_c_value`: 16.1/0.978=16.462167689161554 pg/1C equivalent, with anchor uncertainty excluded. I independently retrieved the [official NCBI dataset report](https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/GCA_050004315.1/dataset_report): GCA_050004315.1, aDesFus1-2.1, total sequence length 16,117,580,045 bp, pooled male/female AMNH A-194068–9, combining SAMN38201181/182. These are different animals from standard specimens 32469/32470. Inference: that assembly record establishes an assembled-length source and pooling provenance; it supplies no measured DNA-content value for these two slides. The selected rounded 16.1-Gbp anchor remains distinct from the exact accession span.

## Visual flags checked against actual files

I viewed the original-resolution source-matched overlays for ordinals 3, 17, 32, 41 and 45 and checked their identities/values directly against original CSV records and masks. The observations in `VISUAL_REVIEW.md` are supported as boundary-review flags, not accepted biological diagnoses.

| Ordinal | Exact source key | CSV line | Area / stored IOD | Independent observation |
|---:|---|---:|---|---|
| 3 | Process_413_green.ome.tiff::tile_y000000_x032768.tiff::2 | 4 | 6083 / 2081.880416 | Selected region surrounded by irregular crowded dark structures; isolation requires review |
| 17 | Process_413_green.ome.tiff::tile_y008192_x045056.tiff::4 | 18 | 3391 / 626.258628 | Mask bbox ends at y=4096 on the actual source edge |
| 32 | Process_414_green.ome.tiff::tile_y000000_x024576.tiff::1 | 33 | 696 / 140.478295 | Small polygon within a substantially larger visible structure |
| 41 | Process_414_green.ome.tiff::tile_y020480_x040960.tiff::1 | 42 | 593 / 160.420607 | Minimum mask y=1; display padding clips, mask does not touch y=0 |
| 45 | Process_414_green.ome.tiff::tile_y028672_x065536.tiff::1 | 46 | 1692 / 484.269010 | Nearly rectangular partial outline with surrounding overlap |

Ordinal 37 also clips display padding but its mask bbox [104,59,178,114] does not touch the tile edge. Among all 51 masks, only ordinal 17 is source-edge touching. This review does not substitute a complete manual accept/reject census for the standards.

## Reviewed file pins

Paths below are relative to this review directory. These identify the final reports I read, including root integration files outside the worker manifest.

| File | SHA256 |
|---|---|
| audit_mask_pixels.py | 18839f8f40a0059cd785bfad627f0cf6dea788deeb2f67572f32c253b7f19ded |
| AUTHOR_IMAGE_CONFIRMATION.json | 3c6057e4a4e85455376176a94f062a37ebae31250f3540d66ade595cbcc7e95a |
| METHODS_UPDATE_PROPOSED.md | efba84ae78e6962436bcd7644df7e404fb560577d152b4a41579df3aa6bb5645 |
| VISUAL_REVIEW.md | fe00385ad43d4f8d78d8fa1c3025eaec6c6efbd2d607f3569efae2c980f8db15 |
| CALIBRATION_DECISION.md | 6b7db9243a54692fec2aa69a4750f4748f892f75b32e60142303a5f7435793b3 |
| CURRENT_EVIDENCE_ADDENDUM.md | fba3a6ab150ca651c255e5206c4585a14c1cc852a5b4ce85cdcbd3c7e3026bd5 |
| FROZEN_STANDARD_QC_INTERFACE.md | 663678978880f7945deaae315c3ab596b1a8b2b2d625c17669ca7b2703fb8e97 |
| PIXEL_AUDIT_EXECUTION.md | d9e65c619684a02f4adbf40aaa5d54a1c2050af9ffcd2b082a5abfc00ef92180 |
| pixel_audit_v1/receipt.json | 1619e4077d98e1c1f36ac4f6737d8c7b508acf50ae914de4ae13de5f29c405b5 |
| pixel_audit_v2/receipt.json | 567178d07d002140c5f7727459b22cdad0ce1da79a80ac6df6698a12a49334ed |
| pixel_audit_v2/FINAL_EVIDENCE_MANIFEST.json | 51759d669cf1b3966b86311357c31dfbe9520ae0c8547e980e6c78f441188ba6 |

Full retained evidence paths: `pixel_audit_v1/objects.json`, `objects.csv`, `controls.json`, `gallery.html`, `receipt.json`, `command.txt/log`; corresponding v2 files; `pixel_audit_v2/validation.json`, `validation_command.txt`, `existing_output_refusal_receipt.json` and `FINAL_EVIDENCE_MANIFEST.json`. Protected manifests are `../submission_preparation_20261013/baseline.json`, `FINAL_PINS.json` and `SUBMISSION_INVENTORY.csv`. Their original schemas, rather than earlier report counts, drove the comparisons.

## Failed observations and remaining limits

- My first independent image-reader command exited **1** at `tifffile.imread` with `TiffFileError: not a TIFF file: header=b'\\x89PNG'`. It mistakenly treated a target PNG as TIFF. The corrected successful command below selects Pillow for PNG and tifffile for TIFF. No audit artifact was changed or regenerated. The worker's separately retained checker failure is in `pixel_audit_v2/validation_failed_v1.log`.
- My first run-map check exited **1** after the manifest/denominator checks because it incorrectly assumed all 51 unknowns had the literal role `target`; the actual schema distinguishes four `current species denominator` rows. The successful scope command below uses all three actual role labels. This was a verifier schema error, not missing run membership.
- Browser retrieval failed for NCBI API/BioC endpoints; publisher retrieval failed after a CDN redirect. Direct UV/urllib retrieval of the official NCBI API succeeded. The attempted article BioC fetch returned **HTTP 429**. I did not independently reproduce the anchor investigator's detailed article/Hardie/historical-original-text retrievals; those bounded provenance/access limits remain in `ANCHOR_EVIDENCE.md`. The dataset provenance and local conditional-anchor metadata were independently checked.
- Native metrics prove parity with the current helper behavior; this is not an independent validation of that helper's photometric model. The observed target full-tile fallback and current float32/percentile normalization are preserved.
- Exact retained-mask IOD can reproduce incomplete masks. No standard nucleus acceptance, repaired boundary, linked-pair eligibility, frozen technical match, scale-linearity measurement or independently validated specimen C-value was established. Original full VSI/ETS sources are unavailable at the recorded paths; original camera response/settings and lifetime preservation remain unresolved.
- Hash checks establish the listed files' present identity against retained baselines. They do not prove provenance of every unlisted file, upstream execution, or historical checkpoint-to-mask attribution.
- No new calibration number, model refit, segmentation, figure generation, HPC run or external edit occurred. The bounded PASS is not absolute-calibration, manuscript-submission or publication approval.

## Exact independent full-pixel command

```bash
UV_PROJECT_ENVIRONMENT=/tmp/desmo-standard-qc-20261013/env uv run --frozen --python 3.13 --with Pillow --with tifffile --with scipy --with scikit-image --with matplotlib --with roifile python - <<'PY'
import sys
sys.dont_write_bytecode=True
import csv,hashlib,html,importlib.util,json,re,subprocess
from collections import Counter
from datetime import datetime,timezone
from pathlib import Path
import numpy as np
import pandas as pd
import tifffile
from PIL import Image
from skimage.segmentation import find_boundaries
R=Path('outputs/standard_calibration_review_20261013').resolve(); P=Path.cwd(); C=Path('/home/jake/Projects/cellprofiler_test')
sys.path.insert(0,str(C/'scripts'))
import audit_iod_measurement_quality as quality
import apply_pair_review_repairs as repair
import run_linked_pair_review as linked
metrics=('iod_inner_intensity_median','iod_ring_intensity_median','iod_local_contrast','iod_local_od_contrast','iod_boundary_grad_mean','iod_edge_sharpness','iod_ring_noise_iqr','iod_contrast_snr')
def sha(p):
 with Path(p).open('rb') as f: return hashlib.file_digest(f,'sha256').hexdigest()
def rows(p):
 with Path(p).open(newline='') as f: return [(i,r) for i,r in enumerate(csv.DictReader(f),2)]
def load(p): return json.loads(Path(p).read_text())
def snapshot(p): return {str(f.relative_to(p)):(sha(f),f.stat().st_size,f.stat().st_mtime_ns) for f in p.rglob('*') if f.is_file()}
print('START',datetime.now(timezone.utc).isoformat(),flush=True)
print('HEADS',{str(p):subprocess.check_output(['git','-C',str(p),'rev-parse','HEAD'],text=True).strip() for p in [P,C,P.with_name('Desmognathus_TE_archive')]},flush=True)
objects=load(R/'pixel_audit_v2/objects.json'); std=rows(C/'output/runs/life_history_pmount_reference_v1/nucleus_measurements/all_measurements.csv'); target=rows(P/'analyses/02_morphology/data/nuclear_iod_objects.csv')
assert [{'line':i,**r} for i,r in std]==load(R/'source_inventory.json')['standard_yolo_candidate_rows']
expected=[('standard',i,r) for i,r in std]+[('target_control',i,r) for process in ['314','315','320','321'] for i,r in [next((i,r) for i,r in target if r['filename']==f'Process_{process}_raw_green.ome.tiff')]]
assert len(objects)==len(expected)==51
metric_max=iod_max=0.; images=[]; clipped=[]; fallback=[]; shapes=Counter()
for obj,(role,line,row) in zip(objects,expected,strict=True):
 assert obj['role']==role and obj['csv_line']==line and obj['source_record']==row
 projected=dict(row)
 if role=='standard':
  for old,new in [('label','nucleus_label'),('area_px','nuc_area_px'),('iod','nuc_iod'),('i_bg','nuc_i_bg'),('mask_path','nucleus_mask_path'),('source_image_path','nucleus_source_image_path')]: projected[new]=row[old]
  projected['nuc_centroid_x']=float(row['centroid_x'])-float(row['tile_x0']); projected['nuc_centroid_y']=float(row['centroid_y'])-float(row['tile_y0']); projected['tile_x0']=projected['tile_y0']=0
 assert projected==obj['metric_input_record']
 path=Path(projected['nucleus_source_image_path']); raw=np.asarray(Image.open(path)) if path.suffix.lower() in ['.png','.jpg','.jpeg'] else tifffile.imread(path)
 labels=tifffile.imread(projected['nucleus_mask_path']); assert raw.dtype==np.uint8
 gray=raw[...,0] if raw.ndim==3 and raw.shape[-1] in [3,4] else np.squeeze(raw)
 assert gray.shape==labels.shape and labels.ndim==2
 shapes[(str(raw.dtype),str(labels.dtype),str(raw.shape))]+=1
 mask=labels==int(projected['nucleus_label']); assert mask.any()
 area=int(np.count_nonzero(mask)); value=float(np.log10(max(float(projected['nuc_i_bg']),1)/np.maximum(gray[mask].astype(np.float64),1)).sum())
 assert area==int(float(projected['nuc_area_px']))==obj['mask_area_px']
 assert f'{value:.6f}'==f"{float(projected['nuc_iod']):.6f}"==obj['recomputed_iod_six_decimals']; assert value==obj['recomputed_iod']
 iod_max=max(iod_max,abs(value-float(projected['nuc_iod'])))
 if role=='standard':
  ys,xs=np.nonzero(mask); assert abs(xs.mean()-projected['nuc_centroid_x'])<=.02 and abs(ys.mean()-projected['nuc_centroid_y'])<=.02
 actual=quality.compute_row_iod_metrics(pd.Series(projected)); assert actual==obj['pixel_metrics']
 if role=='target_control':
  for key in metrics:
   delta=abs(actual[key]-float(row[key])); assert np.isfinite(delta) and delta<=1e-12; metric_max=max(metric_max,delta)
 crop_raw,crop_mask=quality.row_local_masks(pd.Series(projected)); assert int(crop_mask.sum())==area
 if obj['helper_full_tile_fallback']: fallback.append(obj['row_key'])
 ys,xs=np.nonzero(mask); bbox=[int(xs.min()),int(ys.min()),int(xs.max()+1),int(ys.max()+1)]; desired=[bbox[0]-64,bbox[1]-64,bbox[2]+64,bbox[3]+64]
 x0,y0,x1,y1=max(0,desired[0]),max(0,desired[1]),min(gray.shape[1],desired[2]),min(gray.shape[0],desired[3])
 assert bbox==obj['gallery_label_bbox_xyxy'] and [x0,y0,x1,y1]==obj['gallery_bounds_xyxy']
 rawcrop=gray[y0:y1,x0:x1]; maskcrop=mask[y0:y1,x0:x1]; assert int(maskcrop.sum())==area
 overlay=np.repeat(rawcrop[...,None],3,axis=-1); overlay[find_boundaries(maskcrop,mode='inner')]=[255,32,32]
 for name,array in zip(obj['gallery_images'],[rawcrop,maskcrop.astype(np.uint8)*255,overlay],strict=True): assert np.array_equal(np.asarray(Image.open(R/'pixel_audit_v2'/name)),array),name; images.append(name)
 assert sha(obj['raw_path'])==obj['raw_sha256'] and sha(obj['mask_path'])==obj['mask_sha256']
 if desired != [x0,y0,x1,y1]: clipped.append((obj['ordinal'],obj['row_key'],bbox[0]==0 or bbox[1]==0 or bbox[2]==gray.shape[1] or bbox[3]==gray.shape[0]))
 quality.load_row_tile_arrays.cache_clear(); linked.load_raw_array.cache_clear(); linked.load_label_array.cache_clear()
gallery=(R/'pixel_audit_v2/gallery.html').read_text(); assert [html.unescape(s) for s in re.findall(r'<img src="([^"]+)"',gallery)]==images
assert gallery.count('<article>')==51 and len(set(images))==153 and set(images)=={str(p.relative_to(R/'pixel_audit_v2')) for p in (R/'pixel_audit_v2/crops').iterdir()}
print('DIRECT',{'rows':51,'standards':dict(Counter(r['slide_id'] for _,r in std)),'target_lines':[i for role,i,r in expected if role=='target_control'],'area_mismatches':0,'six_decimal_iod_mismatches':0,'max_unrounded_iod_delta':iod_max,'native_metric_values':408,'historical_target_metric_values':32,'max_target_metric_delta':metric_max,'gallery_images_exact':len(images),'gallery_articles':51,'clipped':clipped,'target_full_tile_fallback':fallback,'source_shapes':dict(shapes)},flush=True)
a,b=load(R/'pixel_audit_v1/receipt.json'),load(R/'pixel_audit_v2/receipt.json')
for run,receipt in [('pixel_audit_v1',a),('pixel_audit_v2',b)]:
 assert receipt['exit_status']==0 and receipt['issues']==[] and receipt['inputs_before']==receipt['inputs_after']
 for path,pin in receipt['inputs_before'].items(): assert sha(path)==pin['sha256'],path
 for path,pin in receipt['artifacts'].items(): assert sha(R/run/path)==pin['sha256'],path
 assert receipt['torch_loaded'] is False
assert set(a['artifacts'])==set(b['artifacts'])
for name in a['artifacts']: assert (R/'pixel_audit_v1'/name).read_bytes()==(R/'pixel_audit_v2'/name).read_bytes(),name
assert load(R/'pixel_audit_v1/objects.json')==objects
print('REPEAT',{'source_pins_per_run':len(a['inputs_before']),'actual_identical_artifact_bytes':len(a['artifacts']),'identical_numerical_objects':len(objects),'receipts_are_distinct':(R/'pixel_audit_v1/receipt.json').read_bytes()!=(R/'pixel_audit_v2/receipt.json').read_bytes()},flush=True)
spec=importlib.util.spec_from_file_location('audit',R/'audit_mask_pixels.py'); audit=importlib.util.module_from_spec(spec); spec.loader.exec_module(audit)
for key,val in [('np',np),('pd',pd),('repair',repair),('quality',quality),('linked',linked)]: setattr(audit,key,val)
control=audit.controls(objects[0]); assert control==load(R/'pixel_audit_v2/controls.json') and control['all_rejected']
print('FALSIFICATION',[(x['name'],x['rejected']) for x in control['tests']],flush=True)
guards=[]
for name,out,inputs in [('existing',R/'pixel_audit_v1',()),('symlink',Path('/bin/desmo-verification-output'),()),('input_ancestor',R/'not_created',(R,)),('output_contains_input',R/'not_created',(R/'not_created/input.tiff',)),('missing_parent',R/'absent_parent/not_created',())]:
 try: audit.fresh_path(out,inputs)
 except ValueError as exc: guards.append((name,str(exc))); assert not out.exists() or name=='existing'
 else: raise AssertionError(name)
print('GUARDS',guards,flush=True)
before=snapshot(R/'pixel_audit_v1'); refusal=subprocess.run(['uv','run','--frozen','--python','3.13','--with','Pillow','--with','tifffile','--with','scipy','--with','scikit-image','--with','matplotlib','--with','roifile','python',str(R/'audit_mask_pixels.py'),'--out',str(R/'pixel_audit_v1')],capture_output=True,text=True)
assert refusal.returncode==2 and snapshot(R/'pixel_audit_v1')==before
print('CLI_REFUSAL',{'exit':refusal.returncode,'files_unchanged':len(before),'output':refusal.stdout+refusal.stderr},flush=True)
prep=P/'outputs/submission_preparation_20261013'; baseline=load(prep/'baseline.json'); final=load(prep/'FINAL_PINS.json')
for pin in baseline['protected_figure_run']: assert sha(P/pin['path'])==pin['sha256'],pin['path']
for field in ['final_files_sha256','evidence_sources_sha256']:
 for path,value in final[field].items(): assert sha(P/path)==value,(field,path)
selected=rows(prep/'SUBMISSION_INVENTORY.csv')
for _,row in selected: assert sha(P/row['path'])==row['sha256'],row['path']
print('PROTECTED',{'selected_figure_baseline_pins':len(baseline['protected_figure_run']),'final_files':len(final['final_files_sha256']),'evidence_sources':len(final['evidence_sources_sha256']),'selected_inventory':len(selected)},flush=True)
print('FINISH',datetime.now(timezone.utc).isoformat(),'EXIT 0',flush=True)
PY
```

Observed concluding output: `DIRECT rows=51, area_mismatches=0, six_decimal_iod_mismatches=0, native_metric_values=408, historical_target_metric_values=32, gallery_images_exact=153`; `REPEAT source_pins_per_run=62, actual_identical_artifact_bytes=157, identical_numerical_objects=51`; all six falsifications true; five guards refused; `CLI_REFUSAL exit=2, files_unchanged=161`; `PROTECTED selected_figure_baseline_pins=614, final_files=58, evidence_sources=78, selected_inventory=34`; `EXIT 0`.

## Exact independent scope/manifest command

```bash
UV_PROJECT_ENVIRONMENT=/tmp/desmo-standard-qc-20261013/env uv run --frozen --python 3.13 --with Pillow --with tifffile --with scipy --with scikit-image --with matplotlib --with roifile python - <<'PY'
import sys
sys.dont_write_bytecode=True
import json,hashlib
from pathlib import Path
from datetime import datetime,timezone
import pandas as pd
R=Path('outputs/standard_calibration_review_20261013'); A=Path('/home/jake/Projects/Desmognathus_TE_archive/path_analysis/data/external/derived/image_quality_matched_genome_iod')
def sha(p):
 with Path(p).open('rb') as f: return hashlib.file_digest(f,'sha256').hexdigest()
print('START',datetime.now(timezone.utc).isoformat())
manifest=json.loads((R/'pixel_audit_v2/FINAL_EVIDENCE_MANIFEST.json').read_text())
for x in manifest:
 f=R/x['path']; assert f.stat().st_size==x['bytes'] and sha(f)==x['sha256'],x['path']
print('PIXEL_FINAL_MANIFEST',len(manifest),'verified')
obj=pd.read_csv('analyses/02_morphology/data/nuclear_iod_objects.csv'); f=obj[obj.species=='D. fuscus']; assert not obj.filename.str.contains('Process_413|Process_414',regex=True).any()
assert set(f.filename)=={f'Process_{x}_raw_green.ome.tiff' for x in [314,315,320,321]}
by=f.groupby('filename').nuc_iod.agg(['count','median']); released=pd.read_csv('analyses/02_morphology/data/relative_nuclear_iod_species.csv').set_index('species')
assert abs(by['median'].mean()-released.loc['D. fuscus','iod_equal_image_estimate'])<1e-8
print('TARGET_FUSCUS',by.to_dict('index'),'EQUAL_IMAGE_IOD',by['median'].mean())
std=pd.read_csv('/home/jake/Projects/cellprofiler_test/output/runs/life_history_pmount_reference_v1/nucleus_measurements/all_measurements.csv'); med=std.groupby('slide_id').iod.median()
print('STANDARD_MEDIANS',med.to_dict(),'DIRECTIONAL_PERCENT',(med[414]/med[413]-1)*100)
run=pd.read_csv('outputs/calibration_search_20261013/STAINING_RUN_MAP_AUTHOR_CONFIRMED.csv',dtype=str,keep_default_na=False)
assert len(run)==53 and set(run.staining_run_id)=={'author_confirmed_single_run'} and set(run.run_id_is_internal_label)=={'true'} and set(run.staining_date)=={''}
unknown=run[run.role.isin(['target','current species denominator'])]; standards=run[run.role=='designated standard; QC unresolved']
assert len(unknown)==51 and set(unknown.filename)==set(obj.filename) and len(standards)==2 and set(standards.specimen_id)=={'32469','32470'}
assert set(run.standard_process_for_run)=={'Process_413;Process_414'}
a=json.loads((R/'AUTHOR_IMAGE_CONFIRMATION.json').read_text()); s=json.loads(Path('outputs/calibration_search_20261013/AUTHOR_STAINING_CONFIRMATION.json').read_text())
assert a['answer_verbatim']=='yes these are the 40x brightfield images without oil. they are also the green channel images. they were fixed with the same protocol as the others'
assert s['answer_verbatim']=='One staining run with both standards'
print('AUTHOR_MAP',run.role.value_counts().to_dict(),'VERBATIM_FACTS',a['answer_verbatim'],s['answer_verbatim'])
pin=json.loads(Path('analyses/02_morphology/provenance/iod_source_manifest.json').read_text()); archive=json.loads((A/'manifest.json').read_text()); path=A/'image_quality_matched_nuclei_all_species_sensitivity.csv.gz'
assert sha(path)==pin['matched_panel_sha256'] and archive['feature_scaling']==pin['feature_scaling']
template=pd.read_csv(path); template=template[template.species=='D. apalachicolae'].sort_values('quality_template_id')
assert len(template)==40 and template.quality_template_id.nunique()==40 and set(template.quality_template_id)==set(range(1,41))
deltas={}
for field,params in pin['feature_scaling'].items():
 value=(template[field].clip(params['clip_low'],params['clip_high'])-params['median'])/params['scale_iqr']; delta=float((value-template['z_'+field]).abs().max()); assert delta<1e-12; deltas[field]=delta
anchor=pin['conditional_genome_size_anchor']; assert anchor['reported_assembly_span_gbp']==16.1 and anchor['gbp_per_pg']==.978 and abs(anchor['derived_anchor_pg_1c']-16.1/.978)<1e-12
assert anchor['status']=='conditional_assembly_derived_not_measured_c_value' and anchor['anchor_uncertainty_included'] is False
print('FROZEN_TEMPLATE',{'rows':len(template),'per_image':template.filename.value_counts().to_dict(),'scaling_max_deltas':deltas,'source_sha256':sha(path),'conditional_anchor':anchor['derived_anchor_pg_1c']})
for ordinal in [3,17,32,37,41,45]:
 x=json.loads((R/'pixel_audit_v2/objects.json').read_text())[ordinal-1]
 assert x['source_record']=={str(k):str(v) for k,v in x['source_record'].items()}
 print('VISUAL_SOURCE',ordinal,x['row_key'],'CSV_LINE',x['csv_line'],'AREA',x['mask_area_px'],'IOD',x['stored_iod'],'BBOX',x['gallery_label_bbox_xyxy'])
print('FINISH',datetime.now(timezone.utc).isoformat(),'EXIT 0')
PY
```

Observed concluding output: `PIXEL_FINAL_MANIFEST 333 verified`; `EQUAL_IMAGE_IOD 943.889140125`; `STANDARD_MEDIANS {413: 594.0784, 414: 837.701417}`; `AUTHOR_MAP {'target': 47, 'current species denominator': 4, 'designated standard; QC unresolved': 2}`; `FROZEN_TEMPLATE rows=40`; exact six visual-source rows above; `EXIT 0`.

