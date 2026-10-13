#!/usr/bin/env python3
"""Bounded, read-only retained-mask reconstruction; --out must name a fresh directory."""

from __future__ import annotations

import argparse
import csv
import hashlib
import html
import importlib.metadata
import json
import math
import os
from pathlib import Path
import platform
import subprocess
import sys
import time
import traceback
from datetime import datetime, timezone

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parent
PAPER = ROOT.parents[1]
CELL = Path('/home/jake/Projects/cellprofiler_test')
SCRIPTS = CELL / 'scripts'
STANDARD = CELL / 'output/runs/life_history_pmount_reference_v1/nucleus_measurements/all_measurements.csv'
TARGET = PAPER / 'analyses/02_morphology/data/nuclear_iod_objects.csv'
TARGET_HASH = 'fffc91fa5f53a374711b7b71baf2292733e6185fa7d86bb02dc189a9cfc6d4bf'
HELPER_HASHES = {
    SCRIPTS / 'apply_pair_review_repairs.py': 'ae6e82575b63b50b237c5c85c99e375c14f1636171d905beecf8a76d0eae3327',
    SCRIPTS / 'audit_iod_measurement_quality.py': '41f205ed52ae353a3e836b8b0a1951b256444823ff99561b6189d8ba82cdc696',
}
METRICS = ('iod_inner_intensity_median', 'iod_ring_intensity_median',
           'iod_local_contrast', 'iod_local_od_contrast', 'iod_boundary_grad_mean',
           'iod_edge_sharpness', 'iod_ring_noise_iqr', 'iod_contrast_snr')
TOLERANCE = 1e-12


def clock():
    return datetime.now(timezone.utc).isoformat()


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def require_hash(path, expected):
    actual = digest(path)
    if actual != expected:
        raise ValueError(f'source hash mismatch: {path}: expected {expected}, actual {actual}')


def pin(path):
    path = Path(path)
    stat = path.stat()
    return {'sha256': digest(path), 'bytes': stat.st_size, 'mtime_ns': stat.st_mtime_ns}


def fresh_path(out, inputs=()):
    out = Path(os.path.abspath(out))
    for path in (out, *out.parents):
        if path.is_symlink():
            raise ValueError(f'symlinked output or ancestor: {path}')
    if out.exists():
        raise ValueError(f'output already exists: {out}')
    if not out.parent.is_dir():
        raise ValueError(f'output parent does not exist: {out.parent}')
    for source in inputs:
        source = Path(source).resolve()
        if out == source or out in source.parents or source in out.parents:
            raise ValueError(f'input/output ancestry overlap: {source}')
    return out


def records(path):
    with path.open(newline='') as stream:
        return [(line, row) for line, row in enumerate(csv.DictReader(stream), start=2)]


def rows_from_inventory(inventory):
    require_hash(STANDARD, inventory['pins'][str(STANDARD)]['sha256'])
    require_hash(TARGET, TARGET_HASH)
    for path, expected in HELPER_HASHES.items():
        require_hash(path, expected)
    standard = records(STANDARD)
    expected = inventory['standard_yolo_candidate_rows']
    if len(standard) != 47 or len(expected) != 47:
        raise ValueError('requires exactly 47 standard candidates')
    if [{'line': line, **row} for line, row in standard] != expected:
        raise ValueError('standard CSV rows differ from the inventoried records')
    rows = []
    for line, row in standard:
        key = f"{row['filename']}::{row['tile_name']}::{row['label']}"
        rows.append({'role': 'standard', 'csv_path': str(STANDARD), 'csv_line': line,
                     'row_key': key, 'source_record': row})
    canonical = records(TARGET)
    for process in ('314', '315', '320', '321'):
        line, row = next((line, row) for line, row in canonical
                         if row['filename'] == f'Process_{process}_raw_green.ome.tiff')
        rows.append({'role': 'target_control', 'csv_path': str(TARGET), 'csv_line': line,
                     'row_key': row['review_key'], 'source_record': row})
    if len({r['row_key'] for r in rows}) != 51:
        raise ValueError('duplicate row key')
    return rows


def metric_row(item):
    row = item['source_record']
    if item['role'] == 'target_control':
        return dict(row)
    projected = dict(row)
    for old, new in (('label', 'nucleus_label'), ('area_px', 'nuc_area_px'),
                     ('iod', 'nuc_iod'), ('i_bg', 'nuc_i_bg'),
                     ('mask_path', 'nucleus_mask_path'),
                     ('source_image_path', 'nucleus_source_image_path')):
        projected[new] = row[old]
    projected['nuc_centroid_x'] = float(row['centroid_x']) - float(row['tile_x0'])
    projected['nuc_centroid_y'] = float(row['centroid_y']) - float(row['tile_y0'])
    projected['tile_x0'] = projected['tile_y0'] = 0
    return projected


def validate_label(binary):
    if not binary.any():
        raise ValueError('recorded mask label absent')


def validate_standard_coordinates(projected, binary):
    validate_label(binary)
    coords = np.argwhere(binary)
    cy, cx = coords.mean(axis=0)
    x, y = float(projected['nuc_centroid_x']), float(projected['nuc_centroid_y'])
    if abs(cx - x) > .02 or abs(cy - y) > .02:
        raise ValueError(f'coordinate projection inconsistent: ({x},{y}) vs ({cx},{cy})')
    iy, ix = int(round(y)), int(round(x))
    if not (0 <= iy < binary.shape[0] and 0 <= ix < binary.shape[1] and binary[iy, ix]):
        raise ValueError('projected centroid does not lie on the recorded label')
    return {'full_label_centroid_xy': [float(cx), float(cy)],
            'projected_minus_mask_centroid_xy': [x - float(cx), y - float(cy)]}


def reconstruction(binary, projected):
    validate_label(binary)
    measured = repair.measure_repaired_iod(
        binary, raw_path=projected['nucleus_source_image_path'],
        local_slice=(slice(0, binary.shape[0]), slice(0, binary.shape[1])),
        i_bg=float(projected['nuc_i_bg']))
    stored = float(projected['nuc_iod'])
    area = int(binary.sum())
    return {'mask_area_px': area, 'stored_area_px': int(float(projected['nuc_area_px'])),
            'area_delta_px': area - int(float(projected['nuc_area_px'])),
            'recomputed_iod': measured['iod'], 'stored_iod': stored,
            'iod_delta': measured['iod'] - stored,
            'recomputed_iod_six_decimals': f"{measured['iod']:.6f}",
            'stored_iod_six_decimals': f'{stored:.6f}',
            'iod_six_decimal_match': f"{measured['iod']:.6f}" == f'{stored:.6f}',
            'area_exact_match': area == int(float(projected['nuc_area_px']))}


def helper_crop_diagnostics(projected, binary):
    # Use the actual native initial crop and fallback conditions, without changing the row.
    raw_path, mask_path = projected['nucleus_source_image_path'], projected['nucleus_mask_path']
    raw, labels = quality.load_row_tile_arrays(raw_path, mask_path)
    ox, oy = linked.as_int(projected['tile_x0']), linked.as_int(projected['tile_y0'])
    rx, ry = (0, 0) if linked.is_tile_artifact(Path(raw_path)) else (ox, oy)
    mx, my = linked.mask_origin(Path(mask_path), projected['tile_x0'], projected['tile_y0'], 0, 0)
    cx, cy = linked.as_int(projected['nuc_centroid_x']), linked.as_int(projected['nuc_centroid_y'])
    half = max(32, int(round(math.sqrt(max(float(projected['nuc_area_px']), 1) / math.pi) * 6)))
    initial = linked.crop_global(labels, mx, my, cx, cy, half)
    initial_area = int((initial == int(projected['nucleus_label'])).sum())
    raw_crop, mask_crop = quality.row_local_masks(pd.Series(projected))
    return {'helper_center_int_xy': [cx, cy], 'helper_half_width': half,
            'helper_raw_origin_xy': [rx, ry], 'helper_mask_origin_xy': [mx, my],
            'helper_initial_raw_bounds_unclipped_xyxy': [cx-rx-half, cy-ry-half, cx-rx+half, cy-ry+half],
            'helper_initial_mask_bounds_unclipped_xyxy': [cx-mx-half, cy-my-half, cx-mx+half, cy-my+half],
            'helper_initial_label_pixels': initial_area,
            'helper_full_tile_fallback': initial_area == 0,
            'helper_raw_shape': list(raw_crop.shape), 'helper_mask_shape': list(mask_crop.shape),
            'helper_label_pixels': int(mask_crop.sum()),
            'helper_complete_label_coverage': int(mask_crop.sum()) == int(binary.sum())}


def gallery_crop(item, gray, binary, out, ordinal):
    ys, xs = np.nonzero(binary)
    bbox = [int(xs.min()), int(ys.min()), int(xs.max()) + 1, int(ys.max()) + 1]
    h, w = binary.shape
    desired = [bbox[0]-64, bbox[1]-64, bbox[2]+64, bbox[3]+64]
    bounds = [max(0, desired[0]), max(0, desired[1]), min(w, desired[2]), min(h, desired[3])]
    x0, y0, x1, y1 = bounds
    raw = gray[y0:y1, x0:x1]
    mask = binary[y0:y1, x0:x1]
    rgb = np.repeat(raw[..., None], 3, axis=-1)
    rgb[find_boundaries(mask, mode='inner')] = [255, 32, 32]
    base = f"{ordinal:03d}_Process_{item['source_record']['slide_id']}_label_{metric_row(item)['nucleus_label']}"
    names = [f'crops/{base}_{kind}.png' for kind in ('raw', 'mask', 'overlay')]
    for name, array in zip(names, (raw, mask.astype(np.uint8)*255, rgb), strict=True):
        Image.fromarray(array).save(out / name)
    return {'gallery_images': names, 'gallery_label_bbox_xyxy': bbox,
            'gallery_desired_bounds_xyxy': desired, 'gallery_bounds_xyxy': bounds,
            'gallery_clipped': desired != bounds, 'gallery_pad_px': 64,
            'gallery_label_pixels': int(mask.sum()),
            'gallery_complete_label_coverage': int(mask.sum()) == int(binary.sum()),
            'gallery_mask_touches_tile_edge': bbox[0] == 0 or bbox[1] == 0 or bbox[2] == w or bbox[3] == h}


def controls(item):
    projected = metric_row(item)
    labels = linked.load_label_array(projected['nucleus_mask_path'])
    binary = labels == int(projected['nucleus_label'])
    raw = repair.grayscale_uint8(linked.load_raw_array(projected['nucleus_source_image_path']))
    tests = []

    def rejects(name, call):
        try:
            call()
        except ValueError as exc:
            tests.append({'name': name, 'rejected': True, 'evidence': str(exc)})
        else:
            tests.append({'name': name, 'rejected': False, 'evidence': 'unexpected acceptance'})

    rejects('absent_mask_label', lambda: validate_label(labels == int(labels.max()) + 1))
    wrong = dict(projected)
    wrong['nuc_centroid_x'] = item['source_record']['centroid_x']
    wrong['nuc_centroid_y'] = item['source_record']['centroid_y']
    rejects('unconverted_standard_coordinates', lambda: validate_standard_coordinates(wrong, binary))
    rejects('changed_source_hash', lambda: require_hash(projected['nucleus_source_image_path'], '0'*64))
    changed = dict(projected)
    changed['nuc_i_bg'] = float(projected['nuc_i_bg']) - 1
    result = reconstruction(binary, changed)
    tests.append({'name': 'changed_background', 'rejected': not result['iod_six_decimal_match'],
                  'changed_background': changed['nuc_i_bg'], 'evidence': result})
    ys, xs = np.nonzero(binary & (raw > 0) & (raw < float(projected['nuc_i_bg'])))
    if not len(ys):
        raise ValueError('falsification requires a positive-OD pixel')
    y, x = int(ys[0]), int(xs[0])
    removed = binary.copy()
    removed[y, x] = False
    result = reconstruction(removed, projected)
    tests.append({'name': 'removed_positive_od_pixel',
                  'rejected': not result['area_exact_match'] and not result['iod_six_decimal_match'],
                  'removed_pixel_xy': [x, y], 'removed_raw_value': int(raw[y, x]), 'evidence': result})
    changed = dict(projected)
    changed['nuc_iod'] = float(projected['nuc_iod']) + .000010
    result = reconstruction(binary, changed)
    tests.append({'name': 'changed_stored_iod', 'rejected': not result['iod_six_decimal_match'],
                  'evidence': result})
    return {'row_key': item['row_key'], 'tests': tests,
            'all_rejected': all(test['rejected'] for test in tests)}


def dump(path, value):
    path.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + '\n')


def gallery(out, results):
    cards = []
    for result in results:
        images = ''.join(f'<figure><img src="{html.escape(name)}" alt="{kind}: {html.escape(result["row_key"])}"><figcaption>{kind}</figcaption></figure>'
                         for name, kind in zip(result['gallery_images'], ('RAW', 'MASK', 'OVERLAY'), strict=True))
        coords = result['coordinate_encoding']
        cards.append(f'<article><h2>{result["ordinal"]:03d} {html.escape(result["row_key"])}</h2>'
                     f'<p>{result["role"]}; specimen {html.escape(result["source_record"]["specimen_id"])}; '
                     f'CSV line {result["csv_line"]}; bounds xyxy {result["gallery_bounds_xyxy"]}; '
                     f'clipped={result["gallery_clipped"]}; full-label pixels={result["gallery_label_pixels"]}; '
                     f'helper fallback={result["helper_full_tile_fallback"]}</p>'
                     f'<p>Coordinates: {html.escape(json.dumps(coords, sort_keys=True))}</p>'
                     f'<div class="triplet">{images}</div></article>')
    (out/'gallery.html').write_text('<!doctype html><meta charset="utf-8"><title>Retained mask pixel audit</title>'
        '<style>body{font:16px system-ui;margin:24px;background:#eee}article{background:white;padding:16px;margin:24px 0}'
        'h2{font-size:17px;overflow-wrap:anywhere}.triplet{display:flex;gap:16px;flex-wrap:wrap}figure{margin:0}'
        'img{max-width:440px;height:auto;border:1px solid #bbb}figcaption{font-weight:bold}p{overflow-wrap:anywhere}</style>'
        '<h1>Retained mask pixel audit: all 47 standards and four target controls</h1>'
        '<p>Unadjusted uint8 RAW display, binary recorded label MASK, red inner-boundary OVERLAY. '
        'No IOD selection or manual acceptance. Crops use the full label bounding box plus 64 pixels, clipped at the tile edge. '
        'Gallery crops are independent of helper metric crops. Numerical reconstruction does not qualify calibration.</p>'
        + ''.join(cards) + '\n')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', required=True, type=Path)
    args = parser.parse_args()
    started, monotonic = clock(), time.monotonic()
    try:
        out = fresh_path(args.out)
        inventory = json.loads((ROOT/'source_inventory.json').read_text())
        rows = rows_from_inventory(inventory)
        source_paths = {Path(__file__).resolve(), ROOT/'source_inventory.json', TARGET, STANDARD,
                        ROOT/'AUDIT_DESIGN.md', ROOT/'design_b/DESIGN.md',
                        ROOT/'SOURCE_INVENTORY.md', ROOT/'MEASUREMENT_FLOW.md',
                        ROOT/'AUTHOR_IMAGE_CONFIRMATION.json', PAPER/'pyproject.toml', PAPER/'uv.lock'}
        for item in rows:
            row = metric_row(item)
            source_paths.update((Path(row['nucleus_source_image_path']), Path(row['nucleus_mask_path'])))
        out = fresh_path(out, source_paths)
        before = {str(path): pin(path) for path in sorted(source_paths)}
    except (ValueError, OSError, KeyError, StopIteration) as exc:
        print(f'PREFLIGHT REFUSAL: {exc}', file=sys.stderr)
        return 2
    out.mkdir(exist_ok=False)
    receipt = {'start_utc': started, 'argv': sys.argv, 'cwd': os.getcwd(),
               'uv_project_environment': os.environ.get('UV_PROJECT_ENVIRONMENT'),
               'inputs_before': before, 'metric_absolute_tolerance': TOLERANCE,
               'scope': 'Retained-mask reconstruction and pixel metrics only; no manual acceptance or calibration',
               'exit_status': 1}
    results = []
    try:
        global np, pd, repair, quality, linked, Image, find_boundaries
        sys.path.insert(0, str(SCRIPTS))
        import numpy as np
        import pandas as pd
        import apply_pair_review_repairs as repair
        import audit_iod_measurement_quality as quality
        import run_linked_pair_review as linked
        from PIL import Image
        from skimage.segmentation import find_boundaries

        helpers = {Path(module.__file__).resolve() for module in list(sys.modules.values())
                   if getattr(module, '__file__', None) and Path(module.__file__).resolve().is_relative_to(SCRIPTS)}
        source_paths.update(helpers)
        receipt['inputs_before'].update({str(path): pin(path) for path in sorted(helpers)})
        receipt['helper_source_paths'] = [str(path) for path in sorted(helpers)]
        receipt['python'] = {'version': platform.python_version(), 'executable': sys.executable,
                             'implementation': platform.python_implementation()}
        receipt['packages'] = {name: importlib.metadata.version(name) for name in
                               ('numpy', 'pandas', 'Pillow', 'tifffile', 'scipy', 'scikit-image', 'matplotlib', 'roifile')}
        receipt['heads'] = {str(repo): subprocess.check_output(['git', '-C', str(repo), 'rev-parse', 'HEAD'], text=True).strip()
                            for repo in (PAPER, CELL)}
        (out/'crops').mkdir()
        for ordinal, item in enumerate(rows, start=1):
            row = metric_row(item)
            labels = linked.load_label_array(row['nucleus_mask_path'])
            raw = linked.load_raw_array(row['nucleus_source_image_path'])
            gray = repair.grayscale_uint8(raw)
            if labels.ndim != 2 or gray.shape != labels.shape or raw.dtype != np.uint8:
                raise ValueError(f'invalid retained tile dimensions or depth: {item["row_key"]}')
            binary = labels == int(row['nucleus_label'])
            validate_label(binary)
            result = {**item, 'ordinal': ordinal, 'source_csv_sha256': before[item['csv_path']]['sha256'],
                      'raw_path': row['nucleus_source_image_path'], 'mask_path': row['nucleus_mask_path'],
                      'raw_sha256': before[row['nucleus_source_image_path']]['sha256'],
                      'mask_sha256': before[row['nucleus_mask_path']]['sha256'],
                      'tile_shape': list(labels.shape), 'raw_dtype': str(raw.dtype), 'mask_dtype': str(labels.dtype),
                      'coordinate_encoding': {'mode': 'standard_global_to_local_copy' if item['role']=='standard' else 'target_historical_unchanged',
                          'original_centroid_xy': [item['source_record'].get('centroid_x', item['source_record'].get('nuc_centroid_x')),
                                                   item['source_record'].get('centroid_y', item['source_record'].get('nuc_centroid_y'))],
                          'original_tile_origin_xy': [item['source_record']['tile_x0'], item['source_record']['tile_y0']],
                          'helper_centroid_xy': [row['nuc_centroid_x'], row['nuc_centroid_y']],
                          'helper_tile_origin_xy': [row['tile_x0'], row['tile_y0']]},
                      'metric_input_record': row, 'issues': []}
            if item['role'] == 'standard':
                result.update(validate_standard_coordinates(row, binary))
            result.update(reconstruction(binary, row))
            result.update(helper_crop_diagnostics(row, binary))
            if result['helper_raw_shape'] != result['helper_mask_shape']:
                raise ValueError(f'helper crop shapes disagree: {item["row_key"]}')
            result['pixel_metrics'] = quality.compute_row_iod_metrics(pd.Series(row))
            result['stored_pixel_metric_deltas'] = {}
            if item['role'] == 'target_control':
                result['stored_pixel_metric_deltas'] = {key: result['pixel_metrics'][key] - float(item['source_record'][key]) for key in METRICS}
                for key, delta in result['stored_pixel_metric_deltas'].items():
                    if not math.isfinite(delta) or abs(delta) > TOLERANCE:
                        result['issues'].append(f'historical metric parity failure: {key}: {delta!r}')
            for name in ('iod_six_decimal_match', 'area_exact_match', 'helper_complete_label_coverage'):
                if not result[name]:
                    result['issues'].append(name + ' failed')
            result.update(gallery_crop(item, gray, binary, out, ordinal))
            if not result['gallery_complete_label_coverage']:
                result['issues'].append('gallery_complete_label_coverage failed')
            if not all(math.isfinite(v) for v in result['pixel_metrics'].values()):
                result['issues'].append('nonfinite pixel metric')
            results.append(result)
            print(f'{ordinal:02d}/51 {item["row_key"]}: ' + ('ISSUES '+str(result['issues']) if result['issues'] else 'reconstruction/parity PASS'), flush=True)
            # Bound retained-array caches to the current tile; no source writes.
            quality.load_row_tile_arrays.cache_clear()
            linked.load_raw_array.cache_clear()
            linked.load_label_array.cache_clear()
        falsification = controls(rows[0])
        dump(out/'controls.json', falsification)
        dump(out/'objects.json', results)
        fields = ['ordinal', 'role', 'row_key', 'csv_path', 'csv_line', 'raw_path', 'mask_path',
                  'mask_area_px', 'stored_area_px', 'area_delta_px', 'recomputed_iod', 'stored_iod',
                  'iod_delta', 'iod_six_decimal_match', 'area_exact_match', 'helper_full_tile_fallback',
                  'helper_complete_label_coverage', 'gallery_clipped', 'gallery_complete_label_coverage', 'issues']
        with (out/'objects.csv').open('w', newline='') as stream:
            writer = csv.DictWriter(stream, fieldnames=fields + list(METRICS))
            writer.writeheader()
            for result in results:
                writer.writerow({**{key: result[key] for key in fields}, **result['pixel_metrics']})
        gallery(out, results)
        receipt['counts'] = {'standard': 47, 'target_control': 4, 'gallery_objects': len(results), 'gallery_images': len(results)*3}
        receipt['max_target_metric_absolute_delta'] = max(abs(delta) for result in results for delta in result['stored_pixel_metric_deltas'].values())
        receipt['falsification_all_rejected'] = falsification['all_rejected']
        receipt['torch_loaded'] = any(name == 'torch' or name.startswith('torch.') for name in sys.modules)
        receipt['issues'] = [{'row_key': r['row_key'], 'issues': r['issues']} for r in results if r['issues']]
        receipt['exit_status'] = int(bool(receipt['issues']) or not falsification['all_rejected'] or receipt['torch_loaded'])
    except Exception:
        receipt['execution_exception'] = traceback.format_exc()
        print(receipt['execution_exception'], file=sys.stderr)
        if results:
            dump(out/'partial_objects.json', results)
    finally:
        receipt['inputs_after'] = {str(path): pin(path) for path in sorted(source_paths)}
        receipt['changed_inputs'] = [path for path, value in receipt['inputs_before'].items() if value != receipt['inputs_after'][path]]
        if receipt['changed_inputs']:
            receipt['exit_status'] = 1
        receipt['finish_utc'] = clock()
        receipt['elapsed_seconds'] = time.monotonic() - monotonic
        receipt['artifacts'] = {str(path.relative_to(out)): pin(path) for path in sorted(out.rglob('*')) if path.is_file()}
        dump(out/'receipt.json', receipt)
    print(f"EXIT {receipt['exit_status']}; retained receipt: {out/'receipt.json'}", flush=True)
    return receipt['exit_status']


if __name__ == '__main__':
    sys.exit(main())
