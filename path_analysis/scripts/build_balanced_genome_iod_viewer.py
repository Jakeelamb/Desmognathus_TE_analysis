#!/usr/bin/env python3
"""Build a static mask-validation viewer for balanced genome-IOD pairs."""

from __future__ import annotations

import argparse
import html
import json
import math
from functools import lru_cache
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from PIL import Image


PROJECT_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_SELECTED_PAIRS = (
    PROJECT_ROOT
    / "path_analysis"
    / "data"
    / "external"
    / "derived"
    / "balanced_genome_iod_sensitivity"
    / "balanced_genome_iod_selected_pairs.csv.gz"
)
DEFAULT_OUTPUT_DIR = (
    PROJECT_ROOT
    / "path_analysis"
    / "data"
    / "external"
    / "derived"
    / "balanced_genome_iod_sensitivity"
    / "viewer_balanced_qc_only"
)
DEFAULT_PANEL = "balanced_qc_only"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--selected-pairs-csv", type=Path, default=DEFAULT_SELECTED_PAIRS)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--panel", default=DEFAULT_PANEL)
    parser.add_argument("--thumb-size", type=int, default=384)
    parser.add_argument("--crop-margin", type=int, default=64)
    parser.add_argument("--title", default="Balanced Genome-IOD Mask Viewer")
    parser.add_argument(
        "--default-decision",
        choices=["keep"],
        default=None,
        help="Treat unmarked records as this decision and export an audit row for every record.",
    )
    parser.add_argument(
        "--storage-key",
        default=None,
        help="Browser localStorage key for decisions. Defaults to a panel-specific key.",
    )
    parser.add_argument(
        "--blind-target-metrics",
        action="store_true",
        help="Hide cell area, nucleus area, OD, and IOD during visual review.",
    )
    parser.add_argument(
        "--common-support-only",
        action="store_true",
        help="Show and export only records whose quality_match_status is common_support.",
    )
    parser.add_argument(
        "--sort-by-match-distance",
        action="store_true",
        help="Order records within species from strongest to weakest technical match.",
    )
    parser.add_argument("--max-rows", type=int, default=0, help="Debug limiter; 0 means all rows for the panel.")
    parser.add_argument("--force", action="store_true", help="Regenerate image assets even when they already exist.")
    return parser.parse_args()


def require_exists(path: Path) -> None:
    if not path.exists():
        raise FileNotFoundError(path)


def clean_text(value: Any) -> str:
    if value is None:
        return ""
    try:
        if pd.isna(value):
            return ""
    except Exception:
        pass
    return str(value).strip()


def as_int(value: Any, default: int = 0) -> int:
    try:
        if pd.isna(value):
            return default
    except Exception:
        pass
    try:
        return int(round(float(value)))
    except (TypeError, ValueError):
        return default


def species_slug(species: str) -> str:
    text = species.replace("D. ", "").replace(" ", "_")
    return "".join(ch if ch.isalnum() or ch in {"_", "-"} else "_" for ch in text).strip("_") or "unknown"


@lru_cache(maxsize=8)
def load_raw(path_text: str) -> np.ndarray:
    path = Path(path_text)
    arr = np.asarray(Image.open(path))
    if arr.ndim == 3:
        if arr.shape[0] in {1, 3, 4}:
            arr = arr[0]
        else:
            arr = arr[..., 0]
    if arr.dtype == np.uint16:
        arr = (arr >> 8).astype(np.uint8)
    elif arr.dtype != np.uint8:
        arr = np.clip(arr, 0, 255).astype(np.uint8)
    return arr


@lru_cache(maxsize=8)
def load_mask(path_text: str) -> np.ndarray:
    return np.asarray(Image.open(path_text))


def safe_label_mask(path_text: str, label: int) -> np.ndarray | None:
    if not path_text or label <= 0:
        return None
    path = Path(path_text)
    if not path.exists():
        return None
    arr = load_mask(str(path))
    return np.asarray(arr == label, dtype=bool)


def bbox_for_masks(masks: list[np.ndarray | None]) -> tuple[int, int, int, int] | None:
    boxes: list[tuple[int, int, int, int]] = []
    for mask in masks:
        if mask is None or not np.any(mask):
            continue
        coords = np.argwhere(mask)
        y0, x0 = coords.min(axis=0)
        y1, x1 = coords.max(axis=0) + 1
        boxes.append((int(y0), int(y1), int(x0), int(x1)))
    if not boxes:
        return None
    return (
        min(box[0] for box in boxes),
        max(box[1] for box in boxes),
        min(box[2] for box in boxes),
        max(box[3] for box in boxes),
    )


def square_crop_bounds(
    bbox: tuple[int, int, int, int],
    shape: tuple[int, int],
    *,
    margin: int,
) -> tuple[int, int, int, int]:
    y0, y1, x0, x1 = bbox
    height, width = shape[:2]
    cy = (y0 + y1) / 2.0
    cx = (x0 + x1) / 2.0
    half = int(math.ceil(max(y1 - y0, x1 - x0) / 2.0 + margin))
    half = max(96, half)
    top = max(0, int(round(cy - half)))
    bottom = min(height, int(round(cy + half)))
    left = max(0, int(round(cx - half)))
    right = min(width, int(round(cx + half)))
    if bottom <= top:
        bottom = min(height, top + 1)
    if right <= left:
        right = min(width, left + 1)
    return top, bottom, left, right


def crop_and_resize(arr: np.ndarray, bounds: tuple[int, int, int, int], thumb_size: int, *, nearest: bool) -> Image.Image:
    top, bottom, left, right = bounds
    crop = arr[top:bottom, left:right]
    if crop.dtype == bool:
        img = Image.fromarray((crop.astype(np.uint8) * 255), mode="L")
    else:
        img = Image.fromarray(np.asarray(crop, dtype=np.uint8), mode="L")
    resample = Image.Resampling.NEAREST if nearest else Image.Resampling.BILINEAR
    return img.resize((thumb_size, thumb_size), resample=resample)


def mask_overlay(mask_img: Image.Image, color: tuple[int, int, int], alpha: int) -> Image.Image:
    mask = np.asarray(mask_img.convert("L")) > 0
    rgba = np.zeros((mask.shape[0], mask.shape[1], 4), dtype=np.uint8)
    rgba[mask, 0] = color[0]
    rgba[mask, 1] = color[1]
    rgba[mask, 2] = color[2]
    rgba[mask, 3] = alpha
    return Image.fromarray(rgba, mode="RGBA")


def process_row(row: pd.Series, index: int, output_dir: Path, *, thumb_size: int, crop_margin: int, force: bool) -> dict[str, Any]:
    species = clean_text(row.get("species")) or "unknown"
    slug = species_slug(species)
    asset_dir = output_dir / "assets" / slug
    asset_dir.mkdir(parents=True, exist_ok=True)
    prefix = f"{index:04d}"
    raw_rel = Path("assets") / slug / f"{prefix}_raw.jpg"
    cell_rel = Path("assets") / slug / f"{prefix}_cell.png"
    nucleus_rel = Path("assets") / slug / f"{prefix}_nucleus.png"
    raw_path = output_dir / raw_rel
    cell_path = output_dir / cell_rel
    nucleus_path = output_dir / nucleus_rel

    cell_label = as_int(row.get("mask_label_id"))
    nucleus_label = as_int(row.get("nucleus_label"))

    assets_ready = raw_path.exists() and cell_path.exists() and nucleus_path.exists()
    missing_mask = False
    if force or not assets_ready:
        cell_mask = safe_label_mask(clean_text(row.get("cell_mask_path")), cell_label)
        nucleus_mask = safe_label_mask(clean_text(row.get("nucleus_mask_path")), nucleus_label)
        raw_source = clean_text(row.get("nucleus_source_image_path")) or clean_text(row.get("cell_source_image_path"))
        raw = load_raw(raw_source)
        bbox = bbox_for_masks([cell_mask, nucleus_mask])
        missing_mask = bbox is None
        if bbox is None:
            center = raw.shape[0] // 2, raw.shape[1] // 2
            bbox = (center[0] - 96, center[0] + 96, center[1] - 96, center[1] + 96)
        bounds = square_crop_bounds(bbox, raw.shape, margin=crop_margin)

        raw_img = crop_and_resize(raw, bounds, thumb_size, nearest=False).convert("RGB")
        raw_img.save(raw_path, quality=88, optimize=True)

        if cell_mask is not None:
            cell_img = crop_and_resize(cell_mask, bounds, thumb_size, nearest=True)
        else:
            cell_img = Image.new("L", (thumb_size, thumb_size), 0)
        mask_overlay(cell_img, (255, 80, 40), 96).save(cell_path)

        if nucleus_mask is not None:
            nucleus_img = crop_and_resize(nucleus_mask, bounds, thumb_size, nearest=True)
        else:
            nucleus_img = Image.new("L", (thumb_size, thumb_size), 0)
        mask_overlay(nucleus_img, (25, 190, 255), 128).save(nucleus_path)

    initial_decision = clean_text(row.get("viewer_default_decision")).lower()
    if initial_decision not in {"", "keep", "problem", "unsure"}:
        raise ValueError(f"Unsupported viewer_default_decision: {initial_decision!r}")
    return {
        "id": f"{slug}-{index:04d}",
        "index": index,
        "panel": clean_text(row.get("panel")),
        "species": species,
        "species_slug": slug,
        "rank": as_int(row.get("selection_rank")),
        "review_key": clean_text(row.get("review_key")),
        "filename": clean_text(row.get("filename")),
        "tile_name": clean_text(row.get("tile_name")),
        "cell_label": cell_label,
        "nucleus_label": nucleus_label,
        "cell_area_um2": float(pd.to_numeric(pd.Series([row.get("cell_area_um2")]), errors="coerce").iloc[0]),
        "nucleus_area_um2": float(pd.to_numeric(pd.Series([row.get("nuc_area_um2")]), errors="coerce").iloc[0]),
        "nuc_iod": float(pd.to_numeric(pd.Series([row.get("nuc_iod")]), errors="coerce").iloc[0]),
        "nuc_mean_od": float(pd.to_numeric(pd.Series([row.get("nuc_mean_od")]), errors="coerce").iloc[0]),
        "qc_score": float(pd.to_numeric(pd.Series([row.get("balanced_genome_qc_score")]), errors="coerce").iloc[0]),
        "clarity": float(pd.to_numeric(pd.Series([row.get("clarity_signal")]), errors="coerce").iloc[0]),
        "quality_match_distance": float(
            pd.to_numeric(pd.Series([row.get("quality_match_distance")]), errors="coerce").iloc[0]
        ),
        "quality_match_status": clean_text(row.get("quality_match_status")),
        "initial_decision": initial_decision,
        "review_cohort": clean_text(row.get("review_cohort")),
        "raw_src": raw_rel.as_posix(),
        "cell_overlay_src": cell_rel.as_posix(),
        "nucleus_overlay_src": nucleus_rel.as_posix(),
        "missing_mask": bool(missing_mask),
    }


def build_records(args: argparse.Namespace) -> list[dict[str, Any]]:
    require_exists(args.selected_pairs_csv)
    df = pd.read_csv(args.selected_pairs_csv, low_memory=False)
    df = df.loc[df["panel"].eq(args.panel)].copy()
    if df.empty:
        raise ValueError(f"No rows found for panel {args.panel!r}")
    df = df.sort_values(["species", "selection_rank", "review_key"], kind="mergesort").reset_index(drop=True)
    if args.max_rows > 0:
        df = df.head(args.max_rows).copy()
    records: list[dict[str, Any]] = []
    for idx, row in df.iterrows():
        records.append(
            process_row(
                row,
                idx + 1,
                args.output_dir,
                thumb_size=args.thumb_size,
                crop_margin=args.crop_margin,
                force=args.force,
            )
        )
        if (idx + 1) % 100 == 0:
            print(f"processed {idx + 1}/{len(df)} rows", flush=True)
    return records


def write_index(
    output_dir: Path,
    records: list[dict[str, Any]],
    *,
    panel: str,
    title: str = "Balanced Genome-IOD Mask Viewer",
    storage_key: str | None = None,
    default_decision: str | None = None,
    blind_target_metrics: bool = False,
    common_support_only: bool = False,
    sort_by_match_distance: bool = False,
    auto_fill_target_per_species: int | None = None,
) -> Path:
    ordered_records = list(records)
    if sort_by_match_distance:
        ordered_records.sort(
            key=lambda record: (
                record["species"],
                record["quality_match_distance"]
                if np.isfinite(record["quality_match_distance"])
                else float("inf"),
                record["review_key"],
            )
        )
    review_records = [
        record
        for record in ordered_records
        if not common_support_only or record["quality_match_status"] == "common_support"
    ]
    species = sorted({record["species"] for record in review_records})
    payload = json.dumps(ordered_records, separators=(",", ":"))
    species_payload = json.dumps(species, separators=(",", ":"))
    page_title = html.escape(title)
    storage_key_payload = json.dumps(storage_key or f"mask-viewer-decisions-{panel}-v1")
    export_filename_payload = json.dumps(f"{panel}_decisions.csv")
    default_decision_payload = json.dumps(default_decision or "")
    blind_target_metrics_payload = json.dumps(bool(blind_target_metrics))
    common_support_only_payload = json.dumps(bool(common_support_only))
    auto_fill_target_payload = json.dumps(int(auto_fill_target_per_species or 0))
    index_path = output_dir / "index.html"
    index_path.write_text(
        f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{page_title}</title>
<style>
:root {{
  color-scheme: light;
  --bg: #f5f6f4;
  --panel: #ffffff;
  --ink: #1d2528;
  --muted: #647076;
  --line: #d9dedb;
  --accent: #235789;
  --cell: #ff5028;
  --nucleus: #10a7e2;
}}
* {{ box-sizing: border-box; }}
body {{ margin: 0; font-family: Inter, ui-sans-serif, system-ui, -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif; background: var(--bg); color: var(--ink); }}
header {{ position: sticky; top: 0; z-index: 10; background: rgba(245,246,244,0.96); border-bottom: 1px solid var(--line); backdrop-filter: blur(10px); }}
.toolbar {{ display: grid; grid-template-columns: minmax(220px, 1fr) auto auto; gap: 12px; align-items: center; padding: 12px 16px; }}
.title {{ font-size: 15px; font-weight: 700; }}
.title span {{ color: var(--muted); font-weight: 500; margin-left: 8px; }}
.controls, .toggles {{ display: flex; flex-wrap: wrap; gap: 8px; align-items: center; justify-content: flex-end; }}
select, input[type="search"], button, label.toggle {{ border: 1px solid var(--line); border-radius: 6px; background: #fff; color: var(--ink); height: 34px; padding: 0 10px; font-size: 13px; }}
input[type="search"] {{ min-width: 220px; }}
button {{ cursor: pointer; }}
button.active, .status-btn.active {{ border-color: var(--accent); background: #e8f0f7; }}
label.toggle {{ display: inline-flex; align-items: center; gap: 6px; cursor: pointer; }}
main {{ padding: 14px 16px 24px; }}
.summary {{ display: flex; gap: 10px; align-items: center; flex-wrap: wrap; margin-bottom: 12px; color: var(--muted); font-size: 13px; }}
.grid {{ display: grid; grid-template-columns: repeat(auto-fill, minmax(238px, 1fr)); gap: 12px; }}
.card {{ background: var(--panel); border: 1px solid var(--line); border-radius: 8px; overflow: hidden; box-shadow: 0 1px 2px rgba(20,30,35,0.04); }}
.image-stack {{ position: relative; aspect-ratio: 1 / 1; background: #111; }}
.image-stack img {{ position: absolute; inset: 0; width: 100%; height: 100%; object-fit: cover; image-rendering: auto; }}
.overlay.cell {{ display: var(--show-cell); }}
.overlay.nucleus {{ display: var(--show-nucleus); }}
.card-body {{ padding: 9px; display: grid; gap: 7px; }}
.meta-row {{ display: flex; justify-content: space-between; gap: 8px; font-size: 12px; line-height: 1.25; }}
.meta-row strong {{ font-size: 13px; }}
.muted {{ color: var(--muted); }}
.metrics {{ display: grid; grid-template-columns: repeat(2, minmax(0, 1fr)); gap: 5px; font-size: 11px; color: var(--muted); }}
.metric {{ border: 1px solid var(--line); border-radius: 6px; padding: 4px 5px; min-width: 0; }}
.metric b {{ color: var(--ink); font-weight: 650; }}
.blind-targets .target-metric {{ display: none; }}
.status-row {{ display: grid; grid-template-columns: 1fr 1fr 1fr; gap: 5px; }}
.status-btn {{ height: 28px; border-radius: 6px; font-size: 12px; padding: 0 4px; }}
.card.problem {{ border-color: #c44; }}
.card.keep {{ border-color: #23864b; }}
.card.unsure {{ border-color: #b98b1d; }}
.card.standby {{ opacity: 0.72; }}
.selection-tag {{ border: 1px solid var(--line); border-radius: 999px; padding: 2px 6px; font-size: 10px; color: var(--muted); white-space: nowrap; }}
.selection-tag.auto {{ border-color: #23864b; color: #176238; background: #edf8f1; }}
.selection-tag.standby {{ background: #f0f1ef; }}
.hidden {{ display: none !important; }}
@media (max-width: 780px) {{
  .toolbar {{ grid-template-columns: 1fr; align-items: stretch; }}
  .controls, .toggles {{ justify-content: flex-start; }}
  input[type="search"] {{ min-width: 100%; }}
}}
</style>
</head>
<body>
<header>
  <div class="toolbar">
    <div class="title">{page_title} <span>{html.escape(panel)} | {len(review_records)} pairs</span></div>
    <div class="controls">
      <select id="speciesSelect" aria-label="Species"></select>
      <input id="searchBox" type="search" placeholder="search image, tile, key">
      <button id="exportBtn">Export decisions</button>
    </div>
    <div class="toggles">
      <label class="toggle"><input id="cellToggle" type="checkbox" checked> cell mask</label>
      <label class="toggle"><input id="nucleusToggle" type="checkbox" checked> nucleus mask</label>
      <button id="clearBtn">Clear marks</button>
    </div>
  </div>
</header>
<main>
  <div id="summary" class="summary"></div>
  <div id="grid" class="grid"></div>
</main>
<script>
const RECORDS = {payload};
const SPECIES = {species_payload};
const STORAGE_KEY = {storage_key_payload};
const EXPORT_FILENAME = {export_filename_payload};
const DEFAULT_DECISION = {default_decision_payload};
const EXPORT_ALL_RECORDS = Boolean(DEFAULT_DECISION);
const EXPORT_INITIAL_RECORDS = RECORDS.some(record => record.initial_decision);
const AUTO_FILL_TARGET_PER_SPECIES = {auto_fill_target_payload};
const AUTO_FILL_COHORT = "new_species_unreviewed";
const BLIND_TARGET_METRICS = {blind_target_metrics_payload};
if (BLIND_TARGET_METRICS) document.documentElement.classList.add("blind-targets");
const COMMON_SUPPORT_ONLY = {common_support_only_payload};
const REVIEW_RECORDS = COMMON_SUPPORT_ONLY
  ? RECORDS.filter(record => record.quality_match_status === "common_support")
  : RECORDS;
const state = {{
  species: "all",
  search: "",
  showCell: true,
  showNucleus: true,
  decisions: JSON.parse(localStorage.getItem(STORAGE_KEY) || "{{}}")
}};
const grid = document.getElementById("grid");
const summary = document.getElementById("summary");
const speciesSelect = document.getElementById("speciesSelect");
const searchBox = document.getElementById("searchBox");
const cellToggle = document.getElementById("cellToggle");
const nucleusToggle = document.getElementById("nucleusToggle");

function fmt(value, digits = 2) {{
  const n = Number(value);
  if (!Number.isFinite(n)) return "NA";
  return n.toFixed(digits);
}}
function readBoolParam(params, key, fallback) {{
  if (!params.has(key)) return fallback;
  const value = String(params.get(key) || "").toLowerCase();
  return !["0", "false", "off", "no"].includes(value);
}}
function applyUrlState() {{
  const params = new URLSearchParams(window.location.search);
  const species = params.get("species");
  if (species && SPECIES.includes(species)) state.species = species;
  const search = params.get("search");
  if (search) state.search = search;
  state.showCell = readBoolParam(params, "cellMask", state.showCell);
  state.showNucleus = readBoolParam(params, "nucleusMask", state.showNucleus);
}}
function syncControls() {{
  speciesSelect.value = state.species;
  searchBox.value = state.search;
  cellToggle.checked = state.showCell;
  nucleusToggle.checked = state.showNucleus;
}}
function syncUrl() {{
  const params = new URLSearchParams();
  if (state.species !== "all") params.set("species", state.species);
  if (state.search) params.set("search", state.search);
  if (!state.showCell) params.set("cellMask", "0");
  if (!state.showNucleus) params.set("nucleusMask", "0");
  const query = params.toString();
  const nextUrl = `${{window.location.pathname}}${{query ? `?${{query}}` : ""}}`;
  window.history.replaceState(null, "", nextUrl);
}}
function saveDecisions() {{
  localStorage.setItem(STORAGE_KEY, JSON.stringify(state.decisions));
}}
function setMaskCss() {{
  document.documentElement.style.setProperty("--show-cell", state.showCell ? "block" : "none");
  document.documentElement.style.setProperty("--show-nucleus", state.showNucleus ? "block" : "none");
}}
function speciesOptions() {{
  speciesSelect.innerHTML = `<option value="all">All species</option>` + SPECIES.map(s => `<option value="${{s}}">${{s}}</option>`).join("");
}}
function matches(record) {{
  if (state.species !== "all" && record.species !== state.species) return false;
  if (!state.search) return true;
  const hay = `${{record.species}} ${{record.filename}} ${{record.tile_name}} ${{record.review_key}}`.toLowerCase();
  return hay.includes(state.search.toLowerCase());
}}
function rankValue(record) {{
  const rank = Number(record.rank);
  return Number.isFinite(rank) ? rank : Number.POSITIVE_INFINITY;
}}
function computeAutoSelectedIds() {{
  const selected = new Set();
  if (AUTO_FILL_TARGET_PER_SPECIES <= 0) return selected;
  const groups = new Map();
  REVIEW_RECORDS
    .filter(record => record.review_cohort === AUTO_FILL_COHORT)
    .forEach(record => {{
      if (!groups.has(record.species)) groups.set(record.species, []);
      groups.get(record.species).push(record);
    }});
  groups.forEach(records => {{
    records.sort((a, b) => rankValue(a) - rankValue(b) || a.review_key.localeCompare(b.review_key));
    const explicitKeeps = records.filter(record => state.decisions[record.id]?.decision === "keep");
    explicitKeeps.forEach(record => selected.add(record.id));
    let remaining = Math.max(0, AUTO_FILL_TARGET_PER_SPECIES - explicitKeeps.length);
    records.forEach(record => {{
      if (remaining <= 0 || selected.has(record.id)) return;
      const explicit = state.decisions[record.id]?.decision || "";
      if (explicit === "problem" || explicit === "unsure") return;
      selected.add(record.id);
      remaining -= 1;
    }});
  }});
  return selected;
}}
function resolvedDecision(record, autoSelectedIds = computeAutoSelectedIds()) {{
  const explicit = state.decisions[record.id]?.decision || "";
  if (explicit) return explicit;
  const initial = record.initial_decision || DEFAULT_DECISION || "";
  if (initial) return initial;
  return autoSelectedIds.has(record.id) ? "keep" : "";
}}
function decisionSource(record, autoSelectedIds) {{
  if (state.decisions[record.id]) return "explicit";
  if (record.initial_decision || DEFAULT_DECISION) {{
    return state.decisions[record.id] ? "explicit" : "initial_default";
  }}
  if (autoSelectedIds.has(record.id)) return "auto_fill";
  return "";
}}
function card(record, autoSelectedIds) {{
  const explicitDecision = state.decisions[record.id]?.decision || "";
  const decision = resolvedDecision(record, autoSelectedIds);
  const autoEligible = record.review_cohort === AUTO_FILL_COHORT && AUTO_FILL_TARGET_PER_SPECIES > 0;
  const standby = autoEligible && !decision && !autoSelectedIds.has(record.id);
  const cls = `${{decision ? ` ${{decision}}` : ""}}${{standby ? " standby" : ""}}`;
  const selectionLabel = explicitDecision
    ? `manual ${{explicitDecision}}`
    : record.initial_decision
      ? "previously frozen"
      : autoSelectedIds.has(record.id)
        ? `auto top ${{AUTO_FILL_TARGET_PER_SPECIES}}`
        : standby
          ? "standby replacement"
          : "";
  return `<article class="card${{cls}}" data-id="${{record.id}}">
    <div class="image-stack">
      <img src="${{record.raw_src}}" loading="lazy" alt="">
      <img class="overlay cell" src="${{record.cell_overlay_src}}" loading="lazy" alt="">
      <img class="overlay nucleus" src="${{record.nucleus_overlay_src}}" loading="lazy" alt="">
    </div>
    <div class="card-body">
      <div class="meta-row"><strong>${{record.species}} #${{record.rank}}</strong>${{selectionLabel ? `<span class="selection-tag ${{standby ? "standby" : autoSelectedIds.has(record.id) ? "auto" : ""}}">${{selectionLabel}}</span>` : ""}}</div>
      <div class="meta-row"><span class="muted">${{record.filename.replace("_raw_green.ome.tiff", "")}}</span><span></span></div>
      <div class="meta-row"><span class="muted">${{record.tile_name}}</span><span>C${{record.cell_label}} N${{record.nucleus_label}}</span></div>
      <div class="metrics">
        <div class="metric target-metric">cell <b>${{fmt(record.cell_area_um2, 1)}}</b></div>
        <div class="metric target-metric">nuc <b>${{fmt(record.nucleus_area_um2, 1)}}</b></div>
        <div class="metric target-metric">IOD <b>${{fmt(record.nuc_iod, 1)}}</b></div>
        <div class="metric target-metric">OD <b>${{fmt(record.nuc_mean_od, 3)}}</b></div>
        <div class="metric">QC <b>${{fmt(record.qc_score, 3)}}</b></div>
        <div class="metric">focus <b>${{fmt(record.clarity, 3)}}</b></div>
        <div class="metric">match <b>${{fmt(record.quality_match_distance, 3)}}</b></div>
      </div>
      <div class="status-row">
        <button class="status-btn ${{decision === "keep" ? "active" : ""}}" data-decision="keep">keep</button>
        <button class="status-btn ${{decision === "problem" ? "active" : ""}}" data-decision="problem">problem</button>
        <button class="status-btn ${{decision === "unsure" ? "active" : ""}}" data-decision="unsure">unsure</button>
      </div>
    </div>
  </article>`;
}}
function render() {{
  setMaskCss();
  syncControls();
  syncUrl();
  const rows = REVIEW_RECORDS.filter(matches);
  const autoSelectedIds = computeAutoSelectedIds();
  grid.innerHTML = rows.map(record => card(record, autoSelectedIds)).join("");
  const marked = Object.values(state.decisions).filter(v => v && v.decision).length;
  const problems = Object.values(state.decisions).filter(v => v && v.decision === "problem").length;
  const defaultText = DEFAULT_DECISION ? ` | unmarked: ${{DEFAULT_DECISION}} | ${{problems}} problems` : "";
  const scopedAutoRows = REVIEW_RECORDS.filter(record =>
    record.review_cohort === AUTO_FILL_COHORT && (state.species === "all" || record.species === state.species)
  );
  const autoSpecies = new Set(scopedAutoRows.map(record => record.species));
  const autoSelected = scopedAutoRows.filter(record => resolvedDecision(record, autoSelectedIds) === "keep").length;
  const autoTarget = autoSpecies.size * AUTO_FILL_TARGET_PER_SPECIES;
  const autoText = AUTO_FILL_TARGET_PER_SPECIES > 0 && autoSpecies.size > 0
    ? ` | auto-selected ${{autoSelected}}/${{autoTarget}} target`
    : "";
  summary.textContent = `${{rows.length}} visible | ${{marked}} explicit marks${{defaultText}}${{autoText}} | cell mask: ${{state.showCell ? "on" : "off"}} | nucleus mask: ${{state.showNucleus ? "on" : "off"}}`;
}}
speciesSelect.addEventListener("change", event => {{ state.species = event.target.value; render(); }});
searchBox.addEventListener("input", event => {{ state.search = event.target.value; render(); }});
cellToggle.addEventListener("change", event => {{ state.showCell = event.target.checked; render(); }});
nucleusToggle.addEventListener("change", event => {{ state.showNucleus = event.target.checked; render(); }});
grid.addEventListener("click", event => {{
  const btn = event.target.closest("button[data-decision]");
  if (!btn) return;
  const cardEl = btn.closest(".card");
  const id = cardEl.dataset.id;
  const decision = btn.dataset.decision;
  const record = REVIEW_RECORDS.find(r => r.id === id);
  const initialDecision = record.initial_decision || DEFAULT_DECISION || "";
  if (state.decisions[id]?.decision === decision || (decision === initialDecision && !state.decisions[id])) {{
    delete state.decisions[id];
  }} else {{
    state.decisions[id] = {{
      decision,
      species: record.species,
      review_key: record.review_key,
      filename: record.filename,
      tile_name: record.tile_name,
      rank: record.rank
    }};
  }}
  saveDecisions();
  render();
}});
document.getElementById("clearBtn").addEventListener("click", () => {{
  if (!confirm("Clear local validation marks?")) return;
  state.decisions = {{}};
  saveDecisions();
  render();
}});
document.getElementById("exportBtn").addEventListener("click", () => {{
  const header = ["id","decision","decision_source","species","rank","review_key","filename","tile_name"];
  const lines = [header.join(",")];
  const autoSelectedIds = computeAutoSelectedIds();
  const exportableRecords = AUTO_FILL_TARGET_PER_SPECIES > 0
    ? REVIEW_RECORDS.filter(record => resolvedDecision(record, autoSelectedIds))
    : REVIEW_RECORDS.filter(record => resolvedDecision(record));
  const rowsToExport = (EXPORT_ALL_RECORDS || EXPORT_INITIAL_RECORDS || AUTO_FILL_TARGET_PER_SPECIES > 0)
    ? exportableRecords.map(record => [record.id, resolvedDecision(record, autoSelectedIds), decisionSource(record, autoSelectedIds), record.species, record.rank, record.review_key, record.filename, record.tile_name])
    : Object.entries(state.decisions).map(([id, value]) => [id, value.decision, "explicit", value.species, value.rank, value.review_key, value.filename, value.tile_name]);
  rowsToExport.forEach(row => {{
    const csvRow = row
      .map(v => `"${{String(v ?? "").replaceAll('"', '""')}}"`);
    lines.push(csvRow.join(","));
  }});
  const blob = new Blob([lines.join("\\n") + "\\n"], {{type: "text/csv"}});
  const url = URL.createObjectURL(blob);
  const a = document.createElement("a");
  a.href = url;
  a.download = EXPORT_FILENAME;
  a.click();
  URL.revokeObjectURL(url);
}});
speciesOptions();
applyUrlState();
setMaskCss();
render();
</script>
</body>
</html>
""",
        encoding="utf-8",
    )
    return index_path


def write_summary(output_dir: Path, records: list[dict[str, Any]], panel: str, selected_pairs_csv: Path) -> None:
    summary = {
        "panel": panel,
        "selected_pairs_csv": str(selected_pairs_csv.resolve()),
        "output_dir": str(output_dir.resolve()),
        "index_html": str((output_dir / "index.html").resolve()),
        "n_records": len(records),
        "species_counts": {
            species: int(sum(1 for record in records if record["species"] == species))
            for species in sorted({record["species"] for record in records})
        },
    }
    (output_dir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    records = build_records(args)
    index_path = write_index(
        args.output_dir,
        records,
        panel=args.panel,
        title=args.title,
        storage_key=args.storage_key,
        default_decision=args.default_decision,
        blind_target_metrics=args.blind_target_metrics,
        common_support_only=args.common_support_only,
        sort_by_match_distance=args.sort_by_match_distance,
    )
    write_summary(args.output_dir, records, args.panel, args.selected_pairs_csv)
    print(
        json.dumps(
            {
                "index_html": str(index_path.resolve()),
                "n_records": len(records),
                "panel": args.panel,
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
