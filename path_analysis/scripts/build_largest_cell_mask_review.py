#!/usr/bin/env python3
"""Build and freeze a literal-largest, human-audited linked-cell review set.

The review queue is ranked only by cell area after fixed hard eligibility
gates. A freeze requires 50 human ``keep`` decisions per species and requires
all higher-ranked rows to have an explicit keep or problem decision.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[2]
RUN_ROOT = PROJECT_ROOT.parent / "cellprofiler_test" / "output" / "runs" / "mixed_cellpose_yolo_full_dataset_v1_bgclean"
DEFAULT_CANDIDATES = RUN_ROOT / "pair_review_postrepair_strict" / "triage" / "all_scored_masks.csv"
DEFAULT_OUTPUT_DIR = PROJECT_ROOT / "path_analysis" / "data" / "external" / "derived" / "largest_cell_mask_review"
PANEL = "literal_largest_cell_mask_review"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidates-csv", type=Path, default=DEFAULT_CANDIDATES)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--target-per-species", type=int, default=50)
    parser.add_argument("--review-depth", type=int, default=100, help="Ranked fallbacks per species to render.")
    parser.add_argument("--decisions-csv", type=Path, default=None, help="Viewer export; supplying it attempts a freeze.")
    return parser.parse_args()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def normalize_species(value: Any) -> str:
    text = "" if value is None else str(value).strip()
    if not text or text.lower() == "nan":
        return ""
    if text.startswith("Desmognathus "):
        return "D. " + text.split(" ", 1)[1]
    return text if text.startswith("D. ") else "D. " + text


def bool_series(series: pd.Series) -> pd.Series:
    numeric = pd.to_numeric(series, errors="coerce")
    text = series.fillna("").astype(str).str.strip().str.lower()
    return numeric.fillna(0).ne(0) | text.isin({"true", "t", "yes", "y"})


def prepare_candidates(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(path)
    out = pd.read_csv(path, low_memory=False)
    bool_columns = ["has_cell_match", "physical_pair_ok", "one_to_one_cell", "keep_mask_pair", "keep_strict_core", "cell_edge_touch"]
    numeric_columns = ["cell_area_um2", "nuc_area_um2", "nuc_iod", "nucleus_mask_overlap_label_count", "nucleus_mask_best_overlap_fraction"]
    required = {"species", "review_key", "filename", *bool_columns, *numeric_columns}
    missing = sorted(required - set(out.columns))
    if missing:
        raise ValueError(f"Missing required candidate columns: {missing}")
    out["species"] = out["species"].map(normalize_species)
    for column in bool_columns:
        out[column] = bool_series(out[column])
    for column in numeric_columns:
        out[column] = pd.to_numeric(out[column], errors="coerce")
    source_decision = out.get("decision", pd.Series("", index=out.index)).fillna("").astype(str).str.strip().str.lower()
    shape_suspect = bool_series(out["shape_is_suspect"]) if "shape_is_suspect" in out else pd.Series(False, index=out.index)
    eligible = (
        out["species"].str.len().gt(0)
        & out["review_key"].fillna("").astype(str).str.strip().str.len().gt(0)
        & out["has_cell_match"] & out["physical_pair_ok"] & out["one_to_one_cell"]
        & out["keep_mask_pair"] & out["keep_strict_core"] & ~out["cell_edge_touch"]
        & out["cell_area_um2"].gt(0) & out["nuc_area_um2"].gt(0) & out["nuc_iod"].gt(0)
        & out["nucleus_mask_overlap_label_count"].le(1)
        & out["nucleus_mask_best_overlap_fraction"].ge(0.98)
        & ~shape_suspect & ~source_decision.isin({"discard", "nucleus_only", "maybe"})
    )
    return out.loc[eligible].copy()


def build_review_queue(candidates: pd.DataFrame, *, target_per_species: int, review_depth: int) -> pd.DataFrame:
    if target_per_species < 1 or review_depth < target_per_species:
        raise ValueError("review_depth must be at least a positive target_per_species")
    queues: list[pd.DataFrame] = []
    for species, group in candidates.groupby("species", sort=True):
        ranked = group.sort_values(["cell_area_um2", "review_key"], ascending=[False, True], kind="mergesort").copy()
        if len(ranked) < target_per_species:
            raise ValueError(f"{species} has {len(ranked)} hard-eligible rows; need {target_per_species}.")
        queued = ranked.head(review_depth).copy()
        queued["panel"] = PANEL
        queued["selection_rank"] = np.arange(1, len(queued) + 1)
        queued["selection_metric"] = "cell_area_um2_descending_after_fixed_hard_eligibility"
        queues.append(queued)
    return pd.concat(queues, ignore_index=True) if queues else pd.DataFrame()


def load_decisions(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(path)
    decisions = pd.read_csv(path, low_memory=False)
    required = {"species", "review_key", "decision"}
    missing = sorted(required - set(decisions.columns))
    if missing:
        raise ValueError(f"Decision export is missing columns: {missing}")
    out = decisions[["species", "review_key", "decision"]].copy()
    out["species"] = out["species"].map(normalize_species)
    out["review_key"] = out["review_key"].fillna("").astype(str).str.strip()
    out["decision"] = out["decision"].fillna("").astype(str).str.strip().str.lower()
    out = out.loc[out["species"].str.len().gt(0) & out["review_key"].str.len().gt(0)].drop_duplicates(
        ["species", "review_key"], keep="last"
    )
    invalid = sorted(set(out["decision"]) - {"keep", "problem", "unsure"})
    if invalid:
        raise ValueError(f"Unsupported decisions: {invalid}")
    return out


def freeze_queue(queue: pd.DataFrame, decisions: pd.DataFrame, *, target_per_species: int) -> tuple[pd.DataFrame, pd.DataFrame]:
    # ``all_scored_masks.csv`` has a historical model decision column. It is
    # provenance only here; the exported human decision is the freeze authority.
    merged = queue.drop(columns=["decision"], errors="ignore").merge(
        decisions,
        on=["species", "review_key"],
        how="left",
        validate="one_to_one",
    )
    merged["decision"] = merged["decision"].fillna("")
    frozen_groups: list[pd.DataFrame] = []
    audit_groups: list[pd.DataFrame] = []
    failures: list[str] = []
    for species, group in merged.groupby("species", sort=True):
        group = group.sort_values("selection_rank", kind="mergesort").copy()
        keep = group.loc[group["decision"].eq("keep")]
        if len(keep) < target_per_species:
            failures.append(f"{species}: only {len(keep)} keep decisions; need {target_per_species}")
            continue
        cutoff = int(keep.iloc[target_per_species - 1]["selection_rank"])
        blockers = group.loc[group["selection_rank"].le(cutoff) & ~group["decision"].isin({"keep", "problem"})]
        if not blockers.empty:
            failures.append(f"{species}: ranks through {cutoff} are not fully reviewed: " + ",".join(blockers["selection_rank"].astype(str)))
            continue
        group["freeze_cutoff_rank"] = cutoff
        group["frozen_selected"] = group["review_key"].isin(keep.head(target_per_species)["review_key"])
        frozen_groups.append(group.loc[group["frozen_selected"]].copy())
        audit_groups.append(group)
    if failures:
        raise ValueError("Freeze blocked:\n- " + "\n- ".join(failures))
    frozen = pd.concat(frozen_groups, ignore_index=True).sort_values(["species", "selection_rank"]).reset_index(drop=True)
    audit = pd.concat(audit_groups, ignore_index=True).sort_values(["species", "selection_rank"]).reset_index(drop=True)
    return frozen, audit


def write_queue(queue: pd.DataFrame, *, output_dir: Path, candidates_path: Path, args: argparse.Namespace) -> Path:
    output_dir.mkdir(parents=True, exist_ok=True)
    queue_path = output_dir / "largest_cell_mask_review_queue.csv.gz"
    queue.to_csv(queue_path, index=False)
    manifest = {
        "purpose": "Literal-largest cell-mask review queue; cell area is the only ranking metric after fixed hard eligibility.",
        "candidates_csv": str(candidates_path.resolve()),
        "candidates_sha256": sha256_file(candidates_path),
        "queue_csv": str(queue_path.resolve()),
        "queue_sha256": sha256_file(queue_path),
        "target_per_species": args.target_per_species,
        "review_depth": args.review_depth,
        "species_counts": queue.groupby("species").size().to_dict(),
    }
    (output_dir / "review_queue_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    return queue_path


def write_freeze(frozen: pd.DataFrame, audit: pd.DataFrame, *, output_dir: Path, queue_path: Path, decisions_path: Path, target_per_species: int) -> Path:
    counts = frozen.groupby("species").size()
    if not counts.eq(target_per_species).all():
        raise ValueError("Freeze did not contain exactly the requested count per species.")
    frozen_path = output_dir / "frozen_largest_cell_mask_top50.csv.gz"
    audit_path = output_dir / "frozen_largest_cell_mask_audit.csv"
    frozen.to_csv(frozen_path, index=False)
    audit.to_csv(audit_path, index=False)
    manifest_path = output_dir / "freeze_manifest.json"
    manifest_path.write_text(json.dumps({
        "status": "frozen_after_manual_review",
        "selection_rule": "first 50 manual keeps ranked by cell_area_um2 descending after fixed hard eligibility",
        "queue_csv": str(queue_path.resolve()), "queue_sha256": sha256_file(queue_path),
        "decisions_csv": str(decisions_path.resolve()), "decisions_sha256": sha256_file(decisions_path),
        "frozen_csv": str(frozen_path.resolve()), "frozen_sha256": sha256_file(frozen_path),
        "audit_csv": str(audit_path.resolve()), "target_per_species": target_per_species,
        "species_counts": counts.to_dict(),
    }, indent=2) + "\n", encoding="utf-8")
    return manifest_path


def main() -> None:
    args = parse_args()
    queue = build_review_queue(prepare_candidates(args.candidates_csv), target_per_species=args.target_per_species, review_depth=args.review_depth)
    queue_path = write_queue(queue, output_dir=args.output_dir, candidates_path=args.candidates_csv, args=args)
    result: dict[str, Any] = {"queue_csv": str(queue_path.resolve()), "rows": len(queue), "species": queue["species"].nunique()}
    if args.decisions_csv:
        frozen, audit = freeze_queue(queue, load_decisions(args.decisions_csv), target_per_species=args.target_per_species)
        result["freeze_manifest"] = str(write_freeze(frozen, audit, output_dir=args.output_dir, queue_path=queue_path, decisions_path=args.decisions_csv, target_per_species=args.target_per_species).resolve())
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
