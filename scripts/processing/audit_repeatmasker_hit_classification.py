#!/usr/bin/env python3
"""Audit whether RepeatMasker hits inherited conflicting contig annotations."""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from pathlib import Path
import sys
from typing import Dict

import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))

from scripts.processing.te_classification import classify_repeatmasker_hits


COUNT_FIELDS = (
    "n_hits",
    "aligned_bp",
    "repeatmasker_unmapped_hits",
    "repeatmasker_unmapped_bp",
    "order_comparable_hits",
    "order_comparable_bp",
    "order_mismatch_hits",
    "order_mismatch_bp",
    "superfamily_comparable_hits",
    "superfamily_comparable_bp",
    "superfamily_mismatch_hits",
    "superfamily_mismatch_bp",
)


def _valid_classification(values: pd.Series) -> pd.Series:
    text = values.astype("string").str.strip()
    return values.notna() & text.notna() & ~text.isin(["", "Unclassified"])


def classify_current_combined(chunk: pd.DataFrame) -> pd.DataFrame:
    """Namespace inherited fields, then classify from each native hit label."""

    rename = {}
    for current, namespaced in (
        ("Class", "dnapipete_te_class"),
        ("Order", "dnapipete_order"),
        ("Superfamily", "dnapipete_superfamily"),
    ):
        if namespaced not in chunk.columns and current in chunk.columns:
            rename[current] = namespaced
    return classify_repeatmasker_hits(chunk.rename(columns=rename))


def summarize_chunk(chunk: pd.DataFrame) -> Dict[str, int]:
    """Return additive row and inclusive-base-pair concordance counts."""

    required = {"repeat_class", "query_start", "query_end"}
    missing = sorted(required.difference(chunk.columns))
    if missing:
        raise ValueError(f"RepeatMasker audit chunk is missing columns: {missing}")

    classified = classify_current_combined(chunk)
    start = pd.to_numeric(classified["query_start"], errors="coerce")
    end = pd.to_numeric(classified["query_end"], errors="coerce")
    aligned_bp = (end - start).abs() + 1
    if aligned_bp.isna().any() or (aligned_bp <= 0).any():
        raise ValueError("RepeatMasker audit encountered invalid hit coordinates")

    rm_unmapped = classified["repeatmasker_classification_status"].eq("unmapped")
    rm_order_known = _valid_classification(classified["repeatmasker_order"])
    dna_order_known = _valid_classification(classified["dnapipete_order"])
    order_comparable = rm_order_known & dna_order_known
    order_mismatch = order_comparable & classified[
        "repeatmasker_order"
    ].astype("string").ne(classified["dnapipete_order"].astype("string"))

    rm_superfamily_known = _valid_classification(
        classified["repeatmasker_superfamily"]
    )
    dna_superfamily_known = _valid_classification(classified["dnapipete_superfamily"])
    superfamily_comparable = rm_superfamily_known & dna_superfamily_known
    superfamily_mismatch = superfamily_comparable & classified[
        "repeatmasker_superfamily"
    ].astype("string").ne(classified["dnapipete_superfamily"].astype("string"))

    return {
        "n_hits": int(len(classified)),
        "aligned_bp": int(aligned_bp.sum()),
        "repeatmasker_unmapped_hits": int(rm_unmapped.sum()),
        "repeatmasker_unmapped_bp": int(aligned_bp.loc[rm_unmapped].sum()),
        "order_comparable_hits": int(order_comparable.sum()),
        "order_comparable_bp": int(aligned_bp.loc[order_comparable].sum()),
        "order_mismatch_hits": int(order_mismatch.sum()),
        "order_mismatch_bp": int(aligned_bp.loc[order_mismatch].sum()),
        "superfamily_comparable_hits": int(superfamily_comparable.sum()),
        "superfamily_comparable_bp": int(aligned_bp.loc[superfamily_comparable].sum()),
        "superfamily_mismatch_hits": int(superfamily_mismatch.sum()),
        "superfamily_mismatch_bp": int(aligned_bp.loc[superfamily_mismatch].sum()),
    }


def _rate(numerator: int, denominator: int) -> float:
    return numerator / denominator if denominator else float("nan")


def _portable_source_path(path: Path) -> str:
    try:
        return str(path.resolve().relative_to(PROJECT_ROOT.resolve()))
    except ValueError:
        return str(path)


def audit_file(input_path: Path, chunksize: int = 250_000):
    """Stream the current merged table and return summary plus unmapped labels."""

    totals: Counter = Counter()
    unmapped_hits: Counter = Counter()
    unmapped_bp = defaultdict(int)

    usecols = [
        "repeat_class",
        "Class",
        "Order",
        "Superfamily",
        "query_start",
        "query_end",
    ]
    for chunk in pd.read_csv(
        input_path,
        usecols=usecols,
        chunksize=chunksize,
        low_memory=False,
    ):
        counts = summarize_chunk(chunk)
        totals.update(counts)

        classified = classify_current_combined(chunk)
        start = pd.to_numeric(classified["query_start"], errors="coerce")
        end = pd.to_numeric(classified["query_end"], errors="coerce")
        bp = (end - start).abs() + 1
        unmapped = classified["repeatmasker_classification_status"].eq("unmapped")
        for label, group_index in classified.loc[unmapped].groupby(
            "repeat_class", dropna=False
        ).groups.items():
            label_text = "<missing>" if pd.isna(label) else str(label)
            unmapped_hits[label_text] += len(group_index)
            unmapped_bp[label_text] += int(bp.loc[group_index].sum())

    summary = {field: int(totals[field]) for field in COUNT_FIELDS}
    summary.update(
        {
            "repeatmasker_unmapped_hit_fraction": _rate(
                summary["repeatmasker_unmapped_hits"], summary["n_hits"]
            ),
            "repeatmasker_unmapped_bp_fraction": _rate(
                summary["repeatmasker_unmapped_bp"], summary["aligned_bp"]
            ),
            "order_mismatch_fraction_of_comparable_hits": _rate(
                summary["order_mismatch_hits"], summary["order_comparable_hits"]
            ),
            "order_mismatch_fraction_of_comparable_bp": _rate(
                summary["order_mismatch_bp"], summary["order_comparable_bp"]
            ),
            "superfamily_mismatch_fraction_of_comparable_hits": _rate(
                summary["superfamily_mismatch_hits"],
                summary["superfamily_comparable_hits"],
            ),
            "superfamily_mismatch_fraction_of_comparable_bp": _rate(
                summary["superfamily_mismatch_bp"],
                summary["superfamily_comparable_bp"],
            ),
            "source_path": _portable_source_path(input_path),
            "source_state": "pre_audit_contig_inherited_classification",
        }
    )
    unmapped = pd.DataFrame(
        [
            {
                "repeat_class": label,
                "n_hits": unmapped_hits[label],
                "aligned_bp": unmapped_bp[label],
            }
            for label in sorted(unmapped_hits)
        ]
    )
    return pd.DataFrame([summary]), unmapped


def render_markdown(summary: pd.Series, unmapped: pd.DataFrame) -> str:
    return f"""# RepeatMasker hit-classification audit

The current merged table contains {int(summary.n_hits):,} RepeatMasker hits. Its generic class/order/superfamily columns were inherited from one dnaPipeTE annotation per contig; this audit reclassified each row from its own native `repeat_class` without rewriting the stored table.

- Known-order disagreement: {summary.order_mismatch_fraction_of_comparable_hits:.3%} of comparable hits and {summary.order_mismatch_fraction_of_comparable_bp:.3%} of comparable aligned bp.
- Known-superfamily disagreement: {summary.superfamily_mismatch_fraction_of_comparable_hits:.3%} of comparable hits and {summary.superfamily_mismatch_fraction_of_comparable_bp:.3%} of comparable aligned bp.
- Native RepeatMasker labels not covered by the shared map: {int(summary.repeatmasker_unmapped_hits):,} hits ({summary.repeatmasker_unmapped_hit_fraction:.4%}).
- Aligned length uses inclusive RepeatMasker coordinates: `abs(query_end - query_start) + 1`.

The corrected merge code retains dnaPipeTE classifications under `dnapipete_*`, writes hit-derived fields under `repeatmasker_*`, and keeps generic `Class`, `Order`, and `Superfamily` only as backward-compatible aliases of the hit-derived result. Existing divergence outputs remain pre-audit until regenerated from that corrected table.

Unmapped-label detail is stored in `repeatmasker_unmapped_labels.csv` ({len(unmapped)} labels).
"""


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input",
        type=Path,
        default=PROJECT_ROOT
        / "results"
        / "data"
        / "repeatmasker_detailed_classification_combined.csv",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=PROJECT_ROOT / "plans" / "publication-readiness-deep-audit",
    )
    parser.add_argument("--chunksize", type=int, default=250_000)
    args = parser.parse_args()

    summary, unmapped = audit_file(args.input, chunksize=args.chunksize)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    summary_path = args.output_dir / "repeatmasker_hit_classification_audit.csv"
    unmapped_path = args.output_dir / "repeatmasker_unmapped_labels.csv"
    report_path = args.output_dir / "repeatmasker_hit_classification_audit.md"
    summary.to_csv(summary_path, index=False)
    unmapped.to_csv(unmapped_path, index=False)
    report_path.write_text(
        render_markdown(summary.iloc[0], unmapped), encoding="utf-8"
    )
    print(f"Wrote {summary_path}")
    print(f"Wrote {unmapped_path}")
    print(f"Wrote {report_path}")


if __name__ == "__main__":
    main()
