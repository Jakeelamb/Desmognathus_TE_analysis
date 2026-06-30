#!/usr/bin/env python3
"""
Audit whether core TE and LTR-history variables are dominated by assembly quality.
"""

from __future__ import annotations

import gzip
import hashlib
import math
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
from scipy import stats

import sys

sys.path.insert(0, str(Path(__file__).parent.parent))
from config import PROJECT_ROOT  # noqa: E402


MANIFEST_FILE = PROJECT_ROOT / "results" / "data" / "ltr_age" / "genome_assembly_manifest.csv"
TE_FEATURE_FILE = PROJECT_ROOT / "path_analysis" / "data" / "derived" / "te_path_features.csv"
LTR_HISTORY_FILE = PROJECT_ROOT / "path_analysis" / "data" / "derived" / "ltr_history_features.csv"
OUTPUT_DIR = PROJECT_ROOT / "results" / "data" / "assembly_bias"
REPORT_DIR = PROJECT_ROOT / "results" / "reports"
NOTE_FILE = REPORT_DIR / "ASSEMBLY_BIAS_AUDIT.md"


TARGET_FEATURES = [
    "ltr_line_logratio",
    "order_pielou",
    "weighted_te_divergence_p90",
    "weighted_te_deletions_p90",
    "ectopic_log10_mean_ratio",
    "ltr_history_median_k2p_distance",
    "ltr_history_n_pairs_estimated",
]

ASSEMBLY_METRICS = [
    "assembly_total_bp",
    "assembly_n_contigs",
    "assembly_n50_bp",
    "assembly_n90_bp",
    "assembly_longest_contig_bp",
    "assembly_n_base_fraction",
    "assembly_top10_bp_fraction",
    "assembly_n_contigs_ge_100kb",
    "assembly_local_bytes_gz",
]


def sha256_for_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8192), b""):
            digest.update(chunk)
    return digest.hexdigest()


def canonical_species(value: object) -> str | pd.NA:
    if value is None or pd.isna(value):
        return pd.NA
    text = str(value).strip()
    if not text:
        return pd.NA
    if text.startswith("D. "):
        text = text[3:]
    elif text.startswith("D."):
        text = text[2:]
    return text.strip()


def fasta_lengths(path: Path) -> tuple[list[int], int]:
    lengths: list[int] = []
    seq_len = 0
    n_bases = 0
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="ascii", errors="ignore") as handle:
        for line in handle:
            if line.startswith(">"):
                if seq_len:
                    lengths.append(seq_len)
                seq_len = 0
                continue
            seq = line.strip().upper()
            if not seq:
                continue
            seq_len += len(seq)
            n_bases += seq.count("N")
    if seq_len:
        lengths.append(seq_len)
    return lengths, n_bases


def nx(lengths: Iterable[int], threshold: float) -> tuple[float, float]:
    values = sorted((int(v) for v in lengths if int(v) > 0), reverse=True)
    if not values:
        return math.nan, math.nan
    total = sum(values)
    cutoff = total * threshold
    running = 0
    for i, value in enumerate(values, start=1):
        running += value
        if running >= cutoff:
            return float(value), float(i)
    return float(values[-1]), float(len(values))


def compute_assembly_metrics(path: Path) -> dict[str, float]:
    lengths, n_bases = fasta_lengths(path)
    if not lengths:
        raise ValueError(f"No contigs found in FASTA: {path}")
    total_bp = float(sum(lengths))
    n50_bp, l50_count = nx(lengths, 0.50)
    n90_bp, l90_count = nx(lengths, 0.90)
    top10_bp = float(sum(sorted(lengths, reverse=True)[:10]))
    return {
        "assembly_total_bp": total_bp,
        "assembly_n_contigs": float(len(lengths)),
        "assembly_longest_contig_bp": float(max(lengths)),
        "assembly_mean_contig_bp": float(np.mean(lengths)),
        "assembly_median_contig_bp": float(np.median(lengths)),
        "assembly_n50_bp": n50_bp,
        "assembly_l50_count": l50_count,
        "assembly_n90_bp": n90_bp,
        "assembly_l90_count": l90_count,
        "assembly_n_base_fraction": float(n_bases / total_bp) if total_bp > 0 else math.nan,
        "assembly_top10_bp_fraction": float(top10_bp / total_bp) if total_bp > 0 else math.nan,
        "assembly_n_contigs_ge_100kb": float(sum(1 for x in lengths if x >= 100_000)),
    }


def load_inputs() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    manifest = pd.read_csv(MANIFEST_FILE)
    manifest["species"] = manifest["species"].map(canonical_species)
    te = pd.read_csv(TE_FEATURE_FILE)
    te["species"] = te["species"].map(canonical_species)
    ltr = pd.read_csv(LTR_HISTORY_FILE)
    ltr["species"] = ltr["species"].map(canonical_species)
    return manifest, te, ltr


def build_metrics_table(manifest: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for row in manifest.itertuples(index=False):
        if not row.local_exists:
            continue
        fasta_path = PROJECT_ROOT / str(row.local_path)
        metrics = compute_assembly_metrics(fasta_path)
        metrics.update(
            {
                "species": row.species,
                "gca_accession": row.gca_accession,
                "assembly_local_bytes_gz": float(row.local_bytes) if pd.notna(row.local_bytes) else math.nan,
                "assembly_status": row.status,
                "assembly_manifest_source_id": "repo_genome_assembly_manifest",
                "assembly_manifest_source_path": str(MANIFEST_FILE.relative_to(PROJECT_ROOT)),
                "assembly_manifest_source_sha256": sha256_for_file(MANIFEST_FILE),
            }
        )
        rows.append(metrics)
    return pd.DataFrame(rows).sort_values("species").reset_index(drop=True)


def adjust_pvalues_bh(pvals: pd.Series) -> pd.Series:
    values = pvals.to_numpy(dtype=float)
    n = len(values)
    order = np.argsort(values)
    ranks = np.empty(n, dtype=float)
    ranks[order] = np.arange(1, n + 1, dtype=float)
    adjusted = values * n / ranks
    adjusted[order[::-1]] = np.minimum.accumulate(adjusted[order[::-1]])
    adjusted = np.clip(adjusted, 0, 1)
    return pd.Series(adjusted, index=pvals.index)


def analyze_bias(merged: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    corr_rows = []
    model_rows = []

    for feature in TARGET_FEATURES:
        if feature not in merged.columns:
            continue
        for metric in ASSEMBLY_METRICS:
            if metric not in merged.columns:
                continue
            subset = merged[["species", feature, metric]].dropna().copy()
            n = len(subset)
            if n < 8:
                continue
            rho, pval = stats.spearmanr(subset[metric], subset[feature])
            corr_rows.append(
                {
                    "feature": feature,
                    "assembly_metric": metric,
                    "n": int(n),
                    "spearman_rho": float(rho),
                    "spearman_p": float(pval),
                }
            )

            x = subset[metric].astype(float).to_numpy()
            y = subset[feature].astype(float).to_numpy()
            if np.allclose(np.std(x), 0) or np.allclose(np.std(y), 0):
                continue
            slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)
            model_rows.append(
                {
                    "feature": feature,
                    "assembly_metric": metric,
                    "n": int(n),
                    "slope": float(slope),
                    "intercept": float(intercept),
                    "r_squared": float(r_value**2),
                    "linregress_p": float(p_value),
                    "std_err": float(std_err),
                }
            )

    corr_df = pd.DataFrame(corr_rows)
    if not corr_df.empty:
        corr_df["spearman_p_adjusted_bh"] = adjust_pvalues_bh(corr_df["spearman_p"])
        corr_df["strong_effect"] = corr_df["spearman_rho"].abs().ge(0.6)
        corr_df["nominal_signal"] = corr_df["spearman_p"].lt(0.05)
        corr_df["bh_signal"] = corr_df["spearman_p_adjusted_bh"].lt(0.05)
        corr_df = corr_df.sort_values(
            ["bh_signal", "nominal_signal", "strong_effect", "spearman_p", "feature"],
            ascending=[False, False, False, True, True],
        ).reset_index(drop=True)

    model_df = pd.DataFrame(model_rows)
    if not model_df.empty:
        model_df["linregress_p_adjusted_bh"] = adjust_pvalues_bh(model_df["linregress_p"])
        model_df = model_df.sort_values(
            ["linregress_p_adjusted_bh", "r_squared", "feature"],
            ascending=[True, False, True],
        ).reset_index(drop=True)

    return corr_df, model_df


def write_note(metrics_df: pd.DataFrame, corr_df: pd.DataFrame, model_df: pd.DataFrame) -> None:
    feature_counts = {
        feature: int(metrics_df[feature].notna().sum()) for feature in TARGET_FEATURES if feature in metrics_df.columns
    }
    strongest = corr_df.head(10).copy() if not corr_df.empty else pd.DataFrame()
    strong_hits = corr_df[corr_df["bh_signal"] | corr_df["strong_effect"]].copy() if not corr_df.empty else pd.DataFrame()

    lines = [
        "# Assembly Bias Audit",
        "",
        "This note asks whether the core TE and paired-LTR history variables are",
        "dominated by assembly quality / fragmentation metrics rather than biology.",
        "",
        "## Source Basis",
        "",
        "- `results/data/ltr_age/genome_assembly_manifest.csv`",
        "- `path_analysis/data/derived/te_path_features.csv`",
        "- `path_analysis/data/derived/ltr_history_features.csv`",
        "",
        "## Coverage",
        "",
        f"- Assemblies with local verified FASTA files audited: `{len(metrics_df)}`",
    ]

    for feature, count in feature_counts.items():
        lines.append(f"- `{feature}` non-missing species: `{count}`")

    lines.extend(
        [
            "",
            "## Bottom Line",
            "",
        ]
    )

    bh_hits = corr_df[corr_df["bh_signal"]].copy() if not corr_df.empty else pd.DataFrame()

    if strong_hits.empty:
        lines.append(
            "- No assembly-quality metric reached the current strong-signal screen after BH correction or with an absolute Spearman rho >= 0.6."
        )
        lines.append(
            "- On the present audit, the core TE and paired-LTR history variables do not look dominated by assembly fragmentation or assembly size."
        )
    else:
        if len(bh_hits) == 1 and bh_hits.iloc[0]["feature"] == "ltr_line_logratio":
            row = bh_hits.iloc[0]
            lines.append(
                "- Most core TE and paired-LTR history variables do not show a strong assembly-bias signal in this screen."
            )
            lines.append(
                f"- The main exception is `ltr_line_logratio`, which is negatively associated with "
                f"`{row['assembly_metric']}` (rho `{row['spearman_rho']:.3f}`, BH p `{row['spearman_p_adjusted_bh']:.4g}`)."
            )
            lines.append(
                "- Treat the LTR:LINE balance feature as technically sensitive enough to require extra caution and cross-method validation, but do not downgrade the full TE feature set wholesale."
            )
        else:
            lines.append(
                "- Some assembly-quality relationships were strong enough to warrant caution. Inspect `results/data/assembly_bias/assembly_bias_correlations.csv` before promoting every TE variable equally."
            )

    lines.extend(
        [
            "",
            "## Strongest Screened Associations",
            "",
        ]
    )

    if strongest.empty:
        lines.append("- No analyzable feature/metric pairs were available.")
    else:
        for row in strongest.itertuples(index=False):
            lines.append(
                f"- `{row.feature}` vs `{row.assembly_metric}`: rho `{row.spearman_rho:.3f}`, "
                f"raw p `{row.spearman_p:.4g}`, BH p `{row.spearman_p_adjusted_bh:.4g}`, n `{row.n}`"
            )

    lines.extend(
        [
            "",
            "## Files",
            "",
            "- `results/data/assembly_bias/assembly_quality_metrics.csv`",
            "- `results/data/assembly_bias/assembly_bias_correlations.csv`",
            "- `results/data/assembly_bias/assembly_bias_models.csv`",
        ]
    )

    NOTE_FILE.write_text("\n".join(lines) + "\n", encoding="ascii")


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    REPORT_DIR.mkdir(parents=True, exist_ok=True)

    manifest, te, ltr = load_inputs()
    metrics_df = build_metrics_table(manifest)
    merged = metrics_df.merge(te, on="species", how="left", validate="one_to_one")
    merged = merged.merge(ltr, on="species", how="left", suffixes=("", "_ltr"), validate="one_to_one")
    corr_df, model_df = analyze_bias(merged)

    metrics_df.to_csv(OUTPUT_DIR / "assembly_quality_metrics.csv", index=False)
    corr_df.to_csv(OUTPUT_DIR / "assembly_bias_correlations.csv", index=False)
    model_df.to_csv(OUTPUT_DIR / "assembly_bias_models.csv", index=False)
    write_note(merged, corr_df, model_df)

    print(f"Wrote {OUTPUT_DIR / 'assembly_quality_metrics.csv'}")
    print(f"Wrote {OUTPUT_DIR / 'assembly_bias_correlations.csv'}")
    print(f"Wrote {OUTPUT_DIR / 'assembly_bias_models.csv'}")
    print(f"Wrote {NOTE_FILE}")


if __name__ == "__main__":
    main()
