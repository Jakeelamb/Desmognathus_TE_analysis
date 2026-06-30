#!/usr/bin/env python3
"""
Audit concordance between read-based dnaPipeTE summaries and raw RepeatMasker summaries.
"""

from __future__ import annotations

import hashlib
import math
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

import sys

sys.path.insert(0, str(Path(__file__).parent.parent))
from config import PROJECT_ROOT  # noqa: E402


DNAPIPETE_ORDER_FILE = PROJECT_ROOT / "results" / "data" / "dnaPipeTE_order_breakdown.csv"
DNAPIPETE_SUPERFAMILY_FILE = PROJECT_ROOT / "results" / "data" / "dnaPipeTE_superfamily_breakdown.csv"
REPEATMASKER_RAW_FILE = PROJECT_ROOT / "results" / "data" / "merged_repeatmasker_data.csv"
OUTPUT_DIR = PROJECT_ROOT / "results" / "data" / "cross_method"
REPORT_DIR = PROJECT_ROOT / "results" / "reports"
NOTE_FILE = REPORT_DIR / "CROSS_METHOD_CONCORDANCE_AUDIT.md"


ORDER_COLUMNS = ["DIRS", "Helitron", "LINE", "LTR", "Maverick", "PLE", "SINE", "TIR", "YR"]
SUPERFAMILY_COLUMNS = [
    "Academ",
    "CACTA",
    "Chapaev",
    "DIRS",
    "Dada",
    "EnSpm",
    "Ginger",
    "Gypsy",
    "Helitron",
    "Jockey",
    "L1",
    "MULE",
    "Maverick",
    "Mutator",
    "P",
    "PIF-Harbinger",
    "Penelope",
    "PiggyBac",
    "Tc1-mariner",
    "hAT",
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


def standardize_wide(df: pd.DataFrame) -> pd.DataFrame:
    species_col = df.columns[0]
    df = df.rename(columns={species_col: "species"}).copy()
    df["species"] = df["species"].map(canonical_species)
    return df.dropna(subset=["species"]).drop_duplicates(subset=["species"]).reset_index(drop=True)


def map_repeatmasker_order(repeat_class: object) -> str | None:
    if repeat_class is None or pd.isna(repeat_class):
        return None
    text = str(repeat_class).strip()
    if not text:
        return None
    parts = text.split("/", 1)
    prefix = parts[0]
    prefix_upper = prefix.upper()
    if prefix_upper == "DNA":
        return "TIR"
    mapping = {
        "LTR": "LTR",
        "LINE": "LINE",
        "SINE": "SINE",
        "RC": "Helitron",
        "HELITRON": "Helitron",
        "DIRS": "DIRS",
        "PLE": "PLE",
        "PENELOPE": "PLE",
        "MAVERICK": "Maverick",
        "YR": "YR",
    }
    return mapping.get(prefix_upper)


def map_repeatmasker_superfamily(repeat_class: object) -> str | None:
    if repeat_class is None or pd.isna(repeat_class):
        return None
    text = str(repeat_class).strip()
    if not text or "/" not in text:
        return None
    suffix = text.split("/", 1)[1].strip()
    mapping = {
        "ACADEM": "Academ",
        "CACTA": "CACTA",
        "CHAPAEV": "Chapaev",
        "DIRS": "DIRS",
        "DADA": "Dada",
        "ENSPM": "EnSpm",
        "GINGER": "Ginger",
        "GYPSY": "Gypsy",
        "HELITRON": "Helitron",
        "JOCKEY": "Jockey",
        "L1": "L1",
        "MULE": "MULE",
        "MAVERICK": "Maverick",
        "MUTATOR": "Mutator",
        "P": "P",
        "PIGGYBAC": "PiggyBac",
        "PIF-HARBINGER": "PIF-Harbinger",
        "HARBINGER": "PIF-Harbinger",
        "PENELOPE": "Penelope",
        "TC1-MARINER": "Tc1-mariner",
        "TCMAR-TC1": "Tc1-mariner",
        "MARINER": "Tc1-mariner",
        "HAT": "hAT",
    }
    return mapping.get(suffix.upper())


def build_repeatmasker_summaries() -> tuple[pd.DataFrame, pd.DataFrame]:
    order_totals: dict[tuple[str, str], float] = defaultdict(float)
    superfamily_totals: dict[tuple[str, str], float] = defaultdict(float)

    usecols = ["Desmognathus_Species", "repeat_class", "query_start", "query_end"]
    for chunk in pd.read_csv(REPEATMASKER_RAW_FILE, usecols=usecols, chunksize=250_000):
        chunk["species"] = chunk["Desmognathus_Species"].map(canonical_species)
        chunk["bp"] = (
            pd.to_numeric(chunk["query_end"], errors="coerce")
            - pd.to_numeric(chunk["query_start"], errors="coerce")
            + 1
        )
        chunk = chunk.dropna(subset=["species", "bp"])
        chunk = chunk[chunk["bp"] > 0]
        chunk["order"] = chunk["repeat_class"].map(map_repeatmasker_order)
        chunk["superfamily"] = chunk["repeat_class"].map(map_repeatmasker_superfamily)

        order_grouped = (
            chunk.dropna(subset=["order"])
            .groupby(["species", "order"], dropna=False)["bp"]
            .sum()
            .items()
        )
        for key, value in order_grouped:
            order_totals[key] += float(value)

        superfamily_grouped = (
            chunk.dropna(subset=["superfamily"])
            .groupby(["species", "superfamily"], dropna=False)["bp"]
            .sum()
            .items()
        )
        for key, value in superfamily_grouped:
            superfamily_totals[key] += float(value)

    order_df = pd.DataFrame(
        [{"species": s, "order": o, "bp": bp} for (s, o), bp in order_totals.items()]
    )
    superfamily_df = pd.DataFrame(
        [{"species": s, "superfamily": sf, "bp": bp} for (s, sf), bp in superfamily_totals.items()]
    )

    order_wide = (
        order_df.pivot(index="species", columns="order", values="bp")
        .reindex(columns=ORDER_COLUMNS, fill_value=0.0)
        .fillna(0.0)
    )
    order_wide = order_wide.div(order_wide.sum(axis=1), axis=0) * 100
    order_wide = order_wide.reset_index()

    superfamily_wide = (
        superfamily_df.pivot(index="species", columns="superfamily", values="bp")
        .reindex(columns=SUPERFAMILY_COLUMNS, fill_value=0.0)
        .fillna(0.0)
    )
    mapped_fraction = (
        superfamily_df.groupby("species", dropna=False)["bp"].sum()
        / order_df.groupby("species", dropna=False)["bp"].sum()
    ).rename("repeatmasker_shared_superfamily_fraction")
    superfamily_wide = superfamily_wide.div(superfamily_wide.sum(axis=1).replace(0, np.nan), axis=0) * 100
    superfamily_wide = superfamily_wide.fillna(0.0).reset_index()
    superfamily_wide = superfamily_wide.merge(mapped_fraction.reset_index(), on="species", how="left")

    return order_wide, superfamily_wide


def pielou_evenness(df: pd.DataFrame) -> pd.Series:
    values = df.to_numpy(dtype=float)
    values = np.where(values < 0, np.nan, values)
    row_sums = values.sum(axis=1, keepdims=True)
    probs = np.divide(values, row_sums, out=np.zeros_like(values), where=row_sums > 0)
    with np.errstate(divide="ignore", invalid="ignore"):
        shannon = -(probs * np.log(probs, where=probs > 0, out=np.zeros_like(probs))).sum(axis=1)
    richness = (values > 0).sum(axis=1)
    evenness = np.divide(shannon, np.log(richness), out=np.full_like(shannon, np.nan), where=richness > 1)
    return pd.Series(evenness, index=df.index)


def derive_order_features(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    positive = out[ORDER_COLUMNS].where(out[ORDER_COLUMNS] > 0).stack()
    pseudocount = float(positive.min() / 2) if not positive.empty else 1e-6
    out["ltr_line_logratio"] = np.log((out["LTR"] + pseudocount) / (out["LINE"] + pseudocount))
    retro = out[["LTR", "LINE", "SINE", "PLE"]].sum(axis=1)
    dna = out[["TIR", "Maverick", "Helitron"]].sum(axis=1)
    out["retro_dna_logratio"] = np.log((retro + pseudocount) / (dna + pseudocount))
    out["order_pielou"] = pielou_evenness(out[ORDER_COLUMNS])
    return out


def compare_wide(
    left: pd.DataFrame,
    right: pd.DataFrame,
    columns: list[str],
    level_name: str,
) -> pd.DataFrame:
    rows = []
    merged = left.merge(right, on="species", how="inner", suffixes=("_dna", "_rm"))
    for col in columns:
        left_col = f"{col}_dna"
        right_col = f"{col}_rm"
        subset = merged[["species", left_col, right_col]].dropna()
        if len(subset) < 8:
            continue
        if subset[left_col].nunique() < 2 or subset[right_col].nunique() < 2:
            continue
        spearman_rho, spearman_p = stats.spearmanr(subset[left_col], subset[right_col])
        pearson_r, pearson_p = stats.pearsonr(subset[left_col], subset[right_col])
        rows.append(
            {
                "level": level_name,
                "variable": col,
                "n_species": int(len(subset)),
                "dna_mean": float(subset[left_col].mean()),
                "repeatmasker_mean": float(subset[right_col].mean()),
                "spearman_rho": float(spearman_rho),
                "spearman_p": float(spearman_p),
                "pearson_r": float(pearson_r),
                "pearson_p": float(pearson_p),
            }
        )
    return pd.DataFrame(rows)


def build_species_disagreement(order_dna: pd.DataFrame, order_rm: pd.DataFrame) -> pd.DataFrame:
    merged = order_dna.merge(order_rm, on="species", how="inner", suffixes=("_dna", "_rm"))
    rows = []
    for row in merged.itertuples(index=False):
        dna_vec = np.array([getattr(row, f"{c}_dna") for c in ORDER_COLUMNS], dtype=float)
        rm_vec = np.array([getattr(row, f"{c}_rm") for c in ORDER_COLUMNS], dtype=float)
        bray = np.abs(dna_vec - rm_vec).sum() / (dna_vec.sum() + rm_vec.sum()) if (dna_vec.sum() + rm_vec.sum()) > 0 else math.nan
        rows.append(
            {
                "species": row.species,
                "order_bray_disagreement": float(bray),
                "ltr_line_logratio_abs_diff": float(abs(row.ltr_line_logratio_dna - row.ltr_line_logratio_rm)),
                "order_pielou_abs_diff": float(abs(row.order_pielou_dna - row.order_pielou_rm)),
                "retro_dna_logratio_abs_diff": float(abs(row.retro_dna_logratio_dna - row.retro_dna_logratio_rm)),
            }
        )
    return pd.DataFrame(rows).sort_values(
        ["order_bray_disagreement", "ltr_line_logratio_abs_diff"],
        ascending=[False, False],
    ).reset_index(drop=True)


def write_note(order_df: pd.DataFrame, superfamily_df: pd.DataFrame, disagreement_df: pd.DataFrame) -> None:
    feature_rows = order_df[order_df["variable"].isin(["ltr_line_logratio", "order_pielou", "retro_dna_logratio"])].copy()
    strongest = order_df.sort_values("spearman_p").head(10)
    outliers = disagreement_df.head(5)

    lines = [
        "# Cross-Method Concordance Audit",
        "",
        "This note compares read-based dnaPipeTE summaries against raw RepeatMasker",
        "classification summaries built directly from `merged_repeatmasker_data.csv`.",
        "",
        "## Source Basis",
        "",
        "- `results/data/dnaPipeTE_order_breakdown.csv`",
        "- `results/data/dnaPipeTE_superfamily_breakdown.csv`",
        "- `results/data/merged_repeatmasker_data.csv`",
        "",
        "## Bottom Line",
        "",
    ]

    if not feature_rows.empty and (feature_rows["spearman_rho"].abs() >= 0.6).any():
        lines.append(
            "- The main order-level TE features show at least moderate cross-method agreement, which supports treating them as biological rather than purely method-specific summaries."
        )
    else:
        lines.append(
            "- Cross-method agreement is mixed, so the main TE features should be interpreted with method-specific caution until the strongest disagreements are explained."
        )

    lines.extend(["", "## Strongest Order-Level Concordance", ""])
    for row in strongest.itertuples(index=False):
        lines.append(
            f"- `{row.variable}`: Spearman rho `{row.spearman_rho:.3f}`, Pearson r `{row.pearson_r:.3f}`, n `{row.n_species}`"
        )

    lines.extend(["", "## Highest Species-Level Disagreement", ""])
    for row in outliers.itertuples(index=False):
        lines.append(
            f"- `{row.species}`: order Bray disagreement `{row.order_bray_disagreement:.3f}`, "
            f"|LTR:LINE| diff `{row.ltr_line_logratio_abs_diff:.3f}`, "
            f"|Pielou| diff `{row.order_pielou_abs_diff:.3f}`"
        )

    lines.extend(
        [
            "",
            "## Files",
            "",
            "- `results/data/cross_method/cross_method_order_concordance.csv`",
            "- `results/data/cross_method/cross_method_superfamily_concordance.csv`",
            "- `results/data/cross_method/cross_method_species_disagreement.csv`",
        ]
    )
    NOTE_FILE.write_text("\n".join(lines) + "\n", encoding="ascii")


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    REPORT_DIR.mkdir(parents=True, exist_ok=True)

    order_dna = standardize_wide(pd.read_csv(DNAPIPETE_ORDER_FILE))
    superfamily_dna = standardize_wide(pd.read_csv(DNAPIPETE_SUPERFAMILY_FILE))
    order_rm, superfamily_rm = build_repeatmasker_summaries()

    order_dna = derive_order_features(order_dna)
    order_rm = derive_order_features(order_rm)

    order_compare = compare_wide(
        order_dna,
        order_rm,
        ORDER_COLUMNS + ["ltr_line_logratio", "order_pielou", "retro_dna_logratio"],
        "order",
    )
    superfamily_compare = compare_wide(superfamily_dna, superfamily_rm, SUPERFAMILY_COLUMNS, "superfamily")
    disagreement = build_species_disagreement(order_dna, order_rm)

    order_compare["dna_source_path"] = str(DNAPIPETE_ORDER_FILE.relative_to(PROJECT_ROOT))
    order_compare["repeatmasker_source_path"] = str(REPEATMASKER_RAW_FILE.relative_to(PROJECT_ROOT))
    superfamily_compare["dna_source_path"] = str(DNAPIPETE_SUPERFAMILY_FILE.relative_to(PROJECT_ROOT))
    superfamily_compare["repeatmasker_source_path"] = str(REPEATMASKER_RAW_FILE.relative_to(PROJECT_ROOT))
    disagreement["dna_source_path"] = str(DNAPIPETE_ORDER_FILE.relative_to(PROJECT_ROOT))
    disagreement["repeatmasker_source_path"] = str(REPEATMASKER_RAW_FILE.relative_to(PROJECT_ROOT))

    order_compare.to_csv(OUTPUT_DIR / "cross_method_order_concordance.csv", index=False)
    superfamily_compare.to_csv(OUTPUT_DIR / "cross_method_superfamily_concordance.csv", index=False)
    disagreement.to_csv(OUTPUT_DIR / "cross_method_species_disagreement.csv", index=False)
    write_note(order_compare, superfamily_compare, disagreement)

    print(f"Wrote {OUTPUT_DIR / 'cross_method_order_concordance.csv'}")
    print(f"Wrote {OUTPUT_DIR / 'cross_method_superfamily_concordance.csv'}")
    print(f"Wrote {OUTPUT_DIR / 'cross_method_species_disagreement.csv'}")
    print(f"Wrote {NOTE_FILE}")


if __name__ == "__main__":
    main()
