#!/usr/bin/env python3
"""
LTR age readiness audit.

This script intentionally does not estimate sequence-based LTR insertion ages
from the current local repo state. The previous implementation attempted to
derive insertion ages from RepeatMasker divergence-to-consensus hits on Trinity
assemblies, which is not a defensible substitute for divergence between true
paired 5' and 3' LTR sequences from the same genomic insertion.

The canonical local source of paired-LTR structure is the ectopic master table
written by `scripts/processing/ec.py`. This script audits that paired-element
inventory, records how much age-estimation substrate exists, and reports why
sequence-based age estimation is currently blocked.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict, Iterable, List

import pandas as pd

import sys

sys.path.insert(0, str(Path(__file__).parent.parent))
from config import PROJECT_ROOT, load_lookup_table, paths  # noqa: E402


MASTER_TABLE = paths.results.data / "ectopic_recombination_master.csv"
OUTPUT_DIR = paths.results.data / "ltr_age"
AUDIT_PATH = PROJECT_ROOT / "LTR_AGE_AUDIT.md"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Audit readiness for sequence-based LTR insertion age estimation."
    )
    parser.add_argument(
        "--min-ltr-length",
        type=int,
        default=100,
        help="Minimum LTR length used for readiness counts (default: 100 bp).",
    )
    parser.add_argument(
        "--fail-on-unready",
        action="store_true",
        help="Exit non-zero if genome FASTA files are unavailable for audited species.",
    )
    return parser.parse_args()


def load_master_table() -> pd.DataFrame:
    if not MASTER_TABLE.exists():
        raise FileNotFoundError(f"Missing canonical ectopic master table: {MASTER_TABLE}")

    df = pd.read_csv(MASTER_TABLE, sep="\t")
    required = {
        "species",
        "sequence",
        "element start",
        "element end",
        "lLTR start",
        "lLTR end",
        "lLTR length",
        "rLTR start",
        "rLTR end",
        "rLTR length",
        "Order",
        "Superfamily",
        "Complete",
        "domain_count",
    }
    missing = sorted(required - set(df.columns))
    if missing:
        raise ValueError(f"Ectopic master table is missing required columns: {missing}")
    return df


def find_local_genome_files(gca_values: Iterable[str]) -> Dict[str, str]:
    suffixes = (".fna", ".fa", ".fasta", ".fna.gz", ".fa.gz", ".fasta.gz")
    search_roots = [
        PROJECT_ROOT / "input_data" / "genomes",
        PROJECT_ROOT / "input_data",
        PROJECT_ROOT,
    ]
    found: Dict[str, str] = {}
    for gca in gca_values:
        for root in search_roots:
            if not root.exists():
                continue
            for suffix in suffixes:
                direct = root / f"{gca}{suffix}"
                if direct.exists():
                    found[gca] = str(direct.relative_to(PROJECT_ROOT))
                    break
            if gca in found:
                break
            for suffix in suffixes:
                matches = list(root.glob(f"**/{gca}*{suffix}"))
                if matches:
                    found[gca] = str(matches[0].relative_to(PROJECT_ROOT))
                    break
            if gca in found:
                break
    return found


def build_species_readiness(df: pd.DataFrame, lookup: pd.DataFrame, min_ltr_length: int) -> pd.DataFrame:
    paired = df.copy()
    paired["internal_length"] = paired["rLTR start"] - paired["lLTR end"] - 1
    paired["has_valid_pair_coords"] = (
        paired["lLTR start"].notna()
        & paired["rLTR start"].notna()
        & (paired["lLTR length"] > 0)
        & (paired["rLTR length"] > 0)
        & (paired["internal_length"] >= 0)
    )
    paired["passes_min_ltr_length"] = (
        (paired["lLTR length"] >= min_ltr_length) & (paired["rLTR length"] >= min_ltr_length)
    )
    paired["is_complete_yes"] = paired["Complete"].astype(str).str.lower().eq("yes")
    paired["is_ltr_order"] = paired["Order"].astype(str).eq("LTR")
    paired["has_5plus_domains"] = paired["domain_count"].fillna(0).astype(float) >= 5

    ltr = paired[paired["is_ltr_order"] & paired["has_valid_pair_coords"]].copy()

    gca_to_species = dict(zip(lookup["Genome_Accension"], lookup["Species"]))
    species_to_gca = {species: gca for gca, species in gca_to_species.items()}
    genome_files = find_local_genome_files(species_to_gca.values())

    readiness = (
        ltr.groupby("species")
        .agg(
            n_ltr_pairs=("species", "size"),
            n_complete_yes=("is_complete_yes", "sum"),
            n_with_min_ltr_length=("passes_min_ltr_length", "sum"),
            n_with_5plus_domains=("has_5plus_domains", "sum"),
            median_ltr5_length=("lLTR length", "median"),
            median_ltr3_length=("rLTR length", "median"),
            median_internal_length=("internal_length", "median"),
            gypsy_fraction=("Superfamily", lambda s: (s.astype(str) == "Gypsy").mean()),
        )
        .reset_index()
    )

    readiness["gca_accession"] = readiness["species"].map(species_to_gca)
    readiness["local_genome_fasta_available"] = readiness["gca_accession"].map(genome_files).notna()
    readiness["local_genome_fasta_path"] = readiness["gca_accession"].map(genome_files).fillna("")
    readiness["ready_for_sequence_age_estimation"] = (
        (readiness["n_with_min_ltr_length"] > 0) & readiness["local_genome_fasta_available"]
    )
    readiness = readiness.sort_values(
        ["ready_for_sequence_age_estimation", "n_ltr_pairs"],
        ascending=[False, False],
    ).reset_index(drop=True)
    return readiness


def build_overview(df: pd.DataFrame, readiness: pd.DataFrame, min_ltr_length: int) -> pd.DataFrame:
    paired = df.copy()
    paired["internal_length"] = paired["rLTR start"] - paired["lLTR end"] - 1
    paired["has_valid_pair_coords"] = (
        paired["lLTR start"].notna()
        & paired["rLTR start"].notna()
        & (paired["lLTR length"] > 0)
        & (paired["rLTR length"] > 0)
        & (paired["internal_length"] >= 0)
    )
    paired["passes_min_ltr_length"] = (
        (paired["lLTR length"] >= min_ltr_length) & (paired["rLTR length"] >= min_ltr_length)
    )
    paired["is_complete_yes"] = paired["Complete"].astype(str).str.lower().eq("yes")
    paired["is_ltr_order"] = paired["Order"].astype(str).eq("LTR")
    paired["has_5plus_domains"] = paired["domain_count"].fillna(0).astype(float) >= 5
    ltr = paired[paired["is_ltr_order"] & paired["has_valid_pair_coords"]].copy()

    overview_rows = [
        {"metric": "n_total_master_rows", "value": int(len(df))},
        {"metric": "n_ltr_valid_pairs", "value": int(len(ltr))},
        {"metric": "n_species_with_ltr_pairs", "value": int(ltr["species"].nunique())},
        {"metric": "n_complete_yes_ltr_pairs", "value": int(ltr["is_complete_yes"].sum())},
        {"metric": f"n_ltr_pairs_min_{min_ltr_length}bp", "value": int(ltr["passes_min_ltr_length"].sum())},
        {"metric": "n_ltr_pairs_5plus_domains", "value": int(ltr["has_5plus_domains"].sum())},
        {"metric": "n_ready_species_with_local_genomes", "value": int(readiness["ready_for_sequence_age_estimation"].sum())},
        {"metric": "dominant_superfamily", "value": ltr["Superfamily"].astype(str).value_counts().idxmax()},
        {"metric": "dominant_superfamily_count", "value": int(ltr["Superfamily"].astype(str).value_counts().iloc[0])},
    ]
    return pd.DataFrame(overview_rows)


def build_candidate_inventory(df: pd.DataFrame, min_ltr_length: int) -> pd.DataFrame:
    paired = df.copy()
    paired["internal_length"] = paired["rLTR start"] - paired["lLTR end"] - 1
    paired["has_valid_pair_coords"] = (
        paired["lLTR start"].notna()
        & paired["rLTR start"].notna()
        & (paired["lLTR length"] > 0)
        & (paired["rLTR length"] > 0)
        & (paired["internal_length"] >= 0)
    )
    paired["passes_min_ltr_length"] = (
        (paired["lLTR length"] >= min_ltr_length) & (paired["rLTR length"] >= min_ltr_length)
    )
    paired["is_complete_yes"] = paired["Complete"].astype(str).str.lower().eq("yes")
    paired["has_5plus_domains"] = paired["domain_count"].fillna(0).astype(float) >= 5
    paired = paired[paired["Order"].astype(str).eq("LTR") & paired["has_valid_pair_coords"]].copy()

    columns = [
        "species",
        "sequence",
        "element start",
        "element end",
        "lLTR start",
        "lLTR end",
        "rLTR start",
        "rLTR end",
        "lLTR length",
        "rLTR length",
        "internal_length",
        "Superfamily",
        "Complete",
        "domain_count",
        "passes_min_ltr_length",
        "has_5plus_domains",
        "is_complete_yes",
    ]
    return paired[columns].sort_values(["species", "sequence", "element start"]).reset_index(drop=True)


def render_audit_markdown(
    readiness: pd.DataFrame,
    overview: pd.DataFrame,
    min_ltr_length: int,
) -> str:
    metric_map = dict(zip(overview["metric"], overview["value"]))
    top_species = readiness.head(10)

    lines: List[str] = []
    lines.append("# LTR Age Audit")
    lines.append("")
    lines.append(
        "This note records the audit status of the local LTR insertion-age branch."
    )
    lines.append("")
    lines.append("## Bottom Line")
    lines.append("")
    lines.append(
        "- The previous local approach was not suitable for paper use because it tried to infer insertion age from RepeatMasker divergence-to-consensus hits rather than divergence between true paired 5' and 3' LTR sequences."
    )
    lines.append(
        "- The canonical local paired-element source is the ectopic master table produced by `scripts/processing/ec.py`, not the RepeatMasker adjacency heuristic."
    )
    lines.append(
        "- Local sequence-based age estimation is currently blocked because matching genome FASTA assemblies are not present in the repository workspace."
    )
    lines.append("")
    lines.append("## Current Local Readiness")
    lines.append("")
    lines.append(f"- Valid paired LTR elements in canonical ectopic master table: `{metric_map['n_ltr_valid_pairs']}`")
    lines.append(f"- Species with at least one paired LTR element: `{metric_map['n_species_with_ltr_pairs']}`")
    lines.append(f"- Complete (`Complete = yes`) LTR elements: `{metric_map['n_complete_yes_ltr_pairs']}`")
    lines.append(f"- LTR pairs with both LTRs >= `{min_ltr_length}` bp: `{metric_map[f'n_ltr_pairs_min_{min_ltr_length}bp']}`")
    lines.append(f"- LTR pairs with five or more annotated domains: `{metric_map['n_ltr_pairs_5plus_domains']}`")
    lines.append(
        f"- Dominant superfamily among paired LTR elements: `{metric_map['dominant_superfamily']}` (`{metric_map['dominant_superfamily_count']}` elements)"
    )
    lines.append(
        f"- Species currently ready for true sequence-based age estimation with local genome FASTA present: `{metric_map['n_ready_species_with_local_genomes']}`"
    )
    lines.append("")
    lines.append("## Highest-Coverage Species")
    lines.append("")
    for _, row in top_species.iterrows():
        lines.append(
            f"- `{row['species']}`: `{int(row['n_ltr_pairs'])}` paired LTR elements, "
            f"`{int(row['n_complete_yes'])}` complete, median internal length `{row['median_internal_length']:.1f}` bp"
        )
    lines.append("")
    lines.append("## Interpretation")
    lines.append("")
    lines.append(
        "- The paired-LTR substrate exists locally and is substantial, so the biological question is still tractable."
    )
    lines.append(
        "- What is missing is not paired-element annotation but the sequence-access layer required to compare the 5' and 3' LTRs directly."
    )
    lines.append(
        "- Until genome FASTA assemblies are available locally and wired into a sequence-extraction workflow, this branch should be treated as `blocked for paper-ready age inference`."
    )
    lines.append("")
    return "\n".join(lines) + "\n"


def main() -> int:
    args = parse_args()
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    lookup = load_lookup_table()
    master = load_master_table()

    readiness = build_species_readiness(master, lookup, args.min_ltr_length)
    overview = build_overview(master, readiness, args.min_ltr_length)
    inventory = build_candidate_inventory(master, args.min_ltr_length)

    readiness_path = OUTPUT_DIR / "ltr_age_readiness_by_species.csv"
    overview_path = OUTPUT_DIR / "ltr_age_readiness_overview.csv"
    inventory_path = OUTPUT_DIR / "ltr_age_candidate_inventory.csv"

    readiness.to_csv(readiness_path, index=False)
    overview.to_csv(overview_path, index=False)
    inventory.to_csv(inventory_path, index=False)

    AUDIT_PATH.write_text(
        render_audit_markdown(readiness, overview, args.min_ltr_length),
        encoding="utf-8",
    )

    ready_species = int(readiness["ready_for_sequence_age_estimation"].sum())
    print(f"Wrote species readiness to {readiness_path}")
    print(f"Wrote overview to {overview_path}")
    print(f"Wrote candidate inventory to {inventory_path}")
    print(f"Wrote audit note to {AUDIT_PATH}")
    print(f"Species ready for true local sequence-based age estimation: {ready_species}")

    if args.fail_on_unready and ready_species == 0:
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
