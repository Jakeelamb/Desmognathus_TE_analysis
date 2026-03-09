#!/usr/bin/env python3
"""
LTR age audit and sequence-divergence estimation.

This workflow treats the ectopic master table produced by `scripts/processing/ec.py`
as the canonical paired-LTR annotation source. It first audits paired-LTR readiness
and local genome availability, then optionally estimates true 5'/3' LTR sequence
divergence from the corresponding genome assemblies.

Absolute insertion ages are only reported when an explicit substitution rate is
supplied. By default, the validated output is sequence-based divergence rather than
calibrated age in years.
"""

from __future__ import annotations

import argparse
import gzip
import math
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Optional, Sequence, Tuple

import pandas as pd

import sys

sys.path.insert(0, str(Path(__file__).parent.parent))
from config import PROJECT_ROOT, load_lookup_table, paths  # noqa: E402


MASTER_TABLE = paths.results.data / "ectopic_recombination_master.csv"
OUTPUT_DIR = paths.results.data / "ltr_age"
AUDIT_PATH = PROJECT_ROOT / "LTR_AGE_AUDIT.md"
TRANSITIONS = {
    ("A", "G"),
    ("G", "A"),
    ("C", "T"),
    ("T", "C"),
}
DNA_BASES = {"A", "C", "G", "T"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Audit readiness and estimate sequence-based LTR insertion divergence."
    )
    parser.add_argument(
        "--mode",
        choices=["auto", "audit-only", "estimate"],
        default="auto",
        help=(
            "Run mode: `auto` estimates divergence when local genomes are available, "
            "`audit-only` skips sequence extraction, and `estimate` requires sequence estimation."
        ),
    )
    parser.add_argument(
        "--min-ltr-length",
        type=int,
        default=100,
        help="Minimum LTR length used for readiness counts and estimation candidates (default: 100 bp).",
    )
    parser.add_argument(
        "--min-comparable-sites",
        type=int,
        default=50,
        help="Minimum ungapped A/C/G/T sites required for a valid divergence estimate (default: 50).",
    )
    parser.add_argument(
        "--substitution-rate",
        type=float,
        default=None,
        help=(
            "Optional neutral substitution rate in substitutions/site/year. "
            "If supplied, age_years and age_mya are computed from K2P distance."
        ),
    )
    parser.add_argument(
        "--fail-on-unready",
        action="store_true",
        help="Exit non-zero if no species are locally ready for sequence-based estimation.",
    )
    return parser.parse_args()


def load_master_table() -> pd.DataFrame:
    if not MASTER_TABLE.exists():
        raise FileNotFoundError(f"Missing canonical ectopic master table: {MASTER_TABLE}")

    df = pd.read_csv(MASTER_TABLE, sep="\t")
    required = {
        "element start",
        "element end",
        "sequence",
        "lLTR start",
        "lLTR end",
        "lLTR length",
        "rLTR start",
        "rLTR end",
        "rLTR length",
        "species",
        "Order",
        "Superfamily",
        "Complete",
        "Strand",
        "domain_count",
    }
    missing = sorted(required - set(df.columns))
    if missing:
        raise ValueError(f"Ectopic master table is missing required columns: {missing}")
    return df


def find_local_genome_files(gca_values: Iterable[str]) -> Dict[str, str]:
    suffixes = (".fna", ".fa", ".fasta", ".fna.gz", ".fa.gz", ".fasta.gz")
    search_roots = [
        getattr(paths.input_data, "genomes", PROJECT_ROOT / "input_data" / "genomes"),
        PROJECT_ROOT / "input_data",
        PROJECT_ROOT,
    ]
    found: Dict[str, str] = {}
    for gca in gca_values:
        for root in search_roots:
            root = Path(root)
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


def build_ltr_inventory(df: pd.DataFrame, lookup: pd.DataFrame, min_ltr_length: int) -> pd.DataFrame:
    paired = df.copy()
    numeric_columns = [
        "element start",
        "element end",
        "lLTR start",
        "lLTR end",
        "lLTR length",
        "rLTR start",
        "rLTR end",
        "rLTR length",
        "domain_count",
    ]
    for column in numeric_columns:
        paired[column] = pd.to_numeric(paired[column], errors="coerce")

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
    paired["sequence_root"] = paired["sequence"].astype(str).str.replace(r"_De$", "", regex=True)
    paired["recommended_high_confidence"] = paired["is_complete_yes"] & paired["has_5plus_domains"]

    species_to_gca = dict(zip(lookup["Species"], lookup["Genome_Accension"]))
    genome_files = find_local_genome_files(species_to_gca.values())

    inventory = paired[paired["is_ltr_order"] & paired["has_valid_pair_coords"]].copy()
    inventory["gca_accession"] = inventory["species"].map(species_to_gca).fillna("")
    inventory["local_genome_fasta_path"] = inventory["gca_accession"].map(genome_files).fillna("")
    inventory["local_genome_fasta_available"] = inventory["local_genome_fasta_path"].astype(bool)
    inventory["ready_for_sequence_age_estimation"] = (
        inventory["passes_min_ltr_length"] & inventory["local_genome_fasta_available"]
    )
    return inventory.reset_index(drop=True)


def build_species_readiness(inventory: pd.DataFrame, min_ltr_length: int) -> pd.DataFrame:
    readiness = (
        inventory.groupby("species")
        .agg(
            n_ltr_pairs=("species", "size"),
            n_complete_yes=("is_complete_yes", "sum"),
            n_with_min_ltr_length=("passes_min_ltr_length", "sum"),
            n_with_5plus_domains=("has_5plus_domains", "sum"),
            n_recommended_high_confidence=("recommended_high_confidence", "sum"),
            median_ltr5_length=("lLTR length", "median"),
            median_ltr3_length=("rLTR length", "median"),
            median_internal_length=("internal_length", "median"),
            gypsy_fraction=("Superfamily", lambda s: (s.astype(str) == "Gypsy").mean()),
            gca_accession=("gca_accession", "first"),
            local_genome_fasta_available=("local_genome_fasta_available", "max"),
            local_genome_fasta_path=("local_genome_fasta_path", "first"),
        )
        .reset_index()
    )
    readiness["ready_for_sequence_age_estimation"] = (
        (readiness["n_with_min_ltr_length"] > 0) & readiness["local_genome_fasta_available"]
    )
    readiness = readiness.sort_values(
        ["ready_for_sequence_age_estimation", "n_ltr_pairs"],
        ascending=[False, False],
    ).reset_index(drop=True)
    readiness["min_ltr_length_bp"] = min_ltr_length
    return readiness


def build_overview(
    master: pd.DataFrame,
    inventory: pd.DataFrame,
    readiness: pd.DataFrame,
    min_ltr_length: int,
    pairwise: Optional[pd.DataFrame] = None,
) -> pd.DataFrame:
    overview_rows = [
        {"metric": "n_total_master_rows", "value": int(len(master))},
        {"metric": "n_ltr_valid_pairs", "value": int(len(inventory))},
        {"metric": "n_species_with_ltr_pairs", "value": int(inventory["species"].nunique())},
        {
            "metric": "n_pairs_without_lookup_accession",
            "value": int(inventory["gca_accession"].astype(str).eq("").sum()),
        },
        {
            "metric": "n_species_without_lookup_accession",
            "value": int(inventory.loc[inventory["gca_accession"].astype(str).eq(""), "species"].nunique()),
        },
        {"metric": "n_complete_yes_ltr_pairs", "value": int(inventory["is_complete_yes"].sum())},
        {"metric": f"n_ltr_pairs_min_{min_ltr_length}bp", "value": int(inventory["passes_min_ltr_length"].sum())},
        {"metric": "n_ltr_pairs_5plus_domains", "value": int(inventory["has_5plus_domains"].sum())},
        {
            "metric": "n_ltr_pairs_recommended_high_confidence",
            "value": int(inventory["recommended_high_confidence"].sum()),
        },
        {
            "metric": "n_ready_species_with_local_genomes",
            "value": int(readiness["ready_for_sequence_age_estimation"].sum()),
        },
        {
            "metric": "dominant_superfamily",
            "value": inventory["Superfamily"].astype(str).value_counts().idxmax(),
        },
        {
            "metric": "dominant_superfamily_count",
            "value": int(inventory["Superfamily"].astype(str).value_counts().iloc[0]),
        },
    ]

    if pairwise is not None and not pairwise.empty:
        success = pairwise[pairwise["estimation_status"].eq("estimated")].copy()
        overview_rows.extend(
            [
                {"metric": "n_pairs_attempted_sequence_estimation", "value": int(len(pairwise))},
                {"metric": "n_pairs_successful_sequence_estimation", "value": int(len(success))},
                {"metric": "n_species_with_sequence_estimates", "value": int(success["species"].nunique())},
                {
                    "metric": "n_high_confidence_pairs_successful_sequence_estimation",
                    "value": int(success["recommended_high_confidence"].sum()),
                },
            ]
        )
        if not success.empty:
            overview_rows.extend(
                [
                    {
                        "metric": "median_sequence_p_distance",
                        "value": round(float(success["p_distance"].median()), 6),
                    },
                    {
                        "metric": "median_sequence_k2p_distance",
                        "value": round(float(success["k2p_distance"].median()), 6),
                    },
                    {
                        "metric": "median_sequence_comparable_sites",
                        "value": int(success["comparable_sites"].median()),
                    },
                ]
            )
            if success["age_mya"].notna().any():
                overview_rows.append(
                    {
                        "metric": "median_sequence_age_mya",
                        "value": round(float(success["age_mya"].dropna().median()), 6),
                    }
                )

    return pd.DataFrame(overview_rows)


def iter_fasta_records(path: Path) -> Iterator[Tuple[str, str]]:
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as handle:
        header: Optional[str] = None
        seq_parts: List[str] = []
        for line in handle:
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_parts)
                header = line[1:].strip().split()[0]
                seq_parts = []
            else:
                seq_parts.append(line.strip())
        if header is not None:
            yield header, "".join(seq_parts)


def extract_target_sequences(genome_path: Path, target_ids: Sequence[str]) -> Dict[str, str]:
    remaining = set(target_ids)
    extracted: Dict[str, str] = {}
    if not remaining:
        return extracted
    for header, sequence in iter_fasta_records(genome_path):
        if header in remaining:
            extracted[header] = sequence.upper()
            remaining.remove(header)
        if not remaining:
            break
    return extracted


def load_pairwise_aligner():
    try:
        from Bio.Align import PairwiseAligner
    except ImportError as exc:  # pragma: no cover - environment dependent
        raise ImportError(
            "Sequence-based LTR divergence estimation requires Biopython. "
            "Run this script in the `Dusky` environment or install `biopython`."
        ) from exc

    aligner = PairwiseAligner(mode="global")
    aligner.match_score = 2.0
    aligner.mismatch_score = -3.0
    aligner.open_gap_score = -5.0
    aligner.extend_gap_score = -1.0
    return aligner


def reconstruct_alignment(seq_left: str, seq_right: str, alignment) -> Tuple[str, str]:
    left_parts: List[str] = []
    right_parts: List[str] = []
    coords = alignment.coordinates
    for idx in range(coords.shape[1] - 1):
        left_start = int(coords[0, idx])
        left_end = int(coords[0, idx + 1])
        right_start = int(coords[1, idx])
        right_end = int(coords[1, idx + 1])

        left_step = left_end - left_start
        right_step = right_end - right_start

        if left_step and right_step:
            left_parts.append(seq_left[left_start:left_end])
            right_parts.append(seq_right[right_start:right_end])
        elif left_step:
            left_parts.append(seq_left[left_start:left_end])
            right_parts.append("-" * left_step)
        elif right_step:
            left_parts.append("-" * right_step)
            right_parts.append(seq_right[right_start:right_end])

    return "".join(left_parts), "".join(right_parts)


def compute_alignment_metrics(
    aligned_left: str,
    aligned_right: str,
    min_comparable_sites: int,
    substitution_rate: Optional[float],
) -> Dict[str, float]:
    matches = 0
    mismatches = 0
    gap_columns = 0
    ambiguous_columns = 0
    transitions = 0
    transversions = 0

    for left_base, right_base in zip(aligned_left, aligned_right):
        if left_base == "-" or right_base == "-":
            gap_columns += 1
            continue

        left_base = left_base.upper()
        right_base = right_base.upper()
        if left_base not in DNA_BASES or right_base not in DNA_BASES:
            ambiguous_columns += 1
            continue

        if left_base == right_base:
            matches += 1
        else:
            mismatches += 1
            if (left_base, right_base) in TRANSITIONS:
                transitions += 1
            else:
                transversions += 1

    comparable_sites = matches + mismatches
    p_distance = math.nan
    k2p_distance = math.nan
    age_years = math.nan
    age_mya = math.nan
    status = "estimated"

    if comparable_sites < min_comparable_sites:
        status = "insufficient_comparable_sites"
    elif comparable_sites > 0:
        p_distance = mismatches / comparable_sites
        transition_fraction = transitions / comparable_sites
        transversion_fraction = transversions / comparable_sites
        term_1 = 1 - (2 * transition_fraction) - transversion_fraction
        term_2 = 1 - (2 * transversion_fraction)
        if term_1 > 0 and term_2 > 0:
            k2p_distance = (-0.5 * math.log(term_1)) - (0.25 * math.log(term_2))
            if substitution_rate:
                age_years = k2p_distance / (2 * substitution_rate)
                age_mya = age_years / 1_000_000
        else:
            status = "k2p_out_of_domain"

    return {
        "alignment_length": len(aligned_left),
        "comparable_sites": comparable_sites,
        "matches": matches,
        "mismatches": mismatches,
        "gap_columns": gap_columns,
        "ambiguous_columns": ambiguous_columns,
        "transitions": transitions,
        "transversions": transversions,
        "p_distance": p_distance,
        "k2p_distance": k2p_distance,
        "age_years": age_years,
        "age_mya": age_mya,
        "estimation_status": status,
    }


def estimate_sequence_divergence(
    inventory: pd.DataFrame,
    min_comparable_sites: int,
    substitution_rate: Optional[float],
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    aligner = load_pairwise_aligner()
    result_rows: List[Dict[str, object]] = []
    scaffold_rows: List[Dict[str, object]] = []

    for species, species_rows in inventory.groupby("species", sort=True):
        species_rows = species_rows.copy()
        genome_path_str = str(species_rows["local_genome_fasta_path"].iloc[0]).strip()
        gca_accession = str(species_rows["gca_accession"].iloc[0]).strip()
        target_scaffolds = sorted(set(species_rows["sequence_root"].astype(str)))

        if not genome_path_str:
            scaffold_rows.append(
                {
                    "species": species,
                    "gca_accession": gca_accession,
                    "local_genome_fasta_path": "",
                    "requested_scaffolds": len(target_scaffolds),
                    "extracted_scaffolds": 0,
                    "missing_scaffolds": len(target_scaffolds),
                    "status": "genome_missing",
                }
            )
            for _, row in species_rows.iterrows():
                record = row.to_dict()
                record.update(
                    {
                        "extracted_ltr5_length": math.nan,
                        "extracted_ltr3_length": math.nan,
                        "alignment_length": math.nan,
                        "comparable_sites": math.nan,
                        "matches": math.nan,
                        "mismatches": math.nan,
                        "gap_columns": math.nan,
                        "ambiguous_columns": math.nan,
                        "transitions": math.nan,
                        "transversions": math.nan,
                        "p_distance": math.nan,
                        "k2p_distance": math.nan,
                        "age_years": math.nan,
                        "age_mya": math.nan,
                        "estimation_status": "genome_missing",
                    }
                )
                result_rows.append(record)
            continue

        genome_path = PROJECT_ROOT / genome_path_str
        extracted = extract_target_sequences(genome_path, target_scaffolds)
        scaffold_rows.append(
            {
                "species": species,
                "gca_accession": gca_accession,
                "local_genome_fasta_path": genome_path_str,
                "requested_scaffolds": len(target_scaffolds),
                "extracted_scaffolds": len(extracted),
                "missing_scaffolds": len(target_scaffolds) - len(extracted),
                "status": "complete" if len(extracted) == len(target_scaffolds) else "partial",
            }
        )

        for _, row in species_rows.iterrows():
            record = row.to_dict()
            scaffold_id = str(row["sequence_root"])
            contig = extracted.get(scaffold_id)
            if contig is None:
                record.update(
                    {
                        "extracted_ltr5_length": math.nan,
                        "extracted_ltr3_length": math.nan,
                        "alignment_length": math.nan,
                        "comparable_sites": math.nan,
                        "matches": math.nan,
                        "mismatches": math.nan,
                        "gap_columns": math.nan,
                        "ambiguous_columns": math.nan,
                        "transitions": math.nan,
                        "transversions": math.nan,
                        "p_distance": math.nan,
                        "k2p_distance": math.nan,
                        "age_years": math.nan,
                        "age_mya": math.nan,
                        "estimation_status": "scaffold_missing",
                    }
                )
                result_rows.append(record)
                continue

            left_start = int(row["lLTR start"]) - 1
            left_end = int(row["lLTR end"])
            right_start = int(row["rLTR start"]) - 1
            right_end = int(row["rLTR end"])

            if left_start < 0 or right_start < 0 or left_end > len(contig) or right_end > len(contig):
                record.update(
                    {
                        "extracted_ltr5_length": math.nan,
                        "extracted_ltr3_length": math.nan,
                        "alignment_length": math.nan,
                        "comparable_sites": math.nan,
                        "matches": math.nan,
                        "mismatches": math.nan,
                        "gap_columns": math.nan,
                        "ambiguous_columns": math.nan,
                        "transitions": math.nan,
                        "transversions": math.nan,
                        "p_distance": math.nan,
                        "k2p_distance": math.nan,
                        "age_years": math.nan,
                        "age_mya": math.nan,
                        "estimation_status": "coordinates_out_of_bounds",
                    }
                )
                result_rows.append(record)
                continue

            left_ltr = contig[left_start:left_end]
            right_ltr = contig[right_start:right_end]
            if not left_ltr or not right_ltr:
                record.update(
                    {
                        "extracted_ltr5_length": len(left_ltr),
                        "extracted_ltr3_length": len(right_ltr),
                        "alignment_length": math.nan,
                        "comparable_sites": math.nan,
                        "matches": math.nan,
                        "mismatches": math.nan,
                        "gap_columns": math.nan,
                        "ambiguous_columns": math.nan,
                        "transitions": math.nan,
                        "transversions": math.nan,
                        "p_distance": math.nan,
                        "k2p_distance": math.nan,
                        "age_years": math.nan,
                        "age_mya": math.nan,
                        "estimation_status": "empty_extraction",
                    }
                )
                result_rows.append(record)
                continue

            alignment = aligner.align(left_ltr, right_ltr)[0]
            aligned_left, aligned_right = reconstruct_alignment(left_ltr, right_ltr, alignment)
            metrics = compute_alignment_metrics(
                aligned_left,
                aligned_right,
                min_comparable_sites=min_comparable_sites,
                substitution_rate=substitution_rate,
            )
            record.update(
                {
                    "extracted_ltr5_length": len(left_ltr),
                    "extracted_ltr3_length": len(right_ltr),
                    **metrics,
                }
            )
            result_rows.append(record)

    pairwise = pd.DataFrame(result_rows)
    scaffold_summary = pd.DataFrame(scaffold_rows)

    success = pairwise[pairwise["estimation_status"].eq("estimated")].copy()
    if success.empty:
        species_summary = pd.DataFrame(
            columns=[
                "species",
                "gca_accession",
                "local_genome_fasta_path",
                "n_pairs_attempted",
                "n_pairs_estimated",
                "n_high_confidence_pairs_estimated",
                "median_p_distance",
                "median_k2p_distance",
                "median_comparable_sites",
                "median_age_mya",
            ]
        )
    else:
        species_summary = (
            success.groupby(["species", "gca_accession", "local_genome_fasta_path"])
            .agg(
                n_pairs_attempted=("species", "size"),
                n_pairs_estimated=("species", "size"),
                n_high_confidence_pairs_estimated=("recommended_high_confidence", "sum"),
                median_p_distance=("p_distance", "median"),
                median_k2p_distance=("k2p_distance", "median"),
                median_comparable_sites=("comparable_sites", "median"),
                median_age_mya=("age_mya", "median"),
            )
            .reset_index()
            .sort_values("n_pairs_estimated", ascending=False)
        )

    return pairwise, species_summary, scaffold_summary


def render_audit_markdown(
    readiness: pd.DataFrame,
    overview: pd.DataFrame,
    min_ltr_length: int,
    pairwise: Optional[pd.DataFrame],
    species_summary: Optional[pd.DataFrame],
    substitution_rate: Optional[float],
) -> str:
    metric_map = dict(zip(overview["metric"], overview["value"]))
    top_species = readiness.head(10)
    successful_pairwise = None
    if pairwise is not None and not pairwise.empty:
        successful_pairwise = pairwise[pairwise["estimation_status"].eq("estimated")].copy()

    lines: List[str] = []
    lines.append("# LTR Age Audit")
    lines.append("")
    lines.append("This note records the current status of the local LTR insertion-age branch.")
    lines.append("")
    lines.append("## Bottom Line")
    lines.append("")
    lines.append(
        "- The previous local approach was not suitable for paper use because it tried to infer insertion age from RepeatMasker divergence-to-consensus hits rather than divergence between true paired 5' and 3' LTR sequences."
    )
    lines.append(
        "- The canonical local paired-element source is the ectopic master table produced by `scripts/processing/ec.py`, not the RepeatMasker adjacency heuristic."
    )
    if successful_pairwise is None or successful_pairwise.empty:
        lines.append(
            "- Local sequence-based age estimation is still blocked because matching genome FASTA assemblies are not present or do not yield comparable paired-LTR sequences."
        )
    else:
        lines.append(
            "- Local genome FASTA assemblies are now available, so this branch computes true 5' and 3' LTR sequence divergence directly from assembly coordinates."
        )
        if substitution_rate is None:
            lines.append(
                "- The default local output is sequence divergence, not absolute age in years, because no substitution rate is imposed automatically."
            )
        else:
            lines.append(
                f"- Absolute age estimates are reported using the supplied neutral substitution rate `{substitution_rate}` substitutions/site/year."
            )
    lines.append("")
    lines.append("## Current Local Readiness")
    lines.append("")
    lines.append(f"- Valid paired LTR elements in canonical ectopic master table: `{metric_map['n_ltr_valid_pairs']}`")
    lines.append(f"- Species with at least one paired LTR element: `{metric_map['n_species_with_ltr_pairs']}`")
    lines.append(
        f"- Pairs lacking a lookup-table genome accession: `{metric_map['n_pairs_without_lookup_accession']}` "
        f"across `{metric_map['n_species_without_lookup_accession']}` species bucket(s)"
    )
    lines.append(f"- Complete (`Complete = yes`) LTR elements: `{metric_map['n_complete_yes_ltr_pairs']}`")
    lines.append(f"- LTR pairs with both LTRs >= `{min_ltr_length}` bp: `{metric_map[f'n_ltr_pairs_min_{min_ltr_length}bp']}`")
    lines.append(f"- LTR pairs with five or more annotated domains: `{metric_map['n_ltr_pairs_5plus_domains']}`")
    lines.append(f"- Recommended high-confidence LTR pairs (`Complete = yes` and `5+` domains): `{metric_map['n_ltr_pairs_recommended_high_confidence']}`")
    lines.append(
        f"- Dominant superfamily among paired LTR elements: `{metric_map['dominant_superfamily']}` (`{metric_map['dominant_superfamily_count']}` elements)"
    )
    lines.append(
        f"- Species currently ready for local sequence extraction with genome FASTA present: `{metric_map['n_ready_species_with_local_genomes']}`"
    )
    lines.append("")

    if successful_pairwise is not None and not successful_pairwise.empty:
        lines.append("## Sequence-Based Divergence")
        lines.append("")
        lines.append(
            f"- LTR pairs attempted for sequence estimation: `{metric_map['n_pairs_attempted_sequence_estimation']}`"
        )
        lines.append(
            f"- LTR pairs with successful sequence-based divergence estimates: `{metric_map['n_pairs_successful_sequence_estimation']}`"
        )
        lines.append(
            f"- Species with at least one successful sequence estimate: `{metric_map['n_species_with_sequence_estimates']}`"
        )
        lines.append(
            f"- High-confidence pairs with successful sequence estimates: `{metric_map['n_high_confidence_pairs_successful_sequence_estimation']}`"
        )
        lines.append(
            f"- Median ungapped comparable sites per successful alignment: `{metric_map['median_sequence_comparable_sites']}`"
        )
        lines.append(f"- Median p-distance across successful alignments: `{metric_map['median_sequence_p_distance']}`")
        lines.append(f"- Median K2P distance across successful alignments: `{metric_map['median_sequence_k2p_distance']}`")
        if substitution_rate is not None and "median_sequence_age_mya" in metric_map:
            lines.append(
                f"- Median inferred insertion age across successful alignments: `{metric_map['median_sequence_age_mya']}` Myr"
            )
        lines.append("")

    lines.append("## Highest-Coverage Species")
    lines.append("")
    if species_summary is not None and not species_summary.empty:
        for _, row in species_summary.head(10).iterrows():
            lines.append(
                f"- `{row['species']}`: `{int(row['n_pairs_estimated'])}` successful sequence estimates, "
                f"median K2P `{row['median_k2p_distance']:.4f}`, median comparable sites `{int(row['median_comparable_sites'])}`"
            )
    else:
        for _, row in top_species.iterrows():
            lines.append(
                f"- `{row['species']}`: `{int(row['n_ltr_pairs'])}` paired LTR elements, "
                f"`{int(row['n_complete_yes'])}` complete, median internal length `{row['median_internal_length']:.1f}` bp"
            )
    lines.append("")
    lines.append("## Interpretation")
    lines.append("")
    lines.append(
        "- The paired-LTR substrate exists locally and is substantial, so the biological question is tractable without falling back to divergence-to-consensus heuristics."
    )
    if successful_pairwise is None or successful_pairwise.empty:
        lines.append(
            "- The remaining blocker is sequence access or extraction, not paired-element annotation."
        )
    else:
        lines.append(
            "- The current local branch now supports direct 5'/3' LTR divergence estimation from assembly coordinates, which is the defensible basis for LTR insertion-age inference."
        )
        lines.append(
            "- Absolute age claims should remain conditional on an externally justified substitution rate; the default repo output stays at the divergence level."
        )
    lines.append("")
    return "\n".join(lines) + "\n"


def main() -> int:
    args = parse_args()
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    lookup = load_lookup_table()
    master = load_master_table()
    inventory = build_ltr_inventory(master, lookup, args.min_ltr_length)
    readiness = build_species_readiness(inventory, args.min_ltr_length)

    pairwise: Optional[pd.DataFrame] = None
    species_summary: Optional[pd.DataFrame] = None
    scaffold_summary: Optional[pd.DataFrame] = None

    ready_species = int(readiness["ready_for_sequence_age_estimation"].sum())
    should_estimate = args.mode == "estimate" or (args.mode == "auto" and ready_species > 0)

    if should_estimate:
        pairwise, species_summary, scaffold_summary = estimate_sequence_divergence(
            inventory[inventory["passes_min_ltr_length"]].copy(),
            min_comparable_sites=args.min_comparable_sites,
            substitution_rate=args.substitution_rate,
        )
    elif args.mode == "estimate":
        raise RuntimeError("Sequence-estimation mode requested, but no locally ready species were found.")

    overview = build_overview(
        master=master,
        inventory=inventory,
        readiness=readiness,
        min_ltr_length=args.min_ltr_length,
        pairwise=pairwise,
    )

    readiness_path = OUTPUT_DIR / "ltr_age_readiness_by_species.csv"
    overview_path = OUTPUT_DIR / "ltr_age_readiness_overview.csv"
    inventory_path = OUTPUT_DIR / "ltr_age_candidate_inventory.csv"
    readiness.to_csv(readiness_path, index=False)
    overview.to_csv(overview_path, index=False)
    inventory.to_csv(inventory_path, index=False)

    if pairwise is not None:
        pairwise_path = OUTPUT_DIR / "ltr_age_pairwise_divergence.csv"
        pairwise.to_csv(pairwise_path, index=False)
        print(f"Wrote pairwise divergence to {pairwise_path}")
    if species_summary is not None:
        species_summary_path = OUTPUT_DIR / "ltr_age_species_summary.csv"
        species_summary.to_csv(species_summary_path, index=False)
        print(f"Wrote species summary to {species_summary_path}")
    if scaffold_summary is not None:
        scaffold_summary_path = OUTPUT_DIR / "ltr_age_scaffold_extraction_summary.csv"
        scaffold_summary.to_csv(scaffold_summary_path, index=False)
        print(f"Wrote scaffold extraction summary to {scaffold_summary_path}")

    AUDIT_PATH.write_text(
        render_audit_markdown(
            readiness=readiness,
            overview=overview,
            min_ltr_length=args.min_ltr_length,
            pairwise=pairwise,
            species_summary=species_summary,
            substitution_rate=args.substitution_rate,
        ),
        encoding="utf-8",
    )

    print(f"Wrote species readiness to {readiness_path}")
    print(f"Wrote overview to {overview_path}")
    print(f"Wrote candidate inventory to {inventory_path}")
    print(f"Wrote audit note to {AUDIT_PATH}")
    print(f"Species ready for local sequence extraction: {ready_species}")

    if args.fail_on_unready and ready_species == 0:
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
