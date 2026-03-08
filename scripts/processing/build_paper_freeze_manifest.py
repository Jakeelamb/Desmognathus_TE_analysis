#!/usr/bin/env python3
"""Build a tracked paper-freeze manifest from the current repo state."""

from __future__ import annotations

import argparse
import csv
import hashlib
import subprocess
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Iterable, List, Tuple


PROJECT_ROOT = Path(__file__).resolve().parents[2]
OUTPUT_DIR = PROJECT_ROOT / "paper_freeze"
MANIFEST_PATH = PROJECT_ROOT / "PAPER_FREEZE_MANIFEST.md"

KEY_FILES = [
    "input_data/lookup_table.txt",
    "input_data/phylogeny/desmo900dated_test.tre",
    "results/data/dnaPipeTE_merged_classifications.csv",
    "results/data/dnaPipeTE_order_breakdown.csv",
    "results/data/dnaPipeTE_superfamily_breakdown.csv",
    "results/data/diversity_order_stats.csv",
    "results/data/diversity_superfamily_stats.csv",
    "results/data/ectopic_recombination_filtered_3000bp_5+domains_no_unknown_species.csv",
    "results/tables/pca/te_pca_analysis_manifest.csv",
    "results/tables/pca/te_ppca_analysis_manifest.csv",
    "path_analysis/data/derived/te_path_features.csv",
    "path_analysis/data/derived/te_model_feature_panel.csv",
    "path_analysis/data/derived/path_input_master.csv",
    "path_analysis/data/derived/analysis_panel_summary.csv",
    "path_analysis/data/derived/analysis_species_readiness.csv",
    "path_analysis/data/derived/phylofill_panel_comparison.csv",
]

CANONICAL_SCRIPTS = [
    "scripts/processing/diversity_stats.py",
    "scripts/processing/pca.R",
    "scripts/processing/phylogenetic_pca_analysis.R",
    "scripts/processing/parse_repeatmasker_landscape.py",
    "scripts/processing/pgls_analysis.R",
    "scripts/processing/permanova_analysis.R",
    "scripts/processing/trait_evolution.R",
    "scripts/processing/analyze_phylogenetic_signal.R",
    "scripts/processing/analyze_phylogenetic_correlogram.R",
    "scripts/visualization/plot_all_te_landscapes.R",
    "scripts/R/visualization/te_phylo_landscape.R",
    "path_analysis/scripts/build_analysis_panels.py",
    "path_analysis/scripts/path_model_scaffold.R",
]

PRIMARY_PANEL_FILES = [
    "te_genome_primary_mediumplus_model_ranking.csv",
    "te_genome_primary_strict_body_model_ranking.csv",
    "te_genome_ectopic_primary_mediumplus_model_ranking.csv",
    "te_genome_ectopic_primary_strict_body_model_ranking.csv",
    "te_genome_organismal_primary_mediumplus_model_ranking.csv",
    "te_genome_organismal_primary_strict_body_model_ranking.csv",
    "te_genome_ectopic_organismal_primary_mediumplus_model_ranking.csv",
    "te_genome_ectopic_organismal_primary_strict_body_model_ranking.csv",
]

SUPPLEMENTARY_NOTES = [
    "Compositional PCA and phylogenetic PCA under results/tables/pca and results/figures/pca.",
    "Phylogenetic signal, correlogram, and TE-landscape figures under results/tables/phylogenetic_signal and results/figures/phylo_signal.",
    "Landscape CSVs and landscape/phylo-landscape figures under results/landscapes, results/figures/landscape, and results/figures/phylo_landscape.",
]

PROVISIONAL_NOTES = [
    "Morphology-linked path families remain planning or sensitivity analyses only.",
    "Phylogenetic trait-imputation columns remain sensitivity-only and should not replace observed values.",
    "Order-level all-feature ordinations are supplementary only and should not carry the main biological interpretation.",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--source-commit",
        default=None,
        help="Commit hash for the frozen computational state. Defaults to HEAD.",
    )
    return parser.parse_args()


def run_git(args: List[str]) -> str:
    return subprocess.check_output(
        ["git"] + args,
        cwd=PROJECT_ROOT,
        text=True,
    ).strip()


def sha256_file(path: Path, chunk_size: int = 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while True:
            chunk = handle.read(chunk_size)
            if not chunk:
                break
            digest.update(chunk)
    return digest.hexdigest()


def format_size(num_bytes: int) -> str:
    units = ["B", "KB", "MB", "GB", "TB"]
    value = float(num_bytes)
    for unit in units:
        if value < 1024.0 or unit == units[-1]:
            if unit == "B":
                return f"{int(value)} {unit}"
            return f"{value:.1f} {unit}"
        value /= 1024.0
    return f"{num_bytes} B"


def inventory_group(relative_path: Path) -> str:
    parts = relative_path.parts
    if not parts:
        return ""
    if parts[0] != "results":
        return parts[0]
    if len(parts) == 2:
        return "/".join(parts[:2])
    if len(parts) >= 4 and parts[1] in {"data", "figures", "tables"}:
        return "/".join(parts[:3])
    if len(parts) >= 2:
        return "/".join(parts[:2])
    return "results"


def build_results_inventory() -> Tuple[List[Dict[str, object]], List[Dict[str, object]]]:
    rows: List[Dict[str, object]] = []
    summary: Dict[str, Dict[str, int]] = defaultdict(lambda: {"file_count": 0, "total_bytes": 0})
    for path in sorted((PROJECT_ROOT / "results").rglob("*")):
        if not path.is_file():
            continue
        rel_path = path.relative_to(PROJECT_ROOT)
        group = inventory_group(rel_path)
        size_bytes = path.stat().st_size
        rows.append(
            {
                "group": group,
                "relative_path": rel_path.as_posix(),
                "size_bytes": size_bytes,
                "sha256": sha256_file(path),
            }
        )
        summary[group]["file_count"] += 1
        summary[group]["total_bytes"] += size_bytes

    summary_rows = []
    for group in sorted(summary):
        summary_rows.append(
            {
                "group": group,
                "file_count": summary[group]["file_count"],
                "total_bytes": summary[group]["total_bytes"],
            }
        )
    return rows, summary_rows


def build_key_file_manifest() -> List[Dict[str, object]]:
    rows: List[Dict[str, object]] = []
    for relative in KEY_FILES:
        path = PROJECT_ROOT / relative
        exists = path.exists()
        size_bytes = path.stat().st_size if exists else ""
        sha_value = sha256_file(path) if exists and path.is_file() else ""
        rows.append(
            {
                "relative_path": relative,
                "exists": str(exists).lower(),
                "size_bytes": size_bytes,
                "sha256": sha_value,
            }
        )
    return rows


def read_csv_dicts(path: Path) -> List[Dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def read_panel_counts() -> Dict[str, int]:
    path = PROJECT_ROOT / "path_analysis/data/derived/analysis_panel_summary.csv"
    counts: Dict[str, int] = {}
    for row in read_csv_dicts(path):
        counts[row["panel_name"]] = int(row["n_species"])
    return counts


def read_model_winners() -> List[Dict[str, str]]:
    winners: List[Dict[str, str]] = []
    for filename in PRIMARY_PANEL_FILES:
        path = PROJECT_ROOT / "path_analysis/results" / filename
        rows = read_csv_dicts(path)
        if not rows:
            continue
        top_row = rows[0]
        winners.append(
            {
                "result_file": f"path_analysis/results/{filename}",
                "winner_model": top_row["model"],
                "delta_CICc_to_next": top_row["delta_CICc"],
                "weight": top_row["w"],
            }
        )
    return winners


def write_csv(path: Path, fieldnames: Iterable[str], rows: Iterable[Dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def render_markdown(
    source_commit: str,
    key_rows: List[Dict[str, object]],
    result_rows: List[Dict[str, object]],
    summary_rows: List[Dict[str, object]],
    panel_counts: Dict[str, int],
    winners: List[Dict[str, str]],
) -> str:
    now = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S UTC")
    total_files = len(result_rows)
    total_bytes = sum(int(row["size_bytes"]) for row in result_rows)
    key_present = [row for row in key_rows if row["exists"] == "true"]

    lines: List[str] = []
    lines.append("# Paper Freeze Manifest")
    lines.append("")
    lines.append(
        "This file records the current paper-facing freeze state for the Desmognathus TE project."
    )
    lines.append("")
    lines.append("## Frozen Computational State")
    lines.append("")
    lines.append(f"- Source code commit used for the frozen analysis state: `{source_commit}`")
    lines.append(f"- Manifest generated at: `{now}`")
    lines.append(
        "- Tracked inventory files: `paper_freeze/key_file_manifest.csv`, "
        "`paper_freeze/results_inventory.csv`, and `paper_freeze/results_summary.csv`"
    )
    lines.append(
        f"- Ignored generated outputs inventoried locally: `{total_files}` files totaling `{format_size(total_bytes)}`"
    )
    lines.append("")
    lines.append("## Canonical Entry Points")
    lines.append("")
    for script in CANONICAL_SCRIPTS:
        lines.append(f"- `{script}`")
    lines.append("")
    lines.append("## Analysis Tiers")
    lines.append("")
    lines.append("### Primary Paper Results")
    lines.append("")
    lines.append(
        "These are the most defensible current comparative outputs and should anchor the main text."
    )
    lines.append("")
    primary_panels = [
        ("te_genome_primary_mediumplus", "TE plus genome observed-only medium-plus"),
        ("te_genome_primary_strict_body", "TE plus genome adult-oriented strict-body"),
        ("te_genome_ectopic_primary_mediumplus", "TE plus genome plus ectopic observed-only medium-plus"),
        ("te_genome_ectopic_primary_strict_body", "TE plus genome plus ectopic adult-oriented strict-body"),
        ("te_genome_organismal_primary_mediumplus", "TE plus genome plus organismal observed-only medium-plus"),
        ("te_genome_organismal_primary_strict_body", "TE plus genome plus organismal adult-oriented strict-body"),
        ("te_genome_ectopic_organismal_primary_mediumplus", "TE plus genome plus ectopic plus organismal observed-only medium-plus"),
        ("te_genome_ectopic_organismal_primary_strict_body", "TE plus genome plus ectopic plus organismal adult-oriented strict-body"),
    ]
    for panel_name, label in primary_panels:
        if panel_name in panel_counts:
            lines.append(f"- `{panel_name}`: `{panel_counts[panel_name]}` species. {label}.")
    lines.append("")
    lines.append("Primary winner snapshots:")
    for winner in winners:
        result_file = winner["result_file"]
        model = winner["winner_model"]
        weight = winner["weight"]
        lines.append(f"- `{result_file}` -> `{model}` (weight `{weight}`)")
    lines.append("")
    lines.append("### Supplementary Results")
    lines.append("")
    for note in SUPPLEMENTARY_NOTES:
        lines.append(f"- {note}")
    lines.append("")
    lines.append("### Provisional or Sensitivity-Only Results")
    lines.append("")
    for note in PROVISIONAL_NOTES:
        lines.append(f"- {note}")
    lines.append("")
    lines.append("## Key Frozen Files")
    lines.append("")
    lines.append(
        "These are the most important frozen inputs and downstream tables to reference in methods, supplements, and reproduction notes."
    )
    lines.append("")
    for row in key_present:
        lines.append(
            f"- `{row['relative_path']}`  "
            f"size=`{row['size_bytes']}` bytes  "
            f"sha256=`{row['sha256']}`"
        )
    lines.append("")
    lines.append("## Ignored Results Inventory Summary")
    lines.append("")
    for row in summary_rows:
        lines.append(
            f"- `{row['group']}`: `{row['file_count']}` files, `{format_size(int(row['total_bytes']))}`"
        )
    lines.append("")
    lines.append("## Scope Notes")
    lines.append("")
    lines.append(
        "- `results/` remains git-ignored. The freeze state is therefore represented by the tracked checksum inventories rather than by versioning the generated files themselves."
    )
    lines.append(
        "- Observed trait values remain primary. `phylo_` columns remain explicit sensitivity-only inference."
    )
    lines.append(
        "- Morphology-linked path families remain excluded from the primary paper claims until independent final genome-size estimates are available."
    )
    lines.append("")
    return "\n".join(lines) + "\n"


def main() -> None:
    args = parse_args()
    source_commit = args.source_commit or run_git(["rev-parse", "HEAD"])

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    key_rows = build_key_file_manifest()
    result_rows, summary_rows = build_results_inventory()
    panel_counts = read_panel_counts()
    winners = read_model_winners()

    write_csv(
        OUTPUT_DIR / "key_file_manifest.csv",
        ["relative_path", "exists", "size_bytes", "sha256"],
        key_rows,
    )
    write_csv(
        OUTPUT_DIR / "results_inventory.csv",
        ["group", "relative_path", "size_bytes", "sha256"],
        result_rows,
    )
    write_csv(
        OUTPUT_DIR / "results_summary.csv",
        ["group", "file_count", "total_bytes"],
        summary_rows,
    )
    (OUTPUT_DIR / "source_commit.txt").write_text(source_commit + "\n", encoding="utf-8")
    MANIFEST_PATH.write_text(
        render_markdown(source_commit, key_rows, result_rows, summary_rows, panel_counts, winners),
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
