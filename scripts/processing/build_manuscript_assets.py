#!/usr/bin/env python3
"""Build manuscript-facing summary assets from the frozen primary results."""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Dict, List


PROJECT_ROOT = Path(__file__).resolve().parents[2]
OUTPUT_DIR = PROJECT_ROOT / "paper_freeze"

PRIMARY_RESULT_SPECS = [
    {
        "family_id": "te_genome",
        "panel_id": "te_genome_primary_mediumplus",
        "subset": "mediumplus",
        "tier": "primary",
        "result_role": "main_text",
        "core_claim": "Baseline TE model favors mediated_evenness and supports a stable LTR-balance to TE-evenness relationship.",
    },
    {
        "family_id": "te_genome",
        "panel_id": "te_genome_primary_strict_body",
        "subset": "strict_body",
        "tier": "robustness",
        "result_role": "supplement",
        "core_claim": "Adult-oriented strict-body subset retains the same winner, supporting baseline TE-model robustness.",
    },
    {
        "family_id": "te_genome_ectopic",
        "panel_id": "te_genome_ectopic_primary_mediumplus",
        "subset": "mediumplus",
        "tier": "primary",
        "result_role": "main_text",
        "core_claim": "Adding ectopic support favors ectopic_only before organismal competition is introduced.",
    },
    {
        "family_id": "te_genome_ectopic",
        "panel_id": "te_genome_ectopic_primary_strict_body",
        "subset": "strict_body",
        "tier": "robustness",
        "result_role": "supplement",
        "core_claim": "The ectopic-only result remains the best-supported ectopic model in the strict-body subset.",
    },
    {
        "family_id": "te_genome_organismal",
        "panel_id": "te_genome_organismal_primary_mediumplus",
        "subset": "mediumplus",
        "tier": "primary",
        "result_role": "main_text",
        "core_claim": "Body size is the only organismal covariate that meaningfully improves the TE-genome family.",
    },
    {
        "family_id": "te_genome_organismal",
        "panel_id": "te_genome_organismal_primary_strict_body",
        "subset": "strict_body",
        "tier": "robustness",
        "result_role": "supplement",
        "core_claim": "The body-size signal strengthens in the strict-body subset.",
    },
    {
        "family_id": "te_genome_ectopic_organismal",
        "panel_id": "te_genome_ectopic_organismal_primary_mediumplus",
        "subset": "mediumplus",
        "tier": "primary",
        "result_role": "main_text",
        "core_claim": "When ectopic support and body size compete directly, the best model drops ectopic_index and keeps TE plus body size.",
    },
    {
        "family_id": "te_genome_ectopic_organismal",
        "panel_id": "te_genome_ectopic_organismal_primary_strict_body",
        "subset": "strict_body",
        "tier": "robustness",
        "result_role": "supplement",
        "core_claim": "The same combined-family winner persists in the strict-body subset, with smaller support margins.",
    },
]

MAIN_TEXT_FIGURES = [
    {
        "figure_id": "Figure 1",
        "status": "main_text",
        "asset_path": "path_analysis/results/te_genome_primary_mediumplus_best_model.pdf",
        "purpose": "Baseline path model for the observed-only TE plus genome family.",
        "reason": "Shows the cleanest mechanistic TE result before adding ectopic or organismal competition.",
    },
    {
        "figure_id": "Figure 2",
        "status": "main_text",
        "asset_path": "path_analysis/results/te_genome_ectopic_organismal_primary_mediumplus_best_model.pdf",
        "purpose": "Combined-family best model when ectopic support and body size are allowed to compete.",
        "reason": "Captures the most important updated inference: body size remains while ectopic_index drops.",
    },
    {
        "figure_id": "Figure 3",
        "status": "main_text",
        "asset_path": "results/figures/pca/superfamily_primary_scores_pc1_pc2.png",
        "purpose": "Primary TE ordination at the superfamily level.",
        "reason": "Uses the canonical compositional PCA and provides one visually interpretable summary of TE composition structure.",
    },
]

SUPPLEMENT_FIGURES = [
    {
        "figure_id": "Figure S1",
        "status": "supplement",
        "asset_path": "path_analysis/results/te_genome_primary_strict_body_best_model.pdf",
        "purpose": "Strict-body robustness check for the baseline TE family.",
        "reason": "Confirms winner stability under adult-oriented body-size restriction.",
    },
    {
        "figure_id": "Figure S2",
        "status": "supplement",
        "asset_path": "path_analysis/results/te_genome_ectopic_primary_mediumplus_best_model.pdf",
        "purpose": "Best ectopic-only family before organismal competition.",
        "reason": "Useful to show the intermediate result that ectopic support is strongest only before body-size competition.",
    },
    {
        "figure_id": "Figure S3",
        "status": "supplement",
        "asset_path": "path_analysis/results/te_genome_organismal_primary_mediumplus_best_model.pdf",
        "purpose": "Observed-only TE plus organismal best model.",
        "reason": "Documents the body-size result directly in the organismal family.",
    },
    {
        "figure_id": "Figure S4",
        "status": "supplement",
        "asset_path": "results/figures/pca/superfamily_primary_ppca_phylomorphospace.png",
        "purpose": "Phylogenetic PCA view of the main superfamily ordination.",
        "reason": "Keeps the phylogenetically structured ordination available without making it the main visual.",
    },
    {
        "figure_id": "Figure S5",
        "status": "supplement",
        "asset_path": "results/figures/pca/superfamily_primary_scree_plot.png",
        "purpose": "Variance explained for the primary superfamily PCA.",
        "reason": "Supports the ordination methods and PC1/PC2 interpretation.",
    },
    {
        "figure_id": "Figure S6",
        "status": "supplement",
        "asset_path": "results/figures/phylo_signal/blombergs_k_signal.png",
        "purpose": "Trait-wise Blomberg's K summary from the focused phylogenetic-signal workflow.",
        "reason": "Documents broader TE signal patterns without crowding the main text.",
    },
    {
        "figure_id": "Figure S7",
        "status": "supplement",
        "asset_path": "results/figures/phylo_signal/phylogenetic_correlogram_moran_per_trait.png",
        "purpose": "Distance-binned Moran's I correlogram summary.",
        "reason": "Extends the signal results beyond global K and lambda statistics.",
    },
    {
        "figure_id": "Figure S8",
        "status": "supplement",
        "asset_path": "results/figures/landscape/te_order_distribution.png",
        "purpose": "Across-species TE order landscape summary.",
        "reason": "Provides a compact entry point to the full landscape rebuild without showing dozens of per-species curves.",
    },
    {
        "figure_id": "Figure S9",
        "status": "supplement",
        "asset_path": "results/figures/phylo_landscape/phylogeny_with_landscape.png",
        "purpose": "Phylogeny-linked landscape heatmap summary.",
        "reason": "Connects the landscape rebuild to the comparative tree in a single supplemental graphic.",
    },
]

SUPPLEMENT_TABLES = [
    {
        "table_id": "Table S1",
        "asset_path": "paper_freeze/manuscript_primary_results_summary.csv",
        "purpose": "Primary and robustness path-model winner summary.",
    },
    {
        "table_id": "Table S2",
        "asset_path": "results/tables/pca/te_pca_analysis_manifest.csv",
        "purpose": "Canonical PCA analysis manifest and variance explained.",
    },
    {
        "table_id": "Table S3",
        "asset_path": "results/tables/pca/te_pca_phylogenetic_signal.csv",
        "purpose": "Phylogenetic signal statistics for PCA axes.",
    },
    {
        "table_id": "Table S4",
        "asset_path": "paper_freeze/key_file_manifest.csv",
        "purpose": "Checksummed freeze manifest for critical inputs and derived tables.",
    },
]


def read_csv_rows(path: Path) -> List[Dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def read_panel_counts() -> Dict[str, Dict[str, str]]:
    path = PROJECT_ROOT / "path_analysis/data/derived/analysis_panel_summary.csv"
    return {row["panel_name"]: row for row in read_csv_rows(path)}


def build_primary_summary_rows() -> List[Dict[str, str]]:
    panel_counts = read_panel_counts()
    rows: List[Dict[str, str]] = []
    for spec in PRIMARY_RESULT_SPECS:
        ranking_path = PROJECT_ROOT / "path_analysis/results" / f"{spec['panel_id']}_model_ranking.csv"
        ranking_rows = read_csv_rows(ranking_path)
        top = ranking_rows[0]
        runner = ranking_rows[1] if len(ranking_rows) > 1 else None
        panel_meta = panel_counts[spec["panel_id"]]
        delta_to_runner = ""
        if runner is not None:
            delta_to_runner = f"{float(runner['CICc']) - float(top['CICc']):.6f}"
        rows.append(
            {
                "family_id": spec["family_id"],
                "panel_id": spec["panel_id"],
                "subset": spec["subset"],
                "tier": spec["tier"],
                "result_role": spec["result_role"],
                "n_species": panel_meta["n_species"],
                "winner_model": top["model"],
                "winner_CICc": top["CICc"],
                "winner_weight": top["w"],
                "runner_up_model": runner["model"] if runner else "",
                "runner_up_CICc": runner["CICc"] if runner else "",
                "delta_CICc_to_runner_up": delta_to_runner,
                "best_model_pdf": f"path_analysis/results/{spec['panel_id']}_best_model.pdf",
                "average_model_pdf": f"path_analysis/results/{spec['panel_id']}_average_model.pdf",
                "core_claim": spec["core_claim"],
            }
        )
    return rows


def render_figure_plan(summary_rows: List[Dict[str, str]]) -> str:
    medium_rows = [row for row in summary_rows if row["subset"] == "mediumplus"]
    lines: List[str] = []
    lines.append("# Manuscript Figure Plan")
    lines.append("")
    lines.append(
        "This file translates the frozen TE analysis state into a practical main-text versus supplement asset plan."
    )
    lines.append("")
    lines.append("## Main Text Figures")
    lines.append("")
    for item in MAIN_TEXT_FIGURES:
        lines.append(f"- `{item['figure_id']}`: `{item['asset_path']}`")
        lines.append(f"  Purpose: {item['purpose']}")
        lines.append(f"  Why it belongs in the main text: {item['reason']}")
    lines.append("")
    lines.append("## Supplementary Figures")
    lines.append("")
    for item in SUPPLEMENT_FIGURES:
        lines.append(f"- `{item['figure_id']}`: `{item['asset_path']}`")
        lines.append(f"  Purpose: {item['purpose']}")
        lines.append(f"  Why it belongs in the supplement: {item['reason']}")
    lines.append("")
    lines.append("## Supplementary Tables")
    lines.append("")
    for item in SUPPLEMENT_TABLES:
        lines.append(f"- `{item['table_id']}`: `{item['asset_path']}`")
        lines.append(f"  Purpose: {item['purpose']}")
    lines.append("")
    lines.append("## Recommended Narrative Order")
    lines.append("")
    lines.append(
        "- Start with the baseline TE family and the observed-only medium-plus subset as the clean primary comparative result."
    )
    lines.append(
        "- Then show the combined ectopic-plus-organismal family to make the key updated point: body size persists while ectopic_index drops out."
    )
    lines.append(
        "- Use the superfamily compositional PCA as the single TE ordination figure, but keep it framed as supplementary structure rather than the primary causal result."
    )
    lines.append("")
    lines.append("## Primary Result Snapshot")
    lines.append("")
    for row in medium_rows:
        lines.append(
            f"- `{row['panel_id']}`: winner `{row['winner_model']}` at `n = {row['n_species']}`; {row['core_claim']}"
        )
    lines.append("")
    lines.append("## Exclusions From Main Text")
    lines.append("")
    lines.append(
        "- Morphology-linked families remain sensitivity-only because the genome-size side is still provisional."
    )
    lines.append(
        "- Strict-body reruns should stay supplemental because they are robustness checks on species inclusion, not the primary observed-only dataset."
    )
    lines.append(
        "- Order-level all-feature ordinations and their pPCA variants should stay supplemental because they are less interpretable than the superfamily primary ordination."
    )
    lines.append("")
    return "\n".join(lines) + "\n"


def write_csv(path: Path, fieldnames: List[str], rows: List[Dict[str, str]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    summary_rows = build_primary_summary_rows()
    fieldnames = list(summary_rows[0].keys())
    write_csv(OUTPUT_DIR / "manuscript_primary_results_summary.csv", fieldnames, summary_rows)
    (OUTPUT_DIR / "MANUSCRIPT_FIGURE_PLAN.md").write_text(
        render_figure_plan(summary_rows),
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
