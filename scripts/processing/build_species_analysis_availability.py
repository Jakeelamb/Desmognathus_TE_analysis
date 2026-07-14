#!/usr/bin/env python3
"""Build one species-by-analysis-step data availability table.

The table uses the union of the declared TE34 and Cell21 panels (37 species).
High-level step flags are backed by the exact current analysis tables, not by
assumed panel membership alone.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
OUTPUT = (
    ROOT
    / "results/data/research_review/species_analysis_data_availability_union37_v1.csv"
)


def canonical_species(values: pd.Series) -> pd.Series:
    return (
        values.astype(str)
        .str.strip()
        .str.replace("Desmognathus ", "", regex=False)
        .str.replace("D. ", "", regex=False)
        .str.replace("D.", "", regex=False)
        .str.lower()
    )


def read_species(path: str, column: str = "species") -> pd.DataFrame:
    frame = pd.read_csv(ROOT / path)
    frame[column] = canonical_species(frame[column])
    return frame


def species_set(path: str, column: str = "species") -> set[str]:
    return set(read_species(path, column)[column])


def build_availability() -> pd.DataFrame:
    declared = read_species("path_analysis/data/templates/analysis_species_panels.csv")
    readiness = read_species("path_analysis/data/derived/analysis_species_readiness.csv")
    te_panel = read_species(
        "path_analysis/data/derived/panels/study_te_resource_panel34_v1.csv"
    )

    repeatmasker = read_species(
        "results/data/corrected/te34/repeatmasker_hit_inventory_te34_v1.csv"
    )[["species", "hit_count", "aligned_hit_bp"]]
    dnapipete = read_species(
        "results/data/corrected/te34_replicate_averaged/"
        "dnapipete_mass_accounting_te34_replicate_averaged_v1.csv"
    )[["species", "n_dnapipete_runs", "total_aligned_bases"]]
    diversity_species = species_set(
        "results/data/corrected/te34/te_diversity_mass_sensitivity_te34_v1.csv"
    )
    pca_species = species_set(
        "results/data/corrected/te34/te_pca_scores_te34_v1.csv"
    )

    scatter = read_species(
        "results/data/research_review/te_shannon_genome_size_match_audit_te34_v1.csv",
        "canonical_species",
    ).rename(columns={"canonical_species": "species"})
    scatter = scatter[["species", "included_in_scatter"]].drop_duplicates("species")

    ltr = read_species(
        "results/data/corrected/ectopic_ltr30/ectopic_resource_coverage_ltr30_v1.csv"
    )[[
        "species",
        "resource_status",
        "n_tabout_elements",
        "n_selected_5_6_domain_elements",
        "n_usable_elements",
        "n_source_corrupt_exclusions",
    ]]
    cell = read_species(
        "results/data/corrected/cell21/cell_linked_traits_cell21_v1.csv"
    )[[
        "species",
        "n_specimens_strict",
        "n_images_strict",
        "n_pairs_strict",
        "linked_support_label",
        "iod_n_images",
        "iod_n_specimens",
        "iod_n_nuclei",
        "result_status",
    ]]
    genome = read_species(
        "results/data/research_review/genome_size_estimation/"
        "phylogeny_genome_nucleus_cell_summary.csv"
    )[[
        "species",
        "genome_panel_status",
        "n_genome_nuclei",
        "n_genome_images",
        "n_size_cells",
        "n_size_images",
        "genome_size_pg_fuscus_anchored",
    ]]
    all_finalized_path_species = species_set(
        "results/data/research_review/cell_nucleus_genome_path/"
        "primary_traits_all_finalized.csv"
    )
    integrated = read_species(
        "results/data/corrected/path_analysis/"
        "corrected_path_input_sensitivity_analysis18_v1.csv"
    )
    integrated_summary = (
        integrated.groupby("species", as_index=False)
        .agg(
            step08_input_specifications=("species", "size"),
            step08_ectopic_n_elements=("ectopic_n_elements", "max"),
            step08_terminal_internal_proxy_available=(
                "terminal_internal_ratio_median",
                lambda values: values.notna().any(),
            ),
        )
    )

    table = declared.merge(
        readiness[["species", "has_tree_tip"]], on="species", how="left", validate="one_to_one"
    ).merge(
        te_panel[["species", "te_sra_accession", "te_assembly_accession"]],
        on="species",
        how="left",
        validate="one_to_one",
    )
    for source in (repeatmasker, dnapipete, scatter, ltr, cell, genome, integrated_summary):
        table = table.merge(source, on="species", how="left", validate="one_to_one")

    table.insert(1, "scientific_name", "Desmognathus " + table["species"])
    table = table.rename(
        columns={
            "te_resource_panel34_v1": "panel_te_resource34",
            "cell_linked_panel21_v1": "panel_cell_linked21",
            "integrated_path_panel18_v1": "panel_integrated_path18",
            "hit_count": "step03_repeatmasker_hit_count",
            "aligned_hit_bp": "step03_repeatmasker_aligned_hit_bp",
            "n_dnapipete_runs": "step03_dnapipete_runs",
            "total_aligned_bases": "step03_dnapipete_total_aligned_bases",
            "included_in_scatter": "step03_te_genome_scatter_included",
            "resource_status": "step04_ltr_resource_status",
            "n_tabout_elements": "step04_ltr_tabout_elements",
            "n_selected_5_6_domain_elements": "step04_ltr_selected_elements",
            "n_usable_elements": "step04_ltr_usable_elements",
            "n_source_corrupt_exclusions": "step04_ltr_source_corrupt_exclusions",
            "n_specimens_strict": "step05_cell_specimens",
            "n_images_strict": "step05_cell_images",
            "n_pairs_strict": "step05_cell_nucleus_pairs",
            "linked_support_label": "step05_cell_support",
            "iod_n_images": "step05_relative_iod_images",
            "iod_n_specimens": "step05_relative_iod_specimens",
            "iod_n_nuclei": "step05_relative_iod_nuclei",
            "result_status": "step05_relative_iod_status",
            "genome_panel_status": "step06_genome_panel_status",
            "n_genome_nuclei": "step06_genome_nuclei",
            "n_genome_images": "step06_genome_images",
            "n_size_cells": "step05_finalized_size_cells",
            "n_size_images": "step05_finalized_size_images",
            "genome_size_pg_fuscus_anchored": "step06_genome_size_pg_fuscus_anchored",
        }
    )

    table["step01_phylogeny_tree_trimming"] = (
        table["panel_integrated_path18"] & table["has_tree_tip"].fillna(False)
    )
    table["step02_panel_and_resource_provenance"] = True
    table["step03_te_repeat_analysis"] = (
        table["species"].isin(diversity_species)
        & table["species"].isin(pca_species)
        & table["step03_repeatmasker_hit_count"].notna()
        & table["step03_dnapipete_runs"].notna()
    )
    table["step03_te_diversity_available"] = table["species"].isin(diversity_species)
    table["step03_te_pca_available"] = table["species"].isin(pca_species)
    table["step04_ltr_deletion_footprint"] = table["step04_ltr_usable_elements"].fillna(0).gt(0)
    table["step05_cell_modeling_and_measurement"] = table[
        "step05_finalized_size_cells"
    ].notna()
    table["step06_genome_size_estimation"] = table[
        "step06_genome_size_pg_fuscus_anchored"
    ].notna()
    table["step07_cell_nucleus_genome_path"] = table["species"].isin(
        all_finalized_path_species
    )
    table["step08_integrated_phylogenetic_path"] = table[
        "step08_input_specifications"
    ].notna()
    table["step03_te_genome_scatter_included"] = table[
        "step03_te_genome_scatter_included"
    ].astype("boolean").fillna(False).astype(bool)
    table["step08_terminal_internal_proxy_available"] = table[
        "step08_terminal_internal_proxy_available"
    ].astype("boolean").fillna(False).astype(bool)

    count_columns = [
        "step03_repeatmasker_hit_count",
        "step03_repeatmasker_aligned_hit_bp",
        "step03_dnapipete_runs",
        "step03_dnapipete_total_aligned_bases",
        "step04_ltr_tabout_elements",
        "step04_ltr_selected_elements",
        "step04_ltr_usable_elements",
        "step04_ltr_source_corrupt_exclusions",
        "step05_cell_specimens",
        "step05_cell_images",
        "step05_cell_nucleus_pairs",
        "step05_relative_iod_images",
        "step05_relative_iod_specimens",
        "step05_relative_iod_nuclei",
        "step05_finalized_size_cells",
        "step05_finalized_size_images",
        "step06_genome_nuclei",
        "step06_genome_images",
        "step08_input_specifications",
        "step08_ectopic_n_elements",
    ]
    table[count_columns] = table[count_columns].astype("Int64")

    step_columns = [
        "step01_phylogeny_tree_trimming",
        "step02_panel_and_resource_provenance",
        "step03_te_repeat_analysis",
        "step04_ltr_deletion_footprint",
        "step05_cell_modeling_and_measurement",
        "step06_genome_size_estimation",
        "step07_cell_nucleus_genome_path",
        "step08_integrated_phylogenetic_path",
    ]
    table["analysis_steps_with_data"] = table.apply(
        lambda row: ";".join(column[:6] for column in step_columns if bool(row[column])),
        axis=1,
    )
    table["analysis_steps_missing"] = table.apply(
        lambda row: ";".join(column[:6] for column in step_columns if not bool(row[column])),
        axis=1,
    )

    front = [
        "species",
        "scientific_name",
        "panel_te_resource34",
        "panel_cell_linked21",
        "panel_integrated_path18",
        "decision_basis",
        "te_sra_accession",
        "te_assembly_accession",
        *step_columns,
    ]
    table = table[front + [column for column in table.columns if column not in front]]
    return table.sort_values("species").reset_index(drop=True)


def main() -> None:
    table = build_availability()
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(OUTPUT, index=False)
    print(f"Wrote {OUTPUT.relative_to(ROOT)} ({len(table)} species)")


if __name__ == "__main__":
    main()
