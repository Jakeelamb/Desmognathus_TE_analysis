#!/usr/bin/env python3
"""Build the auditable GBE publication dataset package from frozen artifacts.

This is a packaging step, not an analysis rerun. It copies declared canonical
CSV artifacts, expands small frozen ``.csv.gz`` tables to plain CSV, converts
the focal phylogeny to reviewer-readable node/edge CSVs, and records hashes,
dimensions, release status, and provenance for every exported file.
"""

from __future__ import annotations

import gzip
import hashlib
import json
import shutil
from dataclasses import dataclass
from pathlib import Path
from typing import Callable

import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
PUBLICATION = ROOT / "Publication"
DATASET_DIR = PUBLICATION / "datasets"
TREE_DIR = DATASET_DIR / "trees"


@dataclass(frozen=True)
class DatasetSpec:
    supplement_id: str
    analysis_section: str
    source: str
    filename: str
    title: str
    description: str
    release_status: str
    manuscript_use: str
    transform: Callable[[pd.DataFrame], pd.DataFrame] | None = None


def morphology_species_summary(frame: pd.DataFrame) -> pd.DataFrame:
    columns = [
        "species",
        "n_size_cells",
        "n_size_images",
        "cell_area_um2",
        "cell_area_um2_ci_low",
        "cell_area_um2_ci_high",
        "nucleus_area_um2",
        "nucleus_area_um2_ci_low",
        "nucleus_area_um2_ci_high",
    ]
    return frame[columns].copy()


def species_analysis_availability(frame: pd.DataFrame) -> pd.DataFrame:
    """Normalize panel flags against the frozen panel membership tables.

    The source union table contains two malformed string values in its
    integrated-panel flag and a historical typo in the microscopy column name.
    Publication exports derive all three flags from the frozen panel tables so
    the declared denominators are machine-checkable rather than hand-repaired.
    """
    panel_sources = {
        "Genomic_TE_resource34": "path_analysis/data/derived/panels/study_te_resource_panel34_v1.csv",
        "Microscopy_cell_linked21": "path_analysis/data/derived/panels/study_cell_linked_panel21_v1.csv",
        "Conserved_integrated_path18": "path_analysis/data/derived/panels/study_integrated_path_panel18_v1.csv",
    }
    normalized = frame.copy()
    normalized = normalized.drop(columns=["Microscopyl_cell_linked21"], errors="ignore")
    for column, source in panel_sources.items():
        members = set(pd.read_csv(ROOT / source)["species"].astype(str))
        normalized[column] = normalized["species"].astype(str).isin(members)
    leading = [
        "species",
        "scientific_name",
        "Genomic_TE_resource34",
        "Microscopy_cell_linked21",
        "Conserved_integrated_path18",
    ]
    remaining = [column for column in normalized.columns if column not in leading]
    return normalized[leading + remaining]


SPECS = (
    DatasetSpec("S01", "Study design", "results/data/research_review/species_analysis_data_availability_union37_v1.csv", "Supplementary_Data_S01_species_analysis_availability.csv", "Species-level analysis availability", "Union of species and the analysis steps for which each has data. Panel flags are normalized against the frozen panel membership tables.", "main_candidate", "Defines analysis-specific denominators.", transform=species_analysis_availability),
    DatasetSpec("S02", "Study design", "path_analysis/data/derived/panels/study_te_resource_panel34_v1.csv", "Supplementary_Data_S02_te_resource_panel34.csv", "Genomic TE resource panel", "Frozen accession-explicit 34-species genomic TE panel.", "main_candidate", "Primary TE descriptive denominator."),
    DatasetSpec("S03", "Study design", "path_analysis/data/derived/panels/study_cell_linked_panel21_v1.csv", "Supplementary_Data_S03_cell_linked_panel21.csv", "Linked-cell panel", "Frozen 21-species microscopy panel.", "main_candidate", "Primary cell and nucleus denominator."),
    DatasetSpec("S04", "Study design", "path_analysis/data/derived/panels/study_integrated_path_panel18_v1.csv", "Supplementary_Data_S04_integrated_path_panel18.csv", "Integrated comparative panel", "Exact 18-species overlap used by the corrected comparative sensitivity analyses.", "main_candidate", "Defines the integrated phylogenetic denominator."),
    DatasetSpec("S05", "Genomic data acquisition and processing", "results/data/research_review/assembly_quality_te34_v1.csv", "Supplementary_Data_S05_assembly_quality.csv", "Assembly accessions and quality", "NCBI accession-linked assembly span, N50, and contiguity metadata.", "main_candidate", "Supports genomic-resource acquisition and QC."),
    DatasetSpec("S06", "Genomic data acquisition and processing", "results/data/research_review/dnapipete_quality_te34_v1.csv", "Supplementary_Data_S06_dnapipete_input_quality.csv", "dnaPipeTE input and accounting QC", "Run counts, aligned bases, classification retention, and retry policy.", "main_candidate", "Supports de novo TE analysis QC."),
    DatasetSpec("S07", "De novo TE analysis", "results/data/corrected/te34_replicate_averaged/dnapipete_mass_accounting_te34_replicate_averaged_v1.csv", "Supplementary_Data_S07_dnapipete_mass_accounting.csv", "dnaPipeTE mass accounting", "Species-level retained and unresolved aligned-base mass with equal weighting of retry runs.", "main_candidate", "Primary dnaPipeTE denominator audit."),
    DatasetSpec("S08", "TE diversity", "results/data/corrected/te34_replicate_averaged/te_diversity_mass_sensitivity_te34_replicate_averaged_v1.csv", "Supplementary_Data_S08_te_diversity.csv", "TE diversity summaries", "Order- and superfamily-level diversity under classified-conditional and unresolved-mass-aware compositions.", "main_candidate", "Supports richness, Shannon, Gini-Simpson, Hill, and Pielou results."),
    DatasetSpec("S09", "TE ordination", "results/data/corrected/te34_replicate_averaged/te_pca_scores_te34_replicate_averaged_v1.csv", "Supplementary_Data_S09_te_pca_scores.csv", "TE compositional PCA scores", "CLR-PCA species scores for order and superfamily compositions.", "supplementary", "Supports exploratory ordination."),
    DatasetSpec("S10", "TE ordination", "results/data/research_review/te_pca_variance_te34_v1.csv", "Supplementary_Data_S10_te_pca_variance.csv", "TE PCA explained variance", "Variance and cumulative variance by principal component.", "supplementary", "Supports PCA scree and axis labels."),
    DatasetSpec("S11", "TE ordination", "results/data/research_review/te_pca_superfamily_loadings_te34_v1.csv", "Supplementary_Data_S11_te_pca_loadings.csv", "TE PCA loadings", "Superfamily loadings for each principal component.", "supplementary", "Supports biological interpretation of ordination axes."),
    DatasetSpec("S12", "Repeat landscape", "results/data/corrected/te34_replicate_averaged/repeatmasker_divergence_landscape_te34_replicate_averaged_v1.csv", "Supplementary_Data_S12_repeatmasker_divergence_landscape.csv", "RepeatMasker divergence landscape", "Species- and order-level hit counts and aligned bases in one-percentage-point divergence bins.", "main_candidate", "Primary repeat-landscape input."),
    DatasetSpec("S13", "Repeat landscape", "results/data/corrected/te34_replicate_averaged/repeatmasker_hit_inventory_te34_replicate_averaged_v1.csv", "Supplementary_Data_S13_repeatmasker_species_inventory.csv", "RepeatMasker species inventory", "Species-level hit counts, aligned bases, and detected TE orders.", "main_candidate", "Supports RepeatMasker coverage/accounting."),
    DatasetSpec("S14", "LTR deletion-footprint proxy", "results/data/corrected/ectopic_ltr30/ectopic_element_metrics_ltr30_v1.csv", "Supplementary_Data_S14_ltr_terminal_internal_elements.csv", "Element-level terminal:internal depth metrics", "Zero-aware LTR terminal and internal depth summaries with coverage and provenance for each usable element.", "sensitivity_only", "Exploratory deletion-footprint proxy; not an ectopic-recombination rate."),
    DatasetSpec("S15", "LTR deletion-footprint proxy", "results/data/corrected/ectopic_ltr30/ectopic_species_robustness_ltr30_v1.csv", "Supplementary_Data_S15_ltr_terminal_internal_species_robustness.csv", "Species-level terminal:internal robustness", "Median, geometric mean, bootstrap, trimmed, winsorized, and leave-one-element-out summaries.", "sensitivity_only", "Robustness analysis for the exploratory proxy."),
    DatasetSpec("S16", "LTR deletion-footprint proxy", "results/data/corrected/ectopic_ltr30/ectopic_resource_coverage_ltr30_v1.csv", "Supplementary_Data_S16_ltr_resource_coverage.csv", "LTR resource coverage", "Per-species selected-element and usable-depth-file accounting.", "supplementary", "Defines LTR analysis support."),
    DatasetSpec("S17", "LTR deletion-footprint proxy", "results/data/corrected/ectopic_ltr30/ectopic_excluded_elements_ltr30_v1.csv", "Supplementary_Data_S17_ltr_excluded_elements.csv", "Excluded LTR elements", "Source-corrupt depth files excluded under the frozen policy.", "audit_only", "Reviewer audit trail."),
    DatasetSpec("S18", "Cell and nucleus morphometry", "path_analysis/data/external/derived/largest_cell_mask_review/frozen_largest_cell_mask_top50.csv.gz", "Supplementary_Data_S18_frozen_cell_nucleus_objects.csv", "Frozen reviewed cell-nucleus objects", "The exact 1,050 accepted cell/nucleus pairs: 50 per species.", "main_candidate", "Object-level source for cell and nucleus area estimates."),
    DatasetSpec("S19", "Cell and nucleus morphometry", "results/data/research_review/genome_size_estimation/phylogeny_genome_nucleus_cell_summary.csv", "Supplementary_Data_S19_cell_nucleus_species_estimates.csv", "Cell and nucleus species estimates", "Species medians and conditional bootstrap intervals from the frozen 50-pair panels.", "main_candidate", "Primary morphology species table.", transform=morphology_species_summary),
    DatasetSpec("S20", "Cell and nucleus morphometry", "path_analysis/data/external/derived/largest_cell_mask_review/frozen_largest_cell_mask_audit.csv", "Supplementary_Data_S20_cell_mask_freeze_audit.csv", "Cell-mask selection audit", "Ranked review queue, human decisions, cutoff, and frozen-selection flag.", "audit_only", "Reviewer reconstruction of replacements and exclusions."),
    DatasetSpec("S21", "Nuclear IOD", "path_analysis/data/external/derived/image_quality_matched_genome_iod/image_quality_matched_nuclei_frozen_reviewed.csv.gz", "Supplementary_Data_S21_frozen_nuclear_iod_objects.csv", "Frozen reviewed nuclear-IOD objects", "Exact reviewed nuclei retained by the image-quality-matched panel.", "sensitivity_only", "Object-level source for relative nuclear IOD."),
    DatasetSpec("S22", "Nuclear IOD", "results/data/research_review/genome_size_estimation/species_relative_genome_iod_summary.csv", "Supplementary_Data_S22_relative_nuclear_iod_species.csv", "Relative nuclear-IOD species estimates", "Species/fuscus nuclear-IOD ratios and conditional uncertainty under frozen aggregation choices.", "sensitivity_only", "Approved as relative IOD only; not validated absolute genome size."),
    DatasetSpec("S23", "Nuclear IOD", "results/data/research_review/genome_size_estimation/image_relative_genome_iod_summary.csv", "Supplementary_Data_S23_relative_nuclear_iod_by_image.csv", "Image-level relative nuclear IOD", "Image-specific IOD estimates and ratios used to diagnose between-image variation.", "sensitivity_only", "Image-level assay audit."),
    DatasetSpec("S24", "Nuclear IOD", "results/data/research_review/genome_size_estimation/frozen_quality_balance.csv", "Supplementary_Data_S24_nuclear_iod_quality_balance.csv", "Nuclear-IOD quality balance", "Species-level image-quality and sample-support summary for the frozen panel.", "audit_only", "Assay quality audit."),
    DatasetSpec("S25", "Nuclear IOD", "results/data/research_review/genome_size_estimation/iod_quality_residual_diagnostics.csv", "Supplementary_Data_S25_nuclear_iod_quality_diagnostics.csv", "Nuclear-IOD quality diagnostics", "Residual technical-quality associations used to assess image confounding.", "audit_only", "Assay quality audit."),
    DatasetSpec("S26", "Pairwise comparative analysis", "results/data/research_review/genome_size_estimation/phylogeny_genome_nucleus_cell_correlations.csv", "Supplementary_Data_S26_pairwise_pgls_sensitivity.csv", "Pairwise PGLS sensitivity", "Raw correlations and Pagel-lambda PGLS fits among conditional fuscus-anchored IOD, nucleus area, and cell area.", "sensitivity_only", "Conditional calibration; do not present as independently validated genome size."),
    DatasetSpec("S27", "Phylogeny", "results/data/corrected/phylogeny/phylogeny_tip_crosswalk_analysis18_v1.csv", "Supplementary_Data_S27_phylogeny_tip_crosswalk.csv", "Phylogeny tip crosswalk", "Exact taxon matching, resource identifiers, sample-stream independence, and terminal padding.", "main_candidate", "Defines the comparative tree mapping."),
    DatasetSpec("S28", "Phylogeny", "results/data/corrected/phylogeny/phylogeny_tree_metrics_analysis18_v1.csv", "Supplementary_Data_S28_phylogeny_tree_metrics.csv", "Phylogeny metrics", "Topology, branch-length, ultrametricity, and derivation checks for focal trees.", "main_candidate", "Tree validation."),
    DatasetSpec("S29", "Phylogeny", "results/data/corrected/phylogeny/phylogeny_tree_uncertainty_metrics_analysis18_v1.csv", "Supplementary_Data_S29_phylogeny_tree_uncertainty.csv", "Published-tree uncertainty", "Metrics for the focal, published main, and 200 published bootstrap trees.", "supplementary", "Tree-sensitivity audit."),
    DatasetSpec("S30", "Phylogenetic path sensitivity", "results/data/corrected/path_analysis/corrected_path_anchor_model_comparison_analysis18_v1.csv", "Supplementary_Data_S30_path_model_comparison.csv", "Path-model comparison", "Global-fit, CICc, delta, weights, release gates, and warnings for anchor specifications.", "sensitivity_only", "Exploratory association structures only."),
    DatasetSpec("S31", "Phylogenetic path sensitivity", "results/data/corrected/path_analysis/corrected_path_anchor_edges_analysis18_v1.csv", "Supplementary_Data_S31_path_edge_estimates.csv", "Path edge estimates", "Standardized path coefficients and approximate intervals for top anchor models.", "sensitivity_only", "No causal or absolute-genome-size claim."),
    DatasetSpec("S32", "Phylogenetic path sensitivity", "results/data/corrected/path_analysis/corrected_path_model_stability_analysis18_v1.csv", "Supplementary_Data_S32_path_model_stability.csv", "Path model stability", "Modal models, global-fit rates, and model-weight stability across sensitivity fits.", "sensitivity_only", "Exploratory stability summary."),
    DatasetSpec("S33", "Phylogenetic path sensitivity", "results/data/corrected/path_analysis/corrected_path_tree_sensitivity_rankings_analysis18_v1.csv", "Supplementary_Data_S33_path_tree_sensitivity.csv", "Path tree sensitivity", "Rankings across focal, published, and bootstrap phylogenies.", "sensitivity_only", "Tree uncertainty audit."),
    DatasetSpec("S34", "Phylogenetic path sensitivity", "results/data/corrected/path_analysis/corrected_path_input_sensitivity_analysis18_v1.csv", "Supplementary_Data_S34_path_measurement_sensitivity.csv", "Path measurement sensitivity", "Model results across declared morphology and relative-IOD specifications.", "sensitivity_only", "Measurement-model audit."),
    DatasetSpec("S35", "Phylogenetic path sensitivity", "results/data/corrected/path_analysis/corrected_path_leave_one_out_rankings_analysis18_v1.csv", "Supplementary_Data_S35_path_leave_one_species_out.csv", "Path leave-one-species-out analysis", "Model rankings after omitting each species in turn.", "sensitivity_only", "Influence audit."),
    DatasetSpec("S36", "Phylogenetic path sensitivity", "results/data/corrected/path_analysis/corrected_path_simulation_calibration_analysis18_v1.csv", "Supplementary_Data_S36_path_simulation_calibration.csv", "Path simulation calibration", "Actual-tree false-selection and recovery rates under null and observed-chain simulations.", "audit_only", "Calibration gate; blocks confirmatory causal interpretation."),
    DatasetSpec("S37", "Publication release", "results/data/corrected/path_analysis/publication_release_gate_matrix_analysis18_v1.csv", "Supplementary_Data_S37_publication_release_gates.csv", "Publication release gates", "Evidence-layer status and permitted manuscript language.", "audit_only", "Claim-control table."),
)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def export_csv(spec: DatasetSpec) -> dict[str, object]:
    source = ROOT / spec.source
    if not source.exists():
        raise FileNotFoundError(source)
    destination = DATASET_DIR / spec.filename
    if spec.transform is not None:
        frame = spec.transform(pd.read_csv(source, low_memory=False))
        frame.to_csv(destination, index=False)
    elif source.suffix == ".gz":
        with gzip.open(source, "rb") as source_handle, destination.open("wb") as output_handle:
            shutil.copyfileobj(source_handle, output_handle)
        frame = pd.read_csv(destination, low_memory=False)
    else:
        shutil.copy2(source, destination)
        frame = pd.read_csv(destination, low_memory=False)
    return {
        "supplement_id": spec.supplement_id,
        "analysis_section": spec.analysis_section,
        "filename": spec.filename,
        "file_format": "CSV",
        "title": spec.title,
        "description": spec.description,
        "release_status": spec.release_status,
        "manuscript_use": spec.manuscript_use,
        "source_path": spec.source,
        "source_sha256": sha256_file(source),
        "output_sha256": sha256_file(destination),
        "rows": int(len(frame)),
        "columns": int(len(frame.columns)),
        "bytes": int(destination.stat().st_size),
    }, frame


def export_tree_tables() -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    try:
        from Bio import Phylo
    except ModuleNotFoundError as error:
        raise ModuleNotFoundError("Biopython is required for tree CSV export") from error

    source = ROOT / "results/data/corrected/phylogeny/desmognathus_time_tree_analysis18_v1.nwk"
    tree = Phylo.read(source, "newick")
    node_ids: dict[object, str] = {}
    internal_index = 0
    for clade in tree.find_clades(order="preorder"):
        if clade.is_terminal():
            node_ids[clade] = str(clade.name)
        else:
            internal_index += 1
            node_ids[clade] = f"internal_{internal_index:02d}"
    edge_rows = []
    for parent in tree.find_clades(order="preorder"):
        for child in parent.clades:
            edge_rows.append(
                {
                    "parent_node": node_ids[parent],
                    "child_node": node_ids[child],
                    "child_is_tip": child.is_terminal(),
                    "branch_length_myr": child.branch_length,
                }
            )
    depth = tree.depths()
    node_rows = [
        {
            "node": node_ids[clade],
            "is_tip": clade.is_terminal(),
            "tip_label": clade.name if clade.is_terminal() else "",
            "root_distance_myr": depth[clade],
        }
        for clade in tree.find_clades(order="preorder")
    ]
    edge_frame = pd.DataFrame(edge_rows)
    node_frame = pd.DataFrame(node_rows)
    edge_name = "Supplementary_Data_S38_phylogeny_edges.csv"
    node_name = "Supplementary_Data_S39_phylogeny_nodes.csv"
    edge_frame.to_csv(DATASET_DIR / edge_name, index=False)
    node_frame.to_csv(DATASET_DIR / node_name, index=False)

    non_csv_name = "desmognathus_time_tree_analysis18_v1.nwk"
    shutil.copy2(source, TREE_DIR / non_csv_name)
    bootstrap_source = ROOT / "results/data/corrected/phylogeny/stewart_wiens_2025_bootstrap_analysis18_v1.nex"
    main_source = ROOT / "results/data/corrected/phylogeny/stewart_wiens_2025_main_analysis18_v1.nwk"
    shutil.copy2(bootstrap_source, TREE_DIR / bootstrap_source.name)
    shutil.copy2(main_source, TREE_DIR / main_source.name)

    manifests = []
    for supplement_id, filename, title, frame in [
        ("S38", edge_name, "Focal phylogeny edge list", edge_frame),
        ("S39", node_name, "Focal phylogeny node table", node_frame),
    ]:
        destination = DATASET_DIR / filename
        manifests.append(
            {
                "supplement_id": supplement_id,
                "analysis_section": "Phylogeny",
                "filename": filename,
                "file_format": "CSV",
                "title": title,
                "description": "Reviewer-readable CSV representation of the exact focal Newick tree.",
                "release_status": "main_candidate",
                "manuscript_use": "Tree reconstruction and audit.",
                "source_path": str(source.relative_to(ROOT)),
                "source_sha256": sha256_file(source),
                "output_sha256": sha256_file(destination),
                "rows": int(len(frame)),
                "columns": int(len(frame.columns)),
                "bytes": int(destination.stat().st_size),
            }
        )
    non_csv = []
    for path in [TREE_DIR / non_csv_name, TREE_DIR / main_source.name, TREE_DIR / bootstrap_source.name]:
        non_csv.append(
            {
                "filename": str(path.relative_to(DATASET_DIR)),
                "format": path.suffix.lstrip(".").upper(),
                "reason_not_csv": "Newick/NEXUS preserves standard phylogenetic structure and bootstrap tree sets.",
                "sha256": sha256_file(path),
                "bytes": path.stat().st_size,
            }
        )
    return manifests, non_csv


def main() -> None:
    DATASET_DIR.mkdir(parents=True, exist_ok=True)
    TREE_DIR.mkdir(parents=True, exist_ok=True)
    manifest_rows: list[dict[str, object]] = []
    dictionary_rows: list[dict[str, object]] = []
    for spec in SPECS:
        manifest, frame = export_csv(spec)
        manifest_rows.append(manifest)
        for column in frame.columns:
            dictionary_rows.append(
                {
                    "supplement_id": spec.supplement_id,
                    "filename": spec.filename,
                    "column": column,
                    "pandas_dtype": str(frame[column].dtype),
                    "source_column": column,
                    "definition_status": "source_named; manuscript-critical fields require curated definitions",
                }
            )
    tree_manifests, non_csv = export_tree_tables()
    manifest_rows.extend(tree_manifests)
    manifest = pd.DataFrame(manifest_rows).sort_values("supplement_id")
    dictionary = pd.DataFrame(dictionary_rows).sort_values(["supplement_id", "column"])
    manifest.to_csv(DATASET_DIR / "DATASET_MANIFEST.csv", index=False)
    dictionary.to_csv(DATASET_DIR / "COLUMN_INVENTORY.csv", index=False)
    (DATASET_DIR / "NON_CSV_FILES.json").write_text(
        json.dumps(non_csv, indent=2) + "\n",
        encoding="utf-8",
    )
    summary = {
        "publication_dataset_count": int(len(manifest)),
        "csv_files": int(len(manifest)),
        "non_csv_tree_files": int(len(non_csv)),
        "total_csv_rows": int(manifest["rows"].sum()),
        "total_csv_bytes": int(manifest["bytes"].sum()),
        "release_status_counts": manifest["release_status"].value_counts().to_dict(),
        "manifest_sha256": sha256_file(DATASET_DIR / "DATASET_MANIFEST.csv"),
        "column_inventory_sha256": sha256_file(DATASET_DIR / "COLUMN_INVENTORY.csv"),
    }
    (DATASET_DIR / "release_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
