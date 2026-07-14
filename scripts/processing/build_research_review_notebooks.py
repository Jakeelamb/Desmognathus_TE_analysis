#!/usr/bin/env python3
"""Build the eight canonical collaborator-review notebooks.

Six notebooks are thin review surfaces over frozen, versioned products. The
genome-size and cell/nucleus/genome path notebooks retain their complete
lightweight analysis and presentation code. None launches dnaPipeTE,
RepeatMasker, RepeatModeler, segmentation inference/training, or a new
phylogenetic path-model fit.
"""

from __future__ import annotations

import csv
import hashlib
import json
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable, Optional, Sequence

import nbformat as nbf


ROOT = Path(__file__).resolve().parents[2]
OUTPUT_DIR = ROOT / "notebooks/research_review"
REVIEW_DATA_DIR = ROOT / "results/data/research_review"
REVIEW_FIGURE_DIR = ROOT / "results/figures/research_review"
FROZEN_INPUT_DIR = REVIEW_DATA_DIR / "frozen_inputs"
FROZEN_REGISTRY = REVIEW_DATA_DIR / "frozen_input_registry.csv"
FIGURE_BUILD_COMMAND = (
    "scripts/run_in_dusky.sh Rscript "
    "scripts/processing/build_audited_historical_style_figures.R"
)

NOTEBOOKS = {
    "01_phylogeny_tree_trimming.ipynb": "phylogeny_notebook",
    "02_data_tables_and_provenance.ipynb": "data_tables_notebook",
    "03_repeat_analysis_te34.ipynb": "repeat_notebook",
    "04_ltr_deletion_footprint.ipynb": "ltr_notebook",
    "05_cell_modeling_and_measurement.ipynb": "cell_notebook",
    "06_genome_size_estimation.ipynb": "genome_size_notebook",
    "07_cell_nucleus_genome_path_analysis.ipynb": "cell_nucleus_genome_path_notebook",
    "08_integrated_phylogenetic_path_analysis.ipynb": "integrated_path_notebook",
}

# Small data, tree, and manifest artifacts required to re-execute the review
# notebooks or regenerate their R figures.  Snapshots preserve the original
# repository-relative path below frozen_inputs, making every mapping unambiguous.
FROZEN_DATA_SOURCES = (
    "image_quality_matched_all_species_decisions.csv",
    "image_quality_matched_replacement_review_decisions.csv",
    "input_data/specimen_slides/image_quality_matched_review_extension_20260714_decisions_aeneus_wrighti_orestes_ochrophaeus.csv",
    "input_data/specimen_slides/literal_largest_cell_mask_review_extension_20260714_decisions_aenues_orestes_wrighti.csv",
    "input_data/phylogeny/desmo900dated_test.tre",
    "path_analysis/data/derived/panels/study_cell_linked_panel21_v1.csv",
    "path_analysis/data/derived/panels/study_integrated_path_panel18_v1.csv",
    "path_analysis/data/derived/panels/study_te_resource_panel34_v1.csv",
    "path_analysis/data/external/derived/cellprofiler_final_species_results.csv",
    "path_analysis/data/external/derived/cellprofiler_genome_state_summary.csv",
    "path_analysis/data/external/derived/image_quality_matched_genome_iod/image_quality_matched_nuclei_frozen_reviewed.csv.gz",
    "path_analysis/data/external/derived/image_quality_matched_genome_iod/finalized_quality_matched_nuclei_all_reviewed_species.csv.gz",
    "path_analysis/data/external/derived/image_quality_matched_genome_iod/manifest.json",
    "path_analysis/data/external/derived/largest_cell_mask_review/frozen_largest_cell_mask_top50.csv.gz",
    "path_analysis/data/external/derived/largest_cell_mask_review/finalized_largest_cell_masks_all_reviewed_species.csv.gz",
    "path_analysis/data/external/derived/microscopy_review_finalization_20260714.json",
    "path_analysis/data/derived/source_file_registry.csv",
    "path_analysis/data/templates/source_manifest.csv",
    "results/data/corrected/cell21/cell_linked_traits_cell21_v1.csv",
    "results/data/corrected/microscopy/microscopy_estimand_sensitivity_analysis18_v1.csv",
    "results/data/corrected/microscopy/microscopy_estimator_rank_stability_analysis18_v1.csv",
    "results/data/corrected/microscopy/microscopy_segmentation_validation_analysis18_v1.csv",
    "results/data/corrected/path_analysis/corrected_path_anchor_edges_analysis18_v1.csv",
    "results/data/corrected/path_analysis/corrected_path_anchor_model_comparison_analysis18_v1.csv",
    "results/data/corrected/path_analysis/corrected_path_data_sensitivity_basis_sets_analysis18_v1.csv",
    "results/data/corrected/path_analysis/corrected_path_input_sensitivity_analysis18_v1.csv",
    "results/data/corrected/path_analysis/corrected_path_leave_one_out_rankings_analysis18_v1.csv",
    "results/data/corrected/path_analysis/corrected_path_measurement_specifications_analysis18_v1.csv",
    "results/data/corrected/path_analysis/corrected_path_model_stability_analysis18_v1.csv",
    "results/data/corrected/path_analysis/corrected_path_simulation_calibration_analysis18_v1.csv",
    "results/data/corrected/path_analysis/corrected_path_simulation_edge_calibration_analysis18_v1.csv",
    "results/data/corrected/path_analysis/corrected_path_tree_edge_uncertainty_analysis18_v1.csv",
    "results/data/corrected/path_analysis/publication_release_gate_matrix_analysis18_v1.csv",
    "results/data/corrected/phylogeny/desmognathus_time_tree_analysis18_v1.nwk",
    "results/data/corrected/phylogeny/phylogeny_release_audit_analysis18_v1.manifest.json",
    "results/data/corrected/phylogeny/phylogeny_tip_crosswalk_analysis18_v1.csv",
    "results/data/corrected/phylogeny/phylogeny_topology_difference_analysis18_v1.csv",
    "results/data/corrected/phylogeny/phylogeny_tree_metrics_analysis18_v1.csv",
    "results/data/corrected/phylogeny/phylogeny_tree_uncertainty_metrics_analysis18_v1.csv",
    "results/data/corrected/phylogeny/stewart_wiens_2025_bootstrap_analysis18_v1.nex",
    "results/data/corrected/phylogeny/stewart_wiens_2025_main_analysis18_v1.nwk",
    "results/data/corrected/te34_replicate_averaged/dnapipete_mass_accounting_te34_replicate_averaged_v1.csv",
    "results/data/corrected/te34_replicate_averaged/te34_replicate_averaged_v1.manifest.json",
    "results/data/corrected/te34_replicate_averaged/te_diversity_mass_sensitivity_te34_replicate_averaged_v1.csv",
    "results/data/dnaPipeTE_order_breakdown.csv",
    "results/data/dnaPipeTE_superfamily_breakdown.csv",
    "results/data/corrected/ectopic_ltr30/ectopic_element_metrics_ltr30_v1.csv",
    "results/data/corrected/ectopic_ltr30/ectopic_species_robustness_ltr30_v1.csv",
    "results/data/corrected/ectopic_ltr30/ectopic_excluded_elements_ltr30_v1.csv",
    "results/data/corrected/ectopic_ltr30/ectopic_resource_coverage_ltr30_v1.csv",
    "results/data/corrected/ectopic_ltr30/ectopic_ltr30_v1.manifest.json",
    # Additional inputs read by the lightweight R figure runner.
    "results/data/legacy_descriptive/dnapipete_species_quality_te34_legacy_descriptive_v1.csv",
    "results/data/legacy_descriptive/te_order_composition_te34_legacy_descriptive_v1.csv",
    "results/data/legacy_descriptive/te_superfamily_composition_te34_legacy_descriptive_v1.csv",
    "results/phylogeny/processed_phylogeny.nwk",
)

# Previously generated audit figures used by notebooks 01, 02, 05, and 06.
# Their basenames are unique, so the portable review directory can stay flat.
PORTABLE_FIGURE_SOURCES = (
    "results/figures/corrected/cell21/cell21_cell_nucleus_relative_iod_v1.png",
    "results/figures/corrected/ectopic_ltr30/ectopic_ltr30_coverage_qc_v1.pdf",
    "results/figures/corrected/ectopic_ltr30/ectopic_ltr30_coverage_qc_v1.png",
    "results/figures/corrected/microscopy/microscopy_estimand_sensitivity_analysis18_v1.png",
    "results/figures/corrected/microscopy/microscopy_estimator_rank_stability_analysis18_v1.png",
    "results/figures/corrected/microscopy/microscopy_image_iod_qc_analysis18_v1.png",
    "results/figures/corrected/microscopy/microscopy_iod_decomposition_analysis18_v1.png",
    "results/figures/corrected/microscopy/microscopy_relative_iod_qc_sensitivity_analysis18_v1.png",
    "results/figures/corrected/microscopy/microscopy_segmentation_validation_analysis18_v1.png",
    "results/figures/corrected/microscopy/microscopy_selection_depth_analysis18_v1.png",
    "results/figures/corrected/microscopy/microscopy_support_and_review_analysis18_v1.png",
    "results/figures/corrected/path_analysis/corrected_path_anchor_dag_analysis18_v1.png",
    "results/figures/corrected/path_analysis/corrected_path_anchor_model_comparison_analysis18_v1.png",
    "results/figures/corrected/path_analysis/corrected_path_leave_one_out_influence_analysis18_v1.png",
    "results/figures/corrected/path_analysis/corrected_path_measurement_model_weights_analysis18_v1.png",
    "results/figures/corrected/path_analysis/corrected_path_simulation_calibration_analysis18_v1.png",
    "results/figures/corrected/path_analysis/corrected_path_tree_uncertainty_analysis18_v1.png",
    "results/figures/corrected/path_analysis/publication_release_gate_matrix_analysis18_v1.png",
    "results/figures/corrected/phylogeny/phylogeny_focal_vs_published_topology_analysis18_v1.png",
    "results/figures/corrected/phylogeny/phylogeny_patristic_uncertainty_analysis18_v1.png",
    "results/figures/corrected/phylogeny/phylogeny_processed_scale_comparison_v1.png",
    "results/figures/corrected/phylogeny/phylogeny_published_bootstrap_density_analysis18_v1.png",
    "results/figures/corrected/phylogeny/phylogeny_root_age_uncertainty_analysis18_v1.png",
    "results/figures/presentation_notebooks/cell_models/cell_mask_examples_analysis18_v1.pdf",
    "results/figures/presentation_notebooks/cell_models/cell_mask_examples_analysis18_v1.png",
    "results/figures/presentation_notebooks/cell_models/cell_trait_pairwise_bm_pgls_analysis18_v1.pdf",
    "results/figures/presentation_notebooks/cell_models/cell_trait_pairwise_bm_pgls_analysis18_v1.png",
    "results/figures/presentation_notebooks/cell_models/cell_trait_phylogeny_tracks_analysis18_v1.pdf",
    "results/figures/presentation_notebooks/cell_models/cell_trait_phylogeny_tracks_analysis18_v1.png",
    "results/figures/presentation_notebooks/data_tables/analysis18_pipeline_evidence_streams_v1.pdf",
    "results/figures/presentation_notebooks/data_tables/analysis18_pipeline_evidence_streams_v1.png",
    "results/figures/presentation_notebooks/data_tables/final18_evidence_completeness_analysis18_v1.pdf",
    "results/figures/presentation_notebooks/data_tables/final18_evidence_completeness_analysis18_v1.png",
    "results/figures/presentation_notebooks/phylogeny/phylogeny_source46_vs_final18_trim_analysis18_v1.pdf",
    "results/figures/presentation_notebooks/phylogeny/phylogeny_source46_vs_final18_trim_analysis18_v1.png",
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def git_status(path: str) -> str:
    tracked = subprocess.run(
        ["git", "ls-files", "--error-unmatch", "--", path],
        cwd=ROOT,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
        check=False,
    ).returncode == 0
    if tracked:
        return "tracked"
    ignored = subprocess.run(
        ["git", "check-ignore", "-q", "--", path],
        cwd=ROOT,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
        check=False,
    ).returncode == 0
    return "ignored" if ignored else "untracked"


def existing_registry() -> dict[str, dict[str, str]]:
    if not FROZEN_REGISTRY.exists():
        return {}
    with FROZEN_REGISTRY.open(newline="") as handle:
        return {row["source_path"]: row for row in csv.DictReader(handle)}


def snapshot_artifact(
    source_relative: str,
    frozen_relative: str,
    artifact_kind: str,
    previous: dict[str, dict[str, str]],
    copied_at: str,
) -> dict[str, object]:
    source = ROOT / source_relative
    frozen = ROOT / frozen_relative
    old = previous.get(source_relative, {})

    if source.exists():
        source_hash = sha256(source)
        frozen_hash_before = sha256(frozen) if frozen.exists() else ""
        if frozen_hash_before != source_hash:
            frozen.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(source, frozen)
        frozen_hash = sha256(frozen)
        source_present = True
    elif frozen.exists():
        # A collaborator checkout intentionally contains the frozen artifact even
        # when the ignored upstream result tree is absent.
        frozen_hash = sha256(frozen)
        source_hash = old.get("source_sha256", frozen_hash)
        source_present = False
    else:
        raise FileNotFoundError(
            f"Neither source nor frozen review artifact exists: {source_relative} -> {frozen_relative}"
        )

    if frozen_hash != source_hash:
        raise RuntimeError(f"Frozen artifact hash mismatch for {source_relative}")

    unchanged = (
        old.get("frozen_path") == frozen_relative
        and old.get("source_sha256") == source_hash
        and old.get("frozen_sha256") == frozen_hash
        and bool(old.get("copied_at_utc"))
    )
    return {
        "artifact_kind": artifact_kind,
        "source_path": source_relative,
        "frozen_path": frozen_relative,
        "size_bytes": frozen.stat().st_size,
        "source_sha256": source_hash,
        "frozen_sha256": frozen_hash,
        "copied_at_utc": old["copied_at_utc"] if unchanged else copied_at,
        "source_present_at_build": str(source_present).lower(),
        "source_git_status": git_status(source_relative) if source_present else "source_absent",
    }


def build_frozen_bundle() -> None:
    previous = existing_registry()
    copied_at = datetime.now(timezone.utc).replace(microsecond=0).isoformat()
    rows = []

    for source in FROZEN_DATA_SOURCES:
        frozen = str((FROZEN_INPUT_DIR / source).relative_to(ROOT))
        rows.append(snapshot_artifact(source, frozen, "data_or_tree", previous, copied_at))

    basenames = [Path(source).name for source in PORTABLE_FIGURE_SOURCES]
    if len(basenames) != len(set(basenames)):
        raise RuntimeError("Portable review figure basenames are not unique")
    for source in PORTABLE_FIGURE_SOURCES:
        frozen = str((REVIEW_FIGURE_DIR / Path(source).name).relative_to(ROOT))
        rows.append(snapshot_artifact(source, frozen, "figure", previous, copied_at))

    REVIEW_DATA_DIR.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "artifact_kind",
        "source_path",
        "frozen_path",
        "size_bytes",
        "source_sha256",
        "frozen_sha256",
        "copied_at_utc",
        "source_present_at_build",
        "source_git_status",
    ]
    with FROZEN_REGISTRY.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(sorted(rows, key=lambda row: str(row["source_path"])))

    print(
        f"Verified {len(rows)} portable dependencies in "
        f"{FROZEN_REGISTRY.relative_to(ROOT)}"
    )


def md(text: str):
    return nbf.v4.new_markdown_cell(text.strip())


def code(text: str):
    return nbf.v4.new_code_cell(text.strip())


def common_setup():
    return code(
        f"""
from pathlib import Path
import hashlib
import json

import numpy as np
import pandas as pd
from IPython.display import HTML, IFrame, Image, Markdown, display

ROOT = Path.cwd().resolve()
for candidate in [ROOT, *ROOT.parents]:
    if (candidate / "path_analysis/data/derived").exists() and (candidate / "scripts/processing").exists():
        ROOT = candidate
        break
else:
    raise RuntimeError("Open this notebook from inside the Desmognathus_TE repository")

REVIEW_FIGURE_DIR = ROOT / "results/figures/research_review"
REVIEW_DATA_DIR = ROOT / "results/data/research_review"
FROZEN_INPUT_DIR = REVIEW_DATA_DIR / "frozen_inputs"
FROZEN_REGISTRY_PATH = REVIEW_DATA_DIR / "frozen_input_registry.csv"
FIGURE_BUILD_COMMAND = {FIGURE_BUILD_COMMAND!r}

if not FROZEN_REGISTRY_PATH.exists():
    raise FileNotFoundError(
        f"Portable dependency registry missing: {{FROZEN_REGISTRY_PATH}}. "
        "Rebuild the research-review notebooks first."
    )
frozen_registry = pd.read_csv(FROZEN_REGISTRY_PATH)
required_registry_columns = {{"source_path", "frozen_path", "source_sha256", "frozen_sha256"}}
if not required_registry_columns.issubset(frozen_registry.columns):
    raise RuntimeError("Frozen dependency registry schema is incomplete")
if not frozen_registry.source_sha256.eq(frozen_registry.frozen_sha256).all():
    raise RuntimeError("Frozen dependency registry contains a source/snapshot hash mismatch")
FROZEN_ARTIFACTS = dict(zip(frozen_registry.source_path, frozen_registry.frozen_path))

def file_sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()

missing_snapshots = [path for path in frozen_registry.frozen_path if not (ROOT / path).exists()]
if missing_snapshots:
    raise FileNotFoundError(f"Portable dependency snapshots missing: {{missing_snapshots}}")
corrupt_snapshots = [
    row.frozen_path
    for row in frozen_registry.itertuples()
    if file_sha256(ROOT / row.frozen_path) != row.frozen_sha256
]
if corrupt_snapshots:
    raise RuntimeError(f"Portable dependency snapshot SHA256 mismatch: {{corrupt_snapshots}}")

pd.set_option("display.max_columns", 100)
pd.set_option("display.max_rows", 120)
pd.set_option("display.max_colwidth", 120)

def resolve_artifact(relative):
    path = Path(relative)
    if path.is_absolute():
        return path
    return ROOT / FROZEN_ARTIFACTS.get(str(relative), str(relative))

def read_csv(relative, **kwargs):
    path = resolve_artifact(relative)
    if not path.exists():
        raise FileNotFoundError(path)
    return pd.read_csv(path, **kwargs)

def require_paths(paths, build_command=None):
    missing = [str(path) for path in paths if not Path(path).exists()]
    if missing:
        message = "Missing required frozen review artifacts:\\n- " + "\\n- ".join(missing)
        if build_command:
            message += f"\\n\\nGenerate only the lightweight review outputs with:\\n{{build_command}}"
        raise FileNotFoundError(message)

def show_figure(relative, width=1150):
    path = resolve_artifact(relative)
    require_paths([path])
    display(Image(filename=str(path), width=width))

def show_text(relative, language="text", max_chars=16000):
    path = resolve_artifact(relative)
    require_paths([path])
    text = path.read_text(errors="replace")
    if len(text) > max_chars:
        text = text[:max_chars] + "\\n... [display truncated; open the linked source for the complete file]"
    display(Markdown(f"```{{language}}\\n{{text}}\\n```"))

def artifact_table(relatives):
    rows = []
    for relative in relatives:
        path = resolve_artifact(relative)
        rows.append({{
            "artifact": relative,
            "frozen_artifact": str(path.relative_to(ROOT)),
            "exists": path.exists(),
            "size_mb": round(path.stat().st_size / 1024**2, 3) if path.exists() else np.nan,
        }})
    result = pd.DataFrame(rows)
    if not result.exists.all():
        raise FileNotFoundError(result.loc[~result.exists, "artifact"].tolist())
    return result

print(f"Repository: {{ROOT}}")
print("Mode: frozen-artifact review only; no expensive upstream analysis is run")
print(f"Portable dependency snapshots: {{len(frozen_registry)}}")
"""
    )


def source_cell(paths: Sequence[str], *, show_runner: Optional[str] = None):
    rows = repr([{"role": role, "path": path} for role, path in paths])
    body = f"""
source_files = pd.DataFrame({rows})
source_files["resolved_path"] = source_files.path.map(lambda value: str(resolve_artifact(value).relative_to(ROOT)))
source_files["exists"] = source_files.path.map(lambda value: resolve_artifact(value).exists())
display(source_files)
if not source_files.exists.all():
    raise FileNotFoundError(source_files.loc[~source_files.exists, "path"].tolist())
"""
    if show_runner:
        body += f"\nshow_text({show_runner!r}, language='r')\n"
    return code(body)


def make_notebook(title: str, analysis_id: str, boundary: str, cells: Iterable[object]):
    notebook = nbf.v4.new_notebook()
    notebook["metadata"] = {
        "kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
        "language_info": {"name": "python", "version": "3.8"},
        "desmognathus_analysis_id": analysis_id,
        "artifact_mode": "frozen-read-only",
        "expensive_upstream_tools_executed": False,
    }
    notebook["cells"] = [
        md(
            f"""
# {title}

**{boundary}**

This is a collaborator-review notebook. It loads frozen data, figures, manifests,
and trees and can link to optional local HTML review pages; it does not rerun
expensive bioinformatics, segmentation, model training, or path-model fitting.
Every result is identified by its analysis panel and exact source path.
"""
        ),
        common_setup(),
        *cells,
    ]
    return notebook


def phylogeny_notebook():
    return make_notebook(
        "Phylogenetic tree trimming and validation",
        "path18_v1",
        "Tree review for the integrated path18 panel (n = 18)",
        [
            md(
                """
## Question and decision

How was the time-calibrated tree reduced to the exact integrated panel, and did
the trim preserve the intended topology and branch-time information? The frozen
analysis tree contains exact species matches only. It does not substitute nearby
species, and genomic and microscopy specimens remain independent evidence streams.
"""
            ),
            code(
                """
artifacts = [
    "input_data/phylogeny/desmo900dated_test.tre",
    "results/data/corrected/phylogeny/desmognathus_time_tree_analysis18_v1.nwk",
    "results/data/corrected/phylogeny/stewart_wiens_2025_main_analysis18_v1.nwk",
    "results/data/corrected/phylogeny/stewart_wiens_2025_bootstrap_analysis18_v1.nex",
    "results/data/corrected/phylogeny/phylogeny_tip_crosswalk_analysis18_v1.csv",
    "results/data/corrected/phylogeny/phylogeny_tree_metrics_analysis18_v1.csv",
    "results/data/corrected/phylogeny/phylogeny_release_audit_analysis18_v1.manifest.json",
]
display(artifact_table(artifacts))
panel = read_csv("path_analysis/data/derived/panels/study_integrated_path_panel18_v1.csv")
crosswalk = read_csv("results/data/corrected/phylogeny/phylogeny_tip_crosswalk_analysis18_v1.csv")
metrics = read_csv("results/data/corrected/phylogeny/phylogeny_tree_metrics_analysis18_v1.csv")
assert len(panel) == 18 and panel.species.is_unique
assert len(crosswalk) == 18 and set(crosswalk.species) == set(panel.species)
display(panel[["species", "te_sra_accession", "te_assembly_accession", "decision_basis"]])
display(metrics)
"""
            ),
            md("## Source tree versus exact path18 trim"),
            code(
                """
show_figure("results/figures/presentation_notebooks/phylogeny/phylogeny_source46_vs_final18_trim_analysis18_v1.png")
show_figure("results/figures/corrected/phylogeny/phylogeny_processed_scale_comparison_v1.png")
"""
            ),
            md("## Topology, calibration, and uncertainty checks"),
            code(
                """
for figure in [
    "results/figures/corrected/phylogeny/phylogeny_focal_vs_published_topology_analysis18_v1.png",
    "results/figures/corrected/phylogeny/phylogeny_root_age_uncertainty_analysis18_v1.png",
    "results/figures/corrected/phylogeny/phylogeny_patristic_uncertainty_analysis18_v1.png",
    "results/figures/corrected/phylogeny/phylogeny_published_bootstrap_density_analysis18_v1.png",
]:
    show_figure(figure)
uncertainty = read_csv("results/data/corrected/phylogeny/phylogeny_tree_uncertainty_metrics_analysis18_v1.csv")
topology = read_csv("results/data/corrected/phylogeny/phylogeny_topology_difference_analysis18_v1.csv")
display(uncertainty.describe(include="all").T)
display(topology)
"""
            ),
            md(
                """
## Review boundary

- The 18-tip tree is a deterministic trim used by the integrated analysis, not a
  newly inferred phylogeny.
- The published main tree and 200-tree sample are independent sensitivity surfaces.
- The older processed tree is not treated as an independent phylogenetic hypothesis.
- Tree uncertainty must accompany path coefficients; a single-tree result is not
  the complete inferential result.
"""
            ),
        ],
    )


def data_tables_notebook():
    return make_notebook(
        "Analysis tables, panel boundaries, and provenance",
        "study_panels_v1",
        "TE34 n = 34; Cell21 n = 21; path18 n = 18",
        [
            md(
                """
## Three panels, three denominators

TE descriptions use every vetted genomic resource (TE34), cell descriptions use
every linked-cell species (Cell21), and only the exact intersection enters the
phylogenetic path analysis (path18). Species labels link biological summaries; they
do not imply that the genome and microscopy material came from the same specimen.
"""
            ),
            code(
                """
panel_paths = {
    "TE34": "path_analysis/data/derived/panels/study_te_resource_panel34_v1.csv",
    "Cell21": "path_analysis/data/derived/panels/study_cell_linked_panel21_v1.csv",
    "path18": "path_analysis/data/derived/panels/study_integrated_path_panel18_v1.csv",
}
panels = {name: read_csv(path) for name, path in panel_paths.items()}
expected = {"TE34": 34, "Cell21": 21, "path18": 18}
for name, frame in panels.items():
    assert len(frame) == expected[name], (name, len(frame))
    assert frame.species.is_unique
summary = pd.DataFrame([
    {"panel": name, "n_species": len(frame), "purpose": {
        "TE34": "vetted genomic-resource descriptive analyses",
        "Cell21": "linked-cell morphology and relative nuclear-IOD descriptions",
        "path18": "exact TE/cell/tree intersection for exploratory path analysis",
    }[name]}
    for name, frame in panels.items()
])
display(summary)
display(panels["TE34"])
display(panels["Cell21"])
display(panels["path18"])
"""
            ),
            md("## Inclusion logic and the non-interchangeable sample streams"),
            code(
                """
membership = (
    pd.concat([
        frame[["species"]].assign(**{name: True})
        for name, frame in panels.items()
    ], ignore_index=True)
    .groupby("species", as_index=False).max()
    .fillna(False)
    .sort_values("species")
)
assert "planiceps" not in set(membership.species)
fuscus = panels["TE34"].query("species == 'fuscus'")
assert len(fuscus) == 1
display(membership)
display(fuscus[["species", "te_sra_accession", "te_assembly_accession", "decision_basis"]])
"""
            ),
            md(
                """
### The two *fuscus* genomic resources are separate

The current comparative TE/LTR row uses `SRX20497025 / GCA_032353935.1`,
publicly labeled *planiceps* and accession-specifically reidentified as *fuscus*
by Alex Pyron. The reliable 2025 chromosome-level *fuscus* assembly
`GCA_050004315.1` is a distinct validation resource. It was not used to generate
the current dnaPipeTE, RepeatMasker, PCA, or LTR results and must not be described
as though it was.
"""
            ),
            code(
                """
source_manifest = read_csv("path_analysis/data/templates/source_manifest.csv")
fuscus_resource_rows = source_manifest[
    source_manifest.astype(str).apply(
        lambda row: row.str.contains("GCA_032353935.1|GCA_050004315.1", regex=True).any(),
        axis=1,
    )
]
display(fuscus_resource_rows)
"""
            ),
            md(
                """
## Accession and retry provenance

One TE34 row is one species-level active genomic resource. A dnaPipeTE retry is
not another species or biological replicate. Valid same-species completed runs are
normalized separately and then equal-weight averaged within species before the
species enters descriptive statistics.
"""
            ),
            code(
                """
mass = read_csv("results/data/corrected/te34_replicate_averaged/dnapipete_mass_accounting_te34_replicate_averaged_v1.csv")
assert len(mass) == 34 and mass.species.is_unique
display(mass[[
    "species", "te_sra_accession", "te_assembly_accession", "n_dnapipete_runs",
    "dnapipete_source_labels", "run_aggregation",
]].sort_values(["n_dnapipete_runs", "species"], ascending=[False, True]))
"""
            ),
            md("## Frozen artifact registry and release gates"),
            code(
                """
registry = read_csv("path_analysis/data/derived/source_file_registry.csv")
gates = read_csv("results/data/corrected/path_analysis/publication_release_gate_matrix_analysis18_v1.csv")
current_registry = registry[
    registry.source_id.ne("repo_ectopic_recombination_filtered_3000bp_5plusdomains")
].copy()
display(current_registry[[c for c in ["source_id", "local_path", "sha256", "exists", "origin"] if c in current_registry.columns]])
display(gates[["component", "status", "evidence", "permitted_language"]])
show_figure("results/figures/presentation_notebooks/data_tables/final18_evidence_completeness_analysis18_v1.png")
show_figure("results/figures/presentation_notebooks/data_tables/analysis18_pipeline_evidence_streams_v1.png")
"""
            ),
            source_cell(
                [
                    ("panel contract", "path_analysis/STUDY_SPECIES_PANELS.md"),
                    ("panel builder", "scripts/processing/build_study_species_panels.py"),
                    ("retry-averaged TE34 builder", "scripts/processing/build_te34_replicate_averaged_bundle.py"),
                    ("source manifest", "path_analysis/data/templates/source_manifest.csv"),
                ]
            ),
            md(
                """
## Review boundary

Do not shrink TE34 or Cell21 to the path18 overlap for descriptive figures. Do
not combine retry outputs as extra species. Do not restore *planiceps* as a row;
the current expert-reidentified *fuscus* comparative resource is accession-explicit
in the panel, and the separate chromosome-level validation assembly is not silently
substituted for it.
"""
            ),
        ],
    )


def repeat_notebook():
    figures = [
        "assembly_quality_te34_v1.png",
        "dnapipete_quality_te34_v1.png",
        "te_mass_composition_te34_v1.png",
        "te_mean_order_te34_v1.png",
        "te_mean_superfamily_te34_v1.png",
        "te_diversity_te34_v1.png",
        "te_shannon_cross_taxon_context_te34_v1.png",
        "te_shannon_genome_size_scatter_te18_v1.png",
        "te_pca_species_te34_v1.png",
        "te_pca_scree_te34_v1.png",
        "te_pca_elbow_te34_v1.png",
        "te_pca_silhouette_te34_v1.png",
        "te_pca_clusters_te34_v1.png",
        "te_pca_clusters_phylogeny_te34_v1.png",
    ]
    data = [
        "assembly_quality_te34_v1.csv",
        "dnapipete_quality_te34_v1.csv",
        "te_superfamily_prevalence_te34_v1.csv",
        "te_pca_superfamily_feature_matrix_te34_v1.csv",
        "te_pca_superfamily_loadings_te34_v1.csv",
        "te_pca_variance_te34_v1.csv",
        "te_pca_cluster_metrics_te34_v1.csv",
        "te_pca_scores_clusters_te34_v1.csv",
        "wang_ji_2025_shannon_context_v1.csv",
        "wang_ji_2025_shannon_context_v1.provenance.json",
        "te_shannon_cross_taxon_context_te34_v1.csv",
        "te_shannon_cross_taxon_group_counts_v1.csv",
        "te_top10_superfamily_contributions_te34_v1.csv",
        "te_shannon_genome_size_scatter_te18_v1.csv",
        "te_shannon_genome_size_match_audit_te34_v1.csv",
        "te_shannon_genome_size_stats_te18_v1.txt",
        "audited_historical_style_figures_v1.manifest.json",
    ]
    return make_notebook(
        "Repeat analysis: audited TE34 figures in the historical R style",
        "te34_replicate_averaged_v1",
        "TE34 genomic-resource panel (n = 34); not filtered by microscopy availability",
        [
            md(
                """
## What is being restored

These are regenerated figures from audited current data using the preserved
historical R code and ggplot style. They are not copied screenshots and they are
not extracted historical PDF/PNG panels. The plotting grammar is historical; the
species/resource contract and numerical inputs are current, cleaned, and versioned.

The current retry policy matters: `SRX19952890` and `SRX19952890R2` are two
completed dnaPipeTE runs for the same *D. orestes* resource. Each run is normalized
first, followed by equal-weight within-species averaging. They therefore contribute
one species row, not two observations.
"""
            ),
            code(
                f"""
expected_figures = [REVIEW_FIGURE_DIR / name for name in {figures!r}]
expected_data = [REVIEW_DATA_DIR / name for name in {data!r}]
require_paths(expected_figures + expected_data, build_command=FIGURE_BUILD_COMMAND)
manifest = json.loads((REVIEW_DATA_DIR / "audited_historical_style_figures_v1.manifest.json").read_text())
display(pd.json_normalize(manifest, sep="."))
"""
            ),
            md("## Exact figure-generation source and frozen inputs"),
            source_cell(
                [
                    ("current R review runner", "scripts/processing/build_audited_historical_style_figures.R"),
                    ("preserved historical composition adapter", "scripts/processing/restore_talk_iv_composition_from_historical_r.R"),
                    ("preserved historical PCA adapter", "scripts/processing/restore_talk_iv_pca_from_historical_r.R"),
                    ("TE34 panel", "path_analysis/data/derived/panels/study_te_resource_panel34_v1.csv"),
                    ("path18 exact TE/genome overlap", "path_analysis/data/derived/panels/study_integrated_path_panel18_v1.csv"),
                    ("fuscus-anchored genome-size estimates", "path_analysis/data/external/derived/cellprofiler_final_species_results.csv"),
                    ("genome-estimate calibration and state table", "path_analysis/data/external/derived/cellprofiler_genome_state_summary.csv"),
                    ("retry-averaged mass audit", "results/data/corrected/te34_replicate_averaged/dnapipete_mass_accounting_te34_replicate_averaged_v1.csv"),
                    ("retry-averaged diversity", "results/data/corrected/te34_replicate_averaged/te_diversity_mass_sensitivity_te34_replicate_averaged_v1.csv"),
                    ("run-level order source", "results/data/dnaPipeTE_order_breakdown.csv"),
                    ("run-level superfamily source", "results/data/dnaPipeTE_superfamily_breakdown.csv"),
                    ("audited retry policy and output hashes", "results/data/corrected/te34_replicate_averaged/te34_replicate_averaged_v1.manifest.json"),
                ],
                show_runner="scripts/processing/build_audited_historical_style_figures.R",
            ),
            md(
                """
## 1. Assembly and dnaPipeTE input quality

These panels report assembly contiguity and dnaPipeTE input/component diagnostics.
They are not a substitute for raw-read FastQC/MultiQC base-quality profiles or a
BUSCO completeness analysis; no such missing metric is implied by the word
"quality" here.
"""
            ),
            code(
                """
assembly_qc = read_csv("results/data/research_review/assembly_quality_te34_v1.csv")
run_qc = read_csv("results/data/research_review/dnapipete_quality_te34_v1.csv")
assert assembly_qc.species.nunique() == 34
assert run_qc.species.nunique() == 34
display(assembly_qc)
display(run_qc)
show_figure("results/figures/research_review/assembly_quality_te34_v1.png")
show_figure("results/figures/research_review/dnapipete_quality_te34_v1.png")
"""
            ),
            md("## 2. Mass accounting and mean TE composition"),
            code(
                """
mass = read_csv("results/data/corrected/te34_replicate_averaged/dnapipete_mass_accounting_te34_replicate_averaged_v1.csv")
assert len(mass) == 34 and mass.species.is_unique
orestes = mass.query("species == 'orestes'")
assert len(orestes) == 1
assert "SRX19952890" in orestes.dnapipete_source_labels.iloc[0]
assert "SRX19952890R2" in orestes.dnapipete_source_labels.iloc[0]
assert orestes.n_dnapipete_runs.iloc[0] == 2
display(orestes[["species", "n_dnapipete_runs", "dnapipete_source_labels", "run_aggregation"]])
for figure in [
    "results/figures/research_review/te_mass_composition_te34_v1.png",
    "results/figures/research_review/te_mean_order_te34_v1.png",
    "results/figures/research_review/te_mean_superfamily_te34_v1.png",
]:
    show_figure(figure)
"""
            ),
            md(
                """
## 3. Restored cross-taxon Shannon context and genome-size scatter

This boxplot restores the literal historical `Master_shannon_diversity_vis.R`
grammar using the current audited TE34 top-10-superfamily Shannon values and the
cleaned 87-row published vertebrate composite. The metric is natural-log Shannon
entropy after renormalizing each species' ten most abundant superfamilies to one.

This is **descriptive cross-study context**, not a homogeneous comparative test:
repeat libraries, sequencing depth, genome assemblies versus skims, annotation
protocols, and group sample sizes differ. Source context: [Wang et al. 2025,
*Communications Biology*, DOI 10.1038/s42003-025-08127-3](https://doi.org/10.1038/s42003-025-08127-3)
and [Dryad DOI 10.5061/dryad.zpc866tkv](https://doi.org/10.5061/dryad.zpc866tkv).
The local 87-row historical composite is not an exact export of Wang et al.'s
84-species 2025 analysis—it contains four caecilians and ten salamanders—so it
restores prior project context only.

The historical genome-size scatter is now restored from real, variable values:
the exact 18-species intersection between audited TE34 Shannon and the current
24-species genome-estimate table. The original `GS = 16` placeholder is gone.
Current Desmognathus x values are **fuscus-anchored image-IOD genome-size
estimates**, converted with `1 pg = 0.978 Gb` and shown with their conditional
95% bootstrap intervals. They are estimates, **not direct C-values**, and the
16.36-pg *D. fuscus* anchor does not by itself validate between-slide absolute
calibration. The ten historical salamander points remain descriptive context;
no pooled 28-point inferential test is reported because assay and TE-annotation
workflows differ. Pearson, Spearman, OLS, and Brownian PGLS summaries are
reported only for the 18 current Desmognathus rows and remain exploratory.
Those tests use species point estimates and the single frozen processed time tree;
they do not propagate genome-estimate, Shannon-estimate, or tree uncertainty.
"""
            ),
            code(
                """
cross_taxon = read_csv("results/data/research_review/te_shannon_cross_taxon_context_te34_v1.csv")
group_counts = read_csv("results/data/research_review/te_shannon_cross_taxon_group_counts_v1.csv")
top10 = read_csv("results/data/research_review/te_top10_superfamily_contributions_te34_v1.csv")
wang_context = read_csv("results/data/research_review/wang_ji_2025_shannon_context_v1.csv")
wang_provenance = json.loads((REVIEW_DATA_DIR / "wang_ji_2025_shannon_context_v1.provenance.json").read_text())
assert len(wang_context) == 87
assert len(cross_taxon) == 121
assert cross_taxon.query("group == 'Desmognathus'").species.nunique() == 34
assert len(top10) == 340 and top10.species.nunique() == 34
display(group_counts)
display(pd.json_normalize(wang_provenance, sep="."))
display(top10.groupby("species", as_index=False).agg(
    top10_total=("top10_renormalized_proportion", "sum"),
    shannon_top10=("shannon_component", "sum"),
))
show_figure("results/figures/research_review/te_shannon_cross_taxon_context_te34_v1.png")

genome_scatter = read_csv("results/data/research_review/te_shannon_genome_size_scatter_te18_v1.csv")
genome_match_audit = read_csv("results/data/research_review/te_shannon_genome_size_match_audit_te34_v1.csv")
genome_stats = (REVIEW_DATA_DIR / "te_shannon_genome_size_stats_te18_v1.txt").read_text()
desmo_scatter = genome_scatter.query("group == 'Desmognathus'").copy()
historical_salamanders = genome_scatter.query("group == 'Salamanders'").copy()
assert len(genome_scatter) == 28
assert len(desmo_scatter) == 18 and desmo_scatter.canonical_species.is_unique
assert len(historical_salamanders) == 10
assert desmo_scatter.genome_size_gb.nunique() > 10
assert set(desmo_scatter.is_direct_c_value.astype(str).str.lower()) == {"false"}
assert len(genome_match_audit) == 37
assert genome_match_audit.included_in_scatter.astype(str).str.lower().eq("true").sum() == 18
display(desmo_scatter[[
    "species", "genome_size_gb", "genome_size_ci_low_gb",
    "genome_size_ci_high_gb", "shannon_top10_superfamilies",
    "genome_support_tier", "genome_support_warnings",
]])
display(genome_match_audit)
print(genome_stats)
show_figure("results/figures/research_review/te_shannon_genome_size_scatter_te18_v1.png")
"""
            ),
            md("## 4. Within-TE34 level and denominator sensitivity"),
            code(
                """
diversity = read_csv("results/data/corrected/te34_replicate_averaged/te_diversity_mass_sensitivity_te34_replicate_averaged_v1.csv")
assert diversity.species.nunique() == 34
display(diversity.sort_values(["te_level", "composition_mode", "species"]))
show_figure("results/figures/research_review/te_diversity_te34_v1.png")
"""
            ),
            md(
                """
## 5. Superfamily prevalence and the PCA feature contract

The original 35-species review identified three sparse superfamilies: Chapaev
in 14/35 species (40%), Dada in 10/35 (28.57%), and Ginger in 1/35 (2.86%).
The documented decision was to **retain Chapaev and Dada and exclude Ginger**.

The audited panel now contains 34 genomic resources, so those historical counts
must not be silently relabeled as TE34 results. After resource exclusions and
within-species retry averaging, the corresponding current counts are 12/34,
9/34, and 1/34. The decision itself is unchanged. The table and exact wide
matrix below are generated by the same R runner that performs the PCA, making
the filter executable and reviewable rather than a prose-only convention.
"""
            ),
            code(
                """
prevalence = read_csv("results/data/research_review/te_superfamily_prevalence_te34_v1.csv")
pca_feature_matrix = read_csv("results/data/research_review/te_pca_superfamily_feature_matrix_te34_v1.csv")
pca_loadings = read_csv("results/data/research_review/te_pca_superfamily_loadings_te34_v1.csv")
pca_variance = read_csv("results/data/research_review/te_pca_variance_te34_v1.csv")

policy_features = prevalence.query("superfamily in ['Chapaev', 'Dada', 'Ginger']").copy()
observed_presence = dict(zip(policy_features.superfamily, policy_features.species_present))
observed_decisions = dict(zip(policy_features.superfamily, policy_features.pca_decision))
assert observed_presence == {"Chapaev": 12, "Dada": 9, "Ginger": 1}
assert observed_decisions == {"Chapaev": "retain", "Dada": "retain", "Ginger": "exclude"}
assert policy_features.set_index("superfamily").loc["Ginger", "decision_reason"] == "user_declared_sparse_superfamily"
assert len(pca_feature_matrix) == 34 and pca_feature_matrix.species.is_unique
assert {"Chapaev", "Dada"}.issubset(pca_feature_matrix.columns)
assert "Ginger" not in pca_feature_matrix.columns
assert set(pca_loadings.superfamily) == set(pca_feature_matrix.columns) - {"species"}

display(policy_features[[
    "superfamily", "species_present", "panel_species", "prevalence_percent",
    "pca_decision", "decision_reason",
]])
display(prevalence.sort_values(["pca_decision", "species_present", "superfamily"], ascending=[True, True, True]))
display(Markdown(
    f"**Exact PCA input:** {len(pca_feature_matrix)} species × "
    f"{len(pca_feature_matrix.columns) - 1} retained superfamilies."
))
"""
            ),
            md(
                """
## 6. PCA variance and species scores

This section recreates the historical PCA visual analysis on the audited TE34
composition table. It retains the original R plotting grammar while making the
documented current-panel, deterministic-seed, and silhouette-selected-cluster
adaptations. The current-data variance and coordinates are loaded below;
historical slide percentages are not substituted.
"""
            ),
            code(
                """
scores = read_csv("results/data/research_review/te_pca_scores_clusters_te34_v1.csv")
assert len(scores) == 34 and scores.species.is_unique
display(pca_variance)
display(pca_loadings.sort_values("PC1", key=lambda values: values.abs(), ascending=False).head(10))
display(scores)
show_figure("results/figures/research_review/te_pca_species_te34_v1.png")
show_figure("results/figures/research_review/te_pca_scree_te34_v1.png")
"""
            ),
            md("## 7. Cluster-number diagnostics and the historical k-means view"),
            code(
                """
cluster_metrics = read_csv("results/data/research_review/te_pca_cluster_metrics_te34_v1.csv")
display(cluster_metrics)
cluster_sizes = scores.groupby("cluster", as_index=False).agg(
    n_species=("species", "size"),
    species=("species", lambda values: ", ".join(sorted(values))),
)
display(cluster_sizes)
if cluster_sizes.n_species.min() == 1:
    singleton = cluster_sizes.loc[cluster_sizes.n_species == 1, "species"].iloc[0]
    display(Markdown(
        f"**Current-data result:** maximum silhouette selects a singleton `{singleton}` "
        "outlier versus the remaining species. Treat this as outlier structure, not "
        "evidence for two broad biological TE-composition groups."
    ))
show_figure("results/figures/research_review/te_pca_elbow_te34_v1.png")
show_figure("results/figures/research_review/te_pca_silhouette_te34_v1.png")
show_figure("results/figures/research_review/te_pca_clusters_te34_v1.png")
show_figure("results/figures/research_review/te_pca_clusters_phylogeny_te34_v1.png")
"""
            ),
            md(
                """
## Interpretation and review boundary

- These are current TE34 results regenerated with historical R code/style.
- The PCA feature contract is explicit: Chapaev and Dada are retained; Ginger
  is excluded. Historical 35-species prevalence and current TE34 prevalence are
  reported separately.
- PCA clusters are descriptive. Elbow and silhouette diagnostics must be shown
  beside any selected `k`; the selected partition is not a discovered taxonomy.
- On the current raw-proportion superfamily PCA, maximum silhouette separates
  *D. orestes* as a singleton from the other 33 species. That is an outlier
  diagnostic, not support for two broadly replicated clades. Scaled raw-proportion
  PCA also gives rare categories equal variance, so compositional/phylogenetic PCA
  remains necessary sensitivity analysis.
- Composition is closed and unresolved mass is explicit. Corrected CLR and
  phylogenetic-PCA results remain sensitivity analyses, not replacements for this
  requested historical visual layer.
- Historical presentation rasters are visual references only and are deliberately
  not displayed as current results in this notebook.
"""
            ),
        ],
    )


def ltr_notebook():
    return make_notebook(
        "LTR terminal:internal deletion-footprint analysis",
        "ltr30_v1",
        "All vetted genomic resources with eligible LTR elements (LTR30), not path18",
        [
            md(
                """
## Estimand and denominator

The element-level value is terminal-LTR read depth divided by internal-region read
depth for structurally annotated LTR elements. It is a deletion-footprint and
mapping proxy—not a direct solo-LTR count, an orthologous deletion event, or an
ectopic-recombination rate. The descriptive panel includes every vetted species
with eligible data; it is not reduced to the TE/cell path18 intersection.

The primary ratio includes every expected terminal and internal position, retaining
explicit zero-depth positions. The corrected branch contains **1,086 usable elements
across 30 species**. Two of the original 1,088 selected elements are excluded because
their depth files end in truncated partial lines; both exclusions are preserved in
the source audit. The old 1,088-element nonzero-only table is retained solely for
numerical comparison and is not a primary analysis input.

The violin is regenerated from audited current data using historical R code/style.
No historical PDF/PNG raster is presented as a current result.
"""
            ),
            code(
                """
    expected = [
        REVIEW_FIGURE_DIR / "ltr_terminal_internal_artifact_screened_ltr30_v1.png",
        REVIEW_FIGURE_DIR / "ltr_terminal_internal_violin_ltr30_v1.png",
        REVIEW_FIGURE_DIR / "ltr_terminal_internal_log_ltr30_v1.png",
        REVIEW_FIGURE_DIR / "ectopic_ltr30_coverage_qc_v1.png",
        REVIEW_DATA_DIR / "ltr_mapping_artifact_screen_ltr30_v1.csv",
    REVIEW_DATA_DIR / "ltr_terminal_internal_stats_v1.txt",
    REVIEW_DATA_DIR / "audited_historical_style_figures_v1.manifest.json",
    resolve_artifact("results/data/corrected/ectopic_ltr30/ectopic_element_metrics_ltr30_v1.csv"),
    resolve_artifact("results/data/corrected/ectopic_ltr30/ectopic_species_robustness_ltr30_v1.csv"),
    resolve_artifact("results/data/corrected/ectopic_ltr30/ectopic_excluded_elements_ltr30_v1.csv"),
    resolve_artifact("results/data/corrected/ectopic_ltr30/ectopic_resource_coverage_ltr30_v1.csv"),
    resolve_artifact("results/data/corrected/ectopic_ltr30/ectopic_ltr30_v1.manifest.json"),
]
require_paths(expected, build_command=FIGURE_BUILD_COMMAND)
"""
            ),
            md("## Exact R source and element-level input"),
            source_cell(
                [
                    ("current R review runner", "scripts/processing/build_audited_historical_style_figures.R"),
                    ("corrected zero-aware LTR30 builder", "scripts/processing/build_corrected_ectopic_ltr30.py"),
                    ("preserved historical LTR adapter", "scripts/processing/restore_talk_iv_ectopic_from_historical_r.R"),
                    ("corrected all-position element input", "results/data/corrected/ectopic_ltr30/ectopic_element_metrics_ltr30_v1.csv"),
                    ("species robustness", "results/data/corrected/ectopic_ltr30/ectopic_species_robustness_ltr30_v1.csv"),
                    ("source-corrupt exclusions", "results/data/corrected/ectopic_ltr30/ectopic_excluded_elements_ltr30_v1.csv"),
                    ("resource coverage audit", "results/data/corrected/ectopic_ltr30/ectopic_resource_coverage_ltr30_v1.csv"),
                    ("corrected LTR30 manifest", "results/data/corrected/ectopic_ltr30/ectopic_ltr30_v1.manifest.json"),
                    ("vetted genomic panel", "path_analysis/data/derived/panels/study_te_resource_panel34_v1.csv"),
                    ("statistics", "results/data/research_review/ltr_terminal_internal_stats_v1.txt"),
                ],
                show_runner="scripts/processing/build_audited_historical_style_figures.R",
            ),
            md("## Element coverage and species support"),
            code(
                """
elements = read_csv("results/data/corrected/ectopic_ltr30/ectopic_element_metrics_ltr30_v1.csv")
robustness = read_csv("results/data/corrected/ectopic_ltr30/ectopic_species_robustness_ltr30_v1.csv")
excluded = read_csv("results/data/corrected/ectopic_ltr30/ectopic_excluded_elements_ltr30_v1.csv")
coverage = read_csv("results/data/corrected/ectopic_ltr30/ectopic_resource_coverage_ltr30_v1.csv")
ltr_manifest = json.loads(resolve_artifact(
    "results/data/corrected/ectopic_ltr30/ectopic_ltr30_v1.manifest.json"
).read_text())
vetted = read_csv("path_analysis/data/derived/panels/study_te_resource_panel34_v1.csv")
assert len(elements) == 1086 and elements.species.nunique() == 30
assert set(elements.species).issubset(set(vetted.species))
assert elements.source_depth_status.eq("usable").all()
assert elements.ratio_terminal_internal_all_positions.gt(0).all()
assert np.isfinite(elements.ratio_terminal_internal_all_positions).all()
assert len(excluded) == 2 and excluded.usable.eq(False).all()
assert ltr_manifest["n_selected_elements"] == 1088
assert ltr_manifest["n_usable_elements"] == 1086
assert ltr_manifest["n_source_corrupt_exclusions"] == 2
assert ltr_manifest["n_species_with_usable_elements"] == 30
support = elements.groupby("species", as_index=False).agg(
    n_elements=("element_id", "size"),
    median_terminal_internal=("ratio_terminal_internal_all_positions", "median"),
    mean_terminal_internal=("ratio_terminal_internal_all_positions", "mean"),
    median_reported_position_fraction=("reported_position_fraction", "median"),
)
assert support.species.nunique() == 30, support.species.nunique()
display(support.sort_values("median_terminal_internal"))
display(excluded[[
    "species", "element_id", "depth_file", "reported_position_fraction",
    "partial_final_line", "exclusion_reason", "input_modified",
]])
display(coverage.sort_values(["resource_status", "species"]))
display(robustness.sort_values(["analysis_branch", "species"]))
display(elements[[
    "species", "element_id", "Superfamily", "Complete", "domain_count",
    "mean_depth_terminal_all_positions", "mean_depth_internal_all_positions",
    "ratio_terminal_internal_all_positions", "reported_position_fraction",
    "terminal_positive_coverage_fraction", "internal_positive_coverage_fraction",
]].head(40))
show_figure("results/figures/corrected/ectopic_ltr30/ectopic_ltr30_coverage_qc_v1.png")
"""
            ),
            md(
                """
## Artifact-screened historical-style distribution

Jake's historical analysis intentionally prevented non-biological mapping
pileups from dictating the presentation scale. The corrected implementation keeps
the complete 1,086-element table unchanged and exports a separate row-level screen.
It flags one element only when all three conditions hold: it is above the
within-species 1.5-IQR upper fence, its left and right LTR depths differ by more
than 64-fold, and it contributes over half of its species' summed ratio. This
retains 1,085 elements and all 30 species. It is not blanket IQR pruning.
"""
            ),
            code(
                """
artifact_screen = read_csv("results/data/research_review/ltr_mapping_artifact_screen_ltr30_v1.csv")
flagged_mapping_artifacts = artifact_screen.query("artifact_screen_keep == False")
assert len(artifact_screen) == 1086
assert artifact_screen.artifact_screen_keep.sum() == 1085
assert artifact_screen.query("artifact_screen_keep == True").species.nunique() == 30
assert flagged_mapping_artifacts.element_id.tolist() == ["JAUEJH010597481.1_De_5169_11831"]
display(flagged_mapping_artifacts[[
    "species", "element_id", "ratio_terminal_internal_all_positions",
    "left_right_terminal_log2_imbalance", "species_iqr_upper",
    "fraction_of_species_ratio_sum", "artifact_screen_reason",
]])
show_figure("results/figures/research_review/ltr_terminal_internal_artifact_screened_ltr30_v1.png")
"""
            ),
            md(
                """
## Robust distribution and correctly labeled tests

The unfiltered audit is retained below, followed by the log view. The report prints
raw, log10, and Kruskal-Wallis results for both the complete and artifact-screened
branches. These tests are not interchangeable. The log view exposes skew and
leverage that a raw violin can hide; species-level sample counts remain visible.
The exact report labels are `raw one-way ANOVA`, `log10 one-way ANOVA`, and
`Kruskal-Wallis raw`, with parallel artifact-screened labels.
"""
            ),
            code(
                """
show_figure("results/figures/research_review/ltr_terminal_internal_violin_ltr30_v1.png")
show_figure("results/figures/research_review/ltr_terminal_internal_log_ltr30_v1.png")
show_text("results/data/research_review/ltr_terminal_internal_stats_v1.txt")
"""
            ),
            md(
                """
## Interpretation and review boundary

- State the number of species and eligible elements with every result.
- Use element-level distributions plus species support; do not infer equal
  precision from unequal element counts.
- Element rows are nested within species. The three omnibus tests are descriptive
  diagnostics, not 1,086 independent biological replicates or a substitute for a
  hierarchical species/family model.
- The two truncated source files are excluded without repair or imputation; the
  exact file hashes and exclusion reasons remain in the audit table.
- The one mapping-artifact candidate is excluded only from the screened
  presentation/sensitivity branch. Its complete corrected row remains in the
  unfiltered audit table and raw/log companion figures.
- The raw ANOVA is the historically used descriptive test. The log-scale ANOVA
  and rank-based Kruskal-Wallis result are robust companions, not relabeled raw
  results.
- A terminal:internal depth ratio is mechanistically suggestive, but it does not
  by itself measure the rate of ectopic recombination.
"""
            ),
        ],
    )


def cell_notebook():
    return make_notebook(
        "Cell modeling, mask review, and linked morphometry",
        "cell21_v1",
        "Cell21 linked-cell descriptive panel (n = 21); path18 only where explicitly labeled",
        [
            md(
                """
## Measurement design

Cell and nucleus areas come from linked cell–nucleus pairs after review of
overlaps, broken cells, blur, and irregular objects. The frozen primary morphology
summary is a quality-screened upper-tail estimand: the median of the 50
highest-composite-ranked eligible linked erythrocytes per species and their
corresponding nuclei. It is not a literal sort of the 50 largest masks, and it must
not be called the typical erythrocyte.

Relative nuclear-IOD is an image-density phenotype. It is **not absolute genome
size** and is **not calibrated to pg**. Image quality, staining, illumination, and
nuclear area all affect it.
"""
            ),
            code(
                """
panel = read_csv("path_analysis/data/derived/panels/study_cell_linked_panel21_v1.csv")
traits = read_csv("results/data/corrected/cell21/cell_linked_traits_cell21_v1.csv")
assert len(panel) == 21 and panel.species.is_unique
assert len(traits) == 21 and traits.species.is_unique
assert set(panel.species) == set(traits.species)
display(traits.sort_values("cell_area_um2"))
show_figure("results/figures/corrected/cell21/cell21_cell_nucleus_relative_iod_v1.png")
"""
            ),
            md("## Full mask and model-review pages"),
            code(
                """
CELL_ROOT = ROOT.parent / "cellprofiler_test"
viewer_rows = [
    ("Cell–nucleus linkage", "output/runs/mixed_cellpose_yolo_full_dataset_v1_bgclean/linkage/index.html", "linked masks and assignment QC"),
    ("Tile annotation gallery", "output/tile_annotation_bundle_v1/index.html", "training-domain annotations"),
    ("Raw / overlay / mask review", "output/tile_bootstrap_review_v1/index.html", "segmentation visual QC"),
    ("Reviewed species statistics", "output/runs/mixed_cellpose_yolo_full_dataset_v1_bgclean/linked_species_stats_reviewed/index.html", "historical diagnostic labels; inspect, do not quote pg"),
    ("Historical cell phylogenetic report", "output/runs/mixed_cellpose_yolo_full_dataset_v1_bgclean/phylogenetic_analysis/index.html", "historical diagnostic report"),
]
viewers = pd.DataFrame(viewer_rows, columns=["viewer", "relative_path", "scope"])
viewers["path"] = viewers.relative_path.map(lambda value: CELL_ROOT / value)
viewers["exists"] = viewers.path.map(Path.exists)
display(viewers[["viewer", "relative_path", "scope", "exists"]])
available_viewers = viewers[viewers.exists].copy()

def start_optional_local_viewers():
    # Start local-only HTML viewers interactively; never persist their port.
    if available_viewers.empty:
        display(Markdown(
            "The optional sibling `cellprofiler_test` HTML workspace is not present. "
            "The portable static figures below remain available."
        ))
        return None
    import functools
    import threading
    from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
    global _cell_review_server, _cell_review_thread
    if "_cell_review_server" not in globals():
        handler = functools.partial(SimpleHTTPRequestHandler, directory=str(ROOT.parent))
        _cell_review_server = ThreadingHTTPServer(("127.0.0.1", 0), handler)
        _cell_review_thread = threading.Thread(target=_cell_review_server.serve_forever, daemon=True)
        _cell_review_thread.start()
    base = f"http://127.0.0.1:{_cell_review_server.server_address[1]}/cellprofiler_test"
    live = available_viewers.copy()
    live["url"] = live.relative_path.map(lambda value: f"{base}/{value}")
    display(HTML("".join(
        f'<p><a href="{row.url}" target="_blank"><strong>{row.viewer}</strong></a> — {row.scope}</p>'
        for row in live.itertuples()
    )))
    display(IFrame(src=live.iloc[0].url, width="100%", height=700))
    return _cell_review_server

if available_viewers.empty:
    display(Markdown(
        "**Optional local viewers unavailable.** This does not affect the portable "
        "static audit below."
    ))
else:
    display(Markdown(
        f"**Optional local viewers available: {len(available_viewers)}.** To open them "
        "during an interactive session, run `start_optional_local_viewers()`. The "
        "executed notebook intentionally saves no ephemeral localhost URL."
    ))
"""
            ),
            md("## Static mask examples and segmentation/support diagnostics"),
            code(
                """
show_figure("results/figures/presentation_notebooks/cell_models/cell_mask_examples_analysis18_v1.png")
show_figure("results/figures/corrected/microscopy/microscopy_segmentation_validation_analysis18_v1.png")
show_figure("results/figures/corrected/microscopy/microscopy_support_and_review_analysis18_v1.png")
validation = read_csv("results/data/corrected/microscopy/microscopy_segmentation_validation_analysis18_v1.csv")
display(validation)
"""
            ),
            md("## Largest-50 estimand and selection-depth sensitivity"),
            code(
                """
estimands = read_csv("results/data/corrected/microscopy/microscopy_estimand_sensitivity_analysis18_v1.csv")
stability = read_csv("results/data/corrected/microscopy/microscopy_estimator_rank_stability_analysis18_v1.csv")
display(estimands)
display(stability)
show_figure("results/figures/corrected/microscopy/microscopy_estimand_sensitivity_analysis18_v1.png")
show_figure("results/figures/corrected/microscopy/microscopy_selection_depth_analysis18_v1.png")
show_figure("results/figures/corrected/microscopy/microscopy_estimator_rank_stability_analysis18_v1.png")
"""
            ),
            md(
                """
## Pairwise traits and phylogenetic tracks

The next two frozen figures use the path18 overlap because they require the exact
comparative tree and complete integrated traits. They do not redefine the Cell21
descriptive panel. Raw species points, phylogenetic sensitivity, and standardized
tree tracks are shown together.
"""
            ),
            code(
                """
show_figure("results/figures/presentation_notebooks/cell_models/cell_trait_pairwise_bm_pgls_analysis18_v1.png")
show_figure("results/figures/presentation_notebooks/cell_models/cell_trait_phylogeny_tracks_analysis18_v1.png")
"""
            ),
            md("## Relative nuclear-IOD quality and dependence checks"),
            code(
                """
for figure in [
    "results/figures/corrected/microscopy/microscopy_relative_iod_qc_sensitivity_analysis18_v1.png",
    "results/figures/corrected/microscopy/microscopy_iod_decomposition_analysis18_v1.png",
    "results/figures/corrected/microscopy/microscopy_image_iod_qc_analysis18_v1.png",
]:
    show_figure(figure)
"""
            ),
            source_cell(
                [
                    ("Cell21 frozen builder", "scripts/processing/build_cell21_descriptive_bundle.py"),
                    ("microscopy release audit", "scripts/processing/audit_microscopy_release.py"),
                    ("largest-cell mask review builder", "path_analysis/scripts/build_largest_cell_mask_review.py"),
                    ("CellProfiler provenance audit", "path_analysis/CELLPROFILER_PROVENANCE_AUDIT.md"),
                ]
            ),
            md(
                """
## Interpretation and review boundary

- Cell21 is the descriptive denominator; path18 is used only for explicitly
  comparative plots.
- The selected-50 statistic is a quality-screened upper-tail morphology estimand,
  not a literal largest-50 sort or a species-typical erythrocyte estimate.
- Existing segmentation validation is diagnostic and does not establish held-out
  generalization across all focal species.
- Relative nuclear-IOD cannot be presented as picograms, C-value, or an independent
  validation of genome size.
"""
            ),
        ],
    )


def cell_nucleus_genome_path_notebook():
    return make_notebook(
        "Cell–nucleus–genome path design and identifiability",
        "all_finalized_path_design_v1",
        "All finalized genome-size × nucleus × cell species (n = 24)",
        [
            md(
                """
## What this notebook establishes

This is the design and input audit for the study's central three-trait question.
Genome size is the conditional *D. fuscus*-anchored estimate from Notebook 06;
nucleus and cell areas are medians of the 50 manually vetted paired cells. The
panel contains all 24 finalized species with all three measured traits and uses
no phylogenetic trait filling. Reduced sample sizes and limited-overlap image
quality are retained as interpretation fields, never exclusion rules.

This notebook deliberately separates **identifiability** from model fitting. It
enumerates every labeled acyclic three-node DAG, identifies the Markov-equivalence
classes, and freezes measurement-bootstrap inputs. It does not claim that an
arrow direction has been causally identified.
"""
            ),
            code(
                """
import sys
sys.path.insert(0, str(ROOT / "path_analysis/scripts"))
import build_cell_nucleus_genome_path_notebook as design

design_dir = "results/data/research_review/cell_nucleus_genome_path"
traits = read_csv(f"{design_dir}/primary_traits_all_finalized.csv")
dags = read_csv(f"{design_dir}/all_three_trait_dags.csv")
bootstrap = read_csv(f"{design_dir}/measurement_bootstrap_traits_250.csv")
manifest_path = resolve_artifact(
    f"{design_dir}/cell_nucleus_genome_path_design_manifest.json"
)
manifest = json.loads(manifest_path.read_text())

assert len(traits) >= 3 and traits.species.is_unique
assert len(dags) == 25 and dags.equivalence_class.nunique() == 11
assert bootstrap.bootstrap_replicate.nunique() == 250
assert manifest["actual_phylogenetic_path_models_fitted"] is False
display(traits)
display(pd.DataFrame([manifest]).T.rename(columns={0: "value"}))
"""
            ),
            md("## Exhaustive DAG and equivalence-class audit"),
            code(
                """
class_summary = (
    dags.groupby(["equivalence_class", "basis_claim", "testable_by_dsep"], dropna=False)
    .agg(n_labeled_dags=("dag_id", "size"), min_edges=("n_edges", "min"))
    .reset_index()
    .sort_values(["min_edges", "equivalence_class"])
)
display(class_summary)
display(dags)

mechanism_rows = dags.loc[dags.dag_id.isin(design.USER_MECHANISMS.values())].copy()
mechanism_rows.insert(
    0,
    "proposed_mechanism",
    mechanism_rows.dag_id.map({value: key for key, value in design.USER_MECHANISMS.items()}),
)
display(mechanism_rows[
    ["proposed_mechanism", "dag_id", "equivalence_class", "basis_claim"]
])
assert mechanism_rows.equivalence_class.unique().tolist() == ["nucleus_bridge"]
assert mechanism_rows.basis_claim.unique().tolist() == [
    "genome_size _||_ cell_size | nucleus_size"
]
"""
            ),
            md(
                """
## Central identifiability result

The three proposed mechanisms—genome → nucleus → cell, cell → nucleus →
genome, and nucleus → both genome and cell—have the same skeleton, no collider,
and the same d-separation claim. They are observationally Markov-equivalent.
A cross-sectional covariance matrix, even after phylogenetic correction, can
test whether nucleus size behaves as the supported bridge but cannot orient
those three causal arrow patterns without external temporal, experimental, or
biological direction information.
"""
            ),
            md("## Frozen measurement-bootstrap input"),
            code(
                """
bootstrap_counts = bootstrap.groupby("bootstrap_replicate").agg(
    n_species=("species", "nunique"),
    genome_min=("genome_size_pg", "min"),
    genome_max=("genome_size_pg", "max"),
    nucleus_min=("nucleus_area_um2", "min"),
    nucleus_max=("nucleus_area_um2", "max"),
    cell_min=("cell_area_um2", "min"),
    cell_max=("cell_area_um2", "max"),
)
assert bootstrap_counts.n_species.eq(len(traits)).all()
display(bootstrap_counts.describe().T)
display(bootstrap.head(40))
"""
            ),
            md("## Reference-method audit"),
            code(
                """
show_text(
    "path_analysis/CELL_NUCLEUS_GENOME_PATH_REFERENCE_METHOD_AUDIT.md",
    language="markdown",
)
"""
            ),
            source_cell(
                [
                    ("design/input builder", "path_analysis/scripts/build_cell_nucleus_genome_path_notebook.py"),
                    ("reference-method audit", "path_analysis/CELL_NUCLEUS_GENOME_PATH_REFERENCE_METHOD_AUDIT.md"),
                    ("frozen genome-size builder", "path_analysis/scripts/build_frozen_genome_iod_notebook.py"),
                    ("existing integrated path scaffold", "path_analysis/scripts/path_model_scaffold.R"),
                ]
            ),
            md(
                """
## Interpretation and review boundary

- This notebook establishes the all-finalized-species input panel and candidate-DAG
  identifiability limit; it does not report a fitted causal winner.
- Six saturated DAGs have no d-separation claim and cannot be tested by the
  confirmatory-path procedure.
- The 250 frozen measurement-bootstrap replicates propagate observed image and
  selected-cell sampling uncertainty, but not uncertainty in the 16.36-pg anchor.
- Any fitted follow-up must report global fit, CICc/delta/weights, tree and
  leave-one-species-out sensitivity, and equivalence-class rather than arrow-level
  conclusions.
"""
            ),
        ],
    )


def integrated_path_notebook():
    return make_notebook(
        "Exploratory phylogenetic path analysis",
        "path18_v1",
        "Exact TE34 ∩ Cell21 ∩ tree panel (path18, n = 18)",
        [
            md(
                """
## Inferential boundary

The models compare a small set of biological DAGs on the exact 18-species overlap.
They use relative nuclear-IOD, not absolute genome size. Candidate models were
formalized after historical results existed, so this is exploratory phylogenetic
path analysis—not prospectively confirmatory causal identification.
"""
            ),
            code(
                """
panel = read_csv("path_analysis/data/derived/panels/study_integrated_path_panel18_v1.csv")
inputs = read_csv("results/data/corrected/path_analysis/corrected_path_input_sensitivity_analysis18_v1.csv")
specifications = read_csv("results/data/corrected/path_analysis/corrected_path_measurement_specifications_analysis18_v1.csv")
assert len(panel) == 18 and panel.species.is_unique
assert set(inputs.species) == set(panel.species)
assert "genome_size_pg" not in inputs.columns
display(specifications)
display(inputs.head(36))
"""
            ),
            md("## Release gates before coefficients"),
            code(
                """
gates = read_csv("results/data/corrected/path_analysis/publication_release_gate_matrix_analysis18_v1.csv")
display(gates[["component", "status", "evidence", "permitted_language"]])
show_figure("results/figures/corrected/path_analysis/publication_release_gate_matrix_analysis18_v1.png")
"""
            ),
            md("## Candidate DAGs, global fit, model weights, and edge estimates"),
            code(
                """
models = read_csv("results/data/corrected/path_analysis/corrected_path_anchor_model_comparison_analysis18_v1.csv")
edges = read_csv("results/data/corrected/path_analysis/corrected_path_anchor_edges_analysis18_v1.csv")
basis = read_csv("results/data/corrected/path_analysis/corrected_path_data_sensitivity_basis_sets_analysis18_v1.csv")
display(models)
display(edges)
display(basis[(basis.tree_id == "published_main") & (basis.morphology_estimator == "image_balanced_selected50") & (basis.iod_subset == "image_qc_pass")])
show_figure("results/figures/corrected/path_analysis/corrected_path_anchor_dag_analysis18_v1.png")
show_figure("results/figures/corrected/path_analysis/corrected_path_anchor_model_comparison_analysis18_v1.png")
"""
            ),
            md("## Measurement, tree, and leave-one-species-out sensitivity"),
            code(
                """
for table in [
    "results/data/corrected/path_analysis/corrected_path_model_stability_analysis18_v1.csv",
    "results/data/corrected/path_analysis/corrected_path_tree_edge_uncertainty_analysis18_v1.csv",
    "results/data/corrected/path_analysis/corrected_path_leave_one_out_rankings_analysis18_v1.csv",
]:
    display(read_csv(table))
for figure in [
    "results/figures/corrected/path_analysis/corrected_path_measurement_model_weights_analysis18_v1.png",
    "results/figures/corrected/path_analysis/corrected_path_tree_uncertainty_analysis18_v1.png",
    "results/figures/corrected/path_analysis/corrected_path_leave_one_out_influence_analysis18_v1.png",
]:
    show_figure(figure)
"""
            ),
            md("## Actual-tree calibration"),
            code(
                """
simulation = read_csv("results/data/corrected/path_analysis/corrected_path_simulation_calibration_analysis18_v1.csv")
edge_calibration = read_csv("results/data/corrected/path_analysis/corrected_path_simulation_edge_calibration_analysis18_v1.csv")
display(simulation)
display(edge_calibration)
show_figure("results/figures/corrected/path_analysis/corrected_path_simulation_calibration_analysis18_v1.png")
"""
            ),
            source_cell(
                [
                    ("corrected input builder", "scripts/processing/build_corrected_path_inputs.py"),
                    ("path audit runner", "scripts/processing/audit_corrected_path_models.R"),
                    ("path scaffold", "path_analysis/scripts/path_model_scaffold.R"),
                    ("candidate model record", "path_analysis/CANDIDATE_MODELS.md"),
                ]
            ),
            md(
                """
## Interpretation and review boundary

- Report global fit and CICc/model weights before discussing individual paths.
- Stability across trees, measurement choices, and species omissions is necessary
  but does not establish causality.
- Relative nuclear-IOD algebraically contains nuclear area, so their path is partly
  mechanically coupled and cannot independently validate a nucleotypic mechanism.
- TE evenness and TE log-ratios are derived from the same closed composition.
- At n = 18, simulation calibration and uncertainty are part of the result, not
  optional supplementary detail.
"""
            ),
        ],
    )


BUILDERS = {
    "01_phylogeny_tree_trimming.ipynb": phylogeny_notebook,
    "02_data_tables_and_provenance.ipynb": data_tables_notebook,
    "03_repeat_analysis_te34.ipynb": repeat_notebook,
    "04_ltr_deletion_footprint.ipynb": ltr_notebook,
    "05_cell_modeling_and_measurement.ipynb": cell_notebook,
    "08_integrated_phylogenetic_path_analysis.ipynb": integrated_path_notebook,
}


def build_genome_size_notebook() -> None:
    """Rebuild the lightweight audited IOD products and their canonical notebook."""

    script_dir = ROOT / "path_analysis" / "scripts"
    if str(script_dir) not in sys.path:
        sys.path.insert(0, str(script_dir))
    import build_frozen_genome_iod_notebook as genome_notebook
    import build_frozen_genome_iod_phylogeny_figure as genome_figure

    genome_notebook.build_analysis_outputs(n_bootstrap=2_000)
    genome_figure.build(n_bootstrap=2_000, seed=20260710)
    genome_notebook.write_notebook()
    print(f"Wrote {genome_notebook.NOTEBOOK_PATH.relative_to(ROOT)}")


def build_cell_nucleus_genome_path_outputs() -> None:
    """Freeze the all-finalized-species design and rebuild fitted-result displays."""

    script_dir = ROOT / "path_analysis" / "scripts"
    if str(script_dir) not in sys.path:
        sys.path.insert(0, str(script_dir))
    import build_cell_nucleus_genome_path_notebook as path_design
    import build_cell_nucleus_genome_path_presentation as path_presentation

    path_design.build_outputs(n_bootstrap=250, seed=20260710)
    print(f"Wrote {path_design.OUTPUT_DIR.relative_to(ROOT)}")
    if path_presentation.ANALYSIS_MANIFEST_PATH.exists():
        figure_manifest = path_presentation.build()
        conclusion = json.loads(path_presentation.CONCLUSION_PATH.read_text())
        path_presentation.write_notebook(conclusion)
        if not figure_manifest["release_gates_passed"]:
            raise RuntimeError("All-species path figure release gate did not pass")
        print(f"Wrote {path_presentation.NOTEBOOK_PATH.relative_to(ROOT)}")
    else:
        print(
            "All-species fitted outputs are absent; run the R path analysis and then "
            "build_cell_nucleus_genome_path_presentation.py."
        )


def main() -> None:
    build_frozen_bundle()
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    build_genome_size_notebook()
    build_cell_nucleus_genome_path_outputs()
    for filename, builder in BUILDERS.items():
        output = OUTPUT_DIR / filename
        nbf.write(builder(), output)
        print(f"Wrote {output.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
