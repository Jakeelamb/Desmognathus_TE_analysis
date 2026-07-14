#!/usr/bin/env python3
"""Plot the measured genome-IOD, nucleus, and cell panels beside the time tree."""

from __future__ import annotations

import copy
import hashlib
import json
import sys
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import numpy as np
import pandas as pd
from scipy.optimize import minimize_scalar
from scipy.stats import pearsonr, spearmanr, t as student_t


PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT / "path_analysis" / "scripts"))
import build_frozen_genome_iod_notebook as genome_analysis  # noqa: E402


GENOME_PATH = genome_analysis.FROZEN_PATH
GENOME_SUMMARY_PATH = genome_analysis.SPECIES_SUMMARY_PATH
SOURCE_CELL_PATH = (
    PROJECT_ROOT
    / "path_analysis"
    / "data"
    / "external"
    / "derived"
    / "largest_cell_mask_review"
    / "finalized_largest_cell_masks_all_reviewed_species.csv.gz"
)
SOURCE_TREE_PATH = PROJECT_ROOT / "input_data" / "phylogeny" / "desmo900dated_test.tre"
CELL_PATH = genome_analysis.resolve_portable_artifact(SOURCE_CELL_PATH)
TREE_PATH = genome_analysis.resolve_portable_artifact(SOURCE_TREE_PATH)
OUTPUT_DIR = genome_analysis.REPORT_DIR
FIGURE_DIR = genome_analysis.FIGURE_DIR
PNG_PATH = FIGURE_DIR / "07_measured_phylogeny_genome_nucleus_cell.png"
PDF_PATH = FIGURE_DIR / "07_measured_phylogeny_genome_nucleus_cell.pdf"
PAIRWISE_PNG_PATH = FIGURE_DIR / "08_pairwise_genome_nucleus_cell_relationships.png"
PAIRWISE_PDF_PATH = FIGURE_DIR / "08_pairwise_genome_nucleus_cell_relationships.pdf"
SUMMARY_PATH = OUTPUT_DIR / "phylogeny_genome_nucleus_cell_summary.csv"
CORRELATION_PATH = OUTPUT_DIR / "phylogeny_genome_nucleus_cell_correlations.csv"
MANIFEST_PATH = OUTPUT_DIR / "phylogeny_genome_nucleus_cell_manifest.json"
EXECUTED_NOTEBOOK_PATH = genome_analysis.EXECUTED_NOTEBOOK_PATH
METRIC_ORDER = (
    "genome_size_pg_fuscus_anchored",
    "nucleus_area_um2",
    "cell_area_um2",
)
PAIRWISE_SPECS = (
    (
        "genome_size_vs_nucleus_area",
        "genome_size_pg_fuscus_anchored",
        "nucleus_area_um2",
        "Genome-size estimate (pg; D. fuscus = 16.36)",
        "Nucleus area (µm²)",
        "#4C78A8",
    ),
    (
        "genome_size_vs_cell_area",
        "genome_size_pg_fuscus_anchored",
        "cell_area_um2",
        "Genome-size estimate (pg; D. fuscus = 16.36)",
        "Cell area (µm²)",
        "#D56A4A",
    ),
    (
        "nucleus_area_vs_cell_area",
        "nucleus_area_um2",
        "cell_area_um2",
        "Nucleus area (µm²)",
        "Cell area (µm²)",
        "#765A9A",
    ),
)
METRIC_STYLE = {
    "genome_size_pg_fuscus_anchored": {
        "title": "Genome size",
        "xlabel": "Genome-size estimate (pg)\nD. fuscus anchor = 16.36 pg",
        "color": "#78A6C0",
        "format": ".1f",
    },
    "nucleus_area_um2": {
        "title": "Nucleus area",
        "xlabel": "Nucleus area (µm²)\nmedian of finalized vetted cells",
        "color": "#9A8CB4",
        "format": ".1f",
    },
    "cell_area_um2": {
        "title": "Cell area",
        "xlabel": "Cell area (µm²)\nmedian of finalized vetted cells",
        "color": "#D39768",
        "format": ".1f",
    },
}
PDF_METADATA = {
    "Creator": "Desmognathus_TE frozen microscopy analysis",
    "CreationDate": None,
    "ModDate": None,
}


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def canonical_species(value: Any) -> str:
    text = str(value).strip()
    if text.startswith("D. "):
        return text
    if text.startswith("Desmognathus "):
        return "D. " + text.split(" ", 1)[1]
    return "D. " + text


def bootstrap_equal_image_iod(
    frame: pd.DataFrame,
    *,
    n_bootstrap: int,
    seed: int,
) -> np.ndarray:
    """Draw the mean of image medians by resampling images then nuclei."""
    images = {
        str(filename): group["nuc_iod"].dropna().to_numpy(dtype=float)
        for filename, group in frame.groupby("filename", sort=True)
    }
    images = {name: values for name, values in images.items() if len(values)}
    if not images:
        raise ValueError("No usable image-specific IOD values.")
    names = np.array(sorted(images), dtype=object)
    rng = np.random.default_rng(seed)
    draws = np.empty(n_bootstrap, dtype=float)
    for replicate in range(n_bootstrap):
        sampled_names = rng.choice(names, size=len(names), replace=True)
        medians = []
        for name in sampled_names:
            values = images[str(name)]
            medians.append(
                float(np.median(rng.choice(values, size=len(values), replace=True)))
            )
        draws[replicate] = float(np.mean(medians))
    return draws


def bootstrap_paired_size_medians(
    frame: pd.DataFrame,
    *,
    n_bootstrap: int,
    seed: int,
) -> pd.DataFrame:
    """Bootstrap paired cell and corresponding-nucleus medians."""
    values = frame[["cell_area_um2", "nuc_area_um2"]].to_numpy(dtype=float)
    if len(values) == 0 or not np.isfinite(values).all():
        raise ValueError("Size panel contains no usable paired measurements.")
    rng = np.random.default_rng(seed)
    indices = rng.integers(0, len(values), size=(n_bootstrap, len(values)))
    sampled = values[indices]
    medians = np.median(sampled, axis=1)
    return pd.DataFrame(
        {
            "cell_area_um2": medians[:, 0],
            "nuc_area_um2": medians[:, 1],
        }
    )


def build_figure_data(
    *,
    n_bootstrap: int = 2_000,
    seed: int = 20260710,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    genome = pd.read_csv(GENOME_PATH, low_memory=False)
    genome_summary = pd.read_csv(GENOME_SUMMARY_PATH, low_memory=False)
    cells = pd.read_csv(CELL_PATH, low_memory=False)

    genome_species = set(genome["species"])
    cell_species = set(cells["species"])
    figure_species = sorted(cell_species)
    if not genome_species.issubset(cell_species):
        raise ValueError(
            "Every finalized genome-panel species must have a finalized size panel."
        )
    if not figure_species:
        raise ValueError("Finalized size panel contains no measured species.")
    if cells.groupby("species").size().min() < 1:
        raise ValueError("Every size-panel species must retain a reviewed cell pair.")
    if not cells["review_status"].eq("reviewed_keep").all():
        raise ValueError("Size panel contains a row that is not reviewed keep.")
    if not genome["review_status"].eq("reviewed_keep").all():
        raise ValueError("Genome panel contains a row that is not reviewed keep.")

    genome_points = genome_summary.set_index("species")
    reference_species = genome_analysis.REFERENCE_SPECIES
    if reference_species not in genome_species:
        raise ValueError(f"Genome panel is missing reference species {reference_species}.")
    genome_seed_index = {
        species: index for index, species in enumerate(sorted(genome_species))
    }
    reference_iod = float(
        genome_points.loc[reference_species, "iod_equal_image_estimate"]
    )
    if reference_iod <= 0:
        raise ValueError("Reference-species IOD estimate must be positive.")
    reference_draws = bootstrap_equal_image_iod(
        genome.loc[genome["species"].eq(reference_species)],
        n_bootstrap=n_bootstrap,
        seed=seed + genome_seed_index[reference_species] * 1009,
    )

    summary_rows: list[dict[str, Any]] = []
    draw_frames: list[pd.DataFrame] = []
    for species_index, species in enumerate(figure_species):
        genome_group = genome.loc[genome["species"].eq(species)]
        cell_group = cells.loc[cells["species"].eq(species)]
        size_draws = bootstrap_paired_size_medians(
            cell_group,
            n_bootstrap=n_bootstrap,
            seed=seed + species_index * 2017 + 1,
        )
        metric_draws = {
            "nucleus_area_um2": size_draws["nuc_area_um2"].to_numpy(),
            "cell_area_um2": size_draws["cell_area_um2"].to_numpy(),
        }
        if species in genome_species:
            species_iod_draws = bootstrap_equal_image_iod(
                genome_group,
                n_bootstrap=n_bootstrap,
                seed=seed + genome_seed_index[species] * 1009,
            )
            metric_draws["genome_size_pg_fuscus_anchored"] = (
                species_iod_draws
                / reference_draws
                * genome_analysis.REFERENCE_GENOME_SIZE_PG
            )
            genome_point = float(
                genome_points.loc[species, "iod_equal_image_estimate"]
                / reference_iod
                * genome_analysis.REFERENCE_GENOME_SIZE_PG
            )
            quality_statuses = sorted(
                genome_group["quality_match_status"].dropna().astype(str).unique()
            )
            if quality_statuses == ["common_support"]:
                genome_status = "finalized_included_common_support"
                include_in_primary_genome_analysis = True
            elif quality_statuses == ["limited_overlap"]:
                genome_status = "finalized_included_limited_overlap"
                include_in_primary_genome_analysis = True
            else:
                raise ValueError(
                    f"Unexpected quality-match status for {species}: {quality_statuses}"
                )
        else:
            genome_point = np.nan
            genome_status = "no_finalized_genome_iod"
            include_in_primary_genome_analysis = False
        points = {
            "genome_size_pg_fuscus_anchored": genome_point,
            "nucleus_area_um2": float(cell_group["nuc_area_um2"].median()),
            "cell_area_um2": float(cell_group["cell_area_um2"].median()),
        }
        row: dict[str, Any] = {
            "species": species,
            "genome_panel_status": genome_status,
            "include_in_primary_genome_analysis": include_in_primary_genome_analysis,
            "n_genome_nuclei": int(len(genome_group)),
            "n_genome_images": int(genome_group["filename"].nunique()),
            "n_size_cells": int(len(cell_group)),
            "n_size_images": int(cell_group["filename"].nunique()),
            "size_sample_size_support": (
                "below_top50_target" if len(cell_group) < 50 else "top50_target_met"
            ),
        }
        for metric in METRIC_ORDER:
            row[metric] = points[metric]
            if metric not in metric_draws:
                row[f"{metric}_ci_low"] = np.nan
                row[f"{metric}_ci_high"] = np.nan
                continue
            low, high = np.quantile(metric_draws[metric], [0.025, 0.975])
            row[f"{metric}_ci_low"] = float(low)
            row[f"{metric}_ci_high"] = float(high)
            draw_frames.append(
                pd.DataFrame(
                    {
                        "species": species,
                        "metric": metric,
                        "bootstrap_replicate": np.arange(1, n_bootstrap + 1),
                        "value": metric_draws[metric],
                    }
                )
            )
        summary_rows.append(row)

    summary = pd.DataFrame(summary_rows)
    draws = pd.concat(draw_frames, ignore_index=True)
    correlation_rows = []
    for label, left, right, _x_label, _y_label, _color in PAIRWISE_SPECS:
        analysis_columns = ["species", left, right]
        if "genome_size_pg_fuscus_anchored" in {left, right}:
            complete = summary.loc[
                summary["include_in_primary_genome_analysis"], analysis_columns
            ].dropna()
        else:
            complete = summary[analysis_columns].dropna()
        rho, p_value = spearmanr(complete[left], complete[right])
        pearson_r, pearson_p = pearsonr(complete[left], complete[right])
        pgls = fit_pagel_lambda_pgls(
            complete,
            predictor=left,
            outcome=right,
        )
        correlation_rows.append(
            {
                "comparison": label,
                "left_metric": left,
                "right_metric": right,
                "spearman_rho": float(rho),
                "p_value_descriptive_only": float(p_value),
                "pearson_r": float(pearson_r),
                "pearson_p_value_descriptive_only": float(pearson_p),
                "n_species": int(len(complete)),
                "phylogenetically_corrected": True,
                **pgls,
            }
        )
    return summary, draws, pd.DataFrame(correlation_rows)


def load_measured_tree(species: list[str]) -> Any:
    try:
        from Bio import Phylo
        from Bio.Phylo.BaseTree import Tree
    except ModuleNotFoundError as error:
        raise ModuleNotFoundError(
            "Biopython is required for tree pruning. Run with: uv run --with biopython"
        ) from error

    full_tree = Phylo.read(TREE_PATH, "newick")
    wanted_tips = {value.replace("D. ", "") for value in species}
    available = {terminal.name for terminal in full_tree.get_terminals()}
    missing = sorted(wanted_tips - available)
    if missing:
        raise ValueError(f"Measured species are absent from the time tree: {missing}")
    mrca = full_tree.common_ancestor(
        [terminal for terminal in full_tree.get_terminals() if terminal.name in wanted_tips]
    )
    tree = Tree(root=copy.deepcopy(mrca), rooted=True)
    tree.root.branch_length = 0.0
    for terminal in list(tree.get_terminals()):
        if terminal.name not in wanted_tips:
            tree.prune(terminal)
    observed = {terminal.name for terminal in tree.get_terminals()}
    if observed != wanted_tips:
        raise ValueError("Pruned time tree does not match the measured species set.")
    tip_depths = [tree.distance(tree.root, terminal) for terminal in tree.get_terminals()]
    if max(tip_depths) - min(tip_depths) > 1e-4:
        raise ValueError("Pruned time tree is not ultrametric within tolerance.")
    return tree


def phylogenetic_brownian_covariance(
    tree: Any,
    tip_order: list[str],
) -> np.ndarray:
    """Return the Brownian shared-path covariance for an ordered tip set."""
    terminals = {terminal.name: terminal for terminal in tree.get_terminals()}
    if set(terminals) != set(tip_order):
        raise ValueError("Tree tips do not match the requested covariance order.")
    covariance = np.empty((len(tip_order), len(tip_order)), dtype=float)
    for row_index, row_tip in enumerate(tip_order):
        for column_index, column_tip in enumerate(tip_order):
            if row_index == column_index:
                target = terminals[row_tip]
            else:
                target = tree.common_ancestor(
                    terminals[row_tip], terminals[column_tip]
                )
            covariance[row_index, column_index] = float(
                tree.distance(tree.root, target)
            )
    if not np.isfinite(covariance).all():
        raise ValueError("Phylogenetic covariance contains non-finite values.")
    if np.any(np.diag(covariance) <= 0):
        raise ValueError("Phylogenetic covariance requires positive tip depths.")
    return covariance


def fit_pagel_lambda_pgls(
    frame: pd.DataFrame,
    *,
    predictor: str,
    outcome: str,
) -> dict[str, Any]:
    """Fit ML Pagel-lambda PGLS to log10, sample-standardized traits.

    The likelihood and residual degrees-of-freedom standard error reproduce the
    ``phylolm(..., model = "lambda")`` point fits used by the path analysis.
    """
    data = frame[["species", predictor, outcome]].dropna().copy()
    if len(data) < 4:
        raise ValueError("PGLS requires at least four complete species.")
    if data["species"].duplicated().any():
        raise ValueError("PGLS input contains duplicate species.")
    if (data[[predictor, outcome]] <= 0).any().any():
        raise ValueError("PGLS log10 traits must be strictly positive.")

    tree = load_measured_tree(data["species"].tolist())
    tip_order = [terminal.name for terminal in tree.get_terminals()]
    data["tree_tip"] = data["species"].str.replace("D. ", "", regex=False)
    data = data.set_index("tree_tip").loc[tip_order]
    brownian = phylogenetic_brownian_covariance(tree, tip_order)
    diagonal = np.diag(np.diag(brownian))
    off_diagonal = brownian - diagonal

    predictor_log10 = np.log10(data[predictor].to_numpy(dtype=float))
    outcome_log10 = np.log10(data[outcome].to_numpy(dtype=float))
    predictor_mean = float(predictor_log10.mean())
    predictor_sd = float(predictor_log10.std(ddof=1))
    outcome_mean = float(outcome_log10.mean())
    outcome_sd = float(outcome_log10.std(ddof=1))
    if predictor_sd <= 0 or outcome_sd <= 0:
        raise ValueError("PGLS cannot standardize a zero-variance trait.")
    predictor_z = (predictor_log10 - predictor_mean) / predictor_sd
    outcome_z = (outcome_log10 - outcome_mean) / outcome_sd
    design = np.column_stack([np.ones(len(data)), predictor_z])

    def profile_fit(lambda_value: float) -> dict[str, Any]:
        covariance = diagonal + float(lambda_value) * off_diagonal
        sign, log_determinant = np.linalg.slogdet(covariance)
        if sign <= 0:
            raise ValueError("Pagel-lambda covariance is not positive definite.")
        covariance_inverse_design = np.linalg.solve(covariance, design)
        information = design.T @ covariance_inverse_design
        covariance_inverse_outcome = np.linalg.solve(covariance, outcome_z)
        coefficients = np.linalg.solve(
            information,
            design.T @ covariance_inverse_outcome,
        )
        residual = outcome_z - design @ coefficients
        quadratic = float(residual @ np.linalg.solve(covariance, residual))
        n_species = len(data)
        profile_nll = 0.5 * (
            n_species * np.log(2.0 * np.pi)
            + log_determinant
            + n_species * np.log(quadratic / n_species)
            + n_species
        )
        return {
            "profile_nll": float(profile_nll),
            "coefficients": coefficients,
            "quadratic": quadratic,
            "information": information,
            "covariance": covariance,
        }

    lambda_floor = 1e-7
    optimization = minimize_scalar(
        lambda value: profile_fit(value)["profile_nll"],
        bounds=(lambda_floor, 1.0),
        method="bounded",
        options={"xatol": 1e-13},
    )
    candidates = [
        (lambda_floor, profile_fit(lambda_floor)),
        (float(optimization.x), profile_fit(float(optimization.x))),
        (1.0, profile_fit(1.0)),
    ]
    lambda_estimate, fit = min(
        candidates,
        key=lambda item: item[1]["profile_nll"],
    )
    coefficients = np.asarray(fit["coefficients"], dtype=float)
    degrees_freedom = len(data) - design.shape[1]
    residual_variance = fit["quadratic"] / degrees_freedom
    coefficient_covariance = residual_variance * np.linalg.inv(fit["information"])
    standard_errors = np.sqrt(np.diag(coefficient_covariance))
    t_value = float(coefficients[1] / standard_errors[1])
    p_value = float(2.0 * student_t.sf(abs(t_value), degrees_freedom))
    critical_value = float(student_t.ppf(0.975, degrees_freedom))

    covariance_inverse = np.linalg.inv(fit["covariance"])
    weighted_mean = float(
        (np.ones(len(data)) @ covariance_inverse @ outcome_z)
        / (np.ones(len(data)) @ covariance_inverse @ np.ones(len(data)))
    )
    null_residual = outcome_z - weighted_mean
    null_quadratic = float(null_residual @ covariance_inverse @ null_residual)
    generalized_r_squared = 1.0 - fit["quadratic"] / null_quadratic

    return {
        "pgls_model": "Pagel_lambda_ML",
        "trait_transformation": "log10_then_sample_zscore",
        "pgls_standardized_beta": float(coefficients[1]),
        "pgls_standard_error": float(standard_errors[1]),
        "pgls_ci_low": float(
            coefficients[1] - critical_value * standard_errors[1]
        ),
        "pgls_ci_high": float(
            coefficients[1] + critical_value * standard_errors[1]
        ),
        "pgls_t_value": t_value,
        "pgls_degrees_freedom": int(degrees_freedom),
        "pgls_p_value": p_value,
        "pagel_lambda": float(lambda_estimate),
        "pagel_lambda_at_lower_boundary": bool(lambda_estimate <= 1.01e-7),
        "pgls_generalized_r_squared": float(generalized_r_squared),
        "pgls_intercept_z": float(coefficients[0]),
        "pgls_var_intercept": float(coefficient_covariance[0, 0]),
        "pgls_cov_intercept_slope": float(coefficient_covariance[0, 1]),
        "pgls_var_slope": float(coefficient_covariance[1, 1]),
        "predictor_log10_mean": predictor_mean,
        "predictor_log10_sd": predictor_sd,
        "outcome_log10_mean": outcome_mean,
        "outcome_log10_sd": outcome_sd,
    }


def tree_coordinates(tree: Any) -> tuple[dict[Any, float], dict[Any, float], list[Any]]:
    terminals = list(tree.get_terminals())
    depths = tree.depths()
    max_depth = max(float(depths[terminal]) for terminal in terminals)
    ages = {clade: max_depth - float(depth) for clade, depth in depths.items()}
    y_coordinates: dict[Any, float] = {
        terminal: float(len(terminals) - 1 - index)
        for index, terminal in enumerate(terminals)
    }

    def assign_internal(clade: Any) -> float:
        if clade in y_coordinates:
            return y_coordinates[clade]
        child_values = [assign_internal(child) for child in clade.clades]
        y_coordinates[clade] = float((min(child_values) + max(child_values)) / 2.0)
        return y_coordinates[clade]

    assign_internal(tree.root)
    return ages, y_coordinates, terminals


def draw_tree_axis(
    ax: Any,
    tree: Any,
) -> tuple[list[str], dict[str, float], float]:
    ages, y_coordinates, terminals = tree_coordinates(tree)
    for parent in tree.find_clades(order="preorder"):
        if not parent.clades:
            continue
        child_y = [y_coordinates[child] for child in parent.clades]
        ax.plot(
            [ages[parent], ages[parent]],
            [min(child_y), max(child_y)],
            color="#30343B",
            linewidth=1.2,
        )
        for child in parent.clades:
            ax.plot(
                [ages[parent], ages[child]],
                [y_coordinates[child], y_coordinates[child]],
                color="#30343B",
                linewidth=1.2,
            )
    species_order = []
    y_by_species = {}
    max_age = max(ages.values())
    for terminal in terminals:
        species = canonical_species(terminal.name)
        species_order.append(species)
        y_by_species[species] = y_coordinates[terminal]
        ax.text(
            -0.018 * max_age,
            y_coordinates[terminal],
            species,
            va="center",
            ha="left",
            fontsize=8.6,
            color="#25282D",
            fontstyle="italic",
        )
    ax.set_xlim(max_age * 1.025, -max_age * 0.22)
    ax.set_xticks(np.linspace(0.0, max_age, 5))
    ax.set_ylim(-0.8, len(terminals) - 0.2)
    ax.set_yticks([])
    ax.set_xlabel("Time before present (MYA)")
    ax.set_title("Time-calibrated phylogeny", fontsize=12, fontweight="bold")
    ax.grid(axis="x", alpha=0.16, color="#98A1A8")
    for spine in ["top", "right", "left"]:
        ax.spines[spine].set_visible(False)
    return species_order, y_by_species, max_age


def draw_trait_axis(
    ax: Any,
    *,
    metric: str,
    summary: pd.DataFrame,
    draws: pd.DataFrame,
    species_order: list[str],
    y_by_species: dict[str, float],
) -> None:
    style = METRIC_STYLE[metric]
    available_species = []
    distributions = []
    for species in species_order:
        values = draws.loc[
            draws["species"].eq(species) & draws["metric"].eq(metric), "value"
        ].to_numpy()
        if len(values):
            available_species.append(species)
            distributions.append(values)
    positions = [y_by_species[species] for species in available_species]
    violins = ax.violinplot(
        distributions,
        positions=positions,
        vert=False,
        widths=0.72,
        showmeans=False,
        showmedians=False,
        showextrema=False,
        bw_method=0.25,
    )
    for body in violins["bodies"]:
        body.set_facecolor(style["color"])
        body.set_edgecolor("#52616B")
        body.set_linewidth(0.8)
        body.set_alpha(0.44)

    indexed_summary = summary.set_index("species")
    points = indexed_summary[metric]
    all_values = np.concatenate(distributions)
    lower = float(np.quantile(all_values, 0.002))
    upper = float(np.quantile(all_values, 0.998))
    span = max(upper - lower, abs(upper) * 0.05, 1e-6)
    label_offset = span * 0.018
    for species in species_order:
        value = float(points.loc[species])
        y = y_by_species[species]
        if not np.isfinite(value):
            ax.text(
                0.02,
                y,
                "no finalized IOD estimate",
                transform=ax.get_yaxis_transform(),
                va="center",
                ha="left",
                fontsize=6.6,
                color="#7A7F84",
                fontstyle="italic",
            )
            continue
        limited_overlap = bool(
            metric == "genome_size_pg_fuscus_anchored"
            and indexed_summary.loc[species, "genome_panel_status"]
            == "finalized_included_limited_overlap"
        )
        ax.scatter(
            value,
            y,
            s=30 if limited_overlap else 24,
            marker="D" if limited_overlap else "o",
            facecolor="none" if limited_overlap else "#24282D",
            edgecolor="#9A4D45" if limited_overlap else "#24282D",
            linewidth=1.1 if limited_overlap else 0.6,
            zorder=4,
        )
        ax.text(
            value + label_offset,
            y,
            format(value, style["format"]) + ("†" if limited_overlap else ""),
            va="center",
            ha="left",
            fontsize=7,
            color="#454A50",
        )
    reference = float(points.median())
    ax.axvline(reference, color="#70777E", linestyle="--", linewidth=0.9, alpha=0.65)
    ax.set_xlim(lower - span * 0.06, upper + span * 0.23)
    ax.set_ylim(-0.8, len(species_order) - 0.2)
    ax.set_yticks([])
    ax.set_title(style["title"], fontsize=12, fontweight="bold")
    ax.set_xlabel(style["xlabel"])
    ax.grid(axis="x", alpha=0.16, color="#98A1A8")
    for spine in ["top", "right", "left"]:
        ax.spines[spine].set_visible(False)


def render_figure(
    summary: pd.DataFrame,
    draws: pd.DataFrame,
    correlations: pd.DataFrame,
) -> tuple[list[str], float]:
    tree = load_measured_tree(summary["species"].tolist())
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(
        1,
        4,
        figsize=(21, 12.5),
        gridspec_kw={"width_ratios": [1.55, 1.12, 1.12, 1.12], "wspace": 0.035},
    )
    species_order, y_by_species, max_age = draw_tree_axis(axes[0], tree)
    summary_order = set(summary["species"])
    if set(species_order) != summary_order:
        raise ValueError("Tree order does not contain the exact measured species set.")
    for ax, metric in zip(axes[1:], METRIC_ORDER):
        draw_trait_axis(
            ax,
            metric=metric,
            summary=summary,
            draws=draws,
            species_order=species_order,
            y_by_species=y_by_species,
        )

    correlation_lookup = correlations.set_index("comparison")
    fig.suptitle(
        "Time-calibrated Desmognathus phylogeny with audited bootstrap trait distributions",
        fontsize=16,
        fontweight="bold",
        y=0.978,
    )
    fig.text(
        0.5,
        0.946,
        (
            f"Size panel: {len(summary)} species with all finalized vetted cell–nucleus pairs "
            f"(n = {int(summary['n_size_cells'].min())}–{int(summary['n_size_cells'].max())} per species). "
            f"Genome estimates: {int(summary['genome_size_pg_fuscus_anchored'].notna().sum())} species; "
            f"analysis includes all {int(summary['include_in_primary_genome_analysis'].sum())}. "
            "Hollow diamonds mark limited-overlap quality matches; no phylogenetic fills."
        ),
        ha="center",
        va="center",
        fontsize=9.5,
        color="#4A4F55",
    )
    fig.text(
        0.5,
        0.925,
        (
            "Phylogenetically corrected standardized PGLS β: "
            f"genome–nucleus {correlation_lookup.loc['genome_size_vs_nucleus_area', 'pgls_standardized_beta']:.2f}; "
            f"genome–cell {correlation_lookup.loc['genome_size_vs_cell_area', 'pgls_standardized_beta']:.2f}; "
            f"nucleus–cell {correlation_lookup.loc['nucleus_area_vs_cell_area', 'pgls_standardized_beta']:.2f}. "
            "Pagel's λ was estimated by maximum likelihood for each relationship."
        ),
        ha="center",
        va="center",
        fontsize=9.2,
        color="#4A4F55",
    )
    legend_handles = [
        Patch(
            facecolor="#8EABC0",
            edgecolor="#52616B",
            alpha=0.44,
            label="Bootstrap estimator distribution",
        ),
        Line2D(
            [],
            [],
            marker="o",
            linestyle="",
            color="#24282D",
            label="Current point estimate",
        ),
        Line2D(
            [],
            [],
            linestyle="--",
            color="#70777E",
            label="Across-species median",
        ),
    ]
    fig.legend(
        handles=legend_handles,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.907),
        ncol=3,
        frameon=False,
        fontsize=8.7,
    )
    fig.subplots_adjust(top=0.875, bottom=0.07, left=0.035, right=0.99)
    fig.savefig(PNG_PATH, dpi=220, bbox_inches="tight", facecolor="white")
    fig.savefig(
        PDF_PATH,
        bbox_inches="tight",
        facecolor="white",
        metadata=PDF_METADATA,
    )
    plt.close(fig)
    return species_order, max_age


def render_pairwise_figure(
    summary: pd.DataFrame,
    correlations: pd.DataFrame,
) -> None:
    """Plot all pairwise traits with Pagel-lambda PGLS fits and uncertainty."""
    correlation_lookup = correlations.set_index("comparison")
    fig, axes = plt.subplots(1, 3, figsize=(18, 6.7))
    for panel_index, (
        comparison,
        x_column,
        y_column,
        x_label,
        y_label,
        color,
    ) in enumerate(PAIRWISE_SPECS):
        ax = axes[panel_index]
        columns = [
            "species",
            "include_in_primary_genome_analysis",
            "genome_panel_status",
            x_column,
            f"{x_column}_ci_low",
            f"{x_column}_ci_high",
            y_column,
            f"{y_column}_ci_low",
            f"{y_column}_ci_high",
        ]
        data = summary[columns].dropna().copy()
        x = data[x_column].to_numpy(dtype=float)
        y = data[y_column].to_numpy(dtype=float)
        x_error = np.vstack(
            [
                x - data[f"{x_column}_ci_low"].to_numpy(dtype=float),
                data[f"{x_column}_ci_high"].to_numpy(dtype=float) - x,
            ]
        )
        y_error = np.vstack(
            [
                y - data[f"{y_column}_ci_low"].to_numpy(dtype=float),
                data[f"{y_column}_ci_high"].to_numpy(dtype=float) - y,
            ]
        )
        ax.errorbar(
            x,
            y,
            xerr=x_error,
            yerr=y_error,
            fmt="none",
            ecolor=color,
            elinewidth=1.0,
            alpha=0.28,
            capsize=0,
            zorder=1,
        )
        genome_comparison = "genome_size_pg_fuscus_anchored" in {x_column, y_column}
        primary_mask = (
            data["include_in_primary_genome_analysis"].to_numpy(dtype=bool)
            if genome_comparison
            else np.ones(len(data), dtype=bool)
        )
        ax.scatter(
            x[primary_mask],
            y[primary_mask],
            s=48,
            color=color,
            edgecolor="white",
            linewidth=0.7,
            alpha=0.92,
            zorder=3,
        )
        if genome_comparison and (~primary_mask).any():
            ax.scatter(
                x[~primary_mask],
                y[~primary_mask],
                s=58,
                marker="D",
                facecolors="none",
                edgecolors="#6F7479",
                linewidth=1.2,
                zorder=4,
            )
        stats_row = correlation_lookup.loc[comparison]
        x_line = np.linspace(float(x.min()), float(x.max()), 200)
        predictor_z = (
            np.log10(x_line) - float(stats_row["predictor_log10_mean"])
        ) / float(stats_row["predictor_log10_sd"])
        predicted_z = (
            float(stats_row["pgls_intercept_z"])
            + float(stats_row["pgls_standardized_beta"]) * predictor_z
        )
        design_line = np.column_stack([np.ones(len(x_line)), predictor_z])
        coefficient_covariance = np.array(
            [
                [
                    float(stats_row["pgls_var_intercept"]),
                    float(stats_row["pgls_cov_intercept_slope"]),
                ],
                [
                    float(stats_row["pgls_cov_intercept_slope"]),
                    float(stats_row["pgls_var_slope"]),
                ],
            ]
        )
        mean_variance = np.einsum(
            "ij,jk,ik->i",
            design_line,
            coefficient_covariance,
            design_line,
        )
        model_se = np.sqrt(np.maximum(mean_variance, 0.0))
        critical_value = float(
            student_t.ppf(
                0.975,
                int(stats_row["pgls_degrees_freedom"]),
            )
        )
        outcome_mean = float(stats_row["outcome_log10_mean"])
        outcome_sd = float(stats_row["outcome_log10_sd"])
        y_line = 10.0 ** (outcome_mean + outcome_sd * predicted_z)
        y_line_low = 10.0 ** (
            outcome_mean + outcome_sd * (predicted_z - critical_value * model_se)
        )
        y_line_high = 10.0 ** (
            outcome_mean + outcome_sd * (predicted_z + critical_value * model_se)
        )
        ax.fill_between(
            x_line,
            y_line_low,
            y_line_high,
            color=color,
            alpha=0.12,
            linewidth=0,
            zorder=1,
        )
        ax.plot(
            x_line,
            y_line,
            color=color,
            linewidth=2.0,
            alpha=0.82,
            zorder=2,
        )
        for label_index, row in enumerate(data.itertuples(index=False)):
            x_offset = 4 if label_index % 2 == 0 else -4
            horizontal_alignment = "left" if x_offset > 0 else "right"
            y_offset = 4 + ((label_index % 3) - 1) * 4
            ax.annotate(
                str(row.species).replace("D. ", ""),
                (getattr(row, x_column), getattr(row, y_column)),
                xytext=(x_offset, y_offset),
                textcoords="offset points",
                ha=horizontal_alignment,
                va="center",
                fontsize=7,
                color="#353A40",
                fontstyle="italic",
            )
        lambda_text = (
            "<0.001"
            if float(stats_row["pagel_lambda"]) < 0.001
            else f"{float(stats_row['pagel_lambda']):.2f}"
        )
        ax.set_title(
            (
                f"PGLS βstd = {stats_row['pgls_standardized_beta']:.2f} "
                f"[{stats_row['pgls_ci_low']:.2f}, {stats_row['pgls_ci_high']:.2f}]; "
                f"P = {stats_row['pgls_p_value']:.3f}\n"
                f"Pagel's λ {lambda_text}; n = {int(stats_row['n_species'])}"
            ),
            fontsize=10.0,
            fontweight="bold",
        )
        ax.set_xlabel(x_label)
        ax.set_ylabel(y_label)
        ax.grid(alpha=0.16)
        for spine in ["top", "right"]:
            ax.spines[spine].set_visible(False)

    fig.suptitle(
        "Phylogenetically corrected relationships among genome size, nucleus area, and cell area",
        fontsize=15,
        fontweight="bold",
        y=0.985,
    )
    fig.text(
        0.5,
        0.94,
        (
            "Points are species estimates; bars are 95% bootstrap intervals conditional on "
            "the finalized image and selected-cell panels. Hollow diamonds mark limited-overlap quality matches; "
            "all 24 species are included in genome PGLS fits. Lines are Pagel-λ PGLS fits; ribbons are 95% confidence intervals."
        ),
        ha="center",
        va="center",
        fontsize=9.2,
        color="#4A4F55",
    )
    fig.text(
        0.5,
        0.015,
        (
            "Genome size is calibrated to D. fuscus = 16.36 pg and is scaled from nuclear IOD, which includes nuclear area algebraically. "
            "All λ estimates reached the near-zero boundary, so shared ancestry did not materially change these slopes; PGLS does not establish causal direction."
        ),
        ha="center",
        va="bottom",
        fontsize=9.0,
        color="#4A4F55",
    )
    legend_handles = [
        Line2D(
            [],
            [],
            marker="o",
            linestyle="",
            color="#59636B",
            label="Species estimate with measurement-bootstrap interval",
        ),
        Line2D(
            [],
            [],
            marker="D",
            linestyle="",
            markerfacecolor="none",
            markeredgecolor="#6F7479",
            label="Limited-overlap quality match (included in fit)",
        ),
        Line2D(
            [],
            [],
            color="#59636B",
            linewidth=2.0,
            label="Pagel-λ PGLS fit",
        ),
        Patch(
            facecolor="#59636B",
            alpha=0.12,
            edgecolor="none",
            label="95% model confidence interval",
        ),
    ]
    fig.legend(
        handles=legend_handles,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.905),
        ncol=4,
        frameon=False,
        fontsize=8.5,
    )
    fig.subplots_adjust(top=0.79, bottom=0.15, left=0.06, right=0.99, wspace=0.23)
    fig.savefig(PAIRWISE_PNG_PATH, dpi=220, bbox_inches="tight", facecolor="white")
    fig.savefig(
        PAIRWISE_PDF_PATH,
        bbox_inches="tight",
        facecolor="white",
        metadata=PDF_METADATA,
    )
    plt.close(fig)


def build(*, n_bootstrap: int = 2_000, seed: int = 20260710) -> dict[str, Any]:
    summary, draws, correlations = build_figure_data(
        n_bootstrap=n_bootstrap,
        seed=seed,
    )
    species_order, tree_height = render_figure(summary, draws, correlations)
    render_pairwise_figure(summary, correlations)
    summary["tree_display_order"] = summary["species"].map(
        {species: index + 1 for index, species in enumerate(species_order)}
    )
    summary = summary.sort_values("tree_display_order").reset_index(drop=True)
    summary.to_csv(SUMMARY_PATH, index=False, float_format="%.8f")
    correlations.to_csv(CORRELATION_PATH, index=False, float_format="%.8f")
    manifest = {
        "analysis": "Measured-only time-tree alignment of audited microscopy traits",
        "n_measured_species": int(len(summary)),
        "n_primary_genome_species": int(
            summary["include_in_primary_genome_analysis"].sum()
        ),
        "n_genome_estimate_species": int(
            summary["genome_size_pg_fuscus_anchored"].notna().sum()
        ),
        "limited_overlap_included_species": sorted(
            summary.loc[
                summary["genome_panel_status"].eq(
                    "finalized_included_limited_overlap"
                ),
                "species",
            ].tolist()
        ),
        "all_finalized_species_included_in_genome_models": True,
        "reduced_size_sample_species": {
            str(row.species): int(row.n_size_cells)
            for row in summary.loc[summary["n_size_cells"].lt(50)].itertuples()
        },
        "missing_primary_genome_species": sorted(
            summary.loc[
                summary["genome_size_pg_fuscus_anchored"].isna(), "species"
            ].tolist()
        ),
        "species": species_order,
        "tree_source": str(TREE_PATH.resolve()),
        "tree_source_sha256": sha256_file(TREE_PATH),
        "pruned_tree_height_mya": float(tree_height),
        "genome_source": str(GENOME_PATH.resolve()),
        "genome_source_sha256": sha256_file(GENOME_PATH),
        "size_source": str(CELL_PATH.resolve()),
        "size_source_sha256": sha256_file(CELL_PATH),
        "bootstrap_replicates": int(n_bootstrap),
        "random_seed": int(seed),
        "genome_estimator": (
            "mean of image-specific median nuclear IOD divided by the D. fuscus "
            "equal-image IOD estimate and multiplied by 16.36 pg"
        ),
        "genome_calibration_kind": genome_analysis.CALIBRATION_KIND,
        "genome_reference_species": genome_analysis.REFERENCE_SPECIES,
        "genome_reference_pg": genome_analysis.REFERENCE_GENOME_SIZE_PG,
        "genome_reference_value_status": (
            "existing project convention; exact literature provenance unresolved"
        ),
        "genome_intervals_propagate_reference_image_uncertainty": True,
        "size_estimators": (
            "pooled medians of all finalized manually vetted literal-largest cells "
            "and their corresponding nuclei; observed species counts are retained"
        ),
        "bootstrap_scope": (
            "genome: images then nuclei within images; size: paired resampling "
            "of the finalized selected cell-nucleus rows"
        ),
        "phylogenetic_fills_used": False,
        "absolute_genome_size_claimed": False,
        "conditional_genome_size_estimates_reported": True,
        "correlations_are_phylogenetically_corrected": True,
        "pairwise_phylogenetic_method": (
            "Pagel-lambda maximum-likelihood PGLS on separately log10-transformed "
            "and sample-standardized species traits"
        ),
        "pairwise_phylogenetic_tree": str(TREE_PATH.resolve()),
        "pairwise_complete_case_species": {
            str(row.comparison): int(row.n_species)
            for row in correlations.itertuples(index=False)
        },
        "pairwise_pagel_lambda": {
            str(row.comparison): float(row.pagel_lambda)
            for row in correlations.itertuples(index=False)
        },
        "pairwise_interpretation_limit": (
            "PGLS adjusts pairwise association for modeled shared ancestry but "
            "does not identify causal direction; genome IOD includes nuclear area algebraically"
        ),
        "summary_csv": str(SUMMARY_PATH.resolve()),
        "summary_sha256": sha256_file(SUMMARY_PATH),
        "correlation_csv": str(CORRELATION_PATH.resolve()),
        "correlation_sha256": sha256_file(CORRELATION_PATH),
        "png": str(PNG_PATH.resolve()),
        "png_sha256": sha256_file(PNG_PATH),
        "pdf": str(PDF_PATH.resolve()),
        "pdf_sha256": sha256_file(PDF_PATH),
        "pairwise_png": str(PAIRWISE_PNG_PATH.resolve()),
        "pairwise_png_sha256": sha256_file(PAIRWISE_PNG_PATH),
        "pairwise_pdf": str(PAIRWISE_PDF_PATH.resolve()),
        "pairwise_pdf_sha256": sha256_file(PAIRWISE_PDF_PATH),
    }
    MANIFEST_PATH.write_text(
        json.dumps(manifest, indent=2) + "\n",
        encoding="utf-8",
    )
    return manifest


def main() -> None:
    print(json.dumps(build(), indent=2))


if __name__ == "__main__":
    main()
