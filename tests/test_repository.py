from __future__ import annotations

import hashlib
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]


def run_python(*arguments: str) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [sys.executable, *arguments],
        cwd=ROOT,
        check=False,
        capture_output=True,
        text=True,
    )


def test_compact_top_level_layout() -> None:
    required = {"analyses", "data", "docs", "Publication", "scripts", "tests"}
    forbidden = {
        "input_data",
        "interim",
        "results",
        "path_analysis",
        "notebooks",
        "outputs",
        "plans",
        "Dusky.yml",
        "paths.yaml",
    }
    assert all((ROOT / path).exists() for path in required)
    assert not any((ROOT / path).exists() for path in forbidden)
    assert not (ROOT / "analyses/01_transposable_elements/output").exists()


def test_repository_validator() -> None:
    result = run_python("-m", "scripts.project", "validate")
    assert result.returncode == 0, result.stdout + result.stderr
    assert "TE34=34; Path24=24; overlap=21" in result.stdout
    assert "identity=PASS" in result.stdout


def test_morphology_release_counts_and_species() -> None:
    panel = pd.read_csv(ROOT / "data/identity/path24_panel.csv")
    morphology = pd.read_csv(
        ROOT / "analyses/02_morphology/data/cell_nucleus_objects.csv",
        low_memory=False,
    )
    iod = pd.read_csv(
        ROOT / "analyses/02_morphology/data/nuclear_iod_objects.csv",
        low_memory=False,
    )
    assert (panel["species"].nunique(), len(morphology), len(iod)) == (24, 1152, 805)


def test_expert_reproductive_strategy_coding_and_provenance() -> None:
    traits = pd.read_csv(
        ROOT / "analyses/03_phylogenetic_path/data/organismal_traits.csv",
        dtype=str,
        keep_default_na=False,
    )
    assert (len(traits), traits["species"].nunique()) == (37, 37)
    indexed = traits.set_index("species")

    groups = {
        strategy: set(group["species"])
        for strategy, group in traits.loc[
            traits["reproductive_strategy"].ne("")
        ].groupby("reproductive_strategy")
    }
    assert groups == {
        "direct_development": {"aeneus", "organi", "wrighti"},
        "aquatic_eggs": {
            "amphileucus",
            "aureatus",
            "gvnigeusgwotli",
            "intermedius",
            "kanawha",
            "marmoratus",
            "mavrokoilius",
        },
        "terrestrial_eggs_aquatic_larvae": {
            "abditus",
            "adatsihi",
            "anicetus",
            "apalachicolae",
            "auriculatus",
            "bairdi",
            "balsameus",
            "campi",
            "carolinensis",
            "catahoula",
            "cheaha",
            "conanti",
            "lycos",
            "monticola",
            "ocoee",
            "orestes",
            "pascagoula",
            "perlapsus",
            "santeetlah",
            "tilleyi",
            "valentinei",
            "valtos",
            "welteri",
        },
    }
    missing = traits.loc[traits["reproductive_strategy"].eq(""), "species"]
    assert set(missing) == {"brimleyorum", "folkertsi", "fuscus", "ochrophaeus"}

    source_id = "pyron_personal_communication_reproductive_strategy"
    coded = traits.loc[traits["reproductive_strategy"].ne("")]
    assert set(coded["reproductive_strategy_source_id"]) == {source_id}
    assert set(coded["reproductive_strategy_confidence"]) == {"high"}
    assert coded["organismal_source_ids"].str.split(";").apply(
        lambda source_ids: source_id in source_ids
    ).all()

    # Egg environment is not adult aquaticity: these deliberate cross-codings
    # prevent the two concepts from being collapsed later.
    assert indexed.loc["amphileucus", "aquaticity_index"] == "1.0"
    assert indexed.loc["amphileucus", "reproductive_strategy"] == "aquatic_eggs"
    assert indexed.loc["welteri", "aquaticity_index"] == "2.0"
    assert (
        indexed.loc["welteri", "reproductive_strategy"]
        == "terrestrial_eggs_aquatic_larvae"
    )

    sources = pd.read_csv(
        ROOT / "data/identity/source_manifest.csv",
        dtype=str,
        keep_default_na=False,
    ).set_index("source_id")
    assert source_id in sources.index
    assert sources.loc[source_id, "year"] == ""
    assert "post hoc" in sources.loc[source_id, "notes"]
    assert "Original email date not yet recorded" in sources.loc[source_id, "notes"]


def test_compact_te_inputs_reproduce_released_pca() -> None:
    result = run_python("analyses/01_transposable_elements/recompute_pca.py")
    assert result.returncode == 0, result.stdout + result.stderr
    assert (
        "reproduce feature prevalence, complete CLR, and final superfamily PCA"
        in result.stdout
    )

    data = ROOT / "analyses/01_transposable_elements/data"
    prevalence = pd.read_csv(data / "te_feature_prevalence.csv")
    clr = pd.read_csv(data / "te_pca_clr_matrix.csv")
    retained = prevalence.loc[prevalence["pca_included"]]
    excluded = prevalence.loc[
        prevalence["te_level"].eq("superfamily") & ~prevalence["pca_included"],
        "feature",
    ]
    assert (len(prevalence), len(clr)) == (39, 1156)
    assert retained.groupby("te_level").size().to_dict() == {
        "order": 10,
        "superfamily": 24,
    }
    assert set(excluded) == {"CR1", "Chapaev", "Dada", "Ginger", "Merlin"}
    assert clr.groupby("te_level").size().to_dict() == {
        "order": 340,
        "superfamily": 816,
    }
    assert clr.columns.tolist() == [
        "species",
        "te_level",
        "composition_mode",
        "feature",
        "closed_proportion",
        "clr_value",
    ]

    makefile = (ROOT / "Makefile").read_text()
    assert "\nte-pca:\n" in makefile
    assert "\nreport: figures validate\n" in makefile
    assert "\nfigure-review: figures validate\n" in makefile
    assert "\nte-pca-view:\n" not in makefile
    assert "\nte-structure:\n" not in makefile


def test_final_pca_structure_audit_rejects_discrete_clusters_but_preserves_signal() -> None:
    data = ROOT / "analyses/01_transposable_elements/data"
    manifest = json.loads((data / "te_pca_clustering_analysis_manifest.json").read_text())
    diagnostics = pd.read_csv(data / "te_pca_clustering_diagnostics.csv")
    stability = pd.read_csv(data / "te_pca_clustering_stability.csv")
    assignments = pd.read_csv(data / "te_pca_candidate_cluster_assignments.csv")
    traits = pd.read_csv(data / "te_pca_trait_association_tests.csv").set_index("test_id")
    phylogeny = pd.read_csv(data / "te_pca_phylogenetic_signal_tests.csv").set_index(
        "test_id"
    )

    assert manifest["analysis_status"] == "exploratory_post_hoc_no_discrete_cluster_accepted"
    assert manifest["selected_k"] is None
    assert manifest["clustering"]["accepted_k"] is None
    assert manifest["clustering"]["gap_primary_rule"] == {
        "distance_power": 2,
        "selector": "Tibs2001SEmax",
        "standard_error_multiplier": 1,
        "selected_k": 1,
    }
    assert manifest["pca_geometry"]["species"] == 34
    assert manifest["pca_geometry"]["ubiquitous_superfamilies"] == 24
    assert manifest["pca_geometry"]["clr_rank"] == 23
    assert manifest["pca_geometry"]["maximum_clr_to_all_pc_distance_difference"] < 1e-12

    expected_inputs = {
        "te34_panel",
        "order_composition",
        "superfamily_composition",
        "scores",
        "variance",
        "loadings",
        "clr",
        "prevalence",
        "traits",
        "tree",
        "pca_recompute_script",
        "analysis_script",
        "environment",
        "python_project",
        "python_lock",
    }
    assert set(manifest["inputs"]) == expected_inputs
    for record in manifest["inputs"].values():
        path = ROOT / record["path"]
        assert path.is_file(), path
        assert hashlib.sha256(path.read_bytes()).hexdigest() == record["sha256"]

    primary = diagnostics.loc[diagnostics["representation"].eq("all_23_pcs")].set_index("k")
    assert set(primary.index) == set(range(1, 9))
    assert not diagnostics["accepted_cluster_solution"].any()
    assert primary.loc[2, "cluster_sizes_ascending"] == "9;25"
    assert np.isclose(primary.loc[2, "mean_silhouette"], 0.1789822345)
    assert primary.loc[7, "cluster_sizes_ascending"] == "1;1;2;4;6;9;11"
    assert primary.loc[7, "singleton_count"] == 2
    assert primary.loc[7, "unconstrained_silhouette_best"]

    stability = stability.set_index(["k", "stability_type"])
    assert set(stability.index.get_level_values("k")) == {2}
    assert np.isclose(
        stability.loc[(2, "leave_one_superfamily_out"), "ari_q10"],
        -0.03658536585,
    )
    assert np.isclose(
        stability.loc[(2, "species_subsample_28_of_34"), "ari_median"],
        -0.04177545692,
    )
    assert np.isclose(stability.loc[(2, "dimension_sensitivity"), "ari_median"].max(), 1)

    assert (len(assignments), assignments[["candidate_id", "species"]].duplicated().sum()) == (
        68,
        0,
    )
    assert set(assignments["candidate_status"]) == {"rejected"}
    k2 = assignments.loc[assignments["candidate_id"].eq("singleton_free_k2")]
    assert set(k2.loc[k2["cluster"].eq("C1"), "species"]) == {
        "amphileucus",
        "aureatus",
        "gvnigeusgwotli",
        "intermedius",
        "kanawha",
        "marmoratus",
        "mavrokoilius",
        "organi",
        "wrighti",
    }

    strategy = traits.loc["reproductive_strategy_permanova"]
    assert strategy["n_species"] == 33
    assert strategy["excluded_species"] == "fuscus"
    assert np.isclose(strategy["variance_explained_r2"], 0.262373, atol=1e-6)
    assert strategy["p_value"] == 0.00001
    assert traits.loc["reproductive_strategy_permdisp", "p_value"] > 0.05
    assert traits.loc["size_adjusted_model_reproductive_strategy", "p_value"] == 0.00001
    assert traits.loc["max_svl_log10_permanova", "n_species"] == 34
    assert traits.loc["max_svl_log10_permanova", "excluded_species"] == "none"
    assert traits.loc["size_adjusted_model_log10_max_svl_mm", "n_species"] == 33
    assert traits.loc["size_adjusted_model_log10_max_svl_mm", "excluded_species"] == "fuscus"
    assert traits.loc["size_adjusted_model_log10_max_svl_mm", "p_value"] > 0.05
    assert traits.loc["size_adjusted_model_reproductive_strategy", "n_species"] == 33
    assert traits.loc["size_adjusted_model_reproductive_strategy", "excluded_species"] == "fuscus"
    assert (
        traits.loc["rejected_k2_reproductive_strategy_ari", "statistic"] > 0.8
    )

    k_mult = phylogeny.loc["k_mult_full_24_feature_clr"]
    assert np.isclose(k_mult["statistic"], 0.6581340270)
    assert k_mult["p_value"] == 0.00001
    assert phylogeny.loc["rejected_k2_fitch_transitions", "statistic"] == 2
    assert phylogeny.loc["first_2_pc_k2_patristic_separation", "p_value"] > 0.05


def test_pairwise_release_uses_relative_iod_scale() -> None:
    traits = pd.read_csv(ROOT / "analyses/03_phylogenetic_path/data/path24_traits.csv")
    pairwise = pd.read_csv(ROOT / "analyses/03_phylogenetic_path/data/pairwise_pgls.csv")
    iod_rows = pairwise.loc[pairwise["left_metric"].eq("relative_iod_index")]
    assert set(iod_rows["comparison"]) == {
        "relative_iod_vs_nucleus_area",
        "relative_iod_vs_cell_area",
    }
    assert (
        not pairwise.astype(str)
        .apply(lambda column: column.str.contains("genome_size", regex=False).any())
        .any()
    )
    np.testing.assert_allclose(
        iod_rows["predictor_log10_mean"],
        np.log10(traits["relative_iod_index"]).mean(),
        atol=5e-9,
    )


def test_only_one_active_environment_specification() -> None:
    assert (ROOT / "environment.yml").exists()
    assert not (ROOT / "Dusky.yml").exists()
    assert not (ROOT / "Publication/figure_environment.yml").exists()
    environment = (ROOT / "environment.yml").read_text()
    assert "r-cluster" in environment
    assert "r-rmarkdown=2.31" in environment
    assert "r-vegan" in environment
