from __future__ import annotations

import csv
import hashlib
import importlib.util
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
TE_DATA = ROOT / "analyses/01_transposable_elements/data"
SOURCE_SNAPSHOT = TE_DATA / "wang2025_figure2_te_diversity.csv"
SOURCE_PROVENANCE = TE_DATA / "wang2025_figure2_te_diversity.provenance.json"
COMBINED = TE_DATA / "te_diversity_cross_study.csv"
DIVERSITY = TE_DATA / "te_diversity.csv"
SUPERFAMILY_COMPOSITION = TE_DATA / "te_superfamily_composition.csv"
IOD_SPECIES = ROOT / "analyses/02_morphology/data/relative_nuclear_iod_species.csv"
TE34_PANEL = ROOT / "data/identity/te34_panel.csv"
PATH24_PANEL = ROOT / "data/identity/path24_panel.csv"

PUBLICATION_MANIFEST = ROOT / "Publication/datasets/DATASET_MANIFEST.csv"
NUMBERED_FIGURE_MANIFEST = ROOT / "Publication/figures/FIGURE_MANIFEST.csv"
AUXILIARY_FIGURE_MANIFEST = ROOT / "Publication/figures/AUXILIARY_FIGURE_MANIFEST.csv"
FIGURE_SOURCE = ROOT / "scripts/publication/build_gbe_figures.R"
CANONICAL_DIVERSITY_SOURCE = ROOT / "analyses/01_transposable_elements/recompute_pca.py"
CROSS_STUDY_BUILDER_SOURCE = (
    ROOT / "analyses/01_transposable_elements/build_cross_study_diversity.py"
)


def _load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot import {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


CANONICAL_DIVERSITY = _load_module(CANONICAL_DIVERSITY_SOURCE, "canonical_te_diversity")
CROSS_STUDY_BUILDER = _load_module(CROSS_STUDY_BUILDER_SOURCE, "cross_study_builder")

DATASET_DOI = "10.5061/dryad.zpc866tkv"
ARTICLE_DOI = "10.1038/s42003-025-08127-3"
DRYAD_ARCHIVE_SHA256 = "aff6c3fdcd566a5ec1622ce3b256a38e9c337520a812ef75005a8fd896671dc9"
DRYAD_WORKBOOK_SHA256 = "367a8102fd18f9c17251c442b0a0112af6afc008b51bf4b110f72a98fcc05c70"
SUPERFAMILY_COMPOSITION_SHA256 = "8b3b33eef684a53b474ed00352d924a946ed6820beaa47aaed9699423755c4e1"
DIVERSITY_SHA256 = "3df4c2efc5ec728edf1373c463163850be92d1900c55026dffce346e83489a2a"
IOD_SPECIES_SHA256 = "327542c76109e8a31f26eea88edaf30eb9d756faf6c0f24ebb78f5d3bc3e6699"

SOURCE_COLUMNS = [
    "source_row",
    "english_name",
    "species",
    "family",
    "genome_size_gb",
    "shannon_top10_superfamilies",
    "gini_simpson_top10_superfamilies",
    "vertebrate_class",
    "figure_group",
    "figure_panel",
    "fish_point_shape",
]
COMBINED_COLUMNS = [
    "data_origin",
    "source_record_id",
    "source_doi",
    "source_file",
    "source_file_sha256",
    "metric_scope",
    "scientific_name",
    "common_name",
    "family",
    "comparison_group",
    "genome_size_gb",
    "genome_size_ci_low_gb",
    "genome_size_ci_high_gb",
    "shannon_top10_superfamilies",
    "gini_simpson_top10_superfamilies",
    "figure_inclusion",
]

PUBLISHED_ORIGIN = "wang2025_published"
CURRENT_ORIGIN = "desmognathus_current"
PUBLISHED_SCOPE = "published_top10_superfamilies_repeatmodeler_repeatmasker"
CURRENT_SCOPE = "current_top10_superfamilies_dnapipete_retry_averaged"

AUXILIARY_ID = "Auxiliary_Wang2025_cross_study_te_diversity_reference"
AUXILIARY_BASENAME = "Wang2025_cross_study_TE_diversity_reference"
EXPECTED_AUXILIARY_IDS = [
    "Auxiliary_Wang_style_integrated_te_history_candidate",
    "Auxiliary_Wang_style_diversity_genome_size_candidate",
    "Auxiliary_Path24_development_traits_reference",
    "Auxiliary_Path_DAG_reference",
    AUXILIARY_ID,
]


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _manifest_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def _read_csv_exact(path: Path) -> pd.DataFrame:
    # The Dryad workbook contains more than 15 significant digits for many
    # indices. Pandas' default C-parser conversion can shift those by one ULP.
    return pd.read_csv(path, float_precision="round_trip")


def _semantic_digest(frame: pd.DataFrame, columns: list[str], sort_by: str) -> str:
    records: list[list[object]] = []
    for row in frame.sort_values(sort_by, kind="stable")[columns].itertuples(
        index=False, name=None
    ):
        values: list[object] = []
        for value in row:
            if pd.isna(value):
                values.append(None)
            elif isinstance(value, (float, np.floating)):
                values.append(float(value).hex())
            elif isinstance(value, (int, np.integer)):
                values.append(int(value))
            else:
                values.append(str(value))
        records.append(values)
    payload = json.dumps(records, ensure_ascii=False, separators=(",", ":")).encode()
    return hashlib.sha256(payload).hexdigest()


def _current_top10_metrics() -> pd.DataFrame:
    composition = _read_csv_exact(SUPERFAMILY_COMPOSITION).set_index("species")
    return CANONICAL_DIVERSITY.build_diversity(composition).rename(
        columns={
            "shannon_entropy": "shannon_expected",
            "gini_simpson": "gini_simpson_expected",
        }
    )


def test_cross_study_current_metrics_use_the_canonical_diversity_authority() -> None:
    assert not hasattr(CROSS_STUDY_BUILDER, "top10_metrics")
    assert not hasattr(CROSS_STUDY_BUILDER, "COMPOSITION_INPUT")

    expected = _read_csv_exact(DIVERSITY).set_index("species")
    observed = (
        CROSS_STUDY_BUILDER.build_current_rows()
        .assign(species=lambda frame: frame["scientific_name"].str.removeprefix("Desmognathus "))
        .set_index("species")
    )

    assert observed.index.tolist() == expected.index.tolist()
    assert expected.shape == (34, 3)
    assert expected["n_superfamilies_included"].eq(10).all()
    np.testing.assert_array_equal(
        observed["shannon_top10_superfamilies"].to_numpy(),
        expected["shannon_entropy"].to_numpy(),
    )
    np.testing.assert_array_equal(
        observed["gini_simpson_top10_superfamilies"].to_numpy(),
        expected["gini_simpson"].to_numpy(),
    )

    oracle = _current_top10_metrics().set_index("species")
    np.testing.assert_allclose(expected["shannon_entropy"], oracle["shannon_expected"], atol=1e-12)
    np.testing.assert_allclose(
        expected["gini_simpson"],
        oracle["gini_simpson_expected"],
        atol=1e-12,
    )


def test_wang2025_snapshot_is_the_exact_84_row_dryad_figure2_source() -> None:
    source = _read_csv_exact(SOURCE_SNAPSHOT)

    assert source.shape == (84, 11)
    assert list(source.columns) == SOURCE_COLUMNS
    assert source["source_row"].tolist() == list(range(1, 85))
    assert source["species"].nunique() == 84
    assert (
        not source[
            [
                "english_name",
                "species",
                "family",
                "genome_size_gb",
                "shannon_top10_superfamilies",
                "gini_simpson_top10_superfamilies",
                "vertebrate_class",
                "figure_group",
                "figure_panel",
            ]
        ]
        .isna()
        .any()
        .any()
    )

    assert source.groupby("figure_group").size().to_dict() == {
        "Fishes": 19,
        "Anura": 10,
        "Caudata": 11,
        "Non-avian Reptillia": 16,
        "Aves": 14,
        "Mammalia": 14,
    }
    assert source.groupby("vertebrate_class").size().to_dict() == {
        "Fish": 19,
        "Amphibian": 21,
        "Reptile": 16,
        "Bird": 14,
        "Mammal": 14,
    }
    assert source["figure_panel"].value_counts().to_dict() == {
        "A": 19,
        "B": 10,
        "C": 11,
        "D": 16,
        "E": 14,
        "F": 14,
    }
    shapes = source["fish_point_shape"].dropna().astype(int).value_counts().to_dict()
    assert shapes == {8: 3, 17: 3}

    species = set(source["species"])
    assert "Desmognathus ochrophaeus" in species
    assert "Anas platyrhynchos" in species
    assert "Anas anas" not in species
    assert species.isdisjoint(
        {
            "Geotrypetes seraphini",
            "Ichthyophis bannanicus",
            "Microcaecilia unicolor",
            "Rhinatrema bivittatum",
        }
    )

    # These digests bind every source identity and every published numeric value,
    # while remaining insensitive to CSV byte formatting.
    source_groups = source.copy()
    source_groups["fish_point_shape"] = source_groups["fish_point_shape"].astype("Int64")
    assert (
        _semantic_digest(
            source_groups,
            ["source_row", "english_name", "species", "family"],
            "source_row",
        )
        == "570e43073362ac73eeb94af8d47b4058872982e640aa9ec3ce63ff1b2bab8dfa"
    )
    assert (
        _semantic_digest(
            source,
            [
                "species",
                "genome_size_gb",
                "shannon_top10_superfamilies",
                "gini_simpson_top10_superfamilies",
            ],
            "species",
        )
        == "e6132916b6f903a9b73aa702e98318c73b7b10ae5f67dbf99ea6790b6804d492"
    )
    assert (
        _semantic_digest(
            source_groups,
            [
                "source_row",
                "vertebrate_class",
                "figure_group",
                "figure_panel",
                "fish_point_shape",
            ],
            "source_row",
        )
        == "5c5a1250e8631b087c389d8b2767acbf45693c8201bf996ff173cc7c0604caaa"
    )

    focal = source.set_index("species")
    assert focal.loc["Desmognathus ochrophaeus", "genome_size_gb"] == 15.0
    assert focal.loc["Desmognathus ochrophaeus", "shannon_top10_superfamilies"] == 1.61
    assert focal.loc["Desmognathus ochrophaeus", "gini_simpson_top10_superfamilies"] == 0.71
    # The machine-readable Dryad workbook, not rounded Supplementary Table S4,
    # is the admitted numerical authority.
    assert focal.loc["Ranodon sibiricus", "genome_size_gb"] == 21.0
    assert focal.loc["Pachytriton brevipes", "genome_size_gb"] == 40.0


def test_wang2025_provenance_binds_dois_hashes_license_and_acceptance_boundary() -> None:
    provenance = json.loads(SOURCE_PROVENANCE.read_text(encoding="utf-8"))
    required_keys = {
        "artifact_id",
        "source_dataset_doi",
        "source_article_doi",
        "source_dataset_url",
        "source_archive_name",
        "source_archive_file_id",
        "source_archive_size_bytes",
        "source_archive_sha256",
        "source_member_path",
        "source_member_sha256",
        "source_rows",
        "source_columns",
        "license",
        "retrieved_on",
        "normalization",
        "metric_contract",
        "acceptance_decision",
    }
    assert required_keys <= set(provenance)
    assert provenance["source_dataset_doi"] == DATASET_DOI
    assert provenance["source_article_doi"] == ARTICLE_DOI
    assert DATASET_DOI in provenance["source_dataset_url"]
    assert provenance["source_archive_name"] == "Fig_2.zip"
    assert int(provenance["source_archive_file_id"]) == 4030762
    assert int(provenance["source_archive_size_bytes"]) == 305_298
    assert provenance["source_archive_sha256"] == DRYAD_ARCHIVE_SHA256
    assert provenance["source_member_path"] == "Fig_2/Data_TE_diversity.xlsx"
    assert provenance["source_member_sha256"] == DRYAD_WORKBOOK_SHA256
    assert (int(provenance["source_rows"]), int(provenance["source_columns"])) == (84, 11)
    assert provenance["license"] in {
        "CC0-1.0",
        "https://spdx.org/licenses/CC0-1.0.html",
    }
    assert provenance["retrieved_on"] == "2026-08-12"

    normalization = json.dumps(provenance["normalization"], ensure_ascii=False).casefold()
    metric_contract = json.dumps(provenance["metric_contract"], ensure_ascii=False).casefold()
    acceptance = json.dumps(provenance["acceptance_decision"], ensure_ascii=False).casefold()
    assert "whitespace" in normalization
    assert "round" not in normalization or "no rounding" in normalization
    for phrase in ("natural", "top 10", "gini-simpson", "reclos", "unresolved"):
        assert phrase in metric_contract
    assert "84" in acceptance and "dryad" in acceptance
    assert "87" in acceptance and ("reject" in acceptance or "not" in acceptance)


def test_s41_combines_84_exact_published_and_34_recomputed_current_rows() -> None:
    source = _read_csv_exact(SOURCE_SNAPSHOT)
    combined = _read_csv_exact(COMBINED)

    assert combined.shape == (118, 16)
    assert list(combined.columns) == COMBINED_COLUMNS
    assert combined["data_origin"].value_counts().to_dict() == {
        PUBLISHED_ORIGIN: 84,
        CURRENT_ORIGIN: 34,
    }
    assert combined["source_record_id"].nunique() == 118
    assert combined["figure_inclusion"].value_counts().to_dict() == {
        "published_scatter_and_group": 11,
        "published_group_only": 73,
        "current_scatter_and_group": 21,
        "current_group_only": 13,
    }

    published = combined.loc[combined["data_origin"].eq(PUBLISHED_ORIGIN)].copy()
    current = combined.loc[combined["data_origin"].eq(CURRENT_ORIGIN)].copy()
    assert set(published["source_doi"]) == {DATASET_DOI}
    assert set(published["source_file"]) == {"Fig_2/Data_TE_diversity.xlsx"}
    assert set(published["source_file_sha256"]) == {DRYAD_WORKBOOK_SHA256}
    assert set(published["metric_scope"]) == {PUBLISHED_SCOPE}
    assert set(current["source_doi"]) == {"this_study"}
    assert set(current["metric_scope"]) == {CURRENT_SCOPE}

    expected_published = source.rename(
        columns={
            "species": "scientific_name",
            "english_name": "common_name",
        }
    )[
        [
            "scientific_name",
            "common_name",
            "family",
            "figure_group",
            "fish_point_shape",
            "genome_size_gb",
            "shannon_top10_superfamilies",
            "gini_simpson_top10_superfamilies",
        ]
    ].copy()
    expected_published["comparison_group"] = expected_published["figure_group"].replace(
        {
            "Anura": "Frogs",
            "Caudata": "Published salamanders",
            "Non-avian Reptillia": "Non-avian reptiles",
            "Aves": "Birds",
            "Mammalia": "Mammals",
        }
    )
    fishes = expected_published["figure_group"].eq("Fishes")
    expected_published.loc[
        fishes & expected_published["fish_point_shape"].isna(), "comparison_group"
    ] = "Bony fishes"
    expected_published.loc[
        fishes & expected_published["fish_point_shape"].eq(17), "comparison_group"
    ] = "Cartilaginous fishes"
    expected_published.loc[
        fishes & expected_published["fish_point_shape"].eq(8), "comparison_group"
    ] = "Lobe-finned fishes"
    expected_published = expected_published.drop(columns=["figure_group", "fish_point_shape"])
    observed_published = published[expected_published.columns]
    pd.testing.assert_frame_equal(
        observed_published.sort_values("scientific_name").reset_index(drop=True),
        expected_published.sort_values("scientific_name").reset_index(drop=True),
        check_exact=True,
    )
    published_caudata = set(
        published.loc[
            published["figure_inclusion"].eq("published_scatter_and_group"),
            "scientific_name",
        ]
    )
    assert published_caudata == set(source.loc[source["figure_group"].eq("Caudata"), "species"])
    published_desmognathus = published.loc[
        published["scientific_name"].str.startswith("Desmognathus "),
        "scientific_name",
    ].tolist()
    assert published_desmognathus == ["Desmognathus ochrophaeus"]
    assert "Desmognathus ochrophaeus" not in set(current["scientific_name"])
    assert published["comparison_group"].value_counts().to_dict() == {
        "Bony fishes": 13,
        "Cartilaginous fishes": 3,
        "Lobe-finned fishes": 3,
        "Frogs": 10,
        "Published salamanders": 11,
        "Non-avian reptiles": 16,
        "Birds": 14,
        "Mammals": 14,
    }

    expected_metrics = _current_top10_metrics()
    observed_metrics = current.assign(
        species=current["scientific_name"].str.removeprefix("Desmognathus ")
    ).merge(expected_metrics, on="species", how="outer", validate="one_to_one", indicator=True)
    assert observed_metrics["_merge"].eq("both").all()
    np.testing.assert_allclose(
        observed_metrics["shannon_top10_superfamilies"],
        observed_metrics["shannon_expected"],
        rtol=0,
        atol=1e-12,
    )
    np.testing.assert_allclose(
        observed_metrics["gini_simpson_top10_superfamilies"],
        observed_metrics["gini_simpson_expected"],
        rtol=0,
        atol=1e-12,
    )
    assert set(current["comparison_group"]) == {"Current Desmognathus"}
    assert set(current["family"]) == {"Plethodontidae"}


def test_current_scatter_intervals_are_the_exact_te34_path24_overlap() -> None:
    combined = _read_csv_exact(COMBINED)
    current = combined.loc[combined["data_origin"].eq(CURRENT_ORIGIN)].copy()
    current["species"] = current["scientific_name"].str.removeprefix("Desmognathus ")

    te34 = set(_read_csv_exact(TE34_PANEL)["species"])
    path24 = set(_read_csv_exact(PATH24_PANEL)["species"])
    overlap = te34 & path24
    scatter = current.loc[current["figure_inclusion"].eq("current_scatter_and_group")]
    group_only = current.loc[current["figure_inclusion"].eq("current_group_only")]
    assert set(scatter["species"]) == overlap
    assert set(group_only["species"]) == te34 - overlap
    assert (len(scatter), len(group_only)) == (21, 13)

    iod = _read_csv_exact(IOD_SPECIES).assign(
        species=lambda frame: frame["species"].str.removeprefix("D. ")
    )
    expected = iod.loc[iod["species"].isin(overlap)].assign(
        genome_size_gb_expected=lambda frame: frame["genome_size_pg_fuscus_anchored"] * 0.978,
        genome_size_ci_low_gb_expected=lambda frame: frame["genome_size_pg_ci_low"] * 0.978,
        genome_size_ci_high_gb_expected=lambda frame: frame["genome_size_pg_ci_high"] * 0.978,
    )
    checked = scatter.merge(
        expected[
            [
                "species",
                "genome_size_gb_expected",
                "genome_size_ci_low_gb_expected",
                "genome_size_ci_high_gb_expected",
            ]
        ],
        on="species",
        validate="one_to_one",
    )
    for observed, expected_column in (
        ("genome_size_gb", "genome_size_gb_expected"),
        ("genome_size_ci_low_gb", "genome_size_ci_low_gb_expected"),
        ("genome_size_ci_high_gb", "genome_size_ci_high_gb_expected"),
    ):
        np.testing.assert_allclose(checked[observed], checked[expected_column], rtol=0, atol=1e-10)
    assert checked["genome_size_ci_low_gb"].le(checked["genome_size_gb"]).all()
    assert checked["genome_size_gb"].le(checked["genome_size_ci_high_gb"]).all()
    assert (
        group_only[["genome_size_gb", "genome_size_ci_low_gb", "genome_size_ci_high_gb"]]
        .isna()
        .all()
        .all()
    )

    for values in current["source_file"]:
        assert "analyses/01_transposable_elements/data/te_diversity.csv" in values
        assert "te_superfamily_composition.csv" not in values
        assert "analyses/02_morphology/data/relative_nuclear_iod_species.csv" in values
    for values in current["source_file_sha256"]:
        assert DIVERSITY_SHA256 in values
        assert SUPERFAMILY_COMPOSITION_SHA256 not in values
        assert IOD_SPECIES_SHA256 in values
    assert _sha256(DIVERSITY) == DIVERSITY_SHA256
    assert _sha256(SUPERFAMILY_COMPOSITION) == SUPERFAMILY_COMPOSITION_SHA256
    assert _sha256(IOD_SPECIES) == IOD_SPECIES_SHA256


def test_s41_release_is_manifest_bound_and_contains_no_pooled_statistics() -> None:
    manifest = pd.read_csv(PUBLICATION_MANIFEST)
    assert manifest["supplement_id"].tolist()[-1] == "S43"
    assert "S41" in manifest["supplement_id"].tolist()
    row = manifest.loc[manifest["supplement_id"].eq("S41")]
    assert len(row) == 1
    row = row.iloc[0]
    assert (
        row["source_path"] == "analyses/01_transposable_elements/data/te_diversity_cross_study.csv"
    )
    assert (int(row["rows"]), int(row["columns"])) == (118, 16)
    publication_table = ROOT / "Publication/datasets" / str(row["filename"])
    assert publication_table.read_bytes() == COMBINED.read_bytes()
    assert row["source_sha256"] == _sha256(COMBINED)
    assert row["output_sha256"] == _sha256(publication_table)

    combined = _read_csv_exact(COMBINED)
    assert not any(
        re.search(
            r"(^|_)(p(_?value)?|q(_?value)?|statistic|slope|intercept|r2|rho|tau)($|_)", column
        )
        for column in combined.columns
    )
    release_text = " ".join(
        str(row[field]) for field in ("title", "description", "manuscript_use", "release_status")
    ).casefold()
    for phrase in ("descriptive", "no pooled", "different"):
        assert phrase in release_text


def test_cross_study_reference_is_the_final_unnumbered_auxiliary_without_a_fit() -> None:
    numbered = _manifest_rows(NUMBERED_FIGURE_MANIFEST)
    auxiliary = _manifest_rows(AUXILIARY_FIGURE_MANIFEST)
    assert all(row["basename"] != AUXILIARY_BASENAME for row in numbered)
    assert [row["auxiliary_id"] for row in auxiliary] == EXPECTED_AUXILIARY_IDS

    row = auxiliary[-1]
    assert row["basename"] == AUXILIARY_BASENAME
    assert row["datasets"] == "S41"
    evidence = " ".join(
        row[field] for field in ("title", "release_status", "boundary", "legend", "alt_text")
    ).casefold()
    for phrase in (
        "descriptive",
        "no pooled",
        "repeatmodeler",
        "repeatmasker",
        "dnapipete",
        "different",
    ):
        assert phrase in evidence

    source = FIGURE_SOURCE.read_text(encoding="utf-8")
    basename_position = source.index(AUXILIARY_BASENAME)
    block_start = source.rfind("## ----", 0, basename_position)
    block_end = source.find("## ----", basename_position)
    block = source[block_start : None if block_end < 0 else block_end]
    assert "record_auxiliary_figure(" in block
    assert "record_figure(" not in block
    assert "S41" in block
    assert not re.search(
        r"\b(?:stat_smooth|geom_smooth|lm|glm|gls|pgls|cor\.test)\s*\(",
        block,
    )

    for field in ("pdf", "png", "tiff"):
        asset = ROOT / "Publication/figures" / row[field]
        assert asset.is_file(), asset
        assert asset.stat().st_size > 10_000, asset
