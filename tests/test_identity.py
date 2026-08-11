from __future__ import annotations

import hashlib
import json
import shutil
import subprocess
import sys
from pathlib import Path

import pandas as pd
from Bio import Phylo

ROOT = Path(__file__).resolve().parents[1]
IGNORED_PARTS = {".git", ".pytest_cache", ".ruff_cache", ".venv", "__pycache__"}
IDENTITY_INPUTS = (
    "data/identity/analysis_availability.csv",
    "data/identity/path24_panel.csv",
    "data/identity/source_manifest.csv",
    "data/identity/source_species_aliases.csv",
    "data/identity/species_taxonomy_crosswalk.csv",
    "data/identity/te34_panel.csv",
    "data/identity/trait_registry.csv",
    "analyses/03_phylogenetic_path/data/collaborator_max_svl.csv",
    "analyses/03_phylogenetic_path/data/collaborator_max_svl.xlsx",
    "analyses/03_phylogenetic_path/data/organismal_traits.csv",
    "analyses/03_phylogenetic_path/data/path24_traits.csv",
    "analyses/03_phylogenetic_path/trees/path24_time_tree.nwk",
    "analyses/03_phylogenetic_path/trees/source_time_tree_46.tre",
)


def repository_files() -> set[Path]:
    return {
        path.relative_to(ROOT)
        for path in ROOT.rglob("*")
        if path.is_file() and not IGNORED_PARTS.intersection(path.relative_to(ROOT).parts)
    }


def identity_input_hashes() -> dict[str, str]:
    return {
        relative_path: hashlib.sha256((ROOT / relative_path).read_bytes()).hexdigest()
        for relative_path in IDENTITY_INPUTS
    }


def copy_identity_inputs(destination_root: Path) -> None:
    for relative_path in IDENTITY_INPUTS:
        destination = destination_root / relative_path
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(ROOT / relative_path, destination)


def run_identity() -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [sys.executable, "-m", "scripts.project", "identity"],
        cwd=ROOT,
        check=False,
        capture_output=True,
        text=True,
    )


def audit_errors_for(root: Path) -> tuple[str, ...]:
    code = (
        "import json, sys; "
        "from pathlib import Path; "
        "from scripts.identity import audit_identity; "
        "print(json.dumps(audit_identity(Path(sys.argv[1])).errors))"
    )
    result = subprocess.run(
        [sys.executable, "-c", code, str(root)],
        cwd=ROOT,
        check=False,
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    return tuple(json.loads(result.stdout))


def test_identity_audit_reports_special_resolutions() -> None:
    result = run_identity()
    assert result.returncode == 0, result.stdout + result.stderr
    report = result.stdout.strip()
    assert "orestes-ac=52; orestes-b=50 | exact aliases; max -> orestes=52 mm" in report
    assert "orestes; alternates=orestes_B" in report
    assert "planiceps=48 | source-only row; not transferred to fuscus" in report
    assert "SRX20497025 / GCA_032353935.1" in report
    assert report.endswith(
        "PASS: source aliases, Path24 tree representatives, SVL maxima, "
        "and accession-scoped fuscus identity agree"
    )


def test_svl_aliases_are_exact_and_aggregate_by_registered_maximum() -> None:
    aliases = pd.read_csv(
        ROOT / "data/identity/source_species_aliases.csv",
        dtype=str,
        keep_default_na=False,
    )
    assert set(map(tuple, aliases[["source_id", "source_label", "canonical_species"]].values)) == {
        ("collaborator_desmog_svl", "orestes-ac", "orestes"),
        ("collaborator_desmog_svl", "orestes-b", "orestes"),
    }
    assert not aliases["source_id"].eq("repo_desmognathus_phylogeny_tree").any()
    assert not aliases["source_label"].eq("planiceps").any()

    registry = pd.read_csv(
        ROOT / "data/identity/trait_registry.csv",
        dtype=str,
        keep_default_na=False,
    ).set_index("trait_name")
    assert registry.loc["max_svl_mm", "canonical_species_aggregation"] == "max"

    source = pd.read_csv(
        ROOT / "analyses/03_phylogenetic_path/data/collaborator_max_svl.csv"
    ).set_index("source_taxon")
    traits = pd.read_csv(
        ROOT / "analyses/03_phylogenetic_path/data/organismal_traits.csv"
    ).set_index("species")
    assert source.loc[["orestes-ac", "orestes-b"], "max_svl_mm"].max() == 52
    assert traits.loc["orestes", "max_svl_mm"] == 52
    assert traits.loc["orestes", "max_svl_source_taxon"] == "orestes-ac;orestes-b"


def test_tree_uses_exact_representatives_without_suffix_stripping() -> None:
    source = Phylo.read(
        ROOT / "analyses/03_phylogenetic_path/trees/source_time_tree_46.tre",
        "newick",
    )
    released = Phylo.read(
        ROOT / "analyses/03_phylogenetic_path/trees/path24_time_tree.nwk",
        "newick",
    )
    source_tips = {str(tip.name) for tip in source.get_terminals()}
    released_tips = {str(tip.name) for tip in released.get_terminals()}

    assert {"fuscus", "fuscus_A", "fuscus_E", "orestes", "orestes_B", "planiceps"} <= (
        source_tips
    )
    assert {"fuscus", "orestes"} <= released_tips
    assert {"fuscus_A", "fuscus_E", "orestes_B", "planiceps"}.isdisjoint(released_tips)


def test_identity_audit_rejects_species_to_tip_swaps(tmp_path: Path) -> None:
    copy_identity_inputs(tmp_path)

    traits_path = tmp_path / "analyses/03_phylogenetic_path/data/path24_traits.csv"
    traits = pd.read_csv(traits_path)
    fuscus_tip = traits.loc[traits["species"].eq("D. fuscus"), "tree_tip"].iloc[0]
    orestes_tip = traits.loc[traits["species"].eq("D. orestes"), "tree_tip"].iloc[0]
    traits.loc[traits["species"].eq("D. fuscus"), "tree_tip"] = orestes_tip
    traits.loc[traits["species"].eq("D. orestes"), "tree_tip"] = fuscus_tip
    traits.to_csv(traits_path, index=False)

    errors = audit_errors_for(tmp_path)
    mismatches = [error for error in errors if error.startswith("[TREE_SPECIES_TIP_MISMATCH]")]
    assert len(mismatches) == 2


def test_identity_audit_binds_aliases_to_raw_labels_and_decision(tmp_path: Path) -> None:
    copy_identity_inputs(tmp_path)
    aliases_path = tmp_path / "data/identity/source_species_aliases.csv"
    aliases = pd.read_csv(aliases_path)
    aliases.loc[aliases["source_label"].eq("orestes-ac"), "decision_basis_source_id"] = (
        "collaborator_desmog_svl"
    )
    aliases.loc[aliases["source_label"].eq("orestes-ac"), "decision_date"] = "1900-01-01"
    aliases.loc[len(aliases)] = {
        "source_id": "collaborator_desmog_svl",
        "source_label": "orestes-c",
        "canonical_species": "orestes",
        "decision_basis_source_id": "author_curator_decision_2026_08_10_identity_resolution",
        "decision_date": "2026-08-10",
        "notes": "deliberate invalid test row",
    }
    aliases.to_csv(aliases_path, index=False)

    errors = audit_errors_for(tmp_path)
    assert any(error.startswith("[IDENTITY_ALIAS_MISSING]") for error in errors)
    assert any(error.startswith("[IDENTITY_ALIAS_UNUSED]") for error in errors)


def test_identity_audit_reports_independent_structural_failures(tmp_path: Path) -> None:
    copy_identity_inputs(tmp_path)

    svl_path = tmp_path / "analyses/03_phylogenetic_path/data/collaborator_max_svl.csv"
    svl = pd.read_csv(svl_path)
    svl = pd.concat([svl, svl.loc[svl["source_taxon"].eq("orestes-ac")]], ignore_index=True)
    svl.to_csv(svl_path, index=False)

    panel_path = tmp_path / "data/identity/path24_panel.csv"
    panel = pd.read_csv(panel_path)
    panel = pd.concat([panel, panel.iloc[[0]]], ignore_index=True)
    panel.to_csv(panel_path, index=False)

    manifest_path = tmp_path / "data/identity/source_manifest.csv"
    manifest = pd.read_csv(manifest_path, dtype=str, keep_default_na=False)
    extra = manifest.loc[
        manifest["source_id"].eq("genomic_fuscus_srx20497025_gca032353935_1")
    ].iloc[0].copy()
    extra["source_id"] = "invalid_second_planiceps_to_fuscus"
    extra["sra_accessions"] = "SRX00000000"
    extra["assembly_accession"] = "GCA_000000000.1"
    manifest = pd.concat([manifest, extra.to_frame().T], ignore_index=True)
    manifest.to_csv(manifest_path, index=False)

    errors = audit_errors_for(tmp_path)
    assert any(error.startswith("[SVL_SOURCE_DUPLICATE]") for error in errors)
    assert any(error.startswith("[TREE_PANEL_CONFLICT]") for error in errors)
    assert any(error.startswith("[FUSCUS_SCOPE_BREACH]") for error in errors)


def test_identity_audit_rejects_unregistered_tree_variant_and_fuscus_drift(
    tmp_path: Path,
) -> None:
    copy_identity_inputs(tmp_path)

    tree_path = tmp_path / "analyses/03_phylogenetic_path/trees/source_time_tree_46.tre"
    tree_text = tree_path.read_text()
    assert "abditus:" in tree_text
    tree_path.write_text(tree_text.replace("abditus:", "fuscus_Z:", 1))

    manifest_path = tmp_path / "data/identity/source_manifest.csv"
    manifest = pd.read_csv(manifest_path, dtype=str, keep_default_na=False)
    active = manifest["source_id"].eq("genomic_fuscus_srx20497025_gca032353935_1")
    manifest.loc[active, "analysis_role"] = "invalid_role"
    manifest.to_csv(manifest_path, index=False)

    errors = audit_errors_for(tmp_path)
    assert any(error.startswith("[TREE_TE_TIP_MISSING]") for error in errors)
    assert any(error.startswith("[TREE_UNREGISTERED_VARIANT]") for error in errors)
    assert any(error.startswith("[FUSCUS_RESOURCE_MISMATCH]") for error in errors)


def test_identity_command_is_read_only() -> None:
    before = repository_files()
    hashes_before = identity_input_hashes()
    result = run_identity()
    after = repository_files()
    hashes_after = identity_input_hashes()

    assert result.returncode == 0, result.stdout + result.stderr
    assert "Identity resolution (read-only)" in result.stdout
    assert before == after
    assert hashes_before == hashes_after
