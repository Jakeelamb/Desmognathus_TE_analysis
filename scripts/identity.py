"""Read-only audit of source-label, trait, tree-tip, and accession identity rules."""

from __future__ import annotations

import copy
import math
from dataclasses import dataclass
from itertools import combinations
from pathlib import Path

import pandas as pd
from Bio import Phylo

_ALIAS_COLUMNS = (
    "source_id",
    "source_label",
    "canonical_species",
    "decision_basis_source_id",
    "decision_date",
    "notes",
)
_SVL_SOURCE_ID = "collaborator_desmog_svl"
_SVL_SOURCE_PATH = "analyses/03_phylogenetic_path/data/collaborator_max_svl.xlsx"
_IDENTITY_DECISION_ID = "author_curator_decision_2026_08_10_identity_resolution"
_EXPECTED_UNUSED_SVL_LABELS = frozenset({"imitator", "planiceps"})
_TREE_ALTERNATES = {
    "fuscus": ("fuscus_A", "fuscus_E"),
    "orestes": ("orestes_B",),
}
_TREE_DISTANCE_TOLERANCE = 2e-5


@dataclass(frozen=True)
class IdentityRow:
    """One human-readable identity decision checked by the audit."""

    stream: str
    source: str
    resolution: str


@dataclass(frozen=True)
class IdentityAudit:
    """Deterministic identity report; an empty error tuple means the contract holds."""

    rows: tuple[IdentityRow, ...]
    errors: tuple[str, ...]

    def render(self) -> str:
        lines = ["Identity resolution (read-only)", ""]
        lines.extend(
            f"- {row.stream} | {row.source} | {row.resolution}" for row in self.rows
        )
        lines.append("")
        if self.errors:
            lines.append(f"FAIL: {len(self.errors)} identity contract violation(s)")
            lines.extend(f"- {error}" for error in self.errors)
        else:
            lines.append(
                "PASS: source aliases, Path24 tree representatives, SVL maxima, "
                "and accession-scoped fuscus identity agree"
            )
        return "\n".join(lines)


def _read_csv(root: Path, relative_path: str) -> pd.DataFrame:
    return pd.read_csv(root / relative_path, dtype=str, keep_default_na=False)


def _require_columns(
    frame: pd.DataFrame,
    required: tuple[str, ...],
    label: str,
    errors: list[str],
) -> bool:
    missing = [column for column in required if column not in frame.columns]
    if missing:
        errors.append(f"[IDENTITY_SCHEMA] {label} missing columns: {', '.join(missing)}")
        return False
    return True


def _load_aliases(
    root: Path,
    sources: pd.DataFrame,
    taxonomy: pd.DataFrame,
    errors: list[str],
) -> tuple[pd.DataFrame, dict[tuple[str, str], str]]:
    aliases = _read_csv(root, "data/identity/source_species_aliases.csv")
    if tuple(aliases.columns) != _ALIAS_COLUMNS:
        errors.append(
            "[IDENTITY_SCHEMA] source_species_aliases.csv columns must be exactly: "
            + ", ".join(_ALIAS_COLUMNS)
        )
    if not _require_columns(aliases, _ALIAS_COLUMNS, "source species aliases", errors):
        return aliases, {}

    source_ids = set(sources["source_id"])
    canonical_species = set(taxonomy["current_species"])
    mapping: dict[tuple[str, str], str] = {}
    for row in aliases.to_dict(orient="records"):
        key = (row["source_id"], row["source_label"])
        fields = (
            row["source_id"],
            row["source_label"],
            row["canonical_species"],
            row["decision_basis_source_id"],
            row["decision_date"],
        )
        if any(not value or value != value.strip() for value in fields):
            errors.append(f"[IDENTITY_ALIAS_INVALID] alias fields must be nonblank and exact: {key}")
            continue
        if key in mapping:
            errors.append(f"[IDENTITY_ALIAS_CONFLICT] duplicate source alias: {key}")
            continue
        if row["source_id"] not in source_ids:
            errors.append(f"[IDENTITY_ALIAS_SOURCE_MISSING] unknown source_id: {row['source_id']}")
        if row["canonical_species"] not in canonical_species:
            errors.append(
                "[IDENTITY_CANONICAL_MISSING] unknown canonical species: "
                f"{row['canonical_species']}"
            )
        if row["decision_basis_source_id"] not in source_ids:
            errors.append(
                "[IDENTITY_DECISION_MISSING] unknown decision source: "
                f"{row['decision_basis_source_id']}"
            )
        if row["source_id"] == "repo_desmognathus_phylogeny_tree":
            errors.append(
                "[TREE_ALIAS_FORBIDDEN] tree representatives must be selected exactly, "
                f"not aliased: {row['source_label']}"
            )
        elif row["source_id"] != _SVL_SOURCE_ID:
            errors.append(
                "[IDENTITY_ALIAS_UNSUPPORTED] the active resolver supports only exact "
                f"collaborator-SVL aliases, observed source_id={row['source_id']}"
            )
        mapping[key] = row["canonical_species"]

    required = {
        (_SVL_SOURCE_ID, "orestes-ac"): {
            "canonical_species": "orestes",
            "decision_basis_source_id": _IDENTITY_DECISION_ID,
            "decision_date": "2026-08-10",
        },
        (_SVL_SOURCE_ID, "orestes-b"): {
            "canonical_species": "orestes",
            "decision_basis_source_id": _IDENTITY_DECISION_ID,
            "decision_date": "2026-08-10",
        },
    }
    for key, expected_fields in required.items():
        match = aliases.loc[
            aliases["source_id"].eq(key[0]) & aliases["source_label"].eq(key[1])
        ]
        if len(match) != 1 or any(
            match.iloc[0][column] != value for column, value in expected_fields.items()
        ):
            errors.append(
                "[IDENTITY_ALIAS_MISSING] required exact alias and decision provenance "
                f"are missing or changed: {key}"
            )
    return aliases, mapping


def _audit_svl(
    root: Path,
    taxonomy: pd.DataFrame,
    alias_mapping: dict[tuple[str, str], str],
    rows: list[IdentityRow],
    errors: list[str],
) -> None:
    registry = _read_csv(root, "data/identity/trait_registry.csv")
    if _require_columns(
        registry,
        ("trait_name", "canonical_species_aggregation"),
        "trait registry",
        errors,
    ):
        if registry["trait_name"].duplicated().any():
            errors.append("[TRAIT_RULE_CONFLICT] trait_registry.csv has duplicate trait names")
        rule_rows = registry.loc[registry["trait_name"].eq("max_svl_mm")]
        if len(rule_rows) != 1 or rule_rows.iloc[0]["canonical_species_aggregation"] != "max":
            errors.append("[TRAIT_RULE_MISSING] max_svl_mm must aggregate canonical aliases by max")

    source = _read_csv(
        root,
        "analyses/03_phylogenetic_path/data/collaborator_max_svl.csv",
    )
    source_columns = (
        "source_taxon",
        "max_svl_mm",
        "max_svl_source_id",
        "max_svl_source_path",
    )
    if not _require_columns(source, source_columns, "collaborator SVL CSV", errors):
        return
    source_has_duplicates = source["source_taxon"].duplicated().any()
    if source_has_duplicates:
        errors.append("[SVL_SOURCE_DUPLICATE] collaborator SVL source labels are not unique")

    workbook = pd.read_excel(root / _SVL_SOURCE_PATH, dtype={"Species": str})
    if not _require_columns(workbook, ("Genus", "Species", "maxSVL"), "SVL workbook", errors):
        return

    projection = workbook.loc[
        workbook["Genus"].eq("Desmognathus"), ["Species", "maxSVL"]
    ].reset_index(drop=True)
    source_values = pd.to_numeric(source["max_svl_mm"], errors="coerce")
    projection_values = pd.to_numeric(projection["maxSVL"], errors="coerce")
    projection_matches = (
        projection["Species"].astype(str).tolist() == source["source_taxon"].tolist()
        and projection_values.tolist() == source_values.tolist()
    )
    if not projection_matches:
        errors.append(
            "[SVL_RAW_PROJECTION_MISMATCH] collaborator_max_svl.csv is not the ordered "
            "Genus == Desmognathus workbook projection"
        )
    source_values_invalid = source_values.isna().any() or (source_values <= 0).any()
    if source_values_invalid:
        errors.append("[SVL_VALUE_INVALID] collaborator maximum SVL values must be positive")
    if set(source["max_svl_source_id"]) != {_SVL_SOURCE_ID}:
        errors.append("[SVL_SOURCE_MISMATCH] collaborator CSV source_id changed")
    if set(source["max_svl_source_path"]) != {_SVL_SOURCE_PATH}:
        errors.append("[SVL_SOURCE_MISMATCH] collaborator CSV source path changed")
    raw_source_labels = set(source["source_taxon"])
    unused_aliases = sorted(
        label
        for (source_id, label) in alias_mapping
        if source_id == _SVL_SOURCE_ID and label not in raw_source_labels
    )
    if unused_aliases:
        errors.append(
            "[IDENTITY_ALIAS_UNUSED] collaborator-SVL aliases absent from the raw source: "
            + ";".join(unused_aliases)
        )
    if source_has_duplicates or source_values_invalid:
        return

    traits = _read_csv(
        root,
        "analyses/03_phylogenetic_path/data/organismal_traits.csv",
    )
    trait_columns = (
        "species",
        "max_svl_mm",
        "max_svl_source_taxon",
        "max_svl_source_id",
        "max_svl_source_path",
        "organismal_source_ids",
    )
    if not _require_columns(traits, trait_columns, "organismal traits", errors):
        return
    traits_invalid = traits["species"].duplicated().any() or len(traits) != 37
    if traits_invalid:
        errors.append("[SVL_PANEL_MISMATCH] organismal_traits must contain 37 unique species")
        return

    taxonomy_species = set(taxonomy["current_species"])
    trait_species = set(traits["species"])
    contributions: dict[str, list[tuple[str, float]]] = {}
    unused_labels: set[str] = set()
    for source_row, value in zip(source.to_dict(orient="records"), source_values, strict=True):
        label = source_row["source_taxon"]
        canonical = alias_mapping.get((_SVL_SOURCE_ID, label))
        if canonical is None and label in taxonomy_species:
            canonical = label
        if canonical is None or canonical not in trait_species:
            unused_labels.add(label)
            continue
        contributions.setdefault(canonical, []).append((label, float(value)))

    if unused_labels != _EXPECTED_UNUSED_SVL_LABELS:
        errors.append(
            "[SVL_UNUSED_SOURCE_MISMATCH] expected unused labels imitator;planiceps, observed "
            + ";".join(sorted(unused_labels))
        )

    indexed_traits = traits.set_index("species", drop=False)
    for species in sorted(trait_species):
        observed = contributions.get(species, [])
        if not observed:
            errors.append(f"[SVL_SOURCE_MISSING] no resolved source row for {species}")
            continue
        expected_value = max(value for _, value in observed)
        expected_labels = ";".join(sorted(label for label, _ in observed))
        trait_row = indexed_traits.loc[species]
        try:
            released_value = float(trait_row["max_svl_mm"])
        except (TypeError, ValueError):
            errors.append(f"[SVL_CANONICAL_VALUE_INVALID] nonnumeric SVL for {species}")
            continue
        if not math.isclose(released_value, expected_value, rel_tol=0, abs_tol=1e-12):
            errors.append(
                f"[SVL_CANONICAL_VALUE_MISMATCH] {species}: source max={expected_value:g}; "
                f"organismal_traits={released_value:g}"
            )
        if trait_row["max_svl_source_taxon"] != expected_labels:
            errors.append(
                f"[SVL_PROVENANCE_MISMATCH] {species}: expected source labels {expected_labels}"
            )
        if trait_row["max_svl_source_id"] != _SVL_SOURCE_ID:
            errors.append(f"[SVL_PROVENANCE_MISMATCH] {species}: wrong max_svl_source_id")
        if trait_row["max_svl_source_path"] != _SVL_SOURCE_PATH:
            errors.append(f"[SVL_PROVENANCE_MISMATCH] {species}: wrong max_svl_source_path")
        if _SVL_SOURCE_ID not in trait_row["organismal_source_ids"].split(";"):
            errors.append(f"[SVL_PROVENANCE_MISMATCH] {species}: source absent from source IDs")

    report_labels = {"fuscus", "imitator", "orestes-ac", "orestes-b", "planiceps"}
    missing_report_labels = sorted(report_labels - raw_source_labels)
    if missing_report_labels:
        errors.append(
            "[SVL_SOURCE_MISSING] required audited source labels missing: "
            + ";".join(missing_report_labels)
        )
        return
    source_index = source.assign(_value=source_values).set_index("source_taxon")
    orestes_ac = float(source_index.loc["orestes-ac", "_value"])
    orestes_b = float(source_index.loc["orestes-b", "_value"])
    orestes_max = max(orestes_ac, orestes_b)
    fuscus_svl = float(source_index.loc["fuscus", "_value"])
    rows.extend(
        (
            IdentityRow(
                "SVL",
                f"orestes-ac={orestes_ac:g}; orestes-b={orestes_b:g}",
                f"exact aliases; max -> orestes={orestes_max:g} mm",
            ),
            IdentityRow(
                "SVL",
                f"fuscus={fuscus_svl:g}",
                f"exact source label -> fuscus={fuscus_svl:g} mm",
            ),
            IdentityRow(
                "SVL",
                f"planiceps={float(source_index.loc['planiceps', '_value']):g}",
                "source-only row; not transferred to fuscus",
            ),
            IdentityRow(
                "SVL",
                f"imitator={float(source_index.loc['imitator', '_value']):g}",
                "source-only row outside the union37 organismal table",
            ),
        )
    )

def _tip_names(tree: Phylo.BaseTree.Tree) -> list[str]:
    return [str(tip.name) for tip in tree.get_terminals()]


def _clade_signature(tree: Phylo.BaseTree.Tree) -> frozenset[frozenset[str]]:
    all_tips = frozenset(_tip_names(tree))
    clades = []
    for clade in tree.get_nonterminals(order="preorder"):
        tips = frozenset(str(tip.name) for tip in clade.get_terminals())
        if 1 < len(tips) < len(all_tips):
            clades.append(tips)
    return frozenset(clades)


def _audit_tree(root: Path, rows: list[IdentityRow], errors: list[str]) -> None:
    path_traits = _read_csv(
        root,
        "analyses/03_phylogenetic_path/data/path24_traits.csv",
    )
    path_panel = _read_csv(root, "data/identity/path24_panel.csv")
    te_panel = _read_csv(root, "data/identity/te34_panel.csv")
    if not _require_columns(path_traits, ("species", "tree_tip"), "Path24 traits", errors):
        return
    if not _require_columns(path_panel, ("species",), "Path24 panel", errors):
        return
    if not _require_columns(te_panel, ("species",), "TE34 panel", errors):
        return

    selected = path_traits["tree_tip"].tolist()
    selected_invalid = len(selected) != 24 or len(set(selected)) != 24
    if selected_invalid:
        errors.append("[TREE_REPRESENTATIVE_CONFLICT] Path24 must select 24 unique exact tips")
    traits_invalid = len(path_traits) != 24 or path_traits["species"].duplicated().any()
    if traits_invalid:
        errors.append("[TREE_SPECIES_CONFLICT] Path24 traits must contain 24 unique species rows")
    for record in path_traits[["species", "tree_tip"]].to_dict(orient="records"):
        display_name = record["species"]
        if not display_name.startswith("D. "):
            errors.append(
                f"[TREE_SPECIES_LABEL_INVALID] expected exact 'D. ' prefix: {display_name}"
            )
            continue
        canonical_species = display_name.removeprefix("D. ")
        if canonical_species != record["tree_tip"]:
            errors.append(
                f"[TREE_SPECIES_TIP_MISMATCH] {display_name} must map to exact tip "
                f"{canonical_species}, observed {record['tree_tip']}"
            )
    panel_invalid = len(path_panel) != 24 or path_panel["species"].duplicated().any()
    if panel_invalid:
        errors.append("[TREE_PANEL_CONFLICT] Path24 panel must contain 24 unique species rows")
    if set(selected) != set(path_panel["species"]):
        errors.append("[TREE_PANEL_MISMATCH] path24_traits tree_tip values differ from Path24")
    te_panel_invalid = len(te_panel) != 34 or te_panel["species"].duplicated().any()
    if te_panel_invalid:
        errors.append("[TREE_TE_PANEL_CONFLICT] TE34 panel must contain 34 unique species rows")
    if selected_invalid or traits_invalid or panel_invalid or te_panel_invalid:
        return

    selected_set = set(selected)
    te_species = set(te_panel["species"])

    source = Phylo.read(
        root / "analyses/03_phylogenetic_path/trees/source_time_tree_46.tre",
        "newick",
    )
    released = Phylo.read(
        root / "analyses/03_phylogenetic_path/trees/path24_time_tree.nwk",
        "newick",
    )
    source_tips = _tip_names(source)
    released_tips = _tip_names(released)
    source_invalid = len(source_tips) != 46 or len(set(source_tips)) != 46
    if source_invalid:
        errors.append("[TREE_SOURCE_CONFLICT] source tree must have 46 unique exact tips")
        return
    source_tip_set = set(source_tips)
    missing = sorted(selected_set - source_tip_set)
    if missing:
        errors.append("[TREE_EXACT_TIP_MISSING] selected source tips missing: " + ";".join(missing))
    missing_te_tips = sorted(te_species - source_tip_set)
    if missing_te_tips:
        errors.append(
            "[TREE_TE_TIP_MISSING] TE34 exact source tips missing: " + ";".join(missing_te_tips)
        )
    released_invalid = set(released_tips) != selected_set or len(released_tips) != len(selected)
    if released_invalid:
        errors.append("[TREE_PANEL_MISMATCH] released Path24 tree has the wrong exact tip set")
    allowed_alternates = {
        alternate for alternates in _TREE_ALTERNATES.values() for alternate in alternates
    }
    unregistered_variants = sorted(
        tip
        for tip in source_tips
        if tip not in allowed_alternates
        and any(
            tip.startswith((f"{species}_", f"{species}-"))
            for species in selected_set | te_species
        )
    )
    if unregistered_variants:
        errors.append(
            "[TREE_UNREGISTERED_VARIANT] suffix-like source tips require an explicit decision: "
            + ";".join(unregistered_variants)
        )
    if missing or missing_te_tips or released_invalid:
        return

    pruned = copy.deepcopy(source)
    for tip in list(pruned.get_terminals()):
        if str(tip.name) not in selected_set:
            pruned.prune(tip)

    if _clade_signature(pruned) != _clade_signature(released):
        errors.append("[TREE_PRUNING_MISMATCH] released Path24 topology is not the exact source prune")

    max_pairwise_difference = 0.0
    for left, right in combinations(sorted(selected_set), 2):
        difference = abs(pruned.distance(left, right) - released.distance(left, right))
        max_pairwise_difference = max(max_pairwise_difference, difference)
        if difference > _TREE_DISTANCE_TOLERANCE:
            errors.append(
                f"[TREE_DISTANCE_MISMATCH] {left};{right} differs by {difference:.8g}"
            )
            break
    for tip in sorted(selected_set):
        difference = abs(pruned.distance(tip) - released.distance(tip))
        if difference > _TREE_DISTANCE_TOLERANCE:
            errors.append(f"[TREE_ROOT_DISTANCE_MISMATCH] {tip} differs by {difference:.8g}")
            break

    for representative, alternates in _TREE_ALTERNATES.items():
        absent = [tip for tip in (representative, *alternates) if tip not in source_tips]
        retained_alternates = [tip for tip in alternates if tip in released_tips]
        if absent:
            errors.append(
                f"[TREE_SOURCE_TIP_MISSING] {representative} case missing: {';'.join(absent)}"
            )
        if representative not in released_tips or retained_alternates:
            errors.append(
                f"[TREE_REPRESENTATIVE_MISMATCH] {representative} exact-tip selection changed"
            )
        rows.append(
            IdentityRow(
                "TREE",
                f"{representative}; alternates={';'.join(alternates)}",
                "retain exact unsuffixed representative; prune alternates; no MRCA merge",
            )
        )

    if "planiceps" not in source_tips or "planiceps" in released_tips:
        errors.append("[TREE_PLANICEPS_SCOPE] planiceps must remain a distinct pruned source tip")
    rows.append(
        IdentityRow(
            "TREE",
            "planiceps",
            "distinct source tip pruned; not substituted for fuscus",
        )
    )

    rows.append(
        IdentityRow(
            "TREE",
            "Path24 exact prune",
            f"24 tips; topology preserved; max patristic rounding={max_pairwise_difference:.8g}",
        )
    )


def _audit_fuscus(
    sources: pd.DataFrame,
    taxonomy: pd.DataFrame,
    aliases: pd.DataFrame,
    root: Path,
    rows: list[IdentityRow],
    errors: list[str],
) -> None:
    active_id = "genomic_fuscus_srx20497025_gca032353935_1"
    decision_id = "pyron_personal_communication_2026_07_09"
    active = sources.loc[sources["source_id"].eq(active_id)]
    expected = {
        "source_type": "genomic_resource",
        "evidence_stream": "genomic",
        "sra_accessions": "SRX20497025",
        "assembly_accession": "GCA_032353935.1",
        "biosample_accession": "SAMN34288323",
        "public_taxon_name": "Desmognathus planiceps",
        "analysis_taxon_name": "Desmognathus fuscus",
        "voucher_or_isolate": "IRGN:RAP2245",
        "tree_tip": "fuscus",
        "analysis_role": "current_comparative_te_ltr_ectopic",
        "inclusion_decision": "included_current",
        "decision_basis_source_id": decision_id,
        "decision_date": "2026-07-09",
    }
    if len(active) != 1:
        errors.append(f"[FUSCUS_RESOURCE_MISSING] expected one exact source row: {active_id}")
    else:
        for column, value in expected.items():
            if active.iloc[0][column] != value:
                errors.append(
                    f"[FUSCUS_RESOURCE_MISMATCH] {active_id}.{column} must equal {value}"
                )

    decision = sources.loc[sources["source_id"].eq(decision_id)]
    decision_fields = {
        "source_type": "personal_communication",
        "evidence_stream": "genomic_resource_identity",
        "analysis_taxon_name": "Desmognathus fuscus",
        "tree_tip": "fuscus",
        "analysis_role": "decision_basis",
        "inclusion_decision": "not_applicable",
        "decision_date": "2026-07-09",
    }
    if (
        len(decision) != 1
        or any(decision.iloc[0][column] != value for column, value in decision_fields.items())
        or not all(
            token in decision.iloc[0]["notes"] for token in ("SRX20497025", "GCA_032353935.1")
        )
    ):
        errors.append("[FUSCUS_DECISION_MISSING] expert accession decision is incomplete")

    active_reidentifications = sources.loc[
        sources["public_taxon_name"].eq("Desmognathus planiceps")
        & sources["analysis_taxon_name"].eq("Desmognathus fuscus")
        & sources["inclusion_decision"].eq("included_current")
    ]
    if set(active_reidentifications["source_id"]) != {active_id} or len(
        active_reidentifications
    ) != 1:
        errors.append(
            "[FUSCUS_SCOPE_BREACH] the approved accession pair must be the only active "
            "structured planiceps-to-fuscus mapping"
        )

    for relative_path in (
        "data/identity/te34_panel.csv",
        "data/identity/path24_panel.csv",
        "data/identity/analysis_availability.csv",
    ):
        panel = _read_csv(root, relative_path)
        if not _require_columns(
            panel,
            ("species", "te_sra_accession", "te_assembly_accession"),
            relative_path,
            errors,
        ):
            continue
        fuscus = panel.loc[panel["species"].eq("fuscus")]
        if len(fuscus) != 1 or (
            fuscus.iloc[0]["te_sra_accession"], fuscus.iloc[0]["te_assembly_accession"]
        ) != ("SRX20497025", "GCA_032353935.1"):
            errors.append(f"[FUSCUS_PANEL_MISMATCH] accession pair changed in {relative_path}")

    taxonomy_fuscus = taxonomy.loc[taxonomy["current_species"].eq("fuscus")]
    if len(taxonomy_fuscus) != 1 or taxonomy_fuscus.iloc[0]["search_name_primary"] != (
        "Desmognathus fuscus"
    ):
        errors.append("[FUSCUS_TAXONOMY_MISSING] canonical fuscus taxonomy row is missing")
    executable_mapping_fields = [
        column for column in ("search_name_legacy", "mapping_status") if column in taxonomy.columns
    ]
    if not taxonomy_fuscus.empty and any(
        taxonomy_fuscus.iloc[0][column] for column in executable_mapping_fields
    ):
        errors.append(
            "[FUSCUS_SCOPE_BREACH] accession-specific planiceps mapping belongs only in the source manifest"
        )
    if not aliases.empty and (
        aliases["source_label"].eq("planiceps")
        & aliases["canonical_species"].eq("fuscus")
    ).any():
        errors.append("[FUSCUS_SCOPE_BREACH] planiceps-to-fuscus cannot be a general source alias")

    rows.append(
        IdentityRow(
            "GENOME",
            "SRX20497025 / GCA_032353935.1",
            "public planiceps label -> analysis fuscus for this accession pair only",
        )
    )


def audit_identity(root: Path) -> IdentityAudit:
    """Audit every active identity resolution without changing repository files."""

    root = root.resolve()
    rows: list[IdentityRow] = []
    errors: list[str] = []
    try:
        sources = _read_csv(root, "data/identity/source_manifest.csv")
        taxonomy = _read_csv(root, "data/identity/species_taxonomy_crosswalk.csv")
        source_columns = (
            "source_id",
            "source_type",
            "notes",
            "evidence_stream",
            "sra_accessions",
            "assembly_accession",
            "biosample_accession",
            "public_taxon_name",
            "analysis_taxon_name",
            "voucher_or_isolate",
            "tree_tip",
            "analysis_role",
            "inclusion_decision",
            "decision_basis_source_id",
            "decision_date",
        )
        if not _require_columns(sources, source_columns, "source manifest", errors):
            return IdentityAudit(tuple(rows), tuple(errors))
        if not _require_columns(
            taxonomy,
            ("current_species", "search_name_primary"),
            "taxonomy crosswalk",
            errors,
        ):
            return IdentityAudit(tuple(rows), tuple(errors))
        if sources["source_id"].duplicated().any():
            errors.append("[IDENTITY_SOURCE_CONFLICT] source_manifest source_id values are not unique")
        if taxonomy["current_species"].duplicated().any():
            errors.append("[IDENTITY_CANONICAL_CONFLICT] taxonomy species values are not unique")
        decision = sources.loc[sources["source_id"].eq(_IDENTITY_DECISION_ID)]
        if len(decision) != 1 or any(
            decision.iloc[0][column] != value
            for column, value in {
                "source_type": "author_curator_decision",
                "evidence_stream": "identity_resolution",
                "analysis_role": "decision_basis",
                "decision_date": "2026-08-10",
            }.items()
        ):
            errors.append("[IDENTITY_DECISION_MISSING] author identity-resolution decision is missing")

        aliases, alias_mapping = _load_aliases(root, sources, taxonomy, errors)
        _audit_svl(root, taxonomy, alias_mapping, rows, errors)
        _audit_tree(root, rows, errors)
        _audit_fuscus(sources, taxonomy, aliases, root, rows, errors)
    except (FileNotFoundError, OSError, ValueError, KeyError) as error:
        errors.append(f"[IDENTITY_IO] {type(error).__name__}: {error}")
    return IdentityAudit(tuple(rows), tuple(errors))
