#!/usr/bin/env python3
"""Reconstruct a configured-denominator dnaPipeTE repeat-load sensitivity."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[2]
MASS_INPUT = (
    PROJECT_ROOT
    / "results"
    / "data"
    / "corrected"
    / "dnapipete"
    / "dnapipete_mass_accounting_analysis18_v2.csv"
)
MASS_MANIFEST = MASS_INPUT.with_suffix(".manifest.json")
CLASS_BREAKDOWN = PROJECT_ROOT / "results" / "data" / "dnaPipeTE_class_breakdown.csv"
OUTPUT = (
    PROJECT_ROOT
    / "results"
    / "data"
    / "corrected"
    / "dnapipete"
    / "dnapipete_absolute_load_sensitivity_analysis18_v1.csv"
)
MANIFEST = OUTPUT.with_suffix(".manifest.json")
REPORT = (
    PROJECT_ROOT
    / "plans"
    / "publication-readiness-deep-audit"
    / "dnapipete_absolute_load_sensitivity_analysis18_v1.md"
)

CONFIGURED_GENOME_SIZE_BP = 15_000_000_000
CONFIGURED_GENOME_COVERAGE = 0.1
CONFIGURED_QUANTIFICATION_SAMPLE_BP = int(
    CONFIGURED_GENOME_SIZE_BP * CONFIGURED_GENOME_COVERAGE
)
HISTORICAL_RUN_SCRIPT_COMMIT = "cc64ffb134968250671ff0e3456255063a1aba97"
HISTORICAL_RUN_SCRIPT_BLOB = "738c7cece113ce03b34f7523bf2139431646fac1"
DNAPIPETE_METHOD_URL = "https://pmc.ncbi.nlm.nih.gov/articles/PMC4419797/"
DNAPIPETE_TUTORIAL_URL = "https://tehub.org/tutorials/docs/dnaPipeTE"


def _canonical_species(value: object) -> str:
    text = str(value).strip()
    if text.startswith("D. "):
        text = text[3:]
    elif text.startswith("D."):
        text = text[2:]
    return text.strip()


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _portable(path: Path) -> str:
    return str(path.resolve().relative_to(PROJECT_ROOT.resolve()))


def build_absolute_load(
    mass: pd.DataFrame,
    class_breakdown: pd.DataFrame,
    configured_genome_size_bp: int,
    configured_genome_coverage: float,
) -> pd.DataFrame:
    """Partition repeat-aligned mass against the configured sample denominator."""

    quantification_sample_bp = int(
        configured_genome_size_bp * configured_genome_coverage
    )
    if quantification_sample_bp <= 0:
        raise ValueError("Quantification-sample denominator must be positive")
    required = {
        "species",
        "te_sra_accession",
        "te_assembly_accession",
        "total_aligned_bases",
        "class_retained_aligned_bases",
        "class_unresolved_aligned_bases",
    }
    missing = sorted(required.difference(mass.columns))
    if missing:
        raise ValueError(f"Mass ledger is missing columns: {missing}")

    ledger = mass.copy()
    numeric_mass = [
        "total_aligned_bases",
        "class_retained_aligned_bases",
        "class_unresolved_aligned_bases",
    ]
    ledger[numeric_mass] = ledger[numeric_mass].apply(pd.to_numeric, errors="coerce")
    if ledger[numeric_mass].isna().any().any() or ledger[numeric_mass].lt(0).any().any():
        raise ValueError("Mass ledger contains invalid aligned-base counts")
    if not (
        ledger["class_retained_aligned_bases"]
        + ledger["class_unresolved_aligned_bases"]
    ).equals(ledger["total_aligned_bases"]):
        raise ValueError("Class retained plus unresolved bases do not conserve total mass")

    breakdown = class_breakdown.copy()
    breakdown.index = pd.Index(
        [_canonical_species(value) for value in breakdown.index], name="species"
    )
    if breakdown.index.duplicated().any():
        raise ValueError("Class breakdown has duplicate canonical species")
    breakdown = breakdown.apply(pd.to_numeric, errors="coerce")
    selected = breakdown.reindex(ledger["species"])
    if selected.isna().any().any():
        raise ValueError("Class breakdown is missing analysis species or numeric values")
    row_sums = selected.sum(axis=1)
    if not row_sums.sub(100.0).abs().le(1e-8).all():
        raise ValueError("Class breakdown rows must sum to 100")

    te_columns = [
        column
        for column in selected.columns
        if column.startswith("DNAtransposons")
        or column.startswith("Retrotransposons")
    ]
    if not te_columns or "Other" not in selected or "Unknown" not in selected:
        raise ValueError("Class breakdown lacks the expected TE, Other, or Unknown bins")

    selected = selected.set_axis(ledger.index)
    retained = ledger["class_retained_aligned_bases"].astype(float)
    te_bases = retained * selected[te_columns].sum(axis=1) / 100.0
    other_bases = retained * selected["Other"] / 100.0
    unknown_bases = retained * selected["Unknown"] / 100.0
    unresolved_bases = ledger["class_unresolved_aligned_bases"].astype(float)

    result = ledger[
        ["species", "te_sra_accession", "te_assembly_accession"]
    ].copy()
    result["configured_genome_size_bp"] = configured_genome_size_bp
    result["configured_genome_coverage"] = configured_genome_coverage
    result["configured_quantification_sample_bp"] = quantification_sample_bp
    result["repeat_aligned_bases"] = ledger["total_aligned_bases"].astype(int)
    result["repeat_aligned_fraction"] = (
        result["repeat_aligned_bases"] / quantification_sample_bp
    )
    result["te_classified_bases"] = te_bases
    result["te_classified_fraction"] = te_bases / quantification_sample_bp
    result["other_repeat_bases"] = other_bases
    result["other_repeat_fraction"] = other_bases / quantification_sample_bp
    result["class_unknown_bases"] = unknown_bases
    result["class_unknown_fraction"] = unknown_bases / quantification_sample_bp
    result["class_unresolved_bases"] = unresolved_bases
    result["class_unresolved_fraction"] = unresolved_bases / quantification_sample_bp
    result["class_component_sum_fraction"] = (
        te_bases + other_bases + unknown_bases + unresolved_bases
    ) / quantification_sample_bp
    difference = (
        result["class_component_sum_fraction"] - result["repeat_aligned_fraction"]
    ).abs()
    if not difference.le(1e-10).all():
        raise RuntimeError("Class components do not conserve repeat-aligned fraction")
    if result["repeat_aligned_fraction"].gt(1).any():
        raise RuntimeError("Repeat-aligned mass exceeds the configured sample denominator")
    return result.sort_values("species").reset_index(drop=True)


def _render_report(result: pd.DataFrame) -> str:
    fuscus = result.loc[result["species"].eq("fuscus")].iloc[0]
    return f"""# dnaPipeTE configured-denominator absolute-load sensitivity

Git history recovers the deleted upstream run script at commit `{HISTORICAL_RUN_SCRIPT_COMMIT}` (blob `{HISTORICAL_RUN_SCRIPT_BLOB}`). It configured a 15,000,000,000-bp genome size and 0.1× coverage, implying an independent 1,500,000,000-bp quantification sample for every run. dnaPipeTE documents `aligned_bases` as bases from that independent sample mapped back to annotated repeat contigs.

Dividing conserved repeat-aligned bases by the configured quantification-sample size gives:

- 18 species, exactly the declared final TE/genome panel
- repeat-aligned fraction range: {result.repeat_aligned_fraction.min():.3%}–{result.repeat_aligned_fraction.max():.3%}
- median repeat-aligned fraction: {result.repeat_aligned_fraction.median():.3%}
- *D. fuscus* (`{fuscus.te_sra_accession}` / `{fuscus.te_assembly_accession}`): {fuscus.repeat_aligned_fraction:.3%} repeat-aligned and {fuscus.te_classified_fraction:.3%} assigned to a DNA- or retrotransposon class

This is a **sensitivity estimate**, not yet a confirmatory predictor. The command was recovered from Git history rather than runtime logs; the container version/digest and custom-library checksum are absent; all species used the same configured 15-Gb genome size; only the nuclear-filtered first read mate was supplied; and replicate quantification samples or sampling uncertainty were not retained. The image-IOD genome-size measurements are from independent microscopy specimens and are not substituted into this denominator.

Primary method sources: [Goubert et al. dnaPipeTE paper]({DNAPIPETE_METHOD_URL}) and [author tutorial]({DNAPIPETE_TUTORIAL_URL}).
"""


def main() -> None:
    occupied = [path for path in (OUTPUT, MANIFEST, REPORT) if path.exists()]
    if occupied:
        raise FileExistsError(
            "Non-destructive absolute-load sensitivity requires unused outputs: "
            + ", ".join(str(path) for path in occupied)
        )
    mass = pd.read_csv(MASS_INPUT)
    class_breakdown = pd.read_csv(CLASS_BREAKDOWN, index_col=0)
    result = build_absolute_load(
        mass,
        class_breakdown,
        CONFIGURED_GENOME_SIZE_BP,
        CONFIGURED_GENOME_COVERAGE,
    )

    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    REPORT.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(OUTPUT, index=False)
    REPORT.write_text(_render_report(result), encoding="utf-8")
    manifest = {
        "analysis_scope": "final_te_genome_primary_mediumplus_panel",
        "analysis_status": "sensitivity_only",
        "eligible_for_confirmatory_path_analysis": False,
        "n_species": len(result),
        "mass_input": _portable(MASS_INPUT),
        "mass_input_sha256": _sha256(MASS_INPUT),
        "mass_manifest": _portable(MASS_MANIFEST),
        "mass_manifest_sha256": _sha256(MASS_MANIFEST),
        "class_breakdown_input": _portable(CLASS_BREAKDOWN),
        "class_breakdown_input_sha256": _sha256(CLASS_BREAKDOWN),
        "output": _portable(OUTPUT),
        "output_sha256": _sha256(OUTPUT),
        "report": _portable(REPORT),
        "report_sha256": _sha256(REPORT),
        "historical_run_script_git_commit": HISTORICAL_RUN_SCRIPT_COMMIT,
        "historical_run_script_path_at_commit": "dnaPipeTE.sh",
        "historical_run_script_git_blob": HISTORICAL_RUN_SCRIPT_BLOB,
        "configured_genome_size_bp": CONFIGURED_GENOME_SIZE_BP,
        "configured_genome_coverage": CONFIGURED_GENOME_COVERAGE,
        "configured_quantification_sample_bp": CONFIGURED_QUANTIFICATION_SAMPLE_BP,
        "configured_sample_formula": "genome_size_bp * genome_coverage",
        "input_mate": "nuclear_filtered_R1_only",
        "sample_number": 2,
        "repeatmasker_annotation_threshold": 0.15,
        "repeat_library_historical_path": "/nfs/home/jlamb/TE_libs/dedupe_telib.fasta",
        "repeat_library_sha256_available": False,
        "container_version_or_digest_available": False,
        "runtime_log_available": False,
        "method_sources": [DNAPIPETE_METHOD_URL, DNAPIPETE_TUTORIAL_URL],
        "out_of_panel_species_processed": False,
    }
    MANIFEST.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(f"Wrote {OUTPUT}")
    print(f"Wrote {MANIFEST}")
    print(f"Wrote {REPORT}")


if __name__ == "__main__":
    main()
