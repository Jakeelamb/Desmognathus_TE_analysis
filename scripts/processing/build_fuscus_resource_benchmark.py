#!/usr/bin/env python3
"""Build an accession-labeled benchmark of the two fuscus genome resources.

The current repository and the 2025 publication use different sequencing and
repeat-annotation methods. This script therefore preserves denominators and
marks the comparison diagnostic rather than treating values as interchangeable.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict, List

import pandas as pd


CURRENT_SRA = "SRX20497025"
CURRENT_ASSEMBLY = "GCA_032353935.1"
CURRENT_ASSEMBLY_BP = 652_813_091
VALIDATION_ASSEMBLY = "GCA_050004315.1"
VALIDATION_ASSEMBLY_BP = 16_117_580_045
PAPER_DOI = "https://doi.org/10.1093/g3journal/jkaf157"
ZENODO_RECORD = "https://zenodo.org/records/15255946"


def _load_species_row(path: Path, species: str) -> pd.Series:
    table = pd.read_csv(path, index_col=0)
    if species not in table.index:
        raise ValueError(f"{species} is absent from {path}")
    return table.loc[species]


def _row(
    metric: str,
    current_value: object,
    current_denominator: str,
    published_value: object,
    published_denominator: str,
    interpretation: str,
) -> Dict[str, object]:
    return {
        "metric": metric,
        "current_sra_accession": CURRENT_SRA,
        "current_assembly_accession": CURRENT_ASSEMBLY,
        "current_value": current_value,
        "current_denominator": current_denominator,
        "validation_assembly_accession": VALIDATION_ASSEMBLY,
        "published_value": published_value,
        "published_denominator": published_denominator,
        "directly_comparable": False,
        "interpretation": interpretation,
        "published_source": PAPER_DOI,
    }


def build_benchmark(project_root: Path) -> pd.DataFrame:
    """Return the current-versus-published fuscus diagnostic table."""

    order = _load_species_row(
        project_root / "results" / "data" / "dnaPipeTE_order_breakdown.csv",
        "D.fuscus",
    )
    superfamily = _load_species_row(
        project_root / "results" / "data" / "dnaPipeTE_superfamily_breakdown.csv",
        "D.fuscus",
    )
    absolute_load = pd.read_csv(
        project_root
        / "results"
        / "data"
        / "corrected"
        / "dnapipete"
        / "dnapipete_absolute_load_sensitivity_analysis18_v1.csv"
    )
    absolute_load = absolute_load.loc[absolute_load["species"].eq("fuscus")]
    if len(absolute_load) != 1:
        raise ValueError("Expected exactly one fuscus absolute-load sensitivity row")
    absolute_load = absolute_load.iloc[0]
    if (
        absolute_load["te_sra_accession"] != CURRENT_SRA
        or absolute_load["te_assembly_accession"] != CURRENT_ASSEMBLY
    ):
        raise ValueError("Fuscus absolute-load sensitivity uses an unexpected resource")
    ltr_age = pd.read_csv(
        project_root / "results" / "data" / "ltr_age" / "ltr_age_species_summary.csv"
    )
    ltr_age = ltr_age.loc[ltr_age["species"].eq("D.fuscus")]
    if len(ltr_age) != 1:
        raise ValueError("Expected exactly one D.fuscus LTR-age summary row")
    ltr_age = ltr_age.iloc[0]

    rows: List[Dict[str, object]] = [
        _row(
            "assembly_span_gb",
            CURRENT_ASSEMBLY_BP / 1_000_000_000,
            "assembled_sequence_gb",
            VALIDATION_ASSEMBLY_BP / 1_000_000_000,
            "assembled_sequence_gb",
            "The short-read assembly spans only "
            f"{100 * CURRENT_ASSEMBLY_BP / VALIDATION_ASSEMBLY_BP:.2f}% of the 2025 assembly; "
            "this is an assembly-completeness contrast, not a genome-size estimate comparison.",
        ),
        _row(
            "repeat_aligned_sensitivity",
            100 * float(absolute_load["repeat_aligned_fraction"]),
            "percent_of_configured_1.5Gb_dnaPipeTE_quantification_sample",
            75.0,
            "whole_assembly_percent",
            "The recovered dnaPipeTE denominator permits a sensitivity estimate, but runtime logs, container/library hashes, specimens, and methods differ.",
        ),
        _row(
            "TE_classified_sensitivity",
            100 * float(absolute_load["te_classified_fraction"]),
            "DNA_or_retrotransposon_class_percent_of_configured_1.5Gb_sample",
            75.0,
            "whole_assembly_repeat_percent",
            "This excludes current Other, Unknown, and class-unresolved repeat mass and is not directly equivalent to EarlGrey coverage.",
        ),
        _row(
            "LTR",
            float(order["LTR"]),
            "classified_dnaPipeTE_order_percent",
            36.0,
            "whole_assembly_percent",
            "Both values describe LTRs but use different methods and denominators.",
        ),
        _row(
            "LINE",
            float(order["LINE"]),
            "classified_dnaPipeTE_order_percent",
            15.0,
            "whole_assembly_percent",
            "Numerical similarity does not establish agreement because denominators differ.",
        ),
        _row(
            "Gypsy_Ty3",
            float(superfamily["Gypsy"]),
            "represented_dnaPipeTE_superfamily_percent",
            31.0,
            "published_repeat_coverage_percent",
            "The paper's Ty3 value is repeat coverage; the current value is closed superfamily composition.",
        ),
        _row(
            "LTR_TE_normalized_approx",
            float(order["LTR"]),
            "classified_dnaPipeTE_order_percent",
            48.0,
            "approximate_percent_of_published_TE_coverage_36_over_75",
            "A closer denominator still remains method-dependent and excludes the paper's unclassified fraction.",
        ),
        _row(
            "LINE_TE_normalized_approx",
            float(order["LINE"]),
            "classified_dnaPipeTE_order_percent",
            20.0,
            "approximate_percent_of_published_TE_coverage_15_over_75",
            "A closer denominator still remains method-dependent and excludes the paper's unclassified fraction.",
        ),
        _row(
            "paired_LTR_candidates",
            int(ltr_age["n_pairs_estimated"]),
            "pairs_with_distance_estimate",
            pd.NA,
            "not_reported_as_equivalent_metric",
            "Baseline for rerunning the same paired-LTR procedure on the 2025 assembly.",
        ),
        _row(
            "paired_LTR_high_confidence",
            int(ltr_age["n_high_confidence_pairs_estimated"]),
            "high_confidence_pairs_with_distance_estimate",
            pd.NA,
            "not_reported_as_equivalent_metric",
            "The EarlGrey final GFF alone may not retain the paired 5-prime and 3-prime boundaries needed here.",
        ),
        _row(
            "paired_LTR_median_K2P",
            float(ltr_age["median_k2p_distance"]),
            "paired_LTR_K2P_distance",
            pd.NA,
            "paper_reports_repeat_to_consensus_KIMURA80_instead",
            "Paired-LTR K2P and repeat-to-consensus divergence are not equivalent clocks.",
        ),
    ]
    return pd.DataFrame(rows)


def render_markdown(benchmark: pd.DataFrame) -> str:
    """Render a compact, shareable interpretation of the benchmark."""

    metrics = benchmark.set_index("metric")
    span = metrics.loc["assembly_span_gb"]
    repeat_load = metrics.loc["repeat_aligned_sensitivity"]
    te_load = metrics.loc["TE_classified_sensitivity"]
    ltr = metrics.loc["LTR"]
    line = metrics.loc["LINE"]
    gypsy = metrics.loc["Gypsy_Ty3"]
    pairs = metrics.loc["paired_LTR_candidates"]
    high_confidence = metrics.loc["paired_LTR_high_confidence"]
    k2p = metrics.loc["paired_LTR_median_K2P"]

    return f"""# Fuscus genomic-resource benchmark

**Generated comparison:** `{CURRENT_SRA} / {CURRENT_ASSEMBLY}` versus `{VALIDATION_ASSEMBLY}`
**Interpretation:** diagnostic method-and-resource sensitivity only; no value below authorizes silent substitution.

## What is established now

| Metric | Current repository resource | 2025 published resource | Comparability |
|---|---:|---:|---|
| Assembly span | {span.current_value:.3f} Gb | {span.published_value:.3f} Gb | Same unit, radically different completeness |
| Repeat-aligned load | {repeat_load.current_value:.3f}% of configured 1.5-Gb quantification sample | about {repeat_load.published_value:.0f}% of the assembly | Same percentage scale; different specimens and methods |
| TE-classified load | {te_load.current_value:.3f}% of configured 1.5-Gb quantification sample | about {te_load.published_value:.0f}% whole-assembly repeat coverage | Current value excludes Other/Unknown/unresolved mass |
| LTR | {ltr.current_value:.3f}% of classified dnaPipeTE orders | about {ltr.published_value:.0f}% of the assembly | Different method and denominator |
| LINE | {line.current_value:.3f}% of classified dnaPipeTE orders | about {line.published_value:.0f}% of the assembly | Different method and denominator |
| Gypsy / Ty3 | {gypsy.current_value:.3f}% of represented superfamilies | about {gypsy.published_value:.0f}% of repeat coverage | Different method and denominator |

The old assembly spans only {100 * CURRENT_ASSEMBLY_BP / VALIDATION_ASSEMBLY_BP:.2f}% of the 2025 assembly. That makes the modern resource particularly important for assembly-dependent LTR and structural analyses. It does not invalidate the current shallow-read dnaPipeTE comparison, which should remain a separately named branch. The recovered dnaPipeTE denominator makes the {repeat_load.current_value:.3f}% repeat-aligned value useful as a sensitivity check, not a silent replacement for the published assembly annotation.

The current paired-LTR baseline is {int(pairs.current_value)} estimated pairs, {int(high_confidence.current_value)} high-confidence pairs, and median K2P {k2p.current_value:.6f}.

## Public-data blocker

The [2025 paper]({PAPER_DOI}) states that its EarlGrey annotation is deposited, but the live [Zenodo record]({ZENODO_RECORD}) checked on 2026-07-09 exposes only a 4.76-GB compressed FASTA and a 1.30-MB AGP. The promised filtered-repeat GFF is not present. Downloading the FASTA alone would not reproduce the published repeat classification and would consume about 21 GB when compressed and decompressed copies are retained.

## Minimum asset request

Request the EarlGrey `*_summaryFiles/` directory plus:

- `aDesFus1-2.1.filteredRepeats.gff`;
- `*.highLevelCount.txt`;
- `*_divergence_summary_table.tsv`;
- `*_combined_library.fasta`;
- the raw LTR_FINDER GFF;
- the EarlGrey 4.4.2 command/log and checksums.

The filtered GFF is enough for TE bp, count, length, class/family, and repeat-to-consensus divergence comparisons without the full FASTA. The raw LTR_FINDER GFF is needed to determine whether the repository's paired 5-prime/3-prime LTR K2P metric can be reproduced.

## Decision

Keep `{CURRENT_SRA} / {CURRENT_ASSEMBLY}` as the current comparative branch. Add `{VALIDATION_ASSEMBLY}` as an accession-labeled validation branch once the annotation assets are available. Do not merge either genomic specimen identity into the independently collected microscopy records.
"""


def main() -> None:
    project_root = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=project_root / "plans" / "publication-readiness-deep-audit",
    )
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    benchmark = build_benchmark(project_root)
    csv_path = args.output_dir / "fuscus_resource_benchmark.csv"
    markdown_path = args.output_dir / "fuscus_resource_benchmark.md"
    benchmark.to_csv(csv_path, index=False)
    markdown_path.write_text(render_markdown(benchmark), encoding="utf-8")
    print(f"Wrote {csv_path}")
    print(f"Wrote {markdown_path}")


if __name__ == "__main__":
    main()
