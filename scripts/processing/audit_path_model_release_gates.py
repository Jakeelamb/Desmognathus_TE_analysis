#!/usr/bin/env python3
"""Audit stored phylopath rankings with fail-closed publication gates."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Dict, List

import numpy as np
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[2]
RESULTS_DIR = PROJECT_ROOT / "path_analysis" / "results"
DEFAULT_OUTPUT_DIR = PROJECT_ROOT / "plans" / "publication-readiness-deep-audit"
AUDIT_NAME = "path_model_release_gate_audit"
PHYLOPATH_PAPER = "https://pmc.ncbi.nlm.nih.gov/articles/PMC5923215/"


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def audit_ranking(
    ranking: pd.DataFrame,
    family: str,
    n_species: int,
    alpha: float = 0.05,
    competitive_delta: float = 2.0,
) -> Dict[str, object]:
    """Classify one stored model ranking without choosing a replacement model."""

    required = {"model", "p", "CICc", "delta_CICc", "l", "w"}
    missing = sorted(required.difference(ranking.columns))
    if missing:
        raise ValueError(f"{family} ranking is missing columns: {missing}")
    if ranking.empty:
        raise ValueError(f"{family} ranking is empty")

    table = ranking.copy()
    for column in ["p", "CICc", "delta_CICc", "l", "w"]:
        table[column] = pd.to_numeric(table[column], errors="coerce")
    finite_criterion = np.isfinite(table["CICc"])
    finite_ranking = table[["CICc", "delta_CICc", "l", "w"]].apply(
        lambda column: np.isfinite(column)
    ).all(axis=1)
    global_fit_pass = np.isfinite(table["p"]) & table["p"].ge(alpha)
    eligible = finite_criterion & global_fit_pass

    top = table.iloc[0]
    top_fit_pass = bool(global_fit_pass.iloc[0])
    ranking_complete = bool(finite_ranking.all())
    weights_sum = float(table["w"].sum()) if table["w"].notna().all() else np.nan
    weights_valid = bool(
        ranking_complete and np.isfinite(weights_sum) and abs(weights_sum - 1.0) <= 1e-6
    )
    best_supported_model = pd.NA
    if eligible.any():
        eligible_table = table.loc[eligible]
        best_supported_model = eligible_table.loc[
            eligible_table["CICc"].idxmin(), "model"
        ]

    competitive_supported = (
        global_fit_pass
        & finite_ranking
        & table["delta_CICc"].le(competitive_delta)
    )
    n_competitive = int(competitive_supported.sum())

    if not ranking_complete or not weights_valid:
        release_status = "blocked_nonfinite_ranking"
    elif not eligible.any():
        release_status = "blocked_no_globally_supported_model"
    elif not top_fit_pass:
        release_status = "blocked_top_model_rejected"
    elif n_competitive > 1:
        release_status = "supported_competitive_model_set"
    else:
        release_status = "supported_unique_top_model"

    winner_claim_allowed_by_ranking_gate = (
        release_status == "supported_unique_top_model"
    )
    model_averaging_allowed_by_ranking_gate = (
        release_status == "supported_competitive_model_set"
    )
    return {
        "family": family,
        "n_species": int(n_species),
        "n_models": len(table),
        "n_global_fit_pass": int(global_fit_pass.sum()),
        "n_finite_cicc": int(finite_criterion.sum()),
        "n_competitive_supported_delta_le_2": n_competitive,
        "stored_top_model": str(top["model"]),
        "stored_top_fisher_c_p": float(top["p"]),
        "stored_top_cicc": float(top["CICc"]),
        "stored_top_weight": float(top["w"]),
        "stored_top_global_fit_pass": top_fit_pass,
        "ranking_all_finite": ranking_complete,
        "weights_sum": weights_sum,
        "weights_valid": weights_valid,
        "best_supported_model": best_supported_model,
        "release_status": release_status,
        "winner_claim_allowed_by_ranking_gate": winner_claim_allowed_by_ranking_gate,
        "model_averaging_allowed_by_ranking_gate": model_averaging_allowed_by_ranking_gate,
        "publication_claim_allowed": False,
        "publication_blocker": "historical_pre_audit_inputs_require_corrected_rerun",
        "alpha": alpha,
        "competitive_delta_cicc": competitive_delta,
    }


def _render_markdown(audit: pd.DataFrame) -> str:
    status_counts = audit["release_status"].value_counts()
    lines = [
        "# Phylogenetic path-model release-gate audit",
        "",
        "This audit is non-destructive: it evaluates stored ranking tables and does not rewrite historical model, edge, or figure outputs.",
        "",
        "The gate follows the phylopath interpretation: a significant Fisher-C p-value rejects a causal model; CICc and weights must be finite; and supported models within delta CICc <= 2 form a competitive set rather than a unique winner.",
        "",
        f"Primary method source: [van der Bijl 2018]({PHYLOPATH_PAPER}).",
        "",
        "## Result",
        "",
        f"- Families audited: {len(audit)}.",
        f"- Families passing the narrow unique-top ranking gate: {int(audit.winner_claim_allowed_by_ranking_gate.sum())}.",
        f"- Competitive supported sets requiring model-set or averaging language: {int(audit.model_averaging_allowed_by_ranking_gate.sum())}.",
        f"- Blocked families: {int(audit.release_status.str.startswith('blocked_').sum())}.",
        "- Publication winner claims currently allowed: 0; every stored ranking uses historical/pre-audit inputs.",
    ]
    for status, count in status_counts.items():
        lines.append(f"- `{status}`: {count}.")
    lines.extend(
        [
            "",
            "## Family-level decisions",
            "",
            "| Family | n | Stored top | Fisher-C p | Status | Unique-top ranking gate |",
            "|---|---:|---|---:|---|---|",
        ]
    )
    for row in audit.itertuples(index=False):
        lines.append(
            f"| `{row.family}` | {row.n_species} | `{row.stored_top_model}` | "
            f"{row.stored_top_fisher_c_p:.6g} | `{row.release_status}` | "
            f"{'yes' if row.winner_claim_allowed_by_ranking_gate else 'no'} |"
        )
    lines.extend(
        [
            "",
            "The final column means only that the stored table passes the Fisher-C/finite-CICc/delta-CICc ranking gate. It does not authorize a publication claim from historical inputs. A lower-ranked supported model is reported only as a diagnostic when the stored top model is rejected; this audit does not automatically promote it. Candidate-family definitions and corrected inputs must be frozen before rerunning phylopath.",
            "",
        ]
    )
    return "\n".join(lines)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results-dir", type=Path, default=RESULTS_DIR)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    args = parser.parse_args()

    rankings = sorted(args.results_dir.glob("*_model_ranking.csv"))
    if not rankings:
        raise FileNotFoundError(f"No model-ranking CSVs found in {args.results_dir}")
    rows: List[Dict[str, object]] = []
    source_files = []
    for path in rankings:
        family = path.name[: -len("_model_ranking.csv")]
        input_path = args.results_dir / f"{family}_analysis_input.csv"
        if not input_path.exists():
            raise FileNotFoundError(f"Missing analysis input for {family}: {input_path}")
        n_species = len(pd.read_csv(input_path))
        rows.append(audit_ranking(pd.read_csv(path), family, n_species))
        source_files.extend([path, input_path])
    audit = pd.DataFrame(rows).sort_values("family").reset_index(drop=True)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    csv_path = args.output_dir / f"{AUDIT_NAME}.csv"
    report_path = args.output_dir / f"{AUDIT_NAME}.md"
    manifest_path = args.output_dir / f"{AUDIT_NAME}.manifest.json"
    audit.to_csv(csv_path, index=False)
    report_path.write_text(_render_markdown(audit), encoding="utf-8")
    manifest = {
        "analysis_state": "historical_pre_audit_outputs",
        "mutated_historical_results": False,
        "alpha": 0.05,
        "competitive_delta_cicc": 2.0,
        "n_families": len(audit),
        "n_ranking_gate_unique_top_pass": int(
            audit["winner_claim_allowed_by_ranking_gate"].sum()
        ),
        "n_competitive_sets": int(
            audit["model_averaging_allowed_by_ranking_gate"].sum()
        ),
        "n_publication_claims_allowed": 0,
        "n_blocked_families": int(
            audit["release_status"].str.startswith("blocked_").sum()
        ),
        "audit_csv": str(csv_path.relative_to(PROJECT_ROOT)),
        "audit_csv_sha256": _sha256(csv_path),
        "report": str(report_path.relative_to(PROJECT_ROOT)),
        "report_sha256": _sha256(report_path),
        "method_source": PHYLOPATH_PAPER,
        "source_files": [
            {
                "path": str(path.relative_to(PROJECT_ROOT)),
                "sha256": _sha256(path),
            }
            for path in source_files
        ],
    }
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(f"Wrote {csv_path}")
    print(f"Wrote {report_path}")
    print(f"Wrote {manifest_path}")


if __name__ == "__main__":
    main()
