"""Shared, source-explicit transposable-element classification helpers."""

from __future__ import annotations

from typing import Any, Mapping, Tuple

import pandas as pd

from .dnaPipe import CLASSIFICATION_MAP


Classification = Tuple[Any, Any, Any]

DNAPIPETE_CONTEXT_RENAME = {
    "Source": "SRX_ID",
    "RM_classification": "dnapipete_repeat_class",
    "Class": "dnapipete_te_class",
    "Order": "dnapipete_order",
    "Superfamily": "dnapipete_superfamily",
}


def _present(value: Any) -> bool:
    return isinstance(value, str) and bool(value.strip())


def classify_repeatmasker_hits(
    hits: pd.DataFrame,
    classification_map: Mapping[str, Classification] = CLASSIFICATION_MAP,
) -> pd.DataFrame:
    """Classify every RepeatMasker hit from that hit's own ``repeat_class``.

    Existing columns are retained unchanged so separately sourced dnaPipeTE
    annotations can remain available for disagreement audits.
    """

    if "repeat_class" not in hits.columns:
        raise ValueError("RepeatMasker hit table must contain 'repeat_class'")

    classified = hits.copy()
    mapped = classified["repeat_class"].map(classification_map)

    te_class = mapped.map(
        lambda value: value[0] if isinstance(value, tuple) and _present(value[0]) else pd.NA
    )
    order = mapped.map(
        lambda value: value[1] if isinstance(value, tuple) and _present(value[1]) else pd.NA
    )
    superfamily = mapped.map(
        lambda value: value[2] if isinstance(value, tuple) and _present(value[2]) else pd.NA
    )

    status = pd.Series("unmapped", index=classified.index, dtype="object")
    mapped_label = mapped.notna()
    status.loc[mapped_label & te_class.notna()] = "class_only"
    status.loc[mapped_label & te_class.notna() & order.notna()] = "order_only"
    status.loc[
        mapped_label & te_class.notna() & order.notna() & superfamily.notna()
    ] = "classified"
    status.loc[mapped_label & te_class.isna()] = "mapped_unresolved"

    classified["repeatmasker_te_class"] = te_class.fillna("Unclassified")
    classified["repeatmasker_order"] = order.fillna("Unclassified")
    classified["repeatmasker_superfamily"] = superfamily.fillna("Unclassified")
    classified["repeatmasker_classification_status"] = status
    # Existing downstream readers use these canonical names. They are aliases
    # of the hit-level RepeatMasker result; dnaPipeTE values remain namespaced.
    classified["Class"] = classified["repeatmasker_te_class"]
    classified["Order"] = classified["repeatmasker_order"]
    classified["Superfamily"] = classified["repeatmasker_superfamily"]

    return classified


def prepare_dnapipete_context(dnapipete: pd.DataFrame) -> pd.DataFrame:
    """Namespace and reduce dnaPipeTE annotations to one row per sample/contig."""

    required = {
        "dnaPipeTE_contig_name",
        "Source",
        "RM_classification",
        "Class",
        "Order",
        "Superfamily",
        "hitlength_contiglength",
    }
    missing = sorted(required.difference(dnapipete.columns))
    if missing:
        raise ValueError(f"dnaPipeTE context is missing columns: {missing}")

    context = dnapipete.copy().rename(columns=DNAPIPETE_CONTEXT_RENAME)
    context["hitlength_contiglength"] = pd.to_numeric(
        context["hitlength_contiglength"], errors="coerce"
    )
    context = context.sort_values(
        ["SRX_ID", "dnaPipeTE_contig_name", "hitlength_contiglength"],
        ascending=[True, True, False],
        kind="mergesort",
        na_position="last",
    ).drop_duplicates(["dnaPipeTE_contig_name", "SRX_ID"], keep="first")
    if context.duplicated(["dnaPipeTE_contig_name", "SRX_ID"]).any():
        raise RuntimeError("dnaPipeTE context is not unique by sample and contig")
    return context.reset_index(drop=True).set_index("SRX_ID", drop=False)


def annotate_repeatmasker_hits(
    hits: pd.DataFrame, dnapipete_context: pd.DataFrame
) -> pd.DataFrame:
    """Attach contig context while keeping native hit classification authoritative."""

    required_hits = {"query_name", "SRX_ID", "repeat_class"}
    missing_hits = sorted(required_hits.difference(hits.columns))
    if missing_hits:
        raise ValueError(f"RepeatMasker hits are missing columns: {missing_hits}")
    required_context = {"dnaPipeTE_contig_name", "SRX_ID"}
    missing_context = sorted(required_context.difference(dnapipete_context.columns))
    if missing_context:
        raise ValueError(f"dnaPipeTE context is missing columns: {missing_context}")

    sample_ids = hits["SRX_ID"].dropna().unique()
    available_ids = [sample_id for sample_id in sample_ids if sample_id in dnapipete_context.index]
    if available_ids:
        relevant_context = dnapipete_context.loc[available_ids].reset_index(drop=True)
    else:
        relevant_context = dnapipete_context.iloc[0:0].reset_index(drop=True)
    annotated = pd.merge(
        hits,
        relevant_context,
        left_on=["query_name", "SRX_ID"],
        right_on=["dnaPipeTE_contig_name", "SRX_ID"],
        how="left",
        suffixes=("", "_dnapipete"),
        validate="many_to_one",
        sort=False,
    )
    if len(annotated) != len(hits):
        raise RuntimeError(
            "RepeatMasker/dnaPipeTE merge changed hit count: "
            f"{len(hits)} input rows became {len(annotated)} rows"
        )
    annotated.index = hits.index
    return classify_repeatmasker_hits(annotated)
