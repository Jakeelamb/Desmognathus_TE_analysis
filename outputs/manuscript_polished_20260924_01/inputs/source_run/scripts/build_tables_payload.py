"""Assemble publication tables from copied inputs and the unchanged diversity fits."""
import csv
import json
from pathlib import Path

RUN = Path(__file__).resolve().parents[1]
DATA = RUN / "inputs/datasets"
OUT = RUN / "tables"

def rows(path):
    with Path(path).open(newline="") as handle:
        return list(csv.DictReader(handle))

def supp(stem):
    matches = list(DATA.glob(f"Supplementary_Data_{stem}_*.csv"))
    assert len(matches) == 1, (stem, matches)
    return rows(matches[0])

def animals(data):
    return {r["specimen_id"] for r in data if r["specimen_id"]}

def number(value):
    return float(value) if value not in (None, "", "NA") else None

morph, iod, traits = supp("S18"), supp("S21"), supp("S04")
primary = [r for r in supp("S14") if r["primary_ltr_element"].lower() == "true"]
overlap = supp("S40")
ltr_overlap = rows(RUN / "inputs/ltr_pgls/species_data.csv")
assert len(morph) == 1152 and len(animals(morph)) == 42
assert len(iod) == 805 and len(animals(iod)) == 50
assert len(animals(morph) | animals(iod)) == 51
assert len(primary) == 380 and len({r['species'] for r in primary}) == 30
assert len(traits) == 24 and len(overlap) == 21 and len(ltr_overlap) == 19
coverage = [
    ["TE composition and divergence", 34, None, None, "Species-level repeat summaries", "S02; S08; S12"],
    ["LTR depth distributions", 30, len(primary), None, "Elements; species medians", "S14–S16"],
    ["Cell and nucleus morphology", 24, len(morph), len(animals(morph)), "Selected cell–nucleus pairs", "S18–S19"],
    ["Nuclear optical density", 24, len(iod), len(animals(iod)), "Nuclei; estimated genome size", "S21–S22"],
    ["TE–genome-size comparison", len(overlap), None, None, "Species-level overlap", "S40"],
    ["LTR–genome-size comparison", len(ltr_overlap), 239, None, "Elements summarized to species medians", "LTR PGLS species data"],
]
assert sum(int(r.get("n_elements", r.get("element_n", 0))) for r in ltr_overlap) in (0,239)
result_rows = []
composition = rows(RUN / "inputs/te_pca_trait_association_tests.csv")
for key, label in [("reproductive_strategy_permanova", "Reproductive strategy and TE composition"),
                   ("size_adjusted_model_reproductive_strategy", "Strategy adjusted for recorded body size")]:
    r = next(r for r in composition if r["test_id"] == key)
    result_rows.append({"block":"Composition", "comparison":label, "n_species":int(r["n_species"]),
        "effect_type":"R²", "estimate":number(r["variance_explained_r2"]), "ci_low":None,"ci_high":None,
        "p_value":number(r["p_value"]), "lambda":None, "source":"te_pca_trait_association_tests.csv", "scale":"CLR / Aitchison composition"})
diversity_path = RUN / "figures/figure_3/analysis/model_summary.csv"
diversity = rows(diversity_path)
assert len(diversity) == 2
for r in diversity:
    metric = r["metric"]
    label = "Shannon entropy" if "shannon" in metric.lower() else "Gini–Simpson index"
    result_rows.append({"block":"Diversity", "comparison":label + " versus estimated genome size", "n_species":int(r["n_species"]),
        "effect_type":"Slope per pg/1C", "estimate":number(r["slope"]), "ci_low":number(r["ci_low"]),"ci_high":number(r["ci_high"]),
        "p_value":number(r["p_value"]),"lambda":number(r["lambda"]),"source":"Figure 3 model_summary.csv", "scale":"Untransformed diversity and pg/1C"})
labels = {"genome_size_vs_cell_area":"Cell area versus estimated genome size", "genome_size_vs_nucleus_area":"Nucleus area versus estimated genome size", "nucleus_area_vs_cell_area":"Cell area versus nucleus area"}
cellular = supp("S26")
for key in labels:
    r = next(r for r in cellular if r["comparison"] == key)
    result_rows.append({"block":"Cellular", "comparison":labels[key], "n_species":int(r["n_species"]),
        "effect_type":"Standardized slope", "estimate":number(r["pgls_standardized_beta"]),"ci_low":number(r["pgls_ci_low"]),"ci_high":number(r["pgls_ci_high"]),
        "p_value":number(r["pgls_p_value"]),"lambda":number(r["pagel_lambda"]),"source":"S26", "scale":"Log10 then standardized"})
r = next(r for r in rows(RUN / "inputs/ltr_pgls/model_summary.csv") if r["model"] == "log_log_lambda")
result_rows.append({"block":"LTR", "comparison":"Median LTR ratio versus estimated genome size", "n_species":int(r["n_species"]),
    "effect_type":"Log–log slope", "estimate":number(r["slope"]),"ci_low":number(r["ci_low"]),"ci_high":number(r["ci_high"]),
    "p_value":number(r["p_value"]),"lambda":number(r["lambda"]),"source":"LTR PGLS model_summary.csv", "scale":"Log10 ratio and Gbp"})

OUT.mkdir(exist_ok=True)
with (OUT / "Table_1_data_coverage.csv").open("w", newline="") as handle:
    writer = csv.writer(handle)
    writer.writerow(["Component","Species","Objects","Animals","Measurement","Source"])
    writer.writerows(coverage)
with (OUT / "Table_2_comparative_associations.csv").open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=list(result_rows[0]))
    writer.writeheader(); writer.writerows(result_rows)
payload = {"coverage": coverage,"associations":result_rows,
    "coverage_notes":"Rows overlap. Objects are elements, selected cell–nucleus pairs, or nuclei as indicated. Morphology and IOD animal sets overlap (51 animals in their union). Species-level joins do not imply matching specimens. The four-domain LTR sensitivity retains 989 elements across the same 30 species.",
    "association_notes":"Composition tests use post hoc species-label permutations without phylogenetic adjustment. Regression rows use Pagel-lambda PGLS. Slopes have different scales, specified in the statistic column. Cellular and LTR intervals are t-based. Genome sizes are IOD-derived estimates; morphology describes selected large erythrocytes. Displayed measurement uncertainty is not propagated through these regressions. Both diversity indices summarize each species’ ten most abundant superfamilies. No multiplicity adjustment has been added in this run."}
(OUT / "tables.json").write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n")
print("Prepared 6 coverage rows and",len(result_rows),"association rows from snapshotted sources")
