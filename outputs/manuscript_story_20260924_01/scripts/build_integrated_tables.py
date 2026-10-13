"""Consolidate frozen result exports; never refit statistics or replace old tables."""
from __future__ import annotations
import csv
import json
import math
from pathlib import Path

RUN = Path(__file__).resolve().parents[1]
OUT = RUN / "tables/integrated"
OUT.mkdir(parents=True, exist_ok=True)

def records(path):
    with (RUN / path).open() as handle:
        return list(csv.DictReader(handle))

def num(value):
    if value in (None, "", "NA", "NaN"):
        return None
    x = float(value)
    assert math.isfinite(x)
    return int(x) if x.is_integer() else x

def clean(value):
    return str(value).replace("estimated genome size", "genome size").replace("Estimated genome size", "Genome size")

def sheet(name, title, headers, rows, widths, formats, caption, filename):
    assert rows and all(len(row) == len(headers) for row in rows)
    assert len(widths) == len(headers)
    with (OUT / filename).open("w") as handle:
        writer = csv.writer(handle)
        writer.writerow(headers)
        writer.writerows(rows)
    return dict(name=name, title=title, headers=headers, rows=rows, widths=widths,
                formats=formats, caption=caption, csv=filename)

old = json.loads((RUN / "tables/tables.json").read_text())
coverage = [[clean(x) if isinstance(x, str) else x for x in r] for r in old["coverage"]]
coverage[3][4] = "Nuclei; IOD-derived genome size"
coverage += [["Maximum snout–vent length", 37, None, None, "Species maxima; 24 in cellular analyses", "Frozen organismal traits"],
             ["Integrated TE and cellular traits", 21, None, None, "Exact overlap with SVL and phylogeny", "Figure 11 species join"]]
coverage_note = ("Rows overlap. Morphology: 1,152 selected pairs from 42 animals and 43 images; IOD: 805 nuclei from 50 animals and 51 images. "
                 "The assays share 41 animals (51 in their union), but only 180 exact object identifiers. Seven morphology species and four IOD species have one animal. "
                 "Genome size derives from IOD and a fixed D. fuscus assembly anchor. Four-domain LTR sensitivity retains 989 elements across the same 30 species.")
sheets = [sheet("Data coverage", "Table 1. Data coverage", ["Component", "Species", "Objects", "Animals", "Measurement", "Source"],
                coverage, [38, 10, 12, 11, 48, 30], {"1": "0", "2": "#,##0", "3": "0"}, coverage_note, "Table_1_data_coverage.csv")]

associations = []
for r in old["associations"]:
    associations.append([clean(r["comparison"]), r["n_species"], r["effect_type"], r["estimate"], r["ci_low"],
                         r["ci_high"], 0, r["p_value"], r["lambda"], r["scale"], r["source"]])
body_focal = {("body_to_genome", "b"): "Genome size versus SVL", ("body_to_cell", "b"): "Cell area versus SVL",
              ("body_to_nucleus", "b"): "Nucleus area versus SVL", ("genome_to_cell_body_adjusted", "g"): "Cell area versus genome size, adjusted for SVL",
              ("cell_to_nucleus_body_adjusted", "c"): "Nucleus area versus cell area, adjusted for SVL",
              ("genome_and_cell_to_nucleus_body_adjusted", "g"): "Nucleus area versus genome size, adjusted for cell area and SVL"}
for r in records("figures/figure_8/analysis/regression_coefficients.csv"):
    key = r["relation_id"], r["term"]
    if r["panel"] == "full24" and key in body_focal:
        associations.append([body_focal[key], 24, "Standardized slope", *[num(r[k]) for k in ["standardized_coefficient", "ci_low", "ci_high"]],
                             0, num(r["p_value"]), num(r["fitted_phylogenetic_parameter"]), "Full24-standardized log10 traits", "Figure 8 regression coefficients"])
for r in records("figures/figure_10/analysis/allometry_model_summary.csv"):
    if r["analysis_id"] == "released24":
        associations.append(["Nucleus-on-cell area scaling (test against isometry)", 24, "Log–log slope", num(r["slope"]), num(r["ci_low"]),
                             num(r["ci_high"]), 1, num(r["p_vs_isometry"]), num(r["lambda"]), "Unstandardized log10 areas; P tests slope = 1", "Figure 10 allometry summary"])
te_names = {"genome_body_adjusted": "Genome size versus within-LTR score, adjusted for SVL",
            "cell_te_adjusted": "Cell area versus within-LTR score, adjusted for genome size and SVL",
            "nucleus_te_adjusted": "Nucleus area versus within-LTR score, adjusted for cell area, genome size and SVL"}
te_selected = []
for r in records("figures/figure_11/analysis/regression_coefficients.csv"):
    if r["endpoint"] == "within_ltr" and r["term"] == "t" and r["fit_id"] in te_names:
        te_selected.append(r["fit_id"])
        associations.append([te_names[r["fit_id"]], 21, "Standardized slope", *[num(r[k]) for k in ["estimate", "ci_low", "ci_high"]],
                             0, num(r["p_value"]), num(r["lambda"]), "Logit TE share; full24 cellular scaling and full34 TE scaling", "Figure 11 regression coefficients"])
assert len(te_selected) == 3, te_selected
for r in records("figures/figure_12/analysis/strategy_omnibus_tests.csv"):
    if r["endpoint"] == "within_ltr" and r["svl_adjusted"] == "TRUE":
        associations.append(["Within-LTR low-divergence share and reproductive strategy, adjusted for SVL", 34, "Likelihood ratio (df = 2)",
                             num(r["likelihood_ratio"]), None, None, None, num(r["p_value"]), None,
                             "Logit within-LTR fraction; approximate chi-square ML test", "Figure 12 omnibus tests"])
assoc_note = ("Intervals are 95% coefficient intervals, conditional on fitted phylogenetic parameters. Null specifies the coefficient tested (0 or 1); it is not used for the omnibus likelihood-ratio test. "
              "Composition PERMANOVA rows use ordinary label permutations, not phylogenetic correction. All P values are exploratory and unadjusted for multiplicity. "
              "Within-LTR score is the fraction of LTR aligned-hit bp in divergence bins 0–4, logit transformed. This is not whole-genome TE occupancy or an insertion rate.")
sheets.append(sheet("Focal associations", "Table 2. Focal comparative associations", ["Comparison", "Species", "Statistic", "Estimate", "95% CI low", "95% CI high", "Null", "P", "Lambda", "Scale", "Source"],
                    associations, [56, 9, 24, 12, 13, 13, 9, 13, 11, 47, 34], {"1":"0", "3":"0.0000", "4":"0.0000", "5":"0.0000", "6":"0", "7":"0.00E+00", "8":"0.00E+00"}, assoc_note, "Table_2_focal_associations.csv"))

path_rows = []
for r in records("figures/figure_8/analysis/path_model_comparison.csv"):
    if r["candidate_set"] == "body_adjusted_four":
        path_rows.append(["Four SVL-upstream hypotheses", r["panel"], r["model"], num(r["n_species"]), num(r["rank"]), num(r["CICc"]), num(r["delta_CICc"]), num(r["w"]), num(r["p"]), "Figure 8"])
for r in records("figures/figure_13/analysis/model_comparison.csv"):
    if r["test_order"] == "primary":
        path_rows.append(["Eight body-role hypotheses", r["panel"], r["model"], num(r["n_species"]), num(r["rank"]), num(r["CICc"]), num(r["delta_CICc"]), num(r["w"]), num(r["p"]), "Figure 13; primary test order"])
assert len(path_rows) == 24
for candidate_set in {r[0] for r in path_rows}:
    for cohort in ("full24", "subset22"):
        assert abs(sum(r[7] for r in path_rows if r[0] == candidate_set and r[1] == cohort) - 1) < 1e-10
sheets.append(sheet("Path comparisons", "Table 3. Path-model comparison", ["Candidate set", "Species panel", "Model", "Species", "Rank", "CICc", "Delta CICc", "Weight", "Global-fit P", "Source"],
                    path_rows, [36, 16, 37, 9, 9, 12, 13, 12, 14, 36], {"3":"0", "4":"0", "5":"0.000", "6":"0.000", "7":"0.00E+00", "8":"0.00E+00"},
                    "Compare weights and delta CICc only within the same candidate set and species panel. Subset22 excludes aeneus and wrighti. Weights are relative support among specified hypotheses, not causal probabilities. Reverse test-order scores and equivalent orientations are available with Figure 13.", "Table_3_path_comparisons.csv"))

traits = {r["tree_tip"]: r for r in records("figures/figure_8/analysis/traits_joined.csv")}
organismal = {r["species"]: r for r in records("inputs/body_size/organismal_traits.csv")}
te = {r["species"]: r for r in records("figures/figure_12/analysis/repeat_score_species34.csv")}
support = {(r["tree_tip"], r["assay"]): r for r in records("figures/figure_9/analysis/species_support.csv")}
species = sorted(set(traits) | set(te) | set(organismal))
assert len(species) == 37 and len(set(traits) & set(te)) == 21
species_rows = []
for sp in species:
    t, o, repeat = traits.get(sp, {}), organismal.get(sp, {}), te.get(sp, {})
    m, i = support.get((sp, "Paired morphology"), {}), support.get((sp, "Nuclear IOD"), {})
    species_rows.append(["D. " + sp, num(o.get("max_svl_mm")), num(t.get("genome_size_pg_fuscus_anchored")), num(t.get("cell_area_um2")), num(t.get("nucleus_area_um2")),
                         num(m.get("n_animals")), num(m.get("n_objects")), num(m.get("max_animal_fraction")), num(i.get("n_animals")), num(i.get("n_objects")),
                         num(repeat.get("within_ltr_fraction")), o.get("development_mode") or None, "Yes" if sp in traits and sp in te else "No"])
assert sum(r[6] or 0 for r in species_rows) == 1152 and sum(r[9] or 0 for r in species_rows) == 805
sheets.append(sheet("Species measurements", "Table 4. Species measurements and sampling", ["Species", "Max SVL (mm)", "Genome size (pg/1C)", "Cell area (µm²)", "Nucleus area (µm²)", "Morphology animals", "Paired objects", "Largest animal share", "IOD animals", "IOD nuclei", "Low-divergence share within LTR", "Development mode", "Integrated21"],
                    species_rows, [25, 14, 20, 18, 20, 18, 15, 19, 14, 14, 26, 24, 16], {"1":"0", "2":"0.00", "3":"0.00", "4":"0.00", "5":"0", "6":"0", "7":"0.0%", "8":"0", "9":"0", "10":"0.0%"},
                    "Missing measurements are n.a. Species are alphabetized for lookup. Sources: frozen organismal traits; S04 species traits; Figure 9 animal support; Figure 12 repeat scores. Cell and nucleus areas are the retained-object species medians. Largest animal share is the fraction of paired morphology objects from the most represented animal. Genome size is IOD-derived using a fixed D. fuscus anchor.", "Table_4_species_measurements.csv"))

robust = []
focal_keys = {("genome_to_cell_body_adjusted", "g"): "Genome–cell | SVL", ("cell_to_nucleus_body_adjusted", "c"): "Cell–nucleus | SVL", ("genome_and_cell_to_nucleus_body_adjusted", "g"): "Genome–nucleus | cell + SVL"}
for r in records("figures/figure_9/analysis/point_coefficients.csv"):
    key = r["relation_id"], r["term"]
    if key in focal_keys:
        robust.append([focal_keys[key], r["aggregation"] + "; " + r["panel"], num(r["n_species"]), num(r["estimate"]), num(r["ci_low"]), num(r["ci_high"]), "95% t coefficient interval", num(r["p_value"]), "Figure 9"])
for r in records("figures/figure_9/analysis/coefficient_sensitivity_summary.csv"):
    key = r["relation_id"], r["term"]
    if key in focal_keys and r["mode"] == "hierarchical_bootstrap":
        robust.append([focal_keys[key], "250 hierarchical measurement resamples", 24, num(r["median"]), num(r["q025"]), num(r["q975"]), "Central 95% measurement-sensitivity range", None, "Figure 9"])
for r in records("figures/figure_10/analysis/allometry_model_summary.csv"):
    robust.append(["Nucleus-on-cell unstandardized log slope", r["analysis_id"], num(r["n_species"]), num(r["slope"]), num(r["ci_low"]), num(r["ci_high"]), "95% t interval; P tests slope = 1", num(r["p_vs_isometry"]), "Figure 10"])
assert len(robust) == 19
sheets.append(sheet("Sampling sensitivity", "Table 5. Sampling and allometry sensitivity", ["Relationship", "Analysis", "Species", "Estimate", "Lower", "Upper", "Interval / test", "P", "Source"],
                    robust, [44, 44, 9, 12, 12, 12, 48, 13, 20], {"2":"0", "3":"0.000", "4":"0.000", "5":"0.000", "7":"0.00E+00"},
                    "Figure 9 balances morphology by the median of animal medians; Figure 10 uses the mean of animal medians. Hierarchical resampling preserves cell–nucleus pairing. Assays are resampled independently in Figure 9; single-animal species cannot supply between-animal variance. Measurement-sensitivity ranges condition on sampled animals, selected cells and the fixed tree. Allometry is conditional nucleus-on-cell regression, not an errors-in-variables estimate.", "Table_5_sampling_sensitivity.csv"))

(OUT / "workbook_payload.json").write_text(json.dumps({"sheets": sheets}, indent=2, allow_nan=False) + "\n")
(RUN / "qa/integrated_tables_checks.json").write_text(json.dumps({"passed": True, "sheets": len(sheets), "rows": [len(s["rows"]) for s in sheets], "species_union": len(species), "integrated_overlap": 21, "paired_objects": 1152, "iod_objects": 805, "no_statistics_refitted": True}, indent=2) + "\n")
print(f"Prepared {len(sheets)} tables with {[len(s['rows']) for s in sheets]} rows; previous tables retained.")
