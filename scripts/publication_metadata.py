"""Curated metadata for manuscript-critical supplementary datasets."""

SUPPLEMENT_DEPENDENCIES: dict[str, tuple[str, ...]] = {
    "S01": (
        "analyses/01_transposable_elements/data/ltr_resource_coverage.csv",
        "analyses/01_transposable_elements/recompute_ltr_release.py",
        "analyses/01_transposable_elements/provenance/ltr_terminal_internal_release.json",
        "analyses/02_morphology/data/relative_nuclear_iod_species.csv",
        "analyses/02_morphology/validate_release.py",
        "analyses/02_morphology/provenance/iod_source_manifest.json",
    ),
    "S04": (
        "analyses/02_morphology/data/relative_nuclear_iod_species.csv",
        "analyses/02_morphology/provenance/iod_source_manifest.json",
    ),
    "S08": (
        "analyses/01_transposable_elements/data/te_superfamily_composition.csv",
        "analyses/01_transposable_elements/recompute_pca.py",
    ),
    "S14": (
        "analyses/01_transposable_elements/recompute_ltr_release.py",
        "analyses/01_transposable_elements/provenance/ltr_terminal_internal_release.json",
    ),
    "S15": (
        "analyses/01_transposable_elements/data/ltr_element_metrics.csv",
        "analyses/01_transposable_elements/recompute_ltr_release.py",
        "analyses/01_transposable_elements/provenance/ltr_terminal_internal_release.json",
    ),
    "S16": (
        "analyses/01_transposable_elements/data/ltr_element_metrics.csv",
        "analyses/01_transposable_elements/recompute_ltr_release.py",
        "analyses/01_transposable_elements/provenance/ltr_terminal_internal_release.json",
    ),
    "S17": (
        "analyses/01_transposable_elements/provenance/ltr_terminal_internal_release.json",
    ),
    "S19": (
        "analyses/02_morphology/data/cell_nucleus_objects.csv",
        "analyses/02_morphology/provenance/review_finalization.json",
    ),
    "S22": (
        "analyses/02_morphology/data/nuclear_iod_objects.csv",
        "analyses/02_morphology/data/nuclear_iod_by_image.csv",
        "analyses/02_morphology/validate_release.py",
        "analyses/02_morphology/provenance/iod_source_manifest.json",
        "analyses/02_morphology/provenance/review_finalization.json",
    ),
    "S23": (
        "analyses/02_morphology/data/nuclear_iod_objects.csv",
        "analyses/02_morphology/provenance/review_finalization.json",
    ),
    "S24": (
        "analyses/02_morphology/data/nuclear_iod_objects.csv",
        "analyses/02_morphology/provenance/iod_source_manifest.json",
        "analyses/02_morphology/provenance/review_finalization.json",
    ),
    "S25": (
        "analyses/02_morphology/data/nuclear_iod_objects.csv",
        "analyses/02_morphology/provenance/iod_source_manifest.json",
    ),
    **{
        supplement_id: (
            "analyses/03_phylogenetic_path/data/path24_traits.csv",
            "analyses/03_phylogenetic_path/run_path_analysis.R",
            "analyses/02_morphology/data/relative_nuclear_iod_species.csv",
            "analyses/02_morphology/provenance/iod_source_manifest.json",
        )
        for supplement_id in (
            "S26",
            "S27",
            "S28",
            "S29",
            "S30",
            "S31",
            "S32",
            "S33",
            "S34",
            "S35",
            "S36",
            "S37",
        )
    },
    "S40": (
        "data/identity/te34_panel.csv",
        "data/identity/path24_panel.csv",
        "analyses/01_transposable_elements/data/te_diversity.csv",
        "analyses/02_morphology/data/relative_nuclear_iod_species.csv",
        "analyses/02_morphology/provenance/iod_source_manifest.json",
        "data/identity/species_taxonomy_crosswalk.csv",
    ),
    "S41": (
        "analyses/01_transposable_elements/data/wang2025_figure2_te_diversity.csv",
        "analyses/01_transposable_elements/data/wang2025_figure2_te_diversity.provenance.json",
        "analyses/01_transposable_elements/build_cross_study_diversity.py",
        "analyses/01_transposable_elements/data/te_diversity.csv",
        "analyses/02_morphology/data/relative_nuclear_iod_species.csv",
        "data/identity/te34_panel.csv",
        "data/identity/path24_panel.csv",
    ),
    "S42": (
        "data/identity/te34_panel.csv",
        "analyses/03_phylogenetic_path/data/organismal_traits.csv",
        "analyses/01_transposable_elements/bonsai_comparison/validate_strategy_localization.py",
        "analyses/01_transposable_elements/bonsai_comparison/data/te34_bonsai_strategy_validation_manifest.json",
    ),
    "S43": (
        "analyses/01_transposable_elements/bonsai_comparison/data/te34_bonsai_strategy_label_permutation_tests.csv",
        "analyses/01_transposable_elements/bonsai_comparison/data/te34_bonsai_feature_bootstrap_strategy_stability.csv",
        "analyses/01_transposable_elements/bonsai_comparison/validate_strategy_localization.py",
        "analyses/01_transposable_elements/bonsai_comparison/data/te34_bonsai_strategy_validation_manifest.json",
    ),
}


SUPPLEMENT_METADATA: dict[str, dict[str, str]] = {
    "S01": {
        "description": (
            "Species-level analysis availability, including the exact S22 projection of "
            "the conditional D. fuscus scale derived from the author-selected rounded "
            "16.1-Gbp pooled assembly span at 0.978 Gbp per pg."
        ),
        "manuscript_use": (
            "Availability and audit table only; its projected picogram field is an "
            "assembly-derived sensitivity, not a measured C-value or independently "
            "calibrated absolute genome size."
        ),
    },
    "S04": {
        "title": "Finalized Path24 conditional genome-size traits",
        "description": (
            "Conditional D. fuscus assembly-anchored genome size in pg/1C, nucleus area, "
            "cell area, and support counts for all 24 finalized species."
        ),
        "manuscript_use": (
            "Primary three-trait input. Genome size is the default analysis trait but is "
            "IOD-derived and assembly-anchored, not an independent assay or measured C-value."
        ),
    },
    "S08": {
        "description": (
            "A 34-row, four-column Wang-style TE-diversity table. For each species, the "
            "ten most abundant positive classified superfamilies were selected, reclosed "
            "to one, and summarized with natural-log Shannon and Gini-Simpson indices."
        ),
        "manuscript_use": (
            "Canonical paper-facing Shannon entropy and Gini-Simpson endpoints, matched "
            "to the Wang et al. top-10-superfamily scope. The species-specific top ten "
            "are explicitly reclosed in this study; Wang et al. did not document that "
            "normalization detail. Unresolved aligned-base mass remains only in S06-S07 QC."
        ),
    },
    "S14": {
        "description": (
            "Non-destructive zero-aware audit of 1,086 historically selected usable LTR "
            "elements, including corrected TEsorter domain counts, current LTR/Gypsy "
            "eligibility, regional coverage, IQR decisions, and primary membership."
        ),
    },
    "S15": {
        "description": (
            "Species summaries across declared element branches, including the primary "
            "two-sided IQR-filtered regional-positive-depth branch and robust estimator "
            "diagnostics."
        ),
    },
    "S16": {
        "description": (
            "Per-species LTR resource audit separating historical pipe-split counts from "
            "correct LTR/Gypsy five-domain eligibility, IQR retention, coverage, and primary "
            "support."
        ),
    },
    "S17": {
        "description": (
            "Two source-corrupt rows from the historical pipe-split cohort; both have four "
            "real TEsorter domain annotations and are audit-only under the corrected gate."
        ),
    },
    "S19": {
        "description": (
            "Upper-tail species medians, conditional paired-object bootstrap intervals, and "
            "animal/image support from all 24 finalized reviewed panels."
        ),
        "manuscript_use": (
            "Primary descriptive morphology table; intervals are conditional on observed "
            "animals and images and are not population-level uncertainty."
        ),
    },
    "S21": {
        "description": (
            "All 805 reviewed nuclei retained across the 24-species image-quality-matched "
            "panel as the object-level IOD measurement and calibration substrate."
        ),
        "manuscript_use": (
            "Object-level IOD measurement and calibration audit; not a separate downstream "
            "relative-IOD analysis branch."
        ),
    },
    "S22": {
        "description": (
            "Equal-image means of within-image median nuclear IOD, the retained Path24-median "
            "relative audit index, and the conditional D. fuscus genome-size scale used by "
            "downstream analyses, derived from the selected rounded 16.1-Gbp pooled assembly "
            "span at 0.978 Gbp per pg."
        ),
        "manuscript_use": (
            "IOD fields retain measurement and calibration provenance. Conditional picogram "
            "fields supply the default downstream genome-size trait, but are not a measured "
            "C-value; their intervals exclude uncertainty in the chosen assembly anchor."
        ),
    },
    "S23": {
        "title": "Unnormalized image-level nuclear IOD",
        "description": (
            "Unnormalized image-specific nuclear-IOD summaries used to construct and audit "
            "the equal-image species estimator."
        ),
        "manuscript_use": (
            "Image-level IOD measurement and aggregation audit; values are neither ratios nor "
            "independent genome-size measurements."
        ),
    },
    "S24": {
        "description": (
            "Species-level support and balance against the apalachicolae quality template; "
            "the declared gate is absolute SMD <= 0.10 and KS distance <= 0.25."
        ),
        "manuscript_use": "Underlying IOD support and image-quality audit for genome-size use.",
    },
    "S25": {
        "description": (
            "Within-species-centered associations between image-quality features and log "
            "nuclear IOD, retained as descriptive confounding diagnostics."
        ),
        "manuscript_use": "Underlying IOD assay-quality audit for genome-size use.",
    },
    "S26": {
        "title": "Conditional genome-size pairwise PGLS sensitivity",
        "description": (
            "Raw correlations and Pagel-lambda PGLS fits among conditional D. fuscus "
            "assembly-anchored genome size, nucleus area, and cell area."
        ),
        "manuscript_use": (
            "Conditional genome-size association sensitivity; the trait is IOD-derived and "
            "not an independent genome-size assay."
        ),
    },
    "S27": {
        "title": "Standardized Path24 conditional genome-size traits",
        "description": (
            "Log10-transformed and standardized conditional genome size, nucleus area, and "
            "cell area for the 24-species primary fit."
        ),
        "manuscript_use": (
            "Primary path-analysis trait matrix on the conditional D. fuscus assembly-anchored "
            "genome-size scale; not an independent assay."
        ),
    },
    "S28": {
        "title": "Pairwise conditional-genome-size phylogenetic regressions",
        "description": (
            "Scale-invariant pairwise phylogenetic regressions among conditional genome size, "
            "nucleus area, and cell area."
        ),
        "manuscript_use": (
            "Exploratory associations on the IOD-derived conditional genome-size scale."
        ),
    },
    "S29": {
        "title": "Conditional genome-size evolutionary-model sensitivity",
        "description": (
            "Conditional genome-size path-model rankings under the declared phylogenetic "
            "evolutionary models."
        ),
        "manuscript_use": (
            "Evolutionary-model sensitivity audit; genome size is assembly-anchored and "
            "IOD-derived rather than independently assayed."
        ),
    },
    "S30": {
        "title": "Conditional genome-size path-model comparison",
        "description": (
            "Global fit, CICc, delta, and weights for the finalized 24-species conditional "
            "genome-size analysis."
        ),
        "manuscript_use": (
            "Exploratory association structures only; conditional genome size is not an "
            "independent assay and equivalence classes do not orient a unique causal arrow."
        ),
    },
    "S31": {
        "title": "Conditional genome-size mechanism edge estimates",
        "description": (
            "Standardized coefficients and intervals for the three Markov-equivalent "
            "conditional-genome-size mechanisms."
        ),
        "manuscript_use": (
            "Exploratory edge estimates; no unique causal orientation or independent "
            "genome-size assay claim."
        ),
    },
    "S32": {
        "title": "Conditional genome-size path-model stability",
        "description": (
            "Top-model and competitive-model frequencies across fitted conditional "
            "genome-size sensitivity analyses."
        ),
        "manuscript_use": "Exploratory stability summary on the IOD-derived conditional scale.",
    },
    "S33": {
        "title": "Conditional genome-size path tree sensitivity",
        "description": (
            "Conditional genome-size model rankings across the declared focal and published-tree "
            "sensitivity fits."
        ),
        "manuscript_use": (
            "Tree uncertainty audit; the primary focal-tree fit retains all 24 species and "
            "uses the IOD-derived conditional genome-size trait."
        ),
    },
    "S34": {
        "title": "Conditional genome-size path measurement sensitivity",
        "description": (
            "Conditional genome-size model rankings across 250 measurement-bootstrap replicates "
            "of the finalized 24-species panel."
        ),
        "manuscript_use": "Measurement-model audit for the conditional genome-size analysis.",
    },
    "S35": {
        "title": "Conditional genome-size path leave-one-species-out analysis",
        "description": (
            "Conditional genome-size model rankings after omitting each of the 24 finalized "
            "species in turn."
        ),
        "manuscript_use": "Species-influence audit for the conditional genome-size analysis.",
    },
    "S36": {
        "title": "Conditional genome-size path simulation calibration",
        "description": (
            "Recovery and false-selection summaries for the declared conditional genome-size "
            "simulation regimes."
        ),
        "manuscript_use": "Calibration gate; blocks confirmatory or uniquely causal interpretation.",
    },
    "S37": {
        "title": "Primary conditional genome-size d-separation tests",
        "description": (
            "Basis-set tests for the finalized 24-species conditional genome-size candidate "
            "models."
        ),
        "manuscript_use": (
            "Global-fit and identifiability audit; no independent-assay or unique-causal-arrow "
            "claim."
        ),
    },
    "S40": {
        "analysis_section": "Integrated TE and conditional genome size",
        "title": "TE diversity and conditional genome size overlap",
        "description": (
            "Exact 21-species TE34-by-Path24 overlap joining canonical species-specific "
            "top-10-superfamily Shannon entropy and Gini-Simpson values to conditional "
            "D. fuscus assembly-anchored genome size and intervals."
        ),
        "manuscript_use": (
            "Descriptive integrated sensitivity table only; genome size is IOD-derived "
            "rather than an independent assay or validated absolute C-value."
        ),
    },
    "S41": {
        "analysis_section": "Cross-study TE-diversity context",
        "title": "Wang 2025 and current top-10-superfamily TE diversity",
        "description": (
            "Descriptive cross-study plotting table containing the 84 exact CC0 Dryad "
            "Figure 2 rows from Wang et al. 2025 and 34 current Desmognathus rows "
            "recomputed over each species' ten most abundant positive superfamilies."
        ),
        "release_status": "sensitivity_only",
        "manuscript_use": (
            "Auxiliary descriptive context only. Published RepeatModeler/RepeatMasker and "
            "current dnaPipeTE inputs come from different pipelines; no pooled regression, "
            "hypothesis test, or estimator-equivalence claim is permitted. Current genome "
            "sizes are conditional D. fuscus-anchored IOD-derived sensitivities, not direct "
            "C-values."
        ),
    },
    "S42": {
        "analysis_section": "TE-composition tree representation",
        "title": "TE34 Bonsai tip annotations",
        "description": (
            "Exact 34-species tip annotation table for the pinned Bonsai "
            "TE-composition representation, including expert reproductive-strategy "
            "labels for every species."
        ),
        "release_status": "main_candidate",
        "manuscript_use": (
            "Supplies Figure 2 tip labels and post-hoc reproductive-strategy colors. "
            "The representation is not a species phylogeny, and annotation timing "
            "precludes a confirmatory trait-association claim."
        ),
    },
    "S43": {
        "analysis_section": "TE-composition tree validation",
        "title": "TE34 Bonsai reproductive-strategy validation",
        "description": (
            "Three-row publication summary of fixed-tree strategy localization and "
            "200 bounded feature-resampling Bonsai refits for the TE34 representation."
        ),
        "release_status": "main_candidate",
        "manuscript_use": (
            "Supplies Figure 2 validation annotations and exact-edge stability. "
            "Label tests are post hoc, and coordinate resampling is a feature "
            "perturbation rather than a biological-specimen bootstrap."
        ),
    },
}


COMMON_COLUMN_DEFINITIONS: dict[str, str] = {
    "species": "Canonical analyzed Desmognathus species label.",
    "te_sra_accession": "SRA experiment accession bound to the species TE resource.",
    "te_assembly_accession": "Assembly accession bound to the species TE resource.",
    "genome_size_pg_fuscus_anchored": (
        "Conditional target-to-D. fuscus equal-image IOD ratio multiplied by the selected "
        "16.462167689-pg/1C assembly anchor; not an independent genome-size assay."
    ),
    "genome_size_pg_ci_low": (
        "Lower conditional 95% measurement-bootstrap interval on the D. fuscus-anchored "
        "genome-size scale; anchor-value uncertainty is excluded."
    ),
    "genome_size_pg_ci_high": (
        "Upper conditional 95% measurement-bootstrap interval on the D. fuscus-anchored "
        "genome-size scale; anchor-value uncertainty is excluded."
    ),
    "genome_size_definition": (
        "Declared conditional D. fuscus assembly-anchored, IOD-derived genome-size contract; "
        "not a measured C-value."
    ),
    "path_trait_definition": (
        "Declared path trait contract: conditional D. fuscus assembly-anchored genome size "
        "in pg/1C, IOD-derived and not a measured C-value."
    ),
    "n_genome_size_nuclei": (
        "Number of finalized reviewed nuclei supporting the conditional genome-size estimate."
    ),
    "n_genome_size_images": (
        "Number of distinct microscopy images supporting the conditional genome-size estimate."
    ),
    "n_genome_size_specimens": (
        "Number of distinct animals supporting the conditional genome-size estimate."
    ),
    "n_nuclei": "Number of finalized reviewed nuclei contributing to the row.",
    "n_images": "Number of distinct source images contributing to the row.",
    "n_specimens": "Number of distinct specimen_id values; each value denotes one animal.",
}


SUPPLEMENT_COLUMN_DEFINITIONS: dict[str, dict[str, str]] = {
    "S01": {
        "step06_genome_size_estimation": (
            "Whether the species has a reviewed nuclear-IOD estimate for the named "
            "conditional S22 assembly-derived scaling; not validation of absolute genome size."
        ),
        "step06_genome_panel_status": (
            "Release and image-quality support status for the conditional S22 scaling."
        ),
        "step06_genome_nuclei": (
            "Number of finalized reviewed nuclei contributing to the conditional S22 estimate."
        ),
        "step06_genome_images": (
            "Number of distinct source images contributing to the conditional S22 estimate."
        ),
        "step06_genome_size_pg_fuscus_anchored": (
            "Exact S22 projection of the conditional D. fuscus-anchored pg sensitivity "
            "derived from 16.1 Gbp / 0.978 Gbp per pg; not a measured C-value."
        ),
        "step04_ltr_resource_status": (
            "Corrected LTR/Gypsy five-domain resource status projected from S16."
        ),
        "step04_ltr_tabout_elements": "All candidate elements in the source LTRharvest table.",
        "step04_ltr_selected_elements": (
            "Elements passing the corrected LTR/Gypsy and at-least-five-domain gate."
        ),
        "step04_ltr_usable_elements": (
            "Correctly selected elements with usable zero-aware depth records."
        ),
        "step04_ltr_source_corrupt_exclusions": (
            "Source-corrupt exclusions within the correctly selected cohort."
        ),
    },
    "S04": {
        "genome_size_pg_fuscus_anchored": COMMON_COLUMN_DEFINITIONS[
            "genome_size_pg_fuscus_anchored"
        ],
        "genome_size_pg_ci_low": COMMON_COLUMN_DEFINITIONS["genome_size_pg_ci_low"],
        "genome_size_pg_ci_high": COMMON_COLUMN_DEFINITIONS["genome_size_pg_ci_high"],
        "n_genome_size_nuclei": COMMON_COLUMN_DEFINITIONS["n_genome_size_nuclei"],
        "n_genome_size_images": COMMON_COLUMN_DEFINITIONS["n_genome_size_images"],
        "n_genome_size_specimens": COMMON_COLUMN_DEFINITIONS["n_genome_size_specimens"],
        "genome_size_definition": COMMON_COLUMN_DEFINITIONS["genome_size_definition"],
    },
    "S08": {
        "n_superfamilies_included": (
            "Number of species-specific most-abundant positive classified TE "
            "superfamilies retained before reclosure; fixed at ten."
        ),
        "shannon_entropy": (
            "Natural-log Shannon entropy, -sum(p_i * ln(p_i)), over the species-specific "
            "ten most abundant positive superfamilies after reclosure to one."
        ),
        "gini_simpson": (
            "Gini-Simpson index, 1 - sum(p_i^2), over the species-specific ten most "
            "abundant positive superfamilies after reclosure to one."
        ),
    },
    "S14": {
        "element_id": "Exact sequence-and-coordinate LTR candidate identifier.",
        "domain_count": (
            "Correct number of whitespace-delimited DOMAIN|MODEL annotations in the "
            "TEsorter Domains field."
        ),
        "legacy_pipe_split_domain_count": (
            "Historical erroneous count obtained by splitting the Domains field on the "
            "within-token pipe delimiter; retained only for audit."
        ),
        "domain_gate_keep": (
            "True for an exact TEsorter-joined LTR/Gypsy element at least 3,000 bp long "
            "with at least five correctly parsed domain annotations."
        ),
        "domain_gate_rule_id": "Stable identifier for the corrected TEsorter domain gate.",
        "domain_gate_reason": (
            "Machine-readable eligibility or exclusion reason under the corrected domain gate."
        ),
        "depth_positions_expected": "Number of element positions expected from stored coordinates.",
        "depth_positions_reported": "Number of unique positions present in the depth file.",
        "depth_positions_missing": "Expected positions absent from the depth file.",
        "depth_positions_explicit_zero": "Reported positions whose stored read depth equals zero.",
        "left_ltr_positive_coverage_fraction": (
            "Fraction of expected left-LTR positions with positive read depth."
        ),
        "right_ltr_positive_coverage_fraction": (
            "Fraction of expected right-LTR positions with positive read depth."
        ),
        "terminal_positive_coverage_fraction": (
            "Fraction of expected positions across both terminal LTRs with positive read depth."
        ),
        "internal_positive_coverage_fraction": (
            "Fraction of expected internal-region positions with positive read depth."
        ),
        "mean_depth_left_ltr_all_positions": (
            "Mean left-LTR depth across every expected position, assigning absent positions zero."
        ),
        "mean_depth_right_ltr_all_positions": (
            "Mean right-LTR depth across every expected position, assigning absent positions zero."
        ),
        "mean_depth_terminal_all_positions": (
            "Mean depth across both terminal LTRs and every expected terminal position."
        ),
        "mean_depth_internal_all_positions": (
            "Mean depth across every expected internal-region position."
        ),
        "ratio_terminal_internal_all_positions": (
            "Zero-aware mean terminal depth divided by zero-aware mean internal depth."
        ),
        "log2_ratio_terminal_internal_all_positions": (
            "Base-2 logarithm of ratio_terminal_internal_all_positions."
        ),
        "left_right_terminal_log2_imbalance": (
            "Base-2 logarithm of mean left-LTR depth divided by mean right-LTR depth."
        ),
        "depth_file_sha256": "SHA-256 digest of the preserved source depth file.",
        "source_depth_status": "Frozen usability status of the source depth file.",
        "reported_position_fraction": (
            "Reported unique positions divided by expected positions; distinct from positive depth."
        ),
        "coverage_ge_80pct": (
            "True when both terminal regions and the internal region each have at least 80% "
            "positive-depth coverage; not a reporting-completeness flag."
        ),
        "species_ratio_q1": (
            "Within-species first quartile among correctly selected LTR/Gypsy five-domain "
            "elements; blank for audit-only rows."
        ),
        "species_ratio_q3": (
            "Within-species third quartile among correctly selected LTR/Gypsy five-domain "
            "elements; blank for audit-only rows."
        ),
        "species_ratio_iqr": (
            "Within-species interquartile range among correctly selected elements."
        ),
        "species_ratio_iqr_lower_bound": (
            "Inclusive lower Tukey fence, Q1 minus 1.5 times the within-species IQR."
        ),
        "species_ratio_iqr_upper_bound": (
            "Inclusive upper Tukey fence, Q3 plus 1.5 times the within-species IQR."
        ),
        "iqr_filter_keep": (
            "True when an eligible element ratio lies inclusively within its species' "
            "two-sided 1.5-IQR fences."
        ),
        "iqr_filter_direction": (
            "Whether the row is ineligible, retained, or a lower or upper IQR outlier."
        ),
        "iqr_filter_rule_id": "Stable identifier for the required primary IQR rule.",
        "iqr_filter_reason": "Machine-readable retention or exclusion reason under the IQR rule.",
        "primary_ltr_element": (
            "True when both iqr_filter_keep and coverage_ge_80pct are true; this is the "
            "primary LTR element set."
        ),
        "historical_nonzero_only_ratio": (
            "Historical ratio after conditioning regional means on positive-depth positions only."
        ),
        "nonzero_ratio_reproduction_abs_error": (
            "Absolute difference between recomputed and preserved historical nonzero-only ratios."
        ),
    },
    "S15": {
        "analysis_branch": (
            "Declared element subset; primary_iqr_filtered_coverage_ge_80pct is the paper's "
            "primary LTR branch and other branches are audit or sensitivity summaries."
        ),
        "n_elements": "Number of elements contributing to this species and analysis branch.",
        "ratio_arithmetic_mean": "Arithmetic mean of element terminal:internal ratios.",
        "ratio_median": "Median of element terminal:internal ratios; primary displayed estimator.",
        "ratio_geometric_mean": "Geometric mean of positive element terminal:internal ratios.",
        "ratio_trimmed_mean_10pct": "Arithmetic mean after symmetric 10% trimming.",
        "ratio_winsorized_mean_5pct": "Arithmetic mean after symmetric 5% winsorization.",
        "ratio_q05": "Fifth percentile of element terminal:internal ratios.",
        "ratio_q95": "Ninety-fifth percentile of element terminal:internal ratios.",
        "ratio_max": "Maximum element terminal:internal ratio.",
        "fraction_ratio_gt_1": "Fraction of elements with terminal:internal ratio greater than one.",
        "median_bootstrap_ci_low": "2.5th percentile of 2,000 element-bootstrap medians.",
        "median_bootstrap_ci_high": "97.5th percentile of 2,000 element-bootstrap medians.",
        "geometric_mean_bootstrap_ci_low": (
            "2.5th percentile of 2,000 element-bootstrap geometric means."
        ),
        "geometric_mean_bootstrap_ci_high": (
            "97.5th percentile of 2,000 element-bootstrap geometric means."
        ),
        "max_abs_leave_one_out_mean_shift": (
            "Largest absolute arithmetic-mean shift after omitting one element."
        ),
        "max_abs_leave_one_out_median_shift": (
            "Largest absolute median shift after omitting one element."
        ),
        "largest_element_fraction_of_ratio_sum": (
            "Largest single element ratio divided by the species branch ratio sum."
        ),
        "median_terminal_positive_coverage_fraction": (
            "Median terminal positive-depth fraction across contributing elements."
        ),
        "median_internal_positive_coverage_fraction": (
            "Median internal positive-depth fraction across contributing elements."
        ),
    },
    "S16": {
        "n_historical_pipe_split_5_6_domain_elements": (
            "Elements selected by the archived erroneous pipe-fragment count, before the "
            "TEsorter grammar correction."
        ),
        "n_historical_source_corrupt_elements": (
            "Source-corrupt rows from the historical pipe-split cohort."
        ),
        "n_selected_ltr_gypsy_ge5_domain_elements": (
            "Elements passing the corrected LTR/Gypsy and at-least-five-domain gate."
        ),
        "n_usable_elements": (
            "Correctly selected elements with usable zero-aware depth records."
        ),
        "n_selected_source_corrupt_exclusions": (
            "Source-corrupt exclusions within the correctly selected cohort."
        ),
        "n_coverage_ge_80pct_elements": (
            "Correctly selected elements passing the independent 80% positive-depth gate."
        ),
        "n_iqr_retained_elements": (
            "Correctly selected elements retained by the within-species inclusive two-sided "
            "1.5-IQR rule."
        ),
        "n_primary_iqr_filtered_elements": (
            "Elements retained by both the required IQR rule and the 80% regional "
            "positive-depth gate."
        ),
    },
    "S17": {
        "domain_count": (
            "Correct whitespace-delimited TEsorter domain count; both historical corrupt "
            "rows have four annotations."
        ),
        "legacy_pipe_split_domain_count": (
            "Historical erroneous pipe-fragment count, equal to five for both audit rows."
        ),
        "domain_gate_keep": "False because neither historical corrupt row has five domains.",
        "domain_gate_rule_id": "Stable identifier for the corrected TEsorter domain gate.",
        "domain_gate_reason": "Reason the historical corrupt row is currently ineligible.",
    },
    "S19": {
        "n_size_cells": "Number of finalized reviewed cell-nucleus pairs in the species summary.",
        "n_size_images": "Number of distinct source-image filenames represented.",
        "n_size_specimens": "Number of distinct animals represented.",
        "largest_specimen_n": "Largest number of retained pairs contributed by one animal.",
        "largest_specimen_fraction": (
            "largest_specimen_n divided by n_size_cells; a specimen-concentration diagnostic."
        ),
        "cell_area_um2": "Median cell-mask area among retained upper-tail pairs, in square micrometers.",
        "cell_area_um2_ci_low": "2.5th percentile of the conditional paired-object bootstrap median.",
        "cell_area_um2_ci_high": "97.5th percentile of the conditional paired-object bootstrap median.",
        "nucleus_area_um2": (
            "Median linked nucleus-mask area among retained upper-tail pairs, in square micrometers."
        ),
        "nucleus_area_um2_ci_low": (
            "2.5th percentile of the conditional paired-object bootstrap nucleus-area median."
        ),
        "nucleus_area_um2_ci_high": (
            "97.5th percentile of the conditional paired-object bootstrap nucleus-area median."
        ),
    },
    "S22": {
        "iod_equal_image_estimate": (
            "Arithmetic mean of within-image median per-nucleus integrated optical density."
        ),
        "iod_equal_image_ci_low": "2.5th percentile of the hierarchical image-then-nucleus bootstrap.",
        "iod_equal_image_ci_high": (
            "97.5th percentile of the hierarchical image-then-nucleus bootstrap."
        ),
        "relative_iod_anchor": "Median iod_equal_image_estimate across the complete Path24 panel.",
        "relative_iod_index": "iod_equal_image_estimate divided by relative_iod_anchor.",
        "relative_iod_ci_low": "iod_equal_image_ci_low divided by relative_iod_anchor.",
        "relative_iod_ci_high": "iod_equal_image_ci_high divided by relative_iod_anchor.",
        "iod_ratio_to_fuscus": (
            "Target equal-image IOD divided by the D. fuscus equal-image IOD; a separate "
            "conditional sensitivity ratio, not the primary Path24-median index."
        ),
        "genome_size_pg_fuscus_anchored": (
            "Conditional equal-image IOD ratio to D. fuscus multiplied by the assembly-derived "
            "1C-equivalent anchor (16.1 Gbp / 0.978 Gbp per pg); not a measured C-value or "
            "independently calibrated absolute genome size."
        ),
        "genome_size_pg_ci_low": (
            "2.5th percentile of the hierarchical species-to-fuscus ratio bootstrap on the "
            "conditional pg scale; anchor-value uncertainty is excluded."
        ),
        "genome_size_pg_ci_high": (
            "97.5th percentile of the hierarchical species-to-fuscus ratio bootstrap on the "
            "conditional pg scale; anchor-value uncertainty is excluded."
        ),
        "genome_size_reference_species": (
            "D. fuscus equal-image target estimate used as the conditional IOD denominator; "
            "the recovered standard images are not yet the denominator."
        ),
        "genome_size_reference_pg": (
            "Assembly-derived 1C-equivalent anchor computed as 16.1 Gbp / 0.978 Gbp per pg; "
            "conditional and not a measured C-value."
        ),
        "genome_size_calibration_pg_per_iod": (
            "Conditional assembly-derived reference pg divided by the D. fuscus equal-image "
            "IOD estimate."
        ),
        "genome_size_rank": (
            "Descending rank of the conditional fuscus-anchored pg sensitivity; identical in "
            "order to the underlying equal-image IOD rank."
        ),
        "estimate_support": "Whether one or multiple observed images support the species estimate.",
        "sample_size_support": "Recorded nucleus-count support category.",
        "quality_match_status": "Image-quality common-support or limited-overlap category.",
    },
    "S23": {
        "filename": "Source microscopy image filename.",
        "specimen_id": "Specimen identifier; one value denotes one animal.",
        "image_median_iod": "Median per-nucleus integrated optical density within the image.",
        "image_mean_iod": "Arithmetic mean per-nucleus integrated optical density within the image.",
        "image_median_nucleus_area_um2": "Median nucleus-mask area within the image.",
        "image_median_mean_od": "Median per-nucleus mean optical density within the image.",
        "image_median_match_log_edge_sharpness": (
            "Median log edge-sharpness matching feature within the image."
        ),
        "image_median_match_log_relative_ring_noise": (
            "Median log relative-ring-noise matching feature within the image."
        ),
    },
    "S24": {
        "template_species": "Quality-template species used for balance comparisons.",
        "max_abs_standardized_mean_difference": (
            "Maximum absolute standardized mean difference across the two matching features."
        ),
        "max_ks_distance": "Maximum Kolmogorov-Smirnov distance across the two matching features.",
        "passes_prespecified_balance_thresholds": (
            "True when max absolute SMD <= 0.10 and max KS distance <= 0.25."
        ),
    },
    "S25": {
        "quality_feature": "Image-quality feature tested as a residual technical diagnostic.",
        "within_species_centered_spearman_rho": (
            "Spearman correlation after centering the feature and log IOD within species."
        ),
        "p_value_descriptive_only": "Unadjusted descriptive p-value; not a confirmatory endpoint.",
    },
    "S26": {
        "comparison": (
            "Declared pairwise comparison among genome size, nucleus area, and cell area."
        ),
        "left_metric": (
            "Left-hand trait; genome-size comparisons use genome_size_pg_fuscus_anchored."
        ),
        "right_metric": "Right-hand trait in the declared pairwise comparison.",
    },
    "S27": {
        "genome_size_log10_z": (
            "Sample-standardized log10 conditional D. fuscus-anchored genome size; invariant "
            "to the positive IOD-to-pg scale factor."
        ),
        "genome_size_definition": COMMON_COLUMN_DEFINITIONS["genome_size_definition"],
    },
    "S28": {
        "relation_id": (
            "Declared directional regression identifier using genome, nucleus, and cell labels."
        ),
        "outcome": "Standardized outcome trait; genome_size denotes conditional genome size.",
        "predictor": (
            "Standardized predictor trait for this coefficient row; genome_size denotes "
            "conditional genome size."
        ),
        "path_trait_definition": COMMON_COLUMN_DEFINITIONS["path_trait_definition"],
    },
    "S40": {
        "scientific_name": "Display-form scientific name for the canonical species row.",
        "manual_review_inclusion": (
            "Whether the Path24 species was retained through a documented manual-review "
            "decision below the automatic object-count target."
        ),
        "n_superfamilies_included": (
            "Number of species-specific most-abundant positive classified TE "
            "superfamilies used by canonical S08; fixed at ten."
        ),
        "shannon_entropy": (
            "Natural-log Shannon entropy from the canonical species-specific top-10 "
            "superfamily S08 row after reclosure to one."
        ),
        "gini_simpson": (
            "Gini-Simpson index, 1 - sum(p_i^2), from the canonical species-specific "
            "top-10 superfamily S08 row after reclosure to one."
        ),
        "te_diversity_definition": (
            "Declared TE classification level, denominator, and diversity metric source "
            "for the joined row."
        ),
        "genome_size_pg_fuscus_anchored": (
            "Conditional target-to-D. fuscus equal-image IOD ratio multiplied by the selected "
            "16.462167689-pg/1C assembly anchor."
        ),
        "genome_size_pg_ci_low": (
            "Lower conditional 95% measurement-bootstrap interval for genome size; "
            "anchor-value uncertainty is excluded."
        ),
        "genome_size_pg_ci_high": (
            "Upper conditional 95% measurement-bootstrap interval for genome size; "
            "anchor-value uncertainty is excluded."
        ),
        "estimate_support": (
            "Whether one or multiple observed images support the genome-size estimate."
        ),
        "sample_size_support": (
            "Recorded nucleus-count support category for the genome-size estimate."
        ),
        "quality_match_status": "Image-quality common-support or limited-overlap category.",
        "genome_size_definition": COMMON_COLUMN_DEFINITIONS["genome_size_definition"],
    },
    "S41": {
        "data_origin": (
            "Row origin: exact Wang et al. 2025 Dryad Figure 2 value or current "
            "Desmognathus recomputation."
        ),
        "source_record_id": "Stable row identifier within the declared data origin.",
        "source_doi": "Dataset DOI for published rows or this_study for current rows.",
        "source_file": (
            "Exact published workbook member or semicolon-delimited canonical current "
            "input paths."
        ),
        "source_file_sha256": (
            "SHA-256 of the published workbook member or semicolon-delimited hashes of "
            "the canonical current inputs."
        ),
        "metric_scope": (
            "Method-specific top-ten-superfamily diversity contract; published and current "
            "scopes remain explicitly distinct."
        ),
        "scientific_name": "Scientific name retained from the source or canonical current row.",
        "common_name": "Published English name when supplied; blank for current rows.",
        "family": "Taxonomic family retained from the published source or current genus.",
        "comparison_group": "Display group used by the cross-study distribution panels.",
        "genome_size_gb": (
            "Published genome size in Gbp or current conditional D. fuscus-anchored "
            "IOD-derived sensitivity converted to Gbp at 0.978 Gbp per pg."
        ),
        "genome_size_ci_low_gb": (
            "Lower conditional interval for current genome-size sensitivity; blank for "
            "published and current group-only rows."
        ),
        "genome_size_ci_high_gb": (
            "Upper conditional interval for current genome-size sensitivity; blank for "
            "published and current group-only rows."
        ),
        "shannon_top10_superfamilies": (
            "Natural-log Shannon entropy after retaining the ten most abundant positive "
            "TE superfamilies and reclosing them to one."
        ),
        "gini_simpson_top10_superfamilies": (
            "Gini-Simpson index, 1 - sum(p_i^2), after retaining the ten most abundant "
            "positive TE superfamilies and reclosing them to one."
        ),
        "figure_inclusion": (
            "Declared group-only or scatter-and-group use; no row is used for a pooled fit."
        ),
    },
    "S42": {
        "display_species": "Italic-ready display label for the canonical TE34 species.",
        "reproductive_strategy": "Machine-readable expert reproductive-strategy annotation.",
        "reproductive_strategy_label": "Publication display label for reproductive strategy.",
        "strategy_assigned": "Whether an expert reproductive-strategy label is available.",
        "annotation_timing": (
            "Records that reproductive strategy was supplied after review of the initial "
            "PCA, making the annotation and associated tests post hoc."
        ),
    },
    "S43": {
        "strategy": "Machine-readable reproductive-strategy identifier.",
        "strategy_label": "Publication display label for reproductive strategy.",
        "group_n": "Number of assigned TE34 species in the strategy group.",
        "full_fit_smallest_containing_side_tips": (
            "Fewest tips on an edge-delimited side containing every group member in the "
            "full fitted Bonsai tree."
        ),
        "full_fit_edge_delimited": (
            "Whether one full-fit tree edge exactly delimits all and only the group members."
        ),
        "fixed_tree_smallest_side_p_value": (
            "One-sided fixed-tree label-permutation P value for a containing side this "
            "small or smaller."
        ),
        "fixed_tree_smallest_side_p_value_holm": (
            "Holm-adjusted P value across the nine reported localization tests."
        ),
        "feature_bootstrap_edge_delimited_count": (
            "Number of 200 bounded feature-resampling refits with an exact group edge."
        ),
        "feature_bootstrap_replicates": "Declared feature-resampling replicate count.",
        "feature_bootstrap_edge_delimited_fraction": (
            "Fraction of bounded feature-resampling refits with an exact group edge."
        ),
        "feature_bootstrap_smallest_side_median": (
            "Median smallest containing-side tip count across feature-resampling refits."
        ),
        "global_fitch_observed_steps": "Observed Fitch parsimony steps on the fixed tree.",
        "global_fitch_null_mean_steps": (
            "Mean Fitch parsimony steps across 99,999 fixed-tree label permutations."
        ),
        "global_fitch_p_value": "One-sided global Fitch label-permutation P value.",
        "global_fitch_p_value_holm": (
            "Global Fitch P value after Holm adjustment across nine reported tests."
        ),
    },
}


def column_definition(supplement_id: str, column: str) -> str:
    """Return the curated definition when one exists, otherwise an honest fallback."""

    scoped = SUPPLEMENT_COLUMN_DEFINITIONS.get(supplement_id, {})
    if column in scoped:
        return scoped[column]
    if column in COMMON_COLUMN_DEFINITIONS:
        return COMMON_COLUMN_DEFINITIONS[column]
    return "Source-named field retained for audit; see the canonical component documentation."
