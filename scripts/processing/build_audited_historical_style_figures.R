#!/usr/bin/env Rscript

# Rebuild Jake's historical R/ggplot figure suite on the frozen, audited data.
#
# This script deliberately does NOT embed or copy historical PNGs and does not
# execute dnaPipeTE, RepeatMasker, RepeatModeler, read mapping, or LTR discovery.
# The historical plotting grammar is recovered from:
#   /home/jake/Projects/misc_desmog_scripts/Master_data_vis.R
#   Git 1ad5e8c:te_pca_analysis.R
#   Git 1ad5e8c:scripts/visualization/te_landscape_plots.R
#   Git 1ad5e8c:scripts/visualization/ectopic_recomb_plots.R
#   VS Code history DWDT.R / Git 89c0879:Master_data_vis_PCA.R
#   Git 89c0879:Master_data_vis2_ectopicrecomb.R

options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = FALSE)
script_flag <- grep("^--file=", args, value = TRUE)
script_dir <- if (length(script_flag)) {
  dirname(normalizePath(sub("^--file=", "", script_flag[[1]]), mustWork = FALSE))
} else {
  "scripts/processing"
}
root <- normalizePath(file.path(script_dir, "..", ".."), mustWork = TRUE)

suppressPackageStartupMessages({
  library(ape)
  library(cluster)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(ggtree)
  library(jsonlite)
  library(nlme)
  library(patchwork)
  library(readr)
  library(scales)
  library(tidyr)
  library(viridis)
})

SEED <- 20260710L
CONFIGURED_QUANTIFICATION_SAMPLE_BP <- 1500000000

data_out <- file.path(root, "results", "data", "research_review")
figure_out <- file.path(root, "results", "figures", "research_review")
frozen_input_root <- file.path(data_out, "frozen_inputs")
dir.create(data_out, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_out, recursive = TRUE, showWarnings = FALSE)

WANG_JI_SOURCE_SHA256 <- "4ec3237a6ecd361c7e040b11dccfc9e2ac635cf10b739a4e95920015edccbb0e"
WANG_JI_DOI <- "10.1038/s42003-025-08127-3"
WANG_JI_DRYAD_DOI <- "10.5061/dryad.zpc866tkv"
wang_ji_local_source <- "/home/jake/Projects/Desmognathus_Results/wang_ji_shannondiversity.csv"
wang_ji_snapshot <- file.path(data_out, "wang_ji_2025_shannon_context_v1.csv")
wang_ji_provenance <- file.path(data_out, "wang_ji_2025_shannon_context_v1.provenance.json")

paths <- list(
  panel = file.path(frozen_input_root, "path_analysis", "data", "derived", "panels", "study_te_resource_panel34_v1.csv"),
  integrated_panel = file.path(frozen_input_root, "path_analysis", "data", "derived", "panels", "study_integrated_path_panel18_v1.csv"),
  genome_estimates = file.path(frozen_input_root, "path_analysis", "data", "external", "derived", "cellprofiler_final_species_results.csv"),
  genome_state = file.path(frozen_input_root, "path_analysis", "data", "external", "derived", "cellprofiler_genome_state_summary.csv"),
  retry_manifest = file.path(frozen_input_root, "results", "data", "corrected", "te34_replicate_averaged", "te34_replicate_averaged_v1.manifest.json"),
  mass = file.path(frozen_input_root, "results", "data", "corrected", "te34_replicate_averaged", "dnapipete_mass_accounting_te34_replicate_averaged_v1.csv"),
  quality = file.path(frozen_input_root, "results", "data", "legacy_descriptive", "dnapipete_species_quality_te34_legacy_descriptive_v1.csv"),
  order_composition = file.path(frozen_input_root, "results", "data", "legacy_descriptive", "te_order_composition_te34_legacy_descriptive_v1.csv"),
  superfamily_composition = file.path(frozen_input_root, "results", "data", "legacy_descriptive", "te_superfamily_composition_te34_legacy_descriptive_v1.csv"),
  diversity = file.path(frozen_input_root, "results", "data", "corrected", "te34_replicate_averaged", "te_diversity_mass_sensitivity_te34_replicate_averaged_v1.csv"),
  wang_ji_context = wang_ji_snapshot,
  wang_ji_provenance = wang_ji_provenance,
  ltr_elements = file.path(frozen_input_root, "results", "data", "corrected", "ectopic_ltr30", "ectopic_element_metrics_ltr30_v1.csv"),
  ltr_exclusions = file.path(frozen_input_root, "results", "data", "corrected", "ectopic_ltr30", "ectopic_excluded_elements_ltr30_v1.csv"),
  ltr_manifest = file.path(frozen_input_root, "results", "data", "corrected", "ectopic_ltr30", "ectopic_ltr30_v1.manifest.json"),
  tree = file.path(frozen_input_root, "results", "phylogeny", "processed_phylogeny.nwk")
)

missing_inputs <- setdiff(
  names(paths)[!file.exists(unlist(paths))],
  c("wang_ji_context", "wang_ji_provenance")
)
if (length(missing_inputs)) {
  stop("Missing required frozen input(s): ", paste(missing_inputs, collapse = ", "))
}

canonical_species <- function(value) {
  value <- tolower(trimws(as.character(value)))
  value <- sub("^desmognathus\\s+", "", value)
  value <- sub("^d\\.\\s*", "", value)
  value
}

display_species <- function(value) paste0("D. ", canonical_species(value))

assert_columns <- function(data, required, label) {
  absent <- setdiff(required, names(data))
  if (length(absent)) stop(label, " lacks required column(s): ", paste(absent, collapse = ", "))
}

assert_te34 <- function(data, label) {
  species <- sort(unique(canonical_species(data$species)))
  if (length(species) != 34L || !identical(species, panel_species)) {
    stop(label, " must contain the exact 34-species audited TE-resource panel")
  }
}

sha256 <- function(path) {
  result <- system2("sha256sum", path, stdout = TRUE, stderr = TRUE)
  if (!length(result) || !grepl("^[0-9a-f]{64}", result[[1]])) {
    stop("Could not calculate SHA256 for ", path)
  }
  sub("\\s+.*$", "", result[[1]])
}

if (!file.exists(wang_ji_snapshot) || !file.exists(wang_ji_provenance)) {
  if (!file.exists(wang_ji_local_source)) {
    stop(
      "Missing the portable Wang et al. snapshot and its hash-locked local source: ",
      wang_ji_local_source
    )
  }
  source_sha256 <- sha256(wang_ji_local_source)
  if (!identical(source_sha256, WANG_JI_SOURCE_SHA256)) {
    stop("Wang et al. historical composite no longer matches its audited SHA256")
  }
  wang_raw <- read_csv(wang_ji_local_source, show_col_types = FALSE, name_repair = "minimal")
  names(wang_raw) <- trimws(names(wang_raw))
  wang_raw <- wang_raw %>% mutate(across(where(is.character), trimws))
  assert_columns(wang_raw, c("English_name", "species", "Family", "GS", "SI", "GSI", "Type"), "Wang et al. context")
  wang_clean <- wang_raw %>%
    transmute(
      source_row = row_number(),
      english_name = English_name,
      species = species,
      family = Family,
      genome_size_gb = as.numeric(GS),
      shannon_top10_superfamilies = as.numeric(SI),
      gini_simpson_top10_superfamilies = as.numeric(GSI),
      group = Type,
      source = "Historical published-vertebrate composite used by prior project code"
    )
  if (nrow(wang_clean) != 87L || anyNA(wang_clean)) {
    stop("Cleaned Wang et al. historical composite must contain 87 complete rows")
  }
  write_csv(wang_clean, wang_ji_snapshot)
  write_json(
    list(
      artifact_id = "wang_ji_2025_shannon_context_v1",
      local_source_path = wang_ji_local_source,
      local_source_sha256 = source_sha256,
      local_source_rows = 87L,
      cleaning = "Whitespace-only header/value cleanup and explicit typed column names; published values unchanged.",
      context_reference = paste0("Wang et al. 2025, Communications Biology, DOI ", WANG_JI_DOI),
      published_data_repository = paste0("Dryad DOI ", WANG_JI_DRYAD_DOI),
      published_metric_contract = "Natural-log Shannon entropy over the top 10 most abundant TE superfamilies; TE-superfamily sequence proportions are abundances.",
      source_scope_disclosure = "This 87-row local historical composite is not an exact export of the paper's 84-species comparative analysis; it includes 4 caecilians and 10 salamanders and is retained only as descriptive cross-study context.",
      article_url = paste0("https://doi.org/", WANG_JI_DOI)
    ),
    wang_ji_provenance,
    pretty = TRUE,
    auto_unbox = TRUE
  )
}

wang_context <- read_csv(wang_ji_snapshot, show_col_types = FALSE)
wang_context_provenance <- fromJSON(wang_ji_provenance, simplifyVector = FALSE)
assert_columns(
  wang_context,
  c(
    "source_row", "english_name", "species", "family", "genome_size_gb",
    "shannon_top10_superfamilies", "gini_simpson_top10_superfamilies", "group", "source"
  ),
  "Portable Wang et al. context"
)
if (nrow(wang_context) != 87L || anyNA(wang_context)) {
  stop("Portable Wang et al. context must contain 87 complete rows")
}
if (!identical(wang_context_provenance$local_source_sha256, WANG_JI_SOURCE_SHA256) ||
    !identical(wang_context_provenance$local_source_rows, 87L) ||
    !identical(wang_context_provenance$published_data_repository, paste0("Dryad DOI ", WANG_JI_DRYAD_DOI))) {
  stop("Portable historical vertebrate-context provenance no longer matches its audited contract")
}

relative_path <- function(path) {
  normalized <- normalizePath(path, mustWork = FALSE)
  prefix <- paste0(normalizePath(root, mustWork = TRUE), .Platform$file.sep)
  sub(paste0("^", prefix), "", normalized)
}

theme_publication <- theme_minimal() +
  theme(
    text = element_text(family = "serif", size = 12),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.position = "right"
  )

theme_ltr_publication <- theme_bw() +
  theme(
    text = element_text(family = "serif", size = 12),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
    panel.grid.major = element_line(colour = "grey92"),
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5),
    strip.background = element_blank(),
    strip.text = element_text(size = 14)
  )

figure_specs <- list()
save_figure <- function(plot, stem, width, height) {
  png_path <- file.path(figure_out, paste0(stem, ".png"))
  pdf_path <- file.path(figure_out, paste0(stem, ".pdf"))
  ggsave(png_path, plot = plot, width = width, height = height, units = "in", dpi = 300, bg = "white")
  ggsave(pdf_path, plot = plot, width = width, height = height, units = "in", bg = "white")
  figure_specs[[stem]] <<- list(
    stem = stem,
    png = relative_path(png_path),
    pdf = relative_path(pdf_path),
    width_inches = width,
    height_inches = height,
    png_sha256 = sha256(png_path),
    pdf_sha256 = sha256(pdf_path)
  )
}

# -----------------------------------------------------------------------------
# Panel and retry-policy gates
# -----------------------------------------------------------------------------

panel <- read_csv(paths$panel, show_col_types = FALSE) %>%
  mutate(species = canonical_species(species)) %>%
  arrange(species)
assert_columns(panel, c("species", "te_sra_accession", "te_assembly_accession"), "TE34 panel")
if (nrow(panel) != 34L || anyDuplicated(panel$species)) stop("TE34 panel is not 34 unique species")
panel_species <- sort(panel$species)

if (!identical(panel$te_assembly_accession[panel$species == "fuscus"], "GCA_032353935.1")) {
  stop("Audited fuscus resource must be GCA_032353935.1")
}

retry_manifest <- fromJSON(paths$retry_manifest, simplifyVector = FALSE)
if (!identical(retry_manifest$dnapipete_run_aggregation, "equal_weight_mean_of_run_level_estimates")) {
  stop("TE34 retry manifest no longer records equal-weight run averaging")
}
orestes_labels <- sort(unlist(retry_manifest$orestes_source_labels))
if (!identical(orestes_labels, c("SRX19952890", "SRX19952890R2"))) {
  stop("Unexpected orestes retry labels in audited TE34 manifest")
}

# -----------------------------------------------------------------------------
# Assembly quality: literal historical bar grammar, current audited accessions
# -----------------------------------------------------------------------------

assembly_csv <- file.path(data_out, "assembly_quality_te34_v1.csv")
required_assembly_columns <- c(
  "species", "display_species", "te_assembly_accession", "ncbi_reported_organism",
  "species_mapping_basis", "total_sequence_length_bp", "contig_n50_bp",
  "scaffold_n50_bp", "number_of_scaffolds", "number_of_contigs",
  "ncbi_dataset_report_url", "metadata_retrieved_utc"
)

assembly_cache_valid <- FALSE
if (file.exists(assembly_csv)) {
  cached <- suppressMessages(read_csv(assembly_csv, show_col_types = FALSE))
  if (all(required_assembly_columns %in% names(cached))) {
    cached$species <- canonical_species(cached$species)
    assembly_cache_valid <- nrow(cached) == 34L &&
      identical(sort(cached$species), panel_species) &&
      identical(sort(cached$te_assembly_accession), sort(panel$te_assembly_accession))
  }
}

if (assembly_cache_valid) {
  assembly <- cached %>% arrange(species)
} else {
  message("Resolving current NCBI assembly reports for 34 audited accessions...")
  rows <- lapply(seq_len(nrow(panel)), function(index) {
    record <- panel[index, ]
    accession <- record$te_assembly_accession[[1]]
    endpoint <- paste0(
      "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/",
      accession,
      "/dataset_report"
    )
    payload <- tryCatch(
      fromJSON(endpoint, simplifyVector = FALSE),
      error = function(error) stop("NCBI assembly report failed for ", accession, ": ", conditionMessage(error))
    )
    if (length(payload$reports) != 1L) stop("Expected one NCBI report for ", accession)
    report <- payload$reports[[1]]
    stats <- report$assembly_stats
    reported_name <- report$organism$organism_name
    mapping <- if (record$species[[1]] == "fuscus") {
      if (!identical(reported_name, "Desmognathus planiceps")) {
        stop("Expected the public GCA_032353935.1 report to retain the planiceps label")
      }
      "expert_reidentified_resource: NCBI planiceps row analyzed as fuscus"
    } else {
      "audited_TE34_accession_crosswalk"
    }
    tibble(
      species = record$species[[1]],
      display_species = display_species(record$species[[1]]),
      te_assembly_accession = accession,
      ncbi_reported_organism = reported_name,
      species_mapping_basis = mapping,
      total_sequence_length_bp = as.numeric(stats$total_sequence_length),
      contig_n50_bp = as.numeric(stats$contig_n50),
      scaffold_n50_bp = as.numeric(stats$scaffold_n50),
      number_of_scaffolds = as.numeric(stats$number_of_scaffolds),
      number_of_contigs = as.numeric(stats$number_of_contigs),
      ncbi_dataset_report_url = endpoint,
      metadata_retrieved_utc = format(Sys.time(), tz = "UTC", usetz = TRUE)
    )
  })
  assembly <- bind_rows(rows) %>% arrange(species)
  if (any(!is.finite(as.matrix(assembly[, c(
    "total_sequence_length_bp", "contig_n50_bp", "scaffold_n50_bp",
    "number_of_scaffolds", "number_of_contigs"
  )])))) stop("NCBI returned incomplete assembly metrics")
  write_csv(assembly, assembly_csv)
}

assert_te34(assembly, "Assembly-quality table")

# These four statements deliberately retain Master_data_vis.R's color and
# geom_bar grammar. Only label size/number formatting and 2x2 assembly change.
p_assembly_contig <- ggplot(assembly, aes(x = display_species, y = contig_n50_bp)) +
  geom_bar(stat = "identity", fill = "steelblue", color = "black") +
  labs(title = "Assembly Contig N50", x = "Species", y = "Contig N50 (bp)") +
  scale_y_continuous(labels = label_number(big.mark = ","), expand = expansion(mult = c(0, 0.05))) +
  theme_publication +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7), legend.position = "none")

p_assembly_scaffold <- ggplot(assembly, aes(x = display_species, y = scaffold_n50_bp)) +
  geom_bar(stat = "identity", fill = "steelblue1", color = "black") +
  labs(title = "Assembly Scaffold N50", x = "Species", y = "Scaffold N50 (bp)") +
  scale_y_continuous(labels = label_number(big.mark = ","), expand = expansion(mult = c(0, 0.05))) +
  theme_publication +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7), legend.position = "none")

p_assembly_scaffolds <- ggplot(assembly, aes(x = display_species, y = number_of_scaffolds)) +
  geom_bar(stat = "identity", fill = "darkgreen", color = "black") +
  labs(title = "Assembly Number of Scaffolds", x = "Species", y = "Number of scaffolds") +
  scale_y_continuous(labels = label_number(scale_cut = cut_short_scale()), expand = expansion(mult = c(0, 0.05))) +
  theme_publication +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7), legend.position = "none")

p_assembly_length <- ggplot(assembly, aes(x = display_species, y = total_sequence_length_bp)) +
  geom_bar(stat = "identity", fill = "tomato", color = "black") +
  labs(title = "Assembly Sequence Length", x = "Species", y = "Sequence length (bp)") +
  scale_y_continuous(labels = label_number(scale_cut = cut_short_scale()), expand = expansion(mult = c(0, 0.05))) +
  theme_publication +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7), legend.position = "none")

assembly_plot <- (p_assembly_contig + p_assembly_scaffold) /
  (p_assembly_scaffolds + p_assembly_length) +
  plot_annotation(
    title = "Assembly quality of the audited TE-resource panel",
    subtitle = "TE34; GCA_032353935.1 is the public planiceps-labeled resource expert-reidentified as fuscus",
    theme = theme(text = element_text(family = "serif"), plot.title = element_text(face = "bold", hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))
  )
save_figure(assembly_plot, "assembly_quality_te34_v1", 18, 12)

# -----------------------------------------------------------------------------
# dnaPipeTE input/QC and mass accounting
# -----------------------------------------------------------------------------

quality <- read_csv(paths$quality, show_col_types = FALSE) %>%
  mutate(species = canonical_species(species))
mass <- read_csv(paths$mass, show_col_types = FALSE) %>%
  mutate(species = canonical_species(species))
assert_te34(quality, "dnaPipeTE quality table")
assert_te34(mass, "dnaPipeTE mass-accounting table")
assert_columns(quality, c(
  "species", "n_dnapipete_runs", "run_labels", "n_components", "total_reads",
  "total_aligned_bases", "median_hitlength_contig_ratio", "fraction_components_classified"
), "dnaPipeTE quality table")
assert_columns(mass, c(
  "species", "total_aligned_bases", "order_retained_aligned_bases", "order_unresolved_aligned_bases",
  "order_retained_fraction", "order_unresolved_fraction", "superfamily_retained_fraction",
  "superfamily_unresolved_fraction", "n_dnapipete_runs", "dnapipete_source_labels", "run_aggregation"
), "dnaPipeTE mass table")

quality_export <- panel %>%
  select(species, te_sra_accession, te_assembly_accession) %>%
  left_join(quality, by = "species") %>%
  left_join(
    mass %>% select(
      species,
      mass_total_aligned_bases = total_aligned_bases,
      order_retained_fraction,
      order_unresolved_fraction,
      superfamily_retained_fraction,
      superfamily_unresolved_fraction,
      mass_n_dnapipete_runs = n_dnapipete_runs,
      dnapipete_source_labels,
      run_aggregation
    ),
    by = "species"
  ) %>%
  mutate(
    display_species = display_species(species),
    configured_quantification_sample_bp = CONFIGURED_QUANTIFICATION_SAMPLE_BP,
    repeat_aligned_fraction_of_sample = mass_total_aligned_bases / CONFIGURED_QUANTIFICATION_SAMPLE_BP,
    retry_policy = if_else(n_dnapipete_runs > 1L, "equal_weight_run_mean", "single_active_run")
  ) %>%
  arrange(species)

if (anyNA(quality_export)) stop("dnaPipeTE quality join introduced missing values")
if (max(abs(quality_export$total_aligned_bases - quality_export$mass_total_aligned_bases)) > 1) {
  stop("dnaPipeTE quality and corrected mass ledgers disagree on species-level aligned bases")
}
if (!identical(quality_export$n_dnapipete_runs, quality_export$mass_n_dnapipete_runs)) {
  stop("dnaPipeTE quality and mass ledgers disagree on retry counts")
}
write_csv(quality_export, file.path(data_out, "dnapipete_quality_te34_v1.csv"))

quality_long <- quality_export %>%
  select(display_species, total_reads, n_components, median_hitlength_contig_ratio, fraction_components_classified) %>%
  pivot_longer(-display_species, names_to = "metric", values_to = "value") %>%
  mutate(
    metric = recode(
      metric,
      total_reads = "Reads assigned to components",
      n_components = "dnaPipeTE components",
      median_hitlength_contig_ratio = "Median RepeatMasker-hit / contig length",
      fraction_components_classified = "Fraction of components classified"
    )
  )

dnapipete_quality_plot <- ggplot(quality_long, aes(x = display_species, y = value)) +
  geom_bar(stat = "identity", fill = "steelblue1", color = "black", linewidth = 0.2) +
  facet_wrap(~metric, scales = "free_y", ncol = 2) +
  scale_y_continuous(labels = label_number(scale_cut = cut_short_scale()), expand = expansion(mult = c(0, 0.06))) +
  labs(
    title = "dnaPipeTE input and component quality",
    subtitle = "TE34; same-species retries are equal-weight run means and never additional species",
    x = "Species", y = NULL
  ) +
  theme_publication +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
    legend.position = "none",
    strip.text = element_text(family = "serif", face = "bold")
  )
save_figure(dnapipete_quality_plot, "dnapipete_quality_te34_v1", 16, 10)

species_levels <- quality_export %>%
  arrange(desc(repeat_aligned_fraction_of_sample)) %>%
  pull(display_species)

sample_mass <- quality_export %>%
  transmute(
    display_species = factor(display_species, levels = species_levels),
    denominator = "Configured 1.5-Gb quantification sample",
    Classified = repeat_aligned_fraction_of_sample,
    Unresolved = 1 - repeat_aligned_fraction_of_sample
  ) %>%
  pivot_longer(c(Classified, Unresolved), names_to = "mass_class", values_to = "proportion") %>%
  mutate(mass_class = recode(mass_class, Classified = "Repeat-aligned", Unresolved = "Not aligned to repeat contigs"))

order_mass <- quality_export %>%
  transmute(
    display_species = factor(display_species, levels = species_levels),
    denominator = "Order classification of repeat-aligned mass",
    Classified = order_retained_fraction,
    Unresolved = order_unresolved_fraction
  ) %>%
  pivot_longer(c(Classified, Unresolved), names_to = "mass_class", values_to = "proportion") %>%
  mutate(mass_class = recode(mass_class, Classified = "Classified", Unresolved = "Unresolved"))

superfamily_mass <- quality_export %>%
  transmute(
    display_species = factor(display_species, levels = species_levels),
    denominator = "Superfamily classification of repeat-aligned mass",
    Classified = superfamily_retained_fraction,
    Unresolved = superfamily_unresolved_fraction
  ) %>%
  pivot_longer(c(Classified, Unresolved), names_to = "mass_class", values_to = "proportion") %>%
  mutate(mass_class = recode(mass_class, Classified = "Classified", Unresolved = "Unresolved"))

mass_plot_data <- bind_rows(sample_mass, order_mass, superfamily_mass) %>%
  mutate(
    denominator = factor(denominator, levels = c(
      "Configured 1.5-Gb quantification sample",
      "Order classification of repeat-aligned mass",
      "Superfamily classification of repeat-aligned mass"
    )),
    mass_class = factor(mass_class, levels = c(
      "Repeat-aligned", "Classified", "Unresolved", "Not aligned to repeat contigs"
    ))
  )

te_mass_plot <- ggplot(mass_plot_data, aes(x = display_species, y = proportion, fill = mass_class)) +
  geom_bar(stat = "identity", position = "stack", color = "black", linewidth = 0.15) +
  facet_wrap(~denominator, ncol = 1) +
  scale_fill_manual(values = c(
    "Repeat-aligned" = "tomato", "Classified" = "steelblue1",
    "Unresolved" = "#E69F00", "Not aligned to repeat contigs" = "grey85"
  ), drop = FALSE) +
  scale_y_continuous(labels = label_percent(accuracy = 1), limits = c(0, 1), expand = c(0, 0)) +
  labs(
    title = "dnaPipeTE repeat mass and classification accounting",
    subtitle = "Stacking is corrected; facet titles state the denominator",
    x = "Species", y = "Percentage (%)", fill = "Mass category"
  ) +
  theme_publication +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
    strip.text = element_text(family = "serif", face = "bold")
  )
save_figure(te_mass_plot, "te_mass_composition_te34_v1", 16, 11)

# -----------------------------------------------------------------------------
# Historical mean order and superfamily plots on retry-averaged TE34 composition
# -----------------------------------------------------------------------------

order_comp <- read_csv(paths$order_composition, show_col_types = FALSE) %>%
  mutate(species = canonical_species(species))
superfamily_comp <- read_csv(paths$superfamily_composition, show_col_types = FALSE) %>%
  mutate(species = canonical_species(species))
assert_te34(order_comp, "Order-composition table")
assert_te34(superfamily_comp, "Superfamily-composition table")

order_features <- setdiff(names(order_comp), "species")
superfamily_features <- setdiff(names(superfamily_comp), "species")
if (max(abs(rowSums(order_comp[, order_features]) - 1)) > 1e-8) stop("Order compositions do not close to one")
if (max(abs(rowSums(superfamily_comp[, superfamily_features]) - 1)) > 1e-8) stop("Superfamily compositions do not close to one")

retro_orders <- c("DIRS", "LINE", "LTR", "SINE")
dna_orders <- c("Helitron", "Maverick", "PLE", "TIR", "YR", "Transposon Derivatives")
retro_superfamilies <- c("CR1", "Copia", "DIRS", "Gypsy", "Jockey", "L1", "L2", "Penelope", "tRNA")
dna_superfamilies <- setdiff(superfamily_features, retro_superfamilies)

# Literal structure from 1ad5e8c:scripts/visualization/te_landscape_plots.R.
order_means <- order_comp %>%
  select(-species) %>%
  mutate(RowTotal = rowSums(.)) %>%
  mutate(across(-RowTotal, ~ . / RowTotal * 100)) %>%
  select(-RowTotal) %>%
  summarise(across(everything(), mean)) %>%
  pivot_longer(cols = everything(), names_to = "Category", values_to = "Percentage") %>%
  mutate(Type = case_when(
    Category %in% retro_orders ~ "Retrotransposon",
    Category %in% dna_orders ~ "DNA transposon",
    TRUE ~ "Other"
  )) %>%
  filter(Percentage >= 0.05) %>%
  arrange(desc(Percentage))

order_plot <- ggplot(order_means, aes(x = reorder(Category, -Percentage), y = Percentage, fill = Type)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c(
    "Retrotransposon" = "tomato",
    "DNA transposon" = "steelblue1",
    "Other" = "grey60"
  )) +
  scale_y_continuous(
    breaks = seq(0, ceiling(max(order_means$Percentage) / 10) * 10, by = 10),
    expand = c(0, 0)
  ) +
  labs(
    title = "Mean Transposon Order Distribution Across Species",
    subtitle = "Audited TE34; equal species weight after within-species retry averaging",
    x = "Order", y = "Percentage (%)", fill = NULL
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    text = element_text(family = "serif"),
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5)
  )
save_figure(order_plot, "te_mean_order_te34_v1", 12, 7)

superfamily_means <- superfamily_comp %>%
  select(-species) %>%
  mutate(RowTotal = rowSums(.)) %>%
  mutate(across(-RowTotal, ~ . / RowTotal * 100)) %>%
  select(-RowTotal) %>%
  summarise(across(everything(), mean)) %>%
  pivot_longer(cols = everything(), names_to = "Category", values_to = "Percentage") %>%
  mutate(Type = case_when(
    Category %in% retro_superfamilies ~ "Retrotransposon",
    Category %in% dna_superfamilies ~ "DNA transposon",
    TRUE ~ "Other"
  )) %>%
  filter(Percentage >= 0.5) %>%
  arrange(desc(Percentage))

superfamily_plot <- ggplot(superfamily_means, aes(x = reorder(Category, -Percentage), y = Percentage, fill = Type)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c(
    "Retrotransposon" = "tomato",
    "DNA transposon" = "steelblue1",
    "Other" = "grey60"
  )) +
  scale_y_continuous(
    breaks = seq(0, ceiling(max(superfamily_means$Percentage) / 5) * 5, by = 5),
    minor_breaks = seq(0, ceiling(max(superfamily_means$Percentage)), by = 1),
    expand = c(0, 0)
  ) +
  labs(
    title = "Mean Transposon Superfamily Distribution Across Species",
    subtitle = "Audited TE34; equal species weight after within-species retry averaging",
    x = "Superfamily", y = "Percentage (%)", fill = NULL
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    text = element_text(family = "serif"),
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    panel.grid.minor = element_line(color = "gray90")
  )
save_figure(superfamily_plot, "te_mean_superfamily_te34_v1", 12, 7)

# -----------------------------------------------------------------------------
# Shannon diversity on corrected, retry-averaged TE34 sensitivity table
# -----------------------------------------------------------------------------

diversity <- read_csv(paths$diversity, show_col_types = FALSE) %>%
  mutate(
    species = canonical_species(species),
    display_species = display_species(species),
    Group = case_when(
      te_level == "order" & composition_mode == "classified_conditional" ~ "Order: classified mass",
      te_level == "order" & composition_mode == "mass_aware_unresolved_bin" ~ "Order: mass-aware",
      te_level == "superfamily" & composition_mode == "classified_conditional" ~ "Superfamily: classified mass",
      te_level == "superfamily" & composition_mode == "mass_aware_unresolved_bin" ~ "Superfamily: mass-aware",
      TRUE ~ paste(te_level, composition_mode)
    )
  )
assert_te34(diversity, "Diversity sensitivity table")
if (nrow(diversity) != 136L) stop("Expected 34 species x 2 levels x 2 denominator modes")

# Current denominator/level sensitivity companion; this is intentionally kept
# separate from the restored historical cross-taxon boxplot below.
diversity_plot <- ggplot(diversity, aes(x = Group, y = shannon_entropy, fill = Group)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.12, alpha = 0.65, size = 1.7) +
  labs(
    title = "TE Shannon Diversity Index by Level and Denominator",
    subtitle = "Audited retry-averaged TE34 sensitivity table",
    x = "TE level and mass denominator", y = "Shannon diversity index"
  ) +
  theme_publication +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )
save_figure(diversity_plot, "te_diversity_te34_v1", 12, 7)

# Restore Master_shannon_diversity_vis.R's literal cross-taxon boxplot grammar.
# Current Desmognathus values follow Wang et al. 2025's published metric
# contract: natural-log Shannon entropy over the top 10 most abundant TE
# superfamilies, after renormalizing those ten abundances to sum to one.
desmo_top10_long <- bind_rows(lapply(seq_len(nrow(superfamily_comp)), function(index) {
  species_name <- superfamily_comp$species[[index]]
  values <- as.numeric(superfamily_comp[index, superfamily_features])
  names(values) <- superfamily_features
  ordered <- sort(values[is.finite(values) & values > 0], decreasing = TRUE)
  if (length(ordered) < 10L) stop("Fewer than 10 positive superfamilies for ", species_name)
  selected <- ordered[seq_len(10L)]
  probabilities <- selected / sum(selected)
  tibble(
    species = species_name,
    display_species = display_species(species_name),
    rank = seq_len(10L),
    superfamily = names(selected),
    classified_conditional_proportion = as.numeric(selected),
    top10_renormalized_proportion = as.numeric(probabilities),
    shannon_component = -as.numeric(probabilities * log(probabilities)),
    gini_simpson_squared_component = as.numeric(probabilities^2)
  )
}))

desmo_top10 <- desmo_top10_long %>%
  group_by(species, display_species) %>%
  summarise(
    top10_total_classified_fraction = sum(classified_conditional_proportion),
    shannon_top10_superfamilies = sum(shannon_component),
    gini_simpson_top10_superfamilies = 1 - sum(gini_simpson_squared_component),
    .groups = "drop"
  )
assert_te34(desmo_top10, "Top-10-superfamily Shannon table")
if (nrow(desmo_top10_long) != 340L ||
    any(abs(desmo_top10$shannon_top10_superfamilies) > log(10) + 1e-10) ||
    any(abs(desmo_top10_long %>% group_by(species) %>% summarise(total = sum(top10_renormalized_proportion)) %>% pull(total) - 1) > 1e-10)) {
  stop("Top-10-superfamily diversity contract failed")
}

desmo_context <- desmo_top10 %>%
  transmute(
    source_row = NA_integer_,
    english_name = NA_character_,
    species = paste("Desmognathus", species),
    family = "Plethodontidae",
    genome_size_gb = NA_real_,
    shannon_top10_superfamilies,
    gini_simpson_top10_superfamilies,
    group = "Desmognathus",
    source = "Current audited retry-averaged TE34 superfamily composition"
  )

historical_desmognathus_rows <- grepl("^Desmognathus\\s", wang_context$species, ignore.case = TRUE)
cross_taxon <- bind_rows(wang_context[!historical_desmognathus_rows, ], desmo_context)
cross_taxon_counts <- cross_taxon %>% count(group, source, name = "n_species")
cross_taxon_path <- file.path(data_out, "te_shannon_cross_taxon_context_te34_v1.csv")
cross_taxon_counts_path <- file.path(data_out, "te_shannon_cross_taxon_group_counts_v1.csv")
desmo_top10_long_path <- file.path(data_out, "te_top10_superfamily_contributions_te34_v1.csv")
write_csv(cross_taxon, cross_taxon_path)
write_csv(cross_taxon_counts, cross_taxon_counts_path)
write_csv(desmo_top10_long, desmo_top10_long_path)

cross_taxon_plot <- ggplot(
  cross_taxon,
  aes(x = group, y = shannon_top10_superfamilies, fill = group)
) +
  geom_boxplot() +
  labs(
    title = "Shannon Diversity Index by Group",
    subtitle = "Descriptive cross-study context; library, annotation protocol, and group sample sizes differ",
    x = "Group",
    y = "Shannon Index (SI)"
  ) +
  theme_publication +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )
save_figure(cross_taxon_plot, "te_shannon_cross_taxon_context_te34_v1", 10, 6)

# Restore Master_shannon_diversity_vis.R's historical genome-size/Shannon
# scatter with the current 18-species overlap.  The Desmognathus x values are
# the verified, fuscus-anchored image-IOD estimates carried by the frozen
# CellProfiler bridge; they are estimates with bootstrap intervals, not direct
# flow-cytometric or Feulgen C-value measurements.
genome_estimates <- read_csv(paths$genome_estimates, show_col_types = FALSE) %>%
  mutate(
    species = canonical_species(species),
    across(
      c(
        primary_genome_pg, primary_genome_gb, primary_genome_ci_low_pg,
        primary_genome_ci_high_pg, primary_n_images, primary_n_specimens,
        primary_n_nuclei
      ),
      as.numeric
    )
  )
assert_columns(
  genome_estimates,
  c(
    "species", "primary_genome_pg", "primary_genome_gb",
    "primary_genome_ci_low_pg", "primary_genome_ci_high_pg",
    "primary_n_images", "primary_n_specimens", "primary_n_nuclei",
    "primary_measurement_kind", "primary_support_tier",
    "primary_support_warnings", "result_status", "flag_summary"
  ),
  "Frozen fuscus-anchored genome-size estimates"
)
if (nrow(genome_estimates) != 21L || anyDuplicated(genome_estimates$species)) {
  stop("Frozen genome-size estimate table must contain 21 unique species")
}
if (anyNA(genome_estimates[c(
  "primary_genome_pg", "primary_genome_gb", "primary_genome_ci_low_pg",
  "primary_genome_ci_high_pg"
)]) ||
    any(genome_estimates$primary_genome_ci_low_pg >= genome_estimates$primary_genome_pg) ||
    any(genome_estimates$primary_genome_ci_high_pg <= genome_estimates$primary_genome_pg) ||
    n_distinct(round(genome_estimates$primary_genome_pg, 6)) < 15L ||
    any(abs(genome_estimates$primary_genome_gb - 0.978 * genome_estimates$primary_genome_pg) > 1e-6)) {
  stop("Genome-size estimates failed variability, interval, or pg-to-Gb validation")
}

genome_state <- read_csv(paths$genome_state, show_col_types = FALSE) %>%
  mutate(
    species = canonical_species(species),
    reference_species = display_species(reference_species),
    reference_genome_pg = as.numeric(reference_genome_pg)
  )
assert_columns(
  genome_state,
  c("species", "reference_species", "reference_genome_pg", "genome_pg_estimate"),
  "Frozen genome-state calibration table"
)
calibration_anchors <- genome_state %>%
  distinct(reference_species, reference_genome_pg)
if (nrow(calibration_anchors) != 1L ||
    calibration_anchors$reference_species[[1]] != "D. fuscus" ||
    abs(calibration_anchors$reference_genome_pg[[1]] - 16.36) > 1e-10) {
  stop("Genome-size estimate table must retain the declared D. fuscus = 16.36 pg anchor")
}
fuscus_estimate <- genome_estimates %>% filter(species == "fuscus") %>% pull(primary_genome_pg)
if (length(fuscus_estimate) != 1L || abs(fuscus_estimate - 16.36) > 1e-10) {
  stop("The fuscus anchor row is absent or no longer equals 16.36 pg")
}

integrated_panel <- read_csv(paths$integrated_panel, show_col_types = FALSE) %>%
  mutate(species = canonical_species(species))
if (nrow(integrated_panel) != 18L || anyDuplicated(integrated_panel$species)) {
  stop("Integrated path panel must contain 18 unique species")
}

genome_shannon_match_audit <- full_join(
  desmo_top10 %>% select(species, shannon_top10_superfamilies),
  genome_estimates %>%
    select(species, primary_genome_pg, primary_support_tier, result_status),
  by = "species"
) %>%
  mutate(
    has_te_shannon = is.finite(shannon_top10_superfamilies),
    has_genome_size_estimate = is.finite(primary_genome_pg),
    included_in_scatter = has_te_shannon & has_genome_size_estimate,
    exclusion_reason = case_when(
      included_in_scatter ~ "included_exact_te_genome_overlap",
      !has_te_shannon ~ "no_audited_te34_shannon_value",
      !has_genome_size_estimate ~ "no_current_fuscus_anchored_genome_size_estimate",
      TRUE ~ "excluded_unspecified"
    )
  ) %>%
  rename(canonical_species = species) %>%
  arrange(canonical_species)

included_species <- genome_shannon_match_audit %>%
  filter(included_in_scatter) %>%
  pull(canonical_species)
if (nrow(genome_shannon_match_audit) != 37L ||
    length(included_species) != 18L ||
    !setequal(included_species, integrated_panel$species)) {
  stop("Genome-size/Shannon overlap must equal the exact integrated path18 species set")
}

genome_shannon_match_audit_path <- file.path(
  data_out, "te_shannon_genome_size_match_audit_te34_v1.csv"
)
write_csv(genome_shannon_match_audit, genome_shannon_match_audit_path)

desmo_genome_scatter <- desmo_top10 %>%
  inner_join(genome_estimates, by = "species") %>%
  transmute(
    canonical_species = species,
    species = paste("Desmognathus", species),
    group = "Desmognathus",
    genome_size_gb = primary_genome_gb,
    genome_size_ci_low_gb = 0.978 * primary_genome_ci_low_pg,
    genome_size_ci_high_gb = 0.978 * primary_genome_ci_high_pg,
    genome_size_pg = primary_genome_pg,
    genome_size_ci_low_pg = primary_genome_ci_low_pg,
    genome_size_ci_high_pg = primary_genome_ci_high_pg,
    shannon_top10_superfamilies,
    genome_measurement_kind = "fuscus_anchored_image_iod_genome_size_estimate",
    is_direct_c_value = FALSE,
    genome_support_tier = primary_support_tier,
    genome_support_warnings = coalesce(primary_support_warnings, ""),
    genome_result_status = result_status,
    n_genome_images = primary_n_images,
    n_genome_specimens = primary_n_specimens,
    n_genome_nuclei = primary_n_nuclei,
    genome_source = "Frozen verified CellProfiler linked-nucleus IOD bridge; D. fuscus = 16.36 pg anchor",
    shannon_source = "Current audited retry-averaged TE34 top-10-superfamily composition"
  )

historical_salamander_scatter <- wang_context %>%
  filter(group == "Salamanders") %>%
  transmute(
    canonical_species = tolower(species),
    species,
    group,
    genome_size_gb,
    genome_size_ci_low_gb = NA_real_,
    genome_size_ci_high_gb = NA_real_,
    genome_size_pg = NA_real_,
    genome_size_ci_low_pg = NA_real_,
    genome_size_ci_high_pg = NA_real_,
    shannon_top10_superfamilies,
    genome_measurement_kind = "historical_published_cross_taxon_context",
    is_direct_c_value = NA,
    genome_support_tier = "published_context",
    genome_support_warnings = "Cross-study context; source assays and annotation workflows differ",
    genome_result_status = "historical_context_only",
    n_genome_images = NA_real_,
    n_genome_specimens = NA_real_,
    n_genome_nuclei = NA_real_,
    genome_source = source,
    shannon_source = source
  )

genome_shannon_scatter <- bind_rows(historical_salamander_scatter, desmo_genome_scatter)
if (nrow(genome_shannon_scatter) != 28L ||
    sum(genome_shannon_scatter$group == "Desmognathus") != 18L ||
    sum(genome_shannon_scatter$group == "Salamanders") != 10L ||
    anyNA(genome_shannon_scatter[c("genome_size_gb", "shannon_top10_superfamilies")])) {
  stop("Genome-size/Shannon scatter must contain 18 current Desmognathus and 10 historical salamanders")
}

genome_shannon_scatter_path <- file.path(
  data_out, "te_shannon_genome_size_scatter_te18_v1.csv"
)
write_csv(genome_shannon_scatter, genome_shannon_scatter_path)

# The point/color/lab/theme grammar is the literal historical plot. Horizontal
# intervals are the only visual addition and expose uncertainty in the current
# Desmognathus x estimates.
genome_shannon_scatter_plot <- ggplot(
  genome_shannon_scatter,
  aes(x = genome_size_gb, y = shannon_top10_superfamilies, color = group)
) +
  geom_segment(
    data = desmo_genome_scatter,
    aes(
      x = genome_size_ci_low_gb, xend = genome_size_ci_high_gb,
      y = shannon_top10_superfamilies, yend = shannon_top10_superfamilies,
      color = group
    ),
    inherit.aes = FALSE,
    linewidth = 0.55,
    alpha = 0.45
  ) +
  geom_point(size = 3) +
  labs(
    title = "Shannon Diversity Index for Salamanders and Desmognathus",
    subtitle = "18 fuscus-anchored image-IOD estimates (95% bootstrap intervals) + 10 historical salamander context species",
    x = "Genome Size Estimate (Gb)",
    y = "Shannon Diversity Index",
    color = "Group",
    caption = "Desmognathus: image-IOD estimates conditional on D. fuscus = 16.36 pg; not direct C-values. Historical salamander GS values are cross-study context; no pooled inference."
  ) +
  theme_publication +
  theme(plot.caption = element_text(size = 9, hjust = 0, margin = margin(t = 8)))
save_figure(
  genome_shannon_scatter_plot,
  "te_shannon_genome_size_scatter_te18_v1",
  10,
  6
)

# Desmognathus-only association summaries are exploratory.  The cross-study
# 28-point composite is intentionally not assigned a pooled test because group,
# assay, and repeat-annotation protocol are confounded.
pearson_genome_shannon <- cor.test(
  desmo_genome_scatter$genome_size_gb,
  desmo_genome_scatter$shannon_top10_superfamilies,
  method = "pearson"
)
spearman_genome_shannon <- cor.test(
  desmo_genome_scatter$genome_size_gb,
  desmo_genome_scatter$shannon_top10_superfamilies,
  method = "spearman",
  exact = FALSE
)
ols_genome_shannon <- lm(
  shannon_top10_superfamilies ~ genome_size_gb,
  data = desmo_genome_scatter
)
genome_scatter_tree <- read.tree(paths$tree)
genome_scatter_tree$tip.label <- canonical_species(genome_scatter_tree$tip.label)
genome_scatter_tree <- keep.tip(genome_scatter_tree, desmo_genome_scatter$canonical_species)
pgls_genome_shannon <- gls(
  shannon_top10_superfamilies ~ genome_size_gb,
  data = desmo_genome_scatter,
  correlation = corBrownian(
    value = 1,
    phy = genome_scatter_tree,
    form = ~canonical_species
  ),
  method = "ML"
)
genome_shannon_stats_path <- file.path(
  data_out, "te_shannon_genome_size_stats_te18_v1.txt"
)
writeLines(
  c(
    "Genome-size estimate / Shannon association statistics",
    "Scope: exact 18-species overlap between audited TE34 Shannon and current CellProfiler genome-size estimates",
    "Genome x variable: fuscus-anchored image-IOD-derived estimate converted from pg to Gb with 1 pg = 0.978 Gb",
    "Measurement boundary: these are not direct C-values; uncertainty is conditional on the archived image workflow and D. fuscus = 16.36 pg anchor",
    "Shannon y variable: natural-log entropy of the 10 most abundant TE superfamilies, renormalized within top 10",
    "The 10 historical salamander context points are plotted but excluded from all tests below because study and annotation pipelines differ.",
    "Association tests use species point estimates; they do not propagate genome-estimate or Shannon-estimate uncertainty.",
    "Brownian PGLS uses the single frozen processed TE34 time tree; tree uncertainty is not propagated in this descriptive runner.",
    "",
    "Desmognathus Pearson correlation (non-phylogenetic)",
    capture.output(pearson_genome_shannon),
    "",
    "Desmognathus Spearman correlation (non-phylogenetic)",
    capture.output(spearman_genome_shannon),
    "",
    "Desmognathus ordinary least-squares model",
    capture.output(summary(ols_genome_shannon)),
    "",
    "Desmognathus Brownian-motion PGLS (maximum likelihood)",
    capture.output(summary(pgls_genome_shannon))
  ),
  genome_shannon_stats_path
)

# -----------------------------------------------------------------------------
# Historical PCA and clustering grammar on audited TE34 superfamily composition
# -----------------------------------------------------------------------------

has_variance <- function(x) var(x, na.rm = TRUE) > 0
get_top_features <- function(pca_object, pc_number, top_n = 3L) {
  loadings <- pca_object$rotation[, pc_number]
  ordered <- sort(abs(loadings), decreasing = TRUE)
  paste(names(ordered)[seq_len(min(top_n, length(ordered)))], collapse = ", ")
}

pca_input <- superfamily_comp %>% arrange(match(species, panel_species))
all_superfamily_features <- setdiff(names(pca_input), "species")
explicit_retained_sparse_superfamilies <- c("Chapaev", "Dada")
explicit_excluded_sparse_superfamilies <- c("Ginger")
unresolved_superfamily_features <- intersect("Unknown.TIR", all_superfamily_features)

missing_policy_features <- setdiff(
  c(explicit_retained_sparse_superfamilies, explicit_excluded_sparse_superfamilies),
  all_superfamily_features
)
if (length(missing_policy_features)) {
  stop(
    "TE34 superfamily matrix lacks policy feature(s): ",
    paste(missing_policy_features, collapse = ", ")
  )
}

feature_species_present <- vapply(
  pca_input[, all_superfamily_features, drop = FALSE],
  function(x) sum(is.finite(x) & x > 0),
  integer(1)
)
feature_has_variance <- vapply(
  pca_input[, all_superfamily_features, drop = FALSE],
  has_variance,
  logical(1)
)

pca_features <- all_superfamily_features[
  feature_has_variance &
    !all_superfamily_features %in% unresolved_superfamily_features &
    !all_superfamily_features %in% explicit_excluded_sparse_superfamilies
]
if (!all(explicit_retained_sparse_superfamilies %in% pca_features)) {
  stop("Chapaev and Dada must remain in the TE-superfamily PCA feature set")
}
if (any(explicit_excluded_sparse_superfamilies %in% pca_features)) {
  stop("Ginger must not enter the TE-superfamily PCA feature set")
}

superfamily_prevalence_path <- file.path(data_out, "te_superfamily_prevalence_te34_v1.csv")
superfamily_pca_matrix_path <- file.path(data_out, "te_pca_superfamily_feature_matrix_te34_v1.csv")
superfamily_pca_loadings_path <- file.path(data_out, "te_pca_superfamily_loadings_te34_v1.csv")
superfamily_pca_variance_path <- file.path(data_out, "te_pca_variance_te34_v1.csv")

superfamily_prevalence <- tibble(
  superfamily = all_superfamily_features,
  panel = "TE34",
  panel_species = nrow(pca_input),
  species_present = unname(feature_species_present[all_superfamily_features]),
  prevalence_fraction = species_present / panel_species,
  prevalence_percent = 100 * prevalence_fraction,
  nonzero_variance = unname(feature_has_variance[all_superfamily_features]),
  pca_decision = if_else(superfamily %in% pca_features, "retain", "exclude"),
  decision_reason = case_when(
    superfamily %in% explicit_excluded_sparse_superfamilies ~ "user_declared_sparse_superfamily",
    superfamily %in% unresolved_superfamily_features ~ "unresolved_category",
    !nonzero_variance ~ "zero_variance",
    superfamily %in% explicit_retained_sparse_superfamilies ~ "user_declared_retained_sparse_superfamily",
    TRUE ~ "retained_nonzero_variance"
  )
) %>% arrange(desc(species_present), superfamily)
write_csv(superfamily_prevalence, superfamily_prevalence_path)

pca_feature_matrix <- pca_input %>% select(species, all_of(pca_features))
write_csv(pca_feature_matrix, superfamily_pca_matrix_path)
pca_matrix <- as.matrix(pca_input[, pca_features, drop = FALSE])
rownames(pca_matrix) <- pca_input$species

pca_result <- prcomp(pca_matrix, center = TRUE, scale. = TRUE)
variance_explained <- pca_result$sdev^2 / sum(pca_result$sdev^2)
cumulative_variance <- cumsum(variance_explained)
pca_variance <- tibble(
  PC = seq_along(variance_explained),
  Eigenvalue = pca_result$sdev^2,
  VarExplained = variance_explained,
  CumVarExplained = cumulative_variance
)
write_csv(pca_variance, superfamily_pca_variance_path)

pca_loadings <- as.data.frame(pca_result$rotation, check.names = FALSE) %>%
  mutate(superfamily = rownames(pca_result$rotation), .before = 1)
write_csv(pca_loadings, superfamily_pca_loadings_path)

pc_count_for_clustering <- min(6L, ncol(pca_result$x))
pc_cluster_matrix <- pca_result$x[, seq_len(pc_count_for_clustering), drop = FALSE]
max_k <- min(10L, nrow(pc_cluster_matrix) - 1L)

set.seed(SEED)
kmeans_fits <- lapply(seq_len(max_k), function(k) {
  kmeans(pc_cluster_matrix, centers = k, nstart = 100, iter.max = 100)
})
wss <- vapply(kmeans_fits, function(fit) fit$tot.withinss, numeric(1))
silhouette_scores <- rep(NA_real_, max_k)
for (k in 2:max_k) {
  silhouette_scores[[k]] <- mean(silhouette(kmeans_fits[[k]]$cluster, dist(pc_cluster_matrix))[, 3])
}
selected_k <- which.max(replace(silhouette_scores, is.na(silhouette_scores), -Inf))
selected_fit <- kmeans_fits[[selected_k]]

cluster_metrics <- tibble(
  k = seq_len(max_k),
  wss = wss,
  # k=1 has no mathematical silhouette. Export zero rather than an NA token so
  # the metric column remains machine-numeric; `silhouette_defined` preserves
  # the distinction and selection still considers only k >= 2 above.
  mean_silhouette = replace(silhouette_scores, is.na(silhouette_scores), 0),
  silhouette_defined = seq_len(max_k) >= 2L,
  selected = seq_len(max_k) == selected_k
)
write_csv(cluster_metrics, file.path(data_out, "te_pca_cluster_metrics_te34_v1.csv"))

scores <- as.data.frame(pca_result$x) %>%
  mutate(
    species = rownames(pca_result$x),
    display_species = display_species(species),
    cluster = factor(selected_fit$cluster),
    selected_k = selected_k,
    .before = 1
  )
write_csv(scores, file.path(data_out, "te_pca_scores_clusters_te34_v1.csv"))

# Literal historical PCA species plot with audited labels and no species legend.
pca_species_plot <- ggplot(scores, aes(x = PC1, y = PC2, color = display_species)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_text_repel(aes(label = display_species), size = 3, max.overlaps = Inf, seed = SEED) +
  labs(
    title = "PCA Plot TE Composition by Species",
    subtitle = "Audited TE34 superfamily composition",
    x = paste0(
      "PC1 (", round(100 * variance_explained[1], 1), "% variance)\nTop features: ",
      get_top_features(pca_result, 1)
    ),
    y = paste0(
      "PC2 (", round(100 * variance_explained[2], 1), "% variance)\nTop features: ",
      get_top_features(pca_result, 2)
    )
  ) +
  theme_publication +
  theme(legend.position = "none")
save_figure(pca_species_plot, "te_pca_species_te34_v1", 10, 8)

pca_scree_plot <- ggplot() +
  geom_line(data = pca_variance, aes(x = PC, y = VarExplained), color = "blue") +
  geom_point(data = pca_variance, aes(x = PC, y = VarExplained), color = "blue") +
  geom_line(data = pca_variance, aes(x = PC, y = CumVarExplained), color = "red", linetype = "dashed") +
  geom_hline(yintercept = c(0.8, 0.9), color = "darkgreen", linetype = "dashed") +
  scale_y_continuous(
    name = "Variance Explained",
    labels = label_percent(accuracy = 1),
    sec.axis = sec_axis(~., name = "Cumulative Variance Explained", labels = label_percent(accuracy = 1))
  ) +
  labs(title = "Scree Plot with Cumulative Variance", x = "Principal Component") +
  theme_publication +
  theme(legend.position = "none") +
  annotate(
    "text", x = max(pca_variance$PC) / 2, y = c(0.82, 0.92),
    label = c("80% threshold", "90% threshold"), color = "darkgreen", family = "serif"
  )
save_figure(pca_scree_plot, "te_pca_scree_te34_v1", 10, 8)

pca_elbow_plot <- ggplot(cluster_metrics, aes(x = k, y = wss)) +
  geom_line() +
  geom_point() +
  geom_vline(xintercept = selected_k, linetype = "dashed", color = "darkgreen") +
  scale_x_continuous(breaks = seq_len(max_k)) +
  labs(
    title = "Elbow Plot for Optimal k",
    subtitle = paste0("Dashed line: maximum-silhouette k = ", selected_k),
    x = "Number of Clusters (k)", y = "Total Within-cluster Sum of Squares"
  ) +
  theme_publication +
  theme(legend.position = "none")
save_figure(pca_elbow_plot, "te_pca_elbow_te34_v1", 10, 8)

pca_silhouette_plot <- ggplot(filter(cluster_metrics, k >= 2L), aes(x = k, y = mean_silhouette)) +
  geom_line() +
  geom_point() +
  geom_vline(xintercept = selected_k, linetype = "dashed", color = "darkgreen") +
  scale_x_continuous(breaks = 2:max_k) +
  labs(
    title = "Silhouette Scores",
    subtitle = paste0("Selected k = ", selected_k, " (maximum mean silhouette)"),
    x = "Number of Clusters (k)", y = "Average Silhouette Score"
  ) +
  theme_publication +
  theme(legend.position = "none")
save_figure(pca_silhouette_plot, "te_pca_silhouette_te34_v1", 10, 8)

cluster_counts <- scores %>% count(cluster, name = "cluster_n")
singleton_scores <- scores %>%
  left_join(cluster_counts, by = "cluster") %>%
  filter(cluster_n == 1L)
ellipse_scores <- scores %>% add_count(cluster) %>% filter(n >= 3L)
pca_cluster_plot <- ggplot(scores, aes(x = PC1, y = PC2, color = cluster)) +
  geom_point(size = 3, alpha = 0.6) +
  stat_ellipse(data = ellipse_scores, aes(group = cluster), level = 0.95) +
  geom_text_repel(
    data = singleton_scores, aes(label = display_species),
    seed = SEED, show.legend = FALSE, family = "serif", fontface = "italic"
  ) +
  labs(
    title = "PCA Plot with Clusters",
    subtitle = paste0(
      "Audited TE34; silhouette-selected k = ", selected_k,
      if (nrow(singleton_scores)) paste0("; singleton outlier: ", paste(singleton_scores$display_species, collapse = ", ")) else ""
    ),
    x = paste0("PC1 (", round(100 * variance_explained[1], 1), "%)"),
    y = paste0("PC2 (", round(100 * variance_explained[2], 1), "%)"),
    color = "Cluster"
  ) +
  theme_publication +
  scale_color_viridis_d(direction = -1)
save_figure(pca_cluster_plot, "te_pca_clusters_te34_v1", 10, 8)

# Circular companion requested for the audited cluster assignments.
tree <- read.tree(paths$tree)
tree$tip.label <- canonical_species(tree$tip.label)
if (anyDuplicated(tree$tip.label)) stop("Processed TE34 tree has duplicate canonical tips")
if (!setequal(tree$tip.label, panel_species)) stop("Processed tree does not match exact TE34 panel")
tree <- keep.tip(tree, panel_species)
cluster_lookup <- setNames(as.character(scores$cluster), scores$species)
tip_data <- tibble(
  label = tree$tip.label,
  display_label = display_species(tree$tip.label),
  cluster = factor(cluster_lookup[tree$tip.label], levels = levels(scores$cluster))
)

tree_plot <- ggtree(tree, layout = "circular", size = 0.45) %<+% tip_data +
  geom_tippoint(aes(color = cluster), size = 3, alpha = 0.9) +
  geom_tiplab(aes(label = display_label, color = cluster), size = 2.7, offset = 0.45, show.legend = FALSE) +
  scale_color_viridis_d(direction = -1) +
  labs(
    title = "TE-composition clusters across the Desmognathus phylogeny",
    subtitle = paste0("Audited TE34 time tree; silhouette-selected k = ", selected_k),
    color = "Cluster"
  ) +
  theme_void(base_family = "serif") +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5),
    legend.position = "right",
    plot.margin = margin(15, 30, 15, 30)
  )
save_figure(tree_plot, "te_pca_clusters_phylogeny_te34_v1", 12, 12)

# -----------------------------------------------------------------------------
# Historical unfiltered LTR violin plus readable log-scale companion (LTR30)
# -----------------------------------------------------------------------------

ltr <- read_csv(paths$ltr_elements, show_col_types = FALSE) %>%
  transmute(
    species = canonical_species(species),
    display_species = display_species(species),
    element_id,
    ratio_terminal_internal = as.numeric(ratio_terminal_internal_all_positions),
    left_right_terminal_log2_imbalance = as.numeric(left_right_terminal_log2_imbalance),
    terminal_positive_coverage_fraction = as.numeric(terminal_positive_coverage_fraction),
    internal_positive_coverage_fraction = as.numeric(internal_positive_coverage_fraction),
    domain_count = as.integer(domain_count),
    complete = Complete,
    sequence = sequence
  ) %>%
  filter(is.finite(ratio_terminal_internal), ratio_terminal_internal > 0)

ltr_exclusions <- read_csv(paths$ltr_exclusions, show_col_types = FALSE)
ltr_source_manifest <- fromJSON(paths$ltr_manifest, simplifyVector = FALSE)
if (!identical(ltr_source_manifest$n_selected_elements, 1088L) ||
    !identical(ltr_source_manifest$n_source_corrupt_exclusions, 2L) ||
    !identical(ltr_source_manifest$n_usable_elements, 1086L) ||
    !isTRUE(ltr_source_manifest$zero_depth_positions_retained)) {
  stop("Corrected LTR30 source manifest no longer satisfies the audited 1,088 - 2 = 1,086 contract")
}
if (nrow(ltr_exclusions) != 2L) stop("Expected exactly two proven-truncated LTR depth-file exclusions")

ltr_species <- sort(unique(ltr$species))
if (length(ltr_species) != 30L) stop("Expected the complete historical-style LTR resource table to contain 30 species")
if (!all(ltr_species %in% panel_species)) stop("LTR30 includes species outside audited TE34")
if (nrow(ltr) != 1086L) stop("Unexpected corrected LTR30 element count; audited input may have changed")

# Jake's historical presentation intentionally screened non-biological outliers.
# Keep the complete corrected table untouched, and make the conservative mapping-
# artifact screen a separate, row-auditable branch. Statistical extremeness alone
# is not enough: an element must be an upper within-species 1.5-IQR outlier, show
# >64-fold left/right terminal imbalance, and contribute >50% of its species'
# total element-level ratio.
LTR_TERMINAL_ABS_LOG2_IMBALANCE_THRESHOLD <- 6
ltr_screen <- ltr %>%
  group_by(species) %>%
  mutate(
    species_q1 = quantile(ratio_terminal_internal, 0.25),
    species_q3 = quantile(ratio_terminal_internal, 0.75),
    species_iqr = IQR(ratio_terminal_internal),
    species_iqr_upper = species_q3 + 1.5 * species_iqr,
    fraction_of_species_ratio_sum = ratio_terminal_internal / sum(ratio_terminal_internal),
    upper_within_species_iqr_outlier = ratio_terminal_internal > species_iqr_upper,
    extreme_one_sided_terminal_pileup =
      abs(left_right_terminal_log2_imbalance) > LTR_TERMINAL_ABS_LOG2_IMBALANCE_THRESHOLD,
    dominates_species_ratio_sum = fraction_of_species_ratio_sum > 0.5,
    artifact_screen_keep = !(
      upper_within_species_iqr_outlier &
        extreme_one_sided_terminal_pileup &
        dominates_species_ratio_sum
    ),
    artifact_screen_rule_id = "one_sided_terminal_pileup_v1",
    artifact_screen_reason = if_else(
      artifact_screen_keep,
      "retained",
      "flagged: upper 1.5-IQR outlier; >64-fold left/right terminal imbalance; >50% of species ratio sum"
    )
  ) %>%
  ungroup()

ltr_flagged <- ltr_screen %>% filter(!artifact_screen_keep)
ltr_screened <- ltr_screen %>% filter(artifact_screen_keep)
expected_artifact_element <- "JAUEJH010597481.1_De_5169_11831"
if (nrow(ltr_flagged) != 1L ||
    ltr_flagged$element_id[[1]] != expected_artifact_element ||
    nrow(ltr_screened) != 1085L ||
    n_distinct(ltr_screened$species) != 30L ||
    max(ltr_screened$ratio_terminal_internal) >= 14) {
  stop("Conservative LTR mapping-artifact screen no longer matches its audited one-element contract")
}

ltr_artifact_screen_path <- file.path(data_out, "ltr_mapping_artifact_screen_ltr30_v1.csv")
write_csv(
  ltr_screen %>%
    rename(ratio_terminal_internal_all_positions = ratio_terminal_internal),
  ltr_artifact_screen_path
)

ltr$display_species <- factor(ltr$display_species, levels = display_species(ltr_species))
ltr_screened$display_species <- factor(
  ltr_screened$display_species,
  levels = display_species(ltr_species)
)

create_violin_plot <- function(data, title, subtitle = NULL) {
  ggplot(data, aes(x = display_species, y = ratio_terminal_internal, fill = display_species)) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "red", linewidth = 0.8) +
    geom_violin(trim = FALSE, alpha = 0.8) +
    stat_summary(fun = mean, geom = "crossbar", width = 0.5, color = "black", linewidth = 0.6) +
    labs(
      title = title, subtitle = subtitle,
      x = "Species", y = "Ratio of LTR to Internal"
    ) +
    theme_ltr_publication +
    scale_fill_viridis_d() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 9),
      legend.position = "none"
    )
}

ltr_raw_plot <- create_violin_plot(
  ltr,
  "LTR Elements with 5+ Domains: Historical-style Distribution",
  "LTR30; 1,086 usable elements after two proven-truncated depth files; no IQR deletion"
)
save_figure(ltr_raw_plot, "ltr_terminal_internal_violin_ltr30_v1", 10, 8)

ltr_artifact_screened_plot <- create_violin_plot(
  ltr_screened,
  "LTR Elements with 5+ Domains: Artifact-screened Distribution",
  "LTR30; 1,085 retained elements; one extreme one-sided terminal pileup is preserved in the audit table"
)
save_figure(
  ltr_artifact_screened_plot,
  "ltr_terminal_internal_artifact_screened_ltr30_v1",
  10,
  8
)

ltr_log_plot <- create_violin_plot(
  ltr,
  "LTR Elements with 5+ Domains: Log-scale Distribution",
  "Corrected all-position, zero-aware depth estimator; log axis reveals the central distributions"
) +
  scale_y_log10(
    breaks = c(0.03, 0.1, 0.3, 1, 3, 10, 30),
    labels = label_number(accuracy = 0.01)
  )
save_figure(ltr_log_plot, "ltr_terminal_internal_log_ltr30_v1", 10, 8)

raw_anova <- aov(ratio_terminal_internal ~ species, data = ltr)
log_anova <- aov(log10(ratio_terminal_internal) ~ species, data = ltr)
raw_kruskal <- kruskal.test(ratio_terminal_internal ~ species, data = ltr)
screened_raw_anova <- aov(ratio_terminal_internal ~ species, data = ltr_screened)
screened_log_anova <- aov(log10(ratio_terminal_internal) ~ species, data = ltr_screened)
screened_raw_kruskal <- kruskal.test(ratio_terminal_internal ~ species, data = ltr_screened)
stats_path <- file.path(data_out, "ltr_terminal_internal_stats_v1.txt")
stats_lines <- c(
  "LTR terminal:internal depth-ratio statistics",
  "Analysis scope: LTR30 audited TE-resource subset with available >=3-kb, >=5-domain elements",
  paste0("Elements: ", nrow(ltr)),
  paste0("Species: ", n_distinct(ltr$species)),
  "Selected elements: 1,088; exclusions: 2 depth files proven truncated mid-final-line; usable elements: 1,086",
  "Unfiltered positive finite element audit: 1,086 elements after two source-corrupt exclusions; NO IQR deletion.",
  "Artifact-screened presentation: 1,085 retained elements. It does not apply blanket IQR deletion.",
  "Screen rule: flag only elements that are simultaneously an upper within-species 1.5-IQR outlier, have absolute left/right terminal log2 imbalance > 6 (>64-fold), and contribute >50% of the species ratio sum.",
  paste0("Flagged element: ", expected_artifact_element, "; corrected ratio = ", signif(ltr_flagged$ratio_terminal_internal[[1]], 6), "; left/right log2 imbalance = ", signif(ltr_flagged$left_right_terminal_log2_imbalance[[1]], 6)),
  "Estimator: mean terminal depth across all expected positions / mean internal depth across all expected positions; explicit zero-depth positions are retained.",
  "Interpretation caveat: terminal:internal depth is a deletion footprint proxy, not a directly validated ectopic-recombination rate.",
  "",
  "Raw one-way ANOVA: ratio_terminal_internal ~ species",
  capture.output(summary(raw_anova)),
  "",
  "Log10 one-way ANOVA: log10(ratio_terminal_internal) ~ species",
  capture.output(summary(log_anova)),
  "",
  "Kruskal-Wallis raw ratio: ratio_terminal_internal ~ species",
  capture.output(raw_kruskal),
  "",
  "Artifact-screened raw one-way ANOVA: ratio_terminal_internal ~ species",
  capture.output(summary(screened_raw_anova)),
  "",
  "Artifact-screened log10 one-way ANOVA: log10(ratio_terminal_internal) ~ species",
  capture.output(summary(screened_log_anova)),
  "",
  "Artifact-screened Kruskal-Wallis raw ratio: ratio_terminal_internal ~ species",
  capture.output(screened_raw_kruskal),
  "",
  "Tukey HSD after raw one-way ANOVA",
  capture.output(TukeyHSD(raw_anova)),
  "",
  "Tukey HSD after log10 one-way ANOVA",
  capture.output(TukeyHSD(log_anova))
)
writeLines(stats_lines, stats_path)

# -----------------------------------------------------------------------------
# Machine-auditable manifest
# -----------------------------------------------------------------------------

input_row_counts <- list(
  panel = nrow(panel),
  integrated_panel = nrow(integrated_panel),
  genome_estimates = nrow(genome_estimates),
  genome_state = nrow(genome_state),
  retry_manifest = 1L,
  mass = nrow(mass),
  quality = nrow(quality),
  order_composition = nrow(order_comp),
  superfamily_composition = nrow(superfamily_comp),
  diversity = nrow(diversity),
  wang_ji_context = nrow(wang_context),
  wang_ji_provenance = 1L,
  ltr_elements = nrow(ltr),
  ltr_exclusions = nrow(ltr_exclusions),
  ltr_manifest = 1L,
  tree = length(tree$tip.label)
)

input_manifest <- lapply(names(paths), function(name) {
  list(
    id = name,
    path = relative_path(paths[[name]]),
    sha256 = sha256(paths[[name]]),
    rows_or_tips = input_row_counts[[name]]
  )
})

manifest <- list(
  analysis_id = "audited_historical_style_figures_v1",
  runner = "scripts/processing/build_audited_historical_style_figures.R",
  rebuild_mode = "Current audited and clean frozen data rendered with Jake's historical R and ggplot grammar; no embedded historical PNGs.",
  generated_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  expensive_upstream_tools_executed = FALSE,
  historical_code_provenance = list(
    assembly_and_mass = "/home/jake/Projects/misc_desmog_scripts/Master_data_vis.R",
    pca_git = "1ad5e8c:te_pca_analysis.R",
    pca_original_git = "89c0879:Master_data_vis_PCA.R",
    pca_vscode_history = "DWDT.R",
    composition_git = "1ad5e8c:scripts/visualization/te_landscape_plots.R",
    shannon = "/home/jake/Projects/misc_desmog_scripts/Master_shannon_diversity_vis.R",
    ectopic_git = "1ad5e8c:scripts/visualization/ectopic_recomb_plots.R",
    ectopic_original_git = "89c0879:Master_data_vis2_ectopicrecomb.R"
  ),
  scopes = list(
    te = list(panel = "TE34", n_species = 34L),
    ltr = list(
      panel = "LTR30 available-resource subset",
      n_species = 30L,
      n_elements = nrow(ltr),
      n_artifact_screened_elements = nrow(ltr_screened)
    ),
    cross_taxon_shannon = list(
      role = "descriptive cross-study context only",
      n_historical_context_rows = sum(!historical_desmognathus_rows),
      n_current_desmognathus = nrow(desmo_context),
      n_total = nrow(cross_taxon)
    ),
    genome_size_shannon_scatter = list(
      panel = "exact TE-genome overlap / path18",
      n_current_desmognathus = nrow(desmo_genome_scatter),
      n_historical_salamanders = nrow(historical_salamander_scatter),
      n_total = nrow(genome_shannon_scatter)
    ),
    integrated_path_analysis = "not performed by this figure runner"
  ),
  species_identity = list(
    fuscus_analysis_label = "fuscus",
    fuscus_te_resource = "GCA_032353935.1",
    ncbi_public_label_for_resource = "Desmognathus planiceps",
    mapping_basis = "expert reidentification supplied by Alex Pyron; mapping applies to the genomic resource only"
  ),
  dnapipete_retry_policy = list(
    method_id = "equal_weight_mean_of_run_level_estimates",
    method = "equal-weight mean of run-level estimates within species",
    retries_are_additional_species = FALSE,
    species_with_retry = "orestes",
    source_labels = orestes_labels
  ),
  denominators = list(
    dnapipete_repeat_aligned = "total aligned repeat bases / configured 1,500,000,000-bp quantification sample",
    dnapipete_order_classification = "classified or unresolved aligned bases / total repeat-aligned bases",
    dnapipete_superfamily_classification = "classified or unresolved aligned bases / total repeat-aligned bases",
    mean_order_and_superfamily = "within-species classified-conditional composition, then equal species mean",
    diversity = "both classified_conditional and mass_aware_unresolved_bin modes shown",
    cross_taxon_diversity = "published top-10-superfamily natural-log Shannon contract; top 10 are renormalized to sum to one within species",
    ltr = "all-position mean terminal depth / all-position mean internal depth; explicit zero-depth positions retained"
  ),
  superfamily_pca_filter = list(
    panel = "TE34",
    panel_species = nrow(pca_input),
    data = "dnaPipeTE retry-averaged species-level superfamily relative composition",
    historical_user_decision_context = list(
      panel_species = 35L,
      prevalence = list(
        Chapaev = list(species_present = 14L, prevalence_percent = 40),
        Dada = list(species_present = 10L, prevalence_percent = 28.57),
        Ginger = list(species_present = 1L, prevalence_percent = 2.86)
      ),
      note = "Historical 35-species counts document the original filtering decision; they are not relabeled as TE34 counts."
    ),
    current_audited_prevalence = lapply(
      c(explicit_retained_sparse_superfamilies, explicit_excluded_sparse_superfamilies),
      function(feature) {
        row <- filter(superfamily_prevalence, superfamily == feature)
        list(
          superfamily = feature,
          species_present = row$species_present[[1]],
          panel_species = row$panel_species[[1]],
          prevalence_percent = row$prevalence_percent[[1]],
          pca_decision = row$pca_decision[[1]],
          decision_reason = row$decision_reason[[1]]
        )
      }
    ),
    retained_by_explicit_policy = explicit_retained_sparse_superfamilies,
    excluded_by_explicit_policy = I(explicit_excluded_sparse_superfamilies),
    other_exclusions = unresolved_superfamily_features,
    variance_rule = "exclude zero-variance features; do not prevalence-filter Chapaev or Dada",
    prevalence_audit = list(
      path = relative_path(superfamily_prevalence_path),
      sha256 = sha256(superfamily_prevalence_path),
      rows = nrow(superfamily_prevalence)
    ),
    exact_pca_matrix = list(
      path = relative_path(superfamily_pca_matrix_path),
      sha256 = sha256(superfamily_pca_matrix_path),
      species_rows = nrow(pca_feature_matrix),
      retained_features = length(pca_features)
    )
  ),
  pca_and_clustering = list(
    data = "TE34 retry-averaged superfamily relative composition",
    transform = "historical scaled raw relative composition",
    feature_rule = "retain Chapaev and Dada; exclude Ginger; remove Unknown.TIR if present and all zero-variance features",
    n_features = length(pca_features),
    clustering_axes = pc_count_for_clustering,
    deterministic_seed = SEED,
    k_range = paste0("1-", max_k, "; silhouette evaluated for 2-", max_k),
    selection_rule = "maximum mean silhouette width",
    selected_k = selected_k,
    selected_mean_silhouette = silhouette_scores[[selected_k]]
  ),
  cross_taxon_shannon = list(
    historical_plot_grammar = "literal boxplot statement from Master_shannon_diversity_vis.R",
    current_desmognathus_data = "audited TE34 retry-averaged classified-conditional superfamily composition",
    metric = "natural-log Shannon entropy over each species' 10 most abundant TE superfamilies, renormalized within the top 10",
    context_source = "87-row historical published-vertebrate composite used by prior project code",
    context_source_sha256 = WANG_JI_SOURCE_SHA256,
    published_data_repository = paste0("Dryad DOI ", WANG_JI_DRYAD_DOI),
    context_scope_disclosure = "Not an exact export of Wang et al. 2025's 84-species analysis; includes 4 caecilians and 10 salamanders and is used only as descriptive cross-study context.",
    inference = "none for the pooled cross-study context; Desmognathus-only Pearson, Spearman, OLS, and Brownian PGLS are exported as exploratory summaries",
    genome_size_scatter = list(
      historical_plot_grammar = "literal point/color/lab/theme statement from Master_shannon_diversity_vis.R, with horizontal bootstrap intervals added for current Desmognathus estimates",
      n_current_desmognathus = nrow(desmo_genome_scatter),
      n_historical_salamanders = nrow(historical_salamander_scatter),
      overlap_rule = "exact intersection of current TE34 top-10-superfamily Shannon and the 21-species frozen fuscus-anchored genome estimate table; equals path18",
      calibration_reference_species = "D. fuscus",
      calibration_reference_pg = 16.36,
      pg_to_gb_factor = 0.978,
      measurement_kind = "fuscus_anchored_image_iod_genome_size_estimate",
      absolute_c_value_claimed = FALSE,
      measurement_boundary = "Desmognathus values are fuscus-anchored image-IOD genome-size estimates with conditional bootstrap intervals, not direct C-values; historical salamander values remain descriptive cross-study context.",
      pooled_cross_study_inference = "none",
      data = list(
        scatter = list(path = relative_path(genome_shannon_scatter_path), sha256 = sha256(genome_shannon_scatter_path), rows = nrow(genome_shannon_scatter)),
        match_audit = list(path = relative_path(genome_shannon_match_audit_path), sha256 = sha256(genome_shannon_match_audit_path), rows = nrow(genome_shannon_match_audit)),
        statistics = list(path = relative_path(genome_shannon_stats_path), sha256 = sha256(genome_shannon_stats_path))
      )
    ),
    data_artifacts = list(
      historical_context = list(path = relative_path(wang_ji_snapshot), sha256 = sha256(wang_ji_snapshot), rows = nrow(wang_context)),
      provenance = list(path = relative_path(wang_ji_provenance), sha256 = sha256(wang_ji_provenance), rows = 1L),
      combined_context = list(path = relative_path(cross_taxon_path), sha256 = sha256(cross_taxon_path), rows = nrow(cross_taxon)),
      group_counts = list(path = relative_path(cross_taxon_counts_path), sha256 = sha256(cross_taxon_counts_path), rows = nrow(cross_taxon_counts)),
      desmognathus_top10_contributions = list(path = relative_path(desmo_top10_long_path), sha256 = sha256(desmo_top10_long_path), rows = nrow(desmo_top10_long))
    )
  ),
  ltr_mapping_artifact_screen = list(
    role = "presentation and sensitivity branch; unfiltered corrected input remains intact",
    rule_id = "one_sided_terminal_pileup_v1",
    n_input_elements = nrow(ltr_screen),
    n_flagged_elements = nrow(ltr_flagged),
    n_retained_elements = nrow(ltr_screened),
    n_retained_species = n_distinct(ltr_screened$species),
    left_right_abs_log2_threshold = LTR_TERMINAL_ABS_LOG2_IMBALANCE_THRESHOLD,
    additional_required_conditions = c(
      "upper within-species 1.5-IQR outlier",
      "element contributes more than 50% of species ratio sum"
    ),
    flagged_element_ids = ltr_flagged$element_id,
    unfiltered_input_modified = FALSE,
    audit_table = list(
      path = relative_path(ltr_artifact_screen_path),
      sha256 = sha256(ltr_artifact_screen_path),
      rows = nrow(ltr_screen)
    )
  ),
  methodological_deviations_from_historical_scripts = c(
    "The vetted current TE34 replaces the historical unvetted species set.",
    "dnaPipeTE retries are averaged within species before species-level analysis.",
    "The original sparse-feature decision is explicit and audited on TE34: Chapaev and Dada are retained, while Ginger is excluded.",
    "A deterministic seed and nstart=100 replace unseeded or lightly restarted historical k-means.",
    "k is selected by maximum mean silhouette instead of hard-coding k=5.",
    "Stacked mass bars use mathematically correct stacking and explicitly named denominators.",
    "All PCA species labels are retained with deterministic repulsion for readability.",
    "The discrete viridis cluster palette is reversed so the 33-species majority cluster remains legible on white while the singleton outlier stays visually distinct.",
    "The corrected TE34 denominator-sensitivity panel is retained as an audit companion; the historical cross-taxon boxplot is separately restored with the published top-10-superfamily metric contract.",
    "The assembly panel uses the current audited accessions and explicitly maps the public planiceps row to the expert-reidentified fuscus genomic resource.",
    "The complete 1,086-element corrected LTR audit remains unfiltered and is shown on raw and log scales.",
    "A separate presentation branch screens one row-auditable one-sided terminal pileup only when it also exceeds the within-species upper 1.5-IQR fence and contributes over half the species ratio sum; blanket IQR deletion is not used.",
    "The LTR30 figures use the audited all-position, zero-aware estimator and exclude exactly two proven-truncated depth files.",
    "The restored cross-taxon Shannon boxplot replaces the stale historical Desmognathus rows with current audited TE34 top-10-superfamily values.",
    "The genome-size/Shannon scatter replaces the historical all-16 Desmognathus placeholder with the exact 18-species TE/genome overlap and retains the historical point/color/lab/theme grammar.",
    "Current Desmognathus x values are visibly identified as fuscus-anchored image-IOD genome-size estimates with bootstrap intervals, not direct C-values."
  ),
  inference_limits = c(
    "PCA clustering is descriptive and is not a phylogenetically corrected hypothesis test.",
    "The cross-taxon Shannon boxplot is descriptive context only; heterogeneous repeat libraries, annotation protocols, sampling, and strongly imbalanced group sizes preclude inferential comparison.",
    "The pooled 28-point genome-size/Shannon scatter is descriptive only because current fuscus-anchored image-IOD estimates and historical salamander genome-size values come from non-homogeneous assays and TE annotation workflows.",
    "Desmognathus-only genome-size/Shannon tests are exploratory and conditional on the fuscus calibration anchor, archived image workflow, and selected TE diversity definition.",
    "Genome-size/Shannon association tests use point estimates and one frozen time tree; measurement and tree uncertainty are not propagated.",
    "Element-level LTR ANOVA/Tukey output is historical-style exploratory output and does not remove pseudoreplication or establish an ectopic-recombination rate.",
    "The LTR artifact screen is a conservative mapping-quality sensitivity rule, not proof that every statistically extreme element is non-biological.",
    "No confirmatory phylogenetic path inference is performed in this runner."
  ),
  inputs = input_manifest,
  data_outputs = list(
    relative_path(assembly_csv),
    relative_path(file.path(data_out, "dnapipete_quality_te34_v1.csv")),
    relative_path(file.path(data_out, "te_pca_cluster_metrics_te34_v1.csv")),
    relative_path(file.path(data_out, "te_pca_scores_clusters_te34_v1.csv")),
    relative_path(superfamily_prevalence_path),
    relative_path(superfamily_pca_matrix_path),
    relative_path(superfamily_pca_loadings_path),
    relative_path(superfamily_pca_variance_path),
    relative_path(wang_ji_snapshot),
    relative_path(wang_ji_provenance),
    relative_path(cross_taxon_path),
    relative_path(cross_taxon_counts_path),
    relative_path(desmo_top10_long_path),
    relative_path(genome_shannon_scatter_path),
    relative_path(genome_shannon_match_audit_path),
    relative_path(genome_shannon_stats_path),
    relative_path(ltr_artifact_screen_path),
    relative_path(stats_path)
  ),
  companion_figure_references = list(
    list(
      role = "corrected LTR30 coverage and source-QC companion",
      path = "results/figures/research_review/ectopic_ltr30_coverage_qc_v1.png",
      sha256 = ltr_source_manifest$figures[[1]]$sha256
    )
  ),
  figures = unname(figure_specs)
)

manifest_path <- file.path(data_out, "audited_historical_style_figures_v1.manifest.json")
write_json(manifest, manifest_path, pretty = TRUE, auto_unbox = TRUE, null = "null", digits = NA)

expected_stems <- c(
  "assembly_quality_te34_v1",
  "dnapipete_quality_te34_v1",
  "te_mass_composition_te34_v1",
  "te_mean_order_te34_v1",
  "te_mean_superfamily_te34_v1",
  "te_diversity_te34_v1",
    "te_shannon_cross_taxon_context_te34_v1",
    "te_shannon_genome_size_scatter_te18_v1",
  "te_pca_species_te34_v1",
  "te_pca_scree_te34_v1",
  "te_pca_elbow_te34_v1",
  "te_pca_silhouette_te34_v1",
  "te_pca_clusters_te34_v1",
  "te_pca_clusters_phylogeny_te34_v1",
    "ltr_terminal_internal_violin_ltr30_v1",
    "ltr_terminal_internal_artifact_screened_ltr30_v1",
    "ltr_terminal_internal_log_ltr30_v1"
)
if (!identical(sort(names(figure_specs)), sort(expected_stems))) stop("Figure contract is incomplete")

message("Audited historical-style figure rebuild complete")
message("  TE panel: 34 species")
message("  LTR panel: ", n_distinct(ltr$species), " species / ", nrow(ltr), " elements")
message("  PCA features: ", length(pca_features), "; selected k: ", selected_k)
message("  Figures: ", length(figure_specs), " stems x PNG/PDF")
message("  Manifest: ", relative_path(manifest_path))
