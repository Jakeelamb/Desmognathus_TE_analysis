validate_auxiliary_figure_manifest <- function(manifest) {
  required_columns <- c(
    "auxiliary_id", "title", "basename", "width_mm", "height_mm",
    "final_font_family", "minimum_target_text_pt", "datasets",
    "release_status", "boundary", "legend", "alt_text", "pdf", "png", "tiff"
  )
  if (!all(required_columns %in% names(manifest))) {
    stop("Auxiliary figure manifest is missing required columns: ",
         paste(setdiff(required_columns, names(manifest)), collapse = ", "))
  }
  ids <- c(
    "Auxiliary_Wang_style_integrated_te_history_candidate",
    "Auxiliary_Wang_style_diversity_genome_size_candidate",
    "Auxiliary_Path24_development_traits_reference",
    "Auxiliary_Path_DAG_reference",
    "Auxiliary_Wang2025_cross_study_te_diversity_reference"
  )
  if (!identical(manifest$auxiliary_id, ids)) {
    stop("Auxiliary figure manifest must contain exactly the five ordered IDs")
  }
  retained_titles <- c(
    "Wang-style integrated TE content and proliferation candidate",
    "Wang-style TE diversity and genome-size candidate",
    "Path24 development and phenotype reference",
    "Path DAG equivalence-class reference",
    "Cross-study TE-diversity reference (Wang et al. 2025)"
  )
  current_titles <- retained_titles
  current_titles[c(1, 3)] <- c(
    "Wang-style integrated TE history across all 34 species",
    "Path24 phylogeny and phenotype reference"
  )
  if (identical(manifest$title, retained_titles)) {
    variant <- "retained"
  } else if (identical(manifest$title, current_titles)) {
    variant <- "current"
  } else {
    stop("Auxiliary figure manifest titles do not match a known ordered variant")
  }
  for (column in c("width_mm", "height_mm", "minimum_target_text_pt")) {
    if (any(!is.finite(manifest[[column]])) || any(manifest[[column]] <= 0)) {
      stop("Auxiliary figure manifest requires positive dimensions and text sizes")
    }
  }
  for (column in setdiff(required_columns, c("width_mm", "height_mm", "minimum_target_text_pt"))) {
    if (anyNA(manifest[[column]]) || any(!nzchar(trimws(manifest[[column]])))) {
      stop("Auxiliary figure manifest has empty ", column)
    }
  }
  phrases <- list(
    c("Review candidate, not a numbered manuscript figure",
      "fractions of aligned RepeatMasker hit bp", "not genome occupancy",
      "divergence bins are not insertion ages"),
    c("Review candidate, not a numbered manuscript figure",
      "Reproductive-strategy panels are post-hoc descriptive strata",
      "without a regression because n = 2", "Pagel-lambda PGLS fits"),
    c("Descriptive review only", "D. organi is absent", "not imputed",
      "D. fuscus assembly-anchored pg/1C estimates",
      "not direct independently measured species C-values"),
    c("Not a numbered manuscript figure.", "10 testable equivalence classes",
      "saturated 6-member class is not shown",
      "remains unscored because it has no d-separation claim", "Genome size",
      "D. fuscus assembly-anchored pg/1C",
      "not a direct species C-value measurement",
      "same released genome-size scale is used consistently",
      "do not orient a unique causal arrow"),
    c("Descriptive cross-study context only", "no pooled regression or hypothesis test",
      "10.5061/dryad.zpc866tkv", "10.1038/s42003-025-08127-3", "CC0-1.0",
      "published top-10 TE-superfamily indices", "RepeatModeler/RepeatMasker",
      "reclosed calculations from dnaPipeTE compositions",
      "pipelines are different and not method-harmonized",
      "D. fuscus assembly-anchored estimates", "not direct species C-values")
  )
  if (variant == "retained") {
    phrases[[1]] <- c(phrases[[1]], "Six representatives are selected")
    phrases[[4]] <- c(phrases[[4]], "conditional D. fuscus assembly-anchored pg/1C analysis")
    phrases[[5]] <- c(phrases[[5]], "conditional D. fuscus assembly-anchored estimates")
  } else {
    phrases[[1]] <- c(phrases[[1]], "All 34 TE species are included",
                      "available for 21, with 13 labelled n/a", "intervals omit anchor uncertainty")
    phrases[[3]] <- c(phrases[[3]], "The tree shows topology only")
    phrases[[4]] <- c(phrases[[4]], "D. fuscus assembly-anchored pg/1C scale")
  }
  for (index in seq_along(ids)) {
    evidence <- paste(manifest$boundary[index], manifest$legend[index], manifest$alt_text[index])
    missing <- phrases[[index]][!vapply(phrases[[index]], grepl, logical(1),
                                        x = evidence, fixed = TRUE)]
    if (length(missing)) {
      stop(ids[index], " is missing required review-boundary text: ",
           paste(missing, collapse = "; "))
    }
  }
  list(ids = ids, titles = manifest$title, variant = variant, phrases = unlist(phrases))
}

validate_auxiliary_rendered_html <- function(contract, html) {
  plain <- gsub("<[^>]+>", "", html)
  plain <- gsub("[[:space:]]+", " ", chartr("\u00A0", " ", plain))
  for (value in c(contract$ids, contract$titles, contract$phrases)) {
    if (!grepl(value, plain, fixed = TRUE) && !grepl(value, html, fixed = TRUE)) {
      stop("Rendered auxiliary figures are missing required record or boundary: ", value)
    }
  }
  positions <- vapply(contract$ids, function(id) regexpr(id, html, fixed = TRUE)[[1]], integer(1))
  if (any(positions < 0L) || is.unsorted(positions)) {
    stop("Rendered auxiliary figures do not preserve manifest row order")
  }
}

figure_source_sha256 <- function(path) {
  output <- system2("sha256sum", shQuote(path), stdout = TRUE)
  digest <- strsplit(output[[1]], "[[:space:]]+")[[1]][[1]]
  if (!is.null(attr(output, "status")) || !grepl("^[a-f0-9]{64}$", digest)) {
    stop("Unable to fingerprint report source: ", path)
  }
  digest
}

report_source_fingerprint_html <- function(path, artifact) {
  sprintf('<span data-report-source="%s" data-source-sha256="%s"></span>',
          artifact, figure_source_sha256(path))
}

validate_report_source_fingerprint <- function(source_path, html) {
  artifact <- paste0("Publication/", basename(source_path))
  prefix <- sprintf('<span data-report-source="%s" data-source-sha256="', artifact)
  positions <- gregexpr(prefix, html, fixed = TRUE)[[1]]
  expected <- report_source_fingerprint_html(source_path, artifact)
  if (positions[[1]] < 0L || length(positions) != 1L || !grepl(expected, html, fixed = TRUE)) {
    stop("Rendered report is missing its current source fingerprint")
  }
}
