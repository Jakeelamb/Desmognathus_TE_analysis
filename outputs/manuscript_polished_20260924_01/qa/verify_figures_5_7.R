source('scripts/common.R')

checks <- character()
verify <- function(ok, description) {
  stopifnot(isTRUE(ok))
  checks <<- c(checks, paste('PASS:', description))
}
new_csv <- function(n, file) read_csv(file.path(RUN, 'figures', paste0('figure_', n), 'analysis', file), show_col_types = FALSE)
verify(isTRUE(all.equal(as.data.frame(read_analysis(5, 'model_summary.csv')),
  as.data.frame(new_csv(5, 'frozen_genome_cellular_models.csv')), tolerance = 0)), 'genome-cellular model table exactly equals frozen source')
verify(isTRUE(all.equal(as.data.frame(read_analysis(5, 'fitted_relationships.csv')),
  as.data.frame(new_csv(5, 'frozen_genome_cellular_predictions.csv')), tolerance = 0)), 'genome-cellular predictions exactly equal frozen source')
verify(isTRUE(all.equal(as.data.frame(read_analysis(10, 'allometry_model_summary.csv')),
  as.data.frame(new_csv(5, 'frozen_allometry_models.csv')), tolerance = 0)), 'allometry model table exactly equals frozen source')
verify(isTRUE(all.equal(as.data.frame(read_analysis(10, 'fitted_relationship.csv')),
  as.data.frame(new_csv(5, 'frozen_allometry_predictions.csv')), tolerance = 0)), 'allometry predictions exactly equal frozen source')
verify(isTRUE(all.equal(as.data.frame(read_analysis(7, 'body_size_regressions.csv')),
  as.data.frame(new_csv(7, 'frozen_body_size_models.csv')), tolerance = 0)), 'body-size model table exactly equals frozen source')
verify(isTRUE(all.equal(as.data.frame(read_analysis(7, 'fitted_relationships.csv')),
  as.data.frame(new_csv(7, 'frozen_body_size_predictions.csv')), tolerance = 0)), 'body-size predictions exactly equal frozen source')
original_forest <- read_analysis(9, 'displayed_coefficients.csv')
new_forest <- new_csv(7, 'frozen_conditional_sensitivity_coefficients.csv')
verify(isTRUE(all.equal(as.data.frame(original_forest), as.data.frame(new_forest[,names(original_forest)]),
  tolerance = 0)), 'all 12 conditional coefficients and intervals exactly equal frozen source')
species <- sub('^D\\. ', '', new_csv(5, 'genome_cellular_species.csv')$species)
fractions <- new_csv(5, 'paired_nuclear_fraction_species.csv')
verify(fractions$species[which.max(fractions$y)] == 'anicetus' &&
  fractions$species[which.min(fractions$y)] == 'wrighti',
  'paired nuclear-fraction rows match the other figures: anicetus at top, wrighti at bottom')
for (n in c(5, 7)) {
  folder <- file.path(RUN, 'figures', paste0('figure_', n))
  pdfs <- list.files(folder, pattern = '[.]pdf$', recursive = TRUE, full.names = TRUE)
  verify(length(pdfs) == 5, paste('Figure', n, 'has one collage and four isolated PDF panels'))
  for (pdf in pdfs) {
    fonts <- system2('pdffonts', pdf, stdout = TRUE)
    font_rows <- fonts[-c(1,2)]
    verify(all(grepl('Arial', font_rows)) && all(grepl('yes +yes +yes', font_rows)),
      paste(basename(pdf), 'uses only embedded Arial font subsets'))
  }
  for (tag in letters[1:3]) {
    pdf <- file.path(folder, 'panels', paste0('Figure_', n, tag, '_isolated.pdf'))
    txt <- paste(system2('pdftotext', c('-layout', pdf, '-'), stdout = TRUE), collapse = '\n')
    verify(all(vapply(species, grepl, logical(1), x = txt, fixed = TRUE)),
      paste(basename(pdf), 'contains all 24 species labels'))
  }
  tiffs <- list.files(folder, pattern = '[.]tif$', recursive = TRUE, full.names = TRUE)
  for (tiff in tiffs) {
    info <- system2('identify', c('-format', shQuote('%[colorspace] %x %y'), tiff), stdout = TRUE)
    fields <- strsplit(info, ' +')[[1]]
    verify(length(fields) == 3 && fields[1] == 'CMYK' && as.numeric(fields[2]) >= 300 &&
      fields[2] == fields[3], paste(basename(tiff), 'is CMYK at', fields[2], 'dpi'))
  }
}
write_lines(checks, file.path(RUN, 'qa', 'figures_5_7_export_verification.txt'))
cat(paste(checks, collapse = '\n'), '\n')
