# Inspection exports reuse the collage's ggplots; no assay values are recomputed.
figure1_review_panels <- function(panels, tips) {
  if (!identical(names(panels), letters[1:5]) || nrow(tips) != 21L || anyDuplicated(tips$y)) {
    stop("Figure 1 inspection requires panels a-e and the 21 original species rows")
  }
  result <- lapply(names(panels), function(tag) {
    plot <- panels[[tag]]
    if (tag != "a") {
      plot <- suppressMessages(plot + scale_y_continuous(
        limits = c(0.5, 21.5), breaks = tips$y, labels = tips$display_species,
        expand = c(0, 0)
      )) + theme(axis.text.y = element_text(
        family = resolve_gbe_font(), face = "italic", size = GBE_MIN_TEXT_PT
      ))
    }
    plot + labs(tag = tag) + theme(
      legend.position = "bottom", legend.box = "vertical",
      plot.tag.position = "topleft", plot.margin = margin(8, 12, 8, 8, unit = "pt")
    )
  })
  names(result) <- names(panels)
  result
}

export_figure1_panels <- function(panels, tips, output_dir, collage) {
  plots <- figure1_review_panels(panels, tips)
  panel_dir <- file.path(output_dir, "panels")
  titles <- c(
    "Phylogeny and reproductive strategy", "Estimated genome size",
    "Broad TE composition", "TE superfamilies", "Repeat divergence profiles"
  )
  basenames <- paste0("Figure_1", letters[1:5], "_", c(
    "phylogeny_and_strategy", "conditional_genome_size", "broad_TE_composition",
    "TE_superfamilies", "divergence_profiles"
  ))
  legends <- c(
    "The original topology-only tree and reproductive-strategy symbols for the 21 species shared by the genomic and cellular-trait panels. Species labels retain the collage order. Branch lengths and a time scale are omitted because the source tree's publication, calibration, and time provenance remain unresolved.",
    "Estimated genome sizes (pg/1C*) derived from IOD and anchored to the D. fuscus assembly, with 95% measurement-bootstrap intervals. One common circle symbol is used and each estimate is labelled to one decimal place. The 8-27 pg/1C axis retains every interval. These are not independent species C-values; intervals exclude anchor uncertainty and remain conditional on observed image and animal support. Support cautions are retained in S24 rather than encoded by changing point symbols.",
    "Broad TE-order shares of classified dnaPipeTE aligned bases for the same 21 species. LTR, LINE, TIR, and Other retain the collage's category colours and common 0-100% scale. These are composition fractions within the assay, not genome occupancy.",
    "The same nine abundant dnaPipeTE superfamilies in the same numbered positions and species order. Every row's miniature bars share the original 0-70% classified-aligned-base scale, from the bottom to the top of its frame. The numbered legend identifies every family; the selected subset is not renormalized.",
    "The original 63 LTR, LINE, and TIR RepeatMasker profiles, with every row labelled by species. Each profile retains 0-50% sequence divergence in 1-percentage-point bins and 0-3% of all RepeatMasker-aligned hit bp for that species, from the bottom to the top of its frame. Absent order/bin combinations are zero-filled without subset renormalization. The denominator differs from dnaPipeTE composition. Divergence is not insertion age, and hit-bp fractions are not genome occupancy."
  )
  alt_text <- c(
    "Topology-only phylogeny with 21 italic species labels and reproductive-strategy symbols.",
    "Twenty-one labelled species with uniform circular genome-size points, one-decimal estimate labels, and measurement intervals on an 8-27 pg/1C axis.",
    "Twenty-one labelled horizontal stacked bars show LTR, LINE, TIR, and Other composition shares.",
    "Twenty-one labelled rows each contain nine coloured superfamily bars on a common zero-to-70-percent scale.",
    "Twenty-one labelled rows each contain LTR, LINE, and TIR divergence profiles with shared divergence and abundance scales."
  )
  records <- lapply(seq_along(plots), function(i) {
    save_gbe_figure(plots[[i]], file.path(panel_dir, basenames[[i]]),
                    width_mm = GBE_DOUBLE_COLUMN_MM, height_mm = 190, dpi = 300)
    tibble(
      parent_figure = "Figure_1", panel_id = names(plots)[[i]], title = titles[[i]],
      width_mm = GBE_DOUBLE_COLUMN_MM, height_mm = 190,
      final_font_family = resolve_gbe_font(), minimum_target_text_pt = GBE_MIN_TEXT_PT,
      legend = legends[[i]], alt_text = alt_text[[i]],
      pdf = paste0(basenames[[i]], ".pdf"), png = paste0(basenames[[i]], ".png"),
      tiff = paste0(basenames[[i]], ".tif")
    )
  })
  manifest <- bind_rows(records)
  write_csv(manifest, file.path(panel_dir, "PANEL_MANIFEST.csv"))
  write_lines(c(
    "# Figure 1 isolated panels", "",
    "These inspection exports reuse the collage's plot data, scales, colours, and species order.",
    "Species labels and standalone panel tags are restored; each panel is 180 x 190 mm.", "",
    unlist(lapply(seq_len(nrow(manifest)), function(i) c(
      paste0("## ", manifest$panel_id[[i]], " - ", manifest$title[[i]]), "",
      manifest$legend[[i]], "", paste0("**Alt text:** ", manifest$alt_text[[i]]), ""
    )))
  ), file.path(panel_dir, "PANEL_LEGENDS_AND_ALT_TEXT.md"))

  escape_html <- function(x) {
    x <- gsub("&", "&amp;", x, fixed = TRUE)
    x <- gsub("<", "&lt;", x, fixed = TRUE)
    x <- gsub(">", "&gt;", x, fixed = TRUE)
    gsub('"', "&quot;", x, fixed = TRUE)
  }
  card <- function(id, title, png, pdf, tiff, description, alt) {
    sprintf(paste0(
      '<section id="%s"><h2>%s</h2><p class="downloads">',
      '<a href="%s">PDF</a> · <a href="%s">Full-size PNG</a> · ',
      '<a href="%s">CMYK TIFF</a></p>',
      '<a href="%s"><img src="%s" alt="%s" loading="lazy"></a>',
      '<p class="caption">%s</p></section>'
    ), id, escape_html(title), pdf, png, tiff, png, png, escape_html(alt), escape_html(description))
  }
  cards <- c(
    card("collage", "Full collage", collage$png, collage$pdf, collage$tiff,
         collage$legend, collage$alt_text),
    vapply(seq_len(nrow(manifest)), function(i) card(
      paste0("panel-", manifest$panel_id[[i]]),
      paste0(manifest$panel_id[[i]], " - ", manifest$title[[i]]),
      file.path("panels", manifest$png[[i]]), file.path("panels", manifest$pdf[[i]]),
      file.path("panels", manifest$tiff[[i]]), manifest$legend[[i]], manifest$alt_text[[i]]
    ), character(1))
  )
  write_lines(c(
    '<!doctype html><html lang="en"><head><meta charset="utf-8">',
    '<meta name="viewport" content="width=device-width, initial-scale=1">',
    '<title>Figure 1 - collage and isolated panels</title>',
    '<style>body{margin:0;background:#f3f4f5;color:#222;font:16px/1.5 Arial,sans-serif}',
    'header,main{max-width:1100px;margin:auto;padding:24px}h1{font-size:26px;margin:0 0 8px}',
    'nav{display:flex;flex-wrap:wrap;gap:8px 20px}a{color:#075b8b}section{background:white;',
    'padding:24px;margin:0 0 28px;scroll-margin-top:16px}h2{font-size:21px;margin:0}',
    'img{display:block;width:100%;height:auto;margin:16px auto;max-width:1000px}',
    '.caption{max-width:950px}.downloads{margin:6px 0}a:focus-visible{outline:3px solid #0072b2}',
    '@media(max-width:600px){header,main,section{padding:14px}}',
    '@media print{body{background:white}nav,.downloads{display:none}section{break-before:page}}',
    '</style></head><body><header><h1>Figure 1: collage and isolated panels</h1>',
    '<p>The same 21 species, measurements, colours, and scales. Select an image for its full-size PNG.</p>',
    '<nav aria-label="Figure panels"><a href="#collage">Full collage</a>',
    paste(vapply(seq_len(nrow(manifest)), function(i) sprintf(
      '<a href="#panel-%s">%s - %s</a>', manifest$panel_id[[i]],
      manifest$panel_id[[i]], escape_html(manifest$title[[i]])
    ), character(1)), collapse = "\n"),
    '</nav></header><main>', cards, '</main></body></html>'
  ), file.path(output_dir, "index.html"))
  invisible(manifest)
}
