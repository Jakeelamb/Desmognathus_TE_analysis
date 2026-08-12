# Shared Genome Biology and Evolution figure style for this project. The panel
# grammar and palette follow the supplied salamander-genome reference, while the
# journal-facing sans-serif typography and redundant shapes/linetypes remain.

suppressPackageStartupMessages({
  library(ggplot2)
  library(scales)
})

GBE_FONT <- "Arial"
GBE_FALLBACK_FONT <- "Nimbus Sans"
GBE_BASE_SIZE <- 8
GBE_MIN_TEXT_PT <- 7
GBE_LABEL_SIZE_MM <- GBE_MIN_TEXT_PT / ggplot2::.pt
GBE_SINGLE_COLUMN_MM <- 89
GBE_DOUBLE_COLUMN_MM <- 185

GBE_COLORS <- c(
  blue = "#00BFC4",
  orange = "#F8766D",
  green = "#35B779",
  vermillion = "#C44E52",
  sky = "#31688E",
  purple = "#440154",
  yellow = "#FDE725",
  charcoal = "#1A1A1A",
  gray = "#767676",
  light_gray = "#D6D6D6"
)

resolve_gbe_font <- function() {
  installed <- tryCatch(systemfonts::system_fonts()$family, error = function(e) character())
  if (GBE_FONT %in% installed) return(GBE_FONT)
  if (GBE_FALLBACK_FONT %in% installed) return(GBE_FALLBACK_FONT)
  if ("Helvetica" %in% installed) return("Helvetica")
  if ("Liberation Sans" %in% installed) return("Liberation Sans")
  "sans"
}

theme_gbe <- function(base_size = GBE_BASE_SIZE, base_family = resolve_gbe_font()) {
  theme_classic(base_size = base_size, base_family = base_family) %+replace%
    theme(
      text = element_text(colour = GBE_COLORS[["charcoal"]]),
      axis.title = element_text(size = base_size, face = "plain"),
      axis.text = element_text(size = base_size - 0.5, colour = GBE_COLORS[["charcoal"]]),
      axis.line = element_blank(),
      axis.ticks = element_line(linewidth = 0.3, colour = GBE_COLORS[["charcoal"]]),
      axis.ticks.length = grid::unit(1.5, "mm"),
      panel.background = element_rect(fill = "white", colour = NA),
      panel.border = element_rect(fill = NA, linewidth = 0.4, colour = GBE_COLORS[["charcoal"]]),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.spacing = grid::unit(5, "pt"),
      strip.background = element_rect(fill = "white", linewidth = 0.35, colour = GBE_COLORS[["charcoal"]]),
      strip.text = element_text(size = base_size, face = "bold", hjust = 0.5),
      plot.title = element_blank(),
      plot.subtitle = element_blank(),
      plot.caption = element_text(size = base_size - 1, hjust = 0, colour = GBE_COLORS[["gray"]]),
      plot.tag = element_text(size = base_size + 1, face = "bold"),
      plot.tag.position = c(0, 1),
      legend.position = "bottom",
      legend.justification = "left",
      legend.box.just = "left",
      legend.title = element_text(size = base_size - 0.5, face = "plain"),
      legend.text = element_text(size = base_size - 0.5),
      legend.key.height = grid::unit(3, "mm"),
      legend.key.width = grid::unit(4.5, "mm"),
      legend.spacing.x = grid::unit(2, "pt"),
      legend.margin = margin(1, 0, 0, 0, unit = "pt"),
      plot.background = element_rect(fill = "white", colour = NA),
      plot.margin = margin(3, 4, 3, 4, unit = "pt")
    )
}

scale_color_gbe <- function(..., values = unname(GBE_COLORS[c("purple", "sky", "blue", "green", "yellow", "orange")])) {
  scale_color_manual(values = values, ...)
}

scale_fill_gbe <- function(..., values = unname(GBE_COLORS[c("purple", "sky", "blue", "green", "yellow", "orange")])) {
  scale_fill_manual(values = values, ...)
}

save_gbe_figure <- function(plot, basename, width_mm = GBE_DOUBLE_COLUMN_MM, height_mm, dpi = 300) {
  dir.create(dirname(basename), recursive = TRUE, showWarnings = FALSE)
  width_in <- width_mm / 25.4
  height_in <- height_mm / 25.4
  family <- resolve_gbe_font()

  rgb_pdf <- paste0(basename, ".rgb.pdf")
  final_pdf <- paste0(basename, ".pdf")
  grDevices::cairo_pdf(
    rgb_pdf,
    width = width_in,
    height = height_in,
    family = family,
    onefile = TRUE
  )
  print(plot)
  grDevices::dev.off()

  ghostscript_args <- c(
    "-dSAFER", "-dBATCH", "-dNOPAUSE", "-sDEVICE=pdfwrite",
    "-dPDFSETTINGS=/prepress", "-dEmbedAllFonts=true",
    "-sProcessColorModel=DeviceCMYK", "-sColorConversionStrategy=CMYK",
    paste0("-sOutputFile=", final_pdf), rgb_pdf
  )
  status <- system2("gs", ghostscript_args, stdout = FALSE, stderr = FALSE)
  if (status != 0L || !file.exists(final_pdf)) stop("Ghostscript CMYK PDF conversion failed")
  unlink(rgb_pdf)

  png_path <- paste0(basename, ".png")
  ghostscript_png_args <- c(
    "-dSAFER", "-dBATCH", "-dNOPAUSE", "-sDEVICE=png16m",
    "-dTextAlphaBits=4", "-dGraphicsAlphaBits=4", paste0("-r", dpi),
    paste0("-sOutputFile=", png_path), final_pdf
  )
  status <- system2("gs", ghostscript_png_args, stdout = FALSE, stderr = FALSE)
  if (status != 0L || !file.exists(png_path)) stop("Ghostscript PNG rendering failed")

  tiff_path <- paste0(basename, ".tif")
  status <- system2(
    "magick",
    c(
      png_path, "-colorspace", "CMYK", "-units", "PixelsPerInch",
      "-density", as.character(dpi), "-compress", "LZW", tiff_path
    ),
    stdout = FALSE,
    stderr = FALSE
  )
  if (status != 0L || !file.exists(tiff_path)) stop("ImageMagick CMYK TIFF conversion failed")

  invisible(c(final_pdf, png_path, tiff_path))
}
