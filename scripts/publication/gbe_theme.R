# Shared Genome Biology and Evolution figure style for this project.

suppressPackageStartupMessages({
  library(ggplot2)
  library(scales)
})

GBE_FONT <- "Arial"
GBE_FALLBACK_FONT <- "Nimbus Sans"
GBE_BASE_SIZE <- 8
GBE_SINGLE_COLUMN_MM <- 89
GBE_DOUBLE_COLUMN_MM <- 185

GBE_COLORS <- c(
  blue = "#0072B2",
  orange = "#E69F00",
  green = "#009E73",
  vermillion = "#D55E00",
  sky = "#56B4E9",
  purple = "#CC79A7",
  yellow = "#F0E442",
  charcoal = "#30343B",
  gray = "#7A7F84",
  light_gray = "#D9DEE2"
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
      axis.line = element_line(linewidth = 0.35, colour = GBE_COLORS[["charcoal"]]),
      axis.ticks = element_line(linewidth = 0.3, colour = GBE_COLORS[["charcoal"]]),
      axis.ticks.length = grid::unit(1.5, "mm"),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_line(linewidth = 0.2, colour = "#E7EAEC"),
      strip.background = element_blank(),
      strip.text = element_text(size = base_size, face = "bold", hjust = 0),
      plot.title = element_blank(),
      plot.subtitle = element_blank(),
      plot.caption = element_text(size = base_size - 1, hjust = 0, colour = GBE_COLORS[["gray"]]),
      plot.tag = element_text(size = base_size + 1, face = "bold"),
      plot.tag.position = c(0, 1),
      legend.position = "bottom",
      legend.title = element_text(size = base_size - 0.5, face = "bold"),
      legend.text = element_text(size = base_size - 0.5),
      legend.key.height = grid::unit(3.5, "mm"),
      legend.key.width = grid::unit(5, "mm"),
      plot.margin = margin(4, 5, 4, 5, unit = "pt")
    )
}

scale_color_gbe <- function(..., values = unname(GBE_COLORS[c("blue", "orange", "green", "vermillion", "purple", "sky")])) {
  scale_color_manual(values = values, ...)
}

scale_fill_gbe <- function(..., values = unname(GBE_COLORS[c("blue", "orange", "green", "vermillion", "purple", "sky")])) {
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
