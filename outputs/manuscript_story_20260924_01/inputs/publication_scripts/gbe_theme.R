# Shared Genome Biology and Evolution figure style for this project.
# Arial, restrained colour, and minimal framing are shared across figures.
# The 88/180-mm widths and lower-case panel labels are project conventions;
# GBE requires labels of at least 6 pt, and our working minimum is 7 pt.

suppressPackageStartupMessages({
  library(ggplot2)
  library(scales)
})

GBE_FONT <- "Arial"
GBE_BASE_SIZE <- 8
GBE_MIN_TEXT_PT <- 7
GBE_LABEL_SIZE_MM <- GBE_MIN_TEXT_PT / ggplot2::.pt
GBE_LINE_MM <- 1 / ggplot2::.pt
GBE_SINGLE_COLUMN_MM <- 88
GBE_DOUBLE_COLUMN_MM <- 180

GBE_COLORS <- c(
  blue = "#0072B2",
  orange = "#E69F00",
  green = "#009E73",
  vermillion = "#D55E00",
  sky = "#56B4E9",
  purple = "#6A3D9A",
  yellow = "#F0E442",
  charcoal = "#222222",
  gray = "#6F6F6F",
  light_gray = "#D9D9D9"
)

# Named biological mappings remain stable when figures use different subsets.
GBE_TE_ORDER_COLORS <- c(
  LTR = "#6A3D9A", LINE = "#0072B2", TIR = "#009E73", Other = "#D9D9D9"
)
GBE_TE_SUPERFAMILY_COLORS <- c(
  Gypsy = "#440154", hAT = "#482878", L2 = "#3E4989", DIRS = "#31678E",
  L1 = "#26828E", `PIF-Harbinger` = "#1E9C89", Mutator = "#33B67A",
  `Tc1-mariner` = "#68CC5C", CACTA = "#AEDC30"
)

resolve_gbe_font <- function() {
  installed <- systemfonts::system_fonts()$family
  if (!GBE_FONT %in% installed) {
    stop(
      "Arial is required for GBE figures. Install Arial (regular, bold, italic, ",
      "and bold italic), refresh the font cache, and restart R; font substitution ",
      "is not permitted.",
      call. = FALSE
    )
  }
  GBE_FONT
}

theme_gbe <- function(base_size = GBE_BASE_SIZE, base_family = resolve_gbe_font()) {
  theme_classic(base_size = base_size, base_family = base_family) %+replace%
    theme(
      text = element_text(colour = GBE_COLORS[["charcoal"]]),
      axis.title = element_text(size = base_size, face = "plain"),
      axis.text = element_text(
        size = base_size - 0.5,
        colour = GBE_COLORS[["charcoal"]],
        lineheight = 0.95
      ),
      axis.line = element_line(
        linewidth = GBE_LINE_MM,
        colour = GBE_COLORS[["charcoal"]]
      ),
      axis.ticks = element_line(
        linewidth = GBE_LINE_MM,
        colour = GBE_COLORS[["charcoal"]]
      ),
      axis.ticks.length = grid::unit(1, "mm"),
      panel.background = element_rect(fill = "white", colour = NA),
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.spacing = grid::unit(5, "pt"),
      strip.background = element_rect(fill = "#F4F4F4", colour = NA),
      strip.text = element_text(
        size = base_size,
        face = "bold",
        hjust = 0.5,
        margin = margin(2, 2, 2, 2, unit = "pt")
      ),
      plot.title = element_blank(),
      plot.subtitle = element_blank(),
      plot.caption = element_text(size = base_size - 1, hjust = 0, colour = GBE_COLORS[["gray"]]),
      plot.tag = element_text(size = base_size + 1, face = "bold"),
      plot.tag.position = c(0.005, 0.995),
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

theme_gbe_framed_facets <- function() {
  theme(
    axis.line = element_blank(),
    panel.border = element_rect(
      fill = NA,
      linewidth = GBE_LINE_MM,
      colour = GBE_COLORS[["charcoal"]]
    ),
    strip.background = element_rect(
      fill = "#F4F4F4",
      linewidth = GBE_LINE_MM,
      colour = GBE_COLORS[["charcoal"]]
    )
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

  final_pdf <- paste0(basename, ".pdf")
  grDevices::cairo_pdf(
    final_pdf,
    width = width_in,
    height = height_in,
    family = family,
    onefile = TRUE
  )
  print(plot)
  grDevices::dev.off()

  png_path <- paste0(basename, ".png")
  ghostscript_png_args <- c(
    "-dSAFER", "-dBATCH", "-dNOPAUSE", "-sDEVICE=png16m",
    "-dTextAlphaBits=4", "-dGraphicsAlphaBits=4", paste0("-r", dpi),
    paste0("-sOutputFile=", png_path), final_pdf
  )
  status <- system2("gs", ghostscript_png_args, stdout = FALSE, stderr = FALSE)
  if (status != 0L || !file.exists(png_path)) stop("Ghostscript PNG rendering failed")

  tiff_path <- paste0(basename, ".tif")
  ghostscript_tiff_args <- c(
    "-dSAFER", "-dBATCH", "-dNOPAUSE", "-sDEVICE=tiff32nc",
    "-sCompression=lzw", paste0("-r", dpi),
    paste0("-sOutputFile=", tiff_path), final_pdf
  )
  status <- system2("gs", ghostscript_tiff_args, stdout = FALSE, stderr = FALSE)
  if (status != 0L || !file.exists(tiff_path)) stop("Ghostscript CMYK TIFF rendering failed")

  invisible(c(final_pdf, png_path, tiff_path))
}
