# utils/pub_theme.R
# Sourced by AmAni.summarize.R — relies on the packages (tidyverse, sf,
# patchwork, ggrepel) loaded there and on PUB_SETTINGS from
# utils/pub_settings.R (for save_pub_fig() defaults and base font size).
# Also needs scico (richness palette) and rnaturalearth (map borders).

# Per-session cache for the (slow to fetch) country border geometry
.amani_pub_cache <- new.env(parent = emptyenv())

#
# Shared publication style (P2 of the publication-figures section):
#   theme_amani_pub()   base theme for all publication figures
#   AMANI_* colour sets fixed palettes (taxon group, filter stage, IUCN)
#   organism_palette()  distinct per-organism colours, stable across figures
#   organism_levels()  canonical organism ordering (group, then gamma size)
#   save_pub_fig()     writes a figure as PDF (cairo) + PNG
#   ORGANISM_PHYLOPIC / PHYLOPIC_UUIDS  PhyloPic silhouettes for the P3b
#                       retention figure (rphylopic needed only at plot time)

# --- Base theme ----------------------------------------------------------------

# Neutral theme colours (caption text, facet strips, panel grid)
AMANI_THEME_COLS <- list(
  caption   = "grey40",
  strip_bg  = "grey93",
  grid      = "grey88"
)

# Base font for all publication figures (embedded by cairo_pdf on save).
# PUB_SETTINGS$base_family is overridden here so the family is set once,
# centrally, regardless of the settings profile.
AMANI_FONT <- "Arial"

theme_amani_pub <- function(base_size = PUB_SETTINGS$base_font_size,
                            base_family = AMANI_FONT) {
  theme_minimal(base_size = base_size, base_family = base_family) +
    theme(
      # Compact everything: figures are small on the page
      plot.title       = element_text(face = "bold", size = base_size + 1),
      plot.subtitle    = element_text(size = base_size - 0.5,
                                      margin = margin(b = 4)),
      plot.caption     = element_text(size = base_size - 1.5,
                                      colour = AMANI_THEME_COLS$caption),
      strip.text       = element_text(face = "bold", size = base_size),
      strip.background = element_rect(fill = AMANI_THEME_COLS$strip_bg, colour = NA),
      axis.title       = element_text(size = base_size + 1),
      axis.text        = element_text(size = base_size, colour = "black"),
      axis.ticks       = element_line(linewidth = 0.25),
      axis.ticks.length = unit(1, "mm"),
      panel.grid.major = element_line(linewidth = 0.25,
                                      colour = AMANI_THEME_COLS$grid),
      legend.title     = element_text(size = base_size),
      legend.text      = element_text(size = base_size),
      legend.key.size  = unit(2.5, "mm"),
      legend.box.margin = margin(0, 0, 0, 0),
      panel.grid.minor = element_blank()
    )
}

# --- Fixed palettes ---------------------------------------------------------------

# Coarse taxon groups (colour carries "vertebrate vs invertebrate")
AMANI_TAXON_COLS <- c(
  Vertebrates   = "#0072B2",
  Invertebrates = "#E69F00"
)

scale_fill_taxon <- function(..., drop = FALSE) {
  scale_fill_manual(values = AMANI_TAXON_COLS, drop = drop, ...)
}
scale_color_taxon <- function(..., drop = FALSE) {
  scale_color_manual(values = AMANI_TAXON_COLS, drop = drop, ...)
}

# Before/after filtering stages (matches build_dataset_summary_figure())
AMANI_STAGE_COLS <- c(Before = "grey70", After = "#2c7fb8")

scale_fill_stage <- function(..., drop = FALSE) {
  scale_fill_manual(values = AMANI_STAGE_COLS, drop = drop, ...)
}

# IUCN threat categories (official IUCN red-list colours, as already used
# in summarize_all())
AMANI_IUCN_COLS <- c(
  "NT/LC" = "#60C659",
  "VU"    = "#F9E814",
  "EN"    = "#FC7F3F",
  "CR/EX" = "#D81E05"
)

scale_fill_iucn <- function(..., drop = FALSE) {
  scale_fill_manual(values = AMANI_IUCN_COLS, drop = drop, ...)
}

# --- Per-organism colours -----------------------------------------------------------

# Distinct colours for per-organism lines/points. Hues are spread over the
# full wheel (Okabe-Ita-like saturation) and assigned in the canonical
# organism order, so a given organism keeps its colour in every figure.
organism_palette <- function(levels) {
  if (length(levels) <= 1) return(setNames("#0072B2", levels))
  setNames(scales::hue_pal(h = c(0, 360) + 15, l = 62, c = 80)(length(levels)),
           levels)
}

scale_color_organism <- function(levels, ...) {
  scale_color_manual(values = organism_palette(levels), ...)
}
scale_fill_organism <- function(levels, ...) {
  scale_fill_manual(values = organism_palette(levels), ...)
}

# --- PhyloPic silhouettes (P3b retention figure) ---------------------------------

# Organism -> silhouette taxon. Duplicate organism variants (e.g.
# Ants_combined / Ants_curated) share their taxon's silhouette. Organisms
# absent from this mapping (allCombined, Miridae) are dropped from the figure.
# rphylopic is loaded only by the figure code in AmAni.summarize.R, not here.
ORGANISM_PHYLOPIC <- c(
  cockroaches         = "Blattodea",
  Cockroaches.Curated = "Blattodea",      # Blaberus giganteus
  Antbirds            = "Thamnophilidae",
  Ants_combined       = "Formicidae",
  Ants_curated        = "Formicidae",
  Bees                = "Apidae",
  Bees_combined       = "Apidae",
  Beetles             = "Coleoptera",
  `Butterfl-curated`  = "Papilionoidea",
  Butterflies         = "Papilionoidea",
  DungBeetles         = "Scarabaeidae",
  Finches             = "Fringillidae",
  Flies               = "Diptera",
  Frogs               = "Anura",
  Hemiptera           = "Hemiptera",
  Hymenoptera         = "Hymenoptera",
  Lepidoptera         = "Lepidoptera",     # Actias luna
  Mammals             = "Mammalia",       # Cebus
  Mantodea            = "Mantodea",
  Odonata             = "Odonata",
  `Odonata-combined`  = "Odonata",
  `Odonata-curated`   = "Odonata",
  Squamata            = "Squamata",       # Anolis
  Squamates           = "Squamata",
  Termites            = "Isoptera",
  Tyrannidae          = "Tyrannidae"
)

# Silhouette image UUIDs, picked by hand from candidate grids on phylopic.org
# and hard-coded so the figure renders without network access. Credit the
# contributors in the figure caption:
#   rphylopic::get_attribution(uuid = unname(PHYLOPIC_UUIDS))
PHYLOPIC_UUIDS <- c(
  Anura         = "9e1615e2-4fc8-4967-8a3a-baff830fef51",
  Apidae        = "c8c12aaa-8f2b-4aed-b1b4-d2e4ae1988c3",
  Blattodea     = "5103e59a-304b-454f-b2d5-ded3d64088d2",
  Coleoptera    = "07cbdb20-25a4-4caa-9a2a-f35dd39d2939",
  Diptera       = "60075264-6a3c-43dd-aaa9-d6a52b88556a",
  Formicidae    = "45ccedd3-d5cb-42e8-b158-7e027ab1ff22",
  Fringillidae  = "c0535f8a-a119-44e3-8343-e55a3bc84e58",
  Hemiptera     = "408b5725-338f-48b6-b57f-2317797752cb",
  Hymenoptera   = "23d5cab5-7396-41d4-85d8-747a4524feca",
  Mantodea      = "b673278f-6d7c-4721-8a6e-dd727c1e3bf2",
  Odonata       = "00415ca0-dbca-4312-b871-840fe1c7031b",
  Papilionoidea = "04af22bd-6772-445e-a66f-1009173eea2d",
  Scarabaeidae  = "7cb79685-9ff4-4261-9cc8-0942975a87cb",
  Thamnophilidae = "420c0aff-f49e-472e-99d6-d9dee5b17889",
  Tyrannidae    = "baed3744-96f9-41a4-aa7d-a304ddfbb688",
  Isoptera      = "c44be6d9-6398-4753-84db-9609760fd37c",
  Lepidoptera   = "8229756b-82c3-4a9f-a1c6-e88f958e623e",
  Mammalia      = "156b515d-f25c-4497-b15b-5afb832cc70c",
  Squamata      = "0f1f41e3-8f10-489f-9cbb-7ce8f8a6653e"
)

# --- Canonical organism ordering -------------------------------------------------

# Invertebrates first (larger, more numerous group), then vertebrates; within
# each group, descending full-model gamma richness. Every figure facets or
# orders rows/axes by this so panels line up across figures.
organism_levels <- function(pub) {
  pub$summary |>
    dplyr::left_join(pub$organisms, by = c("organism_dir", "organism")) |>
    dplyr::arrange(
      factor(taxon_group, levels = c("Invertebrates", "Vertebrates")),
      dplyr::desc(gamma_best)) |>
    dplyr::pull(organism)
}

# --- Sequential (richness / map) scales ------------------------------------------

# Grid-cell richness and continuous model outputs. Lajolla (Crameri, via
# scico): perceptually uniform, colourblind-safe, bright high end that reads
# well over maps. Reads from PUB_SETTINGS$richness_palette so the journal
# profile can switch it centrally.

scale_fill_richness <- function(..., palette = PUB_SETTINGS$richness_palette) {
  scico::scale_fill_scico(palette = palette, ...)
}
scale_color_richness <- function(..., palette = PUB_SETTINGS$richness_palette) {
  scico::scale_color_scico(palette = palette, ...)
}

# Bright purple/pink sequential for non-map plots (bar charts, ridge fills,
# panel backgrounds of dot plots, ...). Defaults to viridis "rocket".
scale_fill_bright <- function(..., option = PUB_SETTINGS$bright_palette) {
  ggplot2::scale_fill_viridis_c(option = option, ...)
}
scale_color_bright <- function(..., option = PUB_SETTINGS$bright_palette) {
  ggplot2::scale_color_viridis_c(option = option, ...)
}

# --- Map style --------------------------------------------------------------------

# Borders-over-cells basemap: grid cells are drawn first (fill, no outline),
# then country outlines on top in PUB_SETTINGS$map_border_col so fills stay
# legible while country context reads through. Natural Earth 1:110m, cached
# per session.

amani_country_borders <- function() {
  if (!exists("AMANI_BORDERS", envir = .amani_pub_cache)) {
    assign("AMANI_BORDERS",
           rnaturalearth::ne_countries(scale = 110, returnclass = "sf"),
           envir = .amani_pub_cache)
  }
  get("AMANI_BORDERS", envir = .amani_pub_cache)
}

# All map-layer colours and line widths used by the cell-filter / retained-cell
# maps. Defined once here so every map figure stays consistent; the figure
# code in AmAni.summarize.R only references these.
AMANI_MAP_STYLE <- list(
  land_fill         = "grey95",  # land polygons; water reads as white background
  amazon_fill       = "#c7e9c0",  # Amazon shapefile, slightly darker than land
  border_col        = "grey70",  # country border lines
  border_lwd        = 0.65,      # country border line width
  cell_outline_col  = "black",   # grid-cell outline
  cell_outline_lwd  = 0.1,       # grid-cell outline width
  pct_outline_lwd   = 0.8,       # >= 95th-percentile cell highlight width
  alpha_floor       = 0.2,       # min transparency for alpha-encoded cell fills
  richness_begin    = 0.3,       # low end of the richness palette on maps (its
                                 # bottom is near-white and would read as empty)
  shadow_shift      = 0.08,      # drop-shadow offset, as a fraction of one cell
  shadow_fill       = "grey20",  # drop-shadow fill under highlighted cells
  shadow_alpha      = 0.25       # drop-shadow opacity
)

# Amazon basin (sensu latissimo) outline shapefile
AMAZON_SHP <- "0-CuratingPlantData/data/amazon_shapefile/amazon_sensulatissimo_gmm_v1.shp"

# Load and prepare the Amazon outline for a map: fix invalid geometries,
# reproject to the cells' CRS and crop to the map extent.
amazon_outline <- function(crs, bbox) {
  st_read(AMAZON_SHP, quiet = TRUE) |>
    st_make_valid() |>
    st_transform(crs) |>
    st_crop(st_as_sfc(bbox))
}

# geom to overlay AFTER the cell layer, e.g.:
#   ggplot(cells) + geom_sf(aes(fill = richness), color = NA) +
#     amani_borders() + theme_amani_map()
amani_borders <- function(...) {
  geom_sf(
    data = amani_country_borders(),
    fill = NA,
    colour = PUB_SETTINGS$map_border_col,
    linewidth = PUB_SETTINGS$map_border_lwd,
    ...
  )
}

# Map-flavoured theme: theme_amani_pub() plus map conventions (no gridlines,
# no axes, fixed aspect)
theme_amani_map <- function(...) {
  theme_amani_pub(...) +
    theme(
      panel.grid = element_blank(),
      axis.title = element_blank(),
      axis.text  = element_blank(),
      axis.ticks = element_blank()
    )
}

# --- Saving ---------------------------------------------------------------------------

# Save a figure as vector PDF (cairo) plus an optional 600-dpi PNG preview.
# cairo_pdf keeps all geometry as vector paths and embeds fonts as subsets,
# so the PDF opens directly in Inkscape: strokes remain strokes, and text
# stays editable (Poppler maps the embedded subsets back to text on import).
# Writes to <output_dir>/<PUB_SETTINGS$results_dir>/<filename>.<pdf|png>.
save_pub_fig <- function(plot, filename,
                         dir = file.path(PUB_SETTINGS$output_dir,
                                         PUB_SETTINGS$results_dir),
                         width = PUB_SETTINGS$fig_width_double,
                         height = NULL, units = PUB_SETTINGS$units,
                         dpi = PUB_SETTINGS$dpi,
                         png = PUB_SETTINGS$save_png_preview) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  height <- height %||% (width * 0.6)  # default aspect if not given

  # Journal guardrail (e.g. Nature caps figure height at 247 mm)
  max_h <- PUB_SETTINGS$fig_max_height
  if (units == "mm" && !is.null(max_h) && height > max_h) {
    warning("Figure height ", round(height), " mm exceeds ", max_h,
            " mm (", PUB_SETTINGS$journal, " limit); consider a shorter layout.")
  }

  ggsave(file.path(dir, paste0(filename, ".pdf")), plot = plot,
         width = width, height = height, units = units,
         device = cairo_pdf)
  if (png) {
    ggsave(file.path(dir, paste0(filename, ".png")), plot = plot,
           width = width, height = height, units = units,
           dpi = dpi, bg = "white")
  }

  message("  Wrote: ", file.path(dir, paste0(filename, ".pdf")),
          if (png) paste0(" + .png"),
          " (", width, " x ", round(height), " ", units, ")")
  invisible(file.path(dir, paste0(filename, ".pdf")))
}
