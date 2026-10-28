# utils/pub_settings.R
# Sourced by AmAni.summarize.R — defines PUB_SETTINGS, the single place where
# the publication-figure knobs live (P0 of the publication-figures section).
#
# Everything downstream (collect_pub_data(), the fig_*() builders and
# save_pub_fig()) reads its defaults from this list, so changing a headline
# scenario or a figure width here changes it everywhere.

PUB_SETTINGS <- list(
  # --- Data scope -------------------------------------------------------------
  output_dir = "output",           # organism results root (project-relative)
  include     = NULL,              # NULL = all discovered organisms; otherwise
                                   # a character vector of output-dir names
  group       = NULL,               # NULL, or a name from ORGANISM_GROUPS
                                   # (e.g. "vertebrates"); unioned with include

  # --- Which model run is "the" result ----------------------------------------
  headline = list(
    model    = "Clustered",         # spatial model shown as the headline result
    scenario = "FreqBoundary"       # deforestation frequency scenario ditto
  ),
  # Sensitivity panel (Fig 3b): everything else gets shown compactly
  sensitivity = list(
    models    = c("Neutral", "Clustered", "Shifted"),
    scenarios = c("FreqUnder", "FreqBoundary", "FreqCalibrated", "FreqOver")
  ),

  # --- IUCN categories ---------------------------------------------------------
  iucn_categories = c("NT/LC", "VU", "EN", "CR/EX"),

  # --- Figure geometry ---------------------------------------------------------
  # Journal sizes in mm. "single" = one column (89 mm), "double" = full width
  # (183 mm). Heights are per-figure and passed to save_pub_fig() directly.
  # Nature: 89 mm / 183 mm columns, max height 247 mm, sans-serif (Arial or
  # Helvetica) text at 5-8 pt, line widths 0.25-1 pt.
  journal           = "Nature",
  fig_width_single  = 89,
  fig_width_double  = 183,
  fig_max_height    = 247,
  units             = "mm",
  dpi               = 600,
  base_font_size    = 7,           # pt, for theme_amani_pub()
  base_family       = "sans",      # Arial/Helvetica via the sans alias

  # --- Output formats -----------------------------------------------------------
  # PDF is always written (cairo_pdf: vector strokes, fonts embedded as
  # subsets) so figures open directly in Inkscape with editable text.
  # PNG is an optional preview copy for quick viewing/slides.
  save_png_preview  = TRUE,

  # --- Colour & map style -------------------------------------------------------
  # Sequential palette for grid-cell richness / model outputs (scico/Crameri).
  richness_palette  = "lajolla",
  # Bright sequential for NON-map plots (charts, panel fills): viridis
  # "rocket" runs dark blue -> purple -> magenta/pink. Perceptually uniform
  # and colourblind-safe, unlike Turbo.
  bright_palette    = "rocket",
  # Basemap: cells on the bottom, white country outlines (Natural Earth
  # 1:110m) drawn on top so cell fill stays legible.
  map_border_col   = "white",
  map_border_lwd   = 0.65,

  # --- Output location ----------------------------------------------------------
  # Figures are written to <output_dir>/<results_dir>/ (both PDF and PNG).
  results_dir = "publication_figures"
)
