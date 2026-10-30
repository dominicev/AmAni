# utils/pub_0_settings.R — Publication figures, step 0: settings
#
# Defines PUB_SETTINGS, the single place where the publication-figure knobs
# live. Everything downstream (collect_pub_data(), the fig_*() builders and
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

  # --- Hymenoptera estimate -----------------------------------------------------
  # Fig 2 shows Hymenoptera as ants (Ants_combined), bees (Bees_combined) and
  # "Other Hymenoptera", where
  # Other = (ants + bees) * parasitoid / non-parasitic species. The ratio is
  # the mean over the three paired estimates below (ratios ~13.8-14.2, mean 14.05).
  # Source: Forbes AA, Bagley RK, Beer MA, Hippee AC, Widmayer HA (2018)
  # Quantifying the unquantifiable: why Hymenoptera, not Coleoptera, is the
  # most speciose animal order. BMC Ecology 18:21. doi:10.1186/s12898-018-0176-x
  # Estimated species richness; element i of each vector is the same estimate.
  # Set to NULL to use the modelled Hymenoptera instead.
  hymenoptera_ratios = list(
    nonparasitic = c(79980, 58900, 29450),
    parasitoid   = c(1107487, 833590, 416795),
    source       = "Forbes et al. 2018, BMC Ecology 18:21"
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
