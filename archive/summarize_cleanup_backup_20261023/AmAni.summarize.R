# AmAni.summarize.R — Cross-organism summary of extinction model results
#
# Provides two entry-point functions:
#   summarize_all()        — cross-organism summary (table + figures)
#   summarize_organism()   — single-organism compilation (table + figure)
#
# Reads summary CSVs, gamma richness estimates, and saved plot data from each
# organism's output directory and produces:
#   1. A combined summary table (CSV)
#   2. Comparison figures across organisms (PDF)
#   3. Compilation figures regenerated from saved plot data (PDF)
#
# Designed to be sourced at the end of AmAni.run.R or run standalone.
#
# The function definitions live in utils/ (one script per function); sourcing
# this file loads all of them. As elsewhere in the project, the working
# directory is assumed to be the project root.

library(tidyverse)
library(sf)
library(patchwork)
library(ggrepel)

# --- Load function definitions from utils/ -----------------------------------
# Order matters only for readability: metadata first, then the plot/report
# helpers, then the two entry points and the standalone comparison figures.

# Organism metadata
source("utils/organism_label.R")          # organism_label()
source("utils/organism_groups.R")         # ORGANISM_GROUPS

# Dataset before/after filtering helpers
source("utils/scrape_dataset_summary_stats.R")
source("utils/load_dataset_summary_figure_inputs.R")
source("utils/build_dataset_summary_figure.R")

# Cell-filtering report (per organism)
source("utils/build_cell_filter_map_panels.R")
source("utils/build_completeness_correlation_plot.R")
source("utils/write_cell_filter_report_pdf.R")

# Entry points
source("utils/summarize_all.R")
source("utils/summarize_organism.R")

# Standalone comparison figures
source("utils/plot_spatial_extent_overlap.R")
source("utils/compile_workflow_diagram.R")


# =============================================================================
# PUBLICATION FIGURES (prototype — migrate to utils/ once stable)
# =============================================================================

# --- P0. Settings --------------------------------------------------------------
# PUB_SETTINGS (utils/pub_settings.R): which organisms are included, the
# headline model + scenario (Clustered / FreqBoundary), sensitivity grid,
# IUCN categories, figure widths (89 / 183 mm), output location.

source("utils/pub_settings.R")

# --- P1. Data collection (lifted from summarize_all.R parsing) ----DONE------------
# Five collectors, each returning a tidy tibble tagged with both the raw
# directory name (`organism_dir`) and the display label (`organism`):
#   collect_extinction_summary()    summary_all_models.csv, all models x scenarios
#   collect_organism_summary()     dataset_summary + gamma + sampling_extension
#   collect_dataset_distributions() per-cell five-number summaries (before/after)
#   collect_iucn()                 iucn_threat_classification_data.rds $counts
#   collect_pub_data()             wraps the above into one named list
#
# Parsing logic is lifted from summarize_all() so the two stay consistent;
# once stable, summarize_all() should call these instead of its inline code.
# Like the rest of the project, the working directory is assumed to be the
# project root, so output_dir defaults to "output".

source("utils/pub_data_helpers.R")           # .latest_file, .extract_val, .organism_dirs
source("utils/collect_extinction_summary.R")
source("utils/collect_organism_summary.R")
source("utils/collect_dataset_distributions.R")
source("utils/collect_iucn.R")
source("utils/collect_pub_data.R")

# --- P2. Shared style ------------------------------------------------------------
# theme_amani_pub(), fixed palettes (taxon / stage / IUCN), per-organism
# palette, canonical organism ordering, save_pub_fig() (PDF cairo + PNG).

source("utils/pub_theme.R")
source("utils/retained_cell_site_names.R")   # RETAINED_CELL_SITES
source("utils/label_retained_cells.R")       # label_retained_cells()

# --- P3. Figure 1: Sampling and filtering -------------------------------------
# a. Records / cells / species before vs after (dumbbell, log scale)
# b. % retained (records, cells, species)
# c. Per-cell completeness & richness distributions before vs after
# d. (optional) well-sampled vs super-diverse cells / spatial footprint


collect_retained_cells <- function(output_dir = PUB_SETTINGS$output_dir) {
  files <- Sys.glob(file.path(output_dir, "*", "animal_model", "plot_data",
                              "plotdata_cell_filter_report.rds"))
  names(files) <- basename(dirname(dirname(dirname(files))))
  purrr::imap(files, \(f, org) {
    readRDS(f)$delta |>
      filter(retained) |>
      mutate(organism = org) |>
      select(organism, grid_id, plant_completeness,
             animal_richness, animal_completeness)
  }) |>
    bind_rows()
}

# Per-cell Chao1 completeness (and completeness-weighted richness) from the
# RAW filtered plant counts. The plant_completeness stored in the cell-filter
# reports is computed from the FILLED plant data (when fill_plants = TRUE),
# whose inflated record counts leave no singletons, so it is exactly 1 in
# every cell (see the note at the plant_completeness calculation in
# 02_animalAnalysis.R). Figures that encode plant completeness must use this
# raw version instead.
collect_plant_completeness_raw <- function(
    counts_path = file.path("1-FilteringPlantData", "output",
                            "species_counts4_filtered.rds")) {
  readRDS(counts_path) |>
    st_drop_geometry() |>
    group_by(grid_id) |>
    summarise(
      plant_richness = dplyr::n_distinct(scientificName),
      f1 = sum(count == 1),
      f2 = sum(count == 2),
      .groups = "drop"
    ) |>
    mutate(
      plant_chao1 = plant_richness + if_else(
        f2 > 0, f1^2 / (2 * f2), f1 * (f1 - 1) / 2),
      plant_completeness_raw = plant_richness / plant_chao1,
      # completeness x richness: observed richness discounted by how
      # completely it was sampled
      plant_complete_richness = plant_richness * plant_completeness_raw
    ) |>
    select(grid_id, plant_completeness_raw, plant_complete_richness)
}

# One row per grid cell: geometry from any retaining organism plus
# cross-organism summaries (organism count, first plant stats, summed
# animal richness across retaining organisms) and raw plant completeness.
retained_cell_stats <- function(cells,
                                raw_completeness = collect_plant_completeness_raw()) {
  stats <- cells |>
    st_drop_geometry() |>
    group_by(grid_id) |>
    summarise(n_retained        = n_distinct(organism),
              plant_completeness = dplyr::first(plant_completeness),
              animal_richness    = sum(animal_richness),
              animal_completeness = dplyr::first(animal_completeness),
              .groups = "drop")
  cells |>
    group_by(grid_id) |>
    slice(1) |>
    ungroup() |>
    select(grid_id) |>
    left_join(stats, by = "grid_id") |>
    left_join(raw_completeness, by = "grid_id")
}

fig_retained_cells <- function(retained,
                               raw_completeness = collect_plant_completeness_raw()) {
  retention_freq <- retained_cell_stats(retained, raw_completeness)

  # Zoom to the retained cells: the country-border layer is global, and in
  # EPSG:3857 it would otherwise stretch the panel to the whole world and
  # shrink the cells to specks. Crop it to the same extent so only the
  # Amazonian borders are drawn.
  bbox <- st_bbox(retention_freq)
  borders <- amani_country_borders() |>
    st_transform(st_crs(retention_freq)) |>
    st_crop(st_as_sfc(bbox))
  amazon <- amazon_outline(st_crs(retention_freq), bbox)

  # Land drawn as a light grey fill so water reads as the white background;
  # cells and white country boundaries go on top of it.
  ggplot(retention_freq) +
    geom_sf(data = borders, fill = AMANI_MAP_STYLE$land_fill, colour = NA) +
    geom_sf(data = amazon, fill = AMANI_MAP_STYLE$amazon_fill, colour = NA) +
    geom_sf(data = borders, fill = NA, colour = AMANI_MAP_STYLE$border_col,
            linewidth = AMANI_MAP_STYLE$border_lwd) +
    geom_sf(aes(fill = n_retained), color = AMANI_MAP_STYLE$cell_outline_col,
            linewidth = AMANI_MAP_STYLE$cell_outline_lwd) +
    # begin = 0.3: the palette's low end is near-white and would read as
    # empty on a map of few, scattered cells
    scale_fill_richness(name = "Organisms\nretaining cell",
                        begin = AMANI_MAP_STYLE$richness_begin,
                        breaks = scales::breaks_pretty()) +
    coord_sf(xlim = c(bbox[["xmin"]], bbox[["xmax"]]),
             ylim = c(bbox[["ymin"]], bbox[["ymax"]])) +
    labs(title = "Retained cells across organisms",
         subtitle = paste0("Grid cells retained after filtering; n = ",
                           n_distinct(retained$organism), " organisms, ",
                           nrow(retention_freq), " cells")) +
    theme_amani_map()
}

# Combined per-cell map: fill colour = number of organisms retaining the
# cell, fill transparency = RAW plant completeness x richness (floor 20%,
# log10 scale since the product spans orders of magnitude), and cells at or
# above the 95th percentile of that product get a thicker outline in the same
# colour as their fill.
# (Raw, not the values stored in the cell-filter reports: those are computed
# from the filled plant data and completeness is 1 in every cell — see the
# note in 02_animalAnalysis.R.)
fig_retained_cells_combined <- function(retained,
                                        raw_completeness = collect_plant_completeness_raw()) {
  cell_stats <- retained_cell_stats(retained, raw_completeness)
  q95 <- quantile(cell_stats$plant_complete_richness, 0.95, na.rm = TRUE)

  bbox <- st_bbox(cell_stats)
  borders <- amani_country_borders() |>
    st_transform(st_crs(cell_stats)) |>
    st_crop(st_as_sfc(bbox))
  amazon <- amazon_outline(st_crs(cell_stats), bbox)

  # Drop shadow for the highlighted cells: the same polygons nudged slightly
  # down-right (a fraction of one cell) and filled with translucent grey.
  hi_cells <- filter(cell_stats, plant_complete_richness >= q95)
  cell_bbox <- st_bbox(cell_stats[1, ])
  shadow <- hi_cells
  # sf's affine translation drops the CRS, so re-attach it
  st_geometry(shadow) <- st_set_crs(
    st_geometry(shadow) +
      c(AMANI_MAP_STYLE$shadow_shift * (cell_bbox[["xmax"]] - cell_bbox[["xmin"]]),
        -AMANI_MAP_STYLE$shadow_shift * (cell_bbox[["ymax"]] - cell_bbox[["ymin"]])),
    st_crs(hi_cells))

  ggplot() +
    geom_sf(data = borders, fill = AMANI_MAP_STYLE$land_fill, colour = NA) +
    geom_sf(data = amazon, fill = AMANI_MAP_STYLE$amazon_fill, colour = NA) +
    geom_sf(data = borders, fill = NA, colour = AMANI_MAP_STYLE$border_col,
            linewidth = AMANI_MAP_STYLE$border_lwd) +
    geom_sf(data = cell_stats,
            aes(fill = n_retained, alpha = plant_complete_richness),
            colour = NA) +
    scale_fill_richness(name = "Organisms\nretaining cell",
                        begin = AMANI_MAP_STYLE$richness_begin,
                        breaks = scales::breaks_pretty()) +
    scale_alpha(name = "Plant\ncompleteness x richness",
                range = c(AMANI_MAP_STYLE$alpha_floor, 1), trans = "log10") +
    geom_sf(data = shadow, fill = AMANI_MAP_STYLE$shadow_fill,
            alpha = AMANI_MAP_STYLE$shadow_alpha, colour = NA) +
    geom_sf(data = hi_cells,
            aes(colour = n_retained), fill = NA,
            linewidth = AMANI_MAP_STYLE$pct_outline_lwd) +
    scale_color_richness(begin = AMANI_MAP_STYLE$richness_begin,
                         guide = "none") +
    coord_sf(xlim = c(bbox[["xmin"]], bbox[["xmax"]]),
             ylim = c(bbox[["ymin"]], bbox[["ymax"]])) +
    labs(title = "Retained cells across organisms",
         subtitle = paste0("Colour = organisms retaining cell; transparency = plant ",
                           "completeness x richness (log10, floor 20%); >= 95th percentile (",
                           round(q95), ") outlined; n = ",
                           n_distinct(retained$organism), " organisms, ",
                           nrow(cell_stats), " cells")) +
    theme_amani_map()
}

# 3D version of the combined map (rayshader): column height = number of
# organisms retaining the cell, fill colour = RAW plant completeness x richness
# (log10). Heights are rasterised from the cell polygons (px_per_cell pixels per
# cell) and the flat colour map is draped over them as a texture. With
# file = NULL the interactive rgl window is left open (close with
# rgl::close3d()); with a path, a PNG snapshot is written and the window closed.
# Uses terra::/rayshader:: prefixes: attaching terra masks tidyr::extract.
fig_retained_cells_3d <- function(retained,
                                  raw_completeness = collect_plant_completeness_raw(),
                                  file = NULL,
                                  px_per_cell = 6,
                                  height_per_organism = 7,  # in pixels
                                  phi = 40, theta = -20, zoom = 0.6,
                                  windowsize = c(1100, 850)) {
  cell_stats <- retained_cell_stats(retained, raw_completeness)
  bbox <- st_bbox(cell_stats)
  borders <- amani_country_borders() |>
    st_transform(st_crs(cell_stats)) |>
    st_crop(st_as_sfc(bbox))
  amazon <- amazon_outline(st_crs(cell_stats), bbox)

  cell_bbox <- st_bbox(cell_stats[1, ])
  n_x <- round((bbox[["xmax"]] - bbox[["xmin"]]) / (cell_bbox[["xmax"]] - cell_bbox[["xmin"]]))
  n_y <- round((bbox[["ymax"]] - bbox[["ymin"]]) / (cell_bbox[["ymax"]] - cell_bbox[["ymin"]]))

  # Height matrix: organisms retaining each cell, 0 elsewhere
  template <- terra::rast(
    terra::ext(bbox[["xmin"]], bbox[["xmax"]], bbox[["ymin"]], bbox[["ymax"]]),
    ncols = n_x * px_per_cell, nrows = n_y * px_per_cell,
    crs = st_crs(cell_stats)$wkt)
  height_mat <- terra::rasterize(terra::vect(cell_stats), template,
                                 field = "n_retained", background = 0) |>
    rayshader::raster_to_matrix(verbose = FALSE)

  # Colour texture on exactly the same extent / pixel grid
  texture_gg <- ggplot() +
    geom_sf(data = borders, fill = AMANI_MAP_STYLE$land_fill, colour = NA) +
    geom_sf(data = amazon, fill = AMANI_MAP_STYLE$amazon_fill, colour = NA) +
    geom_sf(data = cell_stats, aes(fill = plant_complete_richness), colour = NA) +
    geom_sf(data = borders, fill = NA, colour = AMANI_MAP_STYLE$border_col,
            linewidth = AMANI_MAP_STYLE$border_lwd) +
    scale_fill_richness(begin = AMANI_MAP_STYLE$richness_begin, trans = "log10") +
    coord_sf(xlim = c(bbox[["xmin"]], bbox[["xmax"]]),
             ylim = c(bbox[["ymin"]], bbox[["ymax"]]), expand = FALSE) +
    theme_void() +
    theme(legend.position = "none")
  texture_file <- tempfile(fileext = ".png")
  on.exit(unlink(texture_file), add = TRUE)
  # rayshader matrices are transposed (rows = x), so width = nrow()
  ggsave(texture_file, texture_gg, width = nrow(height_mat) / 100,
         height = ncol(height_mat) / 100, dpi = 100, bg = "white")

  rayshader::plot_3d(png::readPNG(texture_file), height_mat,
                     zscale = 1 / height_per_organism,
                     solid = FALSE, shadow = FALSE, background = "white",
                     windowsize = windowsize, zoom = zoom, phi = phi, theta = theta)
  rayshader::render_camera(phi = phi, theta = theta, zoom = zoom)
  if (!is.null(file)) {
    # rgl can't draw legends: save one as a separate file next to the snapshot.
    rayshader::render_snapshot(file)
    rgl::close3d()

    legend_gg <- ggplot(cell_stats,
                        aes(plant_complete_richness, plant_complete_richness,
                            colour = plant_complete_richness)) +
      geom_point(alpha = 0) +
      scale_color_richness(begin = AMANI_MAP_STYLE$richness_begin,
                           trans = "log10",
                           breaks = scales::breaks_log(),
                           name = "Plant\ncompleteness x richness") +
      theme_void() +
      theme(legend.text = element_text(size = 10),
            legend.title = element_text(size = 10))
    ggsave(sub("\\.png$", "_legend.png", file), legend_gg,
           width = 2.2, height = 4.5, dpi = 300, bg = "white")
  }
  invisible(file)
}

retained_cells <- collect_retained_cells()
plant_completeness_raw <- collect_plant_completeness_raw()
fig_retained_cells(retained_cells, plant_completeness_raw)
fig_retained_cells_combined(retained_cells, plant_completeness_raw)
# Same map with numbered site labels: top cell per country + next 8 overall
fig_retained_cells_combined(retained_cells, plant_completeness_raw) |>
  label_retained_cells(retained_cell_stats(retained_cells, plant_completeness_raw))
# Interactive 3D window; pass file = "path.png" to save a snapshot instead
fig_retained_cells_3d(retained_cells, plant_completeness_raw, file = "output/retained_cells_3d.png")

# --- P3b. Retention across organisms (PhyloPic violins) ----------------------
# Two-panel figure, no maps: (left) % of records / cells / species retained
# after filtering, (right) the raw number retained (log scale). Each organism
# is a PhyloPic silhouette coloured by taxon group, jittered horizontally
# within its measure; violins show the cross-organism distribution.
#
# rphylopic needs the system library librsvg2-dev (its rsvg dependency).
# Silhouette UUIDs are hard-coded below so the figure renders without network
# access; they were picked by hand from candidate grids on phylopic.org.
# Credit the PhyloPic contributors in the figure caption:
#   rphylopic::get_attribution(uuid = unlist(PHYLOPIC_UUIDS))

library(rphylopic)

fig_sampling_retention <- function(pub,
                                   exclude = "allCombined",
                                   jitter_width = 0.5,
                                   height_pct = 7,   # silhouette size, left panel
                                   height_log = .2,   # silhouette size, right (log) panel
                                   silhouette_alpha = 0.5,
                                   seed = 7406) {
  mapping <- tibble(organism = names(ORGANISM_PHYLOPIC),
                    sil_uuid = unname(PHYLOPIC_UUIDS[ORGANISM_PHYLOPIC]))

  ret <- pub$summary |>
    filter(!organism_dir %in% exclude) |>
    select(organism, starts_with("pre_"), starts_with("post_")) |>
    pivot_longer(-organism, names_to = c("stage", "metric"), names_sep = "_") |>
    pivot_wider(names_from = stage) |>
    mutate(metric = factor(metric, c("records", "cells", "species"),
                           c("Records", "Cells", "Species")),
           pct = 100 * post / pre) |>
    left_join(distinct(pub$organisms, organism, taxon_group), by = "organism") |>
    inner_join(mapping, by = "organism")

  # Fixed (not position_jitter) so the figure is identical every time it is drawn
  set.seed(seed)
  ret$jx <- as.numeric(ret$metric) + runif(nrow(ret), -jitter_width, jitter_width)

  panel <- function(y, ylab, log = FALSE, height) {
    # Violin fill matches the maps' land colour; silhouettes coloured by
    # taxon group (invertebrate / vertebrate) via the shared pub palette
    g <- ggplot(ret, aes(metric, {{ y }})) +
      geom_violin(fill = AMANI_MAP_STYLE$land_fill, colour = NA,
                  scale = "width", width = 0.8) +
      geom_phylopic(aes(x = jx, uuid = sil_uuid, colour = taxon_group),
                    height = height, alpha = silhouette_alpha) +
      scale_color_taxon(name = NULL) +
      labs(x = NULL, y = ylab) +
      theme_amani_pub()
    if (log) {
      # Axis on the right; minor ticks at 2x-9x each decade (classic log ticks).
      # Minor ticks are a guide option here, not a theme element (ggplot2 >= 4).
      g + scale_y_log10(
        position = "right",
        labels = scales::label_comma(),
        minor_breaks = function(lims) {
          decades <- floor(log10(lims[1])):ceiling(log10(lims[2]))
          as.vector(outer(2:9, 10^decades))
        },
        guide = guide_axis(minor.ticks = TRUE))
    } else {
      # Minor ticks every 5%
      g + scale_y_continuous(minor_breaks = function(lims) seq(0, 100, 5),
                             guide = guide_axis(minor.ticks = TRUE))
    }
  }

  # post == 0 (nothing retained) is drawn at 1 on the log panel
  panel(pct, "% retained", height = height_pct) +
    panel(pmax(post, 1), "Number retained (log)", log = TRUE, height = height_log) +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom")
}

fig_sampling_retention(pub)

# --- P3c. Per-cell animal completeness & richness before vs after filtering ----
# Violin counterpart of the P3 item (c): per-cell animal Chao1 completeness and
# species richness distributions, before (all animal cells) vs after (cells
# retained by filtering). Read from the same plotdata_cell_filter_report.rds as
# collect_retained_cells(), but the $animal table, which spans every grid cell
# and carries the `retained` flag, so both stages come from one file. Formatted
# like fig_sampling_retention(): two panels (completeness, richness, both on
# linear scales), violins with scale = "width", and uniformly sized PhyloPic
# silhouettes above each organism's violins (coloured by taxon group) instead
# of x-axis text labels. The two violins overlap at the same x position: the
# retained (After) violin behind, filled by vertebrate vs. arthropod taxon
# group, with the translucent grey Before violin in front so the retained-cell
# distribution reads through it.

collect_cell_distributions <- function(output_dir = PUB_SETTINGS$output_dir) {
  files <- Sys.glob(file.path(output_dir, "*", "animal_model", "plot_data",
                              "plotdata_cell_filter_report.rds"))
  names(files) <- basename(dirname(dirname(dirname(files))))
  purrr::imap(files, \(f, org) {
    readRDS(f)$animal |>
      st_drop_geometry() |>
      mutate(organism_dir = org,
             organism = organism_label(org), .before = 1) |>
      select(organism_dir, organism, grid_id, retained,
             animal_completeness, animal_richness)
  }) |>
    bind_rows()
}

# Coarse clusters within taxon groups so related organisms sit together in
# fig_cell_distributions() (butterflies with butterflies, birds with birds,
# ...). The vector gives the cluster-level panel order; organisms missing from
# the mapping fall back to their own name (own cluster).
ORGANISM_CLUSTER <- c(
  Bees_combined       = "Bees",
  Cockroaches.Curated = "Cockroaches",
  Beetles             = "Beetles",
  Hemiptera           = "Bugs",
  Lepidoptera         = "Lepidoptera",
  Butterflies         = "Lepidoptera",
  `Butterfl-curated`  = "Lepidoptera",
  `Odonata-curated`   = "Odonata",
  `Odonata-combined`  = "Odonata",
  Odonata             = "Odonata",
  Frogs               = "Frogs",
  Finches             = "Birds",
  Antbirds            = "Birds",
  Tyrannidae          = "Birds",
  Mammals             = "Mammals",
  Squamata            = "Squamates"
)

fig_cell_distributions <- function(pub,
                                   cells = collect_cell_distributions(),
                                   exclude = "allCombined",
                                   violin_width = 0.8,
                                   violin_alpha = 0.65,     # Before (front) transparency
                                   silhouette_alpha = 0.9,
                                   silhouette_height = 0.06) {
  # Canonical pub ordering restricted to organisms with cell-filter data: pub
  # also covers organisms without a cell-filter report, and their empty factor
  # slots would push the violins and top-axis silhouettes out of alignment.
  # The canonical order puts invertebrates (arthropods) left, vertebrates right.
  lvls <- organism_levels(pub)
  dist <- cells |>
    filter(!organism_dir %in% exclude) |>
    # Before = every animal cell; After = the retained subset
    (\(x) bind_rows(mutate(x, stage = "Before"),
                    mutate(filter(x, retained), stage = "After")))() |>
    mutate(stage = factor(stage, levels = c("Before", "After"))) |>
    pivot_longer(c(animal_completeness, animal_richness),
                 names_to = "variable", values_to = "value") |>
    filter(organism %in% lvls, !is.na(value)) |>
    left_join(distinct(pub$organisms, organism, taxon_group), by = "organism")

  # Column order, restricted to organisms with cell-filter data (pub also
  # covers organisms without a cell-filter report, whose empty factor slots
  # would push the violins and top-axis silhouettes out of alignment):
  # arthropods left, vertebrates right; within each group, related organisms
  # cluster together via ORGANISM_CLUSTER (butterflies with butterflies, birds
  # with birds, ...), then descending gamma within those clusters
  ord <- dist |>
    distinct(organism, taxon_group) |>
    mutate(cluster = unname(coalesce(ORGANISM_CLUSTER[organism], organism)),
           gamma = pub$summary$gamma_best[match(organism, pub$summary$organism)]) |>
    arrange(factor(taxon_group, levels = c("Invertebrates", "Vertebrates")),
            factor(cluster, levels = unique(ORGANISM_CLUSTER)),
            desc(gamma)) |>
    pull(organism)
  dist <- mutate(dist, organism = factor(organism, ord))

  # One silhouette per organism, placed above its violins on the top axis;
  # x maps to the same factor levels as the violins so the columns line up
  sils <- tibble(organism = levels(dist$organism)) |>
    mutate(sil_uuid = unname(PHYLOPIC_UUIDS[ORGANISM_PHYLOPIC[organism]])) |>
    left_join(distinct(pub$organisms, organism, taxon_group), by = "organism") |>
    filter(!is.na(sil_uuid)) |>
    mutate(organism = factor(organism, levels(dist$organism)))

  panel <- function(var, ylab, minor_breaks = NULL) {
    d <- filter(dist, variable == var)
    # Silhouettes sit at the top of the data range, inside the headroom the
    # y-scale expansion adds above it. geom_phylopic() heights are in native
    # axis units, so express the height as a fraction of the panel's y range
    # to keep silhouettes uniformly sized across panels
    ytop <- max(d$value)
    sil_height <- silhouette_height * diff(range(d$value))

    # N cells per column, stacked and centred under the violins: Before (all
    # cells) in grey on top, After (retained) beneath in the taxon-group
    # colour of the retained violin
    yrange <- diff(range(d$value))
    y_before <- min(d$value) - 0.05 * yrange
    y_after  <- min(d$value) - 0.10 * yrange
    n_cells <- d |>
      count(organism, taxon_group, stage, name = "n") |>
      mutate(label = as.character(n))

    # Overlapped violins, no outlines: the retained (After) violin behind,
    # filled by taxon group; translucent grey Before in front so the
    # retained-cell distribution reads through it
    g <- ggplot(d, aes(organism, value)) +
      geom_violin(data = filter(d, stage == "After"),
                  aes(fill = taxon_group), colour = NA,
                  scale = "width", width = violin_width) +
      geom_violin(data = filter(d, stage == "Before"),
                  fill = AMANI_STAGE_COLS[["Before"]], colour = NA,
                  alpha = violin_alpha, scale = "width", width = violin_width) +
      geom_phylopic(data = mutate(sils, .y = ytop),
                    aes(x = organism, y = .y, uuid = sil_uuid,
                        colour = taxon_group),
                    inherit.aes = FALSE, height = sil_height,
                    alpha = silhouette_alpha) +
      geom_text(data = filter(n_cells, stage == "Before"),
                aes(x = organism, y = y_before, label = label),
                inherit.aes = FALSE, size = 1.7, colour = "grey40") +
      geom_text(data = filter(n_cells, stage == "After"),
                aes(x = organism, y = y_after, label = label,
                    colour = taxon_group),
                inherit.aes = FALSE, size = 1.7) +
      scale_fill_taxon(name = NULL) +
      scale_color_taxon(name = NULL) +
      labs(x = NULL, y = ylab) +
      theme_amani_pub() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      scale_y_continuous(expand = expansion(mult = c(0.15, 0.14)),
                         minor_breaks = minor_breaks,
                         guide = guide_axis(minor.ticks = !is.null(minor_breaks)))
    g
  }

  # Minor ticks every 5% of the completeness scale; richness left at defaults
  panel("animal_completeness", "Chao1 completeness per cell",
        minor_breaks = \(lims) seq(0, 1, 0.05)) +
    panel("animal_richness", "Species richness per cell") +
    plot_layout(guides = "collect") +
    plot_annotation(
      caption = paste("Numbers under each column: grid cells before filtering",
                      "(grey, all cells) and after (colour, retained cells).")) &
    theme(legend.position = "bottom")
}

fig_cell_distributions(pub)

# --- P4. Figure 2: Gamma richness ---------------------------------------------
# a. Ladder: observed -> analysed -> surrogate gamma -> full-model gamma (CIs)
# b. Expansion ratio (gamma / observed) with 1x reference
# c. Composition: observed vs predicted-unsampled share, total annotated
# d. (optional) richness shortfall vs sampling completeness
# e. Map: retained grid cells, coloured by number of organisms retaining each

# Retained cells per organism come from plotdata_cell_filter_report.rds
# (only organisms with that file are included). Carries per-cell plant and
# animal stats from the delta table so figures can encode them.

# --- P5. Figure 3: Extinction and IUCN ----------------------------------------
# a. Species lost (median, IQR, range) as count and % of gamma
# b. Sensitivity across frequency scenarios and spatial models
# c. IUCN category composition (100% stacked), threatened count annotated
# d. (optional) population losses

# --- P6. Build and save -------------------------------------------------------
# pub <- collect_pub_data(...); fig_sampling(pub); fig_gamma(pub);
# fig_extinction(pub); save_pub_fig(...)
