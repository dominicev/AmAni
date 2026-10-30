# utils/pub_3_fig_sampling.R — Publication figures, step 3: Figure 1
# (sampling and filtering)
#
# Uses the collectors in utils/pub_1_collect.R and the style in
# utils/pub_2_style.R. The violin figures need rphylopic (and the system
# library librsvg2-dev for its rsvg dependency); the 3D map needs rayshader,
# terra, png and rgl.
#
#   Retained-cell maps
#     fig_retained_cells()           cells coloured by organisms retaining them
#     fig_retained_cells_combined()  + plant completeness x richness as alpha,
#                                    >= 95th percentile outlined with a shadow
#     fig_retained_cells_3d()        rayshader version of the combined map
#     label_retained_cells()         add numbered site labels to a map
#     RETAINED_CELL_SITES            hand-curated site names per grid cell
#   Retention violins (PhyloPic silhouettes)
#     fig_sampling_retention()       % and number of records/cells/species kept
#     fig_cell_distributions()       per-cell completeness/richness before vs after

# =============================================================================
# Retained-cell maps
# =============================================================================

fig_retained_cells <- function(retained,
                               raw_completeness = collect_plant_completeness_raw()) {
  retention_freq <- retained_cell_stats(retained, raw_completeness)
  ctx <- .map_context(retention_freq)

  # Cells and country boundaries go on top of the land fill
  ggplot(retention_freq) +
    .basemap_fill(ctx) +
    .basemap_lines(ctx) +
    geom_sf(aes(fill = n_retained), color = AMANI_MAP_STYLE$cell_outline_col,
            linewidth = AMANI_MAP_STYLE$cell_outline_lwd) +
    # begin = 0.3: the palette's low end is near-white and would read as
    # empty on a map of few, scattered cells
    scale_fill_richness(name = "Organisms\nretaining cell",
                        begin = AMANI_MAP_STYLE$richness_begin,
                        breaks = scales::breaks_pretty()) +
    .coord_map(ctx) +
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
  ctx <- .map_context(cell_stats)

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
    .basemap_fill(ctx) +
    .basemap_lines(ctx) +
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
    .coord_map(ctx) +
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
# rgl::close3d()); with a path, a PNG snapshot (plus a separate _legend.png)
# is written and the window closed.
# Uses terra::/rayshader:: prefixes: attaching terra masks tidyr::extract.
fig_retained_cells_3d <- function(retained,
                                  raw_completeness = collect_plant_completeness_raw(),
                                  file = NULL,
                                  px_per_cell = 6,
                                  height_per_organism = 7,  # in pixels
                                  phi = 40, theta = -20, zoom = 0.6,
                                  windowsize = c(1500, 850),
                                  labels = FALSE,       # numbered site labels
                                  n_next = 8,
                                  label_top = 6,        # NULL: label every selected cell
                                  label_ranks = c(21, 26),  # extra ranks to always include
                                  # cells sharing one label (ranks listed before the name)
                                  label_groups = list(
                                    "Yasuni National Park" = c(2, 3, 8, 9),
                                    "Rio Urubu State Forest + BDFFP" = c(7, 11)),
                                  label_names = TRUE,  # TRUE: "rank. name" instead of rank only
                                  label_size = 2,
                                  label_lift = 4) {    # label height above column, in organisms
  cell_stats <- retained_cell_stats(retained, raw_completeness)
  ctx <- .map_context(cell_stats)
  bbox <- ctx$bbox

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

  # Colour texture on exactly the same extent / pixel grid; border lines go
  # over the cells here
  texture_gg <- ggplot() +
    .basemap_fill(ctx) +
    geom_sf(data = cell_stats, aes(fill = plant_complete_richness), colour = NA) +
    .basemap_lines(ctx) +
    scale_fill_richness(begin = AMANI_MAP_STYLE$richness_begin, trans = "log10") +
    .coord_map(ctx, expand = FALSE) +
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

  if (labels) {
    lab <- .label_table(cell_stats, n_next)
    # Keep the highest-ranked cells overall (country tops can rank lower), plus
    # any extra ranks and all members of label_groups
    keep <- if (!is.null(label_top)) lab$rank <= label_top else logical(nrow(lab))
    keep <- keep | lab$rank %in% c(label_ranks, unlist(label_groups))
    lab <- lab[keep, , drop = FALSE]

    # Map coordinates -> matrix indices (rows = x west to east, cols = y north to south)
    lab$mx <- (lab$X - bbox[["xmin"]]) / (bbox[["xmax"]] - bbox[["xmin"]]) * nrow(height_mat)
    lab$my <- (bbox[["ymax"]] - lab$Y) / (bbox[["ymax"]] - bbox[["ymin"]]) * ncol(height_mat)
    lab$mx <- pmin(pmax(round(lab$mx), 1), nrow(height_mat))
    lab$my <- pmin(pmax(round(lab$my), 1), ncol(height_mat))
    lab$h <- height_mat[cbind(lab$mx, lab$my)]

    # One label per entity: a single cell, or a group of cells sharing a name
    lab$ent <- as.character(lab$rank)
    lab$grp_name <- NA_character_
    for (g in names(label_groups)) {
      in_g <- lab$rank %in% label_groups[[g]]
      lab$ent[in_g] <- g
      lab$grp_name[in_g] <- g
    }
    ents <- lab |>
      dplyr::group_by(ent) |>
      dplyr::summarise(ax = round(mean(mx)), ay = round(mean(my)), h = max(h),
                       rank_min = min(rank),
                       ranks = paste(sort(rank), collapse = ", "),
                       name = dplyr::first(dplyr::coalesce(grp_name,
                                                           gsub("\n", " ", label))),
                       .groups = "drop") |>
      dplyr::arrange(rank_min) |>
      dplyr::mutate(text = dplyr::case_when(
        !label_names ~ ranks,
        ent %in% names(label_groups) ~ paste0(ranks, ". ", name),
        TRUE ~ name))

    # Stack labels vertically when they sit close to an earlier (higher-ranked) one
    near <- vapply(seq_len(nrow(ents)), function(i) {
      sum(sqrt((ents$ax[seq_len(i - 1)] - ents$ax[i])^2 +
                 (ents$ay[seq_len(i - 1)] - ents$ay[i])^2) < 40)
    }, numeric(1))
    ents$z <- ents$h + label_lift + 5 * near

    zs <- 1 / height_per_organism
    for (i in seq_len(nrow(ents))) {
      # Text only (alpha = 0 hides rayshader's line, which depth-tests away
      # behind neighbouring columns); leader lines are drawn below instead
      rayshader::render_label(height_mat, text = ents$text[i],
                              x = ents$ax[i], y = ents$ay[i],
                              z = ents$z[i], zscale = zs, relativez = FALSE,
                              textsize = label_size, alpha = 0,
                              freetype = FALSE, adjustvec = c(0, 0.5),
                              textcolor = "black")
      members <- lab[lab$ent == ents$ent[i], ]
      ig <- rgl::par3d(ignoreExtent = TRUE)
      for (k in seq_len(nrow(members))) {
        # Tagged "textline" so rayshader's software snapshot renderer draws it
        # (it ignores untagged rgl lines)
        rgl::lines3d(
          cbind(c(members$mx[k], ents$ax[i]) - nrow(height_mat) / 2,
                c(members$h[k], ents$z[i] + 3) / zs,   # +3: snapshot text sits above its anchor
                c(members$my[k], ents$ay[i]) - ncol(height_mat) / 2),
          color = "black", lwd = 2, lit = FALSE, tag = "textline")
      }
      rgl::par3d(ignoreExtent = ig)
    }
  }

  if (!is.null(file)) {
    # rgl can't draw legends: save one as a separate file next to the snapshot.
    rayshader::render_snapshot(file, line_radius = 0.5)
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

# --- Site labels -------------------------------------------------------------------
#
# label_retained_cells() labels (a) the highest plant completeness x richness
# cell in each country and (b) the next `n_next` highest cells overall,
# numbered by their rank among all retained cells (1 = highest). Names come
# from RETAINED_CELL_SITES; cells not logged there are labelled
# "<country> (cell <id>)". Labels are wrapped to two lines and set in Arial
# (must be installed for the PNG/screen device; save PDFs with save_pub_fig()
# so the font is embedded).
#
#   fig_retained_cells_combined(retained_cells, plant_completeness_raw) |>
#     label_retained_cells(retained_cell_stats(retained_cells, plant_completeness_raw))

# Break a label into two lines at the word boundary nearest its midpoint.
.wrap_two_lines <- function(x) {
  vapply(x, function(s) {
    words <- strsplit(s, " ", fixed = TRUE)[[1]]
    if (length(words) < 2) return(s)
    cum <- cumsum(nchar(words) + 1)
    k <- which.min(abs(cum - sum(nchar(words) + 1) / 2))
    paste(paste(words[seq_len(k)], collapse = " "),
          paste(words[-seq_len(k)], collapse = " "), sep = "\n")
  }, character(1), USE.NAMES = FALSE)
}

# Cells to label: top cell per country + next n_next overall.
select_label_cells <- function(cell_stats, n_next = 8) {
  countries <- amani_country_borders() |>
    sf::st_transform(sf::st_crs(cell_stats)) |>
    dplyr::select(country = name)
  pts <- suppressWarnings(sf::st_point_on_surface(cell_stats)) |>
    sf::st_join(countries)

  d <- pts |>
    dplyr::filter(!is.na(plant_complete_richness)) |>
    dplyr::mutate(rank = rank(-plant_complete_richness, ties.method = "min"),
                  X = sf::st_coordinates(geometry)[, 1],
                  Y = sf::st_coordinates(geometry)[, 2]) |>
    sf::st_drop_geometry()

  country_top <- d |>
    dplyr::filter(!is.na(country)) |>
    dplyr::group_by(country) |>
    dplyr::slice_max(plant_complete_richness, n = 1, with_ties = FALSE) |>
    dplyr::ungroup()
  next_top <- d |>
    dplyr::filter(!grid_id %in% country_top$grid_id) |>
    dplyr::arrange(rank) |>
    dplyr::slice_head(n = n_next)

  dplyr::bind_rows(country_top, next_top) |> dplyr::arrange(rank)
}

# Cells to label with their "rank. name" text (shared by the 2D and 3D maps).
.label_table <- function(cell_stats, n_next = 8, sites = RETAINED_CELL_SITES) {
  select_label_cells(cell_stats, n_next) |>
    dplyr::left_join(dplyr::select(sites, grid_id, site_name), by = "grid_id") |>
    dplyr::mutate(
      site_name = dplyr::coalesce(site_name,
                                  paste0(country, " (cell ", grid_id, ")")),
      label = paste0(rank, ". ", .wrap_two_lines(site_name)))
}

label_retained_cells <- function(plot, cell_stats, n_next = 8,
                                 sites = RETAINED_CELL_SITES,
                                 family = AMANI_FONT, size = 2.6) {
  lab <- .label_table(cell_stats, n_next, sites)

  plot +
    ggrepel::geom_label_repel(
      data = lab, aes(X, Y, label = label), family = family, size = size,
      lineheight = 0.9, fill = ggplot2::alpha("white", 0.85),
      min.segment.length = 0, box.padding = 0.5, max.overlaps = Inf,
      seed = 3871, inherit.aes = FALSE)
}

# Parks, reserves and field stations at the highlighted grid cells.
#
# One row per grid_id (EPSG:3857 cells, ~22 km wide). Where a biological
# station sits inside a park/reserve the label names both ("Park - Station").
# Used by label_retained_cells(); cells without
# a row there fall back to "<country> (cell <id>)".
#
#   set         why the cell is logged
#                 richness_top  top plant completeness x richness cells
#                 country_top   highest plant completeness x richness cell in a country
#                 retained_top  most-retained cells (>= 5 organisms)
#                 low_plant     retained by >= 3 organisms but poor plant sampling
#   confidence  how the label was established
#                 user          confirmed by the user on satellite imagery
#                 high/medium/low  matched from coordinates + web/own knowledge, unchecked
#                 unidentified  no site identified yet
#
# Labels are cell-level best matches, not exhaustive lists of what a cell holds.

RETAINED_CELL_SITES <- tibble::tribble(
  ~grid_id, ~set,            ~site_name,                                                              ~confidence,
  # --- Top plant completeness x richness ---------------------------------------
  21604,    "richness_top",  "Saul, Parc amazonien de Guyane",                                        "medium",
  17793,    "richness_top",  "Yasuni National Park, Napo River / Pompeya road",                       "medium",
  17617,    "richness_top",  "Yasuni National Park - PUCE Yasuni Scientific Station",                 "high",
  6916,     "richness_top",  "Tambopata National Reserve - Explorer's Inn / Posada Amazonas",         "medium",
  15008,    "richness_top",  "Amacayacu National Park",                                               "medium",
  18322,    "richness_top",  "Cuyabeno Wildlife Reserve",                                             "medium",
  16292,    "richness_top",  "Rio Urubu State Forest",                                                "user",
  17442,    "richness_top",  "Waorani territory / Yasuni National Park",                              "user",
  17792,    "richness_top",  "Western Yasuni, Napo River lodges",                                     "user",
  7264,     "richness_top",  "Los Amigos Conservation Concession - CICRA Biological Station",         "high",
  16290,    "richness_top",  "Western BDFFP, BR-174 / Ramal do Rio Branquinho",                       "user",
  6724,     "richness_top",  "Machu Picchu Historic Sanctuary",                                       "user",
  15347,    "richness_top",  "Confluence of Amazon and Napo Rivers",                                  "user",
  # --- Top cell in each remaining country (not yet identified) ------------------------
  2191,     "country_top",   "Bosque de Helechos Gigantes",                                           "user",
  19640,    "country_top",   "Waiwoi Village, Rupununi",                                              "user",
  21766,    "country_top",   "Central Suriname Nature Reserve",                                       "user",
  19077,    "country_top",   "Cerro de la Neblina",                                                   "user",
  # --- Most-retained cells ---------------------------------------------------------
  9975,     "retained_top",  "Cristalino State Park - Cristalino Jungle Lodge / RPPN",                "high",
  22489,    "retained_top",  "Kaw-Roura Marshland Reserve - Tresor Regional Nature Reserve",          "high",
  22664,    "retained_top",  "Montsinery, D5 road",                                                   "user",
  22839,    "retained_top",  "Montagne des Singes, Kourou",                                           "medium",
  7786,     "retained_top",  "Manu National Park - Cocha Cashu Biological Station",                   "high",
  14465,    "retained_top",  "Tamshiyacu-Tahuayo Regional Conservation Area",                         "medium",
  17435,    "retained_top",  "Jatun Sacha Biological Station, Napo River",                            "medium",
  17968,    "retained_top",  "Limoncocha Biological Reserve",                                         "medium",
  21959,    "retained_top",  "Nouragues Nature Reserve - Nouragues Research Station",                 "high",
  # --- Retained but poorly sampled plants ----------------------------------------------
  21678,    "low_plant",     "Rio Guamal",                                                            "user",
  7259,     "low_plant",     "Amarakaeri Communal Reserve",                                           "user",
  7083,     "low_plant",     "Amarakaeri Communal Reserve (west)",                                    "user",
  7437,     "low_plant",     "Madre de Dios River lodges near Manu",                                  "user",
  7445,     "low_plant",     "Finca Las Piedras Research Station",                                    "user",
  9588,     "low_plant",     "Balneario Batistao",                                                    "user",
  19730,    "low_plant",     "Finca Bella Vista Eduardo Bravo, La Drymophila, El Cedral, Los Ariscos, Hostal de la Montana", "user",
  12862,    "low_plant",     "Rural area near Casa Alto Rioja",                                       "user",
  18673,    "low_plant",     "Colombia-Ecuador border (unidentified)",                                "unidentified"
)

# =============================================================================
# Retention violins (PhyloPic silhouettes)
# =============================================================================
# Silhouette UUIDs are hard-coded in pub_2_style.R so the figures render
# without network access. Credit the PhyloPic contributors in the caption:
#   rphylopic::get_attribution(uuid = unname(PHYLOPIC_UUIDS))

# Two-panel figure, no maps: (left) % of records / cells / species retained
# after filtering, (right) the raw number retained (log scale). Each organism
# is a PhyloPic silhouette coloured by taxon group, jittered horizontally
# within its measure; violins show the cross-organism distribution.
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

# Per-cell animal Chao1 completeness and species richness, before (all animal
# cells) vs after (cells retained by filtering). Two panels on linear scales,
# violins with scale = "width", and a PhyloPic silhouette above each
# organism's column instead of x-axis text. The two violins overlap at the
# same x position: the retained (After) violin behind, filled by taxon group,
# with the translucent grey Before violin in front so the retained-cell
# distribution reads through it.
fig_cell_distributions <- function(pub,
                                   cells = collect_cell_distributions(),
                                   exclude = "allCombined",
                                   violin_width = 0.8,
                                   violin_alpha = 0.65,     # Before (front) transparency
                                   silhouette_alpha = 0.9,
                                   silhouette_height = 0.06) {
  # Restrict to organisms in pub: their empty factor slots would otherwise
  # push the violins and top-axis silhouettes out of alignment.
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

  # Column order: invertebrates left, vertebrates right; within each group,
  # taxonomic order from ORGANISM_TAXONOMY (Hymenoptera together, Lepidoptera
  # together, Holometabola contiguous, ...). Unclassified organisms go last
  # within their group, by descending gamma
  ord <- dist |>
    distinct(organism, taxon_group) |>
    mutate(tax_rank = match(organism, ORGANISM_TAXONOMY$organism),
           gamma = pub$summary$gamma_best[match(organism, pub$summary$organism)]) |>
    arrange(factor(taxon_group, levels = c("Invertebrates", "Vertebrates")),
            is.na(tax_rank), tax_rank, desc(gamma)) |>
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
    yrange <- diff(range(d$value))
    sil_height <- silhouette_height * yrange

    # N cells per column, stacked and centred under the violins: Before (all
    # cells) in grey on top, After (retained) beneath in the taxon-group
    # colour of the retained violin
    y_before <- min(d$value) - 0.05 * yrange
    y_after  <- min(d$value) - 0.10 * yrange
    n_cells <- d |>
      count(organism, taxon_group, stage, name = "n") |>
      mutate(label = as.character(n))

    ggplot(d, aes(organism, value)) +
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
