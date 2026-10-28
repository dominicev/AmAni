# label_retained_cells.R — Add numbered site labels to the retained-cell map.
#
# Labels (a) the highest plant completeness x richness cell in each country and
# (b) the next `n_next` highest cells overall, numbered by their rank among all
# retained cells (1 = highest). Names come from RETAINED_CELL_SITES
# (utils/retained_cell_site_names.R); cells not logged there are labelled
# "<country> (cell <id>)". Labels are wrapped to two lines and set in Arial.
#
# Usage:
#   cell_stats <- retained_cell_stats(retained_cells, plant_completeness_raw)
#   fig_retained_cells_combined(retained_cells, plant_completeness_raw) |>
#     label_retained_cells(cell_stats)
#
# Requires: sf, dplyr, ggplot2, ggrepel, utils/pub_theme.R (AMANI_FONT,
# amani_country_borders()). Arial must be installed for the PNG/screen device;
# save PDFs with cairo_pdf (save_pub_fig()) so the font is embedded.

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

label_retained_cells <- function(plot, cell_stats, n_next = 8,
                                 sites = RETAINED_CELL_SITES,
                                 family = AMANI_FONT, size = 2.6) {
  lab <- select_label_cells(cell_stats, n_next) |>
    dplyr::left_join(dplyr::select(sites, grid_id, site_name), by = "grid_id") |>
    dplyr::mutate(
      site_name = dplyr::coalesce(site_name,
                                  paste0(country, " (cell ", grid_id, ")")),
      label = paste0(rank, ". ", .wrap_two_lines(site_name)))

  plot +
    ggrepel::geom_label_repel(
      data = lab, aes(X, Y, label = label), family = family, size = size,
      lineheight = 0.9, fill = ggplot2::alpha("white", 0.85),
      min.segment.length = 0, box.padding = 0.5, max.overlaps = Inf,
      seed = 3871, inherit.aes = FALSE)
}
