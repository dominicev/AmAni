# utils/pub_4_fig_gamma_uncertainty.R — Publication figures, step 4b:
# gamma uncertainty
#
#   fig_gamma_uncertainty_treemap()  same layout as fig_gamma_treemap(), but
#                                    tile area ~ gamma uncertainty, i.e. the
#                                    difference between the max and min gamma
#                                    estimates (error bar range). No
#                                    percentages are shown.
#
# Requires utils/pub_4_fig_gamma.R (spiral layout, organism picker, helpers).

# Organisms are picked on the same basis as the gamma figure (`organism_col`)
# so both figures show the same set; tile area uses the range
# `range_cols[2] - range_cols[1]`.
fig_gamma_uncertainty_treemap <- function(pub,
                                          range_cols   = c("gamma_lower", "gamma_upper"),
                                          organism_col = "gamma_lower",
                                          organisms    = NULL,
                                          labels       = TRUE,
                                          phylopics    = TRUE,
                                          sil_max_w    = 0.45,
                                          width_mm     = PUB_SETTINGS$fig_width_double,
                                          min_text_pt  = 4,
                                          max_text_pt  = 14) {
  if (is.null(organisms)) {
    organisms <- gamma_treemap_organisms(pub, organism_col)
  } else if (is.character(organisms)) {
    organisms <- tibble::tibble(organism = organisms, starred = FALSE)
  }

  d <- pub$summary |>
    dplyr::inner_join(organisms, by = "organism") |>
    dplyr::left_join(pub$organisms, by = c("organism_dir", "organism")) |>
    dplyr::mutate(min_est = .data[[range_cols[1]]],
                  max_est = .data[[range_cols[2]]],
                  total   = max_est - min_est) |>
    dplyr::filter(!is.na(total), total > 0) |>
    dplyr::arrange(dplyr::desc(total)) |>
    dplyr::mutate(uuid = unname(PHYLOPIC_UUIDS[ORGANISM_PHYLOPIC[organism]]))
  if (!phylopics) d$uuid <- NA_character_

  max_sz <- max_text_pt / ggplot2::.pt
  min_sz <- min_text_pt / ggplot2::.pt
  fit <- function(s) ifelse(is.na(s) | s < min_sz, NA_real_, pmin(s, max_sz))
  mm <- function(s) ifelse(is.na(s), 0, s / width_mm)

  d <- d |>
    dplyr::bind_cols(.spiral_layout(d$total)) |>
    dplyr::mutate(
      tw = xmax - xmin, rh = ymax - ymin, mid_y = (ymin + ymax) / 2,
      label = paste0(organism, dplyr::if_else(starred, "*", "", missing = "")),
      fill_col = unname(AMANI_TAXON_COLS[taxon_group]),
      # Vertical name along the left edge
      name_size = fit(pmin(0.9 * rh * width_mm / (0.6 * nchar(label)),
                           0.6 * tw * width_mm)),
      name_x = xmin + 0.55 * mm(name_size),
      name_r = xmin + 1.1 * mm(name_size),
      # Silhouette: as tall as the tile allows, clear of the name
      aspect = .phylopic_aspect(uuid),
      sil_h = ifelse(is.na(aspect), 0, pmin(0.9 * rh, sil_max_w * tw / aspect)),
      sil_w = ifelse(is.na(aspect), 0, sil_h * aspect),
      sil_x = pmax(xmin + tw * 0.5 - sil_w / 2, name_r + sil_w / 2),
      sil_r = sil_x + sil_w / 2,
      # Range width + "min-max" in the free space right of the silhouette
      tx0 = pmax(sil_r, name_r),
      tot_x = (tx0 + xmax) / 2,
      tot_lab = scales::comma(round(total)),
      rng_lab = paste0(scales::comma(round(min_est)), "\u2013",
                       scales::comma(round(max_est))),
      tot_size = fit(pmin(0.85 * (xmax - tx0) * width_mm / (0.6 * nchar(tot_lab)),
                          0.3 * rh * width_mm)),
      rng_size = fit(pmin(0.5 * tot_size,
                          0.85 * (xmax - tx0) * width_mm / (0.55 * nchar(rng_lab)))))

  p <- ggplot() +
    geom_rect(data = d,
              aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill_col),
              colour = NA) +
    scale_fill_identity() +
    scale_size_identity() +
    coord_equal(expand = FALSE) +
    theme_void()

  sil <- dplyr::filter(d, !is.na(uuid), sil_h > 0)
  p <- p + purrr::map(seq_len(nrow(sil)), \(i) rphylopic::geom_phylopic(
    data = sil[i, ], aes(x = sil_x, y = mid_y, uuid = uuid),
    height = sil$sil_h[i], fill = "white", alpha = 1, inherit.aes = FALSE))

  if (labels) {
    p <- p +
      geom_text(data = dplyr::filter(d, !is.na(name_size)),
                aes(x = name_x, y = mid_y, label = label, size = name_size),
                angle = 90, family = AMANI_FONT, fontface = "bold",
                colour = "white") +
      geom_text(data = dplyr::filter(d, !is.na(tot_size)),
                aes(x = tot_x, y = mid_y + 0.15 * mm(tot_size),
                    label = tot_lab, size = tot_size),
                family = AMANI_FONT, fontface = "bold", colour = "grey15",
                vjust = 0) +
      geom_text(data = dplyr::filter(d, !is.na(tot_size), !is.na(rng_size)),
                aes(x = tot_x, y = mid_y - 0.15 * mm(tot_size),
                    label = rng_lab, size = rng_size),
                family = AMANI_FONT, colour = "grey30", vjust = 1)
  }

  p + geom_rect(data = d,
                aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                fill = NA, colour = "white", linewidth = 0.5)
}
