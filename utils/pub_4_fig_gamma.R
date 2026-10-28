# utils/pub_4_fig_gamma.R — Publication figures, step 4: Figure 2
# (gamma richness)
#
#   gamma_treemap_organisms()  auto-pick organisms, comparing each umbrella
#                              dataset with its finer-grained alternatives
#   fig_gamma_treemap()        square divided into one rectangle per organism
#                              (area ~ richness), split observed vs unobserved
#   .spiral_layout()           clockwise spiral slice layout used by the treemap

# Clockwise spiral slicing of the rectangle [x0, x0 + w] x [y0, y0 + h]:
# item 1 takes a full-width strip along the top, item 2 a strip down the
# right of what is left, then bottom, then left, and so on inwards. With
# values sorted descending, the largest sits top-left and the order runs
# clockwise round to the bottom left. Areas are exactly proportional to values.
.spiral_layout <- function(values, x0 = 0, y0 = 0, w = 1, h = 1) {
  n <- length(values)
  out <- data.frame(xmin = numeric(n), xmax = numeric(n),
                    ymin = numeric(n), ymax = numeric(n))
  remaining <- sum(values)
  for (i in seq_len(n)) {
    a <- if (i == n) 1 else values[i] / remaining
    side <- c("top", "right", "bottom", "left")[(i - 1) %% 4 + 1]
    if (side == "top") {
      out[i, ] <- c(x0, x0 + w, y0 + h * (1 - a), y0 + h)
      h <- h * (1 - a)
    } else if (side == "right") {
      out[i, ] <- c(x0 + w * (1 - a), x0 + w, y0, y0 + h)
      w <- w * (1 - a)
    } else if (side == "bottom") {
      out[i, ] <- c(x0, x0 + w, y0, y0 + h * a)
      y0 <- y0 + h * a
      h <- h * (1 - a)
    } else {
      out[i, ] <- c(x0, x0 + w * a, y0, y0 + h)
      x0 <- x0 + w * a
      w <- w * (1 - a)
    }
    remaining <- remaining - values[i]
  }
  out
}

# Blend a colour towards light grey to "dull" it for unobserved richness
.dull <- function(col, amount = 0.6) {
  vapply(col, \(cl) grDevices::colorRampPalette(c(cl, "grey85"))(101)[round(amount * 100) + 1],
         character(1), USE.NAMES = FALSE)
}

# --- Organism set -------------------------------------------------------------------

# Base organisms with no overlap (no combined/curated duplicates, no
# allCombined). For umbrella datasets, the richer of the umbrella vs its
# alternatives is kept: if the alternatives together are larger, they replace
# the umbrella and get an asterisk next to their names.
GAMMA_TREEMAP_BASE <- c(
  "Flies", "Hemiptera", "Odonata-curated", "Mantodea", "Termites",
  "Cockroaches.Curated",
  "Frogs", "Mammals", "Squamates", "Antbirds", "Tyrannidae", "Finches")

# Each entry: umbrella dataset first, then its alternatives
GAMMA_TREEMAP_PAIRS <- list(
  c("Hymenoptera", "Ants_curated", "Bees_combined"),
  c("Beetles", "DungBeetles"),
  c("Lepidoptera", "Butterfl-curated")
)

gamma_treemap_organisms <- function(pub, value_col = "gamma_lower") {
  val <- function(org) {
    i <- match(org, pub$summary$organism)
    if (is.na(i)) NA_real_ else pub$summary[[value_col]][i]
  }
  resolve <- function(pair) {
    single <- pair[1]
    alts <- pair[-1]
    s <- val(single)
    alt_vals <- vapply(alts, val, numeric(1))
    alts_ok <- alts[!is.na(alt_vals)]
    p <- sum(alt_vals[!is.na(alt_vals)])
    if (!is.na(s) && (length(alts_ok) == 0 || s >= p))
      return(tibble::tibble(organism = single, starred = FALSE))
    # Alternatives used (umbrella missing, or smaller than the alternatives)
    starred <- !is.na(s) && length(alts_ok) > 0
    tibble::tibble(organism = alts_ok, starred = starred)
  }
  dplyr::bind_rows(
    c(list(tibble::tibble(organism = GAMMA_TREEMAP_BASE, starred = FALSE)),
      purrr::map(GAMMA_TREEMAP_PAIRS, resolve))
  )
}

# --- Treemap ------------------------------------------------------------------------

# Width / height of each PhyloPic silhouette, cached per session so repeated
# renders don't refetch. Falls back to 1 if an image can't be fetched.
.phylopic_aspect_cache <- new.env()
.phylopic_aspect <- function(uuids) {
  vapply(uuids, \(u) {
    if (is.na(u)) return(NA_real_)
    if (is.null(.phylopic_aspect_cache[[u]])) {
      img <- tryCatch(rphylopic::get_phylopic(u, format = "raster", height = 64),
                      error = \(e) NULL)
      .phylopic_aspect_cache[[u]] <- if (is.null(img)) 1 else dim(img)[2] / dim(img)[1]
    }
    .phylopic_aspect_cache[[u]]
  }, numeric(1), USE.NAMES = FALSE)
}

# Square divided into one rectangle per organism, area proportional to
# richness (`value_col`, e.g. "gamma_lower" or "surrogate_best"). Each
# rectangle is split left/right: known (observed) richness on the left in
# full taxon colour, predicted-unobserved on the right, dulled.
#   - organism name: vertical, along the left edge of the known strip
#   - silhouette: pure white, at the far right of the known strip (opposite
#     the name), tile height (capped at `sil_max_w` of the tile width), with
#     that organism's share of the total richness across all tiles over it
#     in large black text
#   - total predicted + "Max:" (upper end of `range_cols`): centred in the
#     space right of the silhouette
# Text is sized to fit its tile at the final width `width_mm`; anything
# that would fall below `min_text_pt` is dropped. `organisms`: tibble from
# gamma_treemap_organisms(), a character vector, or NULL (auto-pick).
fig_gamma_treemap <- function(pub,
                              value_col   = "gamma_lower",
                              range_cols  = c("gamma_lower", "gamma_upper"),
                              organisms   = NULL,
                              labels      = TRUE,
                              phylopics   = TRUE,
                              dull_amount = 0.6,
                              sil_max_w   = 0.45,
                              width_mm    = PUB_SETTINGS$fig_width_double,
                              min_text_pt = 4,
                              max_text_pt = 14) {
  if (is.null(organisms)) {
    organisms <- gamma_treemap_organisms(pub, value_col)
  } else if (is.character(organisms)) {
    organisms <- tibble::tibble(organism = organisms, starred = FALSE)
  }

  d <- pub$summary |>
    dplyr::inner_join(organisms, by = "organism") |>
    dplyr::left_join(pub$organisms, by = c("organism_dir", "organism")) |>
    dplyr::mutate(total = .data[[value_col]],
                  max_pred = .data[[range_cols[2]]]) |>
    dplyr::filter(!is.na(total), total > 0) |>
    dplyr::arrange(dplyr::desc(total)) |>
    dplyr::mutate(observed = pmin(observed_spp, total),
                  unobserved = total - observed,
                  f_obs = observed / total,
                  uuid = unname(PHYLOPIC_UUIDS[ORGANISM_PHYLOPIC[organism]]))
  if (!phylopics) d$uuid <- NA_character_

  # ggplot text size is in mm; 1 data unit = width_mm on the page. Sizes
  # below the readable minimum become NA (label dropped).
  max_sz <- max_text_pt / ggplot2::.pt
  min_sz <- min_text_pt / ggplot2::.pt
  fit <- function(s) ifelse(is.na(s) | s < min_sz, NA_real_, pmin(s, max_sz))
  mm <- function(s) ifelse(is.na(s), 0, s / width_mm)   # text size -> data units

  d <- d |>
    dplyr::bind_cols(.spiral_layout(d$total)) |>
    dplyr::mutate(
      tw = xmax - xmin, rh = ymax - ymin, mid_y = (ymin + ymax) / 2,
      kx1 = xmin + tw * f_obs,                      # known | unobserved border
      kw = kx1 - xmin,
      label = paste0(organism, dplyr::if_else(starred, "*", "", missing = "")),
      fill_known = unname(AMANI_TAXON_COLS[taxon_group]),
      fill_unobs = .dull(fill_known, dull_amount),
      # Vertical name: length fits the tile height, cap height fits the strip
      name_size = fit(pmin(0.9 * rh * width_mm / (0.6 * nchar(label)),
                           0.6 * kw * width_mm)),
      name_x = xmin + 0.55 * mm(name_size),
      name_r = xmin + 1.1 * mm(name_size),           # right edge of the name
      # Silhouette at the far right of the known strip (opposite the name),
      # as tall as the tile allows; pushed right past the known border only
      # if the name would otherwise be covered
      aspect = .phylopic_aspect(uuid),
      sil_h = ifelse(is.na(aspect), 0, pmin(0.9 * rh, sil_max_w * tw / aspect)),
      sil_w = ifelse(is.na(aspect), 0, sil_h * aspect),
      sil_x = pmax(kx1 - sil_w / 2, name_r + sil_w / 2),
      sil_r = sil_x + sil_w / 2,
      # % = this organism's share of the total richness across all tiles
      pct = paste0(round(100 * total / sum(d$total), 1), "%"),
      pct_size = fit(pmin(0.9 * sil_w * width_mm / (0.6 * nchar(pct)),
                          0.4 * rh * width_mm)),
      # Total + Max in the free space right of the silhouette / known strip
      tx0 = pmax(kx1, sil_r, sil_x + mm(pct_size) * 0.35 * nchar(pct)),
      tot_x = (tx0 + xmax) / 2,
      tot_lab = scales::comma(round(total)),
      max_lab = paste0("Max: ", scales::comma(round(max_pred))),
      tot_size = fit(pmin(0.85 * (xmax - tx0) * width_mm / (0.6 * nchar(tot_lab)),
                          0.3 * rh * width_mm)),
      max_size = fit(pmin(0.5 * tot_size,
                          0.85 * (xmax - tx0) * width_mm / (0.55 * nchar(max_lab)))))

  tiles <- dplyr::bind_rows(
    d |> dplyr::transmute(xmin, xmax = kx1, ymin, ymax, fill = fill_known),
    d |> dplyr::filter(kx1 < xmax) |>
      dplyr::transmute(xmin = kx1, xmax, ymin, ymax, fill = fill_unobs))

  p <- ggplot() +
    geom_rect(data = tiles,
              aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill),
              colour = NA) +
    scale_fill_identity() +
    scale_size_identity() +
    coord_equal(expand = FALSE) +
    theme_void()

  # geom_phylopic() takes height only as a fixed parameter, so one layer
  # per organism lets silhouettes scale with their tile
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
      geom_text(data = dplyr::filter(d, !is.na(pct_size)),
                aes(x = sil_x, y = mid_y, label = pct, size = pct_size),
                family = AMANI_FONT, fontface = "bold", colour = "black") +
      geom_text(data = dplyr::filter(d, !is.na(tot_size)),
                aes(x = tot_x, y = mid_y + 0.15 * mm(tot_size),
                    label = tot_lab, size = tot_size),
                family = AMANI_FONT, fontface = "bold", colour = "grey15",
                vjust = 0) +
      geom_text(data = dplyr::filter(d, !is.na(tot_size), !is.na(max_size),
                                     !is.na(max_pred)),
                aes(x = tot_x, y = mid_y - 0.15 * mm(tot_size),
                    label = max_lab, size = max_size),
                family = AMANI_FONT, colour = "grey30", vjust = 1)
  }

  # Organism outlines last so they sit over the silhouettes
  p + geom_rect(data = d,
                aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                fill = NA, colour = "white", linewidth = 0.5)
}

