# Publication Figure 3: extinction and IUCN
#   fig_iucn_stack()   IUCN composition (100% stacked, Figure 2 organisms, mean
#                      species per simulation across all spatial models x
#                      scenarios) with ribbons into a second column that
#                      expands VU, EN, CR and Extinct, and a third column that
#                      splits each of those by organism

AMANI_IUCN_EXT_COL <- "#4A0A02"

# Mean species count per IUCN category for the Figure 2 organisms: each
# organism's counts are summed over models x scenarios and divided by the
# number of simulations, then (if `by` is NULL) summed over organisms. CR/EX
# is split into critically endangered (CR) and completely lost (Ext, 100%
# range loss). Ext is the modelled species loss (pub$extinction$median_spp_lost)
# for the same organism x model x scenario; CR = CR/EX - Ext.
.iucn_pooled_counts <- function(pub, value_col = "gamma_lower", by = NULL) {
  orgs <- gamma_treemap_organisms(pub, value_col = value_col)$organism
  key <- c("organism", "model", "scenario")

  per_org <- pub$iucn |>
    dplyr::filter(organism %in% orgs) |>
    dplyr::left_join(
      dplyr::select(pub$extinction, dplyr::all_of(key), ext = median_spp_lost),
      by = key) |>
    dplyr::mutate(
      ext = dplyr::if_else(iucn_category == "CR/EX",
                           pmin(dplyr::coalesce(ext, 0), count), 0),
      cr  = dplyr::if_else(iucn_category == "CR/EX", count - ext, 0)) |>
    dplyr::summarise(.by = organism,
                     n_sims = dplyr::n_distinct(model, scenario),
                     NT = sum(count[iucn_category == "NT/LC"]) / n_sims,
                     VU = sum(count[iucn_category == "VU"]) / n_sims,
                     EN = sum(count[iucn_category == "EN"]) / n_sims,
                     CR = sum(cr) / n_sims, Ext = sum(ext) / n_sims) |>
    dplyr::select(-n_sims)

  if (is.null(by)) dplyr::summarise(per_org, dplyr::across(-organism, sum)) else per_org
}

# Smooth S-shaped ribbon between two vertical segments
.ribbon_poly <- function(x0, x1, l_lo, l_hi, r_lo, r_hi, id, n = 40) {
  t <- seq(0, 1, length.out = n)
  s <- t^2 * (3 - 2 * t)
  data.frame(x = c(x0 + (x1 - x0) * t, rev(x0 + (x1 - x0) * t)),
             y = c(l_lo + (r_lo - l_lo) * s, rev(l_hi + (r_hi - l_hi) * s)),
             id = id)
}

.fmt_n <- function(n) ifelse(n > 0 & n < 1, "<1", scales::comma(round(n)))
.fmt_pct <- function(x) ifelse(x < 0.01, scales::percent(x, 0.1), scales::percent(x, 1))

# Push label positions apart (bottom to top, then back down from `top`) so
# neighbours are at least `gap` apart
.spread_y <- function(y, gap, top = Inf, bottom = -Inf) {
  o <- order(y); ys <- pmax(y[o], bottom)
  for (i in seq_along(ys)[-1]) ys[i] <- max(ys[i], ys[i - 1] + gap)
  for (i in rev(seq_along(ys)[-length(ys)])) ys[i] <- min(ys[i], ys[i + 1] - gap)
  if (length(ys) && ys[length(ys)] > top) ys <- ys - (ys[length(ys)] - top)
  y[o] <- ys
  y
}

# Layout of stacked organism segments between lo and hi (top = first organism)
.org_seg <- function(d, lo, hi) {
  d <- d[order(d$organism), ]
  d$share <- d$n / sum(d$n)
  d$ymax <- hi - (cumsum(d$share) - d$share) * (hi - lo)
  d$ymin <- hi - cumsum(d$share) * (hi - lo)
  d
}

# Layout: all-species stack in the centre; VU and EN fan out to the left and
# CR and Extinct to the right (each side rescaled), and every expanded bar is
# split by organism on the outer side. A reference bar (organisms in all
# species) at the far right gives the baseline for over/under-representation.
# width_mm: final figure width, used to size text relative to bar height
# (text below min_text_pt moves outside the bar on a leader line).
fig_iucn_stack <- function(pub, value_col = "gamma_lower", gap = 0.04,
                           width_mm = PUB_SETTINGS$fig_width_double,
                           min_text_pt = 5, max_text_pt = 7) {
  # --- Geometry (x units) --------------------------------------------------------
  wo <- 1.6                         # organism bar width
  xlo <- 5.2                        # left organism bars: left edge
  xle <- xlo + wo + 1.8; we <- 0.5  # left expanded bars (VU, EN)
  xc  <- xle + we + 1.8; wc <- 1.3  # centre bar
  xre <- xc + wc + 1.8              # right expanded bars (CR, Ext)
  xro <- xre + we + 1.8             # right organism bars: left edge
  xf  <- xro + wo + 4.1             # reference bar
  x_max <- xf + wo + 5.2
  x_min <- 0.2
  k <- 6                            # one y unit = k x units (coord_fixed)
  mm_x <- width_mm * 0.72 / (x_max - x_min)
  mm_y <- k * mm_x
  max_sz <- max_text_pt / ggplot2::.pt
  min_sz <- min_text_pt / ggplot2::.pt

  # --- Centre column -------------------------------------------------------------
  cnt <- unlist(.iucn_pooled_counts(pub, value_col))
  p <- cnt / sum(cnt)
  cats <- c("NT", "VU", "EN", "CR", "Ext")
  cat_lab <- c(NT = "NT/LC", VU = "VU", EN = "EN", CR = "CR", Ext = "Extinct")
  fills <- c(NT = AMANI_IUCN_COLS[["NT/LC"]], VU = AMANI_IUCN_COLS[["VU"]],
             EN = AMANI_IUCN_COLS[["EN"]], CR = AMANI_IUCN_COLS[["CR/EX"]],
             Ext = AMANI_IUCN_EXT_COL)

  centre <- data.frame(cat = cats, share = p[cats], n = cnt[cats])
  centre$ymin <- cumsum(c(0, head(centre$share, -1)))
  centre$ymax <- cumsum(centre$share)
  centre$fill <- fills[centre$cat]

  # Expanded bars: each side rescaled to fill 0-1 (bottom to top as listed)
  expand_side <- function(side_cats, x) {
    rs <- p[side_cats] / sum(p[side_cats])
    h <- rs * (1 - gap * (length(side_cats) - 1))
    ymax <- cumsum(h + gap) - gap
    data.frame(cat = side_cats, ymin = ymax - h, ymax = ymax, share = p[side_cats],
               n = cnt[side_cats], x = x, fill = fills[side_cats],
               side = if (x < xc) "L" else "R")
  }
  expanded <- rbind(expand_side(c("VU", "EN"), xle), expand_side(c("CR", "Ext"), xre))

  # --- Organism bars ---------------------------------------------------------------
  by_org <- .iucn_pooled_counts(pub, value_col, by = "organism")
  lv <- organism_levels(pub)
  lv <- lv[lv %in% by_org$organism]
  org_long <- by_org |>
    tidyr::pivot_longer(-organism, names_to = "cat", values_to = "n") |>
    dplyr::mutate(organism = factor(organism, levels = lv))
  tax <- pub$organisms |>
    dplyr::distinct(organism, taxon_group) |>
    dplyr::filter(organism %in% lv)
  tax_col <- setNames(AMANI_TAXON_COLS[tax$taxon_group[match(lv, tax$organism)]], lv)

  org_bars <- do.call(rbind, lapply(seq_len(nrow(expanded)), \(i) {
    e <- expanded[i, ]
    d <- .org_seg(org_long[org_long$cat == e$cat, ], e$ymin, e$ymax)
    d$side <- e$side
    d$x <- if (e$side == "L") xlo else xro
    d
  }))
  ref <- org_long |>
    dplyr::summarise(n = sum(n), .by = organism) |>
    .org_seg(0, 1)
  ref$cat <- "All"; ref$side <- "R"; ref$x <- xf

  segs <- rbind(org_bars, ref) |>
    dplyr::mutate(
      fill = tax_col[as.character(organism)],
      mid_y = (ymin + ymax) / 2,
      h_mm = (ymax - ymin) * mm_y,
      lab_in = paste0(.fmt_n(n), "\n(", .fmt_pct(share), ")"),
      # text sized to the segment height (two lines) and bar width
      sz = pmin(0.85 * h_mm / 2.4,
                0.95 * wo * mm_x / (0.6 * pmax(nchar(.fmt_n(n)), nchar(.fmt_pct(share)) + 2)),
                max_sz),
      inside = sz >= min_sz,
      uuid = unname(PHYLOPIC_UUIDS[ORGANISM_PHYLOPIC[as.character(organism)]]),
      aspect = .phylopic_aspect(uuid),
      sil_h = pmin(0.85 * (ymax - ymin), 0.85 * wo / (k * aspect)),
      lab_out = ifelse(side == "L",
                       paste0(organism, "  ", .fmt_n(n), " (", .fmt_pct(share), ")"),
                       paste0(organism, "  ", .fmt_n(n), " (", .fmt_pct(share), ")")))
  inside <- dplyr::filter(segs, inside)
  # Silhouette for every segment above 5% of its column, sized to fit
  sil <- dplyr::filter(segs, !is.na(uuid), share > 0.05, sil_h * mm_y >= 1)

  # Segments too short for text: name + count on a leader line to the outer
  # side. Segments averaging under 0.5 species are left unlabelled.
  outside <- segs |>
    dplyr::filter(!inside, n >= 0.5) |>
    dplyr::mutate(edge = ifelse(side == "L", x, x + wo),
                  lab_x = ifelse(side == "L", edge - 0.9, edge + 0.9)) |>
    dplyr::group_by(cat) |>
    dplyr::mutate(lab_y = .spread_y(mid_y, 1.3 * min_sz / mm_y, top = 0.99, bottom = 0.01)) |>
    dplyr::ungroup()

  # --- Ribbons -------------------------------------------------------------------
  ribbon_set <- \(rows) do.call(rbind, lapply(rows, \(r) {
    e <- expanded[expanded$cat == r$cat, ]; ce <- centre[centre$cat == r$cat, ]
    if (r$kind == "S") {
      if (e$side == "L") .ribbon_poly(xle + we, xc, e$ymin, e$ymax, ce$ymin, ce$ymax, r$cat)
      else .ribbon_poly(xc + wc, xre, ce$ymin, ce$ymax, e$ymin, e$ymax, r$cat)
    } else {
      if (e$side == "L") .ribbon_poly(xlo + wo, xle, e$ymin, e$ymax, e$ymin, e$ymax, r$cat)
      else .ribbon_poly(xre + we, xro, e$ymin, e$ymax, e$ymin, e$ymax, r$cat)
    }
  }))
  rib_s <- ribbon_set(lapply(expanded$cat, \(ct) list(cat = ct, kind = "S")))
  rib_f <- ribbon_set(lapply(expanded$cat, \(ct) list(cat = ct, kind = "F")))
  rib_s$fill <- fills[rib_s$id]; rib_f$fill <- fills[rib_f$id]

  # --- Labels: count, with % of all species below ----------------------------------
  lab_cat <- expanded |>
    dplyr::mutate(label = paste0(cat_lab[cat], "\n", .fmt_n(n), "\n(",
                                 scales::percent(share, 0.1), ")"),
                  lx = ifelse(side == "L", xlo + wo + 0.15, xre + we + 0.15)) |>
    dplyr::group_by(side) |>
    dplyr::mutate(ly = .spread_y((ymin + ymax) / 2, 3.4 * min_sz / mm_y, top = 1, bottom = 0.1)) |>
    dplyr::ungroup()
  lab_nt <- centre[centre$cat == "NT", ] |>
    dplyr::mutate(label = paste0("NT/LC\n", .fmt_n(n), "\n(", .fmt_pct(share), ")"),
                  y = (ymin + ymax) / 2)

  # --- Plot ------------------------------------------------------------------------
  # Fills are literal colours (IUCN and taxon colours share one scale); the
  # legend lists taxon groups only.
  ggplot2::ggplot() +
    ggplot2::geom_polygon(data = rib_s, ggplot2::aes(x, y, group = id, fill = fill),
                          alpha = 0.45) +
    ggplot2::geom_polygon(data = rib_f, ggplot2::aes(x, y, group = id, fill = fill),
                          alpha = 0.25) +
    ggplot2::geom_rect(data = centre, ggplot2::aes(xmin = xc, xmax = xc + wc,
                                                   ymin = ymin, ymax = ymax, fill = fill)) +
    ggplot2::geom_rect(data = expanded, ggplot2::aes(xmin = x, xmax = x + we,
                                                     ymin = ymin, ymax = ymax, fill = fill)) +
    ggplot2::geom_rect(data = segs, ggplot2::aes(xmin = x, xmax = x + wo,
                                                 ymin = ymin, ymax = ymax, fill = fill),
                       colour = "white", linewidth = 0.2) +
    purrr::map(seq_len(nrow(sil)), \(i) rphylopic::geom_phylopic(
      data = sil[i, ], ggplot2::aes(x = x + wo / 2, y = mid_y, uuid = uuid),
      height = sil$sil_h[i], fill = "white", alpha = 1, inherit.aes = FALSE)) +
    ggplot2::geom_text(data = inside,
                       ggplot2::aes(x = x + wo / 2, y = mid_y, label = lab_in, size = sz),
                       colour = "black", fontface = "bold", lineheight = 0.85) +
    ggplot2::geom_segment(data = outside,
                          ggplot2::aes(x = edge, xend = ifelse(side == "L", lab_x + 0.1, lab_x - 0.1),
                                       y = mid_y, yend = lab_y),
                          colour = "grey40", linewidth = 0.2) +
    ggplot2::geom_text(data = outside,
                       ggplot2::aes(x = lab_x, y = lab_y, label = lab_out,
                                    hjust = ifelse(side == "L", 1, 0)),
                       size = min_sz, colour = "black") +
    ggplot2::geom_text(data = lab_cat,
                       ggplot2::aes(x = lx, y = ly, label = label),
                       size = min_sz * 1.1, hjust = 0, lineheight = 0.85) +
    ggplot2::geom_text(data = lab_nt,
                       ggplot2::aes(x = xc + wc / 2, y = y, label = label),
                       size = min_sz * 1.1, lineheight = 0.85) +
    ggplot2::scale_size_identity() +
    ggplot2::scale_fill_identity(guide = "legend", name = "Taxon group",
                                 breaks = unname(AMANI_TAXON_COLS),
                                 labels = names(AMANI_TAXON_COLS)) +
    ggplot2::coord_fixed(ratio = k, expand = FALSE, clip = "off") +
    ggplot2::scale_x_continuous(
      limits = c(x_min, x_max),
      breaks = c(xlo + wo / 2, xle + we / 2, xc + wc / 2, xre + we / 2,
                 xro + wo / 2, xf + wo / 2),
      labels = c("Organisms in\ncategory", "VU, EN\n(rescaled)", "All species",
                 "CR, Extinct\n(rescaled)", "Organisms in\ncategory",
                 "Organisms in\nall species")) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(0.01)) +
    ggplot2::labs(x = NULL, y = NULL) +
    theme_amani_pub() +
    ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                   axis.ticks.y = ggplot2::element_blank(),
                   axis.title.y = ggplot2::element_blank(),
                   panel.grid.major = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank(),
                   legend.position = "bottom")
}


# --- Figure 3a/b: species lost ---------------------------------------------------

# Species lost per organism x model x scenario (Figure 2 organisms) with the
# gamma richness it is expressed against and the taxon group.
.extinction_loss_data <- function(pub, value_col = "gamma_lower") {
  orgs <- gamma_treemap_organisms(pub, value_col = value_col)$organism
  pub$extinction |>
    dplyr::filter(organism %in% orgs) |>
    dplyr::inner_join(
      dplyr::select(pub$summary, organism, gamma = dplyr::all_of(value_col)),
      by = "organism") |>
    dplyr::left_join(dplyr::distinct(pub$organisms, organism, taxon_group),
                     by = "organism") |>
    dplyr::mutate(organism = factor(organism, levels = rev(organism_levels(pub))))
}

# a. Species lost for the headline model x scenario (PUB_SETTINGS$headline):
# median (point), IQR (thick bar) and min-max across replicates (thin line),
# as a count (log scale) and as % of gamma richness.
fig_species_lost <- function(pub, value_col = "gamma_lower",
                             headline = PUB_SETTINGS$headline) {
  d <- .extinction_loss_data(pub, value_col) |>
    dplyr::filter(model == headline$model, scenario == headline$scenario)

  # Zero losses cannot sit on a log axis: counts are floored at 1
  d <- dplyr::mutate(d, dplyr::across(dplyr::ends_with("_spp_lost"),
                                                \(x) pmax(x, 1)))
  d_pct <- .extinction_loss_data(pub, value_col) |>
    dplyr::filter(model == headline$model, scenario == headline$scenario)

  panel <- function(d, scale, xlab, log = FALSE) {
    ggplot2::ggplot(d, ggplot2::aes(y = organism, colour = taxon_group)) +
      ggplot2::geom_linerange(ggplot2::aes(xmin = min_spp_lost / scale,
                                           xmax = max_spp_lost / scale),
                              linewidth = 0.3) +
      ggplot2::geom_linerange(ggplot2::aes(xmin = q25_spp_lost / scale,
                                           xmax = q75_spp_lost / scale),
                              linewidth = 1.4) +
      ggplot2::geom_point(ggplot2::aes(x = median_spp_lost / scale),
                          shape = 21, fill = "white", size = 1.1, stroke = 0.4) +
      scale_color_taxon(name = "Taxon group") +
      (if (log) ggplot2::scale_x_log10(labels = scales::comma)
       else ggplot2::scale_x_continuous(labels = scales::label_percent(scale = 1))) +
      ggplot2::labs(x = xlab, y = NULL) +
      theme_amani_pub()
  }
  a1 <- panel(d, 1, "Count (log scale)", log = TRUE)
  a2 <- panel(d_pct, d_pct$gamma / 100, "% of gamma richness") +
    ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                   axis.ticks.y = ggplot2::element_blank())
  patchwork::wrap_plots(a1, a2, guides = "collect") &
    ggplot2::theme(legend.position = "bottom")
}

# b. Sensitivity: % of gamma richness lost for every spatial model x
# frequency scenario (colour = model, shape = scenario).
fig_extinction_sensitivity <- function(pub, value_col = "gamma_lower",
                                       sensitivity = PUB_SETTINGS$sensitivity,
                                       y_labels = TRUE) {
  d <- .extinction_loss_data(pub, value_col) |>
    dplyr::filter(model %in% sensitivity$models,
                  scenario %in% sensitivity$scenarios) |>
    dplyr::mutate(pct_lost = 100 * median_spp_lost / gamma,
                  model = factor(model, levels = sensitivity$models),
                  scenario = factor(scenario, levels = sensitivity$scenarios))

  ggplot2::ggplot(d, ggplot2::aes(x = pct_lost, y = organism,
                                  colour = model, shape = scenario)) +
    ggplot2::geom_point(size = 1.4, alpha = 0.85,
                        position = ggplot2::position_dodge(width = 0.5)) +
    ggplot2::scale_x_continuous(labels = scales::label_percent(scale = 1)) +
    ggplot2::scale_colour_viridis_d(option = PUB_SETTINGS$bright_palette,
                                    end = 0.8, name = "Spatial model") +
    ggplot2::scale_shape_manual(values = c(16, 17, 15, 18), name = "Scenario") +
    ggplot2::labs(x = "Species lost (% of gamma richness, median)", y = NULL) +
    ggplot2::guides(colour = ggplot2::guide_legend(nrow = 1, order = 1),
                    shape = ggplot2::guide_legend(nrow = 2, order = 2)) +
    theme_amani_pub() +
    ggplot2::theme(legend.position = "bottom", legend.box = "vertical",
                   legend.box.just = "left", legend.margin = ggplot2::margin(0, 0, 0, 0),
                   legend.spacing.y = grid::unit(0, "mm")) +
    (if (y_labels) NULL else
      ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                     axis.ticks.y = ggplot2::element_blank()))
}

# Figure 3 composite: a (species lost) | b (sensitivity) over c (IUCN ribbons)
fig_extinction <- function(pub, value_col = "gamma_lower") {
  # Panel a is itself two plots; wrap_elements() makes it one tagged unit
  a <- patchwork::wrap_elements(full = fig_species_lost(pub, value_col))
  b <- fig_extinction_sensitivity(pub, value_col, y_labels = FALSE)
  c <- fig_iucn_stack(pub, value_col)
  ((a | b) / c) +
    patchwork::plot_layout(heights = c(1, 0.85)) +
    patchwork::plot_annotation(tag_levels = "a")
}
