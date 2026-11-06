# utils/func_report_plots.R — Plot builders for the per-organism reports
#
# Sourced by AmAni.summarize.R, which loads the packages (tidyverse, sf,
# patchwork, ggrepel). Uses organism_label() and
# load_dataset_summary_figure_inputs() from utils/func_core.R.
# Called by summarize_all() / summarize_organism():
#
#   build_dataset_summary_figure()        before/after filtering stats figure
#   build_cell_filter_map_panels()        animal vs plant cell-filter maps
#   build_completeness_correlation_plot() plant vs animal completeness scatter
#   write_cell_filter_report_pdf()        assembles the above into a 2-page PDF
#
# The cell-filter maps are rebuilt from plotdata_cell_filter_report.rds (saved
# by 02_animalAnalysis.R alongside the cell_filter_report_maps_*.pdf), so they
# match that PDF but can be arranged in a single grid.

# --- Dataset before/after summary figure ------------------------------------

#' Build a before/after cell-filtering summary figure.
#'
#' @param thresholds Data frame/tibble with columns `label`, `value` (both
#'   character). One row per filter threshold to display in the text card.
#' @param counts Data frame/tibble with columns `metric` (character),
#'   `before` (numeric), `after` (numeric). Retention is computed as
#'   after / before.
#' @param dists Data frame/tibble with columns `variable`, `stage`
#'   ("Before"/"After"), and the five-number summary `min`, `q1`, `med`,
#'   `q3`, `max`, plus `mean`. One row per variable x stage.
#' @param log_vars Character vector of `variable` names to draw on a log10
#'   y-axis (default richness/occurrences).
#' @param stage_cols Named colour vector for the two stages.
#' @param title,subtitle Overall plot title/subtitle.
#' @return A patchwork object (also printable/ggsave-able), or NULL if there
#'   is no usable data (e.g. counts entirely NA).
build_dataset_summary_figure <- function(
    thresholds,
    counts,
    dists,
    log_vars = c("Richness", "Occurrences"),
    stage_cols = c(Before = "grey70", After = "#2c7fb8"),
    title = "Dataset before vs after cell filtering",
    subtitle = paste(
      "Boxplots show per-cell min/Q1/median/Q3/max; white diamond = mean.",
      "Grey = before, blue = after."
    )) {

  if (is.null(counts) || nrow(counts) == 0 || all(is.na(counts$before))) return(NULL)

  counts <- counts |>
    dplyr::mutate(
      pct_retained = after / before,
      metric = factor(metric, levels = unique(metric))
    )

  dists <- dists |>
    dplyr::mutate(
      stage = factor(stage, levels = c("Before", "After")),
      variable = factor(variable, levels = unique(variable))
    )

  # --- Panel A: filter thresholds (text card) -------------------------------
  p_thresh <- thresholds |>
    dplyr::mutate(y = rev(seq_len(dplyr::n()))) |>
    ggplot() +
    geom_text(aes(x = 0, y = y, label = label), hjust = 0, size = 3.2) +
    geom_text(aes(x = 1, y = y, label = value), hjust = 1, size = 3.2, fontface = "bold") +
    scale_x_continuous(limits = c(0, 1), expand = expansion(mult = 0.02)) +
    scale_y_continuous(limits = c(0, nrow(thresholds) + 1)) +
    labs(title = "Filter thresholds") +
    theme_void() +
    theme(plot.title = element_text(face = "bold", size = 11, hjust = 0))

  # --- Panel B: retention bars (before = full, after = retained portion) -----
  p_counts <- counts |>
    ggplot() +
    geom_rect(aes(xmin = 0, xmax = 1, ymin = as.numeric(metric) - 0.3,
                  ymax = as.numeric(metric) + 0.3),
              fill = stage_cols[["Before"]]) +
    geom_rect(aes(xmin = 0, xmax = pct_retained, ymin = as.numeric(metric) - 0.3,
                  ymax = as.numeric(metric) + 0.3),
              fill = stage_cols[["After"]]) +
    geom_text(aes(x = pmin(pct_retained + 0.02, 0.98), y = as.numeric(metric),
                  label = paste0(before, " to ", after, "  (",
                                 scales::percent(pct_retained, accuracy = 0.1), ")")),
              hjust = 0, size = 3.1) +
    scale_x_continuous(limits = c(0, 1), labels = scales::percent,
                       expand = expansion(mult = c(0, 0.02))) +
    scale_y_continuous(breaks = seq_len(nrow(counts)), labels = levels(counts$metric),
                       limits = c(0.5, nrow(counts) + 0.5)) +
    labs(title = "Retention after filtering", x = "Proportion of full dataset kept", y = NULL) +
    theme_minimal(base_size = 10) +
    theme(
      plot.title = element_text(face = "bold", size = 11),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank()
    )

  # --- Panels C+: five-number-summary boxplots, before vs after --------------
  box_panel <- function(var) {
    d <- dplyr::filter(dists, variable == var)
    p <- ggplot(d, aes(x = stage, fill = stage)) +
      geom_boxplot(
        aes(ymin = min, lower = q1, middle = med, upper = q3, ymax = max),
        stat = "identity", width = 0.55, colour = "grey30"
      ) +
      geom_point(aes(y = mean), shape = 23, size = 2, fill = "white", colour = "grey20") +
      scale_fill_manual(values = stage_cols, guide = "none") +
      labs(title = paste(var, "per cell"), x = NULL, y = NULL) +
      theme_minimal(base_size = 10) +
      theme(
        plot.title = element_text(face = "bold", size = 10),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank()
      )
    if (var %in% log_vars) p <- p + scale_y_log10()
    p
  }

  box_panels <- lapply(levels(dists$variable), box_panel)

  # --- Assemble: top row = thresholds + retention; then boxplots two-per-row -
  top <- p_thresh + p_counts + plot_layout(widths = c(1, 1))
  boxes <- if (length(box_panels) > 0) wrap_plots(box_panels, ncol = 2) else NULL

  fig <- if (!is.null(boxes)) {
    top / boxes + plot_layout(heights = c(1, ceiling(length(box_panels) / 2)))
  } else {
    top
  }

  fig +
    plot_annotation(
      title = title,
      subtitle = subtitle,
      theme = theme(
        plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(size = 9, colour = "grey30")
      )
    )
}

# --- Cell-filter report ---------------------------------------------------------

# Replace plant Chao1 completeness in a cell-filter report with the value from
# the RAW filtered plant counts. Reports saved before this fix computed it from
# the filled plant data (fill_plants = TRUE), which is exactly 1 in every cell.
# Uses collect_plant_completeness_raw() from pub_1_collect.R.
use_raw_plant_completeness <- function(cfr,
    counts_path = file.path(project_root, "1-FilteringPlantData", "output",
                            "species_counts4_filtered.rds")) {
  raw <- collect_plant_completeness_raw(counts_path) |>
    dplyr::select(grid_id, plant_completeness = plant_completeness_raw)

  cfr$plant <- cfr$plant |>
    dplyr::select(-dplyr::any_of("plant_completeness")) |>
    dplyr::left_join(raw, by = "grid_id")
  cfr$delta <- cfr$delta |>
    dplyr::select(-dplyr::any_of(c("plant_completeness", "delta_completeness"))) |>
    dplyr::left_join(raw, by = "grid_id") |>
    dplyr::mutate(delta_completeness = plant_completeness - animal_completeness)
  cfr$plant_completeness_source <- "raw"
  cfr
}

# Seven maps: animal/plant completeness, richness, occurrences, plus the
# plant-minus-animal completeness difference.
build_cell_filter_map_panels <- function(cfr) {

  animal_cells <- cfr$animal
  plant_cells  <- cfr$plant
  delta_cells  <- cfr$delta

  outline_scale <- scale_color_manual(
    values = c(`TRUE` = "black", `FALSE` = "grey70"),
    labels = c(`TRUE` = "Retained", `FALSE` = "Filtered out"),
    name = "Animal cell"
  )

  # Retained/filtered outline drawn as a layer BEHIND the fill, so only a thin
  # halo shows around each cell and the fill colour stays readable. Filtered
  # cells are drawn first and retained cells (outline + fill) on top of them.
  outlined_fill <- function(fill_var) {
    cell_layers <- function(keep) {
      subset_cells <- function(d) d[d$retained %in% keep, ]
      list(
        geom_sf(data = subset_cells, aes(color = retained), fill = NA, linewidth = 0.6),
        geom_sf(data = subset_cells, aes(fill = .data[[fill_var]]), color = NA)
      )
    }
    c(cell_layers(c(FALSE, NA)), cell_layers(TRUE))
  }

  map_theme <- theme_minimal(base_size = 9) +
    theme(plot.title = element_text(face = "bold", size = 9),
          legend.key.size = grid::unit(0.4, "cm"))

  map_pair <- function(an_var, pl_var, title, fill_name, trans = "identity") {
    list(
      ggplot(animal_cells) +
        outlined_fill(an_var) +
        scale_fill_viridis_c(name = fill_name, trans = trans) +
        outline_scale +
        labs(title = paste("Animals:", title)) +
        map_theme,
      ggplot(plant_cells) +
        geom_sf(aes(fill = .data[[pl_var]]), color = NA) +
        scale_fill_viridis_c(name = fill_name, trans = trans) +
        labs(title = paste("Plants:", title)) +
        map_theme
    )
  }

  p_delta <- ggplot(delta_cells) +
    outlined_fill("delta_completeness") +
    scale_fill_gradient2(
      name = "Plant - animal\ncompleteness",
      low = "firebrick", mid = "white", high = "steelblue", midpoint = 0,
      limits = c(-1, 1) * max(abs(delta_cells$delta_completeness), na.rm = TRUE)) +
    outline_scale +
    labs(title = "Delta Chao1 completeness (plant - animal), animal grid cells") +
    map_theme

  c(map_pair("animal_completeness", "plant_completeness",
             "Chao1 completeness", "Completeness"),
    map_pair("animal_richness", "plant_richness",
             "species richness", "Richness"),
    map_pair("animal_occurrences", "plant_occurrences",
             "occurrences (log10)", "Occurrences", trans = "log10"),
    list(p_delta))
}

build_completeness_correlation_plot <- function(cfr) {
  d <- sf::st_drop_geometry(cfr$delta) |>
    dplyr::filter(!is.na(plant_completeness), !is.na(animal_completeness)) |>
    # Retained points drawn last (on top)
    dplyr::arrange(retained)
  if (nrow(d) < 3) return(NULL)

  pearson  <- suppressWarnings(cor(d$plant_completeness, d$animal_completeness))
  spearman <- suppressWarnings(cor(d$plant_completeness, d$animal_completeness,
                                   method = "spearman"))

  ggplot(d, aes(x = plant_completeness, y = animal_completeness)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
    geom_smooth(method = "lm", formula = y ~ x, se = TRUE,
                colour = "black", fill = "grey80", linewidth = 0.6) +
    geom_point(aes(colour = retained, alpha = retained), size = 2.5) +
    scale_colour_manual(values = c(`TRUE` = "#2c7fb8", `FALSE` = "grey60"),
                        labels = c(`TRUE` = "Retained", `FALSE` = "Filtered out"),
                        name = "Animal cell") +
    # Same name/labels as the colour scale so the two legends merge
    scale_alpha_manual(values = c(`TRUE` = 0.8, `FALSE` = 0.4),
                       labels = c(`TRUE` = "Retained", `FALSE` = "Filtered out"),
                       name = "Animal cell") +
    annotate("text", x = -Inf, y = Inf, hjust = -0.05, vjust = 1.4, size = 3.5,
             label = sprintf("Pearson r = %.2f\nSpearman rho = %.2f\nn = %d cells",
                             pearson, spearman, nrow(d))) +
    coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
    labs(title = "Plant vs animal completeness in overlapping grid cells",
         subtitle = "Chao1 completeness per cell; dashed line = 1:1",
         x = "Plant completeness", y = "Animal completeness") +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(face = "bold"))
}

# Page 1: dataset summary figure + cell-filter maps.
# Page 2: plant vs animal Chao1 completeness in overlapping grid cells.
# Writes <org_out>/cell_filter_report_<organism>_<timestamp>.pdf (2 pages).
# Returns the file path, or NULL if the inputs are unavailable.
write_cell_filter_report_pdf <- function(organism, org_am, org_out, timestamp,
                                         page_w = 10, page_h = 16) {

  cfr_path <- file.path(org_am, "plot_data", "plotdata_cell_filter_report.rds")
  if (!file.exists(cfr_path)) {
    message("  No plotdata_cell_filter_report.rds for '", organism,
            "'; skipping cell-filter report PDF.")
    return(NULL)
  }
  cfr <- readRDS(cfr_path)
  if (!identical(cfr$plant_completeness_source, "raw")) {
    cfr <- use_raw_plant_completeness(cfr)
  }

  ds_inputs <- load_dataset_summary_figure_inputs(org_am)
  stats_fig <- if (!is.null(ds_inputs)) {
    build_dataset_summary_figure(
      ds_inputs$thresholds, ds_inputs$counts, ds_inputs$dists,
      title = "Dataset before vs after cell filtering", subtitle = NULL)
  }

  maps <- build_cell_filter_map_panels(cfr)
  # 7 maps in a 2-column grid; the delta map is last, leaving one empty slot
  maps_grid <- wrap_plots(maps, ncol = 2) +
    plot_annotation(tag_levels = "A")
  n_map_rows <- ceiling(length(maps) / 2)

  page1 <- if (!is.null(stats_fig)) {
    (wrap_elements(full = stats_fig) / wrap_elements(full = maps_grid)) +
      plot_layout(heights = c(1.5, n_map_rows))
  } else {
    maps_grid
  }
  page1 <- page1 +
    plot_annotation(
      title = paste0(organism_label(organism), ": cell filtering report"),
      theme = theme(plot.title = element_text(size = 16, face = "bold")))

  p_corr <- build_completeness_correlation_plot(cfr)
  page2 <- if (!is.null(p_corr)) {
    (p_corr / plot_spacer()) + plot_layout(heights = c(1, 0.8)) +
      plot_annotation(
        title = paste0(organism_label(organism), ": plant vs animal completeness"),
        theme = theme(plot.title = element_text(size = 16, face = "bold")))
  }

  out_file <- file.path(org_out,
                        paste0("cell_filter_report_", organism, "_", timestamp, ".pdf"))
  pdf(out_file, width = page_w, height = page_h)
  on.exit(dev.off(), add = TRUE)
  print(page1)
  if (!is.null(page2)) {
    print(page2)
  } else {
    message("  Too few overlapping cells for a completeness correlation; page 2 omitted.")
  }
  message("  Wrote: ", out_file)
  out_file
}
