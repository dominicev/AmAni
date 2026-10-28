# tmp_dataset_summary_figure.R
# Compilation figure summarizing a dataset before vs after cell filtering.
# Source data: output/odonata/animal_model/dataset_summary_*.txt
# Temporary/scratch script.

library(tidyverse)
library(patchwork)

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
#' @return A patchwork object (also printable/ggsave-able).
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

  counts <- counts |>
    mutate(
      pct_retained = after / before,
      metric = factor(metric, levels = unique(metric))
    )

  dists <- dists |>
    mutate(
      stage = factor(stage, levels = c("Before", "After")),
      variable = factor(variable, levels = unique(variable))
    )

  # --- Panel A: filter thresholds (text card) -------------------------------
  p_thresh <- thresholds |>
    mutate(y = rev(seq_len(n()))) |>
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
    d <- filter(dists, variable == var)
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
  boxes <- wrap_plots(box_panels, ncol = 2)

  fig <- top / boxes +
    plot_layout(heights = c(1, ceiling(length(box_panels) / 2))) +
    plot_annotation(
      title = title,
      subtitle = subtitle,
      theme = theme(
        plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(size = 9, colour = "grey30")
      )
    )

  fig
}

# --- Example call with values from dataset_summary_20260902_123826.txt -------
if (identical(environment(), globalenv())) {

  thresholds <- tribble(
    ~label,             ~value,
    "min occurrences",  "4",
    "min completeness", "0.6",
    "min richness",     "20",
    "outlier cutoff",   "135 spp",
    "Track B entry",    "97th pctile (53 spp)"
  )

  counts <- tribble(
    ~metric,      ~before, ~after,
    "Grid cells", 989,     20,
    "Records",    9815,    2342,
    "Species",    674,     303
  )

  dists <- tribble(
    ~variable,      ~stage,   ~min, ~q1,   ~med,  ~q3,  ~max,  ~mean,
    "Richness",     "Before", 1,    1,     2,     7,    132,   6.3,
    "Richness",     "After",  21,   26.5,  33.5,  55,   132,   45.2,
    "Occurrences",  "Before", 1,    1,     2,     8,    524,   9.9,
    "Occurrences",  "After",  40,   51.5,  67,    122,  524,   117.1,
    "Completeness", "Before", 0.1,  0.4,   0.7,   1,    1,     0.7,
    "Completeness", "After",  0.6,  0.6,   0.7,   0.7,  0.9,   0.7
  )

  fig <- build_dataset_summary_figure(thresholds, counts, dists)
  print(fig)

  #out_pdf <- "/home/dominice/Desktop/AmAni/output/odonata/animal_model/dataset_summary_figure.pdf"
  #ggsave(out_pdf, fig, width = 10, height = 10)
  message("Saved: ", out_pdf)
}
