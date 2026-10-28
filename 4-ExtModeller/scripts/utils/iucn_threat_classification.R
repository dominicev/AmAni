#' IUCN Threat Classification Based on Population (Range) Reduction
#'
#' Applies IUCN Red List Criterion A thresholds to species range-reduction
#' output from the extinction model. Categories (mutually exclusive):
#'   Critically Endangered / Extinct (CR/EX): >= 90% loss (includes 100% loss)
#'   Endangered (EN):            >= 70% and < 90% loss
#'   Vulnerable (VU):            >= 50% and < 70% loss
#'   Non-threatened (NT/LC):     < 50% loss
#'
#' Data source: each organism's extinction_model/ directory already contains
#' pre-aggregated replicate summaries — summary_range_reduction_thresholds.csv
#' (median species counts crossing 20/50/75/90/95/99/100% loss, per model
#' [Clustered/Shifted/Neutral] x per scenario [FreqOver/FreqBoundary/
#' FreqCalibrated/FreqUnder]) and summary_range_reduction_single_run.csv
#' (total species modeled, n_species). These exist consistently for all three
#' simulation models, whereas exact per-species pct_range_lost is only stored
#' for the clustered and neutral models.
#'
#' No exact 70% threshold is pre-computed for any model, so the EN/VU
#' boundary is estimated by linear interpolation between the 50% and 75%
#' replicate-median counts: n_reduced_70 ~= n50 + 0.8 * (n75 - n50).

library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)

#' Compute IUCN category counts for each model x scenario
#'
#' @param extinction_dir Path to an organism's extinction_model/ directory
#' @return A list with:
#'   - counts: long-format tibble (model, scenario, iucn_category, count)
#'   - total_species: total number of species modeled
compute_iucn_category_counts <- function(extinction_dir) {

  thresh_path      <- file.path(extinction_dir, "summary_range_reduction_thresholds.csv")
  single_run_path  <- file.path(extinction_dir, "summary_range_reduction_single_run.csv")

  if (!file.exists(thresh_path)) {
    stop("Missing summary_range_reduction_thresholds.csv in: ", extinction_dir)
  }

  thresh <- read.csv(thresh_path, stringsAsFactors = FALSE)

  total_species <- NA_integer_
  if (file.exists(single_run_path)) {
    single <- read.csv(single_run_path, stringsAsFactors = FALSE)
    if ("n_species" %in% names(single) && nrow(single) > 0) {
      total_species <- unique(single$n_species)[1]
    }
  }
  if (is.na(total_species)) {
    warning("Could not determine total_species from summary_range_reduction_single_run.csv; ",
            "NT/LC counts will be NA.")
  }

  category_levels <- c("NT/LC", "VU", "EN", "CR/EX")

  counts <- thresh |>
    mutate(
      # Linear interpolation for the 70% (EN/VU) cutoff between the 50% and
      # 75% replicate-median thresholds (no exact 70% value is precomputed).
      n_reduced_70_est = n_reduced_50_median + 0.8 * (n_reduced_75_median - n_reduced_50_median),
      # CR and fully extinct (100% loss) species are combined into one column.
      `CR/EX` = n_reduced_90_median,
      EN      = pmax(n_reduced_70_est     - n_reduced_90_median, 0),
      VU      = pmax(n_reduced_50_median  - n_reduced_70_est,     0),
      `NT/LC` = pmax(total_species - n_reduced_50_median,          0)
    ) |>
    select(model, scenario, `NT/LC`, VU, EN, `CR/EX`) |>
    pivot_longer(
      cols = c(`NT/LC`, VU, EN, `CR/EX`),
      names_to = "iucn_category",
      values_to = "count"
    ) |>
    mutate(iucn_category = factor(iucn_category, levels = category_levels))

  list(counts = counts, total_species = total_species)
}


#' Plot IUCN threat classification across models, with uncertainty ribbons
#'
#' Draws, for each simulation model (Clustered/Shifted/Neutral):
#'   - a semi-transparent grey ribbon spanning the FreqUnder-FreqOver range
#'     at each threat category (ribbons from different models overlay, so
#'     areas where models agree appear darker)
#'   - short horizontal lines marking the FreqBoundary (solid) and
#'     FreqCalibrated (dashed) best estimates, labeled with species counts
#' The y-axis is log-scaled but labeled with actual species counts.
#'
#' @param classification Output from compute_iucn_category_counts()
#' @param output_dir Directory to save PDF and RDS files
#' @param filename_prefix Prefix for output files (default: "iucn_threat")
#' @param organism_label Optional label for the plot title
#' @return The ggplot object (invisibly)
plot_iucn_threat <- function(classification,
                             output_dir,
                             filename_prefix = "iucn_threat",
                             organism_label = NULL) {

  counts <- classification$counts
  total_species <- classification$total_species

  category_levels <- c("NT/LC", "VU", "EN", "CR/EX")
  x_pos <- setNames(seq_along(category_levels), category_levels)

  counts <- counts |>
    mutate(x = x_pos[as.character(iucn_category)])

  # Floor for the log scale; the real (possibly zero) value is kept for labels
  log_floor <- 0.5

  # Ribbons: FreqUnder-FreqOver span per model, drawn as a rectangle centered
  # on each category (NOT connected across categories) so that models'
  # rectangles overlay directly on top of one another within a category.
  ribbon_df <- counts |>
    filter(scenario %in% c("FreqUnder", "FreqOver")) |>
    select(model, x, scenario, count) |>
    pivot_wider(names_from = scenario, values_from = count) |>
    rename(ymin_raw = FreqUnder, ymax_raw = FreqOver) |>
    mutate(
      ymin = pmax(pmin(ymin_raw, ymax_raw), log_floor),
      ymax = pmax(pmax(ymin_raw, ymax_raw), log_floor),
      xmin = x - 0.32,
      xmax = x + 0.32
    )

  # Best-estimate horizontal lines: FreqBoundary (solid) and FreqCalibrated (dashed).
  # Lines/labels are nudged slightly per model so the three models' ticks
  # don't sit exactly on top of one another and labels stay legible.
  model_nudge <- c("Clustered" = -0.12, "Neutral" = 0, "Shifted" = 0.12)

  best_df <- counts |>
    filter(scenario %in% c("FreqBoundary", "FreqCalibrated")) |>
    mutate(
      count_plot = pmax(count, log_floor),
      x_nudged = x + model_nudge[model],
      xmin = x_nudged - 0.14,
      xmax = x_nudged + 0.14,
      scenario = factor(scenario, levels = c("FreqBoundary", "FreqCalibrated"))
    )

  model_colors <- c(
    "Clustered" = "#1B6CA8",
    "Shifted"   = "#7A4FA3",
    "Neutral"   = "#2E9E5B"
  )

  subtitle_txt <- if (!is.na(total_species)) {
    paste0("Total species modeled: ", scales::comma(total_species))
  } else {
    NULL
  }

  # Only label the FreqBoundary line per model to keep the plot legible;
  # FreqCalibrated remains visible via its dashed linetype.
  label_df <- best_df |> filter(scenario == "FreqBoundary")

  p <- ggplot() +
    geom_rect(
      data = ribbon_df,
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, group = model),
      fill = "grey20", alpha = 0.30, colour = NA
    ) +
    geom_segment(
      data = best_df,
      aes(x = xmin, xend = xmax, y = count_plot, yend = count_plot,
          colour = model, linetype = scenario),
      linewidth = 1
    ) +
    geom_text(
      data = label_df,
      aes(x = x_nudged, y = count_plot,
          label = scales::comma(round(count)), colour = model),
      vjust = -0.6, size = 2.6, show.legend = FALSE
    ) +
    scale_x_continuous(
      breaks = x_pos, labels = names(x_pos),
      expand = expansion(add = 0.9)
    ) +
    scale_y_log10(labels = scales::comma) +
    scale_colour_manual(values = model_colors) +
    scale_linetype_manual(values = c(FreqBoundary = "solid", FreqCalibrated = "dashed")) +
    labs(
      title = organism_label,
      subtitle = subtitle_txt,
      x = "IUCN Threat Category",
      y = "Number of species (log scale)",
      colour = "Model",
      linetype = "Best estimate",
      caption = paste(
        "Grey bands: FreqUnder-FreqOver range per model (darker overlap = higher agreement across models).",
        "Lines: FreqBoundary (solid) / FreqCalibrated (dashed) best estimates. Labels show species counts.",
        sep = "\n"
      )
    ) +
    theme_minimal(base_size = 13) +
    theme(plot.caption = element_text(size = 8, hjust = 0))

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  pdf_path <- file.path(output_dir, paste0(filename_prefix, ".pdf"))
  rds_path <- file.path(output_dir, paste0(filename_prefix, "_data.rds"))

  ggsave(pdf_path, plot = p, width = 8, height = 6)
  saveRDS(list(counts = counts, ribbon_df = ribbon_df, best_df = best_df,
               total_species = total_species), rds_path)

  message("Saved: ", pdf_path)
  message("Saved: ", rds_path)

  invisible(p)
}
