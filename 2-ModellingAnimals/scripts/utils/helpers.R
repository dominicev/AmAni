# helpers.R — Small utility functions used across Stage 2 scripts

rename_grid_id <- function(df) {
  if ("grid_id.x" %in% colnames(df)) {
    df <- dplyr::rename(df, grid_id = grid_id.x)
  }
  # Drop duplicate grid_id columns left by st_join
  if ("grid_id.y" %in% colnames(df)) {
    df <- dplyr::select(df, -grid_id.y)
  }
  df
}


# ── Data-driven filter threshold estimation ──────────────────────────
#
# estimate_filter_thresholds()
#
# Computes empirically-motivated values for min_occurrences and
# min_completeness from the per-cell Chao1 statistics. These can be
# used as defaults or as a starting point for manual tuning.
#
# min_richness is NOT estimated here — it should be set manually based
# on domain knowledge (a conservative estimate of the minimum species
# count one would expect in the least diverse region of the Amazon,
# including deforested areas).
#
# Args:
#   cell_stats  A data frame with columns: grid_id, n_occurrences,
#               n_species, f1, f2, richness_est, completeness (as produced
#               by the cell filter block in ModelAnimal.R).
#   discovery_rate_pct  The discovery rate must drop to this fraction
#               of the baseline rate (measured at 5–15 occurrences)
#               for a cell to be considered adequately sampled. Default
#               0.40 (rate has dropped to 40% of baseline). This scales
#               automatically across taxa with different diversity levels.
#   baseline_occ_range  Occurrence range used to estimate the baseline
#               discovery rate. Default c(5, 15).
#   completeness_gap_quantile  The quantile of the "relative gap"
#               distribution ((Chao1 - S_obs) / S_obs) that must
#               fall below completeness_gap_target. Default 0.50
#               (median).
#   completeness_gap_target  The maximum acceptable relative gap at
#               the chosen quantile. Default 0.50 (i.e., Chao1
#               estimates at most 50% more species than observed).
#
# Returns:
#   A named list with elements:
#     min_occurrences   integer
#     min_completeness  numeric (rounded to 2 decimal places)

estimate_filter_thresholds <- function(
    cell_stats,
    discovery_rate_pct = 0.40,
    baseline_occ_range = c(5, 15),
    completeness_gap_quantile = 0.50,
    completeness_gap_target = 0.50
) {

  # --- min_occurrences: relative discovery rate curve ---
  # The discovery rate (species / occurrences) is computed for each cell.
  # We measure a baseline rate at low occurrence counts, then find the
  # occurrence count where the rate has dropped to discovery_rate_pct of
  # that baseline. This normalises across taxa: hyper-diverse groups need
  # more occurrences, but the *relative* saturation is comparable.
  rate_data <- cell_stats |>
    dplyr::filter(n_occurrences >= 3) |>
    dplyr::mutate(
      log_occ = log10(n_occurrences),
      discovery_rate = n_species / n_occurrences
    )

  # Baseline rate from cells in the reference occurrence range
  baseline_cells <- rate_data |>
    dplyr::filter(n_occurrences >= baseline_occ_range[1],
                  n_occurrences <= baseline_occ_range[2])

  if (nrow(baseline_cells) < 10) {
    warning("Fewer than 10 cells in the baseline occurrence range (",
            baseline_occ_range[1], "–", baseline_occ_range[2],
            "). Estimate may be unreliable.", call. = FALSE)
  }

  baseline_rate <- median(baseline_cells$discovery_rate)
  target_rate <- baseline_rate * discovery_rate_pct

  # Fit loess with fallback for sparse data (too few unique x values or
  # near-singular neighbourhoods cause predLoess to error)
  est_min_occ <- tryCatch({
    rate_fit <- stats::loess(discovery_rate ~ log_occ, data = rate_data, span = 0.3)
    occ_grid <- seq(min(rate_data$log_occ), max(rate_data$log_occ), length.out = 500)
    pred_rate <- stats::predict(rate_fit, newdata = data.frame(log_occ = occ_grid))

    crossing_idx <- which(pred_rate <= target_rate)[1]
    if (!is.na(crossing_idx)) {
      ceiling(10^occ_grid[crossing_idx])
    } else {
      warning("Discovery rate never drops to ", discovery_rate_pct * 100,
              "% of baseline (", round(target_rate, 3),
              "). Returning NA for min_occurrences.", call. = FALSE)
      NA_integer_
    }
  }, error = function(e) {
    # Loess fails with very sparse data — fall back to a simple quantile
    # approach: use the 75th percentile of occurrences as a rough threshold
    fallback <- ceiling(quantile(cell_stats$n_occurrences, 0.75))
    warning("Loess fit failed (", conditionMessage(e),
            "). Falling back to 75th percentile of occurrences: ", fallback,
            call. = FALSE)
    as.integer(fallback)
  })

  # --- min_completeness: relative gap criterion ---
  # Use cells with >= est_min_occ (or >= 30 as fallback) so the completeness
  # distribution isn't dominated by very sparse cells.
  occ_floor <- if (!is.na(est_min_occ)) min(est_min_occ, 30) else 30
  comp_data <- cell_stats |>
    dplyr::filter(n_occurrences >= occ_floor, !is.na(richness_est)) |>
    dplyr::mutate(relative_gap = (richness_est - n_species) / n_species)

  # Sweep completeness values using a narrow window (±0.05) and find the
  # lowest completeness where the chosen quantile of the gap drops below
  # the target.
  comp_sweep <- seq(0.10, 0.95, by = 0.01)
  window <- 0.05
  est_min_comp <- NA_real_
  for (c_val in comp_sweep) {
    cells_in_window <- comp_data |>
      dplyr::filter(completeness >= c_val - window, completeness < c_val + window)
    if (nrow(cells_in_window) < 10) next
    gap_q <- quantile(cells_in_window$relative_gap, completeness_gap_quantile)
    if (gap_q <= completeness_gap_target) {
      est_min_comp <- c_val
      break
    }
  }

  if (is.na(est_min_comp)) {
    warning("Could not find a completeness threshold where the ",
            completeness_gap_quantile * 100, "th percentile of the missing ",
            "fraction drops below ", completeness_gap_target * 100, "%.",
            call. = FALSE)
  }

  result <- list(
    min_occurrences = est_min_occ,
    min_completeness = round(est_min_comp, 2)
  )

  message("  Estimated filter thresholds:")
  message("    min_occurrences  = ", result$min_occurrences,
          " (baseline rate = ", round(baseline_rate, 3),
          ", target = ", round(target_rate, 3),
          " [", discovery_rate_pct * 100, "% of baseline])")
  message("    min_completeness = ", result$min_completeness,
          " (", completeness_gap_quantile * 100, "th pctile of missing fraction < ",
          completeness_gap_target * 100, "%)")

  result
}
