# ─────────────────────────────────────────────────────────────────────────────
# compute_cell_completeness(): abundance-based Chao1 sampling completeness
# per grid cell, calculated from whatever occurrence data is passed in.
#
# Unlike elevation, completeness depends entirely on the occurrence data
# (per-cell record counts and richness), so it must be recalculated whenever
# `occData` changes (e.g. raw vs. filtered vs. interpolated occurrences) —
# it should never be silently reused across different datasets.
#
#   f1 = # species seen once (singletons), f2 = # species seen twice (doubletons)
#   Chao1 = Sobs + f1^2 / (2*f2)         if f2 > 0
#         = Sobs + f1*(f1-1) / 2         if f2 == 0
#   completeness = Sobs / Chao1  (1 = fully inventoried, ->0 = under-sampled)
#
# Returns a tibble with columns grid_id, n_records, completeness.
# ─────────────────────────────────────────────────────────────────────────────
# ─────────────────────────────────────────────────────────────────────────────
# detect_species_col(): identify the column that defines a "species" in an
# occurrence table. Plant data uses `scientificName`; animal data
# (animal_grid) carries `species` / `Morphotype` instead. Auto-detects unless
# the caller specifies one explicitly. `df` should already have geometry
# dropped (a plain data.frame/tibble).
# ─────────────────────────────────────────────────────────────────────────────
detect_species_col <- function(df, species_col = NULL, caller = "detect_species_col()") {
  if (is.null(species_col)) {
    candidates <- c("scientificName", "species", "Morphotype")
    species_col <- candidates[candidates %in% names(df)][1]
    if (is.na(species_col)) {
      stop(caller, ": no species column found. ",
           "Expected one of: scientificName, species, Morphotype. ",
           "Got: ", paste(names(df), collapse = ", "))
    }
  } else if (!species_col %in% names(df)) {
    stop(caller, ": species_col '", species_col, "' not found in occData.")
  }
  species_col
}

compute_cell_completeness <- function(occData, species_col = NULL) {
  tbl <- occData |>
    st_drop_geometry()
  
  species_col <- detect_species_col(tbl, species_col, "compute_cell_completeness()")
  
  # Two input schemas are supported:
  #   * plant data — one row per species-per-cell with an aggregated `count`
  #     column (records for that species in that cell)
  #   * animal data — one row per raw occurrence record and no `count` column
  # For the animal schema we aggregate to species-per-cell record counts first,
  # so the Chao1 f1/f2 (singleton/doubleton species) are computed correctly.
  if (!"count" %in% names(tbl)) {
    tbl <- tbl |>
      group_by(grid_id, .data[[species_col]]) |>
      summarise(count = dplyr::n(), .groups = "drop")
  }
  
  tbl |>
    group_by(grid_id) |>
    summarise(
      n_records = sum(count),
      Sobs      = n_distinct(.data[[species_col]]),
      f1        = sum(count == 1),
      f2        = sum(count == 2),
      .groups   = "drop"
    ) |>
    mutate(
      chao1        = if_else(f2 > 0,
                             Sobs + (f1^2) / (2 * f2),
                             Sobs + (f1 * (f1 - 1)) / 2),
      completeness = if_else(chao1 > 0, Sobs / chao1, NA_real_)
    ) |>
    select(grid_id, n_records, completeness)
}
