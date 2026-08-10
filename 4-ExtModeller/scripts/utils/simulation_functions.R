{
  library(sf)
  library(parallel)
  library(doParallel)
  library(foreach)
}


# =====================================================================
# TOROIDAL SHIFT — spatial null model for grid weights
# =====================================================================
# Shifts weight values on a (nrow × ncol) grid by a random offset,
# wrapping around edges. Preserves all local spatial autocorrelation
# but disassociates the pattern from its real-world geographic location.
#
# grid_weights: named numeric vector (names = grid_id as character)
# nrow_grid, ncol_grid: grid dimensions
# grid_ids: vector of grid_ids in the order produced by st_make_grid
#           (column-major: column 1 bottom→top, then column 2, etc.)
# occupied_ids: grid_ids that contain species (weights on unoccupied
#               grids are zeroed out after shifting)
# shift_row, shift_col: offset amounts. If NULL, sampled randomly.

toroidal_shift_weights <- function(
    grid_weights,
    nrow_grid,
    ncol_grid,
    grid_ids,
    occupied_ids = NULL,
    shift_row = NULL,
    shift_col = NULL
) {
  n <- nrow_grid * ncol_grid
  stopifnot(length(grid_ids) == n)

  # Random shift in each dimension independently
  if (is.null(shift_row)) shift_row <- sample.int(nrow_grid, 1)
  if (is.null(shift_col)) shift_col <- sample.int(ncol_grid, 1)

  # 1. Rectangular weight matrix (column-major, 0-padded for any missing cells)
  # st_make_grid is column-major so grid_ids go down rows then across cols
  orig_weights <- grid_weights[as.character(grid_ids)]
  orig_weights[is.na(orig_weights)] <- 0
  wmat <- matrix(orig_weights, nrow = nrow_grid, ncol = ncol_grid)

  # 2. Toroidal shift — wraps within the original rectangular dimensions
  row_order <- ((seq_len(nrow_grid) - 1 + shift_row) %% nrow_grid) + 1
  col_order <- ((seq_len(ncol_grid) - 1 + shift_col) %% ncol_grid) + 1
  wmat <- wmat[row_order, col_order, drop = FALSE]

  # 3. Convert back to named vector, mask non-occupied cells to 0
  result_weights <- setNames(as.numeric(wmat), as.character(grid_ids))

  if (!is.null(occupied_ids)) {
    occupied_mask <- as.character(grid_ids) %in% as.character(occupied_ids)
    result_weights[!occupied_mask] <- 0
  }

  result_weights
}


# =====================================================================
# RANGE REDUCTION HELPERS (base-R only — safe for parallel workers)
# =====================================================================

# Compute per-species range reduction by comparing surviving grid cells
# to a baseline. Returns a data.frame with one row per species.
#
# baseline_ranges: data.frame with columns scientificName, n_grids_before
#   (precomputed once per dataset to avoid recomputation in each replicate)
# final_table: the surviving occurrence table after deforestation

compute_range_reduction <- function(final_table, baseline_ranges) {
  if (nrow(final_table) > 0) {
    pairs <- unique(final_table[, c("scientificName", "grid_id"), drop = FALSE])
    agg   <- stats::aggregate(grid_id ~ scientificName, data = pairs, FUN = length)
    names(agg)[2] <- "n_grids_after"
  } else {
    agg <- data.frame(scientificName = character(0),
                      n_grids_after  = integer(0),
                      stringsAsFactors = FALSE)
  }

  rr <- merge(baseline_ranges, agg, by = "scientificName", all.x = TRUE)
  rr$n_grids_after[is.na(rr$n_grids_after)] <- 0L
  rr$pct_range_lost <- (rr$n_grids_before - rr$n_grids_after) / rr$n_grids_before
  rr
}

# Count species meeting each range-reduction threshold.
# Returns a named integer vector (n_reduced_20, ..., n_reduced_100).
RANGE_REDUCTION_THRESHOLDS <- c(0.20, 0.50, 0.75, 0.90, 0.95, 0.99, 1.00)
RANGE_REDUCTION_NAMES      <- paste0("n_reduced_", c("20","50","75","90","95","99","100"))

count_range_thresholds <- function(species_rr) {
  counts <- vapply(
    RANGE_REDUCTION_THRESHOLDS,
    function(t) sum(species_rr$pct_range_lost >= t, na.rm = TRUE),
    integer(1)
  )
  setNames(counts, RANGE_REDUCTION_NAMES)
}

# Precompute baseline range sizes (unique grid cells per species) from the
# full occurrence table. Call once before the parallel loop.
precompute_baseline_ranges <- function(sim_tbl) {
  pairs <- unique(sim_tbl[, c("scientificName", "grid_id"), drop = FALSE])
  agg   <- stats::aggregate(grid_id ~ scientificName, data = pairs, FUN = length)
  names(agg)[2] <- "n_grids_before"
  agg
}


# =====================================================================
# SMALL-PATCH THINNING (vectorized)
# =====================================================================
# For each small-patch bin, sample `count` random grid cells and remove
# a fraction (bin_mid / grid_cell_size_km2) of rows from each hit cell.
# Grids hit multiple times across bins accumulate removal probability:
#   p_survive = product of (1 - perc) for each hit
# All removals are computed up front; data.frame is subset ONCE.
#
# grid_weights: optional named numeric vector of sampling weights per
# grid_id. When NULL, sampling is uniform (neutral model). When provided,
# grids are sampled proportional to their weight (clustered model).

run_small_patch_thinning <- function(
    freq_small,
    current_tbl,
    original_tbl,
    grid_cell_size_km2,
    grid_weights = NULL,
    seed = NULL
) {
  if (!is.null(seed)) set.seed(seed)

  freq_small <- as.data.frame(freq_small)
  unique_grids <- unique(current_tbl$grid_id)
  unique_grids <- unique_grids[!is.na(unique_grids)]
  n_grids <- length(unique_grids)

  if (n_grids == 0) {
    return(list(
      final_table = current_tbl,
      n_pop_lost  = nrow(original_tbl) - nrow(current_tbl),
      n_spp_lost  = length(unique(original_tbl$scientificName)) -
                    length(unique(current_tbl$scientificName))
    ))
  }

  # resolve sampling weights for available grids
  if (!is.null(grid_weights)) {
    w <- grid_weights[as.character(unique_grids)]
    w[is.na(w)] <- 0
    # normalize so they sum to 1
    if (sum(w) > 0) w <- w / sum(w) else w <- rep(1 / n_grids, n_grids)
  } else {
    w <- NULL
  }

  # accumulate survival probability per grid cell (start at 1.0)
  grid_survive_prob <- setNames(rep(1.0, n_grids), unique_grids)

  for (i in seq_len(nrow(freq_small))) {
    perc <- min(freq_small$bin_mid[i] / grid_cell_size_km2, 1)
    n_hits <- freq_small$count[i]

    hit_grids <- sample(unique_grids, size = n_hits, replace = TRUE, prob = w)
    hit_counts <- table(hit_grids)

    for (g in names(hit_counts)) {
      grid_survive_prob[g] <- grid_survive_prob[g] * (1 - perc)^hit_counts[g]
    }
  }

  row_grid_ids <- as.character(current_tbl$grid_id)
  row_survive <- grid_survive_prob[row_grid_ids]

  keep <- runif(nrow(current_tbl)) < row_survive
  sim_tbl_final <- current_tbl[keep, , drop = FALSE]

  n_pop_lost <- nrow(original_tbl) - nrow(sim_tbl_final)
  n_spp_lost <- length(unique(original_tbl$scientificName)) -
                length(unique(sim_tbl_final$scientificName))

  list(
    final_table = sim_tbl_final,
    n_pop_lost  = n_pop_lost,
    n_spp_lost  = n_spp_lost
  )
}


# =====================================================================
# LARGE-PATCH DEFORESTATION (vectorized)
# =====================================================================
# grid_weights: optional named numeric vector. When provided, grid cells
# are sampled proportional to their weight rather than uniformly.

run_large_patch_deforestation <- function(
    freq_large,
    sim_tbl,
    grid_cell_size_km2,
    grid_weights = NULL,
    seed = NULL
) {
  if (!is.null(seed)) set.seed(seed)

  freq_large <- as.data.frame(freq_large)

  grid_equiv <- freq_large$count * (freq_large$bin_mid / grid_cell_size_km2)
  n_whole    <- floor(grid_equiv)
  remainder_area <- (grid_equiv - n_whole) * grid_cell_size_km2

  total_grids_to_remove <- sum(n_whole)

  remainderPatches <- data.frame(
    extraPatchSize = remainder_area,
    count          = ifelse(remainder_area > 0, 1L, 0L)
  )

  available_grids <- unique(sim_tbl$grid_id)
  available_grids <- available_grids[!is.na(available_grids)]

  total_grids_to_remove <- min(total_grids_to_remove, length(available_grids))

  # resolve sampling weights
  if (!is.null(grid_weights)) {
    w <- grid_weights[as.character(available_grids)]
    w[is.na(w)] <- 0
    if (sum(w) > 0) w <- w / sum(w) else w <- rep(1 / length(available_grids), length(available_grids))
  } else {
    w <- NULL
  }

  # sample grid cells to deforest (weighted or uniform)
  # note: replace=FALSE with prob requires that we don't ask for more
  # than the number of non-zero-weight grids
  if (!is.null(w)) {
    n_eligible <- sum(w > 0)
    total_grids_to_remove <- min(total_grids_to_remove, n_eligible)
  }
  grids_to_deforest <- sample(
    available_grids, total_grids_to_remove, replace = FALSE, prob = w
  )
  deforestedLarge <- sim_tbl[!(sim_tbl$grid_id %in% grids_to_deforest), , drop = FALSE]

  n_pop_lost <- nrow(sim_tbl) - nrow(deforestedLarge)
  n_spp_lost <- length(unique(sim_tbl$scientificName)) -
                length(unique(deforestedLarge$scientificName))

  list(
    final_table      = deforestedLarge,
    remainderPatches = remainderPatches,
    n_pop_lost       = n_pop_lost,
    n_spp_lost       = n_spp_lost
  )
}


# =====================================================================
# MAIN WRAPPER: simulate one full deforestation–extinction event
# =====================================================================
# grid_weights: passed through to large- and small-patch functions.
# NULL = uniform/neutral; named vector = clustered.

simulateExtinctionFromDeforestation <- function(
    sim_tbl,
    freq_table,
    grid_cell_size_km2,
    grid_weights = NULL,
    seed_large = NULL,
    seed_small = NULL,
    baseline_ranges = NULL
) {
  if (inherits(sim_tbl, "sf")) sim_tbl <- sf::st_drop_geometry(sim_tbl)

  # --- large patches ---
  freq_large <- freq_table[
    freq_table$bin_mid >= grid_cell_size_km2 & freq_table$count > 0, ,
    drop = FALSE
  ]

  out_large <- run_large_patch_deforestation(
    freq_large         = freq_large,
    sim_tbl            = sim_tbl,
    grid_cell_size_km2 = grid_cell_size_km2,
    grid_weights       = grid_weights,
    seed               = seed_large
  )

  deforestedLarge  <- out_large$final_table
  remainderPatches <- out_large$remainderPatches

  # --- fold remainder area into frequency table ---
  freq_mod <- freq_table

  for (i in seq_len(nrow(remainderPatches))) {
    patch_size  <- remainderPatches$extraPatchSize[i]
    patch_count <- remainderPatches$count[i]
    if (is.na(patch_size) || patch_count < 1) next

    bin_idx <- which(
      freq_mod$bin_lower <= patch_size & freq_mod$bin_upper > patch_size
    )
    if (length(bin_idx) == 1) {
      freq_mod$count[bin_idx] <- freq_mod$count[bin_idx] + patch_count
    } else {
      warning(
        "Remainder patch (", round(patch_size, 2), " km2) ",
        "matched ", length(bin_idx), " bins; area lost from simulation."
      )
    }
  }

  # --- small patches ---
  freq_small <- freq_mod[
    freq_mod$bin_mid < grid_cell_size_km2 & freq_mod$count > 0, ,
    drop = FALSE
  ]

  out_small <- run_small_patch_thinning(
    freq_small         = freq_small,
    current_tbl        = deforestedLarge,
    original_tbl       = sim_tbl,
    grid_cell_size_km2 = grid_cell_size_km2,
    grid_weights       = grid_weights,
    seed               = seed_small
  )

  # --- Per-species range reduction ---
  final_tbl <- out_small$final_table

  if (is.null(baseline_ranges)) {
    baseline_ranges <- precompute_baseline_ranges(sim_tbl)
  }

  species_range_reduction <- compute_range_reduction(final_tbl, baseline_ranges)
  threshold_counts        <- count_range_thresholds(species_range_reduction)

  list(
    final_table             = final_tbl,
    freq_table_original     = freq_table,
    freq_table_modified     = freq_mod,
    remainderPatches        = remainderPatches,
    n_pop_lost_large        = out_large$n_pop_lost,
    n_spp_lost_large        = out_large$n_spp_lost,
    n_pop_lost_total        = out_small$n_pop_lost,
    n_spp_lost_total        = out_small$n_spp_lost,
    species_range_reduction = species_range_reduction,
    threshold_counts        = threshold_counts
  )
}


# =====================================================================
# PARALLEL REPLICATE RUNNER with file-based progress tracking
# =====================================================================
# grid_weights: passed through to simulateExtinctionFromDeforestation.
#   NULL = neutral model; named vector = clustered model.
#
# toroidal_params: if provided, a list with elements:
#   nrow_grid, ncol_grid, grid_ids, occupied_ids
# When set, each replicate applies a DIFFERENT random toroidal shift
# to grid_weights before running. This is the "shifted" null model.

SEFD_rep_parallel <- function(
    simAnimalObj,
    freq_table,
    grid_cell_size_km2,
    n_runs = 1,
    grid_weights = NULL,
    toroidal_params = NULL,
    progress_dir = NULL
) {
  sim_tbl <- sf::st_drop_geometry(simAnimalObj)

  # Precompute baseline ranges once (shared across all replicates)
  baseline_ranges <- precompute_baseline_ranges(sim_tbl)

  if (is.null(progress_dir)) {
    progress_dir <- file.path(tempdir(), "sefd_progress")
  }
  if (dir.exists(progress_dir)) unlink(progress_dir, recursive = TRUE)
  dir.create(progress_dir, recursive = TRUE, showWarnings = FALSE)
  message("Progress tracking: ", progress_dir)
  message("Monitor with:  length(list.files('", progress_dir, "'))")

  n_cores_total <- parallel::detectCores(logical = FALSE)
  n_cores_use <- max(1, floor(n_cores_total * 0.7))
  message("Using ", n_cores_use, " of ", n_cores_total, " cores for ", n_runs, " replicates")

  cl <- parallel::makeCluster(n_cores_use)
  on.exit(parallel::stopCluster(cl), add = TRUE)

  doParallel::registerDoParallel(cl)

  # determine which functions/objects to export
  export_fns <- c(
    "simulateExtinctionFromDeforestation",
    "run_large_patch_deforestation",
    "run_small_patch_thinning",
    "compute_range_reduction",
    "count_range_thresholds",
    "precompute_baseline_ranges",
    "RANGE_REDUCTION_THRESHOLDS",
    "RANGE_REDUCTION_NAMES",
    "baseline_ranges"
  )
  if (!is.null(toroidal_params)) {
    export_fns <- c(export_fns, "toroidal_shift_weights")
  }

  extinction_results <- foreach(
    i = seq_len(n_runs),
    .combine   = rbind,
    .export    = export_fns
  ) %dopar% {

    # if toroidal shift requested, apply a random shift per replicate
    weights_i <- grid_weights
    if (!is.null(toroidal_params) && !is.null(grid_weights)) {
      weights_i <- toroidal_shift_weights(
        grid_weights = grid_weights,
        nrow_grid    = toroidal_params$nrow_grid,
        ncol_grid    = toroidal_params$ncol_grid,
        grid_ids     = toroidal_params$grid_ids,
        occupied_ids = toroidal_params$occupied_ids
      )
    }

    out_all <- simulateExtinctionFromDeforestation(
      sim_tbl            = sim_tbl,
      freq_table         = freq_table,
      grid_cell_size_km2 = grid_cell_size_km2,
      grid_weights       = weights_i,
      baseline_ranges    = baseline_ranges
    )

    writeLines(
      as.character(Sys.time()),
      file.path(progress_dir, paste0("run_", i, ".done"))
    )

    tc <- out_all$threshold_counts
    data.frame(
      run_id           = i,
      n_pop_lost_total = out_all$n_pop_lost_total,
      n_spp_lost_total = out_all$n_spp_lost_total,
      n_reduced_20     = unname(tc["n_reduced_20"]),
      n_reduced_50     = unname(tc["n_reduced_50"]),
      n_reduced_75     = unname(tc["n_reduced_75"]),
      n_reduced_90     = unname(tc["n_reduced_90"]),
      n_reduced_95     = unname(tc["n_reduced_95"]),
      n_reduced_99     = unname(tc["n_reduced_99"]),
      n_reduced_100    = unname(tc["n_reduced_100"])
    )
  }

  message("All ", n_runs, " replicates complete.")
  return(extinction_results)
}
