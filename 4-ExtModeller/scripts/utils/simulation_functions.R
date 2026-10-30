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



# DESIGN (data-parallel over rows, replicate loop inside each worker):
#   * The occurrence table is converted once into integer indices
#     (species, grid, species-grid pair), sorted so that every pair's rows
#     are contiguous, and split into one contiguous slice per worker.
#   * Each PSOCK worker receives ONLY its slice (a few MB), never the full
#     table. The same slice is reused for all replicates.
#   * The parent generates the deforestation pattern for each replicate
#     (a survival probability per grid; tiny), broadcasts batches of
#     patterns, and sums the per-species counts the workers return.
#   PSOCK (not fork) is used deliberately: forking the RStudio session
#   can freeze/kill it, and GC in forked children duplicates the parent's
#   memory. Worker memory is bounded by slice size, not table size.


# --- Integer index built once per call ---
# Rows with NA grid_id never survive. Rows with NA grid_id or NA
# scientificName do not count toward range sizes (matching aggregate()'s NA
# handling in precompute_baseline_ranges). NA scientificName is still a
# "species" for the n_spp_lost count (matching unique()).
prepare_sefd_index <- function(sim_tbl) {
  unique_grids <- unique(sim_tbl$grid_id)
  unique_grids <- unique_grids[!is.na(unique_grids)]
  n_grids <- length(unique_grids)

  sp_levels <- unique(sim_tbl$scientificName)
  n_sp      <- length(sp_levels)
  sp_idx    <- match(sim_tbl$scientificName, sp_levels)
  grid_idx  <- match(sim_tbl$grid_id, unique_grids)

  valid <- !is.na(grid_idx) & !is.na(sim_tbl$scientificName)
  # double, not integer, so n_species * n_grids cannot overflow
  key <- (as.numeric(sp_idx) - 1) * n_grids + grid_idx
  key[!valid] <- NA

  # Sort so rows of the same (species, grid) pair are contiguous; NA keys last
  ord      <- order(key, na.last = TRUE)
  key      <- key[ord]
  sp_idx   <- sp_idx[ord]
  grid_idx <- grid_idx[ord]
  valid    <- valid[ord]
  n <- length(key)

  n_valid <- sum(valid)
  if (n_valid > 0) {
    kv <- key[seq_len(n_valid)]
    new_pair <- c(TRUE, kv[-1] != kv[-n_valid])
    pair_id  <- c(cumsum(new_pair), rep(NA_integer_, n - n_valid))
    pair_sp  <- sp_idx[seq_len(n_valid)][new_pair]
  } else {
    pair_id <- rep(NA_integer_, n)
    pair_sp <- integer(0)
  }

  list(
    unique_grids = unique_grids,
    n_grids      = n_grids,
    n_rows       = n,
    n_sp         = n_sp,
    sp_idx       = sp_idx,
    grid_idx     = grid_idx,
    pair_id      = pair_id,   # 1..n_pairs, non-decreasing; NA for invalid rows
    pair_sp      = pair_sp,   # species of each pair
    n_pairs      = length(pair_sp),
    n_before     = tabulate(pair_sp, n_sp)
  )
}


# --- Split the index into contiguous row slices, never splitting a pair ---
split_sefd_index <- function(idx, n_slices) {
  n <- idx$n_rows
  n_slices <- max(1, min(n_slices, n))
  cut_pts <- floor(seq(0, n, length.out = n_slices + 1))
  # move each internal cut forward to the start of the next pair
  for (k in seq_len(n_slices - 1) + 1) {
    j <- cut_pts[k]
    while (j < n && !is.na(idx$pair_id[j + 1]) && !is.na(idx$pair_id[j]) &&
           idx$pair_id[j + 1] == idx$pair_id[j]) j <- j + 1
    cut_pts[k] <- max(j, cut_pts[k - 1])
  }
  slices <- vector("list", n_slices)
  for (k in seq_len(n_slices)) {
    rows <- if (cut_pts[k] < cut_pts[k + 1]) (cut_pts[k] + 1):cut_pts[k + 1] else integer(0)
    pid  <- idx$pair_id[rows]
    pids <- unique(pid[!is.na(pid)])
    slices[[k]] <- list(
      sp      = idx$sp_idx[rows],
      grid    = idx$grid_idx[rows],
      pair    = match(pid, pids),            # local pair id, NA for invalid rows
      pair_sp = idx$pair_sp[pids],           # species of each local pair
      n_pairs = length(pids),
      n_sp    = idx$n_sp
    )
  }
  slices
}


# --- Worker-side functions ---------------------------------------------------
# These are sent to PSOCK workers. They use only base/stats functions, and
# sefd_clean_fn() strips their enclosing environment before sending so that
# nothing large (the index, patterns, ...) is serialized along with them.
sefd_clean_fn <- function(f) { environment(f) <- globalenv(); f }

sefd_set_slice <- function(s) {
  assign(".sefd_slice", s, envir = globalenv())
  length(s$sp)
}

# Evaluate a batch of patterns on this worker's slice.
# `survive` is a numeric vector aligned with unique_grids (0 = deforested).
sefd_worker_eval <- function(survive_list) {
  ch <- get(".sefd_slice", envir = globalenv())
  lapply(survive_list, function(survive) {
    keep <- runif(length(ch$grid)) < survive[ch$grid]
    keep[is.na(keep)] <- FALSE
    p <- ch$pair[keep]
    p <- p[!is.na(p)]
    pair_alive <- tabulate(p, ch$n_pairs) > 0
    list(
      n_kept   = as.numeric(sum(keep)),
      sp_count = tabulate(ch$sp[keep], ch$n_sp),
      n_after  = tabulate(ch$pair_sp[pair_alive], ch$n_sp)
    )
  })
}


# --- Deforestation pattern for one replicate (needs only grids + weights) ---
generate_deforestation_pattern <- function(
    freq_table, unique_grids, grid_cell_size_km2,
    grid_weights = NULL
) {
  n_grids <- length(unique_grids)

  freq_large <- freq_table[
    freq_table$bin_mid >= grid_cell_size_km2 & freq_table$count > 0, ,
    drop = FALSE
  ]

  # --- large patches: which grids to fully remove ---
  grid_equiv <- freq_large$count * (freq_large$bin_mid / grid_cell_size_km2)
  n_whole <- floor(grid_equiv)
  remainder_area <- (grid_equiv - n_whole) * grid_cell_size_km2

  total_grids_to_remove <- sum(n_whole)

  if (!is.null(grid_weights)) {
    w <- grid_weights[as.character(unique_grids)]
    w[is.na(w)] <- 0
    if (sum(w) > 0) w <- w / sum(w)
    else w <- rep(1 / n_grids, n_grids)
    total_grids_to_remove <- min(total_grids_to_remove, sum(w > 0))
  } else {
    w <- NULL
  }
  total_grids_to_remove <- min(total_grids_to_remove, n_grids)

  del_idx <- sample.int(n_grids, total_grids_to_remove, replace = FALSE, prob = w)

  # --- fold remainder into freq table for small patches ---
  freq_mod <- freq_table
  for (i in seq_along(remainder_area)) {
    patch_size <- remainder_area[i]
    if (is.na(patch_size) || patch_size <= 0) next
    bin_idx <- which(freq_mod$bin_lower <= patch_size & freq_mod$bin_upper > patch_size)
    if (length(bin_idx) == 1) {
      freq_mod$count[bin_idx] <- freq_mod$count[bin_idx] + 1L
    }
  }

  # --- small patches: survival probabilities ---
  freq_small <- freq_mod[
    freq_mod$bin_mid < grid_cell_size_km2 & freq_mod$count > 0, ,
    drop = FALSE
  ]

  survive <- rep(1.0, n_grids)
  survive[del_idx] <- 0
  surv_idx <- if (length(del_idx)) seq_len(n_grids)[-del_idx] else seq_len(n_grids)
  n_surviving <- length(surv_idx)

  if (n_surviving > 0 && nrow(freq_small) > 0) {
    if (!is.null(grid_weights)) {
      ws <- grid_weights[as.character(unique_grids[surv_idx])]
      ws[is.na(ws)] <- 0
      if (sum(ws) > 0) ws <- ws / sum(ws)
      else ws <- rep(1 / n_surviving, n_surviving)
    } else {
      ws <- NULL
    }

    s <- survive[surv_idx]
    for (i in seq_len(nrow(freq_small))) {
      perc <- min(freq_small$bin_mid[i] / grid_cell_size_km2, 1)
      hits <- sample.int(n_surviving, size = freq_small$count[i], replace = TRUE, prob = ws)
      s <- s * (1 - perc)^tabulate(hits, n_surviving)
    }
    survive[surv_idx] <- s
  }

  list(
    grids_to_deforest = unique_grids[del_idx],
    survive           = survive
  )
}


# --- Combine worker partial counts for one replicate into a result row ---
sefd_combine <- function(idx, parts) {
  n_kept   <- as.integer(sum(vapply(parts, `[[`, numeric(1), "n_kept")))
  sp_count <- Reduce(`+`, lapply(parts, `[[`, "sp_count"))
  n_after  <- Reduce(`+`, lapply(parts, `[[`, "n_after"))

  in_base  <- idx$n_before > 0
  pct_lost <- (idx$n_before[in_base] - n_after[in_base]) / idx$n_before[in_base]
  tc <- vapply(RANGE_REDUCTION_THRESHOLDS, function(t) sum(pct_lost >= t), integer(1))
  tc <- setNames(tc, RANGE_REDUCTION_NAMES)

  data.frame(
    run_id           = NA_integer_,
    n_pop_lost_total = idx$n_rows - n_kept,
    n_spp_lost_total = idx$n_sp - sum(sp_count > 0),
    n_reduced_20  = unname(tc["n_reduced_20"]),
    n_reduced_50  = unname(tc["n_reduced_50"]),
    n_reduced_75  = unname(tc["n_reduced_75"]),
    n_reduced_90  = unname(tc["n_reduced_90"]),
    n_reduced_95  = unname(tc["n_reduced_95"]),
    n_reduced_99  = unname(tc["n_reduced_99"]),
    n_reduced_100 = unname(tc["n_reduced_100"])
  )
}


# --- Apply a list of patterns using a PSOCK cluster (or sequentially) ---
# Exposed separately from SEFD_rep_parallel so it can be tested with fixed
# patterns. `cl` NULL => run in this process.
sefd_apply_patterns <- function(idx, patterns, cl = NULL, batch_size = 10,
                                progress_dir = NULL) {
  n_runs <- length(patterns)
  survive_list <- lapply(patterns, `[[`, "survive")
  results <- vector("list", n_runs)

  if (is.null(cl)) {
    assign(".sefd_slice", split_sefd_index(idx, 1)[[1]], envir = globalenv())
    on.exit(rm(".sefd_slice", envir = globalenv()), add = TRUE)
  }

  batches <- split(seq_len(n_runs), ceiling(seq_len(n_runs) / batch_size))
  for (b in batches) {
    if (is.null(cl)) {
      per_rep <- lapply(sefd_worker_eval(survive_list[b]), list)
    } else {
      per_node <- parallel::clusterCall(cl, sefd_clean_fn(sefd_worker_eval), survive_list[b])
      # per_node[[node]][[rep]] -> per_rep[[rep]][[node]]
      per_rep <- lapply(seq_along(b), function(j) lapply(per_node, `[[`, j))
    }
    for (j in seq_along(b)) {
      i <- b[j]
      results[[i]] <- sefd_combine(idx, per_rep[[j]])
      results[[i]]$run_id <- i
      if (!is.null(progress_dir)) {
        writeLines(as.character(Sys.time()),
                   file.path(progress_dir, paste0("run_", i, ".done")))
      }
    }
    if (max(b) %% 20 == 0 || max(b) == n_runs) message("  ", max(b), "/", n_runs, " replicates done")
  }
  do.call(rbind, results)
}


# --- Replicate runner (same signature/return as before) ---
# n_cores: defaults to `coresNow` in the global env if set, otherwise 90% of
# logical cores. Each worker holds only n_rows / n_cores rows.
SEFD_rep_parallel <- function(
    simAnimalObj,
    freq_table,
    grid_cell_size_km2,
    n_runs = 1,
    grid_weights = NULL,
    toroidal_params = NULL,
    progress_dir = NULL,
    n_cores = NULL
) {
  t0 <- Sys.time()
  sim_tbl <- if (inherits(simAnimalObj, "sf")) sf::st_drop_geometry(simAnimalObj) else simAnimalObj
  idx <- prepare_sefd_index(sim_tbl)
  rm(sim_tbl)
  message("Indexed ", format(idx$n_rows, big.mark = ","), " rows, ",
          idx$n_sp, " species, ", idx$n_grids, " grids in ",
          round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1), " s")

  freq_table <- as.data.frame(freq_table)

  if (is.null(progress_dir)) progress_dir <- file.path(tempdir(), "sefd_progress")
  if (dir.exists(progress_dir)) unlink(progress_dir, recursive = TRUE)
  dir.create(progress_dir, recursive = TRUE, showWarnings = FALSE)

  # ---- patterns (cheap; generated in the parent) ----
  patterns <- lapply(seq_len(n_runs), function(i) {
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
    generate_deforestation_pattern(
      freq_table         = freq_table,
      unique_grids       = idx$unique_grids,
      grid_cell_size_km2 = grid_cell_size_km2,
      grid_weights       = weights_i
    )
  })

  # ---- workers ----
  if (is.null(n_cores)) {
    n_cores <- if (exists("coresNow", envir = globalenv())) {
      get("coresNow", envir = globalenv())
    } else {
      max(1, floor(parallel::detectCores() * 0.9))
    }
  }
  # no point in more workers than ~50k-row slices
  n_workers <- max(1, min(n_cores, ceiling(idx$n_rows / 5e4)))

  if (n_workers == 1) {
    message("Running ", n_runs, " replicates sequentially")
    res <- sefd_apply_patterns(idx, patterns, cl = NULL, progress_dir = progress_dir)
  } else {
    slices <- split_sefd_index(idx, n_workers)
    message("Running ", n_runs, " replicates on ", n_workers, " PSOCK workers (",
            format(max(lengths(lapply(slices, `[[`, "sp"))), big.mark = ","), " rows max per worker)")
    cl <- parallel::makeCluster(n_workers)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterSetRNGStream(cl, iseed = sample.int(.Machine$integer.max, 1))
    # one slice per worker, sent once
    sent <- unlist(parallel::clusterApply(cl, slices, sefd_clean_fn(sefd_set_slice)))
    rm(slices)
    stopifnot(sum(sent) == idx$n_rows)
    res <- sefd_apply_patterns(idx, patterns, cl = cl, progress_dir = progress_dir)
  }

  message("All ", n_runs, " replicates complete in ",
          round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 2), " min.")
  res
}
