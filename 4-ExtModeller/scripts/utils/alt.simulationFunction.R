# --- NEW: Worker function (tiny memory footprint) ---
generate_deforestation_pattern <- function(
    freq_table, unique_grids, grid_cell_size_km2,
    grid_weights = NULL
) {
  freq_large <- freq_table[
    freq_table$bin_mid >= grid_cell_size_km2 & freq_table$count > 0, ,
    drop = FALSE
  ]

  # --- large patches: which grids to fully remove ---
  grid_equiv <- freq_large$count * (freq_large$bin_mid / grid_cell_size_km2)
  n_whole <- floor(grid_equiv)
  remainder_area <- (grid_equiv - n_whole) * grid_cell_size_km2

  total_grids_to_remove <- sum(n_whole)

  available_grids <- unique_grids
  if (!is.null(grid_weights)) {
    w <- grid_weights[as.character(available_grids)]
    w[is.na(w)] <- 0
    if (sum(w) > 0) w <- w / sum(w)
    else w <- rep(1 / length(available_grids), length(available_grids))
    n_eligible <- sum(w > 0)
    total_grids_to_remove <- min(total_grids_to_remove, n_eligible)
  } else {
    w <- NULL
  }
  total_grids_to_remove <- min(total_grids_to_remove, length(available_grids))

  grids_to_deforest <- sample(
    available_grids, total_grids_to_remove, replace = FALSE, prob = w
  )

  # --- fold remainder into freq table for small patches ---
  freq_mod <- freq_table
  remainderPatches <- data.frame(
    extraPatchSize = remainder_area,
    count = ifelse(remainder_area > 0, 1L, 0L)
  )
  for (i in seq_len(nrow(remainderPatches))) {
    patch_size <- remainderPatches$extraPatchSize[i]
    patch_count <- remainderPatches$count[i]
    if (is.na(patch_size) || patch_count < 1) next
    bin_idx <- which(freq_mod$bin_lower <= patch_size & freq_mod$bin_upper > patch_size)
    if (length(bin_idx) == 1) {
      freq_mod$count[bin_idx] <- freq_mod$count[bin_idx] + patch_count
    }
  }

  # --- small patches: survival probabilities ---
  freq_small <- freq_mod[
    freq_mod$bin_mid < grid_cell_size_km2 & freq_mod$count > 0, ,
    drop = FALSE
  ]

  surviving_grids <- setdiff(available_grids, grids_to_deforest)
  n_surviving <- length(surviving_grids)
  grid_survive_prob <- setNames(rep(1.0, n_surviving), surviving_grids)

  if (n_surviving > 0 && nrow(freq_small) > 0) {
    if (!is.null(grid_weights)) {
      ws <- grid_weights[as.character(surviving_grids)]
      ws[is.na(ws)] <- 0
      if (sum(ws) > 0) ws <- ws / sum(ws)
      else ws <- rep(1 / n_surviving, n_surviving)
    } else {
      ws <- NULL
    }

    for (i in seq_len(nrow(freq_small))) {
      perc <- min(freq_small$bin_mid[i] / grid_cell_size_km2, 1)
      n_hits <- freq_small$count[i]
      hit_grids <- sample(surviving_grids, size = n_hits, replace = TRUE, prob = ws)
      hit_counts <- table(hit_grids)
      for (g in names(hit_counts)) {
        grid_survive_prob[g] <- grid_survive_prob[g] * (1 - perc)^hit_counts[g]
      }
    }
  }

  list(
    grids_to_deforest  = grids_to_deforest,
    grid_survive_prob  = grid_survive_prob
  )
}


# --- NEW: Parent-side application (sequential, no copies sent to workers) ---
apply_deforestation_pattern <- function(sim_tbl, pattern, baseline_ranges) {
  # Remove fully deforested grids
  keep_rows <- !(sim_tbl$grid_id %in% pattern$grids_to_deforest)
  current_tbl <- sim_tbl[keep_rows, , drop = FALSE]

  # Apply small-patch thinning
  row_grid_ids <- as.character(current_tbl$grid_id)
  row_survive <- pattern$grid_survive_prob[row_grid_ids]
  row_survive[is.na(row_survive)] <- 0  # grids not in survive_prob were deforested
  keep <- runif(nrow(current_tbl)) < row_survive
  final_tbl <- current_tbl[keep, , drop = FALSE]

  # Compute results
  species_rr <- compute_range_reduction(final_tbl, baseline_ranges)
  tc <- count_range_thresholds(species_rr)

  data.frame(
    run_id           = NA_integer_,
    n_pop_lost_total = nrow(sim_tbl) - nrow(final_tbl),
    n_spp_lost_total = length(unique(sim_tbl$scientificName)) -
                       length(unique(final_tbl$scientificName)),
    n_reduced_20  = unname(tc["n_reduced_20"]),
    n_reduced_50  = unname(tc["n_reduced_50"]),
    n_reduced_75  = unname(tc["n_reduced_75"]),
    n_reduced_90  = unname(tc["n_reduced_90"]),
    n_reduced_95  = unname(tc["n_reduced_95"]),
    n_reduced_99  = unname(tc["n_reduced_99"]),
    n_reduced_100 = unname(tc["n_reduced_100"])
  )
}


# --- NEW: Replacement for SEFD_rep_parallel ---
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
  baseline_ranges <- precompute_baseline_ranges(sim_tbl)
  unique_grids <- unique(sim_tbl$grid_id)
  unique_grids <- unique_grids[!is.na(unique_grids)]

  if (is.null(progress_dir)) {
    progress_dir <- file.path(tempdir(), "sefd_progress")
  }
  if (dir.exists(progress_dir)) unlink(progress_dir, recursive = TRUE)
  dir.create(progress_dir, recursive = TRUE, showWarnings = FALSE)

  n_cores_total <- parallel::detectCores(logical = FALSE)
  n_cores_use <- max(1, floor(n_cores_total * 0.7))
  message("Using ", n_cores_use, " of ", n_cores_total, " cores for ", n_runs, " replicates")

  # ---- Phase 1: Generate deforestation patterns in parallel (tiny memory) ----
  message("Phase 1: Generating deforestation patterns in parallel...")

  cl <- parallel::makeCluster(n_cores_use)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  doParallel::registerDoParallel(cl)

  export_fns <- "generate_deforestation_pattern"
  if (!is.null(toroidal_params)) {
    export_fns <- c(export_fns, "toroidal_shift_weights")
  }

  patterns <- foreach(
    i = seq_len(n_runs),
    .export = export_fns
  ) %dopar% {
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
      unique_grids       = unique_grids,
      grid_cell_size_km2 = grid_cell_size_km2,
      grid_weights       = weights_i
    )
  }

  parallel::stopCluster(cl)
  on.exit(NULL)  # clear the on.exit since we stopped manually

  # ---- Phase 2: Apply patterns sequentially (big data stays in parent) ----
  message("Phase 2: Applying ", n_runs, " patterns to species data...")

  results_list <- vector("list", n_runs)
  for (i in seq_len(n_runs)) {
    results_list[[i]] <- apply_deforestation_pattern(
      sim_tbl, patterns[[i]], baseline_ranges
    )
    results_list[[i]]$run_id <- i  # overwrite placeholder with actual run number

    writeLines(
      as.character(Sys.time()),
      file.path(progress_dir, paste0("run_", i, ".done"))
    )

    if (i %% 10 == 0) message("  Applied ", i, "/", n_runs, " patterns")
  }

  extinction_results <- do.call(rbind, results_list)

  message("All ", n_runs, " replicates complete.")
  return(extinction_results)
}
