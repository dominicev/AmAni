####ok, finally....now take the scaling factor above (animal_max / max(fitted)) and use it to scale all the plant data.
# # should be the same as above but with a progress bar
# 
###original version of the function that works
# simulateScaledSpeciesRichness <- function(fit_pred, scaling_factor) {
#   
#   #plant_fit_max <- max(fit_pred$fitted)
#   #scaling_factor <- animal_max / plant_fit_max
#   
#   # Apply scaling
#   scaled_plant_distribution <- plant_distribution %>%
#     mutate(n_scaled = n * scaling_factor)
#   
#   # Available grid cells
#   available_grid_ids <- grid_sf$grid_id
#   available_grid_ids <- available_grid_ids[!is.na(available_grid_ids)]
#   
#   # Initialize list to store simulated species assignments
#   sim_species_records <- list()
#   species_counter <- 1
#   
#   ## 🔹 Progress bar setup (species-level)
#   total_species <- sum(round(scaled_plant_distribution$n_scaled), na.rm = TRUE)
#   
#   pb <- progress_bar$new(
#     format = "  Simulating species [:bar] :current/:total (:percent) eta: :eta",
#     total  = total_species,
#     clear  = FALSE,
#     width  = 60
#   )
#   
#   ## 🔹 Simulation loops
#   for (i in seq_len(nrow(scaled_plant_distribution))) {
#     
#     num_sites   <- scaled_plant_distribution$num_grid_cells[i]
#     num_species <- round(scaled_plant_distribution$n_scaled[i])
#     
#     for (j in seq_len(num_species)) {
#       
#       pb$tick()  # ← progress update
#       
#       if (num_sites <= length(available_grid_ids)) {
#         
#         sim_name <- str_glue("SimSpp{str_pad(species_counter, 7, pad = '0')}")
#         
#         assigned_cells <- sample(available_grid_ids, num_sites, replace = FALSE)
#         
#         sim_occurrences <- data.frame(
#           scientificName = rep(sim_name, num_sites),
#           grid_id = assigned_cells,
#           stringsAsFactors = FALSE
#         )
#         
#         sim_species_records[[species_counter]] <- sim_occurrences
#         species_counter <- species_counter + 1
#       }
#     }
#   }
#   
#   # Combine into one large data frame
#   sim_species_df <<- bind_rows(sim_species_records)
#   
#   # Join to grid_sf to get geometry
#   sim_species_sf <<- sim_species_df %>%
#     left_join(grid_sf, by = "grid_id") %>%
#     st_as_sf()
#   
#   ####the line below gives how many species of cockroaches there are under this simulation
#   output<-length(unique(sim_species_sf$scientificName))
#   return(output)
#   
#  #return(sim_species_records)
# }

####new version split into two functions
##inner loop function
# simulate_species_records <- function(
    #     scaled_plant_distribution,
#     available_grid_ids,
#     show_progress = TRUE
# ) {
#   rm(sim_records)
#   sim_species_records <- vector("list", length = sum(round(
#     scaled_plant_distribution$n_scaled
#   ), na.rm = TRUE))
#   
#   species_counter <- 1
#   
#   ## Progress bar (optional)
#   if (show_progress) {
#     total_species <- length(sim_species_records)
#     
#     pb <- progress::progress_bar$new(
#       format = "  Simulating species [:bar] :current/:total (:percent) eta: :eta",
#       total  = total_species,
#       clear  = FALSE,
#       width  = 60
#     )
#   }
#   
#   for (i in seq_len(nrow(scaled_plant_distribution))) {
#     
#     num_sites   <- scaled_plant_distribution$num_grid_cells[i]
#     num_species <- round(scaled_plant_distribution$n_scaled[i])
#     
#     for (j in seq_len(num_species)) {
#       
#       if (show_progress) pb$tick()
#       
#       if (num_sites <= length(available_grid_ids)) {
#         
#         sim_name <- stringr::str_glue(
#           "SimSpp{stringr::str_pad(species_counter, 7, pad = '0')}"
#         )
#         
#         assigned_cells <- sample(
#           available_grid_ids,
#           num_sites,
#           replace = FALSE
#         )
#         
#         sim_species_records[[species_counter]] <- data.frame(
#           scientificName = rep(sim_name, num_sites),
#           grid_id = assigned_cells,
#           stringsAsFactors = FALSE
#         )
#         
#         species_counter <- species_counter + 1
#       }
#     }
#   }
#   
#   sim_records<<-sim_species_records[seq_len(species_counter - 1)] #this might be wrong. if the line below is right then this line is unnecessary but i am keeping it in case it is actually correct
#   return(sim_species_records[seq_len(species_counter - 1)])
# }
#######Same as the above function but with simulated species done with local clustering
#
# Precomputes a sparse adjacency matrix and derives an adjacency list for fast
# contiguous range growth. At each step, the frontier (all 8-connected
# neighbours of the current range not already in the range) is maintained
# incrementally and one frontier cell is sampled proportional to its
# plant-occurrence weight.



build_adjacency <- function(grid_sf, grid_size_m) {
  # Build an 8-connected adjacency list from cell centroid coordinates.
  # Returns a list of length nrow(grid_sf) where each element is an integer
  # vector of row indices of neighbouring cells.
  
  require(Matrix)
  
  cx <- grid_sf$cx
  cy <- grid_sf$cy
  n  <- length(cx)
  thresh <- grid_size_m * 1.5  # captures orthogonal + diagonal neighbours
  
  message("  Building adjacency structure for ", n, " grid cells...")
  
  # Build sparse adjacency matrix in chunks to limit memory from outer()
  chunk_size <- 2000
  from_list <- vector("list", ceiling(n / chunk_size))
  to_list   <- vector("list", ceiling(n / chunk_size))
  
  for (ch in seq_len(ceiling(n / chunk_size))) {
    idx_start <- (ch - 1) * chunk_size + 1
    idx_end   <- min(ch * chunk_size, n)
    chunk_idx <- idx_start:idx_end
    
    dx <- outer(cx[chunk_idx], cx, function(a, b) abs(a - b))
    dy <- outer(cy[chunk_idx], cy, function(a, b) abs(a - b))
    
    is_nb <- dx <= thresh & dy <= thresh
    for (r in seq_along(chunk_idx)) {
      is_nb[r, chunk_idx[r]] <- FALSE
    }
    
    wh <- which(is_nb, arr.ind = TRUE)
    from_list[[ch]] <- chunk_idx[wh[, 1]]
    to_list[[ch]]   <- wh[, 2]
  }
  
  from_vec <- unlist(from_list)
  to_vec   <- unlist(to_list)
  
  adj_mat <- sparseMatrix(i = from_vec, j = to_vec, x = 1, dims = c(n, n))
  
  # Extract adjacency list from column-compressed sparse matrix (fast).
  # @p is 0-based: column i occupies 0-based slots p[i]..p[i+1]-1, i.e. 1-based
  # slots (p[i]+1)..p[i+1]. Empty columns (isolated cells) must yield integer(0).
  adj_list <- vector("list", n)
  for (i in seq_len(n)) {
    p0 <- adj_mat@p[i]
    p1 <- adj_mat@p[i + 1]
    adj_list[[i]] <- if (p1 > p0) adj_mat@i[(p0 + 1):p1] + 1L else integer(0)
  }
  
  message("  Adjacency: ", n, " cells, ~", round(mean(lengths(adj_list)), 1),
          " neighbours/cell.")
  adj_list
}


connected_components <- function(adj_list) {
  # Label the 8-connected components of the adjacency graph. Returns an
  # integer vector of length length(adj_list) giving each cell's component id.
  # A cell's reachable range-growth area equals the size of its component,
  # so this is what determines whether a target range size is achievable.
  n     <- length(adj_list)
  comp  <- integer(n)
  stack <- integer(n)  # preallocated DFS stack; each cell is pushed at most once
  cid   <- 0L
  for (s in seq_len(n)) {
    if (comp[s] != 0L) next
    cid       <- cid + 1L
    sp        <- 1L
    stack[sp] <- s
    comp[s]   <- cid
    while (sp > 0L) {
      v  <- stack[sp]; sp <- sp - 1L
      nb <- adj_list[[v]]
      nb <- nb[comp[nb] == 0L]
      if (length(nb)) {
        comp[nb]                     <- cid
        stack[(sp + 1L):(sp + length(nb))] <- nb
        sp                           <- sp + length(nb)
      }
    }
  }
  comp
}


# ─────────────────────────────────────────────────────────────────────────────
# simulate_species_chunk(): simulate a subset ("chunk") of species jobs.
#
# This is the sequential per-species algorithm previously inlined directly in
# simulate_species_records_clustered()'s double loop, factored out so it can
# be run either once (serial path, one chunk = every job, identical behaviour
# to the pre-parallel implementation) or many times across parallel workers
# (one call per chunk). All arguments are read-only and safe to share
# across workers; no shared mutable state is required because every chunk
# gets its own frontier-growth buffers.
#
# num_sites_vec: integer vector, one target range size (pre-jitter) per
#                species job assigned to this chunk.
# tick_fun:      optional zero-arg callback invoked once per job (used for
#                per-species progress reporting on the serial path only —
#                see simulate_species_records_clustered()).
# ─────────────────────────────────────────────────────────────────────────────
simulate_species_chunk <- function(
    num_sites_vec,
    available_grid_ids,
    all_ids,
    n_cells,
    adj_list,
    w_vec,
    sampling_weights,
    maxCellRange,
    jitter_range,
    size_aware_start,
    avail_reach,        # NULL when size_aware_start = FALSE
    simulate_ranges,
    compactness_gamma = 0,  # frontier weight is multiplied by (# of range neighbours)^gamma;
                            # 0 = off (identical to the pre-compactness algorithm), larger
                            # values fill concavities and give rounder, less porous ranges
    chunk_id = 1L,      # used to make scientificName unique across chunks
    chunk_dir = NULL,   # if non-NULL, write this chunk's records here and return the path
    tick_fun = NULL
) {
  n_jobs  <- length(num_sites_vec)
  # Compact accumulation: one integer vector of assigned grid_ids per species.
  # Far lighter than a list of one data.frame per species, which blew up memory
  # at >100k species.
  cells_list <- vector("list", n_jobs)
  n_filled <- 0L
  
  # ── Frontier exhaustion diagnostics (accumulated within this chunk) ────────
  n_exhausted     <- 0L
  total_shortfall <- 0L
  max_shortfall   <- 0L
  
  # Pre-allocations, reused across every species processed by this chunk —
  # avoids per-species alloc/reset, same rationale as the original single loop.
  visit_gen         <- integer(n_cells)
  frontier_gen      <- integer(n_cells)
  current_gen       <- 0L
  frontier_rows_buf <- integer(n_cells)
  frontier_pos      <- integer(n_cells)
  use_compact       <- compactness_gamma > 0
  # Number of already-in-range neighbours per cell; reset after each species
  nb_count          <- if (use_compact) integer(n_cells) else NULL
  
  for (j in seq_len(n_jobs)) {
    
    if (!is.null(tick_fun)) tick_fun()
    
    num_sites <- num_sites_vec[j]
    
    # ── Per-species range size jitter ────────────────────────────────────
    # Draw target from Poisson(lambda = bin representative) so species
    # within the same bin receive different range sizes, dissolving the
    # point-mass peaks created by assigning all species in a bin the
    # identical median value.
    ##smooth out range sizes larger than 50 cells (which is the same threshold for log binning below)
    # #I added this condition because the binning
    ##behavior for larger ranged species was resulting in
    #  #multi-modal distribution. but doing the poisson jitter
    # #for smaller ranges smooths out the range sizes in an unnatrual way (that doesn't fit the log-normal distribution)
    
    num_sites_j <- if (jitter_range && num_sites > 50) {
      lambda <- min(num_sites, maxCellRange)          # cap lambda before rpois to prevent integer overflow
      max(1L, min(as.integer(round(rpois(1L, lambda))), maxCellRange + 1L))
    } else {
      as.integer(min(num_sites, maxCellRange + 1L))   # also cap the non-jitter path
    }
    
    if (num_sites_j > maxCellRange) next
    
    # Draw starting cell using plant-occurrence weights. When size-aware, keep
    # only candidate starts whose reachable component is large enough to hold the
    # target range; fall back to the largest reachable region if none qualifies.
    if (size_aware_start) {
      eligible <- !is.na(avail_reach) & avail_reach >= num_sites_j
      if (!any(eligible)) eligible <- !is.na(avail_reach) & avail_reach == max(avail_reach, na.rm = TRUE)
      cand_ids <- available_grid_ids[eligible]
      cand_w   <- sampling_weights[eligible]
      starting_cell <- if (length(cand_ids) == 1L) {
        cand_ids
      } else {
        sample(cand_ids, 1, prob = cand_w, replace = FALSE)
      }
    } else {
      starting_cell <- sample(available_grid_ids, 1,
                              prob = sampling_weights, replace = FALSE)
    }
    start_row <- match(starting_cell, all_ids)
    
    if (!simulate_ranges) {
      # Range growth skipped: only the species count (gamma) is needed, so
      # record each species at its single starting cell. This avoids the
      # expensive contiguous frontier growth without changing the species
      # count, which is determined by the num_sites_j <= maxCellRange test
      # above rather than by the realised range.
      assigned_cells <- starting_cell
    } else if (num_sites_j == 1) {
      assigned_cells <- starting_cell
    } else {
      # --- Contiguous weighted growth via adjacency list ---
      # Maintain the frontier incrementally: when a cell is added to the
      # range, remove it from the frontier and add its new neighbors.
      # Each frontier cell is sampled proportional to its plant weight.
      
      range_rows <- integer(num_sites_j)
      range_rows[1] <- start_row
      
      current_gen          <- current_gen + 1L
      visit_gen[start_row] <- current_gen
      actual_size          <- 1L
      flen                 <- 0L
      
      init_nbs <- adj_list[[start_row]]
      init_nbs <- init_nbs[visit_gen[init_nbs] != current_gen]
      flen     <- length(init_nbs)
      if (use_compact) nb_count[init_nbs] <- 1L
      if (flen > 0L) {
        frontier_rows_buf[seq_len(flen)] <- init_nbs
        frontier_gen[init_nbs]           <- current_gen
        frontier_pos[init_nbs]           <- seq_len(flen)
      }
      
      for (step in 2:num_sites_j) {
        if (flen == 0L) break
        
        active <- frontier_rows_buf[seq_len(flen)]
        fw     <- w_vec[active]
        if (use_compact) fw <- fw * nb_count[active]^compactness_gamma
        if (sum(fw) == 0) fw[] <- 1
        
        new_row <- if (flen == 1L) active else sample(active, 1L, prob = fw)
        
        actual_size             <- actual_size + 1L
        range_rows[actual_size] <- new_row
        visit_gen[new_row]      <- current_gen
        
        # O(1) swap-and-pop removal from frontier buffer
        pos <- frontier_pos[new_row]
        if (pos < flen) {
          last_row               <- frontier_rows_buf[flen]
          frontier_rows_buf[pos] <- last_row
          frontier_pos[last_row] <- pos
        }
        frontier_pos[new_row] <- 0L
        frontier_gen[new_row] <- 0L
        flen <- flen - 1L
        
        # Add unvisited non-frontier neighbours (vectorised)
        new_nbs <- adj_list[[new_row]]
        if (use_compact) nb_count[new_nbs] <- nb_count[new_nbs] + 1L
        add_nbs <- new_nbs[
          visit_gen[new_nbs]    != current_gen &
            frontier_gen[new_nbs] != current_gen
        ]
        n_add <- length(add_nbs)
        if (n_add > 0L) {
          new_pos                    <- (flen + 1L):(flen + n_add)
          frontier_rows_buf[new_pos] <- add_nbs
          frontier_pos[add_nbs]      <- new_pos
          frontier_gen[add_nbs]      <- current_gen
          flen                       <- flen + n_add
        }
      }
      if (use_compact) {
        # Reset only the entries touched by this species (neighbours of its range)
        nb_count[unlist(adj_list[range_rows[seq_len(actual_size)]], use.names = FALSE)] <- 0L
      }
      #########################
      # ── Log frontier exhaustion ───────────────────────────────────────────
      if (actual_size < num_sites_j) {
        shortfall       <- num_sites_j - actual_size
        n_exhausted     <- n_exhausted + 1L
        total_shortfall <- total_shortfall + shortfall
        if (shortfall > max_shortfall) max_shortfall <- shortfall
      }
      
      assigned_cells <- all_ids[range_rows[1:actual_size]]
    }
    
    n_filled <- n_filled + 1L
    cells_list[[n_filled]] <- assigned_cells
  }
  
  # ── Assemble this chunk's records as a single compact data.table ──────────
  # Names are made unique per chunk (chunk_id prefix + local species index) so
  # workers need no shared counter and no post-hoc global renaming pass.
  if (n_filled > 0L) {
    cells_list     <- cells_list[seq_len(n_filled)]
    n_cells_per_sp <- lengths(cells_list)
    sim_names      <- sprintf("SimSpp_c%05d_%07d", chunk_id, seq_len(n_filled))
    chunk_dt <- data.table::data.table(
      scientificName = rep(sim_names, n_cells_per_sp),
      grid_id        = unlist(cells_list, use.names = FALSE)
    )
  } else {
    chunk_dt <- data.table::data.table(
      scientificName = character(0),
      grid_id        = integer(0)
    )
  }
  
  # Write to disk (default) and return only the path + summaries, so parallel
  # workers ship back tiny payloads and the caller never holds every chunk's
  # records in memory at once.
  chunk_file <- NULL
  if (!is.null(chunk_dir)) {
    chunk_file <- file.path(chunk_dir, sprintf("chunk_%05d.rds", chunk_id))
    saveRDS(chunk_dt, chunk_file, compress = FALSE)
    chunk_dt <- NULL
  }
  
  list(
    file            = chunk_file,
    records         = chunk_dt,     # NULL when written to disk
    n_species       = n_filled,
    n_exhausted     = n_exhausted,
    total_shortfall = total_shortfall,
    max_shortfall   = max_shortfall
  )
}


simulate_species_records_clustered <- function(
    scaled_plant_distribution,
    available_grid_ids,
    grid_sf,
    adj_list,
    cell_weights,
    sampling_weights = sampling_weights_neutral,
    simulate_ranges = TRUE,  # if FALSE, skip the contiguous range growth and record
    # each species only at its starting cell. Gamma richness
    # and the CI bounds depend solely on species counts (set
    # before range growth), so this is a fast path for them.
    show_progress = TRUE,
    jitter_range  = TRUE,    # if TRUE, Poisson-jitter each species target range size
    size_aware_start = TRUE, # if TRUE, restrict a species' starting cell to connected
                             # components at least as large as its target range size,
                             # so large ranges are not seeded into small isolated regions
    compactness_gamma = get0(".amani_compactnessGamma", envir = .GlobalEnv, ifnotfound = 0),
                             # compactness bias for range growth (see simulate_species_chunk());
                             # 0 restores the original unbiased frontier sampling
    estimateBounds=TRUE,                      # so species within the same bin get different sizes
    # (eliminates point-mass peaks in the output distribution)
    parallel = FALSE,       # if TRUE, distribute species jobs across n_cores workers
                             # via furrr/future on a `multisession` plan (see
                             # fill_all_species_ranges() in fillingRanges.R for the
                             # same pattern). Falls back to serial processing if the
                             # cluster can't be set up. Because each species' range
                             # growth is independent given its starting cell, jobs can
                             # be split across workers without changing the result
                             # distribution — but per-worker RNG streams mean a
                             # parallel run will not reproduce a serial run bit-for-bit
                             # even with the same seed.
    n_cores  = max(1L, parallelly::availableCores() - 1L),
    chunk_target_species = 5000L,  # serial path: species per on-disk chunk
    chunk_dir = NULL               # directory for per-chunk RDS files
) {
  # adj_list:     list of integer vectors — adj_list[[i]] = row indices of
  #               neighbours of cell i in grid_sf.
  # cell_weights: named numeric vector keyed by grid_id — plant-occurrence
  #               weights used to bias frontier cell selection during range growth.
  
  n_cells <- nrow(grid_sf)
  all_ids <- grid_sf$grid_id
  # Evaluate now so the value is fixed before run_chunk() is shipped to workers
  compactness_gamma <- as.numeric(compactness_gamma)
  stopifnot(length(compactness_gamma) == 1L, !is.na(compactness_gamma), compactness_gamma >= 0)
  
  # Weight vector aligned to grid_sf row order
  w_vec <- cell_weights[as.character(all_ids)]
  w_vec[is.na(w_vec)] <- 0
  
  #define this so we don't have to recalculate it inside the loop
  maxCellRange <- length(available_grid_ids)
  avail_reach  <- NULL
  
  # ── Size-aware starting-cell precompute ──────────────────────────────────────
  # A species can only grow within the connected component of its start cell, so
  # the reachable area for a candidate start equals that component's size. Precompute
  # the reachable size for every available start cell once; per species we then keep
  # only starts whose reachable area >= the target range size.
  # size_aware_start requires an adjacency list; when simulate_ranges = FALSE

  # adj_list is NULL (no range growth), so size-aware filtering is meaningless.
  if (size_aware_start && !is.null(adj_list)) {
    max_id            <- max(all_ids)
    id_to_row         <- integer(max_id)
    id_to_row[all_ids] <- seq_len(n_cells) # O(1) grid_id -> row lookup
    comp_id_row    <- connected_components(adj_list)
    comp_sizes     <- tabulate(comp_id_row)
    reach_size_row <- comp_sizes[comp_id_row]              # reachable area per cell (row order)
    avail_reach    <- reach_size_row[id_to_row[available_grid_ids]]  # aligned to available_grid_ids
    message(sprintf(
      "  Size-aware start: %d connected component(s); largest reachable region = %d cells.",
      length(comp_sizes), max(comp_sizes)
    ))
  } else if (size_aware_start && is.null(adj_list)) {
    message("  Size-aware start disabled (no adjacency list; simulate_ranges is FALSE).")
    size_aware_start <- FALSE
  }
  
  # ── Flatten the (bin, species) double loop into one job vector ─────────────
  # Each element is the target range size (pre-jitter) for one species. Jobs
  # stay in the original bin-major order, so with a single chunk (the default,
  # parallel = FALSE) the RNG draw sequence and results are identical to the
  # pre-parallel implementation.
  bin_num_sites   <- scaled_plant_distribution$num_grid_cells
  bin_num_species <- round(scaled_plant_distribution$n_scaled)
  keep <- !is.na(bin_num_sites) & !is.na(bin_num_species) & bin_num_species >= 1
  num_sites_vec <- rep(bin_num_sites[keep], bin_num_species[keep])
  n_jobs <- length(num_sites_vec)
  
  use_parallel <- isTRUE(parallel) && n_cores > 1L && n_jobs > 1L
  
  # ── Chunking strategy ─────────────────────────────────────────────────────
  # Parallel: more chunks than workers, jobs assigned round-robin. num_sites_vec
  # is sorted by ascending range size (bin order), so contiguous blocks would
  # hand the last worker every large, slow-to-grow species; round-robin balances
  # load. Serial: contiguous blocks sized by chunk_target_species so records are
  # streamed to disk one chunk at a time (bounding memory) while preserving the
  # original in-order RNG draw sequence — contiguous + in-order gives identical
  # draws to the pre-chunking single loop.
  if (use_parallel) {
    n_chunks  <- min(n_jobs, max(n_cores * 4L, n_cores))
    chunk_key <- (seq_len(n_jobs) - 1L) %% n_chunks
  } else {
    n_chunks  <- max(1L, ceiling(n_jobs / chunk_target_species))
    chunk_key <- (seq_len(n_jobs) - 1L) %/% chunk_target_species
  }
  job_chunks <- unname(split(num_sites_vec, chunk_key))
  n_chunks   <- length(job_chunks)
  
  # Directory for per-chunk RDS files (created here, cleaned up after combine).
  if (is.null(chunk_dir)) {
    base_dir <- if (exists("project_root", envir = .GlobalEnv)) {
      file.path(get("project_root", envir = .GlobalEnv), "output", "shared")
    } else tempdir()
    chunk_dir <- file.path(
      base_dir, sprintf("sim_chunks_%s", format(Sys.time(), "%Y%m%d_%H%M%S"))
    )
  }
  dir.create(chunk_dir, showWarnings = FALSE, recursive = TRUE)
  on.exit(unlink(chunk_dir, recursive = TRUE), add = TRUE)
  
  run_chunk <- function(chunk_sites, chunk_id, tick_fun = NULL) {
    simulate_species_chunk(
      num_sites_vec       = chunk_sites,
      available_grid_ids  = available_grid_ids,
      all_ids             = all_ids,
      n_cells             = n_cells,
      adj_list            = adj_list,
      w_vec               = w_vec,
      sampling_weights    = sampling_weights,
      maxCellRange        = maxCellRange,
      jitter_range        = jitter_range,
      size_aware_start    = size_aware_start,
      avail_reach         = avail_reach,
      simulate_ranges     = simulate_ranges,
      compactness_gamma   = compactness_gamma,
      chunk_id            = chunk_id,
      chunk_dir           = chunk_dir,
      tick_fun            = tick_fun
    )
  }
  
  chunk_results <- if (use_parallel) {
    tryCatch({
      if (!requireNamespace("furrr", quietly = TRUE) ||
          !requireNamespace("future", quietly = TRUE)) {
        stop("packages 'furrr' and 'future' are required for parallel = TRUE")
      }
      message(sprintf(
        "  Simulating %d species across %d chunk(s) on %d multisession worker(s)...",
        n_jobs, n_chunks, n_cores
      ))
      old_plan <- future::plan(future::multisession, workers = n_cores)
      on.exit(future::plan(old_plan), add = TRUE)
      
      furrr::future_map2(
        job_chunks,
        seq_along(job_chunks),
        run_chunk,
        .options  = furrr::furrr_options(seed = TRUE),
        .progress = isTRUE(show_progress)
      )
    }, error = function(e) {
      warning("simulate_species_records_clustered(): parallel backend failed (",
              conditionMessage(e), "); falling back to serial processing.",
              call. = FALSE)
      NULL
    })
  } else {
    NULL
  }
  
  if (is.null(chunk_results)) {
    # Serial path. With a single chunk this ticks once per species, matching the
    # original per-species progress bar; with several chunks progress is per chunk.
    tick_per_job <- n_chunks == 1L
    if (show_progress) {
      pb <- progress::progress_bar$new(
        format = "  Simulating species [:bar] :current/:total (:percent) eta: :eta",
        total  = if (tick_per_job) n_jobs else n_chunks,
        clear  = FALSE,
        width  = 60
      )
      chunk_results <- lapply(seq_along(job_chunks), function(k) {
        r <- run_chunk(job_chunks[[k]], k,
                       tick_fun = if (tick_per_job) function() pb$tick() else NULL)
        if (!tick_per_job) pb$tick()
        r
      })
    } else {
      chunk_results <- lapply(seq_along(job_chunks),
                              function(k) run_chunk(job_chunks[[k]], k))
    }
  }
  
  # ── Combine chunk results ────────────────────────────────────────────────
  # Read the compact per-chunk data.tables back and stack them into one
  # occurrence table (character name + integer grid_id) — a single compact copy
  # rather than a list of one data.frame per species.
  chunk_files <- vapply(chunk_results, function(r) {
    if (is.null(r$file)) NA_character_ else r$file
  }, character(1))
  chunk_files <- chunk_files[!is.na(chunk_files)]
  
  simulated_occurrence_dt <- data.table::rbindlist(lapply(chunk_files, readRDS))
  
  n_exhausted     <- sum(vapply(chunk_results, `[[`, integer(1), "n_exhausted"))
  total_shortfall <- sum(vapply(chunk_results, `[[`, integer(1), "total_shortfall"))
  max_shortfall   <- if (length(chunk_results)) {
    max(vapply(chunk_results, `[[`, integer(1), "max_shortfall"))
  } else 0L
  n_total <- sum(vapply(chunk_results, `[[`, integer(1), "n_species"))
  
  # Expose the combined occurrence table for downstream export. Replaces the old
  # `sim_records` list-of-data.frames global (see 02_animalAnalysis.R).
  simulated_occurrence_dt <<- simulated_occurrence_dt
  
  # ── Frontier exhaustion summary ───────────────────────────────────────────────
  if (n_exhausted > 0L) {
    message(sprintf(
      "  Frontier exhaustion: %d / %d species (%.1f%%) had ranges truncated.",
      n_exhausted, n_total, 100 * n_exhausted / n_total
    ))
    message(sprintf(
      "  Mean shortfall: %.1f cells | Max shortfall: %d cells | Total cells lost: %d",
      total_shortfall / n_exhausted, max_shortfall, total_shortfall
    ))
  } else {
    message("  Frontier exhaustion: none detected.")
  }
  
  n_total
}




#wrapper function
simulateScaledSpeciesRichness <- function(fit_pred, extent_scaling_factor, sampling_weights = sampling_weights_neutral,
                                          frontier_weights = NULL,  # frontier-growth weights aligned to available_grid_ids; NULL = use sampling_weights
                                          simulate_ranges = if (exists(".amani_simulateRanges", envir = .GlobalEnv)) .amani_simulateRanges else TRUE,
                                          parallel = get0(".amani_simulateParallel", envir = .GlobalEnv, ifnotfound = FALSE),
                                          n_cores  = get0(".amani_simulateNCores",  envir = .GlobalEnv,
                                                          ifnotfound = max(1L, parallelly::availableCores() - 1L))) {
  if (isTRUE(simulate_ranges)){
    print(paste0("Simulating Amazonian species with an extent scaling factor of: ", extent_scaling_factor))}else{message(paste0("Skipping range simulation. Extent scaling factor is: ", extent_scaling_factor))}
  
  
  #   scaling_factor<-0.065 #for testing
  # # Apply scaling
  # scaled_plant_distribution <- plant_distribution %>% #fit_pred should NOT be here. This MUST be plant_distribution, which is the same as fit_pred but for the whole amazon. 
  #   dplyr::mutate(n_scaled = n * scaling_factor)
  
  ###leaving scaled_plant_distribution as it is above will be a problem because it results in us dropping more species with large ranges than expected (e.g., if the scaling factor is 0.5 then we drop all range sizes occupied by <2 species. When the scaling factor is very small (e.g. <0.1) this could result in over 1000 species dropped)
  
  # coerce plant_distribution to data.table without copying if possible
  pd <- data.table::as.data.table(plant_distribution)             # expects column num_grid_cells
  
  # ── Log-spaced binning ──────────────────────────────────────────────────
  # Replace the old hard-edged tiered binning (bin_num_grids) with
  # geometrically spaced bin edges. This avoids the abrupt bin-width
  # transitions at 50/200/600 that created artificial density peaks
  # ("attractor values") in the simulated range-size distribution.
  #
  # Log-spacing naturally gives fine resolution at small range sizes
  # (where most species cluster) and progressively coarser resolution
  # at large range sizes — matching the spirit of the original scheme
  # but without hard discontinuities.
  
  n_log_bins <- 80                                                  # tunable; ~80 gives similar compression to old scheme
  max_cells  <- max(pd$num_grid_cells, na.rm = TRUE)
  
  # keep range sizes 1–50 unbinned (exact), bin everything above
  # using geometrically spaced edges from 50 to max
  unbinned_threshold <- 50
  
  # rows that stay unbinned (small range sizes, high ecological importance)
  pd_small <- pd[num_grid_cells <= unbinned_threshold]
  pd_large <- pd[num_grid_cells >  unbinned_threshold]
  
  if (nrow(pd_large) > 0) {
    # create log-spaced bin edges from the threshold to the max value
    log_edges <- exp(seq(
      log(unbinned_threshold),
      log(max_cells + 1),                                           # +1 so max value falls inside last bin
      length.out = n_log_bins + 1
    ))
    
    # assign each large-range row to a log-spaced bin
    pd_large[, bin_id := findInterval(num_grid_cells, log_edges, rightmost.closed = TRUE)]
    
    # aggregate: median of original values as representative, sum of counts
    pd_large_binned <- pd_large[
      , .(
        num_grid_cells = round(median(num_grid_cells)),
        n = sum(n)
      ),
      by = bin_id
    ][order(num_grid_cells)]
    pd_large_binned[, bin_id := NULL]
    
    # recombine small (exact) and large (log-binned) portions
    plant_distribution_binned <- data.table::rbindlist(
      list(pd_small[, .(num_grid_cells, n)], pd_large_binned)
    )[order(num_grid_cells)]
    
  } else {
    # edge case: all range sizes <= threshold (nothing to bin)
    plant_distribution_binned <- pd_small[, .(num_grid_cells, n)]
  }
  
  # sanity check: total species count is conserved
  stopifnot(abs(sum(plant_distribution_binned$n) - sum(pd$n)) < 0.5)
  
  # ── Apply shape transformation (meanlog shift + sdlog rescaling) ────────
  # If the user selected "meanlog" or "sdlog" in fit_params, the empirical
  # plant distribution's range sizes are warped in log-space before amplitude
  # scaling. This changes the proportion of endemic vs widespread species
  # while preserving the empirical texture of the distribution.
  delta_ml  <- if (exists(".amani_delta_meanlog",     envir = .GlobalEnv)) .amani_delta_meanlog     else 0
  ratio_sd  <- if (exists(".amani_ratio_sdlog",       envir = .GlobalEnv)) .amani_ratio_sdlog       else 1
  ref_ml    <- if (exists(".amani_plant_mdd_meanlog", envir = .GlobalEnv)) .amani_plant_mdd_meanlog else log(2)
  
  if (delta_ml != 0 || ratio_sd != 1) {
    message("  Applying shape transformation: delta_meanlog = ", round(delta_ml, 4),
            ", ratio_sdlog = ", round(ratio_sd, 4))
    plant_distribution_binned <- as.data.frame(plant_distribution_binned)
    plant_distribution_binned$log_x_new <-
      ref_ml + delta_ml +
      (log(plant_distribution_binned$num_grid_cells) - ref_ml) * ratio_sd
    # plant_distribution_binned$num_grid_cells <-
    #   pmax(1L, as.integer(round(exp(plant_distribution_binned$log_x_new))))
    # plant_distribution_binned$log_x_new <- NULL
    # 
    # # Re-aggregate: multiple original bins may map to the same integer after transform
    # plant_distribution_binned <- plant_distribution_binned |>
    #   dplyr::summarise(n = sum(n), .by = num_grid_cells) |>
    #   dplyr::arrange(num_grid_cells)
    
    plant_distribution_binned$num_grid_cells <-
      pmax(1, round(exp(plant_distribution_binned$log_x_new)))  # drop as.integer()
    plant_distribution_binned$log_x_new <- NULL
    
    plant_distribution_binned <- plant_distribution_binned |>
      dplyr::summarise(n = sum(n), .by = num_grid_cells) |>
      dplyr::filter(!is.na(num_grid_cells)) |>          # drop overflow rows
      dplyr::arrange(num_grid_cells)
    
    }
  
  scaled_plant_distribution <- plant_distribution_binned %>%
    dplyr::mutate(n_scaled = ceiling(n * extent_scaling_factor))
  
  # Will hold fullGamma / fullGamma_lwr / fullGamma_upr from the sampling
  # extension model if that path is taken. Stays NULL otherwise, triggering
  # the legacy single-integer return at the end of the function.
  sampling_model_bounds <- NULL

  if (isTRUE(get0(".amani_useSamplingExtension", envir = .GlobalEnv, ifnotfound = FALSE))) {

    # Load the sampling extension model (only once — avoids repeated disk reads
    # across the best-estimate, upper-bound, and lower-bound calls).
    if (!exists("extend_richness_from_sampling_predict", mode = "function")) {
      source("./scripts/utils/SamplingExtensionModel.R", local = FALSE)
    }

    message("Calculating input for sampling-based linear model extension...")
    #compute input parameters for the above function

    occ_data <- animal_grid  # animal_grid is the occurrence data available in global env
    .ssr_comp_method <- if (exists(".amani_completenessMethod", envir = .GlobalEnv)) .amani_completenessMethod else "Chao1"
    if (length(.ssr_comp_method) != 1 || !.ssr_comp_method %in% c("Chao1", "ACE", "Coverage")) .ssr_comp_method <- "Chao1"
    cell_statistics <- compute_cell_stats(occ_data, method = .ssr_comp_method)

    comp <- cell_statistics[[2]]$completeness

    comp_spread            <- max(comp) - min(comp)
    comp_mean              <- mean(comp)
    comp_sd                <- sd(comp)
    # prop_species: overall (regional) proportion of species sampled,
    # estimated via Chao2 using grid cells as incidence sampling units
    # (rather than a per-cell aggregate like median completeness).
    chao2_result            <- compute_chao2_completeness(cell_statistics$spp_abund)
    prop_species            <- chao2_result$completeness
    n_grid_cells           <- length(unique(occ_data$grid_id))
    # weighted.mean over the frequency table gives one mean per species (equivalent
    # to mean(morphotype_grid_counts$num_grid_cells)); mean() on the data.frame would fail
    sample_mean_range_size <- weighted.mean(morphotype_distribution$num_grid_cells,
                                            morphotype_distribution$n)
    interval               <- c("prediction", "confidence")
    level                  <- 0.95

    partialGamma <- sum(scaled_plant_distribution$n_scaled)

    samplingModelResult <- extend_richness_from_sampling(
      partialGamma           = partialGamma,
      comp_spread            = comp_spread,
      comp_mean              = comp_mean,
      prop_species           = prop_species,
      comp_sd                = comp_sd,
      n_grid_cells           = n_grid_cells,
      sample_mean_range_size = sample_mean_range_size,
      delta_meanlog          = delta_ml,
      ratio_sdlog            = ratio_sd,
      fill_Plants            = isTRUE(get0(".amani_fill_plants",
                                          envir = .GlobalEnv, ifnotfound = FALSE)),
      interval               = interval,
      level                  = level
    )

    # Guard: a non-sensical signed_prop_error (1 + error <= 0) makes the model
    # return fullGamma = NA. Left unhandled this turns every n_scaled into NA
    # (zero simulated species) and NA error bounds. Fall back to
    # fullGamma = partialGamma (correction_factor = 1). The lower bound is
    # already floored at partialGamma and the upper bound capped upstream, so
    # normalising fullGamma here keeps every downstream user consistent.
    if (!is.finite(samplingModelResult$fullGamma) ||
        samplingModelResult$fullGamma <= 0) {
      message(sprintf(
        "  Sampling model returned non-sensical fullGamma (%s); defaulting to fullGamma = partialGamma = %d.",
        as.character(samplingModelResult$fullGamma), partialGamma
      ))
      samplingModelResult$fullGamma <- partialGamma
    }

    # Store result together with the input predictors so the caller can write a
    # full parameter report without needing to re-run the computation.
    gammaRichnessFromSamplingModel <<- cbind(
      samplingModelResult,
      data.frame(
        comp_completeness_method = .ssr_comp_method,
        comp_spread            = comp_spread,
        comp_mean              = comp_mean,
        prop_species_used      = prop_species,
        prop_species_method    = "Chao2",
        prop_species_chao2_Sobs  = chao2_result$S_obs,
        prop_species_chao2_Q1    = chao2_result$Q1,
        prop_species_chao2_Q2    = chao2_result$Q2,
        prop_species_chao2_est   = chao2_result$chao2,
        comp_sd                = comp_sd,
        n_grid_cells_raw       = n_grid_cells,
        sample_mean_range_size = sample_mean_range_size,
        delta_meanlog_used     = delta_ml,
        ratio_sdlog_used       = ratio_sd,
        correction_factor      = samplingModelResult$fullGamma / partialGamma
      )
    )
    message(sprintf(
      "  Sampling model: partialGamma = %d -> fullGamma = %.1f (prop_error = %.4f, 95%% PI: [%.1f, %.1f])",
      partialGamma,
      samplingModelResult$fullGamma,
      samplingModelResult$signed_prop_error,
      samplingModelResult$fullGamma_lwr,
      samplingModelResult$fullGamma_upr
    ))
    message(sprintf(
      "  prop_species (Chao2, regional): %.3f (S_obs = %d, Q1 = %d, Q2 = %d, m = %d grid cells, S_Chao2 = %.1f)",
      prop_species, chao2_result$S_obs, chao2_result$Q1, chao2_result$Q2,
      chao2_result$m, chao2_result$chao2
    ))

    # Apply the proportional correction across all bins uniformly so the simulation
    # produces fullGamma species in total.
    correction_factor <- samplingModelResult$fullGamma / partialGamma
    scaled_plant_distribution <- scaled_plant_distribution |>
      dplyr::mutate(n_scaled = ceiling(n_scaled * correction_factor))

    # Capture error bounds so we can propagate them into the return value below.
    # Checks whether the columns actually exist — guards against older function
    # versions that do not return bounds.
    bound_cols <- c("fullGamma", "fullGamma_lwr", "fullGamma_upr")
    if (all(bound_cols %in% names(samplingModelResult)) &&
        !is.na(samplingModelResult$fullGamma_lwr) &&
        !is.na(samplingModelResult$fullGamma_upr)) {
      sampling_model_bounds <- list(
        fullGamma     = samplingModelResult$fullGamma,
        fullGamma_lwr = samplingModelResult$fullGamma_lwr,
        fullGamma_upr = samplingModelResult$fullGamma_upr
      )
    }

  }
  
  
  # Run simulation
  #options(error = NULL)
  
  grid_sf$centroid <- sf::st_centroid(grid_sf$geometry)      # POINT geometry per grid cell
  # extract numeric centroid coordinates for fast filtering
  cent_coords <- sf::st_coordinates(grid_sf$centroid)         # matrix [n x 2]
  
  # drop polygon geometry from the object used in simulation
  # keeps memory footprint smaller when passing to workers
  grid_centroids <- grid_sf[, c("grid_id")]                  # keep only IDs
  grid_centroids$geometry <- grid_sf$centroid                # replace geometry with centroids
  
  # store coordinates directly on grid_sf for O(1) access
  grid_centroids$cx <- cent_coords[, 1]                               # centroid x
  grid_centroids$cy <- cent_coords[, 2]                               # centroid y
  
  
  grid_sub <- grid_centroids[grid_centroids$grid_id %in% available_grid_ids, ]
  grid_sf$centroid <- NULL  # drop centroid column added above
  rm(grid_centroids, cent_coords); gc()
  
  # Load or build adjacency list (8-connected neighbours)
  # Only needed for contiguous range growth; skipped when simulate_ranges = FALSE
  # (e.g. CI-bound calls that need only species counts).
  if (!simulate_ranges) {
    message("  Skipping adjacency list (simulate_ranges = FALSE).")
    adj_list <- NULL
  } else {
    adj_cache <- file.path(project_root, "output", "shared", "adjacency_list_amazon.rds")
    if (file.exists(adj_cache)) {
      message("  Loading cached adjacency list...")
      adj_saved <- readRDS(adj_cache)
      # Verify it matches the current grid
      if (identical(adj_saved$grid_ids, grid_sub$grid_id) &&
          adj_saved$grid_size_m == grid_size_m) {
        adj_list <- adj_saved$adj_list
        message("  Adjacency list loaded (", adj_saved$n_cells, " cells).")
      } else {
        message("  Cached adjacency list doesn't match current grid. Rebuilding...")
        adj_list <- build_adjacency(grid_sub, grid_size_m)
        # Overwrite the stale cache so later runs (and other organisms) reuse it
        saveRDS(list(adj_list = adj_list, grid_ids = grid_sub$grid_id,
                     n_cells = nrow(grid_sub), grid_size_m = grid_size_m), adj_cache)
        message("  Adjacency list cache updated: ", adj_cache)
      }
      rm(adj_saved)
    } else {
      adj_list <- build_adjacency(grid_sub, grid_size_m)
      # Cache for other organisms
      dir.create(dirname(adj_cache), showWarnings = FALSE, recursive = TRUE)
      saveRDS(list(adj_list = adj_list, grid_ids = grid_sub$grid_id,
                   n_cells = nrow(grid_sub), grid_size_m = grid_size_m), adj_cache)
      message("  Adjacency list cached to: ", adj_cache)
    }
  }
  
  # Build named weight vector keyed by grid_id for frontier weighting
  # (smoothed frontier_weights when supplied; start cells still use sampling_weights)
  if (!is.null(frontier_weights)) stopifnot(length(frontier_weights) == length(available_grid_ids))
  cell_weights <- setNames(if (is.null(frontier_weights)) sampling_weights else frontier_weights,
                           as.character(available_grid_ids))
  
  # Runs the simulation, writes per-chunk records to disk, and exposes the
  # combined occurrence table as the global `simulated_occurrence_dt`. Returns
  # the number of simulated species (central gamma) directly — no geometry join
  # is needed just to count, so the previous per-run `sim_species_sf` build
  # (full occurrence table + polygon geometry, held only to be counted) is gone.
  central_gamma <- simulate_species_records_clustered(
    scaled_plant_distribution = scaled_plant_distribution,
    available_grid_ids = available_grid_ids,
    grid_sf = grid_sub,
    adj_list = adj_list,
    cell_weights = cell_weights,
    sampling_weights = sampling_weights,
    simulate_ranges = simulate_ranges,
    estimateBounds = .amani_estimateBounds,
    parallel = parallel,
    n_cores  = n_cores)

  # ── Return: with or without sampling-model error bounds ────────────────────
  # If extend_richness_from_sampling returned fullGamma_lwr / fullGamma_upr,
  # propagate them into the output by scaling the simulated count proportionally.
  # The simulation was run with the central correction_factor (fullGamma /
  # partialGamma), so the bounds are simply scaled by the same ratio relative to
  # that central correction, preserving the stochastic structure of the run.
  # If bounds are unavailable (sampling extension off, or older model version),
  # fall back to returning a plain integer — original behaviour.
  if (!is.null(sampling_model_bounds)) {
    gamma_lwr <- round(
      central_gamma * (sampling_model_bounds$fullGamma_lwr / sampling_model_bounds$fullGamma)
    )
    gamma_upr <- round(
      central_gamma * (sampling_model_bounds$fullGamma_upr / sampling_model_bounds$fullGamma)
    )
    message(sprintf(
      "  Sampling-model bounds applied to simulated count: gamma = %d, lwr = %d, upr = %d",
      central_gamma, gamma_lwr, gamma_upr
    ))
    c(gamma = central_gamma, gamma_lwr = gamma_lwr, gamma_upr = gamma_upr)
  } else {
    # Default: single integer (original behaviour)
    central_gamma
  }
}

########Testing simulateScaledSpeciesRichness with the plant data to see if it gets the right answer#############
# ------------------------------------------------------------
# only run the estimator sanity check if diagnostics are enabled
# this prevents expensive simulations from running in production
# ------------------------------------------------------------
if (exists("diagnostic") && isTRUE(diagnostic)) {
  
  # message to clearly indicate diagnostic mode is active
  print("checking that the estimator is working by replicating the empirical plant richness")
  
  # show the head of the plant spatial-extent distribution as a quick visual sanity check
  head(plant_distribution)
  
  # initialize an empty list to store simulated species records (not strictly needed,
  # but kept here to mirror the structure of the main simulation logic)
  sim_species_recordsRaw <- list()
  
  # initialize species counter (used internally by the simulator)
  species_counter <- 1
  
  # compute the maximum fitted value from the plant extent model
  # this represents the reference scale of the empirical distribution
  plant_fit_max <- max(fit_pred$fitted)
  
  # define a neutral scaling factor (plant scaled to itself)
  # this should reproduce the original plant richness if the estimator is correct
  extent_scaling_factor <- plant_fit_max / plant_fit_max
  
  # apply scaling to the empirical plant distribution
  # n_scaled should equal n under this neutral scaling
  
  scaled_plant_distribution <- plant_distribution %>%
    mutate(n_scaled = n * extent_scaling_factor)
  
  # run the spatial-extent-based richness simulator
  # this should return a richness value comparable to the empirical data
  sim_speciesRaw_df <- simulateScaledSpeciesRichness(
    scaled_plant_distribution,
    extent_scaling_factor = 1,
    simulate_ranges = simulate_ranges 
  )
  
  # check whether simulated richness matches empirical plant richness
  # unique(unlist()) is defensive in case scientificName is a list-column
  is_equal <- sim_speciesRaw_df ==
    length(unique(unlist(species_counts$scientificName)))
  
  # compose a human-readable diagnostic message
  output_text <- if (is_equal) {
    "The spatial extent based richness estimator is working."
  } else {
    "The spatial extent based richness estimator is NOT working."
  }
  
  # generate a timestamp for reproducible diagnostic outputs
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  
  # construct an informative output filename
  output_filename <- paste0(
    "spatEstEst.sanityCheck.",
    timestamp,
    ".txt"
  )
  
  # build the full output path inside the results directory
  path00 <- file.path(output_folder, output_filename)
  
  # write the diagnostic result to disk for later inspection
  writeLines(output_text, con = path00)
  
}


###Good. when i ran it on April 15 and 17 2025 the numbers are exactly the same.

