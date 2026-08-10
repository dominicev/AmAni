# fillingRanges.R
#
# Defines functions for range-filling interpolation of species occurrence
# data across a spatial grid (elevation-buffered expansion, concave-hull
# fill, sampling-completeness pruning of interior absences).
#
# Sourcing this script only attaches packages and defines functions below —
# it does NOT load any data, download elevation tiles, or run any species
# processing. Call the functions explicitly from another script, e.g.:
#
#   source("scripts/fillingRanges.R")
#   occData <- readRDS("output/species_counts4_filtered.rds")
#
#   covs <- build_grid_cell_covariates(occData)
#   filled_ranges_all <- fill_all_species_ranges(occData, covs$grid_unique_cells)
#
#   # single species
#   one_range <- fill_species_range("Drymonia pendula", occData, covs$grid_unique_cells)

library(sf)
library(dplyr)
library(elevatr)
library(terra)


# ─────────────────────────────────────────────────────────────────────────────
# get_grid_unique_cells(): one row per grid cell (grid_id + geometry), reused
# by both the elevation lookup and the per-species range filling.
# ─────────────────────────────────────────────────────────────────────────────
get_grid_unique_cells <- function(occData) {
  occData |>
    group_by(grid_id) |>
    slice(1) |>
    ungroup() |>
    select(grid_id, geometry)
}

# ─────────────────────────────────────────────────────────────────────────────
# get_grid_elevation(): mean elevation (m) per grid cell via AWS terrain
# tiles. This is purely spatial — independent of which occurrences/species
# are loaded — so it's the one covariate that's safe to cache and reuse
# across runs and across different `occData` inputs (as long as the grid
# cells themselves haven't changed).
#
# Args:
#   grid_unique_cells  sf with one row per grid cell (see get_grid_unique_cells)
#   cache_file         path to cache the elevation lookup; NULL disables caching
#   force              recompute even if a cache file exists
#   verbose            print progress messages
#
# Returns a tibble with columns grid_id, elevation_m.
# ─────────────────────────────────────────────────────────────────────────────
get_grid_elevation <- function(grid_unique_cells,
                                cache_file = "output/grid_elevation_cache.RData",
                                force = FALSE,
                                verbose = TRUE) {

  if (!force && !is.null(cache_file) && file.exists(cache_file)) {
    if (verbose) cat("Loading cached elevation from", cache_file, "\n")
    elevation_lookup <- NULL
    load(cache_file)  # restores `elevation_lookup`
    return(elevation_lookup)
  }

  if (verbose) cat("Computing elevation via AWS terrain tiles...\n")

  # elevatr expects geographic (WGS84) coordinates
  grid_wgs84 <- st_transform(grid_unique_cells, 4326)

  # Download AWS Terrain Tiles; z=5 ≈ 4.8 km resolution — appropriate for
  # Amazon-wide coverage. Increase z for finer resolution (slower download).
  dem <- get_aws_terrain(grid_wgs84, z = 5, prj = st_crs(grid_wgs84)$wkt)

  # Mean elevation (metres) within each grid cell polygon
  elev_extracted <- terra::extract(dem, vect(grid_wgs84), fun = mean, na.rm = TRUE)

  elevation_lookup <- tibble::tibble(
    grid_id     = grid_unique_cells$grid_id,
    elevation_m = elev_extracted[, 2]
  )

  if (!is.null(cache_file)) {
    dir.create(dirname(cache_file), showWarnings = FALSE, recursive = TRUE)
    save(elevation_lookup, file = cache_file)
    if (verbose) cat("Saved elevation cache to", cache_file, "\n")
  }

  elevation_lookup
}


# ─────────────────────────────────────────────────────────────────────────────
# build_grid_cell_covariates(): assemble grid_unique_cells with both
# covariates used by fill_species_range() — elevation and sampling
# completeness.
#
# Args:
#   occData               occurrence sf with columns scientificName, grid_id,
#                         count, geometry
#   elevation_cache_file  cache path for the (species-independent) elevation
#                         lookup; NULL disables caching
#   force_elevation       recompute elevation even if a cache exists
#   completeness          how to obtain sampling completeness:
#                           "compute"  (default) — calculate fresh from
#                                      `occData`, so it always reflects
#                                      whichever occurrences were passed in
#                           a file path (character) — load a previously-saved
#                                      completeness table (.rds) instead of
#                                      recomputing; only appropriate if you're
#                                      confident it matches this occData
#                           a data frame — use a completeness table already
#                                      in memory (e.g. reuse across calls)
#   verbose               print progress / summaries
#
# Returns a list with:
#   grid_unique_cells  sf with grid_id, geometry, elevation_m, n_records,
#                      completeness
#   cell_completeness  the completeness table actually used (handy to save
#                      or pass into a later call via `completeness =`)
# ─────────────────────────────────────────────────────────────────────────────
build_grid_cell_covariates <- function(occData,
                                        elevation_cache_file = "output/grid_elevation_cache.RData",
                                        force_elevation = FALSE,
                                        completeness = "compute",
                                        verbose = TRUE) {

  grid_unique_cells <- get_grid_unique_cells(occData)

  elevation_lookup <- get_grid_elevation(
    grid_unique_cells,
    cache_file = elevation_cache_file,
    force = force_elevation,
    verbose = verbose
  )
  grid_unique_cells <- grid_unique_cells |>
    left_join(elevation_lookup, by = "grid_id")

  cell_completeness <- if (identical(completeness, "compute")) {
    if (verbose) cat("Computing sampling completeness from occData...\n")
    compute_cell_completeness(occData)
  } else if (is.character(completeness) && length(completeness) == 1) {
    if (verbose) cat("Loading completeness from", completeness, "\n")
    readRDS(completeness)
  } else if (is.data.frame(completeness)) {
    if (verbose) cat("Using supplied completeness table\n")
    completeness
  } else {
    stop("`completeness` must be \"compute\", a path to a cached .rds file, ",
         "or a data frame with columns grid_id, n_records, completeness")
  }

  grid_unique_cells <- grid_unique_cells |>
    left_join(cell_completeness |> select(grid_id, n_records, completeness),
              by = "grid_id")

  if (verbose) {
    cat("Elevation summary (m):\n");        print(summary(grid_unique_cells$elevation_m))
    cat("Sampling completeness summary:\n"); print(summary(grid_unique_cells$completeness))
  }

  # left_join() can drop the sf class; re-assert it so downstream spatial ops work
  grid_unique_cells <- st_as_sf(grid_unique_cells)

  list(grid_unique_cells = grid_unique_cells, cell_completeness = cell_completeness)
}

# ─────────────────────────────────────────────────────────────────────────────
# build_grid_context(): precompute, once, all the grid-level structures that
# every species reuses. Building these here (instead of inside the per-species
# loop) is the main single-threaded speedup: the fixed grid topology and
# centroids are computed a single time, and per-species work becomes cheap
# integer-index bookkeeping.
#
# Args:
#   grid_unique_cells  sf with grid_id, geometry, elevation_m, completeness
#                      (see build_grid_cell_covariates())
#
# Returns a list with:
#   grid_unique_cells  the input sf (kept for building the final output)
#   n                  number of cells
#   grid_id            plain grid_id vector, aligned to rows 1..n
#   id_to_index        named integer vector mapping grid_id -> row index
#   elevation_m        plain numeric vector of elevation, aligned to rows
#   completeness       plain numeric vector of completeness, aligned to rows
#   adj                st_touches adjacency list over all cells (sgbp)
#   centroids          sf of cell centroids
#   cxy                matrix of centroid coordinates (for bbox pre-filtering)
#   corners            list (length n) of each cell's boundary vertices as a
#                      coordinate matrix — used to build the presence hull from
#                      cell corner points (equivalent to hulling the union of
#                      the cell polygons, but far cheaper)
#   crs                CRS of the grid, for reconstructing hull-input geometry
# ─────────────────────────────────────────────────────────────────────────────
build_grid_context <- function(grid_unique_cells) {
  if (!inherits(grid_unique_cells, "sf")) grid_unique_cells <- st_as_sf(grid_unique_cells)
  n         <- nrow(grid_unique_cells)
  grid_id   <- grid_unique_cells$grid_id
  centroids <- suppressWarnings(st_centroid(grid_unique_cells))

  # Boundary vertices of every cell, split by feature (one st_coordinates call).
  cc      <- st_coordinates(st_geometry(grid_unique_cells))
  feat    <- cc[, ncol(cc)]
  corners <- lapply(split(as.data.frame(cc[, c("X", "Y")]), feat), as.matrix)

  list(
    grid_unique_cells = grid_unique_cells,
    n                 = n,
    grid_id           = grid_id,
    id_to_index       = stats::setNames(seq_len(n), as.character(grid_id)),
    elevation_m       = grid_unique_cells$elevation_m,
    completeness      = grid_unique_cells$completeness,
    adj               = st_touches(grid_unique_cells),
    centroids         = centroids,
    cxy               = st_coordinates(centroids),
    corners           = corners,
    crs               = st_crs(grid_unique_cells)
  )
}

# ─────────────────────────────────────────────────────────────────────────────
# fill_species_range(): build a filled range for a single species.
#
# For one scientificName, this (1) finds the observed grid cells, (2) expands
# into edge-touching cells within an elevation buffer, (3) draws a concave hull
# over observed + buffer cells, (4) fills every grid cell whose centroid falls
# inside that hull and within the elevation band, and (5) prunes well-sampled
# interior cells where the species was never recorded (reliable absences).
#
# Args:
#   species                 scientificName to process (character, length 1)
#   occData                 occurrence sf with columns scientificName, grid_id
#   grid_unique_cells       sf of unique grid cells w/ elevation_m + completeness
#                           (see build_grid_cell_covariates()); only needed when
#                           `grid_context` is not supplied
#   grid_context            optional precomputed grid context from
#                           build_grid_context(); pass once when looping over
#                           many species to avoid recomputing adjacency and
#                           centroids each call. If NULL it is built from
#                           `grid_unique_cells`.
#   cell_idx                optional integer indices (into grid_context) of the
#                           species' observed cells; when supplied, `occData`
#                           is not scanned. Used by fill_all_species_ranges().
#   elev_buffer             metres added above/below the observed elevation range
#   hull_ratio              st_concave_hull ratio (0 = tightest, 1 = convex)
#   completeness_threshold  cells >= this completeness lacking the species are
#                           dropped as reliable absences
#   verbose                 print a one-line per-species summary
#
# Returns an sf of filled grid cells tagged with `scientificName`, `is_observed`,
# and `source` ("observed"/"filled"). If the species has too few cells (or no
# finite elevation data) to build a hull, returns its observed cells unchanged
# (same cells as the input), so no species is lost. Returns NULL only when the
# species has no observed cells at all.
# ─────────────────────────────────────────────────────────────────────────────
fill_species_range <- function(species,
                                occData = NULL,
                                grid_unique_cells = NULL,
                                grid_context = NULL,
                                cell_idx = NULL,
                                elev_buffer = 100,
                                hull_ratio = 0.45,
                                completeness_threshold = 0.5,
                                species_col = NULL,
                                verbose = FALSE) {

  # Precomputed grid structures (adjacency, centroids, covariate vectors).
  # Build them here only when this is a standalone single-species call.
  if (is.null(grid_context)) {
    if (is.null(grid_unique_cells)) {
      stop("Supply either `grid_context` or `grid_unique_cells`.")
    }
    grid_context <- build_grid_context(grid_unique_cells)
  }
  ctx <- grid_context
  guc <- ctx$grid_unique_cells

  # Integer indices (into ctx) of this species' observed cells. When cell_idx
  # is supplied (the fast path used by fill_all_species_ranges) we skip the
  # full scan of occData entirely.
  if (is.null(cell_idx)) {
    if (is.null(occData)) {
      stop("Supply either `cell_idx` or `occData` to locate the species' cells.")
    }
    occ_tbl     <- st_drop_geometry(occData)
    species_col <- detect_species_col(occ_tbl, species_col, "fill_species_range()")
    sp_ids <- occ_tbl |>
      filter(.data[[species_col]] == species) |>
      pull(grid_id) |>
      unique()
    cell_idx <- ctx$id_to_index[as.character(sp_ids)]
  }
  cell_idx <- sort(unique(cell_idx[!is.na(cell_idx)]))

  if (length(cell_idx) == 0) {
    if (verbose) cat(species, "| no observed cells; skipping\n")
    return(NULL)
  }

  # Fallback returned whenever a hull can't be built: the species' observed
  # cells unchanged (same cells as the input), tagged like the normal output.
  observed_only <- guc[cell_idx, ] |>
    mutate(scientificName = species, is_observed = TRUE, source = "observed")

  elev <- ctx$elevation_m
  species_elevation_range <- range(elev[cell_idx], na.rm = TRUE)

  if (!all(is.finite(species_elevation_range))) {
    if (verbose) cat(species, "| no finite elevation data; returning observed cells only\n")
    return(observed_only)
  }

  # Elevation buffer: within `elev_buffer` m of the observed min/max elevation
  elev_low  <- species_elevation_range[1] - elev_buffer
  elev_high <- species_elevation_range[2] + elev_buffer

  # Cells touching any species cell (shared edge/vertex), via the precomputed
  # adjacency list — no per-species spatial predicate.
  touch_idx <- sort(unique(unlist(ctx$adj[cell_idx])))

  # Buffer cells: touching, within the elevation band, not already observed.
  in_band    <- !is.na(elev[touch_idx]) &
                elev[touch_idx] >= elev_low & elev[touch_idx] <= elev_high
  buffer_idx <- setdiff(touch_idx[in_band], cell_idx)

  expanded_idx <- union(cell_idx, buffer_idx)

  # Concave hull over the corner points of the observed + buffer cells.
  # Hulling the cell corners (rather than unioning the cell polygons) is much
  # cheaper and reproduces the union-based hull exactly on our grid, since the
  # union boundary vertices are a subset of these corner points.
  hull_pts <- st_sfc(
    st_multipoint(do.call(rbind, ctx$corners[expanded_idx])),
    crs = ctx$crs
  )
  presence_hull <- tryCatch(
    st_concave_hull(hull_pts, ratio = hull_ratio, allow_holes = FALSE),
    error = function(e) NULL
  )

  if (is.null(presence_hull) || length(presence_hull) == 0 ||
      all(st_is_empty(presence_hull))) {
    if (verbose) cat(species, "| too few cells for a hull; returning observed cells only\n")
    return(observed_only)
  }

  # Grid cells whose centroid falls inside the hull. Restrict the (expensive)
  # st_within test to centroids inside the hull's bounding box first.
  bb   <- st_bbox(presence_hull)
  cxy  <- ctx$cxy
  cand <- which(cxy[, 1] >= bb[["xmin"]] & cxy[, 1] <= bb[["xmax"]] &
                cxy[, 2] >= bb[["ymin"]] & cxy[, 2] <= bb[["ymax"]])

  if (length(cand) == 0) {
    inside_idx <- integer(0)
  } else {
    within_hit <- st_within(ctx$centroids[cand, ], presence_hull,
                            sparse = FALSE)[, 1]
    inside_idx <- cand[within_hit]
  }

  # Same elevation filter used for the buffer expansion
  band       <- !is.na(elev[inside_idx]) &
                elev[inside_idx] >= elev_low & elev[inside_idx] <= elev_high
  filled_idx <- inside_idx[band]

  # Sampling-completeness pruning of interior cells.
  # A hull-filled cell that does NOT contain the species is a "candidate absence".
  # If that cell is well sampled (completeness >= threshold) yet the species was
  # not recorded, we treat it as a reliable true absence and drop it from the hull.
  # If completeness is low (or unknown), the absence is unreliable, so we keep it.
  # Observed species cells are never dropped.
  comp            <- ctx$completeness
  has_species     <- filled_idx %in% cell_idx
  well_sampled    <- !is.na(comp[filled_idx]) & comp[filled_idx] >= completeness_threshold
  drop_as_absence <- !has_species & well_sampled
  n_dropped       <- sum(drop_as_absence)

  # Always retain observed cells (guarantees no species is lost even if a hull
  # vertex centroid falls exactly on the boundary), then build the output once.
  keep <- sort(union(filled_idx[!drop_as_absence], cell_idx))
  keep_observed <- keep %in% cell_idx

  filled_cells <- guc[keep, ] |>
    mutate(
      scientificName = species,
      is_observed    = keep_observed,
      source         = if_else(keep_observed, "observed", "filled")
    )

  if (verbose) {
    cat(species,
        "| observed:", length(cell_idx),
        "| buffer added:", length(buffer_idx),
        "| filled:", nrow(filled_cells),
        "| dropped absences:", n_dropped,
        "| hull ratio:", hull_ratio, "\n")
  }

  filled_cells
}

# ─────────────────────────────────────────────────────────────────────────────
# fill_all_species_ranges(): apply fill_species_range() across every
# scientificName present in occData and combine the results into one sf.
#
# Args:
#   occData             occurrence sf with columns scientificName, grid_id
#   grid_unique_cells   sf w/ elevation_m + completeness (see
#                       build_grid_cell_covariates())
#   ...                 additional arguments forwarded to fill_species_range()
#                       (elev_buffer, hull_ratio, completeness_threshold)
#   parallel            if TRUE, process species in parallel via furrr/future
#                       on a `multisession` plan (works on Windows). The cluster
#                       is created with parallelly::makeClusterPSOCK(), which is
#                       far more robust than base parallel::makeCluster() on
#                       network-synced paths (OneDrive) and locked-down hosts.
#   n_cores             number of worker processes when parallel = TRUE
#   verbose             print per-species progress (serial) or a furrr progress
#                       bar plus a coarse summary (parallel)
#
# Returns one sf combining every species' filled (or observed-only) cells.
# ─────────────────────────────────────────────────────────────────────────────
fill_all_species_ranges <- function(occData, grid_unique_cells, ...,
                                     species_col = NULL,
                                     parallel = FALSE,
                                     n_cores = max(1L, parallelly::availableCores() - 1L),
                                     verbose = TRUE) {
  # Precompute all grid-level structures once; every species reuses them.
  ctx <- build_grid_context(grid_unique_cells)

  # Map each species to the integer indices (into ctx) of its observed cells,
  # computed a single time. This replaces the per-species scan of occData.
  occ_tbl     <- st_drop_geometry(occData)
  species_col <- detect_species_col(occ_tbl, species_col, "fill_all_species_ranges()")
  occ_tab <- occ_tbl |>
    distinct(.data[[species_col]], grid_id)
  occ_tab$idx <- ctx$id_to_index[as.character(occ_tab$grid_id)]
  species_to_cells <- split(occ_tab$idx, occ_tab[[species_col]])

  all_species <- sort(names(species_to_cells))
  dots <- list(...)

  # Per-species worker. Kept as a single closure so the serial and parallel
  # paths run identical code; furrr auto-exports the globals it references
  # (fill_species_range, build_grid_context, ctx, species_to_cells, dots) and
  # the packages they need (sf, dplyr), so there is no manual clusterExport to
  # keep in sync.
  process_one <- function(sp) {
    tryCatch(
      do.call(fill_species_range, c(
        list(species = sp, grid_context = ctx,
             cell_idx = species_to_cells[[sp]], verbose = FALSE),
        dots)),
      error = function(e) NULL
    )
  }

  # Parallelise with furrr/future on a `multisession` plan. Unlike a bare
  # parallel::makeCluster() PSOCK cluster (which fails opaquely with
  # "error reading from connection" on OneDrive-synced paths and locked-down
  # hosts), future's multisession backend uses parallelly::makeClusterPSOCK(),
  # which is specifically hardened against those setup/connection problems.
  # If cluster setup still fails for any reason we fall back to serial rather
  # than aborting the whole run.
  filled_ranges <- if (parallel && n_cores > 1L) {
    tryCatch({
      if (!requireNamespace("furrr", quietly = TRUE) ||
          !requireNamespace("future", quietly = TRUE)) {
        stop("packages 'furrr' and 'future' are required for parallel = TRUE")
      }
      if (verbose) cat("Processing", length(all_species), "species on",
                       n_cores, "multisession workers...\n")

      old_plan <- future::plan(future::multisession, workers = n_cores)
      on.exit(future::plan(old_plan), add = TRUE)

      furrr::future_map(
        all_species,
        process_one,
        .options  = furrr::furrr_options(seed = TRUE),
        .progress = isTRUE(verbose)
      )
    }, error = function(e) {
      warning("fill_all_species_ranges(): parallel backend failed (",
              conditionMessage(e), "); falling back to serial processing.",
              call. = FALSE)
      NULL
    })
  } else {
    NULL
  }

  if (is.null(filled_ranges)) {
    filled_ranges <- lapply(all_species, function(sp) {
      out <- process_one(sp)
      if (verbose && is.null(out)) cat(sp, "| ERROR or no output\n")
      out
    })
  }
  names(filled_ranges) <- all_species

  # Drop species that could not be processed at all, then combine into one sf
  filled_ranges     <- filled_ranges[!vapply(filled_ranges, is.null, logical(1))]
  filled_ranges_all <- bind_rows(filled_ranges)

  if (verbose) {
    cat("Filled ranges for", length(filled_ranges), "of", length(all_species),
        "species |", nrow(filled_ranges_all), "total cells\n")
  }

  filled_ranges_all
}

# ── Example usage (reference only — not executed on source) ─────────────────
# library(rstudioapi); library(ggplot2); library(rnaturalearth)
#
# if (exists(".amani_wd")) {
#   setwd(.amani_wd)
# } else {
#   setwd(dirname(dirname(rstudioapi::getActiveDocumentContext()$path)))
# }
#
# occData <- readRDS("output/species_counts4_filtered.rds")
#
# # Completeness is always computed from `occData` here (the default), so it
# # reflects whichever occurrence table you pass in.
# covs <- build_grid_cell_covariates(occData)
#
# filled_ranges_all <- fill_all_species_ranges(occData, covs$grid_unique_cells)
#
# chosen_species <- "Drymonia pendula"
# one_range <- fill_species_range(chosen_species, occData, covs$grid_unique_cells,
#                                 verbose = TRUE)
#
# dataDir <- paste0(dirname(getwd()), "/0-CuratingPlantData/output")
# amazon_shp_path <- paste0(dirname(dataDir),
#                           "/data/amazon_shapefile/amazon_sensulatissimo_gmm_v1.shp")
# amazon_shape <- st_read(amazon_shp_path, quiet = TRUE)
# amazon_shape <- st_simplify(st_make_valid(amazon_shape), dTolerance = 0.01)
# amazon_3857  <- st_transform(amazon_shape, 3857)
#
# ggplot() +
#   geom_sf(data = amazon_3857, fill = "grey90", color = "grey60", linewidth = 0.3) +
#   geom_sf(data = st_transform(one_range, 3857), aes(fill = source), color = NA) +
#   scale_fill_manual(values = c(observed = "#1a1d1f", filled = "#a8d5b5")) +
#   labs(title = bquote(italic(.(chosen_species)) ~ "— concave hull fill")) +
#   theme_minimal()
