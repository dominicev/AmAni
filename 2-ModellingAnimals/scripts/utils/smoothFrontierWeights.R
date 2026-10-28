# smoothFrontierWeights.R
#
# Spatially smoothed plant-occurrence weights for the frontier step of
# contiguous range growth (simulate_species_chunk()). Raw per-cell plant
# record counts are sharply bimodal, so frontier growth follows corridors of
# well-sampled cells and produces lobed, gappy ranges. Smoothing the surface
# with a Gaussian kernel keeps the large-scale sampling signal while removing
# cell-to-cell spikes.
#
# The weights depend only on the plant data (filled vs unfilled) and the fixed
# Amazon grid, so they are cached in output/shared/ and reused across organisms.


# Gaussian-kernel smoothing on a regular lattice via normalized convolution.
# Cells are mapped to integer (row, col) positions from their centroids; the
# kernel is applied by summing shifted copies of the weight matrix. Dividing by
# the kernel mass over *available* cells only means cells at the Amazon
# boundary (or next to holes in the grid) are not pulled toward zero.
#
# `observed` (default: all TRUE) marks cells whose weight is real data. Cells
# with observed = FALSE contribute nothing to the kernel sums but still receive
# a smoothed value, i.e. they are interpolated from observed neighbours (used
# for grid cells added back by fill_grid_holes()). Returns NaN for a cell with
# no observed cell inside the kernel radius.
smooth_weights_on_grid <- function(weights, cx, cy, grid_size_m,
                                   bandwidth_cells = 2,
                                   radius_cells = ceiling(3 * bandwidth_cells),
                                   observed = rep(TRUE, length(weights))) {
  stopifnot(length(weights) == length(cx), length(cx) == length(cy),
            length(observed) == length(weights),
            !anyNA(weights), bandwidth_cells > 0)

  col <- as.integer(round((cx - min(cx)) / grid_size_m)) + 1L
  row <- as.integer(round((cy - min(cy)) / grid_size_m)) + 1L
  if (anyDuplicated(data.frame(row, col))) {
    stop("smooth_weights_on_grid(): centroids do not map to a regular lattice ",
         "with spacing grid_size_m = ", grid_size_m)
  }

  nr  <- max(row); nc <- max(col)
  idx <- cbind(row, col)
  W <- matrix(0, nr, nc); W[idx] <- weights * observed
  M <- matrix(0, nr, nc); M[idx] <- as.numeric(observed)

  # S[i, j] <- X[i + dr, j + dc], zero outside the matrix
  shift <- function(X, dr, dc) {
    S <- matrix(0, nr, nc)
    src_r <- max(1L, 1L + dr):min(nr, nr + dr)
    src_c <- max(1L, 1L + dc):min(nc, nc + dc)
    if (length(src_r) && length(src_c) && min(src_r) <= max(src_r) && min(src_c) <= max(src_c)) {
      S[src_r - dr, src_c - dc] <- X[src_r, src_c]
    }
    S
  }

  num <- matrix(0, nr, nc)
  den <- matrix(0, nr, nc)
  offsets <- -radius_cells:radius_cells
  for (dr in offsets) {
    for (dc in offsets) {
      d2 <- dr^2 + dc^2
      if (d2 > radius_cells^2) next          # circular kernel footprint
      k   <- exp(-d2 / (2 * bandwidth_cells^2))
      num <- num + k * shift(W, dr, dc)
      den <- den + k * shift(M, dr, dc)
    }
  }

  num[idx] / den[idx]
}


# Replace the weights of filled (unsampled) grid cells with a kernel-weighted
# mean of nearby observed cells. These cells are absent from the plant data
# because of collecting effort, not habitat, so giving them the floor weight
# would recreate the holes this is meant to remove. Observed weights are left
# unchanged. Any filled cell with no observed neighbour in range falls back to
# the median observed weight.
interpolate_filled_weights <- function(weights, grid_ids, grid_sf, grid_size_m,
                                       filled, bandwidth_cells = 1,
                                       radius_cells = ceiling(3 * bandwidth_cells)) {
  stopifnot(length(weights) == length(grid_ids), length(filled) == length(grid_ids))
  if (!any(filled)) return(weights)

  geom <- sf::st_geometry(grid_sf)[match(grid_ids, grid_sf$grid_id)]
  xy   <- sf::st_coordinates(sf::st_centroid(geom))

  smoothed <- smooth_weights_on_grid(
    weights, xy[, 1], xy[, 2], grid_size_m,
    bandwidth_cells = bandwidth_cells, radius_cells = radius_cells,
    observed = !filled
  )
  fallback <- !is.finite(smoothed[filled])
  if (any(fallback)) {
    message("  interpolate_filled_weights(): ", sum(fallback),
            " filled cell(s) had no observed neighbour in range; using median weight.")
    smoothed[filled][fallback] <- stats::median(weights[!filled])
  }
  weights[filled] <- smoothed[filled]
  weights
}


# Load smoothed frontier weights from cache, or build and cache them.
# The cache is reused only if the grid ids, raw weights, and smoothing
# parameters all match; otherwise it is rebuilt. raw_weights is compared in
# full (it is ~10k numbers) so any change to the plant data invalidates it.
load_or_build_frontier_weights <- function(raw_weights,
                                           grid_ids,
                                           grid_sf,
                                           grid_size_m,
                                           fill_plants,
                                           bandwidth_cells = 2,
                                           radius_cells = ceiling(3 * bandwidth_cells),
                                           cache_dir = file.path(project_root, "output", "shared"),
                                           redo = FALSE) {
  stopifnot(length(raw_weights) == length(grid_ids))

  cache_file <- file.path(cache_dir, sprintf(
    "frontier_weights_%s_bw%s_r%s.rds",
    if (isTRUE(fill_plants)) "filledPlants" else "rawPlants",
    format(bandwidth_cells), format(radius_cells)
  ))

  if (file.exists(cache_file) && !isTRUE(redo)) {
    cached <- readRDS(cache_file)
    if (identical(cached$grid_ids, grid_ids) &&
        identical(cached$raw_weights, raw_weights) &&
        identical(cached$grid_size_m, grid_size_m)) {
      message("  Loaded cached smoothed frontier weights: ", cache_file)
      return(cached$weights)
    }
    message("  Cached frontier weights don't match current grid/plant data. Rebuilding...")
  } else if (isTRUE(redo)) {
    message("  Rebuilding smoothed frontier weights (redo = TRUE).")
  } else {
    message("  Building smoothed frontier weights for the first time.")
  }

  geom <- sf::st_geometry(grid_sf)[match(grid_ids, grid_sf$grid_id)]
  if (anyNA(match(grid_ids, grid_sf$grid_id))) {
    stop("load_or_build_frontier_weights(): some grid_ids are missing from grid_sf")
  }
  xy <- sf::st_coordinates(sf::st_centroid(geom))

  weights <- smooth_weights_on_grid(
    raw_weights, xy[, 1], xy[, 2], grid_size_m,
    bandwidth_cells = bandwidth_cells, radius_cells = radius_cells
  )

  dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
  saveRDS(list(weights = weights, raw_weights = raw_weights, grid_ids = grid_ids,
               grid_size_m = grid_size_m, bandwidth_cells = bandwidth_cells,
               radius_cells = radius_cells, fill_plants = isTRUE(fill_plants),
               created = Sys.time()),
          cache_file)
  message("  Smoothed frontier weights cached to: ", cache_file)
  weights
}
