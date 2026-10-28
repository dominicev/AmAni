# fillGridHoles.R
#
# The analysis grid (endemismData) only contains cells that have plant records,
# so the sampled Amazon grid is full of holes (~1/3 of the cells inside the
# basin) that reflect collecting effort rather than habitat. Range growth in
# simulate_species_chunk() walks 8-connected neighbours, so these holes make
# simulated ranges porous and fragment them.
#
# fill_grid_holes() adds the missing cells back: every lattice position whose
# square intersects the Amazon boundary but is absent from grid_sf. New cells
# reuse the existing grid_id scheme (id = (row - 1) * stride + col), so they
# are consistent with Stage 0 ids and cannot collide with existing ones.
# Attribute columns of the new rows are NA, and a logical `grid_filled` column
# marks them (FALSE for original cells).

fill_grid_holes <- function(grid_sf, boundary, grid_size_m, tol_m = 1) {
  stopifnot(inherits(grid_sf, "sf"), "grid_id" %in% names(grid_sf))

  crs_grid <- sf::st_crs(grid_sf)
  boundary <- sf::st_union(sf::st_geometry(sf::st_transform(boundary, crs_grid)))

  ctr <- sf::st_coordinates(sf::st_centroid(sf::st_geometry(grid_sf)))
  x0  <- min(ctr[, 1]); y0 <- min(ctr[, 2])
  col <- as.integer(round((ctr[, 1] - x0) / grid_size_m)) + 1L
  row <- as.integer(round((ctr[, 2] - y0) / grid_size_m)) + 1L

  resid <- max(abs(ctr[, 1] - (x0 + (col - 1) * grid_size_m)),
               abs(ctr[, 2] - (y0 + (row - 1) * grid_size_m)))
  if (resid > tol_m) {
    stop("fill_grid_holes(): centroids do not sit on a regular lattice with ",
         "spacing grid_size_m = ", grid_size_m, " (max residual ", signif(resid, 3), " m).")
  }

  # Recover the id scheme: id = (row - 1) * stride + col
  ids    <- grid_sf$grid_id
  stride <- unique((ids[row > 1L] - col[row > 1L]) / (row[row > 1L] - 1L))
  if (length(stride) != 1L || stride != round(stride) ||
      !all(ids == (row - 1L) * as.integer(stride) + col)) {
    stop("fill_grid_holes(): grid_id is not (row - 1) * stride + col on this lattice; ",
         "cannot generate consistent ids for the missing cells.")
  }
  stride <- as.integer(stride)
  if (max(col) > stride) stop("fill_grid_holes(): lattice wider than the inferred id stride.")

  # Candidate positions: every lattice position inside the existing extent
  cand <- expand.grid(col = seq_len(max(col)), row = seq_len(max(row)))
  cand$grid_id <- (cand$row - 1L) * stride + cand$col
  cand <- cand[!cand$grid_id %in% ids, ]

  cand_xy  <- cbind(x0 + (cand$col - 1) * grid_size_m, y0 + (cand$row - 1) * grid_size_m)
  cand_pts <- sf::st_as_sf(data.frame(x = cand_xy[, 1], y = cand_xy[, 2]),
                           coords = c("x", "y"), crs = crs_grid)
  cand_sq  <- sf::st_geometry(sf::st_buffer(cand_pts, dist = grid_size_m / 2,
                                            endCapStyle = "SQUARE"))
  keep <- lengths(sf::st_intersects(cand_sq, boundary)) > 0
  n_new <- sum(keep)

  grid_sf$grid_filled <- FALSE
  if (n_new == 0L) return(grid_sf)

  new_rows <- grid_sf[rep(NA_integer_, n_new), ]
  new_rows$grid_id     <- cand$grid_id[keep]
  new_rows$grid_filled <- TRUE
  sf::st_geometry(new_rows) <- cand_sq[keep]

  message(sprintf("  fill_grid_holes(): added %d missing cells inside the boundary (%d -> %d cells).",
                  n_new, nrow(grid_sf), nrow(grid_sf) + n_new))
  rbind(grid_sf, new_rows)
}
