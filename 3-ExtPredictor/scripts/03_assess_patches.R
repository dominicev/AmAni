
#library(knitr)
#library(rstudioapi)
#knitr::opts_knit$set(root.dir = dirname(rstudioapi::getActiveDocumentContext()$path))
#setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Set working directory to the project root.
# When run via Rscript, use the --file= argument to locate the script,
# then navigate one level up from scripts/ to the project root.
# When run interactively in RStudio, use rstudioapi instead.
args <- commandArgs(trailingOnly = FALSE)
script_arg <- args[grep("^--file=", args)]

if (length(script_arg) > 0) {
  # Running via Rscript: derive project root from script path
  script_path <- sub("^--file=", "", script_arg)
  setwd(file.path(dirname(normalizePath(script_path)), ".."))
} else if (requireNamespace("rstudioapi", quietly = TRUE) &&
           rstudioapi::isAvailable()) {
  # Running interactively in RStudio: derive from active document
  setwd(file.path(dirname(rstudioapi::getActiveDocumentContext()$path), ".."))
}

message("Working directory set to: ", getwd())


#############Importing patch stats#################
# define the directory where output CSV files live (relative to current working directory)
importDir <- "output/deforestTileOutStats/"                     # folder containing per-seed outputs

getwd()
# list all files in the directory that end with .csv
# full.names = TRUE returns absolute/relative paths, which is safer for reading
# Only match w-prefixed CSVs from the deterministic grid run.
# The output directory also contains ~487 non-w-prefixed CSVs from an
# earlier random sampling run; those have different spatial sampling
# properties and would contaminate the deterministic results.
csv_files <- list.files(
  path       = importDir,                               # directory to scan
  pattern    = "^w[0-9]+\\.patch_metrics\\.csv$",       # deterministic outputs only
  full.names = TRUE                                      # return paths, not just filenames
)

# sanity check: how many CSV files were found
length(csv_files)                                        # total number of CSVs detected

# take only the first 10 CSV files as a test subset
# head() preserves ordering as returned by list.files()
csv_test <- csv_files
#csv_test <-head(csv_files, 500)                          # vector of first N CSV paths, or all paths if you don't want to test

# optional: inspect which files will be read
basename(csv_test)                                       # just the filenames, for clarity

# read CSV files into R, tagging each row with its source window_id
# the window_id is extracted from the filename (e.g. "w1234.patch_metrics.csv" -> 1234)
# and will be used downstream for boundary-based deduplication
csv_data_list <- lapply(csv_test, function(f) {
  d <- read.csv(f)
  d$window_id <- as.integer(gsub("^w|\\.patch_metrics\\.csv$", "", basename(f)))
  d
})

#name the list elements by file name (without directory)
names(csv_data_list) <- basename(csv_test)               # easier downstream referencing

# example check: look at the structure of the first imported table
str(csv_data_list[[1]])                                  # verify columns and types

#############Combining patch stats#################

all_patches <- do.call(
  rbind,                                   # row-bind all elements
  csv_data_list                            # list of patch-metric data.frames
)

head(all_patches)
nrow(all_patches)


#############################################################################
# TRUNCATION DIAGNOSTIC — do large patches hit the window boundary?
#############################################################################
# Before any deduplication, check whether the largest patches are physically
# truncated by the extraction window.  For each of the top N patches (by area),
# load the actual polygon from the .gpkg file and test whether it intersects
# the window's boundary rectangle (derived from the .tif extent).
#
# RESULT (2026-03-20, run on 9,376 boundary-filtered patches):
#   Of the 20 largest patches, only 2 touch their window edge:
#     - Window 2834, patch 3405 (2,707 km²)
#     - Window 1869, patch 2522 (2,674 km²)
#   The remaining 18 (including the largest at 2,911 km²) are fully
#   contained.  Window dimensions (~120 × 120 km) are large enough to
#   contain nearly all deforestation patches.
#
#   Conclusion: truncation affects only a handful of the very largest
#   patches (<0.1% of the dataset) and is negligible for distribution
#   fitting.  A polygon-merge step across adjacent windows is NOT
#   warranted for the current analysis.
#############################################################################

library(terra)
library(sf)

# Read TIF extents for every window that produced output
tif_files_trunc <- list.files(
  "output/deforestTileOutStats/",
  pattern = "^w[0-9]+\\.patches_filtered\\.tif$",
  full.names = TRUE
)

tif_ext_trunc <- lapply(tif_files_trunc, function(f) {
  wid <- as.integer(gsub("^w|\\.patches_filtered\\.tif$", "", basename(f)))
  r <- rast(f)
  e <- ext(r)
  rm(r)
  data.frame(window_id = wid,
             tif_xmin = e$xmin, tif_xmax = e$xmax,
             tif_ymin = e$ymin, tif_ymax = e$ymax)
})
tif_ext_trunc <- do.call(rbind, tif_ext_trunc)
gc(verbose = FALSE)

# Identify the N largest patches (by area) across all windows
n_check <- 20
top_patches <- all_patches[order(-all_patches$area_km2), ]
top_patches <- top_patches[!duplicated(top_patches[, c("centroid_x", "centroid_y")]), ]
top_patches <- head(top_patches, n_check)

# For each, load the polygon and test for edge contact
trunc_results <- lapply(seq_len(nrow(top_patches)), function(i) {
  wid <- top_patches$window_id[i]
  pid <- top_patches$patch_id[i]

  gpkg_path <- file.path("output/deforestTileOutStats",
                          paste0("w", wid, ".patch_polygons.gpkg"))
  if (!file.exists(gpkg_path)) return(NULL)

  polys <- st_read(gpkg_path, quiet = TRUE)
  # gpkg stores the patch label in a column called "patches"
  patch_poly <- polys[polys$patches == pid, ]
  if (nrow(patch_poly) == 0) return(NULL)

  # Build the window boundary rectangle from the TIF extent
  win_ext <- tif_ext_trunc[tif_ext_trunc$window_id == wid, ]
  if (nrow(win_ext) == 0) return(NULL)

  win_box <- st_as_sfc(st_bbox(c(
    xmin = win_ext$tif_xmin, xmax = win_ext$tif_xmax,
    ymin = win_ext$tif_ymin, ymax = win_ext$tif_ymax
  ), crs = st_crs(patch_poly)))
  win_boundary <- st_boundary(win_box)

  touches <- st_intersects(patch_poly, win_boundary, sparse = FALSE)[1, 1]

  data.frame(window_id = wid, patch_id = pid,
             area_km2 = top_patches$area_km2[i],
             touches_edge = touches)
})

trunc_df <- do.call(rbind, trunc_results)

cat("\n=== Truncation check: top", n_check, "largest patches ===\n")
print(trunc_df)
cat("\nPatches touching window edge:", sum(trunc_df$touches_edge),
    "of", nrow(trunc_df), "\n")


# order the combined patch table so that, within each centroid pair,
# the largest patch area comes first; ties will keep the first encountered row
all_patches_ord <- all_patches[
  order(
    all_patches$centroid_x,                # primary grouping key (X)
    all_patches$centroid_y,                # secondary grouping key (Y)
    -all_patches$area_km2                  # sort area descending within centroid
  ),
]

# this avoids carrying unnecessary columns during duplicate detection
centroids_only <- all_patches_ord[
  , c("centroid_x", "centroid_y")          # columns used for spatial identity
]


# create a logical vector marking duplicated centroid pairs
# after ordering, the FIRST row per (centroid_x, centroid_y) is the one we want to keep
dup_centroid <- duplicated(
  all_patches_ord[, c("centroid_x", "centroid_y")]
)


# fromLast = TRUE ensures *all* members of a duplicated group are flagged
print(paste0("Some patches may have been exactly counted more than once. How many times did this happen? ",sum(dup_centroid)))


# remove duplicated centroid entries, keeping only the largest-area patch per centroid
all_UniquePatches <- all_patches_ord[!dup_centroid,  ]
head(all_UniquePatches)

duplicated_patches <- all_patches_ord[
  dup_centroid,                            # select duplicated centroid rows
]

nrow(duplicated_patches)                   # should equal sum(dup_centroid)

dup_groups_logical <- duplicated(
  all_patches_ord[, c("centroid_x", "centroid_y")]
) | duplicated(
  all_patches_ord[, c("centroid_x", "centroid_y")],
  fromLast = TRUE
)

duplicated_groups <- all_patches_ord[
  dup_groups_logical,                      # all rows involved in duplication
]

#order duplicated groups clearly for visual inspection
###View this object if you want to inspect the duplicated patches
duplicated_groups <- duplicated_groups[
  order(
    duplicated_groups$centroid_x,          # group by centroid X
    duplicated_groups$centroid_y,          # then centroid Y
    -duplicated_groups$area_km2             # largest area first within group
  ),
]

#percentage of patches measure > once duplicated...
per_duped<-nrow(duplicated_groups)/nrow(all_patches_ord)

print(paste0("The total proprtion of patches exactly sampled more than once is:", per_duped))


#############################################################################
# BOUNDARY-BASED DEDUPLICATION (using TIF extents)
#############################################################################
# With ~44% overlap between adjacent windows, each point in space is
# covered by 1–4 windows.  A patch straddling window boundaries gets
# measured (with different cropping) in each overlapping window, producing
# near-duplicate rows with slightly different centroids, areas, and
# perimeters.
#
# Strategy: for each window, define a "core zone" — the window extent
# inset by a margin on each side.  Keeping only patches whose centroid
# falls in their source window's core zone removes boundary duplicates.
#
# IMPORTANT: The window-ID-to-grid-position mapping from the original run
# (windows 1–2430) was lost due to a code version change.  Rather than
# reconstructing the grid from expand.grid() (which gives wrong IDs for
# the original run), we read the ACTUAL geographic extent of each window
# directly from its output .tif file.  This bypasses the ID mapping issue
# entirely and works for both the original and resume runs.
#
# Edge windows (those whose extent aligns with the raster boundary) get
# no margin on that side, since there is no adjacent window to overlap.
#############################################################################

library(terra)

# ---- Load raster metadata (for edge detection and margin computation) -----
r <- rast("data/def_00_20.tif")
rast_xmin <- xmin(r)
rast_xmax <- xmax(r)
rast_ymin <- ymin(r)
rast_ymax <- ymax(r)
res_x     <- res(r)[1]    # cell width in degrees (longitude)
res_y     <- res(r)[2]    # cell height in degrees (latitude)
nr        <- nrow(r)
nc        <- ncol(r)
rm(r); gc(verbose = FALSE)

# ---- Compute margin to inset from each window edge -----------------------
# The grid uses 4500-cell spacing with 4000-cell extraction, giving
# ~44% overlap.  The margin is half the overlap in geographic units.
win_grid    <- 4500
win_extract <- 4000
n_windows_row <- ceiling((nr - win_grid) / (win_grid / 2)) + 1   # 55
n_windows_col <- ceiling((nc - win_grid) / (win_grid / 2)) + 1   # 66
step_row <- (nr - win_grid) / (n_windows_row - 1)                # ~2237 px
step_col <- (nc - win_grid) / (n_windows_col - 1)                # ~2244 px
margin_row <- (win_extract - step_row) / 2                        # ~881 px
margin_col <- (win_extract - step_col) / 2                        # ~878 px
margin_x <- margin_col * res_x   # margin in degrees longitude
margin_y <- margin_row * res_y   # margin in degrees latitude

cat("Margin to trim:", round(margin_x, 4), "deg lon,",
    round(margin_y, 4), "deg lat\n")

# ---- Read actual window extents from output TIF files ---------------------
# Each window's .patches_filtered.tif carries its true geographic extent,
# bypassing any window-ID-to-grid-position mapping issues.
tif_files <- list.files(
  "output/deforestTileOutStats/",
  pattern = "^w[0-9]+\\.patches_filtered\\.tif$",
  full.names = TRUE
)

tif_extents <- lapply(tif_files, function(f) {
  wid <- as.integer(gsub("^w|\\.patches_filtered\\.tif$", "", basename(f)))
  r <- rast(f)
  e <- ext(r)
  rm(r)
  data.frame(window_id = wid,
             tif_xmin = e$xmin, tif_xmax = e$xmax,
             tif_ymin = e$ymin, tif_ymax = e$ymax)
})
tif_ext_df <- do.call(rbind, tif_extents)
gc(verbose = FALSE)

cat("TIF extents extracted:", nrow(tif_ext_df), "windows\n")

# ---- Compute core zone for each window from its TIF extent ----------------
# Edge tolerance: if a window edge is within 2 pixels of the raster edge,
# treat it as an edge window (no margin trimmed on that side).
edge_tol_x <- res_x * 2
edge_tol_y <- res_y * 2

tif_ext_df$core_xmin <- ifelse(
  abs(tif_ext_df$tif_xmin - rast_xmin) < edge_tol_x,
  tif_ext_df$tif_xmin, tif_ext_df$tif_xmin + margin_x
)
tif_ext_df$core_xmax <- ifelse(
  abs(tif_ext_df$tif_xmax - rast_xmax) < edge_tol_x,
  tif_ext_df$tif_xmax, tif_ext_df$tif_xmax - margin_x
)
tif_ext_df$core_ymin <- ifelse(
  abs(tif_ext_df$tif_ymin - rast_ymin) < edge_tol_y,
  tif_ext_df$tif_ymin, tif_ext_df$tif_ymin + margin_y
)
tif_ext_df$core_ymax <- ifelse(
  abs(tif_ext_df$tif_ymax - rast_ymax) < edge_tol_y,
  tif_ext_df$tif_ymax, tif_ext_df$tif_ymax - margin_y
)

# ---- Filter patches: keep only those whose centroid is in core zone -------
all_UniquePatches_core <- merge(
  all_UniquePatches,
  tif_ext_df[, c("window_id",
                  "core_xmin", "core_xmax",
                  "core_ymin", "core_ymax")],
  by = "window_id",
  all.x = TRUE
)

in_core <- with(all_UniquePatches_core,
  !is.na(core_xmin) &
  centroid_x >= core_xmin & centroid_x < core_xmax &
  centroid_y >= core_ymin & centroid_y < core_ymax
)

n_core_before <- nrow(all_UniquePatches_core)
all_UniquePatches_core <- all_UniquePatches_core[in_core, ]
n_core_after <- nrow(all_UniquePatches_core)

message(
  "Boundary-based dedup: kept ", n_core_after,
  " of ", n_core_before, " patches (removed ",
  n_core_before - n_core_after, " boundary patches)."
)

# Drop the core-zone columns (no longer needed)
all_UniquePatches_core <- all_UniquePatches_core[
  , !names(all_UniquePatches_core) %in%
      c("core_xmin", "core_xmax", "core_ymin", "core_ymax")
]


#############################################################################
# PROXIMITY-BASED DEDUPLICATION (original method — kept for comparison)
#############################################################################
# The section below is the original proximity-based approach.  It defines
# a per-patch radius (60% of equivalent circle radius), builds a graph of
# "too close" patches, and keeps the largest in each cluster.
#
# With the boundary-based method above, this section is largely redundant
# for deduplication.  It may still be useful for identifying genuinely
# adjacent/overlapping patches (ecological connectivity).  The models at
# the bottom of this script are fit to BOTH:
#   - all_UniquePatches_core  ("boundary-filtered")
#   - all_UniquePatches       ("centroid-only", no boundary or proximity filter)
#############################################################################

# head(all_UniquePatches)
# 
# #these are 60% of the radii if the area was a circle
# withinRadius<-sqrt(all_UniquePatches$area_km2/pi)*.6
# head(withinRadius)
# 
# # load required packages for spatial distance calculations
# library(sf)          # modern spatial framework, efficient distance ops
# library(units)       # ensures distances have explicit units
# 
# #convert patch table to an sf POINT object using centroid lon/lat
# # coords = c("centroid_x","centroid_y") tells sf which columns are spatial
# # crs = 4326 defines WGS84 lon/lat, required for geodesic distances
# patch_sf_ll  <- st_as_sf(
#   all_UniquePatches,                     # data.frame of patch metrics
#   coords = c("centroid_x", "centroid_y"),
#   crs = 4326,
#   remove = FALSE                       # keep original columns
# )
# 
# # project centroids to equal-area CRS in meters (same CRS you used for area)
# patch_sf_m <- st_transform(
#   patch_sf_ll,
#   crs = 5880                            # South America Albers Equal Area (meters)
# )
# 
# # compute pairwise geodesic distances (meters) between all patch centroids
# # by default, sf uses great-circle distances for lon/lat CRS
# # result is a units matrix (n x n)
# dist_mat <- st_distance(
#   patch_sf_m,                            # sf POINT object
#   patch_sf_m                             # distances among all points
# )
# 
# # remove self-distances
# diag(dist_mat) <- NA
# 
# # drop units only after computation
# dist_mat_m <- units::drop_units(dist_mat)
# 
# # ensure withinRadius is in meters
# # if your withinRadius is currently in km, convert once here
# withinRadius_m <- withinRadius * 1000
# #withinRadius_m <- withinRadius         # assume already meters unless you say otherwise
# 
# # test, for each patch i, whether ANY other patch j is within its radius
# # sweep compares each row of dist_mat_m to that patch's radius
# # "<=" means inclusive threshold
# has_neighbor_within_radius <- apply(
#   dist_mat_m <= withinRadius_m,        # logical matrix of distance tests
#   1,                                   # operate row-wise (per patch)
#   any,                                 # TRUE if at least one neighbor qualifies
#   na.rm = TRUE
# )
# 
# # add the logical result back to your patch table
# all_UniquePatches$has_neighbor <- has_neighbor_within_radius
# 
# # OPTIONAL: count how many neighbors each patch has within its radius
# # useful for connectivity or isolation metrics
# all_UniquePatches$n_neighbors <- rowSums(
#   dist_mat_m <= withinRadius_m,
#   na.rm = TRUE
# )
# 
# # OPTIONAL: extract indices of neighboring patches for network-style analyses
# # returns a list where each element contains indices of nearby patches
# neighbor_index_list <- apply(
#   dist_mat_m <= withinRadius_m,
#   1,
#   function(x) which(x)
# )
# 
# # attach neighbor lists if you want explicit IDs
# # converts indices to patch_id values
# all_UniquePatches$neighbor_patch_ids <- lapply(
#   neighbor_index_list,
#   function(idx) all_patches_ord$patch_id[idx]
# )
# 
# # load igraph for clean connected-component handling
# library(igraph)
# 
# # number of patches before filtering
# n_before <- nrow(all_UniquePatches)
# 
# # build logical matrix indicating which patches are "too close"
# # reuse existing distance comparison (no recomputation)
# close_mat <- dist_mat_m <= withinRadius_m   # TRUE = too close
# diag(close_mat) <- FALSE                    # no self-links
# 
# # convert logical adjacency matrix to an edge list
# # which(..., arr.ind = TRUE) gives row/col indices of TRUE entries
# edge_idx <- which(close_mat, arr.ind = TRUE)
# 
# # build an undirected graph where vertices are patch ROW INDICES
# # each edge means "these two patches are too close"
# g <- graph_from_edgelist(
#   as.matrix(edge_idx),
#   directed = FALSE
# )
# 
# # compute connected components (clusters of mutually close patches)
# comp <- components(g)
# 
# # Start by keeping ALL patches. Only patches involved in proximity
# # conflicts (i.e., present as vertices in the graph) need resolution.
# # Previously this was rep(FALSE, ...) which silently dropped every
# # isolated patch — the exact ones we're most confident about.
# keep_patch <- rep(TRUE, n_before)
# 
# # Identify which row indices appear as vertices in the graph.
# # graph_from_edgelist() only creates vertices for patches that have
# # at least one edge, so patches with no proximity conflicts are absent
# # from the graph and should remain TRUE (kept) by default.
# vertices_in_graph <- as.integer(V(g)$name)
# keep_patch[vertices_in_graph] <- FALSE   # reset contested patches; resolve below
# 
# # iterate over each connected component to resolve contested patches
# for (k in seq_len(comp$no)) {
#   
#   # row indices belonging to this cluster
#   members <- which(comp$membership == k)
#   
#   # if cluster has only one patch, keep it
#   if (length(members) == 1) {
#     keep_patch[members] <- TRUE
#     next
#   }
#   
#   # select the patch with the maximum area within the cluster
#   # ties are resolved by keeping the first occurrence
#   best <- members[
#     which.max(all_UniquePatches$area_km2[members])
#   ]
#   
#   # keep only the largest-area patch
#   keep_patch[best] <- TRUE
# }
# 
# # subset the patch table to retained patches
# all_UniquePatches_filtered <- all_UniquePatches[keep_patch, ]
# 
# # count how many patches were removed
# n_after <- nrow(all_UniquePatches_filtered)
# n_removed <- n_before - n_after
# 
# # print a clear summary statement
# message(
#   "Removed ", n_removed, " patches due to proximity filtering (",
#   n_after, " patches retained from ", n_before, ")."
# )
# 

#############Visualizing aggregate patch stats patch stats#################

head(all_UniquePatches_core)

#just plot the small areas
hist(
  all_UniquePatches_core$area_km2,
  breaks = 500,                                  
  xlim = c(0, 200),
  xlab = "Patch area [km²]",
  ylab = "Frequency",
  main = "Patch area distribution (boundary-filtered)"
)

#just plot the large areas
hist(
  all_UniquePatches_core$area_km2[all_UniquePatches_core$area_km2>50],
  breaks = 50,                                  
  #xlim = c(0, 200),
  xlab = "Patch area [km²]",
  ylab = "Frequency",
  main = "Patch area distribution (boundary-filtered, >50 km²)"
)

#compactness
hist(
  all_UniquePatches_core$compactness ,
  breaks = 50,                                  
  #xlim = c(0, 200),
  xlab = "Patch compactness",
  ylab = "Frequency",
  main = "Patch compactness distribution (boundary-filtered)"
)

sum(all_UniquePatches_core$compactness>=0.0000001)/nrow(all_UniquePatches_core)
#I think the measure of compactness (1 is a perfect circle and 0.7 could be considered highly compact) is strongly affected by the high resolution of the original data rasters. This means that a lot of the shapes (I inspected them manually) that I would think are compact have very high perimeters (think about the classic "measuring the coastline of england" problem)
#as such, i think it may not be useful to consider compactness. it will reduce the realism of our projections a bit, but I cannot think of a justifiable compactness filtering method with the data we have

#############################################################################
# FIT AND EXPORT DISTRIBUTION MODELS
#############################################################################
# Refactored into a single function to avoid duplication across filtering
# regimes. Fits power-law and lognormal, runs a formal Vuong test, and
# exports frequency tables + model parameters.
#
# Three filtering regimes (matching stage 4 expectations):
#   - boundary_filtered: core-zone dedup (principled, main analysis)
#   - underfiltered: centroid-only dedup (retains near-dupes from overlap)
#   - overfiltered: boundary-filtered + additional 10 km² area threshold
#############################################################################

library(poweRlaw)

fit_and_export_models <- function(patch_area_vec, suffix, out_dir = "output") {
  # Fit power-law and lognormal distributions, compare via Vuong test, export.
  #
  # Args:
  #   patch_area_vec: numeric vector of patch areas (km²)
  #   suffix: string for output file naming (e.g. "boundary_filtered")
  #   out_dir: output directory (default "output")
  #
  # Returns (invisibly): list with fitted model objects and summary tables

  patch_area <- patch_area_vec[is.finite(patch_area_vec) & patch_area_vec > 0]

  # --- Logarithmic binning ---
  log_breaks <- 10^seq(
    from = floor(log10(min(patch_area))),
    to   = ceiling(log10(max(patch_area))),
    length.out = 25
  )

  patch_hist <- hist(patch_area, breaks = log_breaks, plot = FALSE)

  patch_bins <- data.frame(
    bin_lower = patch_hist$breaks[-length(patch_hist$breaks)],
    bin_upper = patch_hist$breaks[-1],
    count     = patch_hist$counts
  )
  patch_bins$bin_mid       <- sqrt(patch_bins$bin_lower * patch_bins$bin_upper)
  patch_bins$log10_bin_mid <- log10(patch_bins$bin_mid)

  # --- Log-log plot of binned distribution ---
  plot(
    x = patch_bins$bin_lower, y = patch_bins$count,
    log = "xy", type = "b",
    xlab = "Patch area (km², log10 bins)", ylab = "Patch count",
    main = paste("Patch area distribution:", suffix)
  )

  # --- Fit power-law (MLE) ---
  pl_model <- conpl$new(patch_area)
  pl_est   <- estimate_xmin(pl_model)
  pl_model$setXmin(pl_est)

  alpha_hat <- pl_model$getPars()
  xmin_hat  <- pl_model$getXmin()

  cat("\n--- ", suffix, " ---\n")
  cat("Power-law: alpha =", alpha_hat, ", xmin =", xmin_hat, "\n")

  # --- Fit lognormal (same xmin for fair comparison) ---
  ln_model <- conlnorm$new(patch_area)
  ln_model$setXmin(xmin_hat)
  ln_model$setPars(estimate_pars(ln_model))

  ln_pars <- ln_model$getPars()
  cat("Lognormal: meanlog =", ln_pars[1], ", sdlog =", ln_pars[2], "\n")

  # --- Formal model comparison: Vuong likelihood ratio test ---
  vuong <- compare_distributions(pl_model, ln_model)

  favored <- if (vuong$test_statistic > 0) "power_law" else "lognormal"
  cat("Vuong test: statistic =", round(vuong$test_statistic, 3),
      ", p_two_sided =", round(vuong$p_two_sided, 4),
      ", favors:", favored, "\n")

  # --- Visualize CDF with both fits ---
  plot(pl_model, main = paste("CDF with fits:", suffix))
  lines(pl_model, col = "blue", lwd = 2)
  lines(ln_model, col = "red", lwd = 2)
  legend("bottomleft", legend = c("Power-law", "Lognormal"),
         col = c("blue", "red"), lwd = 2, bty = "n")

  # --- Export frequency table ---
  write.csv(
    patch_bins,
    file = file.path(out_dir, paste0("patch_area_frequency_logbins_", suffix, ".csv")),
    row.names = FALSE
  )

  # --- Export lognormal parameters ---
  lognormal_summary <- data.frame(
    model          = "lognormal",
    xmin_km2       = ln_model$getXmin(),
    meanlog        = ln_pars[1],
    sdlog          = ln_pars[2],
    n_patches      = length(patch_area),
    n_above_xmin   = sum(patch_area >= ln_model$getXmin()),
    total_area_km2 = sum(patch_area)
  )
  write.csv(
    lognormal_summary,
    file = file.path(out_dir, paste0("lognormal_model_parameters_", suffix, ".csv")),
    row.names = FALSE
  )

  # --- Export power-law parameters ---
  pl_summary <- data.frame(
    model          = "power_law",
    xmin_km2       = pl_model$getXmin(),
    alpha          = pl_model$getPars(),
    n_patches      = length(patch_area),
    n_above_xmin   = sum(patch_area >= pl_model$getXmin()),
    total_area_km2 = sum(patch_area)
  )
  write.csv(
    pl_summary,
    file = file.path(out_dir, paste0("powerlaw_model_parameters_", suffix, ".csv")),
    row.names = FALSE
  )

  # --- Vuong summary row (for combined export) ---
  vuong_summary <- data.frame(
    regime         = suffix,
    comparison     = "powerlaw_vs_lognormal",
    test_statistic = vuong$test_statistic,
    p_two_sided    = vuong$p_two_sided,
    favors         = favored
  )

  message("Exported ", suffix, " model outputs to ", out_dir)

  invisible(list(
    patch_bins        = patch_bins,
    pl_model          = pl_model,
    ln_model          = ln_model,
    vuong             = vuong_summary,
    lognormal_summary = lognormal_summary,
    pl_summary        = pl_summary
  ))
}


#############################################################################
# CLEAN STALE OUTPUT FILES AND PREPARE REGIMES
#############################################################################

# Remove all previous model parameter / frequency CSVs to avoid stale files
# from old code versions (e.g., "centroid_only" from before the rename)
stale_patterns <- c(
  "^(lognormal|powerlaw)_model_parameters_.*\\.csv$",
  "^patch_area_frequency_logbins_.*\\.csv$",
  "^vuong_test_results\\.csv$"
)
stale_files <- unlist(lapply(stale_patterns, function(p) {
  list.files("output", pattern = p, full.names = TRUE)
}))
if (length(stale_files) > 0) {
  file.remove(stale_files)
  message("Removed ", length(stale_files), " stale output files before re-export.")
}

# Create the overfiltered variant: boundary-filtered with additional 10 km²
# minimum area threshold. This provides a conservative sensitivity check by
# removing the smallest patches most prone to measurement uncertainty.
area_threshold_km2 <- 10
all_UniquePatches_overfiltered <- all_UniquePatches_core[
  all_UniquePatches_core$area_km2 >= area_threshold_km2,
]
message(
  "Overfiltered: ", nrow(all_UniquePatches_overfiltered),
  " patches (>= ", area_threshold_km2, " km²) from ",
  nrow(all_UniquePatches_core), " boundary-filtered patches."
)


#############################################################################
# FIT AND EXPORT ALL THREE REGIMES
#############################################################################

results_boundary <- fit_and_export_models(
  all_UniquePatches_core$area_km2, "boundary_filtered"
)

results_under <- fit_and_export_models(
  all_UniquePatches$area_km2, "underfiltered"
)

results_over <- fit_and_export_models(
  all_UniquePatches_overfiltered$area_km2, "overfiltered"
)

# Combined Vuong test comparison across regimes
vuong_all <- rbind(results_boundary$vuong, results_under$vuong, results_over$vuong)
write.csv(vuong_all, "output/vuong_test_results.csv", row.names = FALSE)
cat("\n=== Vuong test summary ===\n")
print(vuong_all)

message("\nAll models fitted and exported.")

