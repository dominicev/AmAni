{#setwd — navigate to project root (one level above scripts/)
args <- commandArgs(trailingOnly = FALSE)
script_arg <- args[grep("^--file=", args)]

if (length(script_arg) > 0) {
  project_root <- file.path(dirname(normalizePath(sub("^--file=", "", script_arg))), "..")
} else if (requireNamespace("rstudioapi", quietly = TRUE) &&
           rstudioapi::isAvailable()) {
  project_root <- file.path(dirname(rstudioapi::getActiveDocumentContext()$path), "..")
} else {
  stop("Cannot determine project root. Run from RStudio or via Rscript.")
}
setwd(project_root)

#load dependencies
#library(data.table)
library(sf)
#library(utils)
#library(lwgeom)
library(ggplot2)
library(terra)
#library(readr)
#library(tidyr)
#library(stringr)
#library(parzer)
library(dplyr)

library(parallel)     # base R parallel backend
library(doParallel)   # foreach-compatible parallel backend
library(foreach)      # parallel looping

getwd()


#visualizing deforested regions of Amazon based on future projection data from: https://www.raisg.org/en/maps/
dataDir<-"data/raisg-shapefiles/deforestacion_futura_2025/"
shp_path <- "deforestacion_futura_2025.shp"

}
deforest_sf <- st_read(paste0(dataDir, shp_path))

st_geometry_type(deforest_sf)
st_crs(deforest_sf)
plot(st_geometry(deforest_sf))


names(deforest_sf)
glimpse(deforest_sf)

###Variable explanations
# supkm2_s60	  Cell area (km²)
# BosRema_	    Remaining forest area
# pro_m50	      Projected deforested area (baseline, ~50%)
# sup_prom50	  Mean deforested area (km²)
# esc_tend	    Deforestation under trend scenario
# esc_pesc	    Deforestation under pessimistic scenario
# esc_opti	    Deforestation under optimistic scenario

#plot only those deforested under most realistic projection
hist(deforest_sf$esc_tend)

deforest_tend <- deforest_sf |>  filter(esc_tend > quantile(deforest_sf$esc_tend, 0.75))#plot only cells in the 25% most likely to be deforested
plot(st_geometry(deforest_tend))

#plot only those deforested under the most optimistic projection
deforest_opt <- deforest_sf |>  filter(esc_opti > quantile(deforest_sf$esc_opti, 0.75))

plot(st_geometry(deforest_opt))


#density plot by most likely deforested area per cell
deforest_tend <- deforest_sf |>  filter(esc_tend > 0)
ggplot(deforest_tend) +
  geom_sf(aes(fill = esc_tend), color = NA) +
  scale_fill_viridis_c(name = "Deforestation\n(km²)") +
  coord_sf(expand = FALSE)

#density plot by mean deforested area per cell
deforest_tend <- deforest_sf |>  filter(esc_tend > 0)
ggplot(deforest_tend) +
  geom_sf(aes(fill = sup_prom50), color = NA) +
  scale_fill_viridis_c(name = "Deforestation\n(km²)") +
  coord_sf(expand = FALSE)

####
#visualizing deforested regions of Amazon based on total deforestation up to 2020 from: https://www.raisg.org/en/maps/
dataDir<-"data/raisg-shapefiles/deforestacion_2001_2020_total/"
list.files(dataDir)

#deforestAll_2020 <- rast(paste0(dataDir,"deforestacion2001-2020.tif")) #uncomment this section to load, look at, and convert the raw data file.
#plot(deforestAll_2020)

##converting this to polygon data would result in a 70GB+ size object, which is not feasible for us...
#summary(deforestAll_2020)
#head(deforestAll_2020)

# Label contiguous deforestation patches
#class(deforestAll_2020)
#ncell(deforestAll_2020)

#writeRaster(
#  deforestAll_2020 >= 2000,
#  "data/def_00_20.tif",
#  overwrite = TRUE,
#  datatype = "INT1U",
#  gdal = c("COMPRESS=DEFLATE", "TILED=YES")
#)

#import it
def_bin <- rast("data/def_00_20.tif")

#unique(def_bin) #the raster is now binary
plot(def_bin) #it still has all the location data
head(def_bin) #other than 1, it has NA values

#patches_00_20 <- patches(def_bin, directions = 8) # i ran this for 3 days and it only completed at about 40%. I think the amount of time probably increases exponentially the more you run it so i expect that to get to 50% would take more than 3 additional days
#writeRaster(patches_00_20, "patches_2000_2020.tif", overwrite = TRUE)

#I am going to use patches() to
# inspect raster size to decide on a reasonable test window
nrow(def_bin)/1000  # number of rows in the raster
ncol(def_bin)/1000 # number of columns in the raster

source("scripts/utils/calcPatchStat.sub.R")

res<-extract_patch_metrics(
    raster_path = "data/def_00_20.tif",     # SpatRaster: binary raster (1 = deforested, NA = background)
    seed = 85683,                  # numeric or NULL: random seed for reproducibility
    min_area_m2 = 5e6,         # numeric: minimum patch area to retain (5e6 is 50 km^2)
    win_nrow = 4000,           # numeric: the height of the window
    win_ncol = 4000,           # numeric: the width of the window
    out_dir = "output/deforestTileOutStats2/",  # character: directory where outputs will be written
    returnQ = TRUE,            # logical: return results
    apply_morph = FALSE,       # logical: apply morphological operations
    morph_size = 5,            # numeric: size of the erosion kernel
    small_patch_threshold = 50, # numeric: remove small patches before morphology
    method = "random"          # character: "random" or "deterministic" sampling method
    # start_row = NULL,        # numeric: starting row (required if method = "deterministic")
    # start_col = NULL         # numeric: starting column (required if method = "deterministic")
) 
summary(res)
plot(res$raster)

#######do this 10k times
# detect total number of available cores on the machine
n_cores_total <- parallel::detectCores(logical = FALSE)
n_cores_use   <- max(1, floor(n_cores_total * 0.4))

# Detect total system RAM (cross-platform) and pre-compute per-worker limit
total_ram_gb <- tryCatch({
  if (file.exists("/proc/meminfo")) {
    meminfo <- readLines("/proc/meminfo")
    mem_line <- grep("^MemTotal:", meminfo, value = TRUE)
    as.numeric(gsub("[^0-9]", "", mem_line)) / 1024 / 1024
  } else {
    si <- system("systeminfo", intern = TRUE)
    mem_line <- grep("Total Physical Memory", si, value = TRUE)
    as.numeric(gsub("[^0-9]", "", mem_line)) / 1024
  }
}, error = function(e) {
  warning("Could not detect RAM; defaulting to 32 GB")
  32
})
n_cores_for_mem <- max(1, floor(n_cores_total * 0.7))
worker_memmax <- (total_ram_gb / n_cores_for_mem) * 0.7

# start a brand new PSOCK cluster
cl <- parallel::makeCluster(n_cores_use)

# register it for foreach
doParallel::registerDoParallel(cl)

# Calculate grid positions for deterministic sampling with adjusted step size
# to ensure complete coverage with no overhang
win_nrow <- 4500
win_ncol <- 4500

# Get raster dimensions
temp_rast <- rast("data/def_00_20.tif")
nr <- nrow(temp_rast)
nc <- ncol(temp_rast)
rm(temp_rast)
gc()

# Calculate number of windows needed for ~50% overlap
# For 50% overlap, we'd normally use step = window_size / 2
# Number of windows = ceiling((total_size - window_size) / step_size) + 1
n_windows_row <- ceiling((nr - win_nrow) / (win_nrow / 2)) + 1
n_windows_col <- ceiling((nc - win_ncol) / (win_ncol / 2)) + 1

# Calculate adjusted step size to tile perfectly from edge to edge
step_row_adjusted <- (nr - win_nrow) / (n_windows_row - 1)
step_col_adjusted <- (nc - win_ncol) / (n_windows_col - 1)

# Generate starting positions with perfect edge-to-edge coverage
start_rows <- round(seq(1, nr - win_nrow + 1, length.out = n_windows_row))
start_cols <- round(seq(1, nc - win_ncol + 1, length.out = n_windows_col))

# Create all combinations
grid_positions <- expand.grid(start_row = start_rows, start_col = start_cols)
grid_positions$window_id <- 1:nrow(grid_positions)

# Report coverage info
cat("\n=== Deterministic Grid Sampling Configuration ===\n")
cat("Raster dimensions:", nr, "rows x", nc, "cols\n")
cat("Window dimensions:", win_nrow, "rows x", win_ncol, "cols\n")
cat("Adjusted step size:", round(step_row_adjusted, 1), "rows x", 
    round(step_col_adjusted, 1), "cols\n")
cat("Overlap percentage:", round((1 - step_row_adjusted/win_nrow) * 100, 1), 
    "% rows x", round((1 - step_col_adjusted/win_ncol) * 100, 1), "% cols\n")
cat("Number of windows (rows):", length(start_rows), "\n")
cat("Number of windows (cols):", length(start_cols), "\n")
cat("Total windows:", nrow(grid_positions), "\n")
cat("Coverage check - Last window ends at row:", max(start_rows) + win_nrow - 1, 
    "(should be", nr, ")\n")
cat("Coverage check - Last window ends at col:", max(start_cols) + win_ncol - 1, 
    "(should be", nc, ")\n")
cat("===============================================\n\n")

# define paths ONCE in the master so they are captured by closure
func_file <- normalizePath("scripts/utils/calcPatchStat.sub.R")              # function source file
out_dir_i <- normalizePath("output/deforestTileOutStats/", mustWork=FALSE)  # output directory
raster_in <- normalizePath("data/def_00_20.tif")                    # input raster

# Set workers' working directory and export objects they need
abs_project_dir <- getwd()
parallel::clusterExport(cl, "abs_project_dir", envir = environment())
parallel::clusterEvalQ(cl, setwd(abs_project_dir))
parallel::clusterExport(cl, c("grid_positions", "worker_memmax"), envir = environment())

# run extract_patch_metrics() in parallel over grid positions
# each iteration is independent, making this hugely parallel
results <- foreach(
  i = 1:nrow(grid_positions),          # iterate over all grid positions
  .packages = "terra"                  # packages loaded on workers
  # , .errorhandling = "pass"          # continue even if some fail (uncomment to enable)
) %dopar% {
  
  # Apply pre-computed memory limit
  terraOptions(memmax = worker_memmax)
  
  # unset RNG state explicitly
  set.seed(NULL)
  
  # source the function definition on the worker
  source(func_file)                                # defines extract_patch_metrics()
  
  # Get position for this window
  pos <- grid_positions[i, ]
  window_id <- pos$window_id
  
  # skip completed windows
  done_flag_file <- file.path(out_dir_i, paste0("w", window_id, ".done"))
  if (file.exists(done_flag_file)) {
    return(NULL)
  }

  # run the patch extraction and metric calculation
  extract_patch_metrics(
    raster_in,
    output_id   = window_id,
    min_area_m2 = 5e6,
    win_nrow    = 4000,
    win_ncol    = 4000,
    out_dir     = out_dir_i,
    returnQ     = FALSE,
    method      = "deterministic",
    start_row   = pos$start_row,
    start_col   = pos$start_col
  )
  
  # cleanup
  rm(list = ls(), inherits = TRUE)
  gc(FALSE)

}
# stop cluster after completion
parallel::stopCluster(cl)
foreach::registerDoSEQ()


# ============================================================================
# RESUME RUN — Complete remaining deterministic windows
# ============================================================================
# The initial deterministic run (above) processed windows 1–2430 but the full
# grid contains 3630 windows (55 rows × 66 cols).  Windows 2431–3630 were
# never started, and windows 303 and 533 are missing their .done flags.
#
# This block:
#   1. Rebuilds the same deterministic grid used above (win 4500 layout,
#      actual extraction at 4000×4000).
#   2. Identifies every window that lacks a .done flag in the output dir.
#   3. Runs only those remaining windows in parallel.
#
# Machine-adaptive parallelisation:
#   Core count and RAM are detected at runtime so this block is portable
#   across machines without manual edits.
# ============================================================================

library(terra)
library(parallel)
library(doParallel)
library(foreach)

# ---- Detect hardware resources at runtime --------------------------------
# Use physical (not logical/hyper-threaded) cores for heavy IO + compute jobs.
# Reserve ~30 % of cores for OS / RStudio stability.
n_cores_physical <- parallel::detectCores(logical = FALSE)
n_cores_use      <- max(1, floor(n_cores_physical * 0.4))

# Detect total system RAM (in GB).
# Cross-platform: try Linux /proc/meminfo first, then Windows systeminfo.
total_ram_gb <- tryCatch({
  if (file.exists("/proc/meminfo")) {
    # Linux: parse MemTotal from /proc/meminfo (value is in kB)
    meminfo <- readLines("/proc/meminfo")
    mem_line <- grep("^MemTotal:", meminfo, value = TRUE)
    mem_kb <- as.numeric(gsub("[^0-9]", "", mem_line))
    mem_kb / 1024 / 1024
  } else {
    # Windows: parse output of `systeminfo`
    si <- system("systeminfo", intern = TRUE)
    mem_line <- grep("Total Physical Memory", si, value = TRUE)
    mem_mb <- as.numeric(gsub("[^0-9]", "", mem_line))
    mem_mb / 1024
  }
}, error = function(e) {
  warning("Could not detect RAM; defaulting to 32 GB")
  32
})

cat("\n=== Machine Resources ===\n")
cat("Physical cores detected:", n_cores_physical, "\n")
cat("Cores to use:", n_cores_use, "\n")
cat("Total RAM detected:", round(total_ram_gb, 1), "GB\n")
cat("=========================\n\n")

# ---- Rebuild the deterministic grid (must match original layout) ----------
# These parameters MUST match the grid computed in the original run above
# so that window IDs correspond to the same spatial positions.
win_nrow_grid <- 4500
win_ncol_grid <- 4500

temp_rast <- rast("data/def_00_20.tif")
nr <- nrow(temp_rast)
nc <- ncol(temp_rast)
rm(temp_rast); gc()

# Number of windows for ~50 % overlap (same formula as original block)
n_windows_row <- ceiling((nr - win_nrow_grid) / (win_nrow_grid / 2)) + 1
n_windows_col <- ceiling((nc - win_ncol_grid) / (win_ncol_grid / 2)) + 1

# Adjusted step size for perfect edge-to-edge tiling
start_rows <- round(seq(1, nr - win_nrow_grid + 1, length.out = n_windows_row))
start_cols <- round(seq(1, nc - win_ncol_grid + 1, length.out = n_windows_col))

grid_positions <- expand.grid(start_row = start_rows, start_col = start_cols)
grid_positions$window_id <- 1:nrow(grid_positions)

cat("Total grid windows:", nrow(grid_positions), "\n")

# ---- Identify incomplete windows -----------------------------------------
# A window is considered complete only if its .done flag exists.
out_dir_resume <- normalizePath("output/deforestTileOutStats/", mustWork = FALSE)

done_flags <- list.files(out_dir_resume, pattern = "^w[0-9]+\\.done$")
done_ids   <- as.integer(gsub("^w|\\.done$", "", done_flags))

# Windows that still need to be processed
remaining_ids <- setdiff(grid_positions$window_id, done_ids)
grid_remaining <- grid_positions[grid_positions$window_id %in% remaining_ids, ]

cat("Windows already completed:", length(done_ids), "\n")
cat("Windows remaining:", nrow(grid_remaining), "\n")
cat("Remaining IDs include:",
    head(grid_remaining$window_id, 10), "...\n\n")

# ---- Resolve all paths to ABSOLUTE in the master process -----------------
# PSOCK cluster workers do NOT inherit the master's working directory —
# they start in the system default (e.g. $HOME on Linux).  Passing
# absolute paths avoids "file not found" errors on workers.
abs_project_dir <- getwd()
abs_func_file   <- normalizePath("scripts/utils/calcPatchStat.sub.R")
abs_raster_path <- normalizePath("data/def_00_20.tif")
abs_out_dir     <- normalizePath("output/deforestTileOutStats/", mustWork = FALSE)

# Pre-compute per-worker terra memory limit in the master so workers
# don't each need to call systeminfo / read /proc/meminfo
n_cores_for_mem <- max(1, floor(
  parallel::detectCores(logical = FALSE) * 0.7
))
worker_memmax <- (total_ram_gb / n_cores_for_mem) * 0.7
cat("Per-worker terra memmax:", round(worker_memmax, 1), "GB\n\n")

# ---- Launch parallel cluster ---------------------------------------------
cl <- parallel::makeCluster(n_cores_use)
doParallel::registerDoParallel(cl)

# Set workers' working directory to the project folder so that any
# relative paths inside calcPatchStat.sub.R also resolve correctly
parallel::clusterExport(cl, "abs_project_dir", envir = environment())
parallel::clusterEvalQ(cl, setwd(abs_project_dir))

# Export objects that workers need
parallel::clusterExport(
  cl,
  c("grid_remaining", "abs_func_file", "abs_raster_path",
    "abs_out_dir", "worker_memmax"),
  envir = environment()
)

# ---- Parallel foreach over remaining windows ------------------------------
results_resume <- foreach(
  i = seq_len(nrow(grid_remaining)),
  .packages = "terra"
) %dopar% {

  # Apply pre-computed memory limit
  terraOptions(memmax = worker_memmax)

  # Unset any inherited RNG state
  set.seed(NULL)

  # Source the extraction function on this worker (absolute path)
  source(abs_func_file)

  # Look up the grid position for this iteration
  pos       <- grid_remaining[i, ]
  window_id <- pos$window_id

  # Skip-check using absolute output directory
  done_flag_file <- file.path(abs_out_dir, paste0("w", window_id, ".done"))
  if (file.exists(done_flag_file)) {
    return(NULL)
  }

  # Run the patch extraction (all paths absolute)
  extract_patch_metrics(
    raster_path = abs_raster_path,
    output_id   = window_id,
    min_area_m2 = 5e6,
    win_nrow    = 4000,
    win_ncol    = 4000,
    out_dir     = abs_out_dir,
    returnQ     = FALSE,
    method      = "deterministic",
    start_row   = pos$start_row,
    start_col   = pos$start_col
  )

  # Free memory on worker before next iteration
  rm(list = ls(), inherits = TRUE)
  gc(FALSE)
}

# ---- Tear down cluster ----------------------------------------------------
parallel::stopCluster(cl)
foreach::registerDoSEQ()

cat("\n=== Resume run complete ===\n")
cat("Check output/deforestTileOutStats/ for updated outputs.\n")

