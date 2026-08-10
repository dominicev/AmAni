# ============================================================================
# Standalone resume runner for deterministic patch analysis
# Extracted from AmazonHabitatLoss_prep.v2.R (lines 269–463)
#
# This script runs the resume block to process remaining windows.
# It automatically detects the working directory and is cross-platform
# (Linux and Windows).
#
# Usage:
#   Rscript run_resume_deterministic.R
#   OR (for background execution on Linux):
#   nohup Rscript run_resume_deterministic.R > resume_run.log 2>&1 &
#
# ============================================================================

# ---- Set working directory to the script's own folder --------------------
args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", args[grep("^--file=", args)])
# Navigate to project root (one level above scripts/)
if (length(script_path) > 0 && nzchar(script_path)) {
  setwd(file.path(dirname(normalizePath(script_path)), ".."))
}

message("Working directory: ", getwd())
message("Start time: ", Sys.time())

# ---- Load packages -------------------------------------------------------
library(terra)
library(parallel)
library(doParallel)
library(foreach)

# ============================================================================
# RESUME RUN — Complete remaining deterministic windows
# ============================================================================

# ---- Detect hardware resources at runtime --------------------------------
# Use physical (not logical/hyper-threaded) cores for heavy IO + compute jobs.
# Reserve ~30 % of cores for OS stability.
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
# These parameters MUST match the grid computed in the original run
# so that window IDs correspond to the same spatial positions.
win_nrow_grid <- 4500
win_ncol_grid <- 4500

temp_rast <- rast("data/def_00_20.tif")
nr <- nrow(temp_rast)
nc <- ncol(temp_rast)
rm(temp_rast); gc(verbose = FALSE)

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
remaining_ids  <- setdiff(grid_positions$window_id, done_ids)
grid_remaining <- grid_positions[grid_positions$window_id %in% remaining_ids, ]

cat("Windows already completed:", length(done_ids), "\n")
cat("Windows remaining:", nrow(grid_remaining), "\n")
cat("Remaining IDs (sample):",
    head(grid_remaining$window_id, 10), "...\n\n")

if (nrow(grid_remaining) == 0) {
  message("All windows already complete. Exiting.")
  quit(save = "no")
}

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
  .packages = "terra",
  .errorhandling = "pass"
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
  gc(verbose = FALSE)
}

# ---- Tear down cluster ---------------------------------------------------
parallel::stopCluster(cl)
foreach::registerDoSEQ()

message("\n=== Resume run complete ===")
message("End time: ", Sys.time())
message("Check output/deforestTileOutStats/ for updated outputs.")
