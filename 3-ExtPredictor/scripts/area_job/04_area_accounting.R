#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# 04_area_accounting.R
# Compute total deforested area (km²) from the binary raster.
#
# Checkpoints intermediate results to RDS files so expensive steps can be
# skipped on re-run. Safe to interrupt and restart.
#
# Usage (from anywhere):
#   Rscript 04_area_accounting.R
#
# Requirements:
#   - R (>= 4.1)
#   - terra package
#   - def_00_20.tif in 3-ExtPredictor/data/
#
# Output:
#   - output/area_accounting.csv       (final one-row summary)
#   - output/checkpoint_cell_counts.rds (cell count checkpoint)
#   - output/checkpoint_defor_area.rds  (deforested area checkpoint)
#   - output/checkpoint_extent_area.rds (study-area extent checkpoint)
# ---------------------------------------------------------------------------

library(terra)

# --- resolve script location and working directory ---
script_dir <- NULL
if (!interactive()) {
  args <- commandArgs(trailingOnly = FALSE)
  script_path <- normalizePath(
    sub("--file=", "", args[grep("--file=", args)]),
    mustWork = TRUE
  )
  script_dir <- dirname(script_path)
  setwd(script_dir)
} else {
  script_dir <- getwd()
}
cat("Working directory:", getwd(), "\n")

dir.create("output", showWarnings = FALSE)

# --- load raster ---
rast_path <- normalizePath(file.path(script_dir, "..", "..", "data", "def_00_20.tif"))
stopifnot(file.exists(rast_path))

r <- rast(rast_path)
cat("Raster dimensions:", nrow(r), "rows x", ncol(r), "cols\n")
cat("CRS:", crs(r, describe = TRUE)$name, "\n")

# -------------------------------------------------------
# STEP 1: Count deforested and data cells
# -------------------------------------------------------
ckpt1 <- "output/checkpoint_cell_counts.rds"

if (file.exists(ckpt1)) {
  cat("\n[checkpoint] Loading cell counts from", ckpt1, "\n")
  counts <- readRDS(ckpt1)
} else {
  cat("\nStep 1: Counting cells...\n")
  n_def_cells  <- global(r, fun = "sum", na.rm = TRUE)[[1]]
  n_total      <- ncell(r)
  n_na_cells   <- global(is.na(r), fun = "sum")[[1]]
  n_data_cells <- n_total - n_na_cells

  counts <- list(
    n_def_cells  = n_def_cells,
    n_data_cells = n_data_cells,
    n_total      = n_total,
    n_na_cells   = n_na_cells
  )
  saveRDS(counts, ckpt1)
  cat("[checkpoint] Saved cell counts to", ckpt1, "\n")
}

cat("Deforested cells:", format(counts$n_def_cells, big.mark = ","), "\n")
cat("Total non-NA cells:", format(counts$n_data_cells, big.mark = ","), "\n")

# -------------------------------------------------------
# STEP 2: Compute total deforested area (latitude-weighted)
# -------------------------------------------------------
ckpt2 <- "output/checkpoint_defor_area.rds"

if (file.exists(ckpt2)) {
  cat("\n[checkpoint] Loading deforested area from", ckpt2, "\n")
  total_defor_km2 <- readRDS(ckpt2)
} else {
  cat("\nStep 2: Computing cell areas and deforested area sum...\n")
  cat("(this is the slow step)\n")
  t0 <- proc.time()

  cell_area_m2    <- cellSize(r, unit = "m")
  defor_area_m2   <- r * cell_area_m2
  total_defor_m2  <- global(defor_area_m2, fun = "sum", na.rm = TRUE)[[1]]
  total_defor_km2 <- total_defor_m2 / 1e6

  elapsed <- (proc.time() - t0)["elapsed"]
  cat(sprintf("Done in %.1f seconds\n", elapsed))

  saveRDS(total_defor_km2, ckpt2)
  cat("[checkpoint] Saved deforested area to", ckpt2, "\n")
}

cat("Total deforested area:", format(round(total_defor_km2, 1), big.mark = ","), "km2\n")

# -------------------------------------------------------
# STEP 3: Compute total study-area extent (all non-NA cells)
# -------------------------------------------------------
ckpt3 <- "output/checkpoint_extent_area.rds"

if (file.exists(ckpt3)) {
  cat("\n[checkpoint] Loading study-area extent from", ckpt3, "\n")
  total_extent_km2 <- readRDS(ckpt3)
} else {
  cat("\nStep 3: Computing total study-area extent...\n")
  t0 <- proc.time()

  cell_area_m2     <- cellSize(r, unit = "m")
  total_extent_m2  <- global(cell_area_m2 * !is.na(r), fun = "sum", na.rm = TRUE)[[1]]
  total_extent_km2 <- total_extent_m2 / 1e6

  elapsed <- (proc.time() - t0)["elapsed"]
  cat(sprintf("Done in %.1f seconds\n", elapsed))

  saveRDS(total_extent_km2, ckpt3)
  cat("[checkpoint] Saved study-area extent to", ckpt3, "\n")
}

cat("Total study-area extent:", format(round(total_extent_km2, 1), big.mark = ","), "km2\n")

# -------------------------------------------------------
# REPORT & SAVE
# -------------------------------------------------------
cat("\n=== AREA ACCOUNTING ===\n")
cat("Total deforested area:  ", format(round(total_defor_km2, 1), big.mark = ","), "km2\n")
cat("Total study-area extent:", format(round(total_extent_km2, 1), big.mark = ","), "km2\n")
cat(sprintf("Deforested fraction:     %.2f%%\n", 100 * total_defor_km2 / total_extent_km2))

out <- data.frame(
  n_deforested_cells = counts$n_def_cells,
  n_data_cells       = counts$n_data_cells,
  total_defor_km2    = total_defor_km2,
  total_extent_km2   = total_extent_km2,
  defor_fraction     = total_defor_km2 / total_extent_km2
)

write.csv(out, "output/area_accounting.csv", row.names = FALSE)
cat("\nSaved to output/area_accounting.csv\n")
cat("Done.\n")
