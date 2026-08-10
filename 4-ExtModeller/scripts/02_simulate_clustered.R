# =====================================================================
# 02_simulate_clustered.R — Spatially clustered deforestation model
# =====================================================================
# Same overall pipeline as 01_simulate_neutral.R, but deforestation
# patches are placed preferentially in grid cells that are closer to
# observed deforestation (from the RAISG raster centroids).
#
# The key difference: grid_weights — a named numeric vector that biases
# which grid cells get deforested. Grids near known deforestation
# hotspots have higher weight.
# =====================================================================

{
  # set working directory to project root (4-ExtModeller/)
  library(knitr)
  library(rstudioapi)
  if (exists(".amani_wd")) {
    setwd(.amani_wd)
  } else {
    project_root <- dirname(dirname(rstudioapi::getActiveDocumentContext()$path))
    knitr::opts_knit$set(root.dir = project_root)
    setwd(project_root)
  }

  library(sf)
  library(terra)
  library(ggplot2)
  library(dplyr)
  library(parallel)
  library(doParallel)
  library(foreach)

  getwd()
}

source("scripts/utils/simulation_functions.R")

# Output directory (overridden by manager when running organism-specific pipelines)
if (exists(".amani_output_dir")) {
  outdir <- .amani_output_dir
} else {
  outdir <- "output"
}
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)


############################
# IMPORT DATA (same as neutral script)
############################

extModelDir <- paste0(dirname(getwd()), "/3-ExtPredictor/output")

# model parameter CSVs
model_files <- list.files(
  path = extModelDir, pattern = "model_parameters_.*\\.csv$", full.names = TRUE
)
model_tables <- lapply(model_files, read.csv)
names(model_tables) <- tools::file_path_sans_ext(basename(model_files))

# frequency tables
freq_files <- list.files(
  path = extModelDir, pattern = "patch_area_frequency_logbins_.*\\.csv$", full.names = TRUE
)
freq_tables <- lapply(freq_files, read.csv)
names(freq_tables) <- tools::file_path_sans_ext(basename(freq_files))

# organize
lognormal_over  <- model_tables$lognormal_model_parameters_overfiltered
lognormal_under <- model_tables$lognormal_model_parameters_underfiltered
freq_over       <- freq_tables$patch_area_frequency_logbins_overfiltered
freq_under      <- freq_tables$patch_area_frequency_logbins_underfiltered
freq_boundary   <- freq_tables$patch_area_frequency_logbins_boundary_filtered

# total deforested area (from terra::cellSize computation)
total_deforested_km2 <- 454509.4


############################
# IMPORT ANIMAL DATA & GRID
############################

if (exists(".amani_animal_model_dir")) {
  animalModelDir <- .amani_animal_model_dir
} else {
  animalModelDir <- paste0(dirname(getwd()), "/2-ModellingAnimals/output/madre de dios analysis results")
}
sim_files <- list.files(
  path = animalModelDir, pattern = "^simulatedOccurencesForAnimal.*\\.rds$", full.names = TRUE
)
simAnimalObj <- readRDS(sim_files[1])

dataDir <- paste0(dirname(getwd()), "/1-FilteringPlantData/output")
occ_sf <- readRDS(paste0(dataDir, "/species_counts4_filtered.rds"))

grid_size_m <- 22360
grid <- st_make_grid(occ_sf, cellsize = grid_size_m, what = "polygons")
grid_sf <- st_sf(grid_id = seq_along(grid), geometry = grid)

if (st_crs(simAnimalObj) != st_crs(grid_sf)) {
  simAnimalObj <- st_transform(simAnimalObj, st_crs(grid_sf))
}

grid_sf$grid_area_m2  <- as.numeric(st_area(grid_sf))
grid_sf$grid_area_km2 <- grid_sf$grid_area_m2 / 1e6
grid_cell_size_km2 <- round(median(grid_sf$grid_area_km2))
stopifnot(length(grid_cell_size_km2) == 1)

cat("Grid cell size:", grid_cell_size_km2, "km2\n")
cat("Total grid cells:", nrow(grid_sf), "\n")
cat("Total species:", length(unique(simAnimalObj$scientificName)), "\n")


############################
# COMPUTE SPATIAL DEFORESTATION WEIGHTS
############################
# Strategy: use the patch centroids from the boundary-filtered extraction
# (3-ExtPredictor output) to estimate where deforestation is concentrated.
# For each grid cell, count how many patch centroids fall within it (or
# nearby). Grids with more/larger nearby patches get higher weight.
#
# This makes deforestation spatially realistic: grid cells in the "arc of
# deforestation" are much more likely to be hit than remote interior cells.

# load patch centroids from the boundary-filtered CSV files
patch_csv_dir <- paste0(dirname(getwd()), "/3-ExtPredictor/output/deforestTileOutStats")
patch_csvs <- list.files(
  path = patch_csv_dir,
  pattern = "^w[0-9]+\\.patch_metrics\\.csv$",
  full.names = TRUE
)

# read all patch CSVs and combine
all_patch_data <- do.call(rbind, lapply(patch_csvs, read.csv))
cat("Loaded", nrow(all_patch_data), "raw patch records\n")

# create sf from centroids (lon/lat WGS84, despite extraction using EPSG:5880
# for area calculations, the centroids were written in geographic coordinates)
patch_pts <- st_as_sf(
  all_patch_data,
  coords = c("centroid_x", "centroid_y"),
  crs = 4326
)
patch_pts <- st_transform(patch_pts, st_crs(grid_sf))

# spatial join: assign each patch centroid to a grid cell
patch_in_grid <- st_join(patch_pts, grid_sf, join = st_within)

# compute weight per grid cell: sum of patch areas within each cell
# this weights by total deforested area, not just count
grid_defor_weight <- patch_in_grid |>
  st_drop_geometry() |>
  group_by(grid_id) |>
  summarise(
    n_patches      = n(),
    total_area_km2 = sum(area_km2, na.rm = TRUE),
    .groups = "drop"
  )

# merge onto grid and fill zeros
grid_sf <- grid_sf |>
  left_join(grid_defor_weight, by = "grid_id")
grid_sf$total_area_km2[is.na(grid_sf$total_area_km2)] <- 0
grid_sf$n_patches[is.na(grid_sf$n_patches)] <- 0

# build the named weight vector
# use total_area_km2 as the weight (grids with more deforested area
# are more likely to be deforested further)
# add a small floor so every grid has nonzero probability
weight_floor <- 0.01 * mean(grid_sf$total_area_km2[grid_sf$total_area_km2 > 0])
grid_weights <- setNames(
  grid_sf$total_area_km2 + weight_floor,
  grid_sf$grid_id
)

cat("\nDeforestation weight summary:\n")
cat("Grids with deforestation:", sum(grid_sf$total_area_km2 > 0), "\n")
cat("Grids without:", sum(grid_sf$total_area_km2 == 0), "\n")
cat("Max weight:", round(max(grid_weights), 1), "\n")
cat("Min weight (floor):", round(min(grid_weights), 4), "\n")
cat("Weight ratio (max/min):", round(max(grid_weights) / min(grid_weights), 0), "\n")


############################
# TOROIDAL SHIFT SETUP
############################
# Build the parameter list needed by toroidal_shift_weights().
# The grid is 148 rows x 175 cols (column-major from st_make_grid).

nrow_grid <- round((st_bbox(grid_sf)["ymax"] - st_bbox(grid_sf)["ymin"]) / grid_size_m)
ncol_grid <- round((st_bbox(grid_sf)["xmax"] - st_bbox(grid_sf)["xmin"]) / grid_size_m)
cat("Grid dimensions for toroidal shift:", nrow_grid, "rows x", ncol_grid, "cols\n")
stopifnot(nrow_grid * ncol_grid == nrow(grid_sf))

# grid_ids in the order st_make_grid produces them (1:n, column-major)
grid_ids_ordered <- grid_sf$grid_id

# occupied grids = those that contain at least one species record
occupied_ids <- unique(simAnimalObj$grid_id)
occupied_ids <- occupied_ids[!is.na(occupied_ids)]

toroidal_params <- list(
  nrow_grid    = nrow_grid,
  ncol_grid    = ncol_grid,
  grid_ids     = grid_ids_ordered,
  occupied_ids = occupied_ids
)

# quick visual test of one shifted weight surface
set.seed(3841)
shifted_test <- toroidal_shift_weights(
  grid_weights, nrow_grid, ncol_grid, grid_ids_ordered, occupied_ids
)
grid_sf$shifted_test <- shifted_test[as.character(grid_sf$grid_id)]

cat("Shifted weight summary:\n")
cat("  Nonzero grids:", sum(shifted_test > 0), "\n")
cat("  Max:", round(max(shifted_test), 1), "\n")


############################
# VISUALIZE DEFORESTATION WEIGHTS
############################

ggplot(grid_sf) +
  geom_sf(aes(fill = total_area_km2), color = NA) +
  scale_fill_viridis_c(
    option = "inferno", direction = -1,
    na.value = "grey90",
    name = "Deforested area (km²)"
  ) +
  theme_minimal() +
  labs(
    title = "Deforestation intensity per grid cell (observed)",
    subtitle = "Used as spatial weights for clustered simulation"
  )

ggplot(grid_sf) +
  geom_sf(aes(fill = shifted_test), color = NA) +
  scale_fill_viridis_c(
    option = "inferno", direction = -1,
    na.value = "grey90",
    name = "Shifted weight"
  ) +
  theme_minimal() +
  labs(
    title = "Toroidal-shifted deforestation weights (one example)",
    subtitle = "Same spatial pattern, dissociated from real-world location"
  )


############################
# SINGLE TEST SIMULATION (clustered)
############################

out_test <- simulateExtinctionFromDeforestation(
  sim_tbl            = simAnimalObj,
  freq_table         = freq_over,
  grid_cell_size_km2 = grid_cell_size_km2,
  grid_weights       = grid_weights
)

cat("\nClustered test (freq_over):\n")
cat("Pop lost:", out_test$n_pop_lost_total, "\n")
cat("Spp lost:", out_test$n_spp_lost_total, "\n")


############################
# FULL PARALLEL RUNS (clustered)
############################

#######FREQ OVER#####

results_file_over <- list.files(
  path = outdir, pattern = "^clustered_FreqOver_.*\\.RData$", full.names = TRUE
)

if (length(results_file_over) > 0) {
  cat("Found existing clustered FreqOver results:\n")
  cat(paste(" ", results_file_over, collapse = "\n"), "\n")
  cat("Loading most recent file.\n")
  cat("To re-run, delete the file and re-source this section.\n")
  load(results_file_over[length(results_file_over)])
} else {
  clustered_results_FreqOver <- SEFD_rep_parallel(
    simAnimalObj,
    freq_table = freq_over,
    grid_cell_size_km2,
    n_runs = 100,
    grid_weights = grid_weights
  )
  save(
    clustered_results_FreqOver,
    file = paste0(
      outdir, "/clustered_FreqOver_",
      format(Sys.time(), "%Y%m%d_%H%M%S"),
      ".RData"
    )
  )
}


#######FREQ UNDER#####

results_file_under <- list.files(
  path = outdir, pattern = "^clustered_FreqUnder_.*\\.RData$", full.names = TRUE
)

if (length(results_file_under) > 0) {
  cat("Found existing clustered FreqUnder results:\n")
  cat(paste(" ", results_file_under, collapse = "\n"), "\n")
  cat("Loading most recent file.\n")
  load(results_file_under[length(results_file_under)])
} else {
  clustered_results_FreqUnder <- SEFD_rep_parallel(
    simAnimalObj,
    freq_table = freq_under,
    grid_cell_size_km2,
    n_runs = 100,
    grid_weights = grid_weights
  )
  save(
    clustered_results_FreqUnder,
    file = paste0(
      outdir, "/clustered_FreqUnder_",
      format(Sys.time(), "%Y%m%d_%H%M%S"),
      ".RData"
    )
  )
}


#######FREQ BOUNDARY#####

results_file_boundary <- list.files(
  path = outdir, pattern = "^clustered_FreqBoundary_.*\\.RData$", full.names = TRUE
)

if (length(results_file_boundary) > 0) {
  cat("Found existing clustered FreqBoundary results:\n")
  cat(paste(" ", results_file_boundary, collapse = "\n"), "\n")
  cat("Loading most recent file.\n")
  load(results_file_boundary[length(results_file_boundary)])
} else {
  clustered_results_FreqBoundary <- SEFD_rep_parallel(
    simAnimalObj,
    freq_table = freq_boundary,
    grid_cell_size_km2,
    n_runs = 100,
    grid_weights = grid_weights
  )
  save(
    clustered_results_FreqBoundary,
    file = paste0(
      outdir, "/clustered_FreqBoundary_",
      format(Sys.time(), "%Y%m%d_%H%M%S"),
      ".RData"
    )
  )
}


############################
# COMPARE ALL THREE CLUSTERED SCENARIOS
############################

n_spp_total <- length(unique(simAnimalObj$scientificName))

clustered_compare <- data.frame(
  scenario = c("FreqOver (conservative)", "FreqBoundary (middle)", "FreqUnder (liberal)"),
  total_patch_area_km2 = c(
    sum(freq_over$bin_mid * freq_over$count),
    sum(freq_boundary$bin_mid * freq_boundary$count),
    sum(freq_under$bin_mid * freq_under$count)
  ),
  median_spp_lost = c(
    median(clustered_results_FreqOver$n_spp_lost_total),
    median(clustered_results_FreqBoundary$n_spp_lost_total),
    median(clustered_results_FreqUnder$n_spp_lost_total)
  ),
  sd_spp_lost = round(c(
    sd(clustered_results_FreqOver$n_spp_lost_total),
    sd(clustered_results_FreqBoundary$n_spp_lost_total),
    sd(clustered_results_FreqUnder$n_spp_lost_total)
  ), 1),
  median_pop_lost = c(
    median(clustered_results_FreqOver$n_pop_lost_total),
    median(clustered_results_FreqBoundary$n_pop_lost_total),
    median(clustered_results_FreqUnder$n_pop_lost_total)
  ),
  pct_spp_lost = round(100 * c(
    median(clustered_results_FreqOver$n_spp_lost_total),
    median(clustered_results_FreqBoundary$n_spp_lost_total),
    median(clustered_results_FreqUnder$n_spp_lost_total)
  ) / n_spp_total, 1)
)

cat("\n=== CLUSTERED DEFORESTATION RESULTS ===\n")
clustered_compare


############################
# TOROIDAL-SHIFTED RUNS (null model)
############################
# Each replicate uses a different random toroidal shift of the weight
# surface. This preserves the spatial autocorrelation of the deforestation
# pattern but places it at a random location relative to species ranges.

#######FREQ OVER (shifted)#####

results_file_shifted_over <- list.files(
  path = outdir, pattern = "^shifted_FreqOver_.*\\.RData$", full.names = TRUE
)

if (length(results_file_shifted_over) > 0) {
  cat("Found existing shifted FreqOver results:\n")
  cat(paste(" ", results_file_shifted_over, collapse = "\n"), "\n")
  cat("Loading most recent file.\n")
  load(results_file_shifted_over[length(results_file_shifted_over)])
} else {
  shifted_results_FreqOver <- SEFD_rep_parallel(
    simAnimalObj,
    freq_table = freq_over,
    grid_cell_size_km2,
    n_runs = 100,
    grid_weights = grid_weights,
    toroidal_params = toroidal_params
  )
  save(
    shifted_results_FreqOver,
    file = paste0(
      outdir, "/shifted_FreqOver_",
      format(Sys.time(), "%Y%m%d_%H%M%S"),
      ".RData"
    )
  )
}

#######FREQ UNDER (shifted)#####

results_file_shifted_under <- list.files(
  path = outdir, pattern = "^shifted_FreqUnder_.*\\.RData$", full.names = TRUE
)

if (length(results_file_shifted_under) > 0) {
  cat("Found existing shifted FreqUnder results:\n")
  cat(paste(" ", results_file_shifted_under, collapse = "\n"), "\n")
  cat("Loading most recent file.\n")
  load(results_file_shifted_under[length(results_file_shifted_under)])
} else {
  shifted_results_FreqUnder <- SEFD_rep_parallel(
    simAnimalObj,
    freq_table = freq_under,
    grid_cell_size_km2,
    n_runs = 100,
    grid_weights = grid_weights,
    toroidal_params = toroidal_params
  )
  save(
    shifted_results_FreqUnder,
    file = paste0(
      outdir, "/shifted_FreqUnder_",
      format(Sys.time(), "%Y%m%d_%H%M%S"),
      ".RData"
    )
  )
}

#######FREQ BOUNDARY (shifted)#####

results_file_shifted_boundary <- list.files(
  path = outdir, pattern = "^shifted_FreqBoundary_.*\\.RData$", full.names = TRUE
)

if (length(results_file_shifted_boundary) > 0) {
  cat("Found existing shifted FreqBoundary results:\n")
  cat(paste(" ", results_file_shifted_boundary, collapse = "\n"), "\n")
  cat("Loading most recent file.\n")
  load(results_file_shifted_boundary[length(results_file_shifted_boundary)])
} else {
  shifted_results_FreqBoundary <- SEFD_rep_parallel(
    simAnimalObj,
    freq_table = freq_boundary,
    grid_cell_size_km2,
    n_runs = 100,
    grid_weights = grid_weights,
    toroidal_params = toroidal_params
  )
  save(
    shifted_results_FreqBoundary,
    file = paste0(
      outdir, "/shifted_FreqBoundary_",
      format(Sys.time(), "%Y%m%d_%H%M%S"),
      ".RData"
    )
  )
}


############################
# COMPARE: CLUSTERED VS SHIFTED VS NEUTRAL
############################

# helper to build one row of the comparison table
make_row <- function(label, res) {
  data.frame(
    model           = label,
    median_spp_lost = median(res$n_spp_lost_total),
    sd_spp_lost     = round(sd(res$n_spp_lost_total), 1),
    median_pop_lost = median(res$n_pop_lost_total),
    pct_spp_lost    = round(100 * median(res$n_spp_lost_total) / n_spp_total, 1)
  )
}

all_compare <- rbind(
  make_row("Clustered - FreqOver",     clustered_results_FreqOver),
  make_row("Shifted - FreqOver",       shifted_results_FreqOver),
  make_row("Clustered - FreqBoundary", clustered_results_FreqBoundary),
  make_row("Shifted - FreqBoundary",   shifted_results_FreqBoundary),
  make_row("Clustered - FreqUnder",    clustered_results_FreqUnder),
  make_row("Shifted - FreqUnder",      shifted_results_FreqUnder)
)

cat("\n=== CLUSTERED vs SHIFTED (toroidal null model) ===\n")
all_compare


############################
# INVESTIGATE EXTINCT SPECIES
############################

sim_tbl <- sf::st_drop_geometry(simAnimalObj)

set.seed(7293)
out_over <- simulateExtinctionFromDeforestation(
  sim_tbl = sim_tbl, freq_table = freq_over,
  grid_cell_size_km2 = grid_cell_size_km2, grid_weights = grid_weights
)
out_boundary <- simulateExtinctionFromDeforestation(
  sim_tbl = sim_tbl, freq_table = freq_boundary,
  grid_cell_size_km2 = grid_cell_size_km2, grid_weights = grid_weights
)
out_under <- simulateExtinctionFromDeforestation(
  sim_tbl = sim_tbl, freq_table = freq_under,
  grid_cell_size_km2 = grid_cell_size_km2, grid_weights = grid_weights
)

range_sizes <- sim_tbl |>
  group_by(scientificName) |>
  summarise(
    n_grids = n_distinct(grid_id),
    n_pops  = n(),
    .groups = "drop"
  )

surv_over     <- unique(out_over$final_table$scientificName)
surv_boundary <- unique(out_boundary$final_table$scientificName)
surv_under    <- unique(out_under$final_table$scientificName)

range_sizes$extinct_over     <- !(range_sizes$scientificName %in% surv_over)
range_sizes$extinct_boundary <- !(range_sizes$scientificName %in% surv_boundary)
range_sizes$extinct_under    <- !(range_sizes$scientificName %in% surv_under)

cat("\n=== Range sizes of extinct vs surviving species (clustered) ===\n")
for (scen in c("over", "boundary", "under")) {
  col <- paste0("extinct_", scen)
  cat("\n--- Freq", scen, "---\n")
  cat("Extinct (n =", sum(range_sizes[[col]]), "):\n")
  print(summary(range_sizes$n_grids[range_sizes[[col]]]))
  cat("Surviving:\n")
  print(summary(range_sizes$n_grids[!range_sizes[[col]]]))
}

n_single_grid <- sum(range_sizes$n_grids == 1)
cat("\n=== Single-grid endemic species ===\n")
cat("Total:", n_single_grid, "of", nrow(range_sizes),
    "(", round(100 * n_single_grid / nrow(range_sizes), 1), "% )\n\n")

cat("% of single-grid spp extinct per scenario:\n")
cat("  FreqOver:    ",
    round(100 * sum(range_sizes$extinct_over & range_sizes$n_grids == 1) / n_single_grid, 1), "%\n")
cat("  FreqBoundary:",
    round(100 * sum(range_sizes$extinct_boundary & range_sizes$n_grids == 1) / n_single_grid, 1), "%\n")
cat("  FreqUnder:   ",
    round(100 * sum(range_sizes$extinct_under & range_sizes$n_grids == 1) / n_single_grid, 1), "%\n")

cat("\nMulti-grid species that went extinct:\n")
cat("  FreqOver:    ", sum(range_sizes$extinct_over & range_sizes$n_grids > 1), "\n")
cat("  FreqBoundary:", sum(range_sizes$extinct_boundary & range_sizes$n_grids > 1), "\n")
cat("  FreqUnder:   ", sum(range_sizes$extinct_under & range_sizes$n_grids > 1), "\n")


############################
# PER-SPECIES RANGE REDUCTION (single detailed clustered runs)
############################

baseline_ranges <- precompute_baseline_ranges(sim_tbl)

rr_over     <- compute_range_reduction(out_over$final_table,     baseline_ranges)
rr_boundary <- compute_range_reduction(out_boundary$final_table, baseline_ranges)
rr_under    <- compute_range_reduction(out_under$final_table,    baseline_ranges)

rr_over$scenario     <- "FreqOver"
rr_boundary$scenario <- "FreqBoundary"
rr_under$scenario    <- "FreqUnder"

clustered_range_reduction <- rbind(rr_over, rr_boundary, rr_under)

write.csv(clustered_range_reduction,
          paste0(outdir, "/clustered_species_range_reduction.csv"),
          row.names = FALSE)
cat("Saved output/clustered_species_range_reduction.csv\n")

cat("\n=== Range Reduction Thresholds (clustered, single run) ===\n")
for (scen in c("FreqOver", "FreqBoundary", "FreqUnder")) {
  rr_sub <- clustered_range_reduction[clustered_range_reduction$scenario == scen, ]
  tc <- count_range_thresholds(rr_sub)
  cat("\n", scen, ":\n")
  for (j in seq_along(tc)) {
    cat("  ", names(tc)[j], ": ", tc[j], "\n")
  }
}


############################
# ENDEMISM-BASED WEIGHTS
############################
# Derive spatial weights from the input plant occurrence data (occ_sf) rather
# than from observed deforestation. Weighted endemism (WE) per cell is the sum
# of 1/range_size across all species present, where range_size is the number of
# grid cells a species occupies. Cells with many range-restricted species score
# highest and are deforested preferentially in the clustered-endemism model.

plant_ranges <- st_drop_geometry(occ_sf) |>
  group_by(scientificName) |>
  summarise(range_size = n_distinct(grid_id), .groups = "drop")

endemism_per_cell <- st_drop_geometry(occ_sf) |>
  left_join(plant_ranges, by = "scientificName") |>
  group_by(grid_id) |>
  summarise(weighted_endemism = sum(1 / range_size), .groups = "drop")

grid_sf <- grid_sf |>
  left_join(endemism_per_cell, by = "grid_id")
grid_sf$weighted_endemism[is.na(grid_sf$weighted_endemism)] <- 0

endemism_floor <- 0.01 * mean(grid_sf$weighted_endemism[grid_sf$weighted_endemism > 0])
grid_weights_endemism <- setNames(
  grid_sf$weighted_endemism + endemism_floor,
  grid_sf$grid_id
)

cat("\nEndemism weight summary:\n")
cat("  Cells with endemism data:", sum(grid_sf$weighted_endemism > 0), "\n")
cat("  Max WE:", round(max(grid_weights_endemism), 3), "\n")
cat("  Weight ratio (max/min):", round(max(grid_weights_endemism) / min(grid_weights_endemism), 0), "\n")

toroidal_params_endemism <- list(
  nrow_grid    = nrow_grid,
  ncol_grid    = ncol_grid,
  grid_ids     = grid_ids_ordered,
  occupied_ids = occupied_ids
)


############################
# FULL PARALLEL RUNS (clustered-endemism)
############################

#######FREQ OVER#####

results_file_endemism_clustered_over <- list.files(
  path = outdir, pattern = "^endemism_clustered_FreqOver_.*\\.RData$", full.names = TRUE
)
if (length(results_file_endemism_clustered_over) > 0) {
  cat("Loading existing endemism_clustered FreqOver results.\n")
  load(results_file_endemism_clustered_over[length(results_file_endemism_clustered_over)])
} else {
  endemism_clustered_FreqOver <- SEFD_rep_parallel(
    simAnimalObj, freq_over, grid_cell_size_km2,
    n_runs       = 100,
    grid_weights = grid_weights_endemism
  )
  save(endemism_clustered_FreqOver,
       file = paste0(outdir, "/endemism_clustered_FreqOver_",
                     format(Sys.time(), "%Y%m%d_%H%M%S"), ".RData"))
}

#######FREQ BOUNDARY#####

results_file_endemism_clustered_boundary <- list.files(
  path = outdir, pattern = "^endemism_clustered_FreqBoundary_.*\\.RData$", full.names = TRUE
)
if (length(results_file_endemism_clustered_boundary) > 0) {
  cat("Loading existing endemism_clustered FreqBoundary results.\n")
  load(results_file_endemism_clustered_boundary[length(results_file_endemism_clustered_boundary)])
} else {
  endemism_clustered_FreqBoundary <- SEFD_rep_parallel(
    simAnimalObj, freq_boundary, grid_cell_size_km2,
    n_runs       = 100,
    grid_weights = grid_weights_endemism
  )
  save(endemism_clustered_FreqBoundary,
       file = paste0(outdir, "/endemism_clustered_FreqBoundary_",
                     format(Sys.time(), "%Y%m%d_%H%M%S"), ".RData"))
}

#######FREQ UNDER#####

results_file_endemism_clustered_under <- list.files(
  path = outdir, pattern = "^endemism_clustered_FreqUnder_.*\\.RData$", full.names = TRUE
)
if (length(results_file_endemism_clustered_under) > 0) {
  cat("Loading existing endemism_clustered FreqUnder results.\n")
  load(results_file_endemism_clustered_under[length(results_file_endemism_clustered_under)])
} else {
  endemism_clustered_FreqUnder <- SEFD_rep_parallel(
    simAnimalObj, freq_under, grid_cell_size_km2,
    n_runs       = 100,
    grid_weights = grid_weights_endemism
  )
  save(endemism_clustered_FreqUnder,
       file = paste0(outdir, "/endemism_clustered_FreqUnder_",
                     format(Sys.time(), "%Y%m%d_%H%M%S"), ".RData"))
}


############################
# FULL PARALLEL RUNS (shifted-endemism)
############################

#######FREQ OVER#####

results_file_endemism_shifted_over <- list.files(
  path = outdir, pattern = "^endemism_shifted_FreqOver_.*\\.RData$", full.names = TRUE
)
if (length(results_file_endemism_shifted_over) > 0) {
  cat("Loading existing endemism_shifted FreqOver results.\n")
  load(results_file_endemism_shifted_over[length(results_file_endemism_shifted_over)])
} else {
  endemism_shifted_FreqOver <- SEFD_rep_parallel(
    simAnimalObj, freq_over, grid_cell_size_km2,
    n_runs          = 100,
    grid_weights    = grid_weights_endemism,
    toroidal_params = toroidal_params_endemism
  )
  save(endemism_shifted_FreqOver,
       file = paste0(outdir, "/endemism_shifted_FreqOver_",
                     format(Sys.time(), "%Y%m%d_%H%M%S"), ".RData"))
}

#######FREQ BOUNDARY#####

results_file_endemism_shifted_boundary <- list.files(
  path = outdir, pattern = "^endemism_shifted_FreqBoundary_.*\\.RData$", full.names = TRUE
)
if (length(results_file_endemism_shifted_boundary) > 0) {
  cat("Loading existing endemism_shifted FreqBoundary results.\n")
  load(results_file_endemism_shifted_boundary[length(results_file_endemism_shifted_boundary)])
} else {
  endemism_shifted_FreqBoundary <- SEFD_rep_parallel(
    simAnimalObj, freq_boundary, grid_cell_size_km2,
    n_runs          = 100,
    grid_weights    = grid_weights_endemism,
    toroidal_params = toroidal_params_endemism
  )
  save(endemism_shifted_FreqBoundary,
       file = paste0(outdir, "/endemism_shifted_FreqBoundary_",
                     format(Sys.time(), "%Y%m%d_%H%M%S"), ".RData"))
}

#######FREQ UNDER#####

results_file_endemism_shifted_under <- list.files(
  path = outdir, pattern = "^endemism_shifted_FreqUnder_.*\\.RData$", full.names = TRUE
)
if (length(results_file_endemism_shifted_under) > 0) {
  cat("Loading existing endemism_shifted FreqUnder results.\n")
  load(results_file_endemism_shifted_under[length(results_file_endemism_shifted_under)])
} else {
  endemism_shifted_FreqUnder <- SEFD_rep_parallel(
    simAnimalObj, freq_under, grid_cell_size_km2,
    n_runs          = 100,
    grid_weights    = grid_weights_endemism,
    toroidal_params = toroidal_params_endemism
  )
  save(endemism_shifted_FreqUnder,
       file = paste0(outdir, "/endemism_shifted_FreqUnder_",
                     format(Sys.time(), "%Y%m%d_%H%M%S"), ".RData"))
}


############################
# FREQ-CALIBRATED FREQUENCY TABLE
############################
# Scale FreqOver patch counts so the total simulated area matches the raster-
# derived total deforested area. The patch size distribution shape is preserved;
# only the overall volume of deforestation is adjusted.

freq_calibrated_area <- sum(freq_over$bin_mid * freq_over$count)
freq_calibrated_scale <- total_deforested_km2 / freq_calibrated_area
freq_calibrated       <- freq_over
freq_calibrated$count <- round(freq_over$count * freq_calibrated_scale)

cat("\nFreqCalibrated summary:\n")
cat("  Scale factor:          ", round(freq_calibrated_scale, 4), "\n")
cat("  Original area (km2):   ", round(freq_calibrated_area, 1), "\n")
cat("  Target area (km2):     ", total_deforested_km2, "\n")
cat("  Calibrated area (km2): ", round(sum(freq_calibrated$bin_mid * freq_calibrated$count), 1), "\n")


############################
# FULL PARALLEL RUNS (FreqCalibrated — all models)
############################

#######NEUTRAL#####

results_file_neutral_calibrated <- list.files(
  path = outdir, pattern = "^neutral_FreqCalibrated_.*\\.RData$", full.names = TRUE
)
if (length(results_file_neutral_calibrated) > 0) {
  cat("Loading existing neutral FreqCalibrated results.\n")
  load(results_file_neutral_calibrated[length(results_file_neutral_calibrated)])
} else {
  neutral_results_FreqCalibrated <- SEFD_rep_parallel(
    simAnimalObj, freq_calibrated, grid_cell_size_km2,
    n_runs       = 100,
    grid_weights = NULL
  )
  save(neutral_results_FreqCalibrated,
       file = paste0(outdir, "/neutral_FreqCalibrated_",
                     format(Sys.time(), "%Y%m%d_%H%M%S"), ".RData"))
}

#######CLUSTERED#####

results_file_clustered_calibrated <- list.files(
  path = outdir, pattern = "^clustered_FreqCalibrated_.*\\.RData$", full.names = TRUE
)
if (length(results_file_clustered_calibrated) > 0) {
  cat("Loading existing clustered FreqCalibrated results.\n")
  load(results_file_clustered_calibrated[length(results_file_clustered_calibrated)])
} else {
  clustered_results_FreqCalibrated <- SEFD_rep_parallel(
    simAnimalObj, freq_calibrated, grid_cell_size_km2,
    n_runs       = 100,
    grid_weights = grid_weights
  )
  save(clustered_results_FreqCalibrated,
       file = paste0(outdir, "/clustered_FreqCalibrated_",
                     format(Sys.time(), "%Y%m%d_%H%M%S"), ".RData"))
}

#######SHIFTED#####

results_file_shifted_calibrated <- list.files(
  path = outdir, pattern = "^shifted_FreqCalibrated_.*\\.RData$", full.names = TRUE
)
if (length(results_file_shifted_calibrated) > 0) {
  cat("Loading existing shifted FreqCalibrated results.\n")
  load(results_file_shifted_calibrated[length(results_file_shifted_calibrated)])
} else {
  shifted_results_FreqCalibrated <- SEFD_rep_parallel(
    simAnimalObj, freq_calibrated, grid_cell_size_km2,
    n_runs          = 100,
    grid_weights    = grid_weights,
    toroidal_params = toroidal_params
  )
  save(shifted_results_FreqCalibrated,
       file = paste0(outdir, "/shifted_FreqCalibrated_",
                     format(Sys.time(), "%Y%m%d_%H%M%S"), ".RData"))
}

#######CLUSTERED-ENDEMISM#####

results_file_endemism_clustered_calibrated <- list.files(
  path = outdir, pattern = "^endemism_clustered_FreqCalibrated_.*\\.RData$", full.names = TRUE
)
if (length(results_file_endemism_clustered_calibrated) > 0) {
  cat("Loading existing endemism_clustered FreqCalibrated results.\n")
  load(results_file_endemism_clustered_calibrated[length(results_file_endemism_clustered_calibrated)])
} else {
  endemism_clustered_FreqCalibrated <- SEFD_rep_parallel(
    simAnimalObj, freq_calibrated, grid_cell_size_km2,
    n_runs       = 100,
    grid_weights = grid_weights_endemism
  )
  save(endemism_clustered_FreqCalibrated,
       file = paste0(outdir, "/endemism_clustered_FreqCalibrated_",
                     format(Sys.time(), "%Y%m%d_%H%M%S"), ".RData"))
}

#######SHIFTED-ENDEMISM#####

results_file_endemism_shifted_calibrated <- list.files(
  path = outdir, pattern = "^endemism_shifted_FreqCalibrated_.*\\.RData$", full.names = TRUE
)
if (length(results_file_endemism_shifted_calibrated) > 0) {
  cat("Loading existing endemism_shifted FreqCalibrated results.\n")
  load(results_file_endemism_shifted_calibrated[length(results_file_endemism_shifted_calibrated)])
} else {
  endemism_shifted_FreqCalibrated <- SEFD_rep_parallel(
    simAnimalObj, freq_calibrated, grid_cell_size_km2,
    n_runs          = 100,
    grid_weights    = grid_weights_endemism,
    toroidal_params = toroidal_params_endemism
  )
  save(endemism_shifted_FreqCalibrated,
       file = paste0(outdir, "/endemism_shifted_FreqCalibrated_",
                     format(Sys.time(), "%Y%m%d_%H%M%S"), ".RData"))
}


############################
# EXPORT ALL RESULTS TO CSV
############################

n_spp_total <- length(unique(sim_tbl$scientificName))

make_summary_row <- function(model, scenario, res) {
  data.frame(
    model           = model,
    scenario        = scenario,
    median_spp_lost = median(res$n_spp_lost_total),
    mean_spp_lost   = round(mean(res$n_spp_lost_total), 1),
    sd_spp_lost     = round(sd(res$n_spp_lost_total), 1),
    q25_spp_lost    = unname(quantile(res$n_spp_lost_total, 0.25)),
    q75_spp_lost    = unname(quantile(res$n_spp_lost_total, 0.75)),
    median_pop_lost = median(res$n_pop_lost_total),
    mean_pop_lost   = round(mean(res$n_pop_lost_total), 1)
  )
}

tag <- function(df, scenario) { df$scenario <- scenario; df }

# 1. Master summary — all models × all scenarios
summary_all_models <- rbind(
  make_summary_row("Neutral",            "FreqOver",       extinction_results_FreqOver),
  make_summary_row("Neutral",            "FreqBoundary",   extinction_results_FreqBoundary),
  make_summary_row("Neutral",            "FreqUnder",      extinction_results_FreqUnder),
  make_summary_row("Neutral",            "FreqCalibrated", neutral_results_FreqCalibrated),
  make_summary_row("Clustered",          "FreqOver",       clustered_results_FreqOver),
  make_summary_row("Clustered",          "FreqBoundary",   clustered_results_FreqBoundary),
  make_summary_row("Clustered",          "FreqUnder",      clustered_results_FreqUnder),
  make_summary_row("Clustered",          "FreqCalibrated", clustered_results_FreqCalibrated),
  make_summary_row("Shifted",            "FreqOver",       shifted_results_FreqOver),
  make_summary_row("Shifted",            "FreqBoundary",   shifted_results_FreqBoundary),
  make_summary_row("Shifted",            "FreqUnder",      shifted_results_FreqUnder),
  make_summary_row("Shifted",            "FreqCalibrated", shifted_results_FreqCalibrated),
  make_summary_row("Clustered_endemism", "FreqOver",       endemism_clustered_FreqOver),
  make_summary_row("Clustered_endemism", "FreqBoundary",   endemism_clustered_FreqBoundary),
  make_summary_row("Clustered_endemism", "FreqUnder",      endemism_clustered_FreqUnder),
  make_summary_row("Clustered_endemism", "FreqCalibrated", endemism_clustered_FreqCalibrated),
  make_summary_row("Shifted_endemism",   "FreqOver",       endemism_shifted_FreqOver),
  make_summary_row("Shifted_endemism",   "FreqBoundary",   endemism_shifted_FreqBoundary),
  make_summary_row("Shifted_endemism",   "FreqUnder",      endemism_shifted_FreqUnder),
  make_summary_row("Shifted_endemism",   "FreqCalibrated", endemism_shifted_FreqCalibrated)
)
write.csv(summary_all_models, paste0(outdir, "/summary_all_models.csv"), row.names = FALSE)
cat("Saved output/summary_all_models.csv\n")
print(summary_all_models, row.names = FALSE)

# 2. Raw replicates — neutral
write.csv(
  rbind(tag(extinction_results_FreqOver,     "FreqOver"),
        tag(extinction_results_FreqBoundary, "FreqBoundary"),
        tag(extinction_results_FreqUnder,    "FreqUnder"),
        tag(neutral_results_FreqCalibrated,  "FreqCalibrated")),
  paste0(outdir, "/neutral_all_replicates.csv"), row.names = FALSE
)
cat("Saved output/neutral_all_replicates.csv\n")

# 3. Raw replicates — clustered
write.csv(
  rbind(tag(clustered_results_FreqOver,       "FreqOver"),
        tag(clustered_results_FreqBoundary,   "FreqBoundary"),
        tag(clustered_results_FreqUnder,      "FreqUnder"),
        tag(clustered_results_FreqCalibrated, "FreqCalibrated")),
  paste0(outdir, "/clustered_all_replicates.csv"), row.names = FALSE
)
cat("Saved output/clustered_all_replicates.csv\n")

# 4. Raw replicates — shifted
write.csv(
  rbind(tag(shifted_results_FreqOver,       "FreqOver"),
        tag(shifted_results_FreqBoundary,   "FreqBoundary"),
        tag(shifted_results_FreqUnder,      "FreqUnder"),
        tag(shifted_results_FreqCalibrated, "FreqCalibrated")),
  paste0(outdir, "/shifted_all_replicates.csv"), row.names = FALSE
)
cat("Saved output/shifted_all_replicates.csv\n")

# 5. Raw replicates — clustered-endemism
write.csv(
  rbind(tag(endemism_clustered_FreqOver,       "FreqOver"),
        tag(endemism_clustered_FreqBoundary,   "FreqBoundary"),
        tag(endemism_clustered_FreqUnder,      "FreqUnder"),
        tag(endemism_clustered_FreqCalibrated, "FreqCalibrated")),
  paste0(outdir, "/endemism_clustered_all_replicates.csv"), row.names = FALSE
)
cat("Saved output/endemism_clustered_all_replicates.csv\n")

# 6. Raw replicates — shifted-endemism
write.csv(
  rbind(tag(endemism_shifted_FreqOver,       "FreqOver"),
        tag(endemism_shifted_FreqBoundary,   "FreqBoundary"),
        tag(endemism_shifted_FreqUnder,      "FreqUnder"),
        tag(endemism_shifted_FreqCalibrated, "FreqCalibrated")),
  paste0(outdir, "/endemism_shifted_all_replicates.csv"), row.names = FALSE
)
cat("Saved output/endemism_shifted_all_replicates.csv\n")

# 7. Species extinction status (from single detailed clustered run)
write.csv(range_sizes, paste0(outdir, "/clustered_species_extinction_status.csv"), row.names = FALSE)
cat("Saved output/clustered_species_extinction_status.csv\n")


############################
# DONE / STATUS
############################
# [DONE] Clustered deforestation model (FreqOver, FreqBoundary, FreqUnder, FreqCalibrated)
# [DONE] Toroidal-shifted null model (FreqOver, FreqBoundary, FreqUnder, FreqCalibrated)
# [DONE] Clustered-endemism model (FreqOver, FreqBoundary, FreqUnder, FreqCalibrated)
# [DONE] Shifted-endemism null model (FreqOver, FreqBoundary, FreqUnder, FreqCalibrated)
# [DONE] Investigate which species go extinct (range-restricted vs widespread)
#
# See NOTES.txt in the project root for outstanding items and next steps.
