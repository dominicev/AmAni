# =====================================================================
# 03_visualize_deforestation_maps.R — Diagnostic spatial visualization
# =====================================================================
# Runs a few single-replicate simulations for each deforestation model
# (neutral, clustered, shifted) and maps per-cell deforestation intensity
# overlaid on species distributions. Useful for visually comparing how
# the three placement models distribute habitat loss across the grid.
#
# Prerequisites: run the data-loading sections of 01_simulate_neutral.R
# or 02_simulate_clustered.R first (through weight computation), OR
# source this script standalone — it loads everything it needs.
# =====================================================================

{
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
}

source("scripts/utils/simulation_functions.R")

# Output directory (overridden by manager when running organism-specific pipelines)
if (exists(".amani_output_dir")) {
  outdir <- .amani_output_dir
} else {
  outdir <- "output"
}
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)


# =====================================================================
# DATA LOADING (mirrors 02_simulate_clustered.R)
# =====================================================================

extModelDir <- paste0(dirname(getwd()), "/3-ExtPredictor/output")

freq_files <- list.files(
  path = extModelDir, pattern = "patch_area_frequency_logbins_.*\\.csv$",
  full.names = TRUE
)
freq_tables <- lapply(freq_files, read.csv)
names(freq_tables) <- tools::file_path_sans_ext(basename(freq_files))
freq_over <- freq_tables$patch_area_frequency_logbins_overfiltered

# Animal data & grid
if (exists(".amani_animal_model_dir")) {
  animalModelDir <- .amani_animal_model_dir
} else {
  animalModelDir <- paste0(dirname(getwd()), "/2-ModellingAnimals/output/madre de dios analysis results")
}
sim_files <- list.files(
  path = animalModelDir, pattern = "^simulatedOccurencesForAnimal\\.rds$",
  full.names = TRUE
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

# Deforestation weights
patch_csv_dir <- paste0(dirname(getwd()), "/3-ExtPredictor/output/deforestTileOutStats")
patch_csvs <- list.files(
  path = patch_csv_dir,
  pattern = "^w[0-9]+\\.patch_metrics\\.csv$",
  full.names = TRUE
)

all_patch_data <- do.call(rbind, lapply(patch_csvs, read.csv))
patch_pts <- st_as_sf(all_patch_data, coords = c("centroid_x", "centroid_y"), crs = 4326)
patch_pts <- st_transform(patch_pts, st_crs(grid_sf))
patch_in_grid <- st_join(patch_pts, grid_sf, join = st_within)

grid_defor_weight <- patch_in_grid |>
  st_drop_geometry() |>
  group_by(grid_id) |>
  summarise(
    n_patches = n(),
    total_area_km2 = sum(area_km2, na.rm = TRUE),
    .groups = "drop"
  )

grid_sf <- grid_sf |>
  left_join(grid_defor_weight, by = "grid_id")
grid_sf$total_area_km2[is.na(grid_sf$total_area_km2)] <- 0
grid_sf$n_patches[is.na(grid_sf$n_patches)] <- 0

weight_floor <- 0.01 * mean(grid_sf$total_area_km2[grid_sf$total_area_km2 > 0])
grid_weights <- setNames(
  grid_sf$total_area_km2 + weight_floor,
  grid_sf$grid_id
)

# Toroidal shift parameters
nrow_grid <- round((st_bbox(grid_sf)["ymax"] - st_bbox(grid_sf)["ymin"]) / grid_size_m)
ncol_grid <- round((st_bbox(grid_sf)["xmax"] - st_bbox(grid_sf)["xmin"]) / grid_size_m)
occupied_ids <- unique(simAnimalObj$grid_id)
occupied_ids <- occupied_ids[!is.na(occupied_ids)]
grid_ids_ordered <- grid_sf$grid_id

sim_tbl <- sf::st_drop_geometry(simAnimalObj)

cat("Data loaded.\n")
cat("Grid:", nrow_grid, "rows x", ncol_grid, "cols,", nrow(grid_sf), "cells\n")
cat("Species:", length(unique(sim_tbl$scientificName)),
    "| Records:", nrow(sim_tbl),
    "| Occupied cells:", length(occupied_ids), "\n")


# =====================================================================
# DIAGNOSTIC SIMULATION WRAPPER
# =====================================================================
# Runs one simulation and returns per-cell deforestation intensity
# (0 = untouched, 1 = fully deforested). Does NOT modify source files.

simulate_with_cell_info <- function(
    sim_tbl, freq_table, grid_cell_size_km2, grid_weights = NULL
) {
  if (inherits(sim_tbl, "sf")) sim_tbl <- sf::st_drop_geometry(sim_tbl)

  # --- Large patches ---
  freq_large <- freq_table[
    freq_table$bin_mid >= grid_cell_size_km2 & freq_table$count > 0, , drop = FALSE
  ]

  available_grids <- unique(sim_tbl$grid_id)
  available_grids <- available_grids[!is.na(available_grids)]

  large_deforested_grids <- character(0)
  if (nrow(freq_large) > 0) {
    grid_equiv <- freq_large$count * (freq_large$bin_mid / grid_cell_size_km2)
    n_whole <- floor(grid_equiv)
    total_grids_to_remove <- sum(n_whole)

    if (!is.null(grid_weights)) {
      w <- grid_weights[as.character(available_grids)]
      w[is.na(w)] <- 0
      if (sum(w) > 0) w <- w / sum(w) else w <- rep(1 / length(available_grids), length(available_grids))
      n_eligible <- sum(w > 0)
      total_grids_to_remove <- min(total_grids_to_remove, n_eligible)
    } else {
      w <- NULL
    }
    total_grids_to_remove <- min(total_grids_to_remove, length(available_grids))

    if (total_grids_to_remove > 0) {
      large_deforested_grids <- as.character(
        sample(available_grids, total_grids_to_remove, replace = FALSE, prob = w)
      )
    }
  }

  current_tbl <- sim_tbl[!(sim_tbl$grid_id %in% large_deforested_grids), , drop = FALSE]

  # --- Small patches: compute per-cell survival probability ---
  freq_small <- freq_table[
    freq_table$bin_mid < grid_cell_size_km2 & freq_table$count > 0, , drop = FALSE
  ]

  unique_grids <- unique(current_tbl$grid_id)
  unique_grids <- unique_grids[!is.na(unique_grids)]
  n_grids <- length(unique_grids)

  if (!is.null(grid_weights)) {
    w <- grid_weights[as.character(unique_grids)]
    w[is.na(w)] <- 0
    if (sum(w) > 0) w <- w / sum(w) else w <- rep(1 / n_grids, n_grids)
  } else {
    w <- NULL
  }

  grid_survive_prob <- setNames(rep(1.0, n_grids), unique_grids)

  for (i in seq_len(nrow(freq_small))) {
    perc <- min(freq_small$bin_mid[i] / grid_cell_size_km2, 1)
    n_hits <- freq_small$count[i]
    hit_grids <- sample(unique_grids, size = n_hits, replace = TRUE, prob = w)
    hit_counts <- table(hit_grids)
    for (g in names(hit_counts)) {
      grid_survive_prob[g] <- grid_survive_prob[g] * (1 - perc)^hit_counts[g]
    }
  }

  # Combine: large patches get intensity = 1, small patches = 1 - survive_prob
  all_grids <- unique(as.character(sim_tbl$grid_id))
  all_grids <- all_grids[!is.na(all_grids)]

  defor_intensity <- setNames(rep(0, length(all_grids)), all_grids)
  for (g in names(grid_survive_prob)) {
    defor_intensity[g] <- 1 - grid_survive_prob[g]
  }
  defor_intensity[large_deforested_grids] <- 1.0

  data.frame(
    grid_id = as.integer(names(defor_intensity)),
    defor_intensity = as.numeric(defor_intensity),
    large_patch_removed = names(defor_intensity) %in% large_deforested_grids
  )
}


# =====================================================================
# RUN DIAGNOSTIC SIMULATIONS
# =====================================================================

n_reps <- 3

set.seed(4817)

# Neutral (uniform random — no grid_weights)
neutral_cells <- lapply(seq_len(n_reps), function(i) {
  res <- simulate_with_cell_info(sim_tbl, freq_over, grid_cell_size_km2, grid_weights = NULL)
  res$replicate <- i
  res$model <- "Neutral"
  res
})

# Clustered (real deforestation weights)
clustered_cells <- lapply(seq_len(n_reps), function(i) {
  res <- simulate_with_cell_info(sim_tbl, freq_over, grid_cell_size_km2, grid_weights)
  res$replicate <- i
  res$model <- "Clustered"
  res
})

# Shifted (toroidally shifted weights — spatial null)
shifted_cells <- lapply(seq_len(n_reps), function(i) {
  shifted_w <- toroidal_shift_weights(
    grid_weights, nrow_grid, ncol_grid, grid_ids_ordered, occupied_ids
  )
  res <- simulate_with_cell_info(sim_tbl, freq_over, grid_cell_size_km2, shifted_w)
  res$replicate <- i
  res$model <- "Shifted"
  res
})

all_cell_data <- do.call(rbind, c(neutral_cells, clustered_cells, shifted_cells))
all_cell_data$model <- factor(all_cell_data$model, levels = c("Neutral", "Shifted", "Clustered"))

cat("Simulations complete:", nrow(all_cell_data), "rows\n")


# =====================================================================
# PREPARE PLOT DATA
# =====================================================================

# Species richness per grid cell (background layer)
spp_richness <- sim_tbl |>
  group_by(grid_id) |>
  summarise(n_species = n_distinct(scientificName), .groups = "drop")

grid_occupied <- grid_sf |>
  inner_join(spp_richness, by = "grid_id")

# Merge deforestation intensity with grid geometry
all_cell_sf <- all_cell_data |>
  left_join(
    grid_sf |> select(grid_id) |> st_drop_geometry() |> cbind(geometry = st_geometry(grid_sf)),
    by = "grid_id"
  )
# Safer merge via grid_sf directly
all_cell_sf <- all_cell_data |>
  left_join(grid_sf |> select(grid_id, geometry), by = "grid_id") |>
  st_as_sf()

# Only plot cells with meaningful deforestation
defor_cells <- all_cell_sf |> filter(defor_intensity > 0.001)


# =====================================================================
# FACETED MAP: Neutral vs Shifted vs Clustered
# =====================================================================

p <- ggplot() +
  geom_sf(data = grid_occupied, fill = "grey85", color = "grey70", linewidth = 0.05) +
  geom_sf(data = defor_cells, aes(fill = defor_intensity), color = NA) +
  scale_fill_gradient(
    low = "#FEE08B", high = "#D73027",
    name = "Deforestation\nintensity",
    limits = c(0, 1)
  ) +
  facet_grid(model ~ replicate,
             labeller = labeller(replicate = function(x) paste("Rep", x))) +
  theme_minimal(base_size = 10) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank()
  ) +
  labs(
    title = "Simulated deforestation across three placement models (freq_over)",
    subtitle = "Grey = occupied cells; color = deforestation intensity (yellow → red)"
  )

ggsave(paste0(outdir, "/diagnostic_deforestation_maps.pdf"), p, width = 12, height = 10)


# =====================================================================
# DETERMINISTIC TRANSFORM DIAGNOSTIC — shifted model only
# =====================================================================
# One simulation per transformation type, each with a fixed half-grid
# shift where applicable, so results are fully reproducible and the
# effect of each transform can be assessed in isolation.

half_row <- floor(nrow_grid / 2)
half_col <- floor(ncol_grid / 2)

transform_cases <- list(
  list(label = "Reference\n(no shift)",      shift_row = 0,        shift_col = 0       ),
  list(label = "Shift rows\n(quarter)",       shift_row = floor(nrow_grid / 4), shift_col = 0       ),
  list(label = "Shift rows\n(half)",          shift_row = half_row, shift_col = 0       ),
  list(label = "Shift rows\n(three-quarter)", shift_row = floor(nrow_grid * 3/4), shift_col = 0     ),
  list(label = "Shift cols\n(quarter)",       shift_row = 0,        shift_col = floor(ncol_grid / 4)),
  list(label = "Shift cols\n(half)",          shift_row = 0,        shift_col = half_col            ),
  list(label = "Shift cols\n(three-quarter)", shift_row = 0,        shift_col = floor(ncol_grid * 3/4)),
  list(label = "Shift rows + cols\n(quarter)",shift_row = floor(nrow_grid / 4), shift_col = floor(ncol_grid / 4)),
  list(label = "Shift rows + cols\n(half)",   shift_row = half_row, shift_col = half_col)
)

transform_cells <- do.call(rbind, lapply(transform_cases, function(tc) {
  w <- toroidal_shift_weights(
    grid_weights, nrow_grid, ncol_grid, grid_ids_ordered, occupied_ids,
    shift_row = tc$shift_row,
    shift_col = tc$shift_col
  )
  res <- simulate_with_cell_info(sim_tbl, freq_over, grid_cell_size_km2, w)
  res$transform <- tc$label
  res
}))

transform_cells$transform <- factor(transform_cells$transform,
                                    levels = sapply(transform_cases, `[[`, "label"))

transform_sf <- transform_cells |>
  left_join(grid_sf |> select(grid_id, geometry), by = "grid_id") |>
  st_as_sf() |>
  filter(defor_intensity > 0.001)

p_transforms <- ggplot() +
  geom_sf(data = grid_occupied, fill = "grey85", color = "grey70", linewidth = 0.05) +
  geom_sf(data = transform_sf, aes(fill = defor_intensity), color = NA) +
  scale_fill_gradient(
    low = "#FEE08B", high = "#D73027",
    name = "Deforestation\nintensity",
    limits = c(0, 1)
  ) +
  facet_wrap(~ transform, ncol = 3) +
  theme_minimal(base_size = 10) +
  theme(
    axis.text  = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank()
  ) +
  labs(
    title    = "Shifted model — deterministic transform diagnostic",
    subtitle = "Each panel isolates one transformation; shift = half grid dimensions"
  )

ggsave(paste0(outdir, "/diagnostic_shifted_transforms.pdf"), p_transforms, width = 15, height = 15)
