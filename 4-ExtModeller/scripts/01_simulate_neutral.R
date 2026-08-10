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
  
  #load dependencies
  #library(data.table)
  library(sf)
  library(terra)
  #library(utils)
  #library(lwgeom)
  library(ggplot2)
  library(terra)
  #library(readr)
  library(tidyr)
  #library(stringr)
  #library(parzer)
  library(dplyr)
  library(scales)    # for nicer color scaling
  
  library(parallel)     # base R parallel backend
  library(doParallel)   # foreach-compatible parallel backend
  library(foreach)      # parallel looping
  library(stats)
  
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

extModelDir<-paste0(dirname(getwd()),"/3-ExtPredictor/output")

list.files(paste0(dirname(getwd()),"/3-ExtPredictor/output"))


############################
# IMPORT MODEL PARAMETER CSVs
############################

# list all CSV files that encode model parameters explicitly
model_files <- list.files(
  path       = extModelDir,                          # directory to search
  pattern    = "model_parameters_.*\\.csv$",          # only model parameter CSVs
  full.names = TRUE                                   # return absolute paths
)

# read all model parameter tables into a named list
# lapply keeps this fast and memory-light
model_tables <- lapply(
  model_files,                                        # vector of CSV paths
  read.csv                                            # base R CSV reader
)

# name each list element using the filename (without extension)
names(model_tables) <- tools::file_path_sans_ext(
  basename(model_files)
)

# quick structural check of imported models
str(model_tables)


############################
# IMPORT PATCH AREA FREQUENCY TABLES
############################

# identify frequency-bin CSV files
freq_files <- list.files(
  path       = extModelDir,
  pattern    = "patch_area_frequency_logbins_.*\\.csv$",
  full.names = TRUE
)

# read frequency tables
freq_tables <- lapply(
  freq_files,
  read.csv
)

# name by filtering regime (overfiltered / underfiltered)
names(freq_tables) <- tools::file_path_sans_ext(
  basename(freq_files)
)

# verify bin structure
lapply(freq_tables, head)



############################
# ORGANIZE OBJECTS FOR DOWNSTREAM USE
############################

# extract commonly used objects explicitly for clarity
lognormal_over  <- model_tables$lognormal_model_parameters_overfiltered
lognormal_under <- model_tables$lognormal_model_parameters_underfiltered

powerlaw_over   <- model_tables$powerlaw_model_parameters_overfiltered
powerlaw_under  <- model_tables$powerlaw_model_parameters_underfiltered

freq_over       <- freq_tables$patch_area_frequency_logbins_overfiltered
freq_under      <- freq_tables$patch_area_frequency_logbins_underfiltered


############################
# BASIC SANITY CHECKS
############################

# confirm lognormal parameters are sane
lognormal_over
lognormal_under


############################
# DEFINE A COMMON X-AXIS
############################

# define a common patch-area range (km²) for comparison
# using log-spaced values ensures smooth curves on log axes
x_vals <- 10^seq(
  from = floor(log10(min(c(freq_over$bin_lower, freq_under$bin_lower)))),  # smallest area
  to   = ceiling(log10(max(c(freq_over$bin_upper, freq_under$bin_upper)))), # largest area
  length.out = 500                                                          # smooth curves
)

############################
# EXTRACT LOGNORMAL PARAMETERS
############################

# pull parameters for the overfiltered model
meanlog_over <- lognormal_over$meanlog        # mean of log(area)
sdlog_over   <- lognormal_over$sdlog          # sd of log(area)
xmin_over    <- lognormal_over$xmin_km2       # lower cutoff

# pull parameters for the underfiltered model
meanlog_under <- lognormal_under$meanlog
sdlog_under   <- lognormal_under$sdlog
xmin_under    <- lognormal_under$xmin_km2

############################
# COMPUTE LOGNORMAL PDF VALUES
############################

# compute probability density for overfiltered model
pdf_over <- dlnorm(
  x_vals,                                     # patch area values
  meanlog = meanlog_over,
  sdlog   = sdlog_over
)

# compute probability density for underfiltered model
pdf_under <- dlnorm(
  x_vals,
  meanlog = meanlog_under,
  sdlog   = sdlog_under
)

# set densities to NA below xmin to reflect fitting domain
pdf_over[x_vals < xmin_over]   <- NA
pdf_under[x_vals < xmin_under] <- NA

############################
# PLOT: MODEL VS MODEL (PDF)
############################

# open a clean plotting window
plot(
  x_vals, pdf_over,
  log  = "xy",                                 # log–log axes
  type = "l",                                  # line plot
  lwd  = 3,
  col  = "blue",
  xlab = "Patch area (km²)",
  ylab = "Probability density",
  main = "Lognormal patch-size models: over- vs under-filtered"
)

# add underfiltered model
lines(
  x_vals, pdf_under,
  lwd = 3,
  col = "red"
)

# add vertical lines showing xmin cutoffs
abline(
  v   = xmin_over,
  col = "black",
  lty = 2
)

abline(
  v   = xmin_under,
  col = "black",
  lty = 2
)

# add legend
legend(
  "topright",
  legend = c(
    "Overfiltered model",
    "Underfiltered model",
    "Overfiltered xmin",
    "Underfiltered xmin"
  ),
  col = c("blue", "red", "black", "black"),
  lwd = c(3, 3, 1, 1),
  lty = c(1, 1, 2, 2),
  bty = "n"
)

###############Import animal data#######################

if (exists(".amani_animal_model_dir")) {
  animalModelDir <- .amani_animal_model_dir
} else {
  animalModelDir <- paste0(dirname(getwd()), "/2-ModellingAnimals/output/madre de dios analysis results")
}
sim_files <- list.files(
  path       = animalModelDir,                         # directory containing simulation outputs
  pattern    = "^simulatedOccurencesForAnimal.*\\.rds$",# regex: files starting with name, ending in .csv
  full.names = TRUE                                   # return full paths (required for read.csv)
)

simAnimalObj <- readRDS(sim_files[1])                  # path to .rds file

geom_types <- st_geometry_type(
  simAnimalObj,              # sf object with geometry column
  by_geometry = TRUE         # one geometry type per row
)
# quick overview: how many of each geometry type?
table(geom_types)

# confirm the coordinate reference system (this is critical for grid alignment)
st_crs(simAnimalObj)                              # prints CRS; note EPSG and units

###import/make the grid

dataDir<-paste0(dirname(getwd()),"/1-FilteringPlantData/output")

occ_sf<-readRDS(paste0(dataDir,"/species_counts4_filtered.rds"))

# Define grid size in meters (22.36 km ≈ 22360 meters)
grid_size_m <- 22360

grid <- st_make_grid(occ_sf, cellsize = grid_size_m, what = "polygons")

# Convert grid to an sf object
grid_sf <- st_sf(grid_id = 1:length(grid), geometry = grid)


st_crs(grid_sf)                                   # CRS of the grid
all.equal(st_crs(simAnimalObj), st_crs(grid_sf))  # sanity check CRS match

grid_areas_m2 <- st_area(grid_sf)                 # area of each grid cell in m²
summary(as.numeric(grid_areas_m2))                 # should be nearly constant 

# convert grid cell area to km² for reporting
grid_areas_km2 <- as.numeric(grid_areas_m2) / 1e6 # m² → km²
summary(grid_areas_km2)                            # expected grid size in km²
simAnimal_sf<-simAnimalObj

# ensure CRS match so spatial join is valid
if (st_crs(simAnimal_sf) != st_crs(grid_sf)) {
  simAnimal_sf <- st_transform(simAnimal_sf, st_crs(grid_sf))  # reproject points only
}


# remove rows with missing species names
# NA scientificName cannot contribute to species richness
simAnimal_tbl <- simAnimal_sf |> 
  filter(!is.na(scientificName))            # keep only valid species names

# count number of unique species per grid cell
# n_distinct() ensures multiple occurrences of the same species in a grid
# are only counted once
species_richness_by_grid <- simAnimal_tbl |> 
  group_by(grid_id) |>                      # group by grid cell ID
  summarise(
    species_richness = n_distinct(scientificName),  # unique species count per grid
    .groups = "drop"                        # return an ungrouped data.frame
  )


head(species_richness_by_grid)

#replace NA richness with zero (grids with no animals)
species_richness_by_grid$species_richness[
  is.na(species_richness_by_grid$species_richness)
] <- 0

# join species richness back onto the GRID POLYGONS
# this restores the correct geometry for mapping
grid_sf_with_richness <- grid_sf |>
  left_join(
    species_richness_by_grid |> 
      st_drop_geometry(),          # drop point geometry before join
    by = "grid_id"
  )


p_richness <- ggplot() +
  geom_sf(
    data = grid_sf_with_richness,            # grid with richness already joined
    aes(fill = species_richness),             # color by richness
    color = NA                                # remove polygon borders for cleaner look
  ) +
  scale_fill_viridis_c(
    option = "magma",                         # perceptually uniform color scale
    direction = -1,                           # darker = higher richness
    na.value = "white",                       # grids with no data
    name = "Species richness"                 # legend title
  ) +
  theme_minimal() +                           # clean background
  labs(
    title = "Animal Species Richness per Grid Cell",
    subtitle = "Unique scientificName values per grid_id",
    caption = "Simulated animal occurrences"
  )

# print the plot
print(p_richness)


#OK, now the next step is to build simulated pools of deforested amazon that match our min amd max models

########deforested Amazon: high simulation########
#this gives the probability of drawing patches of various sizes
patchModel<-lognormal_under

#the total empirical area deforested when this model was estimated
totalDeforestedArea<-lognormal_under$total_area_km2

# --- Total deforested area ---
# Computed via terra::cellSize on def_00_20.tif (latitude-weighted cell areas).
# The binary raster contains 514,868,669 deforested cells (all non-NA cells).
# Source: 3-ExtPredictor/scripts/area_job/output/area_accounting.csv

total_deforested_km2 <- 454509.4

# totalDeforestedArea (from lognormal model fit) only counts patches >= 50 km2,
# so it represents a subset of total_deforested_km2. The gap is mostly small
# patches below the 50 km2 extraction threshold, plus minor truncation of
# large patches at window edges.

patchModel

# compute area of each grid cell in square meters
grid_sf$grid_area_m2 <- as.numeric(
  st_area(grid_sf)           # polygon area in m²
)

# convert grid cell area to km²
grid_sf$grid_area_km2 <- grid_sf$grid_area_m2 / 1e6
grid_cell_size_km2 <- round(median(grid_sf$grid_area_km2))
stopifnot(length(grid_cell_size_km2) == 1)

#use the frequency table calculated before to draw random patchs

freq_over[22,]
freq_over[5,]
simAnimalObj

###below is a single, full simulation
out_all <- simulateExtinctionFromDeforestation(
     sim_tbl            = simAnimalObj,
     freq_table         = freq_over,
     grid_cell_size_km2 = grid_cell_size_km2
   )

   out_all$n_pop_lost_total
   out_all$n_spp_lost_total
   out_all$final_table

   
########FULL RUN WITH MULTIPLE REPLICATES#######
#######FREQ OVER#####

results_file_over <- list.files(
  path = outdir, pattern = "^neutral_FreqOver_.*\\.RData$", full.names = TRUE
)

if (length(results_file_over) > 0) {
  cat("Found existing FreqOver results:\n")
  cat(paste(" ", results_file_over, collapse = "\n"), "\n")
  cat("Loading most recent file.\n")
  cat("To re-run, delete the file and re-source this section.\n")
  load(results_file_over[length(results_file_over)])
} else {
  extinction_results_FreqOver <- SEFD_rep_parallel(
    simAnimalObj,
    freq_table = freq_over,
    grid_cell_size_km2,
    n_runs = 100
  )
  save(
    extinction_results_FreqOver,
    file = paste0(
      outdir, "/neutral_FreqOver_",
      format(Sys.time(), "%Y%m%d_%H%M%S"),
      ".RData"
    )
  )
}

   # quick diagnostics
   str(extinction_results_FreqOver)
   summary(extinction_results_FreqOver)
   
   # OPTIONAL quick diagnostics
    hist(extinction_results_FreqOver$n_spp_lost_total, breaks = 30,
         main = "Species extinctions across parallel simulations",
         xlab = "Number of species extinct")
   
    hist(extinction_results_FreqOver$n_pop_lost_total, breaks = 30,
         main = "Population extinctions across parallel simulations",
         xlab = "Number of populations extinct")
   
############################
# SANITY CHECK: back-of-envelope expected population losses
############################
# Compare simulation output against a simple linear expectation.
# For each bin, estimate how many grid-cell equivalents are deforested,
# then multiply by the mean number of population records per grid cell.
# This ignores cumulative depletion, so it overestimates slightly.

    popsizes <- table(simAnimalObj$grid_id)
    mean_pops_per_grid <- mean(popsizes)
    total_pops <- nrow(sf::st_drop_geometry(simAnimalObj))

    # large patches (>= grid cell): remove whole grid cells
    freq_large <- freq_over[
      freq_over$bin_mid >= grid_cell_size_km2 & freq_over$count > 0, ,
      drop = FALSE
    ]

    expected_lost_large <- sum(
      (freq_large$bin_mid / grid_cell_size_km2) *
        freq_large$count *
        mean_pops_per_grid
    )

    # small patches (< grid cell): fractional removal per grid cell
    freq_small <- freq_over[
      freq_over$bin_mid < grid_cell_size_km2 & freq_over$count > 0, ,
      drop = FALSE
    ]

    expected_lost_small <- sum(
      (freq_small$bin_mid / grid_cell_size_km2) *
        freq_small$count *
        mean_pops_per_grid
    )

    # total area deforested by each regime
    area_large_km2 <- sum(freq_large$bin_mid * freq_large$count)
    area_small_km2 <- sum(freq_small$bin_mid * freq_small$count)

    cat("--- Back-of-envelope sanity check ---\n")
    cat("Mean pops per grid cell:", round(mean_pops_per_grid, 2), "\n")
    cat("Expected pops lost (large patches):", round(expected_lost_large), "\n")
    cat("Expected pops lost (small patches):", round(expected_lost_small), "\n")
    cat("Expected pops lost (total):", round(expected_lost_large + expected_lost_small), "\n")
    cat("As % of total pops:", round(100 * (expected_lost_large + expected_lost_small) / total_pops, 1), "%\n")
    cat("Total area large patches:", round(area_large_km2), "km2\n")
    cat("Total area small patches:", round(area_small_km2), "km2\n")
    cat("Ratio large/small area:", round(area_large_km2 / area_small_km2, 2), "\n")
#######FREQ UNDER#####

results_file_under <- list.files(
  path = outdir, pattern = "^neutral_FreqUnder_.*\\.RData$", full.names = TRUE
)

if (length(results_file_under) > 0) {
  cat("Found existing FreqUnder results:\n")
  cat(paste(" ", results_file_under, collapse = "\n"), "\n")
  cat("Loading most recent file.\n")
  cat("To re-run, delete the file and re-source this section.\n")
  load(results_file_under[length(results_file_under)])
} else {
  extinction_results_FreqUnder <- SEFD_rep_parallel(
    simAnimalObj,
    freq_table = freq_under,
    grid_cell_size_km2,
    n_runs = 100
  )
  save(
    extinction_results_FreqUnder,
    file = paste0(
      outdir, "/neutral_FreqUnder_",
      format(Sys.time(), "%Y%m%d_%H%M%S"),
      ".RData"
    )
  )
}


#######FREQ BOUNDARY (middle scenario)#####

freq_boundary <- freq_tables$patch_area_frequency_logbins_boundary_filtered

results_file_boundary <- list.files(
  path = outdir, pattern = "^neutral_FreqBoundary_.*\\.RData$", full.names = TRUE
)

if (length(results_file_boundary) > 0) {
  cat("Found existing FreqBoundary results:\n")
  cat(paste(" ", results_file_boundary, collapse = "\n"), "\n")
  cat("Loading most recent file.\n")
  cat("To re-run, delete the file and re-source this section.\n")
  load(results_file_boundary[length(results_file_boundary)])
} else {
  extinction_results_FreqBoundary <- SEFD_rep_parallel(
    simAnimalObj,
    freq_table = freq_boundary,
    grid_cell_size_km2,
    n_runs = 100
  )
  save(
    extinction_results_FreqBoundary,
    file = paste0(
      outdir, "/neutral_FreqBoundary_",
      format(Sys.time(), "%Y%m%d_%H%M%S"),
      ".RData"
    )
  )
}


############################
# COMPARE ALL THREE SCENARIOS
############################

n_spp_total <- length(unique(simAnimalObj$scientificName))

results_compare <- data.frame(
  scenario = c("FreqOver (conservative)", "FreqBoundary (middle)", "FreqUnder (liberal)"),
  total_patch_area_km2 = c(
    sum(freq_over$bin_mid * freq_over$count),
    sum(freq_boundary$bin_mid * freq_boundary$count),
    sum(freq_under$bin_mid * freq_under$count)
  ),
  n_patches = c(sum(freq_over$count), sum(freq_boundary$count), sum(freq_under$count)),
  median_spp_lost = c(
    median(extinction_results_FreqOver$n_spp_lost_total),
    median(extinction_results_FreqBoundary$n_spp_lost_total),
    median(extinction_results_FreqUnder$n_spp_lost_total)
  ),
  sd_spp_lost = round(c(
    sd(extinction_results_FreqOver$n_spp_lost_total),
    sd(extinction_results_FreqBoundary$n_spp_lost_total),
    sd(extinction_results_FreqUnder$n_spp_lost_total)
  ), 1),
  median_pop_lost = c(
    median(extinction_results_FreqOver$n_pop_lost_total),
    median(extinction_results_FreqBoundary$n_pop_lost_total),
    median(extinction_results_FreqUnder$n_pop_lost_total)
  ),
  pct_spp_lost = round(100 * c(
    median(extinction_results_FreqOver$n_spp_lost_total),
    median(extinction_results_FreqBoundary$n_spp_lost_total),
    median(extinction_results_FreqUnder$n_spp_lost_total)
  ) / n_spp_total, 1)
)

results_compare


############################
# INVESTIGATE EXTINCT SPECIES: RANGE-RESTRICTED VS WIDESPREAD
############################
# Run one detailed simulation per scenario to identify which species go
# extinct, then compare their range sizes to surviving species.

sim_tbl <- sf::st_drop_geometry(simAnimalObj)

set.seed(4817)
out_over <- simulateExtinctionFromDeforestation(
  sim_tbl = sim_tbl, freq_table = freq_over, grid_cell_size_km2 = grid_cell_size_km2
)
out_boundary <- simulateExtinctionFromDeforestation(
  sim_tbl = sim_tbl, freq_table = freq_boundary, grid_cell_size_km2 = grid_cell_size_km2
)
out_under <- simulateExtinctionFromDeforestation(
  sim_tbl = sim_tbl, freq_table = freq_under, grid_cell_size_km2 = grid_cell_size_km2
)

# range size = number of distinct grid cells a species occupies
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

# summary: range sizes of extinct vs surviving species
cat("=== Range size (n_grids) of extinct vs surviving species ===\n")
for (scen in c("over", "boundary", "under")) {
  col <- paste0("extinct_", scen)
  cat("\n--- Freq", scen, "---\n")
  cat("Extinct:\n")
  print(summary(range_sizes$n_grids[range_sizes[[col]]]))
  cat("Surviving:\n")
  print(summary(range_sizes$n_grids[!range_sizes[[col]]]))
}

# single-grid endemic analysis
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

cat("\n% of extinct spp that are single-grid:\n")
cat("  FreqOver:    ",
    round(100 * mean(range_sizes$n_grids[range_sizes$extinct_over] == 1), 1), "%\n")
cat("  FreqBoundary:",
    round(100 * mean(range_sizes$n_grids[range_sizes$extinct_boundary] == 1), 1), "%\n")
cat("  FreqUnder:   ",
    round(100 * mean(range_sizes$n_grids[range_sizes$extinct_under] == 1), 1), "%\n")

# multi-grid extinctions
cat("\nMulti-grid species that went extinct:\n")
cat("  FreqOver:    ", sum(range_sizes$extinct_over & range_sizes$n_grids > 1), "\n")
cat("  FreqBoundary:", sum(range_sizes$extinct_boundary & range_sizes$n_grids > 1), "\n")
cat("  FreqUnder:   ", sum(range_sizes$extinct_under & range_sizes$n_grids > 1), "\n")

multi_extinct_under <- range_sizes |>
  filter(extinct_under & n_grids > 1)
if (nrow(multi_extinct_under) > 0) {
  cat("  FreqUnder multi-grid range sizes:", sort(multi_extinct_under$n_grids), "\n")
}


############################
# PER-SPECIES RANGE REDUCTION (single detailed runs)
############################
# For each scenario, compute what % of each species' range was lost.
# Includes species that survived (partial reduction) and those that
# went fully extinct (100% reduction).

baseline_ranges <- precompute_baseline_ranges(sim_tbl)

rr_over     <- compute_range_reduction(out_over$final_table,     baseline_ranges)
rr_boundary <- compute_range_reduction(out_boundary$final_table, baseline_ranges)
rr_under    <- compute_range_reduction(out_under$final_table,    baseline_ranges)

rr_over$scenario     <- "FreqOver"
rr_boundary$scenario <- "FreqBoundary"
rr_under$scenario    <- "FreqUnder"

neutral_range_reduction <- rbind(rr_over, rr_boundary, rr_under)

write.csv(neutral_range_reduction,
          paste0(outdir, "/neutral_species_range_reduction.csv"),
          row.names = FALSE)
cat("Saved output/neutral_species_range_reduction.csv\n")

# Summary: number of species at each reduction threshold
cat("\n=== Range Reduction Thresholds (neutral, single run) ===\n")
for (scen in c("FreqOver", "FreqBoundary", "FreqUnder")) {
  rr_sub <- neutral_range_reduction[neutral_range_reduction$scenario == scen, ]
  tc <- count_range_thresholds(rr_sub)
  cat("\n", scen, ":\n")
  for (j in seq_along(tc)) {
    cat("  ", names(tc)[j], ": ", tc[j], "\n")
  }
}


############################
# EXPORT RESULTS TO CSV
############################

# 1. Scenario comparison summary
write.csv(results_compare, paste0(outdir, "/neutral_scenario_comparison.csv"), row.names = FALSE)
cat("Saved output/neutral_scenario_comparison.csv\n")

# 2. Raw replicate results (all 100 runs per scenario)
extinction_results_FreqOver$scenario <- "FreqOver"
extinction_results_FreqBoundary$scenario <- "FreqBoundary"
extinction_results_FreqUnder$scenario <- "FreqUnder"

neutral_all_replicates <- rbind(
  extinction_results_FreqOver,
  extinction_results_FreqBoundary,
  extinction_results_FreqUnder
)
write.csv(neutral_all_replicates, paste0(outdir, "/neutral_all_replicates.csv"), row.names = FALSE)
cat("Saved output/neutral_all_replicates.csv\n")

# 3. Species extinction status (from single detailed run)
write.csv(
  range_sizes,
  paste0(outdir, "/neutral_species_extinction_status.csv"),
  row.names = FALSE
)
cat("Saved output/neutral_species_extinction_status.csv\n")


############################
# DONE / STATUS
############################
# [DONE] Incorporate spatially clustered deforestation placement (see 02_simulate_clustered.R)
# [DONE] Compare neutral vs clustered vs shifted (toroidal null model)
# [DONE] Investigate which species go extinct (range-restricted vs widespread)
#
# See NOTES.txt in the project root for outstanding items and next steps.
