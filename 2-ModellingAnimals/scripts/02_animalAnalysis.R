print("Starting spatial extent analysis....")

# Working directory should already be set to project root (2-ModellingAnimals/)
# by the calling script (01_CockroachAnalysis.R, 01_GBIFAnalysis.R, etc.)

SEEDnumber<-2024299

###START: LOAD LIBRARIES####
library(progress)
library(data.table)
library(sf)
#library(lwgeom)
library(ggplot2)
library(terra)
library(tidyr)
#library(readr)
library(stringr)
#library(parzer)
#library(MASS)
library(dplyr)
#library(readxl)
#library(vegan)
#library(scales)
#library(glue)
#library(purrr)
library(minpack.lm)


###Output results folder####
if (exists(".amani_output_dir")) {
  output_folder <- .amani_output_dir
} else {
  output_folder <- "output/madre de dios analysis results"
}

# Ensure plot_data/ subdirectory exists for saving plot data as RDS
plot_data_dir <- file.path(output_folder, "plot_data")
dir.create(plot_data_dir, showWarnings = FALSE, recursive = TRUE)

#### Write dataset summary report (pre- and post-filter) ####
# This captures key stats about the data before and after cell filtering,
# which is essential for evaluating whether downstream simulation estimates
# are realistic.

if (exists(".amani_filter_summary", envir = .GlobalEnv)) {
  fs <- .amani_filter_summary
  pre <- fs$pre_filter
  post <- fs$post_filter
  thr <- fs$thresholds

  fmt_summary <- function(s) {
    paste0("min=", round(s["Min."], 1), ", Q1=", round(s["1st Qu."], 1),
           ", med=", round(s["Median"], 1), ", mean=", round(s["Mean"], 1),
           ", Q3=", round(s["3rd Qu."], 1), ", max=", round(s["Max."], 1))
  }

  dataset_report <- paste0(
    "=======================================================\n",
    " Dataset Summary: Before and After Cell Filtering\n",
    " Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n",
    "=======================================================\n\n",
    "--- Filter thresholds ---\n",
    "  min_occurrences:   ", thr$min_occurrences, "\n",
    "  min_completeness:  ", round(thr$min_completeness, 3), "\n",
    "  min_richness:      ", thr$min_richness, "\n",
    "  outlier cutoff:    ", thr$outlier_source, "\n",
    "  Track B entry:     97th percentile (", round(thr$upper_fence_97), " spp)\n",
    "  Filter applied:    ", ifelse(fs$filter_applied, "YES", "NO (< 30 cells)"), "\n\n",
    "--- Before filtering ---\n",
    "  Grid cells:        ", pre$n_cells, "\n",
    "  Total records:     ", pre$n_records, "\n",
    "  Total species:     ", pre$n_species, "\n",
    "  Richness/cell:     ", fmt_summary(pre$richness_summary), "\n",
    "  Occurrences/cell:  ", fmt_summary(pre$occ_summary), "\n",
    "  Completeness/cell: ", fmt_summary(pre$comp_summary), "\n\n",
    "--- After filtering ---\n",
    "  Grid cells:        ", post$n_cells,
    " (", post$n_well_sampled, " well-sampled + ",
    post$n_super_diverse, " super-diverse; ",
    post$n_outliers, " outliers dropped)\n",
    "  Total records:     ", post$n_records,
    " (", round(post$n_records / pre$n_records * 100, 1), "% of original)\n",
    "  Total species:     ", post$n_species,
    " (", round(post$n_species / pre$n_species * 100, 1), "% of original)\n",
    "  Richness/cell:     ", fmt_summary(post$richness_summary), "\n",
    "  Occurrences/cell:  ", fmt_summary(post$occ_summary), "\n",
    "  Completeness/cell: ", fmt_summary(post$comp_summary), "\n\n",
    "--- Species loss from filtering ---\n",
    "  Species retained:  ", post$n_species, " / ", pre$n_species, "\n",
    "  Species lost:      ", pre$n_species - post$n_species,
    " (", round((1 - post$n_species / pre$n_species) * 100, 1), "%)\n\n",
    "=======================================================\n"
  )

  # Print to console
  cat(dataset_report)

  # Write to file
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  ds_report_path <- file.path(output_folder, paste0("dataset_summary_", timestamp, ".txt"))
  writeLines(dataset_report, con = ds_report_path)
  message("Dataset summary written to: ", ds_report_path)

} else if (exists("animal_grid", envir = .GlobalEnv)) {
  # Filtering was disabled (e.g. min_occurrences = 0) — write a minimal summary
  # with the basic counts so every run produces a dataset_summary_*.txt.
  ag_tbl <- sf::st_drop_geometry(animal_grid) |>
    dplyr::filter(!is.na(grid_id))

  n_cells_raw   <- dplyr::n_distinct(ag_tbl$grid_id)
  n_records_raw <- nrow(ag_tbl)
  n_species_raw <- dplyr::n_distinct(ag_tbl$Morphotype)

  dataset_report <- paste0(
    "=======================================================\n",
    " Dataset Summary (no cell filtering applied)\n",
    " Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n",
    "=======================================================\n\n",
    "  Filter applied:   NO (min_occurrences = 0 or filtering disabled)\n\n",
    "--- Dataset counts ---\n",
    "  Grid cells:       ", n_cells_raw, "\n",
    "  Total records:    ", n_records_raw, "\n",
    "  Total species:    ", n_species_raw, "\n\n",
    "=======================================================\n"
  )

  cat(dataset_report)

  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  ds_report_path <- file.path(output_folder, paste0("dataset_summary_", timestamp, ".txt"))
  writeLines(dataset_report, con = ds_report_path)
  message("Dataset summary written to: ", ds_report_path)
}

################Plant import and setup###################
###step 1 make the grid

dataDir<-paste0(dirname(getwd()),"/1-FilteringPlantData/output")
#dataDir<-paste0(dirname(dirname(getwd())),"/1-FilteringPlantData/output")


species_counts<-occ_sf#readRDS(paste0(dataDir,"/species_counts4_filtered.rds"))
#species_counts<-readRDS(paste0(dataDir,"/species_counts4_filtered.rds"))
#endemismData<-readRDS(paste0(dataDir,"/FINAL_allEndemism_FullyFiltered.rds"))
 rm(occ_sf)
# Extract the grid from the endemism data to preserve Stage 0 grid_ids.
# Do NOT rebuild via st_make_grid() — different bounding box produces
# different grid_ids that break downstream joins with species_counts.
grid_size_m <- 22360
grid_sf <- endemismData


timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
output_filename <- paste0("rawEndemismGridplot.", timestamp, ".pdf")
output_path <- file.path(output_folder, output_filename)

# Create the PDF plot
pdf(output_path, width = 8, height = 6)
p<-ggplot() +
  geom_sf(data = grid_sf, aes(fill = raw_endemism), color = NA) +
  scale_fill_viridis_c(option = "plasma", na.value = "white") +
  theme_minimal() +
  labs(title = "Species-Level Endemism (500 sq km grid)", fill = "Raw Endemism")
print(p)
dev.off()

timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
output_filename <- paste0("weightedEndemismGridplot.", timestamp, ".pdf")
output_path <- file.path(output_folder, output_filename)
pdf(output_path, width = 8, height = 6)

p<-ggplot() +
  geom_sf(data = grid_sf, aes(fill = weighted_endemism), color = NA) +
  scale_fill_viridis_c(option = "plasma", na.value = "white") +
  theme_minimal() +
  labs(title = "Species-Level Endemism (500 sq km grid)", fill = "Weighted Endemism")
print(p)
dev.off()


#########simulate gamma richness based on distribution of all endemism (not just single cell endemics)#########

#####Now, get an estimate of the distribution of all cell-presences for each species
#In other words, simulating endemic species only is simulating species that are present in only 1 cell.
#If I can now take presence-absence data for each species, scale it based on the animal data, I can simulate similar datasets and just count the number of total species

#First, take the animal data and plot the distribution of presences per cell.

###basic abundance histogram (commented out — crashes RStudio when Plots pane is small)
# animal_grid %>%
#   count(Morphotype) %>%
#   ggplot(aes(x = reorder(Morphotype, n), y = n)) +
#   geom_col(fill = "steelblue") +
#   coord_flip() +
#   labs(
#     title = "Morphotype Frequency in animal_grid",
#     x = "Morphotype",
#     y = "Number of Occurrences"
#   ) +
#   theme_minimal()

# Now plot the distribution of morphotypes and number of grids cells they occur in
# Step 1: Count unique grid cells per morphotype
#file_path<-"C:\\Users\\dominice\\OneDrive - University of Illinois - Urbana\\Projects\\Peruvian Blattodea\\Cockroaches in a coal mine\\AmAni\\2-ModellingAnimals\\data\\DarwinCoreAnimalOccurences.8-1-2026.xlsx"
#animal_grid<-

animal_grid<-rename_grid_id(animal_grid) #this function is defined in rawEndemismCalculator.R

morphotype_grid_counts <- animal_grid %>%
  distinct(Morphotype, grid_id) %>%
  count(Morphotype, name = "num_grid_cells")

# Step 2: Count how many morphotypes have each grid cell count
morphotype_distribution <- morphotype_grid_counts %>%
  count(num_grid_cells)

# Step 3: Plot
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
output_filename <- paste0("dist.morphotype.spat.Spread.", timestamp, ".pdf")
output_path <- file.path(output_folder, output_filename)

print("Printing histogram of animal spatial extent frequencies...")

# Create the PDF plot
pdf(output_path, width = 8, height = 6)


p<-ggplot(morphotype_distribution, aes(x = num_grid_cells, y = n)) +
  geom_col(fill = "lightblue") +
  labs(
    title = "Distribution of Morphotypes by Grid Cell Count",
    x = "Number of Grid Cells",
    y = "Number of Morphotypes"
  ) +
  theme_minimal()
print(p)
dev.off()

###now do the same with the plant data:
###1: all the plant data
#class(species_counts)

setDF(species_counts)
species_counts <- st_as_sf(
  species_counts,
  sf_column_name = "geometry",
  crs = 3857
)


# Step 1: Count unique grid cells per plant species

#plant_grid_counts <- species_counts %>%
#  distinct(scientificName, grid_id) %>%
#  count(scientificName, name = "num_grid_cells")
length(unique(species_counts$scientificName))

plant_grid_counts <- species_counts %>%
  sf::st_drop_geometry() %>%
  distinct(scientificName, grid_id) %>%
  count(scientificName, name = "num_grid_cells")
nrow(plant_grid_counts)

species_counts |>
  sf::st_drop_geometry() |>
  dplyr::distinct(scientificName) |>
  nrow()
data.table::uniqueN(species_counts$scientificName)


# Step 2: Count how many times each grid cell count occurs in the list
plant_distribution <- plant_grid_counts %>%
  count(num_grid_cells)
head(plant_distribution)

sum(plant_distribution$n)

plant_grid_counts

# Step 3: Plot
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
output_filename <- paste0("hist_plantSpp_spatExt.", timestamp, ".pdf")
output_path <- file.path(output_folder, output_filename)

# Open PDF device
pdf(output_path, width = 8, height = 6)

p<-ggplot(plant_distribution, aes(x = num_grid_cells, y = n)) +
  geom_col(fill = "lightblue") +
  labs(
    title = "Distribution of Plant Species by Grid Cell Count",
    x = "Number of Grid Cells",
    y = "Number of Plant Species"
  ) +
  theme_minimal()
print(p)
dev.off()

###2: just the madre de dios the plant data

compute_plant_richness <- function(
    species_counts,
    redoeverything = FALSE,  # logical flag: TRUE = force recompute, FALSE = reuse if available
    fill_plants = FALSE      # which plant-data type this run uses (keys the cache)
) {
  # Cache file path (saved alongside the plant data so it's shared across organisms)
  cache_file <- paste0(dirname(getwd()), "/1-FilteringPlantData/output/plant_richness_cache.rds")
  cache_file_Filled <- paste0(dirname(getwd()), "/1-FilteringPlantData/output/filled_plant_rich_cache.rds")
  # Richness differs between filled and unfilled plant data, so key the cache by type.
  cache_file <- if (isTRUE(fill_plants)) cache_file_Filled else cache_file
  
  # 1. Check in-memory cache (fastest)
  if (exists("plant_richness", envir = .GlobalEnv) && !redoeverything &&
      identical(get0(".plant_richness_fill_state", envir = .GlobalEnv, ifnotfound = NA), isTRUE(fill_plants))) {
    message("plant_richness already exists in memory (matching fill_plants=", isTRUE(fill_plants), ") — reusing.")
    return(get("plant_richness", envir = .GlobalEnv))
  }

  # 2. Check on-disk cache
  if (file.exists(cache_file) && !redoeverything) {
    message("Loading plant_richness from cache: ", cache_file)
    plant_richness <- readRDS(cache_file)
    assign(".plant_richness_fill_state", isTRUE(fill_plants), envir = .GlobalEnv)
    return(plant_richness)
  }

  # 3. Compute from scratch
  if (redoeverything) {
    message("plant_richness — recomputing because redoeverything = TRUE.")
  } else {
    message("plant_richness — computing for the first time.")
  }
  message("Summarizing the plant data (fill_plants=", isTRUE(fill_plants), ") for a manipulation (this may take a while)...")

  plant_richness <- species_counts %>%
    dplyr::group_by(grid_id) %>%
    dplyr::summarize(
      plant_richness = dplyr::n_distinct(scientificName),
      .groups = "drop"
    )

  message("Done summarizing. Saving cache to: ", cache_file)
  saveRDS(plant_richness, cache_file)
  assign(".plant_richness_fill_state", isTRUE(fill_plants), envir = .GlobalEnv)

  plant_richness
}

plant_richness <- compute_plant_richness(
  species_counts,
  redoeverything = redoeverything,
  fill_plants = isTRUE(get0(".amani_fill_plants", envir = .GlobalEnv, ifnotfound = FALSE))
)

# Recompute plant single-cell endemism within the overlap region only,
# matching the approach in rawEndemismCalculator.R.
# Use grid_id filtering rather than st_crop to avoid sf/data.table class issues.
overlap_grid_ids_02 <- grid_sf$grid_id[
  as.logical(st_intersects(grid_sf, st_as_sfc(overlap_bbox), sparse = FALSE))
]

local_plant_endem_02 <- species_counts %>%
  sf::st_drop_geometry() %>%
  filter(grid_id %in% overlap_grid_ids_02) %>%
  group_by(scientificName) %>%
  filter(n_distinct(grid_id) == 1) %>%
  ungroup() %>%
  group_by(grid_id) %>%
  summarise(raw_endemism = n_distinct(scientificName), .groups = "drop")

# Select only grid_id before joining to avoid column name collisions
# (grid_sf may already have a raw_endemism column from upstream)
plant_richness_grid_narrow <- grid_sf %>%
  dplyr::select(grid_id, geometry) %>%
  st_crop(overlap_bbox) %>%
  left_join(local_plant_endem_02, by = "grid_id") %>%
  filter(!is.na(raw_endemism), raw_endemism > 0) %>%
  dplyr::select(grid_id, geometry, raw_endemism)

plant_richness_overlap <- st_crop(plant_richness, overlap_bbox)


plant_grid_joined <- st_join(grid_sf, plant_richness_overlap %>% dplyr::select(-grid_id))
plant_grid_joined <- rename_grid_id(plant_grid_joined)
plant_grid_joined <- st_crop(plant_grid_joined, overlap_bbox)

plant_richness_grid <- plant_grid_joined %>%
  group_by(grid_id)%>%
  filter(plant_richness > 0) %>%
  summarise(
    plant_richness = sum(plant_richness, na.rm = TRUE),
    geometry = st_union(geometry)
  ) %>%
  ungroup()


print("Plotting plant richness in area of animal overlap...")
#####plot plant richness in the madre de dios######
#Step 4: Create the plant plot
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
output_filename <- paste0("plant_richness_map_overlap", timestamp, ".pdf")
output_path <- file.path(output_folder, output_filename)

# Open PDF device
pdf(output_path, width = 8, height = 6)

# Plant Richness Plot (Now Properly Cropped)
p<-ggplot(data = plant_richness_grid) +
  geom_sf(aes(fill = plant_richness), color = "black", size = 0.1) +
  scale_fill_viridis_c(option = "magma", direction = -1, name = "Plant Richness") +
  theme_minimal() +
  theme(legend.position = "right") +
  labs(
    title = "Plant Richness in Area of Animal Overlap",
    caption = "Data: Aggregated plant richness per Amazon grid cell"
  )

print(p)

# Close the device
dev.off()

######get the distribution
colnames(grid_sf)[1]<-"grid_id"
#species_counts_joined<-st_join(grid_sf, species_counts, by = "grid_id")

species_counts_joined <- grid_sf %>%
  left_join(
    sf::st_drop_geometry(species_counts),
    by = "grid_id"
  )

#if RawEndemism Calculator was skipped we need to calculate some things here
if(.amani_skipRawEndemism){ overlap_bbox <- st_bbox(animal_endemism_filtered) }

species_counts_overlap <- st_crop(species_counts_joined, overlap_bbox)

filtered_species_counts <- species_counts_overlap %>%
  filter(grid_id %in% animal_grid$grid_id)



plant_grid_counts <- as.data.table(sf::st_drop_geometry(filtered_species_counts))[
  , .(num_grid_cells = uniqueN(grid_id)),
  by = scientificName
]



#slower alternative to above
#plant_grid_counts <- filtered_species_counts %>%
#  distinct(scientificName , grid_id) %>%
#  count(scientificName , name = "num_grid_cells")



# Step 2: Count how many morphotypes have each grid cell count
plantSpp_distribution <- plant_grid_counts %>%
  count(num_grid_cells)

# Guard: lognormal fitting requires enough distinct range-size values to be
# numerically stable. With <=3 unique num_grid_cells values the 3-parameter
# NLS converges to pathological solutions (e.g. meanlog ~50), making dlnorm
# return ~0 and causing the extent scaling factor to overflow.
n_distinct_cells <- nrow(plantSpp_distribution)
if (n_distinct_cells < 4) {
  stop(
    "Insufficient plant range-size variation in the animal overlap area: only ",
    n_distinct_cells, " distinct num_grid_cells value(s) found (minimum 4 required). ",
    "The overlap between the plant occurrence data and the animal grid is too small ",
    "to fit a log-normal distribution reliably. Check that 'species_counts' covers ",
    "the full study extent and that the animal grid has adequate spatial overlap ",
    "with the plant data."
  )
}

# Free large intermediates no longer needed (~1.5 GB)
rm(species_counts_joined, species_counts_overlap, filtered_species_counts)
gc()

#####plot spatial extent of plant richness in the area of overlap######
print("Plotting spatial extent of plant richness in the area of animal overlap...")

# Step 3: Plot
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
output_filename <- paste0("dist.PlantSpp.MDD.spat.Spread.", timestamp, ".pdf")
output_path <- file.path(output_folder, output_filename)

# Create the PDF plot
pdf(output_path, width = 8, height = 6)
 
p<-ggplot(plantSpp_distribution, aes(x = num_grid_cells, y = n)) +
  geom_col(fill = "lightblue") +
  labs(
    title = "Distribution of Plant Species by Grid Cell Count",
    x = "Number of Grid Cells",
    y = "Number of Plant Species"
  ) +
  theme_minimal()
print(p)
dev.off()

#####compare spatial extent of plant richness and animal richness in the area of overlap######

#ok now i want to compare the two distributions

# Add identifier columns
plantSpp_distribution$group <- "Plant species"
morphotype_distribution$group <- "Animal morphotypes"

# Bind into one dataframe
combined_distribution <- bind_rows(plantSpp_distribution, morphotype_distribution)

timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
output_filename <- paste0("combined_SpatExtents.", timestamp, ".pdf")
output_path <- file.path(output_folder, output_filename)

# Open PDF device
pdf(output_path, width = 8, height = 6)
p<-ggplot(combined_distribution, aes(x = factor(num_grid_cells), y = n, fill = group)) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  scale_fill_viridis_d(name = NULL, option = "magma", direction = -1) +
  labs(
    title = "Distribution of Species by Grid Cell Occupancy",
    x = "Number of Grid Cells",
    y = "Number of Species / Morphotypes"
  ) +
  theme_minimal() +
  theme(legend.position = "top")
print(p)
dev.off()

# Save combined spatial extent distributions for cross-organism summaries
saveRDS(morphotype_distribution,
        file.path(plot_data_dir, "plotdata_animal_extent_distribution.rds"))
saveRDS(combined_distribution,
        file.path(plot_data_dir, "plotdata_spatial_extent_distributions.rds"))

#####using plant overlap spatial extents to extrapolate on animal spatial extents########
####now fit a log-normal distribution to the plant data

# Data prep
plant_df <- plantSpp_distribution %>%
  arrange(num_grid_cells)

# Fit log-normal: convert x to numeric (meanlog, sdlog estimated manually or with nls)
x_vals <- plant_df$num_grid_cells
y_vals <- plant_df$n


# Estimate starting values from weighted moments of the distribution
weights_norm <- y_vals / sum(y_vals)
start_meanlog <- sum(weights_norm * log(x_vals + 1))
start_sdlog   <- sqrt(sum(weights_norm * (log(x_vals + 1) - start_meanlog)^2))

# Use nls to fit a log-normal curve: f(x) = a * dlnorm(x, meanlog, sdlog)
# nlsLM (Levenberg-Marquardt) is more robust than base nls for this kind of fit
lognorm_fit <- nlsLM(
  y_vals ~ a * dlnorm(x_vals, meanlog, sdlog),
  start = list(a = max(y_vals), meanlog = start_meanlog, sdlog = start_sdlog)
)


# Predict values from fit
x_pred <- seq(min(x_vals), max(x_vals), by = 1) #this used to be b=0.1 but I see no reason to do that other than for aesthics
fit_pred <- data.frame(
  num_grid_cells = x_pred,
  fitted = predict(lognorm_fit, newdata = data.frame(x_vals = x_pred))
)
colnames(fit_pred)[[2]]<-"n"

###now rescale that model distribution to match the animal data

# Get max value from morphotype data
animal_max <- max(morphotype_distribution$n, na.rm = TRUE)

# ── Shape-shift parameters ──────────────────────────────────────────────
# Resolve which log-normal parameters to estimate from the animal data.
# fit_params is a character vector: subset of c("amplitude", "meanlog", "sdlog").
# - "amplitude" (always on): scale the total species count (existing behaviour).
# - "meanlog": shift the distribution's central tendency (range sizes).
# - "sdlog":   rescale the distribution's spread (endemic-vs-widespread ratio).
# Parameters NOT listed are inherited unmodified from the plant fit.

fit_params <- if (exists(".amani_fit_params", envir = .GlobalEnv)) {
  .amani_fit_params
} else {
  c("amplitude")
}

# Plant log-normal parameters (already fitted above)
plant_mdd_meanlog <- coef(lognorm_fit)[["meanlog"]]
plant_mdd_sdlog   <- coef(lognorm_fit)[["sdlog"]]

# Fit log-normal to animal spatial extent distribution
animal_df_for_fit <- morphotype_distribution |>
  dplyr::filter(n > 0) |>
  dplyr::arrange(num_grid_cells)

lognorm_fit_animal <- tryCatch({
  nls(
    n ~ a * dlnorm(num_grid_cells, meanlog, sdlog),
    data = animal_df_for_fit,
    start = list(a = max(animal_df_for_fit$n),
                 meanlog = plant_mdd_meanlog,
                 sdlog   = plant_mdd_sdlog),
    control = nls.control(maxiter = 200, warnOnly = TRUE)
  )
}, error = function(e) {
  warning("Could not fit log-normal to animal MDD data: ", e$message,
          "\n  Falling back to amplitude-only scaling.", call. = FALSE)
  NULL
})

if (!is.null(lognorm_fit_animal) && !lognorm_fit_animal$convInfo$isConv) {
  warning("Animal log-normal fit did not converge. Falling back to amplitude-only.")
  lognorm_fit_animal <- NULL
}

# Compute shape deltas (only used for params selected by the user)
if (!is.null(lognorm_fit_animal)) {
  animal_mdd_meanlog <- coef(lognorm_fit_animal)[["meanlog"]]
  animal_mdd_sdlog   <- coef(lognorm_fit_animal)[["sdlog"]]
  delta_meanlog <- animal_mdd_meanlog - plant_mdd_meanlog
  ratio_sdlog   <- animal_mdd_sdlog / plant_mdd_sdlog
} else {
  delta_meanlog <- 0
  ratio_sdlog   <- 1
}

# Override deltas for params the user did NOT select
if (!"meanlog" %in% fit_params) delta_meanlog <- 0
if (!"sdlog"   %in% fit_params) ratio_sdlog   <- 1

# Push shape-shift parameters to global env for simulateScaledSpeciesRichness()
.amani_delta_meanlog    <<- delta_meanlog
.amani_ratio_sdlog      <<- ratio_sdlog
.amani_plant_mdd_meanlog <<- plant_mdd_meanlog

message("  Shape-shift parameters (fit_params = [",
        paste(fit_params, collapse = ", "), "]):")
message("    Plant MDD meanlog = ", round(plant_mdd_meanlog, 4),
        ",  sdlog = ", round(plant_mdd_sdlog, 4))
if (!is.null(lognorm_fit_animal)) {
  message("    Animal MDD meanlog = ", round(animal_mdd_meanlog, 4),
          ",  sdlog = ", round(animal_mdd_sdlog, 4))
}
message("    delta_meanlog = ", round(delta_meanlog, 4),
        " (", ifelse("meanlog" %in% fit_params, "ACTIVE", "off — using plant shape"), ")")
message("    ratio_sdlog   = ", round(ratio_sdlog, 4),
        " (", ifelse("sdlog" %in% fit_params, "ACTIVE", "off — using plant shape"), ")")

# Generate the shape-shifted plant curve for plotting
# (applies the same transform that simulateScaledSpeciesRichness will use)
if (delta_meanlog != 0 || ratio_sdlog != 1) {
  fit_pred_shifted <- fit_pred |>
    dplyr::mutate(
      log_x_new = plant_mdd_meanlog + delta_meanlog +
                  (log(num_grid_cells) - plant_mdd_meanlog) * ratio_sdlog,
      num_grid_cells_shifted = pmax(1, round(exp(log_x_new)))
    ) |>
    dplyr::summarise(n = sum(n), .by = num_grid_cells_shifted) |>
    dplyr::rename(num_grid_cells = num_grid_cells_shifted) |>
    dplyr::arrange(num_grid_cells)
  fit_pred_for_scale <- fit_pred_shifted
} else {
  fit_pred_for_scale <- fit_pred
}

# Scale the fitted (and optionally shape-shifted) curve to match animal amplitude
fit_pred <- fit_pred |>
  dplyr::mutate(fitted_scaled = n * (animal_max / max(n)))

if (exists("fit_pred_shifted")) {
  fit_pred_shifted <- fit_pred_shifted |>
    dplyr::mutate(fitted_scaled = n * (animal_max / max(n)))
}




####now plot it

# Add raw animal data for reference
animal_df <- morphotype_distribution %>%
  mutate(group = "Animal morphotypes")

timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
output_filename <- paste0("animal_Spat_Ext_model.", timestamp, ".pdf")
output_path <- file.path(output_folder, output_filename)

# Open PDF device
pdf(output_path, width = 8, height = 6)

p <- ggplot() +
  geom_bar(data = animal_df, aes(x = factor(num_grid_cells), y = n, fill = group),
           stat = "identity", position = "identity", alpha = 0.7, color = "black") +
  geom_line(data = fit_pred, aes(x = num_grid_cells, y = fitted_scaled),
            color = "darkred", linewidth = 1.3) +
  scale_fill_viridis_d(option = "magma", direction = -1, name = NULL) +
  theme_minimal() +
  theme(legend.position = "top")

# If shape-shifting was applied, overlay the shifted curve and annotate
if (exists("fit_pred_shifted")) {
  p <- p +
    geom_line(data = fit_pred_shifted,
              aes(x = num_grid_cells, y = fitted_scaled),
              color = "steelblue", linewidth = 1.3, linetype = "dashed") +
    labs(
      title = "Animal Spatial Extent: Amplitude-Only vs Shape-Shifted Fit",
      x = "Number of Grid Cells",
      y = "Number of Morphotypes",
      caption = paste0("Red = amplitude-only (plant shape); ",
                       "Blue dashed = shape-shifted (",
                       paste(fit_params, collapse = "+"), ")")
    )
} else {
  p <- p +
    labs(
      title = "Animal Morphotype Grid Distribution with Rescaled Log-normal Fit",
      x = "Number of Grid Cells",
      y = "Number of Morphotypes",
      caption = "Log-normal curve originally fit to plant data and rescaled to animal amplitude"
    )
}

print(p)
dev.off()




# Available grid cells — restrict to cells within the Amazon boundary
# (grid_sf_clipped contains only cells that intersect the Amazon polygon)
available_grid_ids <- grid_sf$grid_id[grid_sf$grid_id %in% grid_sf_clipped$grid_id]
available_grid_ids <- available_grid_ids[!is.na(available_grid_ids)]
head(available_grid_ids)
length(available_grid_ids)
length(unique(available_grid_ids))
#calculate a weight matrix of the available_grid_ids to use for drawing the starting cell in the simulation.
#this adds a spatially explicit element to the model
species_counts$grid_id

# drop geometry if species_counts is sf, to avoid overhead
# this ensures we are only working with attributes
species_counts_dt <- data.table::as.data.table(
  sf::st_drop_geometry(species_counts)
)

# count how many times each grid_id appears in species_counts
# this gives frequency per grid cell
grid_freq <- species_counts_dt[
  , .N,                                  # .N = count of rows
  by = grid_id                           # group by grid cell
]


# ensure all available_grid_ids are represented
# grids with zero occurrences get frequency 0
grid_freq_full <- data.table::data.table(
  grid_id = available_grid_ids           # all possible grid IDs
)[
  grid_freq,                             # join observed frequencies
  on = "grid_id"
][
  is.na(N), N := 0                       # replace NA with 0 for missing grids
]

nrow(grid_freq)
length(available_grid_ids)

# extract the weights vector in the same order as available_grid_ids
# this is directly usable as prob = in sample()
sampling_weights <- grid_freq_full[
  match(available_grid_ids, grid_id),    # enforce identical ordering
  N                                      # frequencies as weights
]


# replace any NA values with the minimum weight
# this allows simulated species to be put anywhere in the grid, just with low probability
sampling_weights[is.na(sampling_weights)] <-  0  # ensures no NA probabilities
sampling_weights<-sampling_weights+1 #change the probabilities so the simulated species can go anywhere in the grid, at least with low probability

# optional sanity check: confirm no NAs remain
anyNA(sampling_weights)                                # should return FALSE

# create a neutral (uniform) sampling weight vector
# every grid cell has equal probability
sampling_weights_neutral <- rep(
  1,
  length(available_grid_ids)        # one weight per grid cell
)

# optional sanity check
length(sampling_weights_neutral) == length(available_grid_ids)


####ok, finally....now take the scaling factor above (animal_max / max(fitted)) and use it to scale all the plant data.
# # should be the same as above but with a progress bar
# 
###original version of the function that works
# simulateScaledSpeciesRichness <- function(fit_pred, scaling_factor) {
#   
#   #plant_fit_max <- max(fit_pred$fitted)
#   #scaling_factor <- animal_max / plant_fit_max
#   
#   # Apply scaling
#   scaled_plant_distribution <- plant_distribution %>%
#     mutate(n_scaled = n * scaling_factor)
#   
#   # Available grid cells
#   available_grid_ids <- grid_sf$grid_id
#   available_grid_ids <- available_grid_ids[!is.na(available_grid_ids)]
#   
#   # Initialize list to store simulated species assignments
#   sim_species_records <- list()
#   species_counter <- 1
#   
#   ## 🔹 Progress bar setup (species-level)
#   total_species <- sum(round(scaled_plant_distribution$n_scaled), na.rm = TRUE)
#   
#   pb <- progress_bar$new(
#     format = "  Simulating species [:bar] :current/:total (:percent) eta: :eta",
#     total  = total_species,
#     clear  = FALSE,
#     width  = 60
#   )
#   
#   ## 🔹 Simulation loops
#   for (i in seq_len(nrow(scaled_plant_distribution))) {
#     
#     num_sites   <- scaled_plant_distribution$num_grid_cells[i]
#     num_species <- round(scaled_plant_distribution$n_scaled[i])
#     
#     for (j in seq_len(num_species)) {
#       
#       pb$tick()  # ← progress update
#       
#       if (num_sites <= length(available_grid_ids)) {
#         
#         sim_name <- str_glue("SimSpp{str_pad(species_counter, 7, pad = '0')}")
#         
#         assigned_cells <- sample(available_grid_ids, num_sites, replace = FALSE)
#         
#         sim_occurrences <- data.frame(
#           scientificName = rep(sim_name, num_sites),
#           grid_id = assigned_cells,
#           stringsAsFactors = FALSE
#         )
#         
#         sim_species_records[[species_counter]] <- sim_occurrences
#         species_counter <- species_counter + 1
#       }
#     }
#   }
#   
#   # Combine into one large data frame
#   sim_species_df <<- bind_rows(sim_species_records)
#   
#   # Join to grid_sf to get geometry
#   sim_species_sf <<- sim_species_df %>%
#     left_join(grid_sf, by = "grid_id") %>%
#     st_as_sf()
#   
#   ####the line below gives how many species of cockroaches there are under this simulation
#   output<-length(unique(sim_species_sf$scientificName))
#   return(output)
#   
#  #return(sim_species_records)
# }

####new version split into two functions
##inner loop function
# simulate_species_records <- function(
#     scaled_plant_distribution,
#     available_grid_ids,
#     show_progress = TRUE
# ) {
#   rm(sim_records)
#   sim_species_records <- vector("list", length = sum(round(
#     scaled_plant_distribution$n_scaled
#   ), na.rm = TRUE))
#   
#   species_counter <- 1
#   
#   ## Progress bar (optional)
#   if (show_progress) {
#     total_species <- length(sim_species_records)
#     
#     pb <- progress::progress_bar$new(
#       format = "  Simulating species [:bar] :current/:total (:percent) eta: :eta",
#       total  = total_species,
#       clear  = FALSE,
#       width  = 60
#     )
#   }
#   
#   for (i in seq_len(nrow(scaled_plant_distribution))) {
#     
#     num_sites   <- scaled_plant_distribution$num_grid_cells[i]
#     num_species <- round(scaled_plant_distribution$n_scaled[i])
#     
#     for (j in seq_len(num_species)) {
#       
#       if (show_progress) pb$tick()
#       
#       if (num_sites <= length(available_grid_ids)) {
#         
#         sim_name <- stringr::str_glue(
#           "SimSpp{stringr::str_pad(species_counter, 7, pad = '0')}"
#         )
#         
#         assigned_cells <- sample(
#           available_grid_ids,
#           num_sites,
#           replace = FALSE
#         )
#         
#         sim_species_records[[species_counter]] <- data.frame(
#           scientificName = rep(sim_name, num_sites),
#           grid_id = assigned_cells,
#           stringsAsFactors = FALSE
#         )
#         
#         species_counter <- species_counter + 1
#       }
#     }
#   }
#   
#   sim_records<<-sim_species_records[seq_len(species_counter - 1)] #this might be wrong. if the line below is right then this line is unnecessary but i am keeping it in case it is actually correct
#   return(sim_species_records[seq_len(species_counter - 1)])
# }
#######Same as the above function but with simulated species done with local clustering
#
# Precomputes a sparse adjacency matrix and derives an adjacency list for fast
# contiguous range growth. At each step, the frontier (all 8-connected
# neighbours of the current range not already in the range) is maintained
# incrementally and one frontier cell is sampled proportional to its
# plant-occurrence weight.

build_adjacency <- function(grid_sf, grid_size_m) {
  # Build an 8-connected adjacency list from cell centroid coordinates.
  # Returns a list of length nrow(grid_sf) where each element is an integer
  # vector of row indices of neighbouring cells.

  require(Matrix)

  cx <- grid_sf$cx
  cy <- grid_sf$cy
  n  <- length(cx)
  thresh <- grid_size_m * 1.5  # captures orthogonal + diagonal neighbours

  message("  Building adjacency structure for ", n, " grid cells...")

  # Build sparse adjacency matrix in chunks to limit memory from outer()
  chunk_size <- 2000
  from_list <- vector("list", ceiling(n / chunk_size))
  to_list   <- vector("list", ceiling(n / chunk_size))

  for (ch in seq_len(ceiling(n / chunk_size))) {
    idx_start <- (ch - 1) * chunk_size + 1
    idx_end   <- min(ch * chunk_size, n)
    chunk_idx <- idx_start:idx_end

    dx <- outer(cx[chunk_idx], cx, function(a, b) abs(a - b))
    dy <- outer(cy[chunk_idx], cy, function(a, b) abs(a - b))

    is_nb <- dx <= thresh & dy <= thresh
    for (r in seq_along(chunk_idx)) {
      is_nb[r, chunk_idx[r]] <- FALSE
    }

    wh <- which(is_nb, arr.ind = TRUE)
    from_list[[ch]] <- chunk_idx[wh[, 1]]
    to_list[[ch]]   <- wh[, 2]
  }

  from_vec <- unlist(from_list)
  to_vec   <- unlist(to_list)

  adj_mat <- sparseMatrix(i = from_vec, j = to_vec, x = 1, dims = c(n, n))

  # Extract adjacency list from column-compressed sparse matrix (fast)
  adj_list <- vector("list", n)
  for (i in seq_len(n)) {
    adj_list[[i]] <- adj_mat@i[adj_mat@p[i]:(adj_mat@p[i + 1] - 1)] + 1L
  }

  message("  Adjacency: ", n, " cells, ~", round(mean(lengths(adj_list)), 1),
          " neighbours/cell.")
  adj_list
}


simulate_species_records_clustered <- function(
    scaled_plant_distribution,
    available_grid_ids,
    grid_sf,
    adj_list,
    cell_weights,
    sampling_weights = sampling_weights_neutral,
    show_progress = TRUE,
    jitter_range  = TRUE,    # if TRUE, Poisson-jitter each species target range size
    estimateBounds=TRUE                       # so species within the same bin get different sizes
                            # (eliminates point-mass peaks in the output distribution)
) {
  # adj_list:     list of integer vectors — adj_list[[i]] = row indices of
  #               neighbours of cell i in grid_sf.
  # cell_weights: named numeric vector keyed by grid_id — plant-occurrence
  #               weights used to bias frontier cell selection during range growth.

  rm(sim_records)
  sim_species_records <- vector("list", length = sum(round(
    scaled_plant_distribution$n_scaled
  ), na.rm = TRUE))

  species_counter <- 1
  n_cells <- nrow(grid_sf)
  all_ids <- grid_sf$grid_id

  # Weight vector aligned to grid_sf row order
  w_vec <- cell_weights[as.character(all_ids)]
  w_vec[is.na(w_vec)] <- 0

  # ── Frontier exhaustion diagnostics ────────────────────────────────────────────
  # Track species where contiguous growth hits the reachable frontier before
  # the target range size is met. Silently truncated species create artificial
  # peaks at the maximum connected-region size.
  n_exhausted     <- 0L   # species whose frontier ran out early
  total_shortfall <- 0L   # cumulative cells short across all truncated species
  max_shortfall   <- 0L   # worst-case single-species shortfall

  ## Progress bar
  if (show_progress) {
    total_species <- length(sim_species_records)
    pb <- progress::progress_bar$new(
      format = "  Simulating species [:bar] :current/:total (:percent) eta: :eta",
      total  = total_species,
      clear  = FALSE,
      width  = 60
    )
  }
  
  #define this so we don't have to recalculate it inside the loop
  maxCellRange<-length(available_grid_ids)
  # Pre-allocations outside species loops — avoids per-species alloc/reset
  visit_gen         <- integer(n_cells)   # in-range stamp: cell in range when visit_gen[r] == current_gen
  frontier_gen      <- integer(n_cells)   # frontier stamp: cell in frontier when frontier_gen[r] == current_gen
  current_gen       <- 0L
  frontier_rows_buf <- integer(n_cells)   # active frontier row indices in positions 1:flen
  frontier_pos      <- integer(n_cells)   # frontier_pos[r] = position of r in frontier_rows_buf
  max_id            <- max(all_ids)
  id_to_row         <- integer(max_id)
  id_to_row[all_ids] <- seq_len(n_cells) # O(1) grid_id -> row lookup
  #all_sim_names     <- sprintf("SimSpp%07d", seq_len(length(sim_species_records)))
  
  for (i in seq_len(nrow(scaled_plant_distribution))) {

    num_sites   <- scaled_plant_distribution$num_grid_cells[i]
    num_species <- round(scaled_plant_distribution$n_scaled[i])
    if (is.na(num_sites) || is.na(num_species) || num_species < 1) next
    if (is.na(num_sites) || is.na(num_species) || num_species < 1) next
    
    for (j in seq_len(num_species)) {

      if (show_progress) pb$tick()

      # ── Per-species range size jitter ────────────────────────────────────
      # Draw target from Poisson(lambda = bin representative) so species
      # within the same bin receive different range sizes, dissolving the
      # point-mass peaks created by assigning all species in a bin the
      # identical median value.
      ##smooth out range sizes larger than 50 cells (which is the same threshold for log binning below)
      # #I added this condition because the binning
      ##behavior for larger ranged species was resulting in
      #  #multi-modal distribution. but doing the poisson jitter
      # #for smaller ranges smooths out the range sizes in an unnatrual way (that doesn't fit the log-normal distribution)
       
      num_sites_j <- if (jitter_range && num_sites > 50) {
        lambda <- min(num_sites, maxCellRange)          # cap lambda before rpois to prevent integer overflow
        max(1L, min(as.integer(round(rpois(1L, lambda))), maxCellRange + 1L))
      } else {
        as.integer(min(num_sites, maxCellRange + 1L))   # also cap the non-jitter path
      }
      
      
      if (num_sites_j <= maxCellRange) {

        sim_name <- stringr::str_glue(
          "SimSpp{stringr::str_pad(species_counter, 7, pad = '0')}"
        )

        # Draw starting cell using plant-occurrence weights
        starting_cell <- sample(available_grid_ids, 1,
                                prob = sampling_weights, replace = FALSE)
        start_row <- match(starting_cell, all_ids)

        if (num_sites_j == 1) {
          assigned_cells <- starting_cell
        } else {
          # --- Contiguous weighted growth via adjacency list ---
          # Maintain the frontier incrementally: when a cell is added to the
          # range, remove it from the frontier and add its new neighbors.
          # Each frontier cell is sampled proportional to its plant weight.

          range_rows <- integer(num_sites_j)
          range_rows[1] <- start_row
          #########################
          # in_range <- logical(n_cells)
          # in_range[start_row] <- TRUE
          # actual_size <- 1L
          # 
          # # Initialize frontier from neighbors of starting cell
          # frontier_set <- adj_list[[start_row]]
          # frontier_set <- frontier_set[!in_range[frontier_set]]

          ###########################The code block below replaces the commented-out version above
          current_gen          <- current_gen + 1L
          visit_gen[start_row] <- current_gen
          actual_size          <- 1L
          flen                 <- 0L
          
          init_nbs <- adj_list[[start_row]]
          init_nbs <- init_nbs[visit_gen[init_nbs] != current_gen]
          flen     <- length(init_nbs)
          if (flen > 0L) {
            frontier_rows_buf[seq_len(flen)] <- init_nbs
            frontier_gen[init_nbs]           <- current_gen
            frontier_pos[init_nbs]           <- seq_len(flen)
          }
          #########################
          
          #########################
          # for (step in 2:num_sites_j) {
          #   if (length(frontier_set) == 0) break
          # 
          #   fw <- w_vec[frontier_set]
          #   if (sum(fw) == 0) fw <- rep(1, length(frontier_set))
          # 
          #   new_row <- if (length(frontier_set) == 1) {
          #     frontier_set
          #   } else {
          #     sample(frontier_set, 1, prob = fw)
          #   }
          # 
          #   actual_size <- actual_size + 1L
          #   range_rows[actual_size] <- new_row
          #   in_range[new_row] <- TRUE
          # 
          #   # Update frontier: remove new_row, add its new neighbours
          #   frontier_set <- frontier_set[frontier_set != new_row]
          #   new_nbs <- adj_list[[new_row]]
          #   new_nbs <- new_nbs[!in_range[new_nbs]]
          #   frontier_set <- unique(c(frontier_set, new_nbs))
          # }
          ###########################The code block below replaces the commented-out version above
          for (step in 2:num_sites_j) {
            if (flen == 0L) break
            
            active <- frontier_rows_buf[seq_len(flen)]
            fw     <- w_vec[active]
            if (sum(fw) == 0) fw[] <- 1
            
            new_row <- if (flen == 1L) active else sample(active, 1L, prob = fw)
            
            actual_size             <- actual_size + 1L
            range_rows[actual_size] <- new_row
            visit_gen[new_row]      <- current_gen
            
            # O(1) swap-and-pop removal from frontier buffer
            pos <- frontier_pos[new_row]
            if (pos < flen) {
              last_row               <- frontier_rows_buf[flen]
              frontier_rows_buf[pos] <- last_row
              frontier_pos[last_row] <- pos
            }
            frontier_pos[new_row] <- 0L
            frontier_gen[new_row] <- 0L
            flen <- flen - 1L
            
            # Add unvisited non-frontier neighbours (vectorised)
            new_nbs <- adj_list[[new_row]]
            add_nbs <- new_nbs[
              visit_gen[new_nbs]    != current_gen &
                frontier_gen[new_nbs] != current_gen
            ]
            n_add <- length(add_nbs)
            if (n_add > 0L) {
              new_pos                    <- (flen + 1L):(flen + n_add)
              frontier_rows_buf[new_pos] <- add_nbs
              frontier_pos[add_nbs]      <- new_pos
              frontier_gen[add_nbs]      <- current_gen
              flen                       <- flen + n_add
            }
          }
          #########################
          # ── Log frontier exhaustion ───────────────────────────────────────────
          if (actual_size < num_sites_j) {
            shortfall       <- num_sites_j - actual_size
            n_exhausted     <- n_exhausted + 1L
            total_shortfall <- total_shortfall + shortfall
            if (shortfall > max_shortfall) max_shortfall <- shortfall
          }

          assigned_cells <- all_ids[range_rows[1:actual_size]]
        }

        sim_species_records[[species_counter]] <- data.frame(
          scientificName = rep(sim_name, length(assigned_cells)),
          grid_id = assigned_cells,
          stringsAsFactors = FALSE
        )

        species_counter <- species_counter + 1
      }
    }
  }

  sim_records <<- sim_species_records[seq_len(species_counter - 1)]

  # ── Frontier exhaustion summary ───────────────────────────────────────────────
  n_total <- species_counter - 1L
  if (n_exhausted > 0L) {
    message(sprintf(
      "  Frontier exhaustion: %d / %d species (%.1f%%) had ranges truncated.",
      n_exhausted, n_total, 100 * n_exhausted / n_total
    ))
    message(sprintf(
      "  Mean shortfall: %.1f cells | Max shortfall: %d cells | Total cells lost: %d",
      total_shortfall / n_exhausted, max_shortfall, total_shortfall
    ))
  } else {
    message("  Frontier exhaustion: none detected.")
  }
}




#wrapper function
simulateScaledSpeciesRichness <- function(fit_pred, extent_scaling_factor, sampling_weights = sampling_weights_neutral) {

  print(paste0("Simulating Amazonian species with an extent scaling factor of: ", extent_scaling_factor))
  
  
#   scaling_factor<-0.065
# # Apply scaling
# scaled_plant_distribution <- plant_distribution %>% #fit_pred should NOT be here. This MUST be plant_distribution, which is the same as fit_pred but for the whole amazon. 
#   dplyr::mutate(n_scaled = n * scaling_factor)

###leaving scaled_plant_distribution as it is above will be a problem because it results in us dropping more species with large ranges than expected (e.g., if the scaling factor is 0.5 then we drop all range sizes occupied by <2 species. When the scaling factor is very small (e.g. <0.1) this could result in over 1000 species dropped)

# coerce plant_distribution to data.table without copying if possible
pd <- data.table::as.data.table(plant_distribution)             # expects column num_grid_cells

# ── Log-spaced binning ──────────────────────────────────────────────────
# Replace the old hard-edged tiered binning (bin_num_grids) with
# geometrically spaced bin edges. This avoids the abrupt bin-width
# transitions at 50/200/600 that created artificial density peaks
# ("attractor values") in the simulated range-size distribution.
#
# Log-spacing naturally gives fine resolution at small range sizes
# (where most species cluster) and progressively coarser resolution
# at large range sizes — matching the spirit of the original scheme
# but without hard discontinuities.

n_log_bins <- 80                                                  # tunable; ~80 gives similar
                                                                  # compression to old scheme
max_cells  <- max(pd$num_grid_cells, na.rm = TRUE)

# keep range sizes 1–50 unbinned (exact), bin everything above
# using geometrically spaced edges from 50 to max
unbinned_threshold <- 50

# rows that stay unbinned (small range sizes, high ecological importance)
pd_small <- pd[num_grid_cells <= unbinned_threshold]
pd_large <- pd[num_grid_cells >  unbinned_threshold]

if (nrow(pd_large) > 0) {
  # create log-spaced bin edges from the threshold to the max value
  log_edges <- exp(seq(
    log(unbinned_threshold),
    log(max_cells + 1),                                           # +1 so max value falls inside last bin
    length.out = n_log_bins + 1
  ))

  # assign each large-range row to a log-spaced bin
  pd_large[, bin_id := findInterval(num_grid_cells, log_edges, rightmost.closed = TRUE)]

  # aggregate: median of original values as representative, sum of counts
  pd_large_binned <- pd_large[
    , .(
      num_grid_cells = round(median(num_grid_cells)),
      n = sum(n)
    ),
    by = bin_id
  ][order(num_grid_cells)]
  pd_large_binned[, bin_id := NULL]

  # recombine small (exact) and large (log-binned) portions
  plant_distribution_binned <- data.table::rbindlist(
    list(pd_small[, .(num_grid_cells, n)], pd_large_binned)
  )[order(num_grid_cells)]

} else {
  # edge case: all range sizes <= threshold (nothing to bin)
  plant_distribution_binned <- pd_small[, .(num_grid_cells, n)]
}

# sanity check: total species count is conserved
stopifnot(abs(sum(plant_distribution_binned$n) - sum(pd$n)) < 0.5)

# ── Apply shape transformation (meanlog shift + sdlog rescaling) ────────
# If the user selected "meanlog" or "sdlog" in fit_params, the empirical
# plant distribution's range sizes are warped in log-space before amplitude
# scaling. This changes the proportion of endemic vs widespread species
# while preserving the empirical texture of the distribution.
delta_ml  <- if (exists(".amani_delta_meanlog",     envir = .GlobalEnv)) .amani_delta_meanlog     else 0
ratio_sd  <- if (exists(".amani_ratio_sdlog",       envir = .GlobalEnv)) .amani_ratio_sdlog       else 1
ref_ml    <- if (exists(".amani_plant_mdd_meanlog", envir = .GlobalEnv)) .amani_plant_mdd_meanlog else log(2)

if (delta_ml != 0 || ratio_sd != 1) {
  message("  Applying shape transformation: delta_meanlog = ", round(delta_ml, 4),
          ", ratio_sdlog = ", round(ratio_sd, 4))
  plant_distribution_binned <- as.data.frame(plant_distribution_binned)
  plant_distribution_binned$log_x_new <-
    ref_ml + delta_ml +
    (log(plant_distribution_binned$num_grid_cells) - ref_ml) * ratio_sd
  # plant_distribution_binned$num_grid_cells <-
  #   pmax(1L, as.integer(round(exp(plant_distribution_binned$log_x_new))))
  # plant_distribution_binned$log_x_new <- NULL
  # 
  # # Re-aggregate: multiple original bins may map to the same integer after transform
  # plant_distribution_binned <- plant_distribution_binned |>
  #   dplyr::summarise(n = sum(n), .by = num_grid_cells) |>
  #   dplyr::arrange(num_grid_cells)
  
  plant_distribution_binned$num_grid_cells <-
    pmax(1, round(exp(plant_distribution_binned$log_x_new)))  # drop as.integer()
  plant_distribution_binned$log_x_new <- NULL
  
  plant_distribution_binned <- plant_distribution_binned |>
    dplyr::summarise(n = sum(n), .by = num_grid_cells) |>
    dplyr::filter(!is.na(num_grid_cells)) |>          # drop overflow rows
    dplyr::arrange(num_grid_cells)
  
  
  
  
}

scaled_plant_distribution <- plant_distribution_binned %>%
     dplyr::mutate(n_scaled = ceiling(n * extent_scaling_factor))


# Run simulation
#options(error = NULL)

grid_sf$centroid <- sf::st_centroid(grid_sf$geometry)      # POINT geometry per grid cell
# extract numeric centroid coordinates for fast filtering
cent_coords <- sf::st_coordinates(grid_sf$centroid)         # matrix [n x 2]

# drop polygon geometry from the object used in simulation
# keeps memory footprint smaller when passing to workers
grid_centroids <- grid_sf[, c("grid_id")]                  # keep only IDs
grid_centroids$geometry <- grid_sf$centroid                # replace geometry with centroids

# store coordinates directly on grid_sf for O(1) access
grid_centroids$cx <- cent_coords[, 1]                               # centroid x
grid_centroids$cy <- cent_coords[, 2]                               # centroid y


grid_sub <- grid_centroids[grid_centroids$grid_id %in% available_grid_ids, ]
grid_sf$centroid <- NULL  # drop centroid column added above
rm(grid_centroids, cent_coords); gc()

# Load or build adjacency list (8-connected neighbours)
adj_cache <- file.path(project_root, "output", "shared", "adjacency_list_amazon.rds")
if (file.exists(adj_cache)) {
  message("  Loading cached adjacency list...")
  adj_saved <- readRDS(adj_cache)
  # Verify it matches the current grid
  if (identical(adj_saved$grid_ids, grid_sub$grid_id) &&
      adj_saved$grid_size_m == grid_size_m) {
    adj_list <- adj_saved$adj_list
    message("  Adjacency list loaded (", adj_saved$n_cells, " cells).")
  } else {
    message("  Cached adjacency list doesn't match current grid. Rebuilding...")
    adj_list <- build_adjacency(grid_sub, grid_size_m)
  }
  rm(adj_saved)
} else {
  adj_list <- build_adjacency(grid_sub, grid_size_m)
  # Cache for other organisms
  dir.create(dirname(adj_cache), showWarnings = FALSE, recursive = TRUE)
  saveRDS(list(adj_list = adj_list, grid_ids = grid_sub$grid_id,
               n_cells = nrow(grid_sub), grid_size_m = grid_size_m), adj_cache)
  message("  Adjacency list cached to: ", adj_cache)
}

# Build named weight vector keyed by grid_id for frontier weighting
cell_weights <- setNames(sampling_weights, as.character(available_grid_ids))

simulate_species_records_clustered(
  scaled_plant_distribution = scaled_plant_distribution,
  available_grid_ids = available_grid_ids,
  grid_sf = grid_sub,
  adj_list = adj_list,
  cell_weights = cell_weights,
  sampling_weights = sampling_weights,
  estimateBounds = .amani_estimateBounds)

# Combine into one large data frame
sim_species_df <- dplyr::bind_rows(sim_records)

# Join to grid_sf to get geometry
sim_species_sf <- sim_species_df %>%
  dplyr::left_join(grid_sf, by = "grid_id") %>%
  sf::st_as_sf()

# Number of simulated species
length(unique(sim_species_sf$scientificName))
}



####the line below gives how many species of animals there are under this simulation
####This is just a test, so it's commented out

#print("Starting a simulation...this will take long if the animal has a lot of data...")
#plant_fit_max <- max(fit_pred$fitted)
#simulateScaledSpeciesRichness(fit_pred, scaling_factor = (animal_max / plant_fit_max))

###there are no error ranges on this because I am using a single scaling factor. 


###now i am going to repeat this using the plant sp. range extent distribution unmodified. if it yields roughly the same spp. richness as the raw plant data then the process is working correctly.


########Testing simulateScaledSpeciesRichness with the plant data to see if it gets the right answer#############
# ------------------------------------------------------------
# only run the estimator sanity check if diagnostics are enabled
# this prevents expensive simulations from running in production
# ------------------------------------------------------------
if (exists("diagnostic") && isTRUE(diagnostic)) {
  
  # message to clearly indicate diagnostic mode is active
  print("checking that the estimator is working by replicating the empirical plant richness")
  
  # show the head of the plant spatial-extent distribution as a quick visual sanity check
  head(plant_distribution)
  
  # initialize an empty list to store simulated species records (not strictly needed,
  # but kept here to mirror the structure of the main simulation logic)
  sim_species_recordsRaw <- list()
  
  # initialize species counter (used internally by the simulator)
  species_counter <- 1
  
  # compute the maximum fitted value from the plant extent model
  # this represents the reference scale of the empirical distribution
  plant_fit_max <- max(fit_pred$fitted)
  
  # define a neutral scaling factor (plant scaled to itself)
  # this should reproduce the original plant richness if the estimator is correct
  extent_scaling_factor <- plant_fit_max / plant_fit_max
  
  # apply scaling to the empirical plant distribution
  # n_scaled should equal n under this neutral scaling
  
  scaled_plant_distribution <- plant_distribution %>%
    mutate(n_scaled = n * extent_scaling_factor)
  
  # run the spatial-extent-based richness simulator
  # this should return a richness value comparable to the empirical data
  sim_speciesRaw_df <- simulateScaledSpeciesRichness(
    scaled_plant_distribution,
    extent_scaling_factor = 1
  )
  
  # check whether simulated richness matches empirical plant richness
  # unique(unlist()) is defensive in case scientificName is a list-column
  is_equal <- sim_speciesRaw_df ==
    length(unique(unlist(species_counts$scientificName)))
  
  # compose a human-readable diagnostic message
  output_text <- if (is_equal) {
    "The spatial extent based richness estimator is working."
  } else {
    "The spatial extent based richness estimator is NOT working."
  }
  
  # generate a timestamp for reproducible diagnostic outputs
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  
  # construct an informative output filename
  output_filename <- paste0(
    "spatEstEst.sanityCheck.",
    timestamp,
    ".txt"
  )
  
  # build the full output path inside the results directory
  path00 <- file.path(output_folder, output_filename)
  
  # write the diagnostic result to disk for later inspection
  writeLines(output_text, con = path00)
  
}


###Good. when i ran it on April 15 and 17 2025 the numbers are exactly the same.

###ok now I am going to use the simulated species extents to do a randomization simulation to get 95% CIs

############Estimate gamma from total spatial extents##################
head(species_counts)

# # Sum total abundance per grid cell
# grid_abundance <- species_counts %>%
#   group_by(grid_id) %>%
#   summarise(
#     total_abundance = sum(count),
#     geometry = st_geometry(geometry)[1],  # Keep one representative geometry
#     .groups = "drop"
#   ) %>%
#   st_as_sf()
# 
# 
# set.seed(SEEDnumber)  # for reproducibility
# 
# # 1. Setup
# 
# simulated_records <- list()  # initialize an empty list
# record_counter <- 1
# 
# head(sim_species_sf)
# # Loop through all grid cells
# for (i in 1:nrow(grid_abundance)) {
#   
#   this_abundance <- grid_abundance$total_abundance[i]
#   this_grid_id <- grid_abundance$grid_id[i]
#   
#   # Skip if abundance is zero or NA (just in case)
#   if (is.na(this_abundance) || this_abundance == 0) next
#   
#   # Sample species names
#   sampled_species <- sample(sim_species_sf$scientificName, size = this_abundance, replace = TRUE)
#   
#   # Create a data.frame of simulated occurrences
#   df <- data.frame(
#     grid_id = rep(this_grid_id, this_abundance),
#     scientificName = sampled_species,
#     stringsAsFactors = FALSE
#   )
#   
#   # Store in list
#   simulated_records[[record_counter]] <- df
#   record_counter <- record_counter + 1
# }

plant_fit_max <- max(fit_pred$n)


###allows specifying a custom scaling factor (mostly for troubleshooting)
if ( is.null( customScale)) {
 
  extent_sf_here <- (animal_max / plant_fit_max)
  
} else {
  
  extent_sf_here <-customScale
  
}


# run spatially explicit or neutral simulation depending on the flag
if (spatial == "explicit") {
    # spatially explicit simulation (uses empirical / distance-weighted structure)
  gammaRichnessSPRIch <- simulateScaledSpeciesRichness(
    fit_pred,
    extent_scaling_factor = extent_sf_here,
    sampling_weights = sampling_weights
  )
  
} else {
  
  # non-spatial / neutral simulation (uniform sampling weights)
  
  gammaRichnessSPRIch <- simulateScaledSpeciesRichness(
    fit_pred,
    extent_scaling_factor = extent_sf_here,
    sampling_weights = sampling_weights_neutral
  )
  
}

# Combine all rows into one data frame
simulated_occurrence_df <- do.call(rbind, sim_records)
head(simulated_occurrence_df)

###export this reslt, as it is the main output of this model, and can be used downstream to estimate other things
print("Joining simulated species with spatial geometry for exportation...")

geom_tbl <- as.data.table(species_counts)[
  , .SD[1], by = grid_id
][, .(grid_id, geometry)]


simulated_occurrence_sf_TOEXPORT <-
  as.data.table(simulated_occurrence_df)[
    geom_tbl,
    on = "grid_id"
  ] |>
  st_as_sf()
head(simulated_occurrence_sf_TOEXPORT)

timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
output_filename <- paste0("simulatedOccurencesForAnimal.", timestamp, ".csv")
output_path <- file.path(output_folder, output_filename)
write.csv(simulated_occurrence_sf_TOEXPORT, output_path, row.names = FALSE)
saveRDS(simulated_occurrence_sf_TOEXPORT, paste0(output_folder,"/simulatedOccurencesForAnimal.rds"))

print("...Simulated species with spatial geometry exported to CSV and RData (see...simulatedOccurencesForAnimal...files)")


gammaRichnessEst<-length(unique(simulated_occurrence_df$scientificName))
print(gammaRichnessEst)


###########Get 90% CI's for the extent-based gamma diversity estimator#################
# Original data
x_vals <- plant_df$num_grid_cells
y_vals <- plant_df$n

# Define x range for prediction
x_pred <- seq(min(x_vals), max(x_vals), by = 0.1)

# Define a function to fit model and return predictions
fit_lognorm <- function(df) {
  x_vals <- df$num_grid_cells
  y_vals <- df$n
  
  weights_norm <- y_vals / sum(y_vals)
  start_meanlog <- sum(weights_norm * log(x_vals + 1))
  start_sdlog   <- sqrt(sum(weights_norm * (log(x_vals + 1) - start_meanlog)^2))
  
  tryCatch({
    fit <- nlsLM(
      y_vals ~ a * dlnorm(x_vals, meanlog, sdlog),
      start = list(a = max(y_vals), meanlog = start_meanlog, sdlog = start_sdlog)
    )
    preds <- predict(fit, newdata = data.frame(x_vals = x_pred))
    return(preds)
  }, error = function(e) {
    return(rep(NA, length(x_pred)))
  })
}

# Bootstrapping: resample and fit multiple times
set.seed(SEEDnumber+1)
n_boot <- 200


# Expand plant_df into individual species occurrences
species_expanded <- plant_df %>%
  uncount(weights = n)  # creates one row per species, based on n

# Check the distribution
table(species_expanded$num_grid_cells)

set.seed(SEEDnumber+2)
n_species <- nrow(species_expanded)




boot_preds <- replicate(n_boot, {

  # Sample species' range sizes from empirical distribution
  boot_sample <- sample(
    species_expanded$num_grid_cells,
    size = ceiling(n_species*sample(c(0.95, 1.05), size = 1)), #to increase the spread of the bootstrapped data, have it draw either 5% less or 5% more species than what's actually in the sample
    replace = TRUE
  )

  # Tabulate into a new frequency table
  boot_df <- as.data.frame(table(boot_sample)) %>%
    rename(num_grid_cells = boot_sample, n = Freq) %>%
    mutate(num_grid_cells = as.integer(as.character(num_grid_cells)),
           group = "Plant species")
  
  # Optional: sort to match original
  boot_df <- boot_df %>% arrange(num_grid_cells)
  
  
  
  fit_lognorm(boot_df)
}, simplify = "matrix")

# Calculate summary statistics across bootstraps
ci_df <- data.frame(
  num_grid_cells = x_pred,
  lower = apply(boot_preds, 1, function(x) quantile(x, 0.001, na.rm = TRUE)),
  upper = apply(boot_preds, 1, function(x) quantile(x, 0.999, na.rm = TRUE)),
  median = apply(boot_preds, 1, function(x) median(x, na.rm = TRUE))

)

x_valsLower <- ci_df$num_grid_cells
y_valsLower <- ci_df$lower

x_valsUpper <- ci_df$num_grid_cells
y_valsUpper <- ci_df$upper

# Fit log-normal model to lower CI
weights_norm_lower <- y_valsLower / sum(y_valsLower)
start_meanlog_lower <- sum(weights_norm_lower * log(x_valsLower + 1))
start_sdlog_lower   <- sqrt(sum(weights_norm_lower * (log(x_valsLower + 1) - start_meanlog_lower)^2))

lognorm_lower_fit <- nlsLM(
  y_valsLower ~ a * dlnorm(x_valsLower, meanlog, sdlog),
    start = list(a = max(y_valsLower), meanlog = start_meanlog_lower, sdlog = start_sdlog_lower)
  )

# lognorm_lower_fit <- nls(
#   y_valsLower ~ a * dlnorm(x_valsLower, meanlog, sdlog),
#   start = list(a = max(y_valsLower), meanlog = mean(log(x_valsLower)), sdlog = sd(log(x_valsLower))),
#   control = nls.control(maxiter = 200, warnOnly = TRUE)
# )

# Fit log-normal model to upper CI
weights_norm_Upper <- y_valsUpper / sum(y_valsUpper)
start_meanlog_Upper <- sum(weights_norm_Upper * log(x_valsUpper + 1))
start_sdlog_Upper   <- sqrt(sum(weights_norm_Upper * (log(x_valsUpper + 1) - start_meanlog_Upper)^2))

lognorm_Upper_fit <- nlsLM(
  y_valsUpper ~ a * dlnorm(x_valsUpper, meanlog, sdlog),
  start = list(a = max(y_valsUpper), meanlog = start_meanlog_Upper, sdlog = start_sdlog_Upper)
)

# lognorm_Upper_fit <- nls(
#   y_valsUpper ~ a * dlnorm(x_valsUpper, meanlog, sdlog),
#   start = list(a = max(y_valsUpper), meanlog = mean(log(x_valsUpper)), sdlog = sd(log(x_valsUpper))),
#   control = nls.control(maxiter = 200, warnOnly = TRUE)
# )


# Predict over range
x_pred <- ci_df$num_grid_cells
y_pred <- predict(lognorm_lower_fit, newdata = data.frame(x_vals = x_pred))

y_predUpper <- predict(lognorm_Upper_fit, newdata = data.frame(x_vals = x_pred))

print("Exporting the plant extent model...")

timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
output_filename <- paste0("plantPatExtentModel.MDD.", timestamp, ".pdf")
output_path <- file.path(output_folder, output_filename)

# Open PDF device
pdf(output_path, width = 8, height = 6)

# Plot
p<-ggplot(ci_df, aes(x = num_grid_cells)) +
  geom_line(aes(y = y_pred), color = "steelblue", size = 1.2, linetype = "dashed") +
  geom_line(aes(y = y_predUpper), color = "steelblue", size = 1.2, linetype = "dashed") +
  labs(
    title = "Log-normal Fit to Lower Confidence Interval",
    x = "Number of Grid Cells",
    y = "Lower 95% CI",
    caption = "Dashed lines = fitted log-normal curve of upper and lower CI"
  ) +
  geom_point(data = plant_df, aes(x = num_grid_cells, y = n), size = 3, color = "darkgreen") +
  geom_ribbon(data = ci_df, aes(x = num_grid_cells, ymin = lower, ymax = upper), fill = "darkgreen", alpha = 0.2) +
  geom_line(data = ci_df, aes(x = num_grid_cells, y = median), color = "darkgreen", size = 1.2) +
  scale_x_log10() +
  scale_y_log10() +
  labs(
    title = "Log-normal Fit to Plant Species Distribution with Bootstrapped CIs",
    x = "Number of Grid Cells (log scale)",
    y = "Number of Plant Species (log scale)"
  ) +
  theme_minimal()

print(p)

dev.off()

saveRDS(list(plant_df = plant_df, ci_df = ci_df, y_pred = y_pred, y_predUpper = y_predUpper),
        file.path(plot_data_dir, "plotdata_plant_extent_model.rds"))

fit_predUpper <- data.frame(
  num_grid_cells = x_pred,
  fitted = predict(lognorm_Upper_fit, newdata = data.frame(x_vals = x_pred))
)

fit_predLower <- data.frame(
  num_grid_cells = x_pred,
  fitted = predict(lognorm_lower_fit, newdata = data.frame(x_vals = x_pred))
)

#######CALCULATE (animal) gamma richness ranges based on species richness scaling########
#NOTE: We already calculated the best estimate above

print("Calculating gamma richness of animal based on species richness scaling, with error ranges.")

plant_fit_max <- max(fit_predUpper$fitted)
colnames(fit_predUpper)[2]<-"n"

if ( !.amani_estimateBounds){
  upperEstSpRich<-9999999999
  
}else
{
###allows specifying a custom scaling factor (mostly for troubleshooting)
if ( is.null( customScale)) {
  
  extent_sf_here <- (animal_max / plant_fit_max)
  
  if (spatial == "explicit") {
    # spatially explicit simulation (uses empirical / distance-weighted structure)
    upperEstSpRich <- simulateScaledSpeciesRichness(
      fit_predUpper,
      extent_scaling_factor = extent_sf_here,
      sampling_weights = sampling_weights
    )
    
  } else {
    
    # non-spatial / neutral simulation (uniform sampling weights)
    
    upperEstSpRich <- simulateScaledSpeciesRichness(
      fit_predUpper,
      extent_scaling_factor = extent_sf_here,
      sampling_weights = sampling_weights_neutral
    )
    
  }
  
  
  
} else {
  
  extent_sf_here <-customScale
  upperEstSpRich <-999999999
  
}
}





plant_fit_max <- max(fit_predLower$fitted)
###allows specifying a custom scaling factor (mostly for troubleshooting)


if ( !.amani_estimateBounds){
  lowerEstSpRich<-0
  
}else
{

if ( is.null( customScale)) {
  
  extent_sf_here <- (animal_max / plant_fit_max)
  colnames(fit_predLower)[2]<-"n"
  
  if (spatial == "explicit") {
    # spatially explicit simulation (uses empirical / distance-weighted structure)
    lowerEstSpRich <- simulateScaledSpeciesRichness(
      fit_predLower,
      extent_scaling_factor = extent_sf_here,
      sampling_weights = sampling_weights
    )
    
  } else {
    
    # non-spatial / neutral simulation (uniform sampling weights)
    
    lowerEstSpRich <- simulateScaledSpeciesRichness(
      fit_predLower,
      extent_scaling_factor = extent_sf_here,
      sampling_weights = sampling_weights_neutral
    )
    
  }
  
  
} else {
  
  extent_sf_here <-customScale
  lowerEstSpRich<-0
  
}}


output_text <- paste0(
  "Gamma Richness Estimate (spatial extent scaling)\n",
  "================================================\n\n",
  "Fit parameters:    [", paste(fit_params, collapse = ", "), "]\n",
  "Extent scaling factor: ", round(animal_max / max(fit_pred$n), 6), "\n",
  "delta_meanlog:     ", round(.amani_delta_meanlog, 6),
    ifelse("meanlog" %in% fit_params, " (active)", " (off)"), "\n",
  "ratio_sdlog:       ", round(.amani_ratio_sdlog, 6),
    ifelse("sdlog" %in% fit_params, " (active)", " (off)"), "\n",
  "Best estimate:     ", gammaRichnessEst, " species\n",
  "Upper bound:       ", upperEstSpRich, " species\n",
  "Lower bound:       ", lowerEstSpRich, " species\n\n",
  "Description: Projected total animal species richness of the Amazon ",
  "based on a log-normal distribution of plant spatial extents rescaled ",
  "to match observed animal data. Bounds are based on bootstrap resampling ",
  "at 95% and 105% of the original dataset size.\n\n",
  "Shape-shifting: When fit_params includes 'meanlog' or 'sdlog', the ",
  "empirical plant range-size distribution is warped in log-space to match ",
  "the animal's MDD spatial extent shape before amplitude scaling."
)

timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
output_path <- file.path(output_folder, paste0("gamma_richnessfromSPRICH.", timestamp, ".txt"))
writeLines(output_text, con = output_path)

#### Optional: filter plant cells by Chao1 sampling completeness ########
# When .amani_filter_plant_completeness is set, only plant grid cells with
# Chao1 completeness >= that threshold are used in the correlation analyses
# below. This reduces sample size but strengthens the signal by removing
# poorly-sampled plant cells. The filter applies ONLY to the correlations,
# not to the rest of the pipeline.

#NOTE: Should this be incorporated earlier in the workflow? 


filter_plant_comp <- if (exists(".amani_filter_plant_completeness", envir = .GlobalEnv)) {
  .amani_filter_plant_completeness
} else {
  NULL
}

if (!is.null(filter_plant_comp) && is.numeric(filter_plant_comp)) {
  message("Computing plant Chao1 completeness for correlation filtering (threshold = ",
          filter_plant_comp, ")...")

  plant_chao1_stats <- species_counts |>
    sf::st_drop_geometry() |>
    dplyr::group_by(grid_id) |>
    dplyr::summarise(
      plant_n_occ     = sum(count),
      plant_n_species = dplyr::n_distinct(scientificName),
      f1 = sum(count == 1),
      f2 = sum(count == 2),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      chao1 = plant_n_species + dplyr::if_else(
        f2 > 0, f1^2 / (2 * f2), f1 * (f1 - 1) / 2
      ),
      plant_completeness = plant_n_species / chao1
    )

  well_sampled_plant_ids <- plant_chao1_stats$grid_id[
    plant_chao1_stats$plant_completeness >= filter_plant_comp
  ]
  message("  Plant cells passing filter: ", length(well_sampled_plant_ids),
          " / ", nrow(plant_chao1_stats))
} else {
  well_sampled_plant_ids <- NULL
}

#### plotting animal endemism against plant endemism in the area of overlap###############
#### This is important be because the correlation between the animal and plant data is one of the assumptions of this analysis. The stronger those correlations, the more meaningful the result is.
plant_richness_grid_narrow
animal_endemism_filtered

# Ensure both layers are in the same CRS
plant_richness_grid_narrow <- st_transform(plant_richness_grid_narrow, crs = st_crs(animal_endemism_filtered))

# Filter to plant grid cells that overlap with any cockroach points
plant_richness_grid_narrow_filtered <- st_filter(plant_richness_grid_narrow, animal_endemism_filtered)

# Optionally restrict to well-sampled plant cells (with safety guard)
plant_filter_used_endemism <- FALSE
if (!is.null(well_sampled_plant_ids)) {
  candidate <- plant_richness_grid_narrow_filtered |>
    dplyr::filter(grid_id %in% well_sampled_plant_ids)
  if (nrow(candidate) < 15) {
    warning("  Plant completeness filter would leave only ", nrow(candidate),
            " cells for endemism correlation (< 15). Skipping plant filter.",
            call. = FALSE)
  } else {
    plant_richness_grid_narrow_filtered <- candidate
    plant_filter_used_endemism <- TRUE
    message("  Endemism correlation: ", nrow(plant_richness_grid_narrow_filtered),
            " cells after plant completeness filter")
  }
}

# Perform spatial join by ID
endemism_joined <- st_join(plant_richness_grid_narrow_filtered, animal_endemism_filtered %>% dplyr::select(-grid_id))
endemism_joined <- rename_grid_id(endemism_joined)

# Drop incomplete cases (cells where either variable is NA)
endemism_joined <- endemism_joined %>%
  filter(!is.na(raw_endemism), !is.na(animal_endemism_count))

endemism_model <- lm(animal_endemism_count ~ raw_endemism, data = endemism_joined)

timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")

outfile <- paste0(output_folder,paste0("/end_model_validation_", timestamp, ".txt"))

sink(outfile)
cat("Plant completeness filter: ", ifelse(!plant_filter_used_endemism, "OFF",
    paste0(">= ", filter_plant_comp)), "\n\n")
print(summary(endemism_model))
sink()


#### plotting animal richness against plant richness###############
#### This is important be because the correlation between the animal and plant data is one of the assumptions of this analysis. The stronger those correlations, the more meaningful the result is.

animal_richness_sf <- animal_grid %>%
  group_by(grid_id) %>%
  summarise(
    morphotype_richness = n_distinct(Morphotype),
    geometry = geometry[1]  # retain one geometry per grid_id
  ) %>%
  sf::st_as_sf()

# Ensure both layers are in the same CRS
plant_richness_grid <- st_transform(plant_richness_grid, crs = st_crs(animal_grid))

# Filter to plant grid cells that overlap with any cockroach points
plant_filtered <- st_filter(plant_richness_grid, animal_grid)

# Optionally restrict to well-sampled plant cells (with safety guard)
plant_filter_used_richness <- FALSE
if (!is.null(well_sampled_plant_ids)) {
  candidate <- plant_filtered |>
    dplyr::filter(grid_id %in% well_sampled_plant_ids)
  if (nrow(candidate) < 15) {
    warning("  Plant completeness filter would leave only ", nrow(candidate),
            " cells for richness correlation (< 15). Skipping plant filter.",
            call. = FALSE)
  } else {
    plant_filtered <- candidate
    plant_filter_used_richness <- TRUE
    message("  Richness correlation: ", nrow(plant_filtered),
            " cells after plant completeness filter")
  }
}

#plot overlap of data
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
output_filename <- paste0("animalPointsOnplantOverlap.", timestamp, ".pdf")
output_path <- file.path(output_folder, output_filename)

# Open PDF device
pdf(output_path, width = 8, height = 6)
p<-ggplot() +
  geom_sf(data = plant_richness_grid, fill = "grey90", color = "black") +
  geom_sf(data = plant_filtered, fill = "green", alpha = 0.5) +
  geom_sf(data = animal_grid, color = "red", size = 0.5) +
  labs(title = "Filtered Plant Richness Grid Overlapping Animal Points")
print(p)
dev.off()

# Perform spatial join by ID
richness_joined <- st_join(plant_filtered, animal_richness_sf %>% dplyr::select(-grid_id))
richness_joined <- rename_grid_id(richness_joined)

# Drop incomplete cases
richness_joined <- richness_joined %>%
  filter(!is.na(plant_richness), !is.na(morphotype_richness))

richness_model <- lm(morphotype_richness ~ plant_richness, data = richness_joined)

timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")

outfile <- paste0(output_folder,paste0("/rich_model_validation_", timestamp, ".txt"))

sink(outfile)
cat("Plant completeness filter: ", ifelse(!plant_filter_used_richness, "OFF",
    paste0(">= ", filter_plant_comp)), "\n\n")
print(summary(richness_model))
sink()


print("Plotting correlation among raw endemism")

# Plot raw endemism correlation
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
output_filename <- paste0("plant_vs_insect.RichnessCorr.", timestamp, ".pdf")
output_path <- file.path(output_folder, output_filename)

# Open PDF device
pdf(output_path, width = 8, height = 6)



# Extract R² value
r2 <- summary(richness_model)$r.squared
r2_text <- sprintf("R² = %.3f", r2)  # Format to 3 decimal places

p<-ggplot(richness_joined, aes(x = plant_richness, y = morphotype_richness)) +
  geom_point(color = "darkgreen", size = 3, alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  theme_minimal() +
  labs(
    title = "Plant vs. Animal Richness per Grid Cell",
    x = "Plant Richness",
    y = "Animal Morphotype Richness",
    caption = paste("Each point represents a grid cell —", r2_text)
  )
print(p)
# Close the device
dev.off()

saveRDS(sf::st_drop_geometry(richness_joined)[, c("plant_richness", "morphotype_richness")],
        file.path(plot_data_dir, "plotdata_richness_correlation.rds"))

# ── Correlation summary report ───────────────────────────────────────
# Write a single file summarizing the cell counts, filter settings,
# and correlation statistics for both models. This serves as a quick
# reference for assessing model reliability.
# ─────────────────────────────────────────────────────────────────────

endemism_summ <- summary(endemism_model)
richness_summ <- summary(richness_model)

n_animal_cells <- length(unique(animal_grid$grid_id))
n_endemism_cells <- nrow(endemism_joined)
n_richness_cells <- nrow(richness_joined)

correlation_report <- paste0(
  "=======================================================\n",
  " Plant–Animal Correlation Summary\n",
  " Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n",
  "=======================================================\n\n",
  "--- Filter settings ---\n",
  "Animal cell filter:            ", ifelse(n_animal_cells < nrow(grid_sf), "ON", "OFF"), "\n",
  "  Animal cells used:           ", n_animal_cells, "\n",
  "Plant completeness filter:     ",
    ifelse(is.null(filter_plant_comp), "OFF", paste0(">= ", filter_plant_comp)), "\n",
  "  Applied to endemism corr:    ", ifelse(plant_filter_used_endemism, "YES", "NO (< 15 cells or OFF)"), "\n",
  "  Applied to richness corr:    ", ifelse(plant_filter_used_richness, "YES", "NO (< 15 cells or OFF)"), "\n\n",
  "--- Endemism correlation (animal endemism ~ plant raw endemism) ---\n",
  "  Cells:       ", n_endemism_cells, "\n",
  "  R²:          ", round(endemism_summ$r.squared, 4), "\n",
  "  Adj. R²:     ", round(endemism_summ$adj.r.squared, 4), "\n",
  "  F-statistic: ", round(endemism_summ$fstatistic[1], 2),
  " (df = ", endemism_summ$fstatistic[2], ", ", endemism_summ$fstatistic[3], ")\n",
  "  p-value:     ", format.pval(pf(endemism_summ$fstatistic[1],
                                     endemism_summ$fstatistic[2],
                                     endemism_summ$fstatistic[3],
                                     lower.tail = FALSE), digits = 4), "\n",
  "  Slope:       ", round(coef(endemism_model)[2], 4),
  " (SE = ", round(endemism_summ$coefficients[2, 2], 4), ")\n\n",
  "--- Richness correlation (animal richness ~ plant richness) ---\n",
  "  Cells:       ", n_richness_cells, "\n",
  "  R²:          ", round(richness_summ$r.squared, 4), "\n",
  "  Adj. R²:     ", round(richness_summ$adj.r.squared, 4), "\n",
  "  F-statistic: ", round(richness_summ$fstatistic[1], 2),
  " (df = ", richness_summ$fstatistic[2], ", ", richness_summ$fstatistic[3], ")\n",
  "  p-value:     ", format.pval(pf(richness_summ$fstatistic[1],
                                     richness_summ$fstatistic[2],
                                     richness_summ$fstatistic[3],
                                     lower.tail = FALSE), digits = 4), "\n",
  "  Slope:       ", round(coef(richness_model)[2], 4),
  " (SE = ", round(richness_summ$coefficients[2, 2], 4), ")\n\n",
  "=======================================================\n"
)

# Print to console
cat(correlation_report)

# Write to file
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
corr_report_path <- file.path(output_folder, paste0("correlation_summary_", timestamp, ".txt"))
writeLines(correlation_report, con = corr_report_path)
message("Correlation summary written to: ", corr_report_path)



# ------------------------------------------------------------
# plot simulated species richness per grid cell
# this uses the exported simulated occurrence sf object
# ------------------------------------------------------------

# drop geometry to avoid unnecessary spatial overhead during counting
# and count distinct simulated species per grid cell
sim_richness_by_grid <- simulated_occurrence_sf_TOEXPORT |>
  sf::st_drop_geometry() |>                     # remove geometry for fast tabulation
  dplyr::group_by(grid_id) |>                   # group by grid cell
  dplyr::summarise(
    sim_richness = dplyr::n_distinct(scientificName),  # richness = unique simulated species
    .groups = "drop"
  )

# join the simulated richness back to the grid geometry
# this ensures every grid cell has a polygon for plotting
sim_richness_grid_sf <- grid_sf |>
  dplyr::left_join(
    sim_richness_by_grid,
    by = "grid_id"                              # spatial key
  )

# replace NA richness values (cells with no simulated species) with zero
sim_richness_grid_sf$sim_richness[
  is.na(sim_richness_grid_sf$sim_richness)
] <- 0

# quick sanity check on the richness distribution
summary(sim_richness_grid_sf$sim_richness)

# ------------------------------------------------------------
# plot the simulated richness map
# ------------------------------------------------------------

p_sim_richness <- ggplot2::ggplot() +
  ggplot2::geom_sf(
    data = sim_richness_grid_sf,                # grid polygons with simulated richness
    ggplot2::aes(fill = sim_richness),           # color by richness
    color = NA                                  # remove borders for cleaner map
  ) +
  ggplot2::scale_fill_viridis_c(
    option = "magma",                           # perceptually uniform palette
    direction = -1,
    name = "Simulated richness",
    na.value = "white"
  ) +
  ggplot2::theme_minimal() +
  ggplot2::labs(
    title = "Simulated Species Richness per Grid Cell",
    caption = "Richness computed from simulated species occurrences"
  )

# Bare print commented out — crashes RStudio when Plots pane is small.
# The plot is saved to PDF below.
# print(p_sim_richness)

# ------------------------------------------------------------
#save to file
# ------------------------------------------------------------

timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
out_file <- file.path(
  output_folder,
  paste0("simulated_richness_map_", timestamp, ".pdf")
)
pdf(out_file, width = 8, height = 6)
print(p_sim_richness)
dev.off()

saveRDS(sim_richness_grid_sf[, c("grid_id", "sim_richness", "geometry")],
        file.path(plot_data_dir, "plotdata_simulated_richness_map.rds"))

print("Done.")

