print("Starting raw endemism calculator script...")

# Working directory should already be set to project root (2-ModellingAnimals/)
# by the calling script (01_CockroachAnalysis.R, 01_GBIFAnalysis.R, etc.)

SEEDnumber<-2024299

###START: LOAD LIBRARIES####
#library(progress)
library(data.table)
library(sf)
#library(lwgeom)
library(ggplot2)
#library(terra)
#library(tidyr)
#library(readr)
#library(stringr)
#library(parzer)
#library(MASS)
library(dplyr)
#library(readxl)
#library(vegan)
#library(scales)
#library(glue)
#library(purrr)


#this function will be useful to rename spatial opject grid id headers, because they get changed a lot
rename_grid_id <- function(df) {
  if ("grid_id.x" %in% colnames(df)) {
    df <- dplyr::rename(df, grid_id = grid_id.x)
  }
  df
}



print("...setting up some variables and importing data...")
###Output results folder####
if (exists(".amani_output_dir")) {
  output_folder <- .amani_output_dir
} else {
  output_folder <- "output/analysis_results"
}

# Create plot_data/ subdirectory for saving plot data as RDS
plot_data_dir <- file.path(output_folder, "plot_data")
dir.create(plot_data_dir, showWarnings = FALSE, recursive = TRUE)

################Plant import and setup###################
###step 1 make the grid

dataDir<-paste0(dirname(getwd()),"/1-FilteringPlantData/output")

species_counts<-occ_sf#readRDS(paste0(dataDir,"/species_counts4_filtered.rds"))
species_counts<-rename_grid_id(species_counts)

#endemismData<-#readRDS(paste0(dataDir,"/FINAL_allEndemism_FullyFiltered.rds"))
endemismData<-rename_grid_id(endemismData)

# Extract the grid from the endemism data to preserve Stage 0 grid_ids.
# Do NOT rebuild via st_make_grid() — different bounding box produces
# different grid_ids that break downstream joins with species_counts.
grid_size_m <- 22360
grid_sf <- endemismData

# Load Amazon boundary for clipping plots
amazon_shp <- st_read(
  paste0(dirname(getwd()), "/0-CuratingPlantData/data/amazon_shapefile/amazon_sensulatissimo_gmm_v1.shp"),
  quiet = TRUE
) |> st_transform(st_crs(grid_sf))

# Clip grid to Amazon boundary for plotting
grid_sf_clipped <- st_intersection(grid_sf, amazon_shp)

# Save plant endemism plot data for cross-organism compilation
saveRDS(grid_sf_clipped[, c("grid_id", "raw_endemism", "geometry")],
        file.path(plot_data_dir, "plotdata_plant_raw_endemism.rds"))
saveRDS(grid_sf_clipped[, c("grid_id", "weighted_endemism", "geometry")],
        file.path(plot_data_dir, "plotdata_plant_weighted_endemism.rds"))

# Bare plots commented out to avoid crashing RStudio when Plots pane is small.
# These are saved to PDF downstream anyway.
# ggplot() +
#   geom_sf(data = grid_sf_clipped, aes(fill = raw_endemism), color = NA) +
#   scale_fill_viridis_c(option = "plasma", na.value = "white") +
#   theme_minimal() +
#   labs(title = "Species-Level Endemism (500 sq km grid)", fill = "Raw Endemism")
#
# ggplot() +
#   geom_sf(data = grid_sf_clipped, aes(fill = weighted_endemism), color = NA) +
#   scale_fill_viridis_c(option = "plasma", na.value = "white") +
#   theme_minimal() +
#   labs(title = "Species-Level Endemism (500 sq km grid)", fill = "Weighted Endemism")

########Plot and output animal endemism and richness in area of overlap##########
print("Plotting and outputting animal endemism and richness in area of overlap...")

timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
output_filename <- paste0("animal_richness_map_", timestamp, ".pdf")
output_path <- file.path(output_folder, output_filename)

# Open PDF device
pdf(output_path, width = 8, height = 6)

# Clip to Amazon boundary before plotting
animal_grid_data_clipped <- st_intersection(
  animal_grid_data %>% filter(animal_richness > 0),
  amazon_shp
)

p<-ggplot(data = animal_grid_data_clipped) +
  geom_sf(aes(fill = animal_richness), color = "black", size = 0.1) +
  scale_fill_viridis_c(option = "magma", direction = -1, name = "Species Richness") +
  theme_minimal() +
  theme(legend.position = "right") +
  labs(title = "Animal Richness (Amazon basin)",
       caption = "Data: Animal occurrences mapped to Amazon grid cells")
print(p)

# Close the device
dev.off()

saveRDS(animal_grid_data_clipped[, c("grid_id", "animal_richness", "geometry")],
        file.path(plot_data_dir, "plotdata_animal_richness_map.rds"))

#Now do it for endemism
#Step 1: Compute Endemism per Grid Cell
#
animal_grid<-rename_grid_id(animal_grid)

animal_endemism <- animal_grid %>%
  group_by(Morphotype) %>%  
  mutate(total_occurrences = n()) %>%  # Count total occurrences per species
  ungroup() %>%
  filter(total_occurrences == 1) %>%  # Keep only species that appear in ONE grid
  group_by(grid_id) %>%
  summarize(animal_endemism_count = n(), .groups = "drop")  # Count unique species per grid cell


animal_endemism<-rename_grid_id(animal_endemism)
grid_sf<-rename_grid_id(grid_sf)
#Step 2: Merge Endemism Data with Grid
animal_endemism_grid <- st_join(grid_sf %>% dplyr::select(grid_id, geometry),
                               animal_endemism %>% dplyr::select(-grid_id),
                               join = st_intersects) %>%
  mutate(animal_endemism_count = replace_na(animal_endemism_count, 0))


timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
output_filename <- paste0("animal_endemism_map_", timestamp, ".pdf")
output_path <- file.path(output_folder, output_filename)

# Open PDF device
pdf(output_path, width = 8, height = 6)

#Step 3: Plot Endemism Map (clipped to Amazon)
animal_endemism_clipped <- st_intersection(
  animal_endemism_grid %>% filter(animal_endemism_count > 0),
  amazon_shp
)
p<-ggplot(data = animal_endemism_clipped) +
  geom_sf(aes(fill = animal_endemism_count), color = "black", size = 0.1) +
  scale_fill_viridis_c(option = "plasma", direction = -1, name = "Endemism") +
  theme_minimal() +
  theme(legend.position = "right") +
  labs(title = "Animal Endemism (Amazon basin)",
       caption = "Data: Single-cell endemic species per grid cell")
print(p)
# Close the device
dev.off()

saveRDS(animal_endemism_clipped[, c("grid_id", "animal_endemism_count", "geometry")],
        file.path(plot_data_dir, "plotdata_animal_endemism_map.rds"))

#########plot plant endemism in the area of overlap#########
###now get the plant map
print("Plotting plant endemism and richness in the same area as the animals...")
#Step 1: Filter the Grid to Match Animal Data

#animal_endemism_grid <- animal_endemism_grid %>%
#  rename(grid_id.x = grid_id.x.x) %>%
#  dplyr::select(-grid_id.x.y)  # Remove duplicate if unnecessary

animal_endemism_grid<-rename_grid_id(animal_endemism_grid)

animal_endemism_filtered <- animal_endemism_grid %>%
  filter(animal_endemism_count > 0)

# Check whether there are enough single-cell endemics overlapping with
# plant endemism data to run the correlation and scaling analysis.
n_endemic_cells <- nrow(animal_endemism_filtered)

if (n_endemic_cells > 0) {
  # Check whether any plant occurrence data exists in the animal overlap area.
  # Endemism will be recomputed locally, so we just need plant records present.
  overlap_bbox <- st_bbox(animal_endemism_filtered)
  overlap_test_ids <- grid_sf$grid_id[
    as.logical(st_intersects(grid_sf, st_as_sfc(overlap_bbox), sparse = FALSE))
  ]
  n_plant_in_overlap <- sum(
    sf::st_drop_geometry(species_counts)$grid_id %in% overlap_test_ids
  )
  n_overlap <- if (n_plant_in_overlap > 0) n_endemic_cells else 0
} else {
  n_overlap <- 0
}

can_run_endemism <- n_endemic_cells > 0 && n_overlap >= 2

if (!can_run_endemism) {
  warning("Insufficient endemism data for correlation analysis ",
          "(", n_endemic_cells, " endemic cells, ", n_overlap, " with plant overlap). ",
          "Skipping endemism-plant scaling. ",
          "Using occupied grid cells for spatial extent instead.")

  animal_endemism_filtered <- animal_endemism_grid %>%
    filter(grid_id %in% unique(animal_grid$grid_id))

  overlap_bbox <- st_bbox(animal_endemism_filtered)

  # Set scaled_results so ModelAnimal.R conditional logic works
  scaled_results <- NA

} else {

overlap_bbox <- st_bbox(animal_endemism_filtered)  # Get spatial extent of animal data

# Recompute plant single-cell endemism using ONLY the overlap region.
# A plant species is "endemic" here if it appears in exactly one grid cell
# within the animal data extent, regardless of its Amazon-wide range.
# This gives a much more meaningful comparison for localised animal datasets.

# Use grid_id filtering rather than st_crop to avoid sf/data.table class issues.
# Identify which grid cells fall inside the overlap bbox.
overlap_grid_ids <- grid_sf$grid_id[
  as.logical(st_intersects(grid_sf, st_as_sfc(overlap_bbox), sparse = FALSE))
]

local_plant_endemism <- species_counts %>%
  sf::st_drop_geometry() %>%
  filter(grid_id %in% overlap_grid_ids) %>%
  group_by(scientificName) %>%
  filter(n_distinct(grid_id) == 1) %>%     # single-cell endemics in overlap
  ungroup() %>%
  group_by(grid_id) %>%
  summarise(raw_endemism_local = n_distinct(scientificName), .groups = "drop")

# Join local endemism back onto the grid geometry in the overlap
plant_endemism_grid <- grid_sf %>%
  st_crop(overlap_bbox) %>%
  left_join(local_plant_endemism, by = "grid_id") %>%
  filter(!is.na(raw_endemism_local), raw_endemism_local > 0) %>%
  dplyr::select(grid_id, geometry, raw_endemism = raw_endemism_local)

message("  Local plant endemism computed over ", nrow(local_plant_endemism),
        " cells in overlap region (", sum(local_plant_endemism$raw_endemism_local),
        " single-cell endemic plant species).")

#Step 4: Create the plant plot
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
output_filename <- paste0("plant_endemism_map_overlap", timestamp, ".pdf")
output_path <- file.path(output_folder, output_filename)

# Open PDF device
pdf(output_path, width = 8, height = 6)
# Plant Endemism Plot (Now Properly Cropped)
p<-ggplot(data = plant_endemism_grid) +
  geom_sf(aes(fill = raw_endemism), color = "black", size = 0.1) +
  scale_fill_viridis_c(option = "magma", direction = -1, name = "Plant Endemism\n(local)") +
  theme_minimal() +
  theme(legend.position = "right") +
  labs(title = "Plant Endemism in Area of Animal Overlap",
       caption = "Endemism recomputed within overlap region only")

print(p)
# Close the device
dev.off()


# Spatial join to align plant and cockroach endemism by grid cell
combined_endemism <- st_join(animal_endemism_filtered,
                             plant_endemism_grid %>% dplyr::select(-grid_id),
                             join = st_equals)

# Note: grid cell filtering by min_occurrences is now applied early in
# ModelAnimal.R (to animal_grid and animal_grid_data) so all downstream
# plots and analyses use only well-sampled cells.

# Drop rows where either variable is NA — these are animal-endemic cells
# without single-cell endemic plants (raw_endemism is NA when a cell has
# zero single-cell endemic plant species in the upstream data).
combined_endemism_complete <- combined_endemism %>%
  filter(!is.na(raw_endemism), !is.na(animal_endemism_count))

n_complete <- nrow(combined_endemism_complete)
n_total    <- nrow(combined_endemism)

if (n_complete < n_total) {
  message("  Endemism correlation: ", n_total - n_complete, " of ", n_total,
          " animal-endemic cells dropped (no single-cell endemic plants).")
}

#########get correlation between cockroach and plant endemism#########

if (n_complete < 3) {
  warning("Only ", n_complete, " cells with both animal and plant single-cell endemics. ",
          "Endemism scaling requires >= 3 complete cases. Skipping.",
          call. = FALSE)

  # Fallback global outputs so downstream scripts don't crash
  endemism_r2 <<- NA_real_
  scaled_results <- rep(NA_real_, 1000)
  scaled_endemism_median <<- NA_real_
  scaled_endemism_ci <<- c(NA_real_, NA_real_)

  # Write empty/NA plot-data RDS files so summarize_organism() has something to read
  saveRDS(sf::st_drop_geometry(combined_endemism_complete)[, c("raw_endemism", "animal_endemism_count")],
          file.path(plot_data_dir, "plotdata_endemism_correlation.rds"))
  saveRDS(data.frame(endemism_ratio = numeric(0)),
          file.path(plot_data_dir, "plotdata_endemism_ratios.rds"))
  saveRDS(dplyr::mutate(grid_sf_clipped[, c("grid_id", "geometry")],
                        scaled_plant_endemism = NA_real_),
          file.path(plot_data_dir, "plotdata_scaled_endemism_map.rds"))
  saveRDS(list(results = scaled_results, med = NA_real_, mn = NA_real_,
               ci_low = NA_real_, ci_high = NA_real_),
          file.path(plot_data_dir, "plotdata_simulated_endemism_dist.rds"))

} else {

# Plot raw endemism correlation
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
output_filename <- paste0("plant_vs_insect.corr", timestamp, ".pdf")
output_path <- file.path(output_folder, output_filename)

# Open PDF device
pdf(output_path, width = 8, height = 6)


# Filter to positive values for log-log regression
combined_endemism_log <- combined_endemism_complete %>%
  filter(raw_endemism > 0, animal_endemism_count > 0)
n_log <- nrow(combined_endemism_log)

# Fit log-log linear model (power-law relationship)
lm_model <- lm(log(animal_endemism_count) ~ log(raw_endemism), data = combined_endemism_log)

# Extract R² value
r2 <- summary(lm_model)$r.squared
endemism_r2 <<- r2
r2_text <- sprintf("R² = %.3f  (n = %d)", r2, n_log)
print("Outputting the log-log correlation between the two organisms' raw endemism...")
# Plot on log-log scale
p<-ggplot(combined_endemism_log, aes(x = raw_endemism, y = animal_endemism_count)) +
  geom_point(alpha = 0.7, color = "darkgreen") +
  geom_smooth(method = "lm", se = FALSE, color = "red") +
  scale_x_log10() +
  scale_y_log10() +
  theme_minimal() +
  labs(
    x = "Plant Raw Endemism (log scale)",
    y = "Animal Raw Endemism (log scale)",
    title = "Correlation Between Plant and Animal Endemism (log-log)",
    caption = paste("Each point represents a grid cell —", r2_text)
  )
print(p)

# Close the device
dev.off()

saveRDS(sf::st_drop_geometry(combined_endemism_log)[, c("raw_endemism", "animal_endemism_count")],
        file.path(plot_data_dir, "plotdata_endemism_correlation.rds"))

##########Scale plant endemism to match cockroach endemism########

#Now, use the ratio of animal endemism to plant endemism to generate a distribution of endemism ratios
#The endemism ratios will be used to scale the plant endemism data appropriately for the animal
#we are going to use a log-normal distribution because that is good for ratios
#I also tried with a gamma distribution, but then it seems like most of the values are basically 0, which I think is not realistic

# Step 1: Clean and compute the ratio per cell
# 
print("Doing a rescaling of plant data to match animal parameters...")

ratio_data <- combined_endemism %>%
  filter(animal_endemism_count > 0, raw_endemism > 0) %>%
  mutate(endo_ratio = animal_endemism_count / raw_endemism)

#ratio_data$endo_ratio are the ratios of roach end to plant endemism

# Step 2: Get log-transformed mean and sd
log_ratios <- log(ratio_data$endo_ratio)
log_mean <- mean(log_ratios)
log_sd   <- sd(log_ratios)

# Step 3: Sample endemism ratios from a log-normal distribution
set.seed(SEEDnumber-1)
log_endemism_ratios <- rlnorm(1000, meanlog = log_mean, sdlog = log_sd)

# Convert back to real values
endemism_ratios <- exp(log_endemism_ratios)
# endemism_ratios are a simulated distribution of ratios between animal endemism and plant endemism



# Optional: convert to a dataframe for plotting
endemism_ratio_df <- data.frame(endemism_ratio = endemism_ratios)

endemism_ratio_df_under2 <- endemism_ratio_df %>%
  filter(endemism_ratio < 2)

print("Plotting distribution of endemism ratios...")
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
output_filename <- paste0("distEndemismRatios.", timestamp, ".pdf")
output_path <- file.path(output_folder, output_filename)

# Open PDF device
pdf(output_path, width = 8, height = 6)
p<-ggplot(endemism_ratio_df_under2, aes(x = endemism_ratio)) +
  geom_histogram(bins = 30, fill = "darkgreen", color = "black") +
  theme_minimal() +
  labs(
    title = "Endemism Ratios < 2 (Animal / Plant)",
    subtitle = "LogNorm-distributed values sampled from observed cell ratios",
    x = "Endemism Ratio",
    y = "Frequency"
  )

print(p)
# Close the device
dev.off()

saveRDS(endemism_ratio_df_under2, file.path(plot_data_dir, "plotdata_endemism_ratios.rds"))

print("Getting total Amazonian endemism based on raw endemism per grid cell...")
#########Get total Amazonian animal endemism#########
######now use the median value to recalculate the world grid for animals
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
output_filename <- paste0("animalEnd.Amaz.medianLogNorm.", timestamp, ".pdf")
output_path <- file.path(output_folder, output_filename)

pdf(output_path, width = 8, height = 6)
ratio_median<-median(endemism_ratios)

#ratio_median is the simulated median ratio between animal endemism and plant endemism


# Apply multiplicative scaling to the data based on the median value
grid_sf <- grid_sf %>%
  mutate(scaled_plant_endemism = raw_endemism * ratio_median) #we multiply these because the ratio is animal/plant*plant so the final result is in "animal units"

p<-ggplot() +
  geom_sf(data = grid_sf, aes(fill = scaled_plant_endemism), color = "black") +
  scale_fill_viridis_c(option = "plasma", na.value = "white") +
  theme_minimal() +
  labs(title = "Scaled Plant Endemism (Adjusted by Animal Offset)",
       fill = "Scaled Endemism")

print(p)
# Close the device
dev.off()

saveRDS(grid_sf[!is.na(grid_sf$scaled_plant_endemism), c("grid_id", "scaled_plant_endemism", "geometry")],
        file.path(plot_data_dir, "plotdata_scaled_endemism_map.rds"))

###this value should be equal to the total number of predicted endemic species
totalAmazonianScaledEndemicSpp<-sum(
  grid_sf$scaled_plant_endemism[!is.na(grid_sf$scaled_plant_endemism)]
)

print("Calculating error ranges...")
#########Get error bounds for animal endemism scaling#########
##now i need to use the endemism ratios to simulate new data, multiple times since we are drawing from a distribution, and then find the total number of species


# Simulate species_counts using scaled endemism and sampled abundance values
# Apply the model slope to scale plant endemism
# Step 1: Apply the model slope to scale raw plant endemism

abundance_pool <- species_counts$count
head(species_counts)


###now do it with the randomized endemism ratios to get error bounds
# Initialize a vector to store total richness per run

##step 1. make the above a function
#scalingFactor<-ratio_median###set this for testing


# scaleEndemismValues<-function(endemismRatio){
#   
#   grid_sf <- grid_sf %>%
#     mutate(scaled_plant_endemism = raw_endemism * endemismRatio)  # multiply each cell's raw endemism by the slope
#   #totalAmazonianScaledEndemicSpp<-
#   #  sum(unique(grid_sf$scaled_plant_endemism)[!is.na(unique(grid_sf$scaled_plant_endemism))])
#   #sum(unique(grid_sf$raw_endemism)[!is.na(unique(grid_sf$raw_endemism))])
#   
#   grid_sf <- grid_sf %>%
#     mutate(n_species = round(scaled_plant_endemism))
#  
#    # if(max(grid_sf$n_species)>10000000){ message("Error: too many predicted species....skipping for...") 
#    #  message(paste0(toString(max(grid_sf$n_species)), " max species"))
#    #  numberOfEndemicsInAmazon <-9999999
#    #  }else{
#    #    message(paste0(toString(max(grid_sf$n_species)), " max species"))
#   
#   # Step 2: Simulate species occurrences based on the scaled endemism
#   simulated_species_counts <- grid_sf %>%
#     filter(n_species > 0) %>%
#     rowwise() %>%
#     mutate(
#       species_list = list(paste0(
#         #grid_id.x, 
#         grid_id,
#         "_Plant_species_", str_pad(1:n_species, 12, pad = "0")))
#     ) %>%
#     unnest(species_list) %>%
#     rename(scientificName = species_list) %>%
#     mutate(count = sample(abundance_pool, n(), replace = TRUE)) %>%
#     dplyr::select(
#       #grid_id.x, 
#       grid_id,
#       scientificName, count, geometry)
#   
#   numberOfEndemicsInAmazon<-length(unique(simulated_species_counts$scientificName))
#   
#   return(numberOfEndemicsInAmazon)
# }


scaleEndemismValues2 <- function(endemismRatio) {
  
  n_species <- round(grid_sf$raw_endemism * endemismRatio)
  
  sum(n_species[n_species > 0], na.rm = TRUE)
}


median(endemism_ratios)
min(endemism_ratios)

length(unique(species_counts$scientificName))

scaled_results <- lapply(endemism_ratios, scaleEndemismValues2)
scaled_results<-unlist(scaled_results)

#scaleEndemismValues(ratio_median)

# Calculate statistics
med <- median(scaled_results)
mn <- mean(scaled_results)
ci_low <- quantile(scaled_results, 0.05)
ci_high <- quantile(scaled_results, 0.95)

# Expose to global env for model_summary.rds
scaled_endemism_median <<- med
scaled_endemism_ci <<- c(ci_low, ci_high)


print("Outputting final raw endemism calculation...")
#########plot the output of endemism simulation#########


# Generate timestamped filename and path
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
output_filename <- paste0("dist.Sim.endemSpp.", timestamp, ".pdf")
output_path <- file.path(output_folder, output_filename)

# Create the PDF plot
# Zoom x-axis to bulk of the distribution (99th percentile + 20% padding)
x_upper <- quantile(scaled_results, 0.99, na.rm = TRUE) * 1.2

pdf(output_path, width = 8, height = 6)
p<-ggplot(data.frame(richness = scaled_results), aes(x = richness)) +
  geom_histogram(bins = 30, fill = "skyblue", color = "black") +
  geom_vline(xintercept = med, color = "red", linetype = "dashed", linewidth = 1) +
  geom_vline(xintercept = mn, color = "blue", linetype = "dotted", linewidth = 1) +
  geom_vline(xintercept = ci_low, color = "darkgreen", linetype = "solid", linewidth = 1) +
  geom_vline(xintercept = ci_high, color = "darkgreen", linetype = "solid", linewidth = 1) +
  annotate("text", x = med, y = Inf, label = paste0("Median: ", round(med, 2)), vjust = 5.5, color = "red") +
  annotate("text", x = mn, y = Inf, label = paste0("Mean: ", round(mn, 2)), vjust = 10.5, color = "blue") +
  annotate("text", x = ci_low, y = Inf, label = paste0("5% CI: ", round(ci_low, 2)), vjust = 2.5, color = "darkgreen") +
  annotate("text", x = ci_high, y = Inf, label = paste0("95% CI: ", round(ci_high, 2)), vjust = 3.5, color = "darkgreen") +
  scale_x_continuous(limits = c(0, x_upper)) +
  theme_minimal() +
  labs(
    title = "Simulated total endemic species (endemism-ratio method)",
    subtitle = "Using 1000 random endemism ratios from 95% CI",
    x = "Total simulated endemic species",
    y = "Frequency"
  )
print(p)
dev.off()

saveRDS(list(results = scaled_results, med = med, mn = mn,
             ci_low = ci_low, ci_high = ci_high),
        file.path(plot_data_dir, "plotdata_simulated_endemism_dist.rds"))

# Write key results to a text file
txt_filename <- paste0("endemism_ratio_simulation_results_", timestamp, ".txt")
txt_path <- file.path(output_folder, txt_filename)
writeLines(c(
  "=== Endemism-Ratio Simulation Results ===",
  "",
  "Method: Raw endemism extrapolation using per-cell animal/plant endemism ratios.",
  "This estimate scales plant single-cell endemism by observed animal/plant ratios",
  "and sums across the Amazon grid. It does NOT use the spatial-extent model.",
  "",
  paste0("Median endemism ratio:     ", round(ratio_median, 6)),
  paste0("Number of ratio samples:   ", length(endemism_ratios)),
  "",
  "--- Simulated total endemic species ---",
  paste0("Median:              ", round(med, 2)),
  paste0("Mean:                ", round(mn, 2)),
  paste0("5th percentile:      ", round(ci_low, 2)),
  paste0("95th percentile:     ", round(ci_high, 2)),
  paste0("Point estimate:      ", round(totalAmazonianScaledEndemicSpp, 2)),
  "",
  paste0("Timestamp: ", Sys.time())
), txt_path)
message("  Wrote endemism-ratio simulation results to: ", txt_path)

quantile(scaled_results, probs = c(0.40, 0.60))

} # end if (n_complete >= 3)

} # end if (can_run_endemism)

