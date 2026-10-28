# set working directory to project root (2-ModellingAnimals/)
library(rstudioapi)
if (exists(".amani_wd")) {
  setwd(.amani_wd)
} else {
  project_root <- dirname(dirname(rstudioapi::getActiveDocumentContext()$path))
  setwd(project_root)
}

###START: LOAD LIBRARIES####

library(data.table)
library(sf)
library(utils)
library(lwgeom)
library(ggplot2)
library(terra)
library(tidyr)  # For replace_na
library(readr)
library(tidyr)
library(stringr)
library(parzer)
library(MASS)
library(dplyr)
library(readxl)
library(vegan)
library(scales)
library(glue)
library(purrr)


###Output results folder####
output_folder <- "data/finches"

################Plant import and setup###################
###step 1 make the grid
occ_sf<-readRDS("../0-CuratingPlantData/output/amazon_filtered_data.occ_sf.rds")

# Define grid size in meters (22.36 km ≈ 22360 meters)
grid_size_m <- 22360

grid <- st_make_grid(occ_sf, cellsize = grid_size_m, what = "polygons")

# Convert grid to an sf object
grid_sf <- st_sf(grid_id = 1:length(grid), geometry = grid)

species_counts<-readRDS("../0-CuratingPlantData/output/species_counts3.rds")

species_endemism<-readRDS("../0-CuratingPlantData/output/species_endemism2.rds")

raw_endemism<-readRDS("../0-CuratingPlantData/output/raw_species_endemism2.rds")

# Merge grid with endemism data
grid_sf <- st_join(grid_sf, raw_endemism, by = "grid_id")

# Plot plant endemism ###this takes about 5 minutes and I saved it already (April 8 2025)
#ggplot() +
#  geom_sf(data = grid_sf, aes(fill = raw_endemism), color = "black") +
#  scale_fill_viridis_c(option = "plasma", na.value = "white") +
#  theme_minimal() +
#  labs(title = "Species-Level Endemism (500 sq km grid)", fill = "Endemism")

######COCKROACH DATA IMPORT AND SETUP#######
# Define the file path
file_path <- "data/finches/0011407-250402121839773.csv"

# Load the GBIF CSV data
roach_data <- read_tsv(file_path)
names(roach_data)
# Select and rename the relevant columns
selected_roach_data <- roach_data %>%
  dplyr::select(
    Morphotype = scientificName,       # Or verbatimScientificName if needed
    Event_Date = eventDate,
    verB_Lat = decimalLatitude,
    verb_Long = decimalLongitude,
    Locality = locality,
    Issue = issue                      # Add issue column here
  )

clean_roach_data <- selected_roach_data %>%
  filter(!grepl("TAXON_MATCH_FUZZY", Issue))

selected_roach_data <- clean_roach_data %>%
  mutate(
    decimal_lat = parse_lat(verB_Lat),
    decimal_lon = parse_lon(verb_Long)
  )

# Remove rows with missing coordinates
selected_roach_data <- selected_roach_data %>%
  filter(!is.na(decimal_lat) & !is.na(decimal_lon))

###filter out juvenile-only morphotypes
#adult_morphotypes <- selected_roach_data %>%
#  filter(Life_Stage == "Adult") %>%
#  distinct(Morphotype) %>%
#  pull(Morphotype)  # Extract as a vector

#selected_roach_data2 <- selected_roach_data %>%
#  filter(Morphotype %in% adult_morphotypes)


# Convert cockroach occurrence data into an sf object
roach_sf <- selected_roach_data %>%
  st_as_sf(coords = c("decimal_lon", "decimal_lat"), crs = 4326) %>%  # WGS84
  st_transform(st_crs(grid_sf))  # Match CRS with grid

# Perform a spatial join: Assign each cockroach occurrence to a grid cell
roach_grid <- st_join(roach_sf, grid_sf, join = st_within)

roach_endemism <- roach_grid %>%
  group_by(grid_id.x) %>%
  summarize(cockroach_richness = n_distinct(Morphotype), .groups = "drop")

roach_grid_data <- st_join(grid_sf %>% dplyr::select(grid_id.x, geometry),  # Keep only grid_id.x and geometry
                           roach_endemism %>% dplyr::select(grid_id.x, cockroach_richness),  
                           join = st_intersects) %>%
  mutate(cockroach_richness = replace_na(cockroach_richness, 0))  # Replace NAs with 0


########Plot and output cockroach endemism and richness in madre de dios##########


timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
output_filename <- paste0("cockroach_richness_map_", timestamp, ".pdf")
output_path <- file.path(output_folder, output_filename)

# Open PDF device
pdf(output_path, width = 8, height = 6)

# Plot
ggplot(data = roach_grid_data %>% filter(cockroach_richness > 0)) +
  geom_sf(aes(fill = cockroach_richness), color = "black", size = 0.1) +
  scale_fill_viridis_c(option = "magma", direction = -1, name = "Cockroach Richness") +
  theme_minimal() +
  theme(legend.position = "right") +
  labs(title = "Cockroach Richness in Madre de Dios, Peru",
       caption = "Data: Cockroach occurrences mapped to Amazon grid cells")

# Close the device
dev.off()


#Now do it for endemism
#Step 1: Compute Endemism per Grid Cell
roach_endemism <- roach_grid %>%
  group_by(Morphotype) %>%  
  mutate(total_occurrences = n()) %>%  # Count total occurrences per species
  ungroup() %>%
  filter(total_occurrences == 1) %>%  # Keep only species that appear in ONE grid
  group_by(grid_id.x) %>%
  summarize(cockroach_endemism = n(), .groups = "drop")  # Count unique species per grid cell
#Step 2: Merge Endemism Data with Grid
roach_endemism_grid <- st_join(grid_sf %>% dplyr::select(grid_id.x, geometry), 
                               roach_endemism %>% dplyr::select(grid_id.x, cockroach_endemism),
                               join = st_intersects) %>%
  mutate(cockroach_endemism = replace_na(cockroach_endemism, 0))  # Replace NA with 0


timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
output_filename <- paste0("cockroach_endemism_map_", timestamp, ".pdf")
output_path <- file.path(output_folder, output_filename)

# Open PDF device
pdf(output_path, width = 8, height = 6)

#Step 3: Plot Endemism Map
ggplot(data = roach_endemism_grid %>% filter(cockroach_endemism > 0)) +  # Remove empty cells
  geom_sf(aes(fill = cockroach_endemism), color = "black", size = 0.1) +
  scale_fill_viridis_c(option = "plasma", direction = -1, name = "Cockroach Endemism") +
  theme_minimal() +
  theme(legend.position = "right") +
  labs(title = "Cockroach Endemism in Madre de Dios, Peru",
       caption = "Data: Unique species per Amazon grid cell")

# Close the device
dev.off()



#########plot plant endemism in the madre de dios#########
###now get the plant map

#Step 1: Filter the Grid to Match Cockroach Data

roach_endemism_grid <- roach_endemism_grid %>%
  rename(grid_id.x = grid_id.x.x) %>%
  dplyr::select(-grid_id.x.y)  # Remove duplicate if unnecessary

roach_endemism_filtered <- roach_endemism_grid %>%
  filter(cockroach_endemism > 0)


madre_de_dios_bbox <- st_bbox(roach_endemism_filtered)  # Get spatial extent of cockroach data


plant_endemism_grid <- st_intersection(grid_sf, roach_endemism_filtered) %>%
  filter(raw_endemism > 0) %>%
  dplyr::select(grid_id.x, geometry, raw_endemism)


#Step 4: Create the plant plot
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
output_filename <- paste0("plant_endemism_map_MDDios", timestamp, ".pdf")
output_path <- file.path(output_folder, output_filename)

# Open PDF device
pdf(output_path, width = 8, height = 6)
# Plant Endemism Plot (Now Properly Cropped)
ggplot(data = plant_endemism_grid) +
  geom_sf(aes(fill = raw_endemism), color = "black", size = 0.1) +
  scale_fill_viridis_c(option = "magma", direction = -1, name = "Plant Endemism") +
  theme_minimal() +
  theme(legend.position = "right") +
  labs(title = "Plant Endemism in Madre de Dios, Peru",
      caption = "Data: Raw plant endemism per Amazon grid cell")

# Close the device
dev.off()


# Spatial join to align plant and cockroach endemism by grid cell
combined_endemism <- st_join(roach_endemism_filtered, plant_endemism_grid, join = st_equals)

#########get correlation between cockroach and plant endemism#########

# Plot raw endemism correlation
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
output_filename <- paste0("plant_vs_insect.corr", timestamp, ".pdf")
output_path <- file.path(output_folder, output_filename)

# Open PDF device
pdf(output_path, width = 8, height = 6)


# Fit the linear model
lm_model <- lm(cockroach_endemism ~ raw_endemism, data = combined_endemism)

# Extract R² value
r2 <- summary(lm_model)$r.squared
r2_text <- sprintf("R² = %.3f", r2)  # Format to 3 decimal places

# Plot with R² in caption
ggplot(combined_endemism, aes(x = raw_endemism, y = cockroach_endemism)) +
  geom_point(alpha = 0.7, color = "darkgreen") +
  geom_smooth(method = "lm", se = FALSE, color = "red") +
  theme_minimal() +
  labs(
    x = "Plant Raw Endemism",
    y = "Cockroach Raw Endemism",
    title = "Correlation Between Plant and Cockroach Endemism",
    caption = paste("Each point represents a grid cell —", r2_text)
  )


# Close the device
dev.off()



##########Scale plant endemism to match cockroach endemism########

#Now, use the ratio of cockroach endemism to plant endemism to generate a number of scaling factors
#The scaling factors will be used to scale the plant endemism data appropriately for cockroaches
#we are going to use a log-normal distribution because that is good for ratios
#I also tried with a gamma distribution, but then it seems like most of the values are basically 0, which I think is not realistic

# Step 1: Clean and compute the ratio per cell
ratio_data <- combined_endemism %>%
  filter(cockroach_endemism > 0, raw_endemism > 0) %>%
  mutate(endo_ratio = cockroach_endemism / raw_endemism)

#ratio_data$endo_ratio are the ratios of roach end to plant endemism

# Step 2: Get log-transformed mean and sd
log_ratios <- log(ratio_data$endo_ratio)
log_mean <- mean(log_ratios)
log_sd   <- sd(log_ratios)

# Step 3: Sample scaling factors from a log-normal distribution
set.seed(123)
log_scaling_factors <- rlnorm(1000, meanlog = log_mean, sdlog = log_sd)

# Convert back to real (scaling) values
scaling_factors <- exp(log_scaling_factors)
#scaling_factors are a simulated distribution of ratios between cockroach endemism and plant endemism



# Optional: convert to a dataframe for plotting
scaling_df <- data.frame(scaling_factor = scaling_factors)

scaling_df_under2 <- scaling_df %>%
  filter(scaling_factor < 2)


# Plot raw endemism correlation
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
output_filename <- paste0("distScalingFactors.corr", timestamp, ".pdf")
output_path <- file.path(output_folder, output_filename)

# Open PDF device
pdf(output_path, width = 8, height = 6)
ggplot(scaling_df_under2, aes(x = scaling_factor)) +
  geom_histogram(bins = 30, fill = "darkgreen", color = "black") +
  theme_minimal() +
  labs(
    title = "Scaling Factors < 2",
    subtitle = "Subset of LogNorm-Distributed Values",
    x = "Scaling Factor",
    y = "Frequency"
  )
# Close the device
dev.off()


#########Get total Amazoniian cockroach endemism#########
######now use the median value to recalculate the world grid for cockroaches
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
output_filename <- paste0("roachEnd.Amaz.medianLogNorm.", timestamp, ".pdf")
output_path <- file.path(output_folder, output_filename)

pdf(output_path, width = 8, height = 6)
ratio_median<-median(scaling_factors)

#ratio_median is the simulated mean ratio between cockroach endemism and plant endemism


# Apply multiplicative scaling to the data based on the median value
grid_sf <- grid_sf %>%
  mutate(scaled_plant_endemism = raw_endemism * ratio_median) #we multiply these because the ratio is cockroach/plant*plant so the final result is in "cockroach units"

ggplot() +
  geom_sf(data = grid_sf, aes(fill = scaled_plant_endemism), color = "black") +
  scale_fill_viridis_c(option = "plasma", na.value = "white") +
  theme_minimal() +
  labs(title = "Scaled Plant Endemism (Adjusted by Cockroach Offset)",
       fill = "Scaled Endemism")

# Close the device
dev.off()

###this value should be equal to the total number of predicted endemic species
totalAmazonianScaledEndemicSpp<-sum(
  grid_sf$scaled_plant_endemism[!is.na(grid_sf$scaled_plant_endemism)]
                                    )


#########Get error bounds for cockroach endemism scaling#########
##now i need to use th scaling factors to simulate new data, multiple times since we are drawing the scaling factor from a distribution, and then find the total number of species


# Simulate species_counts using scaled endemism and sampled abundance values
# Apply the model slope to scale plant endemism
# Step 1: Apply the model slope to scale raw plant endemism

abundance_pool <- species_counts$count
head(species_counts)


###now do it with the randomized scaling factors to get error bounds
# Initialize a vector to store total richness per run

##step 1. make the above a function
#scalingFactor<-ratio_median###set this for testing


scaleEndemismValues<-function(scalingFactor){
  
  grid_sf <- grid_sf %>%
    mutate(scaled_plant_endemism = raw_endemism * scalingFactor)  # multiply each cell's raw endemism by the slope
  totalAmazonianScaledEndemicSpp<-
    sum(unique(grid_sf$scaled_plant_endemism)[!is.na(unique(grid_sf$scaled_plant_endemism))])
  sum(unique(grid_sf$raw_endemism)[!is.na(unique(grid_sf$raw_endemism))])

  grid_sf <- grid_sf %>%
    mutate(n_species = round(scaled_plant_endemism))
  
  
  # Step 2: Simulate species occurrences based on the scaled endemism
  simulated_species_counts <- grid_sf %>%
    filter(n_species > 0) %>%
    rowwise() %>%
    mutate(
      species_list = list(paste0(grid_id.x, "_Plant_species_", str_pad(1:n_species, 10, pad = "0")))
    ) %>%
    unnest(species_list) %>%
    rename(scientificName = species_list) %>%
    mutate(count = sample(abundance_pool, n(), replace = TRUE)) %>%
    dplyr::select(grid_id.x, scientificName, count, geometry)
  
  
  numberOfEndemicsInAmazon<-length(unique(simulated_species_counts$scientificName))
  
  return(numberOfEndemicsInAmazon)
}


median(scaling_factors)
min(scaling_factors)

length(unique(species_counts$scientificName))
scaleEndemismValues(1)

scaleEndemismValues(5)
scaled_results <- lapply(scaling_factors, scaleEndemismValues)
scaled_results<-unlist(scaled_results)

scaleEndemismValues(ratio_median)

raw_endemism
# Calculate statistics
med <- median(scaled_results)
mn <- mean(scaled_results)
ci_low <- quantile(scaled_results, 0.05)
ci_high <- quantile(scaled_results, 0.95)

#########plot the output of endemism simulation#########


# Generate timestamped filename and path
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
output_filename <- paste0("dist.Sim.endemSpp.", timestamp, ".pdf")
output_path <- file.path(output_folder, output_filename)

# Create the PDF plot
pdf(output_path, width = 8, height = 6)
ggplot(data.frame(richness = scaled_results), aes(x = richness)) +
  geom_histogram(bins = 30, fill = "skyblue", color = "black") +
  geom_vline(xintercept = med, color = "red", linetype = "dashed", size = 1) +
  geom_vline(xintercept = mn, color = "blue", linetype = "dotted", size = 1) +
  geom_vline(xintercept = ci_low, color = "darkgreen", linetype = "solid", size = 1) +
  geom_vline(xintercept = ci_high, color = "darkgreen", linetype = "solid", size = 1) +
  annotate("text", x = med, y = Inf, label = paste0("Median: ", round(med, 2)), vjust = 5.5, color = "red") +
  annotate("text", x = mn, y = Inf, label = paste0("Mean: ", round(mn, 2)), vjust = 10.5, color = "blue") +
  annotate("text", x = ci_low, y = Inf, label = paste0("5% CI: ", round(ci_low, 2)), vjust = 2.5, color = "darkgreen") +
  annotate("text", x = ci_high, y = Inf, label = paste0("95% CI: ", round(ci_high, 2)), vjust = 3.5, color = "darkgreen") +
  scale_x_continuous(limits = c(0, 300000)) +
  theme_minimal() +
  labs(
    title = "Simulated total endemic species",
    subtitle = "Using 100 random scaling factors from 95% CI",
    x = "Total simulated endemic species",
    y = "Frequency"
  )
dev.off()




quantile(scaled_results, probs = c(0.40, 0.60))

#########simulate gamma richness based on distribution of all endemism (not just since cell endemics)#########

#####Now, get an estimate of the distribution of all cell-presences for each species
#In other words, simulating endemic species only is simulating species that are present in only 1 cell.
#If I can now take presence-absence data for each species, scale it based on the cockroach data, I can simulate similar datasets and just count the number of total species

#First, take the cockroach data and plot the distribution of presences per cell.

###basic abundance histogram
roach_grid %>%
  count(Morphotype) %>%
  ggplot(aes(x = reorder(Morphotype, n), y = n)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  labs(
    title = "Morphotype Frequency in roach_grid",
    x = "Morphotype",
    y = "Number of Occurrences"
  ) +
  theme_minimal()

# Now plot the distribution of morphotypes and number of grids cells they occur in
# Step 1: Count unique grid cells per morphotype
morphotype_grid_counts <- roach_grid %>%
  distinct(Morphotype, grid_id.x) %>%
  count(Morphotype, name = "num_grid_cells")

# Step 2: Count how many morphotypes have each grid cell count
morphotype_distribution <- morphotype_grid_counts %>%
  count(num_grid_cells)

# Step 3: Plot
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
output_filename <- paste0("dist.morphotype.spat.Spread.", timestamp, ".pdf")
output_path <- file.path(output_folder, output_filename)

# Create the PDF plot
pdf(output_path, width = 8, height = 6)


ggplot(morphotype_distribution, aes(x = num_grid_cells, y = n)) +
  geom_col(fill = "lightblue") +
  labs(
    title = "Distribution of Morphotypes by Grid Cell Count",
    x = "Number of Grid Cells",
    y = "Number of Morphotypes"
  ) +
  theme_minimal()
dev.off()

###now do the same with the plant data:
###1: all the plant data

# Step 1: Count unique grid cells per plant species
plant_grid_counts <- species_counts %>%
  distinct(scientificName, grid_id) %>%
  count(scientificName, name = "num_grid_cells")

# Step 2: Count how many species have each grid cell count
plant_distribution <- plant_grid_counts %>%
  count(num_grid_cells)

# Step 3: Plot
ggplot(plant_distribution, aes(x = num_grid_cells, y = n)) +
  geom_col(fill = "lightblue") +
  labs(
    title = "Distribution of Plant Species by Grid Cell Count",
    x = "Number of Grid Cells",
    y = "Number of Plant Species"
  ) +
  theme_minimal()



###2: just the madre de dios the plant data


plant_richness <- species_counts %>%
  group_by(grid_id) %>%
  summarize(plant_richness = n_distinct(scientificName), .groups = "drop")

plant_richness_grid_narrow <- st_intersection(grid_sf, roach_endemism_filtered)

# Extract raw_endemism from list column safely
plant_richness_grid_narrow <- plant_richness_grid_narrow %>%
  mutate(raw_endemism = purrr::map_dbl(raw_endemism.x, ~ as.numeric(.x))) %>%
  filter(raw_endemism.x > 0) %>%
  dplyr::select(grid_id = grid_id.x, geometry, raw_endemism.x)

plant_richness_madre_de_dios <- st_crop(plant_richness, madre_de_dios_bbox)


                             
                             
plant_grid_joined <- st_join(grid_sf, plant_richness_madre_de_dios, by = "grid_id")
plant_grid_joined <- st_crop(plant_grid_joined, madre_de_dios_bbox)

plant_richness_grid <- plant_grid_joined %>%
  group_by(grid_id.x)%>%
  filter(plant_richness > 0) %>%
  summarise(
    plant_richness = sum(plant_richness, na.rm = TRUE),
    geometry = st_union(geometry)
  ) %>%
  ungroup()

#####plot plant richness in the madre de dios######
#Step 4: Create the plant plot
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
output_filename <- paste0("plant_richness_map_MDDios", timestamp, ".pdf")
output_path <- file.path(output_folder, output_filename)

# Open PDF device
pdf(output_path, width = 8, height = 6)
# Plant Endemism Plot (Now Properly Cropped)


ggplot(data = plant_richness_grid) +
  geom_sf(aes(fill = plant_richness), color = "black", size = 0.1) +
  scale_fill_viridis_c(option = "magma", direction = -1, name = "Plant Richness") +
  theme_minimal() +
  theme(legend.position = "right") +
  labs(
    title = "Plant Richness in Madre de Dios, Peru",
    caption = "Data: Aggregated plant richness per Amazon grid cell"
  )



# Close the device
dev.off()

######get the distribution
species_counts_joined<-st_join(grid_sf, species_counts, by = "grid_id")
species_counts_MDD <- st_crop(species_counts_joined, madre_de_dios_bbox)

filtered_species_counts <- species_counts_MDD %>%
  filter(grid_id %in% roach_grid$grid_id.x)



plant_grid_counts <- filtered_species_counts %>%
  distinct(scientificName , grid_id.x) %>%
  count(scientificName , name = "num_grid_cells")

# Step 2: Count how many morphotypes have each grid cell count
plantSpp_distribution <- plant_grid_counts %>%
  count(num_grid_cells)



#####plot spatial extent of plant richness in the madre de dios######

# Step 3: Plot
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
output_filename <- paste0("dist.PlantSpp.MDD.spat.Spread.", timestamp, ".pdf")
output_path <- file.path(output_folder, output_filename)

# Create the PDF plot
pdf(output_path, width = 8, height = 6)
 
ggplot(plantSpp_distribution, aes(x = num_grid_cells, y = n)) +
  geom_col(fill = "lightblue") +
  labs(
    title = "Distribution of Plant Species by Grid Cell Count",
    x = "Number of Grid Cells",
    y = "Number of Plant Species"
  ) +
  theme_minimal()
dev.off()

#####compate spatial extent of plant richness and cockroach richness in the madre de dios######

#ok now i want to compare the two distributions

# Add identifier columns
plantSpp_distribution$group <- "Plant species"
morphotype_distribution$group <- "Roach morphotypes"

# Bind into one dataframe
combined_distribution <- bind_rows(plantSpp_distribution, morphotype_distribution)

ggplot(combined_distribution, aes(x = factor(num_grid_cells), y = n, fill = group)) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  scale_fill_viridis_d(name = NULL, option = "magma", direction = -1) +
  labs(
    title = "Distribution of Species by Grid Cell Occupancy",
    x = "Number of Grid Cells",
    y = "Number of Species / Morphotypes"
  ) +
  theme_minimal() +
  theme(legend.position = "top")

#####using plant MDD spatial extents to extrapolate on cockroach spatial extents########
####now fit a log-normal distirbution to the plant data

# Data prep
plant_df <- plantSpp_distribution %>%
  arrange(num_grid_cells)

# Fit log-normal: convert x to numeric (meanlog, sdlog estimated manually or with nls)
x_vals <- plant_df$num_grid_cells
y_vals <- plant_df$n

# Use nls to fit a log-normal curve: f(x) = a * dlnorm(x, meanlog, sdlog)
lognorm_fit <- nls(
  y_vals ~ a * dlnorm(x_vals, meanlog, sdlog),
  start = list(a = max(y_vals), meanlog = log(2), sdlog = 0.5)
)

# Predict values from fit
x_pred <- seq(min(x_vals), max(x_vals), by = 0.1)
fit_pred <- data.frame(
  num_grid_cells = x_pred,
  fitted = predict(lognorm_fit, newdata = data.frame(x_vals = x_pred))
)

ggplot() +
  geom_point(data = plant_df, aes(x = num_grid_cells, y = n), size = 3, color = "darkgreen") +
  geom_line(data = fit_pred, aes(x = num_grid_cells, y = fitted), color = "darkgreen", size = 1.2) +
  labs(
    title = "Log-normal Fit to Plant Species Distribution",
    x = "Number of Grid Cells",
    y = "Number of Plant Species",
    caption = "Fitted using nls and dlnorm"
  ) +
  theme_minimal()


###now rescale that model distribution to match the cockroach data

# Get max value from morphotype data
roach_max <- max(morphotype_distribution$n, na.rm = TRUE)

# Scale the fitted curve to match roach amplitude
fit_pred <- fit_pred %>%
  mutate(fitted_scaled = fitted * (roach_max / max(fitted)))




####now plot it

# Add raw roach data for reference
roach_df <- morphotype_distribution %>%
  mutate(group = "Roach morphotypes")

ggplot() +
  geom_bar(data = roach_df, aes(x = factor(num_grid_cells), y = n, fill = group),
           stat = "identity", position = "identity", alpha = 0.7, color = "black") +
  geom_line(data = fit_pred, aes(x = num_grid_cells, y = fitted_scaled),
            color = "darkred", size = 1.3) +
  scale_fill_viridis_d(option = "magma", direction = -1, name = NULL) +
  labs(
    title = "Roach Morphotype Grid Distribution with Rescaled Log-normal Fit",
    x = "Number of Grid Cells",
    y = "Number of Morphotypes",
    caption = "Log-normal curve originally fit to plant data and rescaled to roach amplitude"
  ) +
  theme_minimal() +
  theme(legend.position = "top")



####ok, finally....now take the scaling factor above (roach_max / max(fitted)) and use it to scale all the plant data. I originally did this out of a function, but made it into a function because i replicate it later

grid_sf <- grid_sf %>%
  dplyr::select(-grid_id)

grid_sf <- grid_sf %>%
  rename(grid_id = grid_id.x)

simulateScaledSpeciesRichness<-function(fit_pred){
  

  
  
  plant_fit_max <- max(fit_pred$fitted)
  scaling_factor <- roach_max / plant_fit_max
  
  # 2. Apply scaling to plant_distribution
  scaled_plant_distribution <- plant_distribution %>%
    mutate(n_scaled = n * scaling_factor)
  
  
  #ggplot(scaled_plant_distribution, aes(x = num_grid_cells, y = n)) +
  #  geom_col(fill = "lightblue") +
  #  labs(
  #    title = "Distribution of Simulated Cockroach Species by Grid Cell Count for Amazon",
  #    x = "Number of Grid Cells",
  #    y = "Number of Cockroach Species"
  #  ) +
  #  theme_minimal()
  
  # Create a reference list of all available grid cells
  available_grid_ids <- grid_sf$grid_id
  available_grid_ids <- available_grid_ids[!is.na(available_grid_ids)]
  
  # Initialize list to store simulated species assignments
  sim_species_records <- list()
  
  species_counter <- 1
  
  for (i in 1:nrow(scaled_plant_distribution)) {
    
    num_sites <- scaled_plant_distribution$num_grid_cells[i]
    num_species <- round(scaled_plant_distribution$n_scaled[i])
    
    for (j in 1:num_species) {
      
      if (num_sites <= length(available_grid_ids)) {
        
        # ✔ Move name generation here
        sim_name <- str_glue("SimSpp{str_pad(species_counter, 7, pad = '0')}")
        
        assigned_cells <- sample(available_grid_ids, num_sites, replace = FALSE)
        
        sim_occurrences <- data.frame(
          scientificName = rep(sim_name, num_sites),
          grid_id = assigned_cells,
          stringsAsFactors = FALSE
        )
        
        sim_species_records[[species_counter]] <- sim_occurrences
        species_counter <- species_counter + 1
      }
    }
  }
  
  
  
  # Combine into one large data frame
  sim_species_df <- bind_rows(sim_species_records)
  
  # Join to grid_sf to get geometry
  sim_species_sf <<- sim_species_df %>%
    left_join(grid_sf, by = "grid_id") %>%
    st_as_sf()
  
  ####the line below gives how many species of cockroaches there are under this simulation
  output<-length(unique(sim_species_sf$scientificName))
  return(output)
  
}



####the line below gives how many species of cockroaches there are under this simulation
simulateScaledSpeciesRichness(fit_pred)

###there are no error ranges on this because I am using a single scaling factor. 

###what to do next?
  #a. repeat the above simulation without a scaling factor to see if you get roughly the same number of species as in the original plant data.<<DONE. it worked
  #b. think about a way to get a distribution of scaling factors so you can get CI's for this estimate <<Tried it, and it worked in theory, but in practice the result wasn't useful think about other ways to get error estimates.
  #c. plot madre-de-dios cockroach richness against plant MDD richness


###now i am going to repeat this using the plant sp. range extent distribution unmodified. if it yields roughly the same spp. richness as the raw plant data then the process is working correctly.


########Testing simulateScaledSpeciesRichness with the plant data to see if it gets the right answer#############

head(plant_distribution)

# Initialize list to store simulated species assignments
sim_species_recordsRaw <- list()

species_counter <- 1

plant_fit_max <- max(fit_pred$fitted)
scaling_factor <- roach_max / plant_fit_max

# 2. Apply scaling to plant_distribution
scaled_plant_distribution <- plant_distribution %>%
  mutate(n_scaled = n * scaling_factor)

available_grid_ids <- grid_sf$grid_id
available_grid_ids <- available_grid_ids[!is.na(available_grid_ids)]


for (i in 1:nrow(scaled_plant_distribution)) {
  
  num_sites <- scaled_plant_distribution$num_grid_cells[i]
  num_species <- round(scaled_plant_distribution$n[i])
  
  for (j in 1:num_species) {
    
    if (num_sites <= length(available_grid_ids)) {
      
      sim_name <- str_glue("SimSpp{str_pad(species_counter, 7, pad = '0')}")
      
      assigned_cells <- sample(available_grid_ids, num_sites, replace = FALSE)
      
      sim_occurrences <- data.frame(
        scientificName = rep(sim_name, num_sites),
        grid_id = assigned_cells,
        stringsAsFactors = FALSE
      )
      
      sim_species_recordsRaw[[species_counter]] <- sim_occurrences
      species_counter <- species_counter + 1
    }
  }
}


# Combine into one large data frame
sim_speciesRaw_df <- bind_rows(sim_species_recordsRaw)

#unique(sim_speciesRaw_df$scientificName)

# Join to grid_sf to get geometry
sim_speciesRaw_sf <- sim_speciesRaw_df %>%
  left_join(grid_sf, by = "grid_id") %>%
  st_as_sf()

head(sim_speciesRaw_df)
####the line below gives how many species of cockroaches there are under this simulation
# 1. Check if the richness estimators agree
is_equal <- length(unique(sim_speciesRaw_sf$scientificName)) == length(unique(species_counts$scientificName))

# 2. Compose the message based on the result
output_text <- if (is_equal) {
  "The spatial extent based richness estimator is working."
} else {
  "The spatial extent based richness estimator is NOT working."
}

# 3. Create a timestamped output filename
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
output_filename <- paste0("spatEstEst.sanityCheck.", timestamp, ".txt")

# 4. Define the full file path (change 'output_folder' if needed)
path00 <- file.path(output_folder, output_filename)

# 5. Write the result to file
writeLines(output_text, con = path00)


###Good. when i ran it on April 15 and 17 the numbers are exactly the same.

###ok now I am going to use the simulated species extents to do a randomization simulation to get 95% CIs
###i did this below one way that didn't really work. but i did it another way elsewhere. this can be ignored


############estimate gamma from total spatial extents##################
head(species_counts)

# Sum total abundance per grid cell
grid_abundance <- species_counts %>%
  group_by(grid_id) %>%
  summarise(
    total_abundance = sum(count),
    geometry = st_geometry(geometry)[1],  # Keep one representative geometry
    .groups = "drop"
  ) %>%
  st_as_sf()


set.seed(2024299)  # for reproducibility

# 1. Setup

simulated_records <- list()  # initialize an empty list
record_counter <- 1

# Loop through all grid cells
for (i in 1:nrow(grid_abundance)) {
  
  this_abundance <- grid_abundance$total_abundance[i]
  this_grid_id <- grid_abundance$grid_id[i]
  
  # Skip if abundance is zero or NA (just in case)
  if (is.na(this_abundance) || this_abundance == 0) next
  
  # Sample species names
  sampled_species <- sample(sim_species_sf$scientificName, size = this_abundance, replace = TRUE)
  
  # Create a data.frame of simulated occurrences
  df <- data.frame(
    grid_id = rep(this_grid_id, this_abundance),
    scientificName = sampled_species,
    stringsAsFactors = FALSE
  )
  
  # Store in list
  simulated_records[[record_counter]] <- df
  record_counter <- record_counter + 1
}

# Combine all rows into one data frame
simulated_occurrence_df <- do.call(rbind, simulated_records)
head(simulated_occurrence_df)

gammaRichness<-length(unique(simulated_occurrence_df$scientificName))
print(gammaRichness)

# Create the text to output
output_text <- paste0("gammaRichness: ", gammaRichness, "\n",
                      "Description: This damma diversity estimate is the projected total (cockroach) species richness of the amazon across based on a log-normal distribution of plant richness rescaled to match (cockroach) data in the madre-de-dios",
                      "In this case, the result is ", gammaRichness, ".")

timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
output_filename <- paste0("gamma_richnessfromSPRICH.", timestamp, ".txt")
path00<-file.path(output_folder, output_filename)
# Write to file in current working directory
writeLines(output_text, con =path00 )

############Redundant: SKIP##################

#one of the problems that occurred when i first ran the code below is that i was getting the same richness value every time despite random sampling. it was because the abundance was so high that i was maxing out the species pool in every simulation. The abundance was so high is because (i think) some species are extremely common within each grid cell. so we make the correction below to fix it.
#sim_species_unique <- sim_species_sf %>%
#  distinct(scientificName, grid_id, .keep_all = TRUE)

#dupes <- sim_species_sf %>%
#  count(scientificName, grid_id) %>%
#  filter(n > 1)

# nrow(dupes)  # how many (species, grid) combos are duplicated?
# 
# 
# length(sim_species_unique$scientificName)
# length(sim_species_sf$scientificName)
####the above works. now let's make it a function, and run it 100 times

simulate_gamma_richness <- function() {
  sim_species_list <- unique(sim_species_sf$scientificName)
  simulated_records <- list()
  record_counter <- 1
  
  for (i in 1:nrow(grid_abundance)) {
    this_abundance <- grid_abundance$total_abundance[i]
    this_grid_id <- grid_abundance$grid_id[i]
    
    if (is.na(this_abundance) || this_abundance == 0) next
    
    sampled_species <- sample(sim_species_unique$scientificName, size = ceiling(this_abundance), replace = TRUE)

    
    df <- data.frame(
      grid_id = rep(this_grid_id, ceiling(this_abundance)),
      scientificName = sampled_species,
      stringsAsFactors = FALSE
    )
    simulated_records[[record_counter]] <- df
    record_counter <- record_counter + 1
  }

  all_occurrences <- do.call(rbind, simulated_records)
  thisRun=length(unique(all_occurrences$scientificName))  # gamma richness
  print(thisRun)
  thisRun
}

sum(grid_abundance$total_abundance)*0.05
length(sim_species_sf$scientificName)

set.seed(4498572)  # for reproducibility

#####ALTHOUGH THIS WORKS, IT IS NOT USEFUL, SO I AM COMMENTING IT OUT SO IT ISN'T RUN EVERYTIME WE RUN THE FULL NOTEBOOK
#gamma_richness_values <- replicate(100, simulate_gamma_richness())
#table(table(sim_species_sf$scientificName))

#ggplot(data.frame(gamma = gamma_richness_values), aes(x = gamma)) +
#  geom_histogram(bins = 30, fill = "skyblue", color = "black") +
#  theme_minimal() +
#  labs(
#    title = "Simulated Gamma Richness Distribution (100 runs)",
#    x = "Gamma Richness (unique species)",
#    y = "Frequency"
#  )

###I AM GOING TO TRY A DIFFERENT WAY: by doing the bootstrapping before making the model



###########Get 90% CI's for the extent-based gamma diversity estimator#################
# Original data
x_vals <- plant_df$num_grid_cells
y_vals <- plant_df$n

# Define x range for prediction
x_pred <- seq(min(x_vals), max(x_vals), by = 0.1)

# Define a function to fit model and return predictions
fit_lognorm <- function(df) {
  tryCatch({
    fit <- nls(
      n ~ a * dlnorm(num_grid_cells, meanlog, sdlog),
      data = df,
      start = list(a = max(df$n), meanlog = log(2), sdlog = 0.5),
      control = nls.control(maxiter = 100)
    )
    preds <- predict(fit, newdata = data.frame(num_grid_cells = x_pred))
    return(preds)
  }, error = function(e) {
    return(rep(NA, length(x_pred)))
  })
}

# Bootstrapping: resample and fit multiple times
set.seed(123)
n_boot <- 200


# Expand plant_df into individual species occurrences
species_expanded <- plant_df %>%
  uncount(weights = n)  # creates one row per species, based on n

# Check the distribution
table(species_expanded$num_grid_cells)

set.seed(123)
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
lognorm_lower_fit <- nls(
  y_valsLower ~ a * dlnorm(x_valsLower, meanlog, sdlog),
  start = list(a = max(y_valsLower), meanlog = log(2), sdlog = 0.5),
  control = nls.control(maxiter = 200, warnOnly = TRUE)
)

# Fit log-normal model to lower CI
lognorm_Upper_fit <- nls(
  y_valsUpper ~ a * dlnorm(x_valsUpper, meanlog, sdlog),
  start = list(a = max(y_valsUpper), meanlog = log(2), sdlog = 0.5),
  control = nls.control(maxiter = 200, warnOnly = TRUE)
)


# Predict over range
x_pred <- ci_df$num_grid_cells
y_pred <- predict(lognorm_lower_fit, newdata = data.frame(x_vals = x_pred))

y_predUpper <- predict(lognorm_Upper_fit, newdata = data.frame(x_vals = x_pred))



timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
output_filename <- paste0("plantPatExtentModel.MDD.", timestamp, ".pdf")
output_path <- file.path(output_folder, output_filename)

# Open PDF device
pdf(output_path, width = 8, height = 6)

# Plot
ggplot(ci_df, aes(x = num_grid_cells)) +
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
  labs(
    title = "Log-normal Fit to Plant Species Distribution with Bootstrapped CIs",
    x = "Number of Grid Cells",
    y = "Number of Plant Species"
  ) +
  theme_minimal()
dev.off()

fit_predUpper <- data.frame(
  num_grid_cells = x_pred,
  fitted = predict(lognorm_Upper_fit, newdata = data.frame(x_vals = x_pred))
)

fit_predLower <- data.frame(
  num_grid_cells = x_pred,
  fitted = predict(lognorm_lower_fit, newdata = data.frame(x_vals = x_pred))
)

#######CALCULATE (cockroach) gamma richness based on species richness scaling########
gammaRichnessSPRIch<-simulateScaledSpeciesRichness(fit_pred)##best estimate
upperEstSpRich<-simulateScaledSpeciesRichness(fit_predUpper)###upper estimate
lowerEstSpRich<-simulateScaledSpeciesRichness(fit_predLower)###lower estimate

output_text <- paste0("gammaRichness: ", gammaRichnessSPRIch, "\n",
                      "Description: This damma diversity estimate is the projected total (cockroach) species richness of the amazon across based on a log-normal distribution of plant richness rescaled to match (cockroach) data in the madre-de-dios",
                      "In this case, the result is ", gammaRichness, ".\n\n\nThe estimates range from ",upperEstSpRich, " to ", lowerEstSpRich, ". This is based on bootstrap resampling at 95% and 105% of the original dataset size" )

timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
output_filename <- paste0("gamma_richnessfromSPRICH.", timestamp, ".txt")
path00<-file.path(output_folder, output_filename)
# Write to file in current working directory
writeLines(output_text, con =path00 )


#### plotting madre-de-dios cockroach richness against plant MDD richness###############

ggplot(data = plant_richness_grid) +
  geom_sf(aes(fill = plant_richness), color = "black", size = 0.1) +
  scale_fill_viridis_c(option = "magma", direction = -1, name = "Plant Richness") +
  theme_minimal() +
  theme(legend.position = "right") +
  labs(
    title = "Plant Richness in Madre de Dios, Peru",
    caption = "Data: Aggregated plant richness per Amazon grid cell"
  )



roach_richness <- roach_grid %>%
  group_by(grid_id.x) %>%
  summarise(
    morphotype_richness = n_distinct(Morphotype),
    geometry = geometry[1]  # retain one geometry per grid_id
  ) %>%
  sf::st_as_sf()


# Ensure both layers are in the same CRS
plant_richness_grid <- st_transform(plant_richness_grid, crs = st_crs(roach_grid))

# Filter to plant grid cells that overlap with any cockroach points
plant_filtered <- st_filter(plant_richness_grid, roach_grid)

ggplot() +
  geom_sf(data = plant_richness_grid, fill = "grey90", color = "black") +
  geom_sf(data = plant_filtered, fill = "green", alpha = 0.5) +
  geom_sf(data = roach_grid, color = "red", size = 0.5) +
  labs(title = "Filtered Plant Richness Grid Overlapping Roach Points")

# First ensure plant_filtered and roach_richness have matching grid IDs
names(plant_filtered)
names(roach_richness)

# Adjust column names if needed
roach_richness <- roach_richness %>%
  rename(grid_id = grid_id.x)

# Perform spatial join by ID
richness_joined <- st_join(plant_filtered, roach_richness, by = "grid_id")



model <- lm(morphotype_richness ~ plant_richness, data = richness_joined)
summary(model)





# Plot raw endemism correlation
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
output_filename <- paste0("plant_vs_insect.RichnessCorr.", timestamp, ".pdf")
output_path <- file.path(output_folder, output_filename)

# Open PDF device
pdf(output_path, width = 8, height = 6)



# Extract R² value
r2 <- summary(model)$r.squared
r2_text <- sprintf("R² = %.3f", r2)  # Format to 3 decimal places



ggplot(richness_joined, aes(x = plant_richness, y = morphotype_richness)) +
  geom_point(color = "darkgreen", size = 3, alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  theme_minimal() +
  labs(
    title = "Plant vs. Cockroach Richness per Grid Cell",
    x = "Plant Richness",
    y = "Roach Morphotype Richness",
    caption = paste("Each point represents a grid cell —", r2_text)
  )

# Close the device
dev.off()





###(NOT)PROBLEM: THE RAW ENDEMISM PLOT EARLY IN THIS NOTEBOOK MAKES IT SEEM LIKE THERE ARE ONLY A FEW DOZEN ENDEMIC PLANTS. BUT THE HISTOGRAM (FROM SPECIES COUNTS) MAKES IT SEEM LIKE THERE ARE >1000. WHAT IS THE SOURCE OF THIS DISPARTIY? i THINK IT MIGHT BE THAT THE EARLY PLOT WAS CREATED AFTER REMOVING ALL NON-ENDEMICS BUT THE LATTER ONE WAS DONE USING THE FULL SPECIES LIST. however, I NEED TO DOUBLE-CHECK IF THIS IS THE CASE. IF IT IS INDEED THE CASE THEN I THINK THE DATA IS OK, IT'S JUST THAT EACH ONE RELIES ON DIFFERENT ASSUMPTIONS FOR THE COCKROACHES. IN THE FIRST ONE (WITH FEW ENDEMIC PLANTS) THE ASSUMPTION OF THAT MODEL IS THAT THE COCKROACHES ARE ALL TRULY ENDEMIC. IN THE LATTER ONE, THE ASSUMPTION OF THAT MODEL IS THAT WEIGHTED ENDEMISM IN THE COCKROACHES HAS THE SAME DISTRIBUTION AS THE WEIGHTED ENDEMISM IN PLANTS. SO THAT SPATIALLY NARROW ENDEMISM (I.E., ENDEMISM WHERE YOU DON'T CONSIDER SPP THAT OCCUR OUTSIDE OF THE AREA CONSIDERED) WOULD HAVE THE SAME DISTRIBUTION IN PLANTS AND COCKROACHES.
#ANSWER: the raw_endemism object was indeed created after filtering out all species that occur in more than 1 grid cell. so, it's ok.






