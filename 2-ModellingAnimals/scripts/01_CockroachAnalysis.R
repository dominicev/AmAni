# set working directory to project root (2-ModellingAnimals/)
library(rstudioapi)
if (exists(".amani_wd")) {
  setwd(.amani_wd)
} else {
  project_root <- dirname(dirname(rstudioapi::getActiveDocumentContext()$path))
  setwd(project_root)
}

library(readxl)
library(dplyr)
library(tidyr)
library(parzer)
library(sf)
#library(readr)

###Output results folder####
if (exists(".amani_output_dir")) {
  output_folder <- .amani_output_dir
} else {
  output_folder <- "output/analysis_results"
}

################Plant import and setup###################
###step 1 make the grid

dataDir<-paste0(dirname(getwd()),"/1-FilteringPlantData/output")

#occ_sf<-readRDS(paste0(dataDir,"/species_counts4_filtered.rds"))

# Extract the grid from the plant data to preserve Stage 0 grid_ids.
# Do NOT rebuild via st_make_grid() — different bounding box produces
# different grid_ids that break downstream joins with species_counts.
grid_size_m <- 22360
grid_sf <- occ_sf |>
  dplyr::group_by(grid_id) |>
  dplyr::slice(1) |>
  dplyr::ungroup() |>
  dplyr::select(grid_id, geometry)

######ANIMAL DATA IMPORT AND SETUP#######
if (!exists(".amani_data_file")) {
  stop("Set .amani_data_file to the path of the Darwin Core Excel file before running.", call. = FALSE)
}
file_path <- .amani_data_file

# Load the data, skipping the first row to correctly set headers
animal_raw_data <- read_excel(file_path, sheet = 1, skip = 1)
selected_animal_data <- animal_raw_data %>%
  dplyr::select(
    Morphotype = 47,
    Sex = 50,
    Life_Stage = 51,
    Event_Date = 56,
    verB_Lat = 54,
    verb_Long = 55,
    Locality = 53
  )

selected_animal_data <- selected_animal_data %>%
  mutate(
    decimal_lat = parse_lat(verB_Lat),
    decimal_lon = parse_lon(verb_Long)
  )

# Remove rows with missing coordinates
selected_animal_data <- selected_animal_data %>%
  filter(!is.na(decimal_lat) & !is.na(decimal_lon))

###filter out juvenile-only morphotypes
adult_morphotypes <- selected_animal_data %>%
  filter(Life_Stage == "Adult") %>%
  distinct(Morphotype) %>%
  pull(Morphotype)  # Extract as a vector

selected_animal_data2 <- selected_animal_data %>%
  filter(Morphotype %in% adult_morphotypes)


# Convert animal occurrence data into an sf object
animal_sf <- selected_animal_data2 %>%
  st_as_sf(coords = c("decimal_lon", "decimal_lat"), crs = 4326) %>%  # WGS84
  st_transform(st_crs(grid_sf))  # Match CRS with grid

# Perform a spatial join: Assign each animal occurrence to a grid cell
animal_grid <- st_join(animal_sf, grid_sf, join = st_within)

animal_endemism <- animal_grid %>%
  group_by(grid_id) %>%
  summarize(animal_richness = n_distinct(Morphotype), .groups = "drop")

animal_grid_data <- grid_sf |>
  dplyr::select(grid_id, geometry) |>
  dplyr::left_join(
    animal_endemism |> sf::st_drop_geometry() |> dplyr::select(grid_id, animal_richness),
    by = "grid_id"
  ) |>
  mutate(animal_richness = replace_na(animal_richness, 0))

source("scripts/utils/ModelAnimal.R")

modelAnimal(animal_grid_data, 
            redoeverything = ifelse(exists(".amani_redo"), .amani_redo, TRUE), 
            fill_animals = if (exists(".amani_fillAnimals")) .amani_fillAnimals else FALSE,
            alphaRichnessMethod = if (exists(".amani_alphaRichnessMethod")) .amani_alphaRichnessMethod else "empirical",
            fit_params = .amani_fit_params)

#source("scripts/utils/rawEndemismCalculator.R")
#source("scripts/animalAnalysis.R")

#the spatial clustering works
