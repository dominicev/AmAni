# 01_GBIFAnalysis.R — Entry point for standard GBIF occurrence CSV downloads
# Expects .amani_data_file to be set by the manager before sourcing.

library(data.table)
library(dplyr)
library(tidyr)
library(sf)

if (exists(".amani_wd")) {
  setwd(.amani_wd)
} else {
  project_root <- dirname(dirname(rstudioapi::getActiveDocumentContext()$path))
  setwd(project_root)
}

###Output results folder####
if (exists(".amani_output_dir")) {
  output_folder <- .amani_output_dir
} else {
  output_folder <- "output"
}

################Plant import and setup###################

dataDir <- paste0(dirname(getwd()), "/1-FilteringPlantData/output")

#occ_sf <- readRDS(paste0(dataDir, "/species_counts4_filtered.rds"))

# Extract the grid from the plant data to preserve Stage 0 grid_ids.
# Do NOT rebuild via st_make_grid() — different bounding box produces
# different grid_ids that break downstream joins with species_counts.
grid_size_m <- 22360
grid_sf <- occ_sf |>
  dplyr::group_by(grid_id) |>
  dplyr::slice(1) |>
  dplyr::ungroup() |>
  dplyr::select(grid_id, geometry)

################GBIF DATA IMPORT AND SETUP###################

if (!exists(".amani_data_file")) {
  stop("Set .amani_data_file to the path of the GBIF CSV before running.", call. = FALSE)
}

message("Reading GBIF data from: ", .amani_data_file)
ext <- tolower(tools::file_ext(.amani_data_file))
# Read CSV directly from a zip archive if needed — streams the file
# through unzip to stdout so nothing is extracted to disk.
# 
# 
#FUTURE NOTE: Make fread_auto a stand-alone function in a helper file. different version of it are called in this script and the cleaning script
fread_auto <- function(path, ...) {
  file_sep <- if (ext == "tsv") "\t" else ","
  if (tolower(tools::file_ext(path)) == "zip") {
    # Find the CSV inside the archive
    #csv_name <- grep("\\.csv$", unzip(path, list = TRUE)$Name, value = TRUE)
    csv_name <- grep("\\.(csv|tsv)$", unzip(path, list = TRUE)$Name, ignore.case = TRUE, value = TRUE) #fix to add TSV
    
    # Check inner extension for zipped files
    inner_ext <- tolower(tools::file_ext(csv_name))
    inner_sep <- if (inner_ext == "tsv") "\t" else ","
    
    if (length(csv_name) == 0) stop("No CSV found inside: ", path, call. = FALSE)
    if (length(csv_name) > 1) {
      message("  Multiple CSVs in archive; using: ", csv_name[1])
      csv_name <- csv_name[1]
    }
    cmd <- sprintf('unzip -p "%s" "%s"', normalizePath(path, winslash = "/"), csv_name)
    data.table::fread(cmd = cmd, ...)
  } else {
    data.table::fread(file = path, sep = file_sep, ...)  
  }
}
#Peek at the first row of the file to see what headers actually exist
file_headers <- names(data.table::fread(.amani_data_file, nrows = 0))
print(paste0("There are ", toString( length(file_headers)), " header columns in the raw data file."))


# Define the columns your model ideally wants
desired_columns <- c("gbifID", "species", "taxonRank", "scientificName",
                     "decimalLatitude", "decimalLongitude",
                     "coordinateUncertaintyInMeters", "occurrenceStatus", "issue",
                     "establishmentMeans", "year")

#Only select the columns that are actually present in the file
columns_to_keep <- intersect(desired_columns, file_headers)

#Read the file safely
gbif_raw <- fread_auto(
  .amani_data_file,
  select = columns_to_keep,
  quote = "\"" ####NOTE: This line might cause errors....change it to make it flexible based on the status of commas in the actual file headers
)

# Fall back to scientificName when no species column is available
# (mirrors the logic in 00_cleanGBIF.R)
if (!"species" %in% names(gbif_raw) && "scientificName" %in% names(gbif_raw)) {
  message("  'species' column not found; using scientificName as species")
  data.table::setnames(gbif_raw, "scientificName", "species")
} else if ("species" %in% names(gbif_raw) && "scientificName" %in% names(gbif_raw)) {
  # Fill blank/NA species values from scientificName
  gbif_raw[is.na(species) | species == "", species := scientificName]
}

# Coordinates may arrive as character columns (e.g. from curated CSVs), and
# some rows contain non-decimal formats (DMS strings, blanks, literal "NA")
# that silently become NA on coercion. Catch that here — before filtering —
# so bad rows are dropped with a clear message instead of failing later
# inside st_as_sf().
gbif_raw[, decimalLatitude := suppressWarnings(as.numeric(decimalLatitude))]
gbif_raw[, decimalLongitude := suppressWarnings(as.numeric(decimalLongitude))]

n_bad_coords <- sum(is.na(gbif_raw$decimalLatitude) | is.na(gbif_raw$decimalLongitude))
if (n_bad_coords > 0) {
  message(
    "  Dropping ", n_bad_coords, " record(s) with non-numeric or missing coordinates ",
    "(e.g. DMS-format strings, blanks, literal \"NA\") that could not be coerced to decimal degrees."
  )
}

# Basic filters
gbif_clean <- gbif_raw[
  occurrenceStatus == "PRESENT" &
    !is.na(decimalLatitude) &
    !is.na(decimalLongitude) &
    species != ""
]

message("  Records after filtering: ", nrow(gbif_clean),
        " (", uniqueN(gbif_clean$species), " species)")

# Convert to sf and reproject to match grid
animal_sf <- st_as_sf(
  gbif_clean,
  coords = c("decimalLongitude", "decimalLatitude"),
  crs = 4326
) |>
  st_transform(st_crs(grid_sf))

# Spatial join: keep only occurrences falling within the Amazon grid
animal_grid <- st_join(animal_sf, grid_sf, join = st_within) |>
  filter(!is.na(grid_id))

# Map "species" to "Morphotype" so downstream scripts work unchanged
animal_grid$Morphotype <- animal_grid$species

message("  Records within Amazon grid: ", nrow(animal_grid),
        " (", n_distinct(animal_grid$Morphotype), " species)")

# Richness per grid cell
animal_endemism <- animal_grid |>
  group_by(grid_id) |>
  summarize(animal_richness = n_distinct(Morphotype), .groups = "drop")

animal_grid_data <- grid_sf |>
  dplyr::select(grid_id, geometry) |>
  dplyr::left_join(
    animal_endemism |> sf::st_drop_geometry() |> dplyr::select(grid_id, animal_richness),
    by = "grid_id"
  ) |>
  mutate(animal_richness = replace_na(animal_richness, 0))

################RUN DOWNSTREAM PIPELINE###################

source("scripts/utils/ModelAnimal.R")

modelAnimal(animal_grid_data,
            redoeverything = ifelse(exists(".amani_redo"), .amani_redo, TRUE), 
            fill_animals = if (exists(".amani_fillAnimals")) .amani_fillAnimals else FALSE,
            alphaRichnessMethod = if (exists(".amani_alphaRichnessMethod")) .amani_alphaRichnessMethod else "empirical",
            fit_params = .amani_fit_params, estimateBounds = .amani_estimateBounds)
