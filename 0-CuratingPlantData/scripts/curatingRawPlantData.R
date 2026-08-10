# set working directory to project root (0-CuratingPlantData/)
library(rstudioapi)
if (exists(".amani_wd")) {
  setwd(.amani_wd)
} else {
  project_root <- dirname(dirname(rstudioapi::getActiveDocumentContext()$path))
  setwd(project_root)
}

#load dependencies
library(data.table)
library(sf)
library(lwgeom)
library(stringr)
library(dplyr)
library(purrr)
library(WorldFlora)
library(digest)
library(tidyterra)
library(ggplot2)
library(elevatr)
library(terra)

zip_file <- file.path("data", "0017480-250310093411724.zip")

data <- fread(
  zip_file,
  select = c(
    "scientificName",
    "decimalLatitude",
    "decimalLongitude",
    "coordinateUncertaintyInMeters",
    "coordinatePrecision",
    "elevation",
    "elevationAccuracy",
    "countryCode",
    "stateProvince",
    "eventDate",
    "year",
    "occurrenceStatus"
  )
)

message("--- Stage 0: Curating Plant Data ---")
message("Raw GBIF records loaded: ", format(nrow(data), big.mark = ","))

# --- Data quality filters (consistent with 00_cleanGBIF.R for animals) ---

# Remove absence records
if ("occurrenceStatus" %in% names(data)) {
  n_before <- nrow(data)
  data <- data[is.na(occurrenceStatus) | toupper(occurrenceStatus) != "ABSENT"]
  message("Presence records only: ", format(nrow(data), big.mark = ","),
          " (", n_before - nrow(data), " absence records removed)")
}

# Remove implausible collection years
if ("year" %in% names(data)) {
  n_before <- nrow(data)
  current_year <- as.integer(format(Sys.Date(), "%Y"))
  data <- data[is.na(year) | (year >= 1800 & year <= current_year)]
  message("Year filter (>= 1800, <= ", current_year, "): ",
          format(nrow(data), big.mark = ","),
          " (", n_before - nrow(data), " removed)")
}

# Remove records with coordinate uncertainty > 10 km (NAs kept)
if ("coordinateUncertaintyInMeters" %in% names(data)) {
  n_before <- nrow(data)
  data <- data[is.na(coordinateUncertaintyInMeters) | coordinateUncertaintyInMeters <= 10000]
  message("Coordinate uncertainty filter (<= 10 km, NAs kept): ",
          format(nrow(data), big.mark = ","),
          " (", n_before - nrow(data), " removed)")
}

# Remove records with missing coordinates
n_before <- nrow(data)
data <- data[!is.na(decimalLatitude) & !is.na(decimalLongitude)]
message("NA coordinates removed: ", n_before - nrow(data),
        " | Remaining: ", format(nrow(data), big.mark = ","))

# Remove exact spatial duplicates (same scientificName + lat + lon)
n_before <- nrow(data)
data <- unique(data, by = c("scientificName", "decimalLatitude", "decimalLongitude"))
message("Exact spatial duplicates removed: ",
        format(n_before - nrow(data), big.mark = ","),
        " | Remaining: ", format(nrow(data), big.mark = ","))
#head(data)
saveRDS(data, "output/filtered_gbif_data2.rds")

#########filter occurences to only those in the Amazon basin

# Download amazon shapefile if not already present
dest_dir <- "data/amazon_shapefile"
shp_path <- file.path(dest_dir, "amazon_sensulatissimo_gmm_v1.shp")

if (!file.exists(shp_path)) {
  base_url <- "https://github.com/gamamo/AmazonBasinLimits/raw/main/"
  files <- c("amazon_sensulatissimo_gmm_v1.shp",
             "amazon_sensulatissimo_gmm_v1.shx",
             "amazon_sensulatissimo_gmm_v1.dbf",
             "amazon_sensulatissimo_gmm_v1.prj")
  if (!dir.exists(dest_dir)) dir.create(dest_dir)
  for (file in files) {
    download.file(paste0(base_url, file), destfile = file.path(dest_dir, file), mode = "wb")
  }
  message("Amazon shapefile downloaded.")
} else {
  message("Amazon shapefile already present, skipping download.")
}

amazon_shape <- st_read(shp_path, quiet = TRUE)
amazon_shape <- st_make_valid(amazon_shape)
amazon_shape <- st_simplify(amazon_shape, dTolerance = 0.01)



# 2. Download elevation data matching the shapefile's extent
# z = 5 or 6 provides a good balance between resolution and download speed for large areas like the Amazon
elevation_raster <- get_elev_raster(amazon_shape, z = 5, src = "aws")

# 3. Convert to a terra SpatRaster object for easier manipulation
el_terra <- rast(elevation_raster)

# 4. Crop and mask the elevation data precisely to the shapefile boundary
el_masked <- crop(el_terra, amazon_shape)
el_masked <- mask(el_masked, vect(amazon_shape))

# 5. Get the min and max elevation values (ignoring NA values outside the boundary)
elev_values <- values(el_masked, na.rm = TRUE)
min_elev <- min(elev_values)
max_elev <- max(elev_values)

# Print results
cat("Minimum Elevation:", round(min_elev, 2), "meters\n")
cat("Maximum Elevation:", round(max_elev, 2), "meters\n")

ggplot() +
  # Plot the topographic raster layer
  geom_spatraster(data = el_masked) +
  # Use a terrain color palette
  scale_fill_terrain_c(name = "Elevation (m)", na.value = "transparent", guide_colorbar(theme = 
                      theme(axis.text.x = element_text(angle = 90, hjust = 1)))) +
  # Overlay the shapefile border
  geom_sf(data = amazon_shape, fill = NA, color = "black", linewidth = 0.5) +
  # Add titles and labels
  labs(
    title = "Topographical Map of the Amazon Basin",
    subtitle = paste0("Elevational Range: ", min_elev, "m to ", max_elev, "m"),
    x = "Longitude",
    y = "Latitude"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

rm(el_masked, el_terra, elevation_raster)

# Convert GBIF data to sf object
data_sf <- st_as_sf(data, coords = c("decimalLongitude", "decimalLatitude"), crs = 4326)

# [Speedup C] Crop to Amazon bounding box first (fast bbox test),
# then filter against the actual polygon (on the smaller set)
data_cropped <- st_crop(data_sf, st_bbox(amazon_shape))
data_filtered <- st_filter(data_cropped, amazon_shape)
rm(data_cropped)
message("Records after Amazon spatial filter: ", format(nrow(data_filtered), big.mark = ","))

saveRDS(data_filtered, "output/amazon_filtered_data2.rds")

#########Generate endemism model##########

# [Speedup B] Build grid in EPSG:3857, then transform just the grid
# polygons (~20k) back to 4326 for the spatial join — avoids reprojecting
# millions of individual points.
# EPSG:3857 (Web Mercator) is used consistently across all pipeline stages
# for the species grid. See stages_0_1_issues.txt #4 for rationale.

grid_size_m <- 22360  # 22.36 km ≈ 500 km² cells

amazon_3857 <- st_transform(amazon_shape, 3857)
grid <- st_make_grid(amazon_3857, cellsize = grid_size_m, what = "polygons")
grid_sf <- st_sf(grid_id = 1:length(grid), geometry = grid, crs = 3857)

# Transform grid to 4326 for joining with occurrence data (which stays in 4326)
grid_sf_4326 <- st_transform(grid_sf, 4326)

# Spatial join: assign grid_id to each occurrence
message("Spatial join: assigning grid cells to occurrences...")
occ_grid <- st_join(data_filtered, grid_sf_4326, left = FALSE)

# [Speedup A] Drop geometry and switch to data.table for all aggregation
# and name cleaning — geometry is only needed for final output.
occ_dt <- as.data.table(st_drop_geometry(occ_grid))
rm(occ_grid, data_filtered, data_sf, grid_sf_4326)
gc()

#(a) Compute Species Occurrences per Grid
species_dt <- occ_dt[, .(count = .N), by = .(grid_id, scientificName)]
rm(occ_dt)

message("Grid cells with data: ", format(uniqueN(species_dt$grid_id), big.mark = ","),
        " | Species x cell rows: ", format(nrow(species_dt), big.mark = ","),
        " | Unique names: ", format(uniqueN(species_dt$scientificName), big.mark = ","))

# Name cleaning (data.table + stringr)
species_dt[, c("genus", "epithet") := {
  parts <- str_split(scientificName, " ")
  list(
    vapply(parts, `[`, character(1), 1),
    vapply(parts, function(x) if (length(x) > 1) x[2] else NA_character_, character(1))
  )
}]

# Filter invalid epithets
species_dt[str_detect(epithet, "\\."), epithet := NA_character_]
species_dt[str_detect(epithet, "^[A-Z]"), epithet := NA_character_]
species_dt[str_detect(epithet, "^[^a-zA-Z]"), epithet := NA_character_]
species_dt[, cleaned_name := fifelse(!is.na(epithet), paste(genus, epithet), NA_character_)]

species_dt <- species_dt[!is.na(cleaned_name)]
species_dt[, scientificName := cleaned_name]
species_dt <- species_dt[, .(grid_id, scientificName, count)]

message("After name cleaning: ",
        format(nrow(species_dt), big.mark = ","), " rows | ",
        format(uniqueN(species_dt$scientificName), big.mark = ","), " species | ",
        format(uniqueN(species_dt$grid_id), big.mark = ","), " grid cells")

# [Speedup D] Cache WFO matching — only re-run if the species list changed
species_vec <- unique(species_dt$scientificName)
wfo_cache_path <- "output/wfo_data.rds"
species_hash_path <- "output/wfo_species_hash.rds"

run_wfo <- TRUE
if (file.exists(wfo_cache_path) && file.exists(species_hash_path)) {
  current_hash <- digest::digest(sort(species_vec), algo = "xxhash64")
  cached_hash  <- readRDS(species_hash_path)
  if (identical(current_hash, cached_hash)) {
    message("WFO cache hit — species list unchanged, loading cached matches.")
    wfo_best <- readRDS(wfo_cache_path)
    run_wfo <- FALSE
  } else {
    message("WFO cache miss — species list changed, re-matching.")
  }
}

if (run_wfo) {
  wfo_matches <- WFO.match(species_vec,
                           WFO.file = "data/_DwC_backbone_R/classification.csv")

  # Keep one best match per input name, preferring Accepted over Synonym
  wfo_best <- wfo_matches %>%
    filter(!is.na(scientificName)) %>%
    mutate(
      status_rank = case_when(
        taxonomicStatus == "Accepted" ~ 1L,
        taxonomicStatus == "Synonym"  ~ 2L,
        TRUE                          ~ 3L
      )
    ) %>%
    group_by(spec.name.ORIG) %>%
    arrange(status_rank) %>%
    slice(1) %>%
    ungroup() %>%
    select(-status_rank)

  saveRDS(wfo_best, wfo_cache_path)
  saveRDS(digest::digest(sort(species_vec), algo = "xxhash64"), species_hash_path)
}

n_matched   <- sum(!is.na(wfo_best$scientificName))
n_accepted  <- sum(wfo_best$taxonomicStatus == "Accepted", na.rm = TRUE)
n_synonym   <- sum(wfo_best$taxonomicStatus == "Synonym", na.rm = TRUE)
message("WFO matching: ", format(n_matched, big.mark = ","), " names matched (",
        n_accepted, " accepted, ", n_synonym, " synonyms, ",
        n_matched - n_accepted - n_synonym, " other)")

# Apply corrected names
wfo_lookup <- setDT(wfo_best[, c("spec.name.ORIG", "scientificName")])
setnames(wfo_lookup, c("spec.name.ORIG", "scientificName"), c("old_name", "corrected_name"))
species_dt <- merge(species_dt, wfo_lookup, by.x = "scientificName", by.y = "old_name", all.x = TRUE)
species_dt[!is.na(corrected_name), scientificName := corrected_name]
species_dt[, corrected_name := NULL]

# Re-aggregate counts — name corrections may have merged distinct names
species_dt <- species_dt[, .(count = sum(count)), by = .(grid_id, scientificName)]

message("After WFO correction + re-aggregation: ",
        format(nrow(species_dt), big.mark = ","), " rows | ",
        format(uniqueN(species_dt$scientificName), big.mark = ","), " species | ",
        format(uniqueN(species_dt$grid_id), big.mark = ","), " grid cells")

# Rejoin grid polygon geometry (in 3857) for the final sf output
grid_ids_with_data <- unique(species_dt$grid_id)
grid_geom <- grid_sf[grid_sf$grid_id %in% grid_ids_with_data, ]

species_counts_final <- merge(grid_geom, species_dt, by = "grid_id", all.x = FALSE)
species_counts_final <- st_as_sf(species_counts_final)

saveRDS(species_counts_final, "output/species_counts3.rds")



# --- QA visualization ---#(b) Calculate Weighted Endemism
# data.table for speed; geometry not needed
range_sizes <- species_dt[, .(range_size = uniqueN(grid_id)), by = scientificName]
species_dt_e <- merge(species_dt, range_sizes, by = "scientificName")
species_dt_e[, endemism_weight := 1 / range_size]
weighted_endemism <- species_dt_e[, .(weighted_endemism = sum(endemism_weight)), by = grid_id]

saveRDS(as.data.frame(weighted_endemism), "output/species_endemism2.rds")
message("Weighted endemism: ", nrow(weighted_endemism), " grid cells, ",
        "range ", round(min(weighted_endemism$weighted_endemism), 2),
        " - ", round(max(weighted_endemism$weighted_endemism), 2))

#(c) Calculate Raw Endemism (count of single-cell species per grid cell)
endemics_dt <- range_sizes[range_size == 1]
endemic_records <- species_dt[scientificName %in% endemics_dt$scientificName]
raw_endemism <- endemic_records[, .(raw_endemism = uniqueN(scientificName)), by = grid_id]

saveRDS(as.data.frame(raw_endemism), "output/raw_species_endemism2.rds")
message("Raw endemism: ", nrow(raw_endemism), " grid cells with endemics, ",
        nrow(endemics_dt), " single-cell species")
message("--- Stage 0 complete ---")
# Outlier detection and plotting for weighted + raw endemism.
# Uses the same threshold as Stage 1 (Q95 + 3*IQR).
# Kept here for reference if quick QA of Stage 0 outputs is needed.
#
# grid_sf <- st_join(grid_sf, species_endemism, by = "grid_id")
# q95  <- quantile(grid_sf$weighted_endemism, 0.95, na.rm = TRUE)
# iqr  <- IQR(grid_sf$weighted_endemism, na.rm = TRUE)
# upper_cutoff <- q95 + 3 * iqr
# grid_sf <- grid_sf |>
#   dplyr::mutate(weight_endemism_filt = dplyr::if_else(
#     weighted_endemism > upper_cutoff, NA_real_, as.numeric(weighted_endemism)))
# ggplot() +
#   geom_sf(data = grid_sf, aes(fill = weight_endemism_filt), color = "black") +
#   scale_fill_viridis_c(option = "plasma", na.value = "white") +
#   theme_minimal() +
#   labs(title = "Weighted Endemism (QA check)", fill = "Endemism")
