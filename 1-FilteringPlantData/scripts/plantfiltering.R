# set working directory to project root (1-FilteringPlantData/)
library(rstudioapi)
if (exists(".amani_wd")) {
  setwd(.amani_wd)
} else {
  project_root <- dirname(dirname(rstudioapi::getActiveDocumentContext()$path))
  setwd(project_root)
}

library(sf)
library(ggplot2)
library(dplyr)
library(svglite)
library(rnaturalearth)

###set directory for loading data from upstream stage
dataDir<-paste0(dirname(getwd()),"/0-CuratingPlantData/output")

message("--- Stage 1: Filtering Plant Data ---")

###load species occurences
species_counts<-readRDS(paste0(dataDir,"/species_counts3.rds"))

###load weighted endemism
species_endemism<-readRDS(paste0(dataDir,"/species_endemism2.rds"))

message("Loaded from Stage 0: ",
        format(nrow(species_counts), big.mark = ","), " rows | ",
        format(n_distinct(species_counts$scientificName), big.mark = ","), " species | ",
        format(n_distinct(species_counts$grid_id), big.mark = ","), " grid cells")

grid_sf<-species_counts

# --- Filter non-native species using GBIF establishmentMeans ---
# Reads establishmentMeans from the raw GBIF download. Species tagged
# "introduced" (n > 4 records) with zero "native" records are removed
# entirely — all occurrences of that species are dropped, even if individual
# records lacked the tag. A whitelist protects native species sometimes
# recorded in cultivation contexts.
introduced_cache <- "output/introduced_species.rds"
gbif_zip <- paste0(dirname(dataDir), "/data/0017480-250310093411724.zip")

if (file.exists(introduced_cache)) {
  message("Loading cached introduced species list...")
  introduced_final <- readRDS(introduced_cache)
} else {
  message("Reading establishmentMeans from GBIF zip (slow, ~30s)...")
  em_dt <- data.table::fread(
    cmd = paste0("unzip -p '", gbif_zip, "' '*.csv'"),
    select = c("scientificName", "establishmentMeans"),
    sep = "\t"
  )
  em_dt <- em_dt[establishmentMeans != ""]

  # Extract binomial (first two words) to match WFO-corrected names
  em_dt[, binomial := sub("^(\\S+\\s+\\S+).*", "\\1", scientificName)]

  em_summary <- em_dt[,
    .(n_native = sum(tolower(establishmentMeans) == "native"),
      n_introduced = sum(tolower(establishmentMeans) %in%
                           c("introduced", "managed", "invasive"))),
    by = binomial]

  native_whitelist <- c("Psidium guajava", "Arachis hypogaea",
                        "Capsicum frutescens")
  introduced_strict <- em_summary[n_introduced > 4 & n_native == 0, binomial]
  introduced_strict <- setdiff(introduced_strict, native_whitelist)

  # Map through WFO name corrections
  wfo_cache <- paste0(dataDir, "/wfo_data.rds")
  if (exists("wfo_lookup")) {
    wfo_dt <- data.table::as.data.table(wfo_lookup)
  } else if (file.exists(wfo_cache)) {
    wfo_dt <- data.table::as.data.table(readRDS(wfo_cache)$lookup)
  } else {
    stop("WFO lookup not found — run Stage 0 first or ensure wfo_data.rds exists")
  }
  corrected_from_wfo <- wfo_dt[old_name %in% introduced_strict, corrected_name]
  direct_matches <- intersect(introduced_strict, unique(grid_sf$scientificName))
  introduced_final <- unique(c(corrected_from_wfo, direct_matches))

  saveRDS(introduced_final, introduced_cache)
  rm(em_dt, em_summary)
  gc()
}

n_spp_before <- n_distinct(grid_sf$scientificName)
n_rows_before <- nrow(grid_sf)
grid_sf <- grid_sf |> dplyr::filter(!scientificName %in% introduced_final)
species_counts <- species_counts |> dplyr::filter(!scientificName %in% introduced_final)
message("Non-native filter: removed ",
        n_spp_before - n_distinct(grid_sf$scientificName), " introduced species, ",
        format(n_rows_before - nrow(grid_sf), big.mark = ","), " rows | ",
        format(n_distinct(grid_sf$scientificName), big.mark = ","), " species remain")

# ============================================================
# CELL-LEVEL FILTERING
# Three complementary filters are applied, then unioned:
#   1. Edge-buffer filter — endemism outlier candidates that
#      fall within direction-specific buffers around the Amazon
#      boundary (ecological transition zones).
#   2. City/urban filter — cells overlapping Natural Earth
#      urban extent polygons (cultivated / botanical garden risk).
#   3. Extreme outlier filter — cells with raw_endemism > 1000,
#      regardless of location (herbarium georeferencing artifacts).
# ============================================================

# --- Load Amazon boundary for edge-buffer calculations ---
amazon_shp_path <- paste0(dirname(dataDir),
                          "/data/amazon_shapefile/amazon_sensulatissimo_gmm_v1.shp")
amazon_shape <- st_read(amazon_shp_path, quiet = TRUE)
amazon_shape <- st_make_valid(amazon_shape)
amazon_shape <- st_simplify(amazon_shape, dTolerance = 0.01)
amazon_3857  <- st_transform(amazon_shape, 3857)

# --- One row per grid cell (reused by multiple filters) ---
grid_sf <- st_as_sf(grid_sf, sf_column_name = "geometry", crs = 3857)
grid_unique <- grid_sf |>
  dplyr::group_by(grid_id) |>
  dplyr::slice(1) |>
  dplyr::ungroup()

cell_width_m <- diff(st_bbox(grid_unique[1, ])[c("xmin", "xmax")])

# ---- 1. Edge-buffer filter ----
# Identify endemism outlier candidates with a broad threshold (Q95 + 10*IQR,
# union of raw and weighted endemism). Then keep only candidates that lie
# within direction-specific buffers around the Amazon boundary — the idea
# is that edge cells in ecological transitions (Andes, coast, cerrado) are
# most likely to be sampling artifacts rather than genuine hotspots.

iqr_multiplier <- 10

# Weighted endemism — cell-level values from Stage 0
cell_level_we <- species_endemism |> st_drop_geometry()
q95_we  <- quantile(cell_level_we$weighted_endemism, 0.95, na.rm = TRUE)
iqr_we  <- IQR(cell_level_we$weighted_endemism, na.rm = TRUE)
upper_cutoff_we <- q95_we + iqr_multiplier * iqr_we

# Raw endemism — cell-level values from Stage 0
raw_endemism <- readRDS(paste0(dataDir, "/raw_species_endemism2.rds"))
cell_level_re <- raw_endemism |> st_drop_geometry()
q95_re  <- quantile(cell_level_re$raw_endemism, 0.95, na.rm = TRUE)
iqr_re  <- IQR(cell_level_re$raw_endemism, na.rm = TRUE)
upper_cutoff_re <- q95_re + iqr_multiplier * iqr_re

message("Endemism outlier thresholds (Q95 + ", iqr_multiplier, "*IQR):")
message("  Weighted endemism: ", round(upper_cutoff_we, 2))
message("  Raw endemism: ", round(upper_cutoff_re, 2))

filtered_ids_we <- cell_level_we |>
  dplyr::filter(weighted_endemism > upper_cutoff_we) |>
  dplyr::pull(grid_id)
filtered_ids_re <- cell_level_re |>
  dplyr::filter(raw_endemism > upper_cutoff_re) |>
  dplyr::pull(grid_id)

endemism_candidate_ids <- union(filtered_ids_we, filtered_ids_re)
message("Endemism outlier candidates (union): ", length(endemism_candidate_ids))

# Compute distance from each candidate to the Amazon boundary and classify
# which edge it is nearest to.
candidate_cells <- grid_unique |>
  dplyr::filter(grid_id %in% endemism_candidate_ids) |>
  st_transform(3857)

amazon_boundary_3857 <- st_boundary(st_union(amazon_3857))
candidate_centroids <- st_centroid(candidate_cells)
dist_to_boundary <- as.numeric(
  st_distance(candidate_centroids, amazon_boundary_3857)
)
candidate_cells$dist_to_boundary_cells <- dist_to_boundary / cell_width_m

# Find nearest boundary point for each candidate to classify edge direction
candidate_centroids_ll <- st_transform(candidate_centroids, 4326)
amazon_boundary_ll <- st_boundary(st_union(st_transform(amazon_3857, 4326)))
nearest_lines <- st_nearest_points(candidate_centroids_ll, amazon_boundary_ll)
boundary_pts_coords <- st_coordinates(st_cast(nearest_lines, "POINT"))
# Every other point is the boundary endpoint (odd = centroid, even = boundary)
idx_bnd <- seq(2, nrow(boundary_pts_coords), by = 2)
bnd_lon <- boundary_pts_coords[idx_bnd, 1]
bnd_lat <- boundary_pts_coords[idx_bnd, 2]

# Classify edge: western Andes, northern/coastal, or southern (cerrado)
candidate_cells$edge <- dplyr::case_when(
  bnd_lon < -67                          ~ "west_andes",
  bnd_lat > 3 | (bnd_lon > -55 & bnd_lat > -2) ~ "north_coast",
  bnd_lat < -8                           ~ "south",
  TRUE                                   ~ "other"
)

# Direction-specific buffer thresholds (in cell widths):
#   West (Andes):   5 cells (~110 km) — steep Andes-to-lowland gradient
#   North (coast):  1 cell  (~22 km)  — tight coastal edge (Guianas)
#   South (cerrado):4 cells (~90 km)  — diffuse cerrado-Amazon ecotone
edge_buffer_thresholds <- c(west_andes = 5, north_coast = 1, south = 4, other = 2)
candidate_cells$buffer_threshold <- edge_buffer_thresholds[candidate_cells$edge]
candidate_cells$within_buffer <- candidate_cells$dist_to_boundary_cells <=
  candidate_cells$buffer_threshold

edge_buffer_ids <- candidate_cells$grid_id[candidate_cells$within_buffer]

message("Edge-buffer filter: ", length(edge_buffer_ids), " of ",
        length(endemism_candidate_ids), " candidates within buffer")

# ---- 2. City/urban filter ----
# Remove cells overlapping Natural Earth urban extent polygons.
urban_cache_path <- "output/urban_areas_amazon_3857.rds"

if (file.exists(urban_cache_path)) {
  message("Loading cached urban area polygons...")
  urban_crop <- readRDS(urban_cache_path)
} else {
  study_bbox <- st_bbox(st_transform(grid_unique, 4326))

  message("Downloading Natural Earth urban areas layer...")
  urban_sf <- ne_download(
    scale = 10, type = "urban_areas",
    category = "cultural", returnclass = "sf"
  )
  urban_sf <- st_make_valid(urban_sf)
  urban_crop <- st_crop(urban_sf, study_bbox)
  rm(urban_sf); gc()

  urban_crop <- st_transform(urban_crop, 3857)
  saveRDS(urban_crop, urban_cache_path)
  message("Urban polygons cached to ", urban_cache_path)
}
message("Urban polygons in study area: ", nrow(urban_crop))

grid_unique_3857 <- grid_unique |> st_transform(3857)
hits <- st_intersects(grid_unique_3857, urban_crop, sparse = TRUE)
grid_ids_cities <- unique(grid_unique_3857$grid_id[lengths(hits) > 0])
rm(urban_crop); gc()

message("Cells flagged — urban area polygons: ", length(grid_ids_cities))

# ---- 3. Extreme outlier filter ----
# Cells with raw_endemism > 1000 are removed regardless of location.
# These are herbarium georeferencing artifacts (e.g., specimens pinned
# to an institution's coordinates rather than collection sites).
extreme_endemism_threshold <- 1000
extreme_outlier_ids <- cell_level_re |>
  dplyr::filter(raw_endemism > extreme_endemism_threshold) |>
  dplyr::pull(grid_id)

message("Cells flagged — extreme outliers (raw_endemism > ",
        extreme_endemism_threshold, "): ", length(extreme_outlier_ids))

# ---- Combine all filters ----
bad_grid_ids <- unique(c(edge_buffer_ids, grid_ids_cities, extreme_outlier_ids))
message("Total unique cells to remove: ", length(bad_grid_ids),
        " (edge-buffer: ", length(edge_buffer_ids),
        ", cities: ", length(grid_ids_cities),
        ", extreme: ", length(extreme_outlier_ids), ")")

# Save removal log for reproducibility
removal_log <- dplyr::bind_rows(
  candidate_cells |>
    st_drop_geometry() |>
    dplyr::filter(within_buffer) |>
    dplyr::transmute(grid_id, reason = paste0("edge_buffer_", edge),
                     edge, dist_cells = round(dist_to_boundary_cells, 1)),
  tibble::tibble(grid_id = setdiff(grid_ids_cities, edge_buffer_ids),
                 reason = "city_urban", edge = NA_character_,
                 dist_cells = NA_real_),
  tibble::tibble(grid_id = extreme_outlier_ids,
                 reason = "extreme_outlier", edge = NA_character_,
                 dist_cells = NA_real_)
)
saveRDS(removal_log, "output/removed_cells_log.rds")
message("Removal log saved: output/removed_cells_log.rds")

species_counts_filt <- species_counts |>
  dplyr::filter(!grid_id %in% bad_grid_ids)

n_rows_removed  <- nrow(species_counts) - nrow(species_counts_filt)
n_cells_remaining <- n_distinct(species_counts_filt$grid_id)
n_spp_remaining  <- n_distinct(species_counts_filt$scientificName)
message("After cell removal: ",
        format(nrow(species_counts_filt), big.mark = ","), " rows (",
        format(n_rows_removed, big.mark = ","), " removed) | ",
        format(n_spp_remaining, big.mark = ","), " species | ",
        format(n_cells_remaining, big.mark = ","), " grid cells")

saveRDS(species_counts_filt, "output/species_counts4_filtered.rds")

# ============================================================
# RANGE-FILLING INTERPOLATION
# fillingRanges.R only defines functions when sourced (no data loading /
# computation runs automatically). Build grid-cell covariates and run the
# interpolation explicitly on species_counts_filt — the occurrence table
# being interpolated here — so sampling completeness is always calculated
# from these same (filtered) occurrences, not reused from an unrelated cache.
# ============================================================
message("--- Range-filling interpolation (fillingRanges.R) ---")
source("scripts/computeCompleteness.R")
source("scripts/fillingRanges.R")


species_counts_filt<-readRDS("output/species_counts4_filtered.rds")

grid_covariates <- build_grid_cell_covariates(species_counts_filt)
filled_ranges_all <- fill_all_species_ranges(species_counts_filt,
                                             grid_covariates$grid_unique_cells)

species_counts_interp <- filled_ranges_all
saveRDS(species_counts_interp, "output/species_c4_filt_interpolated.rds")

# Species with too few cells to build a concave hull are dropped by
# fill_species_range(); flag how many, since this affects raw endemism
# (which counts single-cell species).
n_spp_filt   <- n_distinct(species_counts_filt$scientificName)
n_spp_interp <- n_distinct(species_counts_interp$scientificName)
message("Interpolated occurrences saved: output/species_c4_filt_interpolated.rds")
message("  ", format(nrow(species_counts_interp), big.mark = ","), " (species, cell) presences | ",
        format(n_spp_interp, big.mark = ","), " species (",
        format(n_spp_filt - n_spp_interp, big.mark = ","),
        " dropped — too few cells to interpolate)")



## Recalculate endemism from the filtered species counts so values reflect
## the post-removal dataset (removing cells can change which species are
## single-cell endemics and alters range sizes for weighted endemism).

# Weighted endemism (recalculated)
weighted_endemism_filt <- species_counts_filt |>
  dplyr::group_by(scientificName) |>
  dplyr::mutate(range_size = n_distinct(grid_id)) |>
  dplyr::ungroup() |>
  dplyr::mutate(endemism_weight = 1 / range_size) |>
  dplyr::group_by(grid_id) |>
  dplyr::summarise(weighted_endemism = sum(endemism_weight, na.rm = TRUE),
                   .groups = "drop")

# Raw endemism (recalculated)
endemics_filt <- species_counts_filt |>
  dplyr::group_by(scientificName) |>
  dplyr::filter(n_distinct(grid_id) == 1) |>
  dplyr::ungroup()

raw_endemism_filt <- endemics_filt |>
  dplyr::group_by(grid_id) |>
  dplyr::summarise(raw_endemism = n_distinct(scientificName), .groups = "drop")

# Save individual filtered endemism layers
grid_geom_filt <- species_counts_filt |>
  dplyr::group_by(grid_id) |>
  dplyr::slice(1) |>
  dplyr::ungroup() |>
  dplyr::select(grid_id, geometry)

grid_sf_RE_filtered <- grid_geom_filt |>
  dplyr::left_join(st_drop_geometry(raw_endemism_filt), by = "grid_id")

grid_sf_SE_filtered <- grid_geom_filt |>
  dplyr::left_join(st_drop_geometry(weighted_endemism_filt), by = "grid_id")

saveRDS(grid_sf_RE_filtered, "output/filteredRaw_species_endemism.rds")
saveRDS(grid_sf_SE_filtered, "output/filteredWeighted_species_endemism.rds")

# Consolidate into single output: one row per grid cell with both metrics
allEndemismAndOccurenceFiltered <- grid_geom_filt |>
  dplyr::left_join(st_drop_geometry(weighted_endemism_filt), by = "grid_id") |>
  dplyr::left_join(st_drop_geometry(raw_endemism_filt), by = "grid_id")

saveRDS(allEndemismAndOccurenceFiltered, "output/FINAL_allEndemism_FullyFiltered.rds")

message("Recalculated endemism on filtered data:")
message("  Weighted endemism: ", nrow(weighted_endemism_filt), " cells, range ",
        round(min(weighted_endemism_filt$weighted_endemism, na.rm = TRUE), 2), " - ",
        round(max(weighted_endemism_filt$weighted_endemism, na.rm = TRUE), 2))
message("  Raw endemism: ", nrow(raw_endemism_filt), " cells with endemics, ",
        n_distinct(endemics_filt$scientificName), " single-cell species")
message("Final output: ", nrow(allEndemismAndOccurenceFiltered), " grid cells in ",
        "FINAL_allEndemism_FullyFiltered.rds")
message("--- Stage 1 complete ---")

####plot each endemism dataset

###weighted plot
pWeight <- ggplot(grid_sf_SE_filtered) +
  geom_sf(aes(color = weighted_endemism), size = 0.6) +
  scale_color_viridis_c(option = "plasma") +
  theme_minimal()

ggsave(
  filename = "output/weighted_endemism.svg",
  plot = pWeight,
  device = "svg",
  width = 8, height = 6, units = "in"
)

###raw plot
pRaw <- ggplot(grid_sf_RE_filtered) +
  geom_sf(aes(color = raw_endemism), size = 0.6) +
  scale_color_viridis_c(option = "plasma") +
  theme_minimal()

ggsave(
  filename = "output/raw_endemism.svg",
  plot = pRaw,
  device = "svg",
  width = 8, height = 6, units = "in"
)

# ============================================================
# INTERPOLATED ENDEMISM
# Repeat the endemism recalculation + maps on the range-filled
# (interpolated) occurrences so the two can be compared side by side.
# Mirrors the filtered block above with `_interp` variables and
# distinct output filenames.
# ============================================================
message("--- Interpolated endemism recalculation + maps ---")

# Weighted endemism (interpolated)
weighted_endemism_interp <- species_counts_interp |>
  dplyr::group_by(scientificName) |>
  dplyr::mutate(range_size = n_distinct(grid_id)) |>
  dplyr::ungroup() |>
  dplyr::mutate(endemism_weight = 1 / range_size) |>
  dplyr::group_by(grid_id) |>
  dplyr::summarise(weighted_endemism = sum(endemism_weight, na.rm = TRUE),
                   .groups = "drop")

# Raw endemism (interpolated)
endemics_interp <- species_counts_interp |>
  dplyr::group_by(scientificName) |>
  dplyr::filter(n_distinct(grid_id) == 1) |>
  dplyr::ungroup()

raw_endemism_interp <- endemics_interp |>
  dplyr::group_by(grid_id) |>
  dplyr::summarise(raw_endemism = n_distinct(scientificName), .groups = "drop")

# One row per grid cell (interpolated)
grid_geom_interp <- species_counts_interp |>
  dplyr::group_by(grid_id) |>
  dplyr::slice(1) |>
  dplyr::ungroup() |>
  dplyr::select(grid_id, geometry)

grid_sf_RE_interp <- grid_geom_interp |>
  dplyr::left_join(st_drop_geometry(raw_endemism_interp), by = "grid_id")

grid_sf_SE_interp <- grid_geom_interp |>
  dplyr::left_join(st_drop_geometry(weighted_endemism_interp), by = "grid_id")

saveRDS(grid_sf_RE_interp, "output/interp_filteredRaw_species_endemism.rds")
saveRDS(grid_sf_SE_interp, "output/interp_filteredWeighted_species_endemism.rds")

# Consolidate into single output: one row per grid cell with both metrics
allEndemismAndOccurenceInterp <- grid_geom_interp |>
  dplyr::left_join(st_drop_geometry(weighted_endemism_interp), by = "grid_id") |>
  dplyr::left_join(st_drop_geometry(raw_endemism_interp), by = "grid_id")

saveRDS(allEndemismAndOccurenceInterp, "output/INTERP_allEndemism_FullyFiltered.rds")

message("Recalculated endemism on interpolated data:")
message("  Weighted endemism: ", nrow(weighted_endemism_interp), " cells, range ",
        round(min(weighted_endemism_interp$weighted_endemism, na.rm = TRUE), 2), " - ",
        round(max(weighted_endemism_interp$weighted_endemism, na.rm = TRUE), 2))
message("  Raw endemism: ", nrow(raw_endemism_interp), " cells with endemics, ",
        n_distinct(endemics_interp$scientificName), " single-cell species")

####plot each interpolated endemism dataset

###weighted plot (interpolated)
pWeight_interp <- ggplot(grid_sf_SE_interp) +
  geom_sf(aes(color = weighted_endemism), size = 0.6) +
  scale_color_viridis_c(option = "plasma") +
  theme_minimal()

ggsave(
  filename = "output/interp_weighted_endemism.svg",
  plot = pWeight_interp,
  device = "svg",
  width = 8, height = 6, units = "in"
)

###raw plot (interpolated)
pRaw_interp <- ggplot(grid_sf_RE_interp) +
  geom_sf(aes(color = raw_endemism), size = 0.6) +
  scale_color_viridis_c(option = "plasma") +
  theme_minimal()

ggsave(
  filename = "output/interp_raw_endemism.svg",
  plot = pRaw_interp,
  device = "svg",
  width = 8, height = 6, units = "in"
)

