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
  sf::st_drop_geometry() %>%
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

#sum(plant_distribution$n)

#plant_grid_counts

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

#### Cell filter report: maps of ALL original cells (animals and plants) ####
# Colours cells by Chao1 completeness, richness and occurrence count, plus the
# plant-minus-animal completeness difference. Retained animal cells are outlined.
if (exists(".amani_filter_summary", envir = .GlobalEnv) &&
    !is.null(.amani_filter_summary$cell_table)) {

  animal_cells <- .amani_filter_summary$cell_table

  # Plant per-cell stats for every plant cell (Chao1, same formula as the
  # optional plant completeness filter below). Occurrence counts live in
  # `count` for the raw plant data and `n_records` for the filled/interpolated
  # data, so pick whichever column exists.
  #
  # NOTE (downstream usage): when fill_plants = TRUE this table is computed
  # from the FILLED/interpolated plant data, whose inflated record counts
  # leave ~no singletons per cell, so plant_completeness collapses to exactly
  # 1 in every cell. Any DOWNSTREAM use of plant completeness (figures,
  # summaries, the delta_completeness panel below) that is meant to measure
  # actual sampling quality should instead be calculated from the RAW
  # filtered counts (1-FilteringPlantData/output/species_counts4_filtered.rds,
  # same formula, keyed by grid_id) — NOT from the filled value computed here.
  plant_tbl <- species_counts |> sf::st_drop_geometry()
  cnt_col <- if ("count" %in% names(plant_tbl)) "count" else "n_records"
  plant_cells <- plant_tbl |>
    dplyr::group_by(grid_id) |>
    dplyr::summarise(
      plant_richness    = dplyr::n_distinct(scientificName),
      plant_occurrences = sum(.data[[cnt_col]]),
      f1 = sum(.data[[cnt_col]] == 1),
      f2 = sum(.data[[cnt_col]] == 2),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      plant_chao1 = plant_richness + dplyr::if_else(
        f2 > 0, f1^2 / (2 * f2), f1 * (f1 - 1) / 2
      ),
      plant_completeness = plant_richness / plant_chao1
    ) |>
    dplyr::select(-f1, -f2, -plant_chao1)

  # Polygon geometry for plant cells; use grid polygons for animal cells too
  # if the stored geometry is not polygonal
  plant_cells_sf <- grid_sf |>
    dplyr::select(grid_id) |>
    dplyr::inner_join(plant_cells, by = "grid_id") |>
    sf::st_transform(sf::st_crs(animal_cells))

  if (!all(sf::st_geometry_type(animal_cells) %in% c("POLYGON", "MULTIPOLYGON"))) {
    poly <- sf::st_geometry(grid_sf)[match(animal_cells$grid_id, grid_sf$grid_id)]
    animal_cells <- sf::st_set_geometry(animal_cells,
                                        sf::st_transform(poly, sf::st_crs(animal_cells)))
  }

  # Plant-minus-animal completeness: negative = plants less complete
  delta_cells <- animal_cells |>
    dplyr::inner_join(sf::st_drop_geometry(plant_cells), by = "grid_id") |>
    dplyr::mutate(delta_completeness = plant_completeness - animal_completeness)

  saveRDS(list(animal = animal_cells, plant = plant_cells_sf, delta = delta_cells),
          file.path(plot_data_dir, "plotdata_cell_filter_report.rds"))

  outline_scale <- scale_color_manual(
    values = c(`TRUE` = "black", `FALSE` = "grey70"),
    labels = c(`TRUE` = "Retained", `FALSE` = "Filtered out"),
    name = "Animal cell"
  )

  map_pair <- function(an_var, pl_var, title, fill_name, trans = "identity") {
    pa <- ggplot(animal_cells) +
      geom_sf(aes(fill = .data[[an_var]], color = retained), linewidth = 0.4) +
      scale_fill_viridis_c(name = fill_name, trans = trans) +
      outline_scale +
      labs(title = paste("Animals:", title)) +
      theme_minimal()
    pp <- ggplot(plant_cells_sf) +
      geom_sf(aes(fill = .data[[pl_var]]), color = NA) +
      scale_fill_viridis_c(name = fill_name, trans = trans) +
      labs(title = paste("Plants:", title)) +
      theme_minimal()
    patchwork::wrap_plots(pa, pp, nrow = 1)
  }

  p_delta <- ggplot(delta_cells) +
    geom_sf(aes(fill = delta_completeness, color = retained), linewidth = 0.4) +
    scale_fill_gradient2(name = "Plant - animal\ncompleteness",
                         low = "firebrick", mid = "white", high = "steelblue",
                         midpoint = 0,
                         limits = c(-1, 1) * max(abs(delta_cells$delta_completeness), na.rm = TRUE)) +
    outline_scale +
    labs(title = "Delta Chao1 completeness (plant - animal), animal grid cells") +
    theme_minimal()

  ts_cfr <- format(Sys.time(), "%Y%m%d_%H%M%S")
  pdf(file.path(output_folder, paste0("cell_filter_report_maps_", ts_cfr, ".pdf")),
      width = 14, height = 7)
  print(map_pair("animal_completeness", "plant_completeness",
                 "Chao1 completeness", "Completeness"))
  print(map_pair("animal_richness", "plant_richness",
                 "species richness", "Richness"))
  print(map_pair("animal_occurrences", "plant_occurrences",
                 "number of occurrences (log10)", "Occurrences", trans = "log10"))
  print(p_delta)
  dev.off()
  message("Cell filter report maps written to: ", output_folder)
}

# Filter species_counts directly by animal grid IDs — no need to join grid_sf
# geometry or st_crop since species_counts already has grid_id and we
# immediately drop geometry anyway. This avoids building ~1.5 GB intermediates.
species_counts_joined  <- NULL  # keep downstream rm() calls safe
species_counts_overlap <- NULL
filtered_species_counts <- NULL

plant_grid_counts <- as.data.table(sf::st_drop_geometry(species_counts))[
  grid_id %in% animal_grid$grid_id,
  .(num_grid_cells = uniqueN(grid_id)),
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
print("Done with that heavy plot...")

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
#
# The grid only has cells with plant records, so ~1/3 of the cells inside the
# basin are missing (a sampling artefact, not habitat) and fragment simulated
# ranges. Add them back so contiguous range growth sees a gap-free grid.
# Set .amani_fillGridHoles <- FALSE to use the original holey grid.
source(file.path(project_root, "2-ModellingAnimals/scripts/utils/fillGridHoles.R"))
if (isTRUE(get0(".amani_fillGridHoles", envir = .GlobalEnv, ifnotfound = TRUE))) {
  amazon_boundary <- get0("amazon_shp", envir = .GlobalEnv, ifnotfound = NULL)
  if (is.null(amazon_boundary)) {
    amazon_boundary <- sf::st_read(
      file.path(dirname(getwd()),
                "0-CuratingPlantData/data/amazon_shapefile/amazon_sensulatissimo_gmm_v1.shp"),
      quiet = TRUE
    )
  }
  grid_sf <- fill_grid_holes(grid_sf, amazon_boundary, grid_size_m)
} else {
  grid_sf$grid_filled <- FALSE
}
available_grid_ids <- grid_sf$grid_id[
  grid_sf$grid_id %in% c(grid_sf_clipped$grid_id, grid_sf$grid_id[grid_sf$grid_filled])
]
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


# Smoothed frontier weights for contiguous range growth. Raw plant counts are
# sharply bimodal, so weighting the frontier by them produces lobed, gappy
# ranges; a Gaussian-smoothed surface keeps the regional sampling signal
# without the cell-level spikes. Start-cell selection still uses the raw
# sampling_weights. Plant-only, so cached in output/shared/ across organisms.
# Set .amani_frontierBandwidth <- 0 to restore the old raw-weight frontier.
source(file.path(project_root, "2-ModellingAnimals/scripts/utils/smoothFrontierWeights.R"))

# Filled (previously missing) cells have no plant records, so the code above
# gave them the floor weight. Replace it with a local interpolation from
# sampled neighbours so they are not treated as unsuitable habitat.
filled_cells <- grid_sf$grid_filled[match(available_grid_ids, grid_sf$grid_id)]
sampling_weights <- interpolate_filled_weights(
  weights = sampling_weights, grid_ids = available_grid_ids,
  grid_sf = grid_sf, grid_size_m = grid_size_m, filled = filled_cells
)

frontier_bandwidth <- get0(".amani_frontierBandwidth", envir = .GlobalEnv, ifnotfound = 2)
frontier_weights <- if (is.numeric(frontier_bandwidth) && frontier_bandwidth > 0) {
  load_or_build_frontier_weights(
    raw_weights     = sampling_weights,
    grid_ids        = available_grid_ids,
    grid_sf         = grid_sf,
    grid_size_m     = grid_size_m,
    fill_plants     = isTRUE(get0(".amani_fill_plants", envir = .GlobalEnv, ifnotfound = FALSE)),
    bandwidth_cells = frontier_bandwidth,
    redo            = isTRUE(get0("redoeverything", envir = .GlobalEnv, ifnotfound = FALSE))
  )
} else {
  NULL
}


source(file.path(project_root, "2-ModellingAnimals/scripts/utils/SimulateScaledRichness.R"))


# Range-simulation switch for the top-level calls below. Set by run_animal_model()
# via .amani_simulateRanges; defaults to TRUE for standalone runs of this script.
# When FALSE, species ranges are not grown (only species counts are produced),
# which is all that gamma richness and the CI bounds require.
simulate_ranges <- if (exists(".amani_simulateRanges", envir = .GlobalEnv)) {
  # TRUE, "BestOnly", and "SkipUpper" all enable range simulation for the

  # best estimate; the upper/lower CI runs are selectively skipped downstream.
  isTRUE(.amani_simulateRanges) ||
    identical(.amani_simulateRanges, "BestOnly") ||
    identical(.amani_simulateRanges, "SkipUpper")
} else {
  TRUE
}


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
# NOTE (future): When _amani_useSamplingExtension = TRUE, simulateScaledSpeciesRichness()
# now returns a named vector c(gamma, gamma_lwr, gamma_upr) derived from the
# extend_richness_from_sampling() prediction interval. Consider replacing the
# separate upper/lower simulation runs below (upperEstSpRich / lowerEstSpRich,
# which each run the full simulation again with fit_predUpper / fit_predLower)
# with these model-derived bounds, which are cheaper to obtain and may represent
# a conceptually distinct source of uncertainty (sampling model error vs.
# extent-scaling curve uncertainty). Decide which error source is appropriate
# to report before making this change.
if (spatial == "explicit") {
    # spatially explicit simulation (uses empirical / distance-weighted structure)
  gammaRichnessSPRIch <- simulateScaledSpeciesRichness(
    fit_pred,
    extent_scaling_factor = extent_sf_here,
    sampling_weights = sampling_weights,
    frontier_weights = frontier_weights,
    simulate_ranges = simulate_ranges 
  )
  
} else {
  
  # non-spatial / neutral simulation (uniform sampling weights)
  
  gammaRichnessSPRIch <- simulateScaledSpeciesRichness(
    fit_pred,
    extent_scaling_factor = extent_sf_here,
    sampling_weights = sampling_weights_neutral,
    simulate_ranges = simulate_ranges
  )
  
}

# Determine which error-bound mode to apply throughout this section.
# "full"   -> run separate upper/lower simulations (extent-curve bootstrap CI);
#             report gamma_upr of upperEstSpRich and gamma_lwr of lowerEstSpRich.
# "narrow" -> skip upper/lower simulations; report gamma_upr / gamma_lwr from
#             the best-estimate run only (sampling-model prediction interval).
gammaErrorType <- get0(".amani_gammaErrorType", envir = .GlobalEnv, ifnotfound = "full")
message("  gammaErrorType = '", gammaErrorType, "'")

# Extract narrow bounds from the best-estimate simulation return value.
# Only available when .amani_useSamplingExtension = TRUE, which causes
# simulateScaledSpeciesRichness() to return a named c(gamma, gamma_lwr, gamma_upr).
# Falls back to NA when the sampling extension is off (plain-integer return).
if (is.numeric(gammaRichnessSPRIch) && !is.null(names(gammaRichnessSPRIch))) {
  gamma_lwr_narrow <- unname(gammaRichnessSPRIch["gamma_lwr"])
  gamma_upr_narrow <- unname(gammaRichnessSPRIch["gamma_upr"])
} else {
  gamma_lwr_narrow <- NA_real_
  gamma_upr_narrow <- NA_real_
}

# Capture the sampling-extension-model result from THIS (best-estimate) run only.
# simulateScaledSpeciesRichness() also invokes the sampling extension model
# internally for the upper/lower bound re-simulations further below, and each
# call overwrites the global `gammaRichnessFromSamplingModel` via <<-. Without
# snapshotting it here, anything read from that global later in the script
# (e.g. in the report) would silently reflect the lower-bound run's fit curve
# instead of the best estimate's, which is the bug reported downstream.
gammaRichnessFromSamplingModel_best <- if (
  isTRUE(get0(".amani_useSamplingExtension", envir = .GlobalEnv, ifnotfound = FALSE)) &&
  exists("gammaRichnessFromSamplingModel", envir = .GlobalEnv)
) {
  gammaRichnessFromSamplingModel
} else {
  NULL
}

# The simulation now exposes the combined occurrence table directly as a compact
# data.table (simulated_occurrence_dt), replacing the old list-of-data.frames
# `sim_records` that was rbind-ed here (which duplicated the whole table in RAM).
simulated_occurrence_df <- as.data.frame(simulated_occurrence_dt)
head(simulated_occurrence_df)

###export this result, as it is the main output of this model, and can be used downstream to estimate other things
print("Joining simulated species with spatial geometry for exportation...")

### Speed-up: build geometry lookup from compact grid_sf (10k rows) instead
### of deduplicating the full species_counts table (~16.8M rows).
needed_ids <- unique(simulated_occurrence_df$grid_id)

if (all(needed_ids %in% grid_sf$grid_id)) {
  geom_tbl <- as.data.table(grid_sf)[, .(grid_id, geometry)]
} else {
  warning(
    "grid_sf does not cover all grid_ids in simulated_occurrence_df; ",
    "falling back to slower lookup derived from species_counts."
  )
  first_idx <- which(!duplicated(species_counts$grid_id))
  geom_tbl <- data.table::data.table(
    grid_id  = species_counts$grid_id[first_idx],
    geometry = sf::st_geometry(species_counts)[first_idx]
  )
}


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


###########Get CI's for the extent-based gamma diversity estimator surrogate model#################
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
  lower = apply(boot_preds, 1, function(x) quantile(x, 0.025, na.rm = TRUE)),
  upper = apply(boot_preds, 1, function(x) quantile(x, 0.975, na.rm = TRUE)),
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
    title = "Log-normal Fit to Confidence Interval Bounds",
    x = "Number of Grid Cells",
    y = "95% CI",
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

if ( !.amani_estimateBounds || gammaErrorType == "narrow"){
  # Skip the full upper-bound simulation: either bounds are disabled globally,
  # or "narrow" mode uses the sampling-model prediction interval instead.
  upperEstSpRich <- NA_real_
  
} else if (identical(.amani_simulateRanges, "BestOnly") || identical(.amani_simulateRanges, "SkipUpper")) {
  # BestOnly: ranges were only simulated for the best estimate; skip both CI re-runs.
  # SkipUpper: skip the upper-bound re-simulation only (lower still runs below).
  # Both fall back to the sampling-model prediction upper interval.
  upperEstSpRich <- gamma_upr_narrow
  message("  ", .amani_simulateRanges, ": skipping upper-bound re-simulation; upperEstSpRich <- gamma_upr_narrow = ", gamma_upr_narrow)
  
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
      sampling_weights = sampling_weights,
      frontier_weights = frontier_weights,
      simulate_ranges = simulate_ranges   # only the species count is needed for the CI bound
    )
    
  } else {
    
    # non-spatial / neutral simulation (uniform sampling weights)
    
    upperEstSpRich <- simulateScaledSpeciesRichness(
      fit_predUpper,
      extent_scaling_factor = extent_sf_here,
      sampling_weights = sampling_weights_neutral,
      simulate_ranges = simulate_ranges   # only the species count is needed for the CI bound
    )
    
  }
  
  
  
} else {
  
  extent_sf_here <-customScale
  upperEstSpRich <-999999999
  
}
}





plant_fit_max <- max(fit_predLower$fitted)
###allows specifying a custom scaling factor (mostly for troubleshooting)


if ( !.amani_estimateBounds || gammaErrorType == "narrow"){
  # Skip the full lower-bound simulation: either bounds are disabled globally,
  # or "narrow" mode uses the sampling-model prediction interval instead.
  lowerEstSpRich <- NA_real_
  
} else if (identical(.amani_simulateRanges, "BestOnly")) {
  # Range simulation was run for the best estimate only; skip the lower-bound
  # re-simulation and fall back to the sampling-model prediction lower interval.
  lowerEstSpRich <- gamma_lwr_narrow
  message("  BestOnly: skipping lower-bound re-simulation; lowerEstSpRich <- gamma_lwr_narrow = ", gamma_lwr_narrow)
  
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
      sampling_weights = sampling_weights,
      frontier_weights = frontier_weights,
      simulate_ranges = simulate_ranges    # only the species count is needed for the CI bound
    )
    
  } else {
    
    # non-spatial / neutral simulation (uniform sampling weights)
    
    lowerEstSpRich <- simulateScaledSpeciesRichness(
      fit_predLower,
      extent_scaling_factor = extent_sf_here,
      sampling_weights = sampling_weights_neutral,
      simulate_ranges = simulate_ranges    # only the species count is needed for the CI bound
    )
    
  }
  
  
} else {
  
  extent_sf_here <-customScale
  lowerEstSpRich<-0
  
}}

# Surrogate-only bounds: total species count from the plant CI curves
# *without* the sampling-extension correction.  Uses column position [[2]]
# because the upper frame is renamed "n" on line 1170 while the lower frame
# may still be called "fitted" if bounds were skipped.
partialGamma_upr <- round(sum(fit_predUpper[[2]] * (animal_max / max(fit_predUpper[[2]]))))
partialGamma_lwr <- round(sum(fit_predLower[[2]] * (animal_max / max(fit_predLower[[2]]))))

# Resolve the final upper/lower bounds to report, based on gammaErrorType.
# "full"   -> gamma_upr of the upperEstSpRich run  +  gamma_lwr of the lowerEstSpRich run
#             (widest CI: both extent-curve bootstrap uncertainty AND sampling-model PI)
# "narrow" -> gamma_upr and gamma_lwr from the best-estimate run only
#             (cheaper; captures sampling-model prediction interval without re-running simulation)
if (gammaErrorType == "full") {
  gamma_upr_report <- if (is.numeric(upperEstSpRich) && !is.null(names(upperEstSpRich))) {
    unname(upperEstSpRich["gamma_upr"])
  } else {
    upperEstSpRich   # plain-integer fallback (sampling extension off)
  }
  gamma_lwr_report <- if (is.numeric(lowerEstSpRich) && !is.null(names(lowerEstSpRich))) {
    unname(lowerEstSpRich["gamma_lwr"])
  } else {
    lowerEstSpRich   # plain-integer fallback
  }
} else {
  # "narrow": pull bounds directly from best-estimate simulation
  gamma_upr_report <- gamma_upr_narrow
  gamma_lwr_report <- gamma_lwr_narrow
}

# "Best estimate from plant surrogate model" should be the surrogate-only
# figure (partialGamma), i.e. before the sampling-extension correction is
# applied — not gammaRichnessEst, which (when the extension is active) is the
# stochastic realization of the *corrected* (post-extension) distribution.
# Using gammaRichnessEst here mislabeled the corrected estimate as the
# surrogate-only one, and — because it was compared against a stale/overwritten
# gammaRichnessFromSamplingModel (see snapshot above) — could appear
# inconsistent with the reported bounds.
gammaRichnessEst_report <- if (!is.null(gammaRichnessFromSamplingModel_best)) {
  gammaRichnessFromSamplingModel_best$partialGamma
} else {
  gammaRichnessEst
}

output_text <- paste0(
  "Gamma Richness Estimate (spatial extent scaling)\n",
  "================================================\n\n",
  "Fit parameters:    [", paste(fit_params, collapse = ", "), "]\n",
  "Plant range filling...",  if (isTRUE(get0(".amani_fill_plants", envir = .GlobalEnv, ifnotfound = FALSE))){"was used."}else{"was not used."},
  "Amazonian Animal gamma richness calculated from: ", if (isTRUE(get0(".amani_useSamplingExtension", envir = .GlobalEnv, ifnotfound = FALSE))) {"plant surrogate model + animal sampling extension model.\n"}else{"plant surrogate model\n"},
  "Extent scaling factor: ", round(animal_max / max(fit_pred$n), 6), "\n",
  "delta_meanlog:     ", round(.amani_delta_meanlog, 6),
    ifelse("meanlog" %in% fit_params, " (active)", " (off)"), "\n",
  "ratio_sdlog:       ", round(.amani_ratio_sdlog, 6),
    ifelse("sdlog" %in% fit_params, " (active)", " (off)"), "\n",
  "Best estimate from plant surrogate model:     ", gammaRichnessEst_report, " species\n",
  "Upper bound (plant CI):  ", partialGamma_upr, " species\n",
  "Lower bound (plant CI):  ", partialGamma_lwr, " species\n\n",
  
  if (!is.null(gammaRichnessFromSamplingModel_best)) {
    smp <- gammaRichnessFromSamplingModel_best
    paste0(
      "Expected proportional richness shortfall of surrogate model: ", smp$signed_prop_error, ". \n",
      "Best estimate from total model (surrogate + extension if used):     ", smp$fullGamma, " species. \n",
      "Upper bound:       ", smp$fullGamma_upr, " species\n",
      "Lower bound:       ", smp$fullGamma_lwr, " species\n\n"
    )
  } else {
    "Sampling extension model was not used.\n\n"
  },
  
  
  "Description of surrogate model: Projected total animal species richness of the Amazon ",
  "based on a log-normal distribution of plant spatial extents rescaled ",
  "to match observed animal data. ",
  if (gammaErrorType == "full") {
    "Bounds (full mode): re-ran simulation on bootstrap CI curves of the plant extent model (95%/105% resampling).\n\n"
  } else {
    "Bounds (narrow mode): gamma_upr and gamma_lwr from the sampling-model prediction interval of the best-estimate run only.\n\n"
  },
  "Shape-shifting: When fit_params includes 'meanlog' or 'sdlog', the ",
  "empirical plant range-size distribution is warped in log-space to match ",
  "the animal data's Amazonian spatial extent shape before amplitude scaling.",
  
  "Description of animal sampling extension model: Surrogate model shows predictable richness shortfalls in simulation, based mostly on incomplete sampling of species (also endemism of the sites being sampled and sampling evenness across sites). The surrogate model is a multi-dimensional linear model (parameterized in simulation) that takes the predicted gamma of the surrogate model and extends it to the predicted actual gamma. The particular model used varies depending on if the fill_plants option was used."
)

timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
output_path <- file.path(output_folder, paste0("gamma_richnessfromSPRICH.", timestamp, ".txt"))
writeLines(output_text, con = output_path)

#### Sampling extension model parameter report ####
# Written only when .amani_useSamplingExtension is TRUE and the model ran
# (i.e., gammaRichnessFromSamplingModel was populated by simulateScaledSpeciesRichness).
if (!is.null(gammaRichnessFromSamplingModel_best)) {

  se <- gammaRichnessFromSamplingModel_best  # shorthand; snapshot from the best-estimate run

  sampling_ext_report <- paste0(
    "=======================================================\n",
    " Sampling Extension Model — Parameter Report\n",
    " Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n",
    "=======================================================\n\n",

    "--- Model configuration ---\n",
    "  fill_Plants:       ", se$fill_Plants, "\n",
    "    (selects which sub-model was used: plant ranges filled vs. empirical)\n\n",

    "--- Input predictors (supplied to the linear model) ---\n",
    "  comp_spread:            ", round(se$comp_spread, 6),
    "  (max - min per-cell ", se$comp_completeness_method, " completeness)\n",
    "  comp_mean:              ", round(se$comp_mean, 6),
    "  (mean per-cell ", se$comp_completeness_method, " completeness)\n",
    "  prop_species (proxy):   ", round(se$prop_species_used, 6),
    "  (", se$prop_species_method, " regional completeness across grid cells; proxy for prop. of true species observed)\n",
    "    Chao2 inputs: S_obs = ", se$prop_species_chao2_Sobs,
    ", Q1 = ", se$prop_species_chao2_Q1,
    ", Q2 = ", se$prop_species_chao2_Q2,
    ", S_Chao2 = ", round(se$prop_species_chao2_est, 2), "\n",
    "  comp_sd:                ", round(se$comp_sd, 6),
    "  (SD of per-cell ", se$comp_completeness_method, " completeness)\n",
    "  n_grid_cells (raw):     ", se$n_grid_cells_raw,
    "  (number of occupied grid cells; log-transformed internally)\n",
    "  n_grid_cells (log):     ", round(log(se$n_grid_cells_raw), 6), "\n",
    "  sample_mean_range_size: ", round(se$sample_mean_range_size, 6),
    "  (weighted mean range size across animal morphotypes, in grid cells)\n",
    "  delta_meanlog:          ", round(se$delta_meanlog_used, 6),
    "  (animal - plant log-normal meanlog shift; ", ifelse("meanlog" %in% fit_params, "active", "off — using plant shape"), ")\n",
    "  ratio_sdlog:            ", round(se$ratio_sdlog_used, 6),
    "  (animal / plant log-normal sdlog ratio; ", ifelse("sdlog" %in% fit_params, "active", "off — using plant shape"), ")\n\n",

    "--- Model output ---\n",
    "  partialGamma:           ", se$partialGamma,
    "  (total species from plant surrogate model, before correction)\n",
    "  abs_prop_error:         ", round(se$signed_prop_error, 6),
    "  (predicted proportional richness shortfall of surrogate model)\n",
    "  correction_factor:      ", round(se$correction_factor, 6),
    "  (= 1 + abs_prop_error; multiplied uniformly across range-size bins)\n",
    "  fullGamma (best):       ", round(se$fullGamma, 2),
    "  species  (= partialGamma * correction_factor)\n",
    "  fullGamma (lower 95%):  ", round(se$fullGamma_lwr, 2), " species\n",
    "  fullGamma (upper 95%):  ", round(se$fullGamma_upr, 2), " species\n\n",

    "--- Interpretation ---\n",
    "  The linear model predicts that this dataset's surrogate gamma underestimates\n",
    "  true richness by ~", round(se$signed_prop_error * 100, 1), "% given the observed sampling\n",
    "  characteristics. Applying correction_factor = ", round(se$correction_factor, 4),
    " to each range-size bin\n",
    "  before simulation yields the corrected fullGamma above.\n\n",

    "=======================================================\n"
  )

  timestamp_se <- format(Sys.time(), "%Y%m%d_%H%M%S")
  se_report_path <- file.path(output_folder,
                              paste0("sampling_extension_report_", timestamp_se, ".txt"))
  writeLines(sampling_ext_report, con = se_report_path)
  message("Sampling extension report written to: ", se_report_path)
}

#### Optional: filter plant cells by Chao1 sampling completeness ########
# When .amani_filter_plant_completeness is set, only plant grid cells with
# Chao1 completeness >= that threshold are used in the correlation analyses
# below. This reduces sample size but strengthens the signal by removing
# poorly-sampled plant cells. The filter applies ONLY to the correlations,
# not to the rest of the pipeline.
#
# NOTE: completeness here is computed from the RAW filtered counts
# (species_counts4_filtered.rds), not species_counts, because species_counts is
# the FILLED plant data when fill_plants = TRUE (completeness would be ~1).
# The plant_completeness column in the cell-filter report above still uses
# species_counts.

#NOTE: Should this be incorporated earlier in the workflow? 


filter_plant_comp <- if (exists(".amani_filter_plant_completeness", envir = .GlobalEnv)) {
  .amani_filter_plant_completeness
} else {
  NULL
}

if (!is.null(filter_plant_comp) && is.numeric(filter_plant_comp)) {
  message("Computing plant Chao1 completeness for correlation filtering (threshold = ",
          filter_plant_comp, ")...")

  # Occurrence counts live in `count` (raw plant data) or `n_records`
  # (filled/interpolated data) — pick whichever exists.
  # Always use RAW (unfilled) counts: with fill_plants = TRUE, species_counts
  # is interpolated and completeness would be ~1 everywhere.
  raw_plant_counts <- readRDS(file.path(dataDir, "species_counts4_filtered.rds"))
  plant_tbl2 <- raw_plant_counts |> sf::st_drop_geometry()
  rm(raw_plant_counts)
  cnt_col2 <- if ("count" %in% names(plant_tbl2)) "count" else "n_records"
  plant_chao1_stats <- plant_tbl2 |>
    dplyr::group_by(grid_id) |>
    dplyr::summarise(
      plant_n_occ     = sum(.data[[cnt_col2]]),
      plant_n_species = dplyr::n_distinct(scientificName),
      f1 = sum(.data[[cnt_col2]] == 1),
      f2 = sum(.data[[cnt_col2]] == 2),
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

endemism_model_ok <- nrow(endemism_joined) >= 3
if (!endemism_model_ok) {
  warning("Only ", nrow(endemism_joined),
          " cells available for endemism correlation — model will be degenerate.",
          call. = FALSE)
}
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
  "Animal cell filter:            ", ifelse(n_animal_cells < sum(!grid_sf$grid_filled), "ON", "OFF"), "\n",
  "  Animal cells used:           ", n_animal_cells, "\n",
  "Plant completeness filter:     ",
    ifelse(is.null(filter_plant_comp), "OFF", paste0(">= ", filter_plant_comp)), "\n",
  "  Applied to endemism corr:    ", ifelse(plant_filter_used_endemism, "YES", "NO (< 15 cells or OFF)"), "\n",
  "  Applied to richness corr:    ", ifelse(plant_filter_used_richness, "YES", "NO (< 15 cells or OFF)"), "\n\n",
  "--- Endemism correlation (animal endemism ~ plant raw endemism) ---\n",
  "  Cells:       ", n_endemism_cells, "\n",
  if (endemism_model_ok && !is.null(endemism_summ$fstatistic)) {
    paste0(
      "  R²:          ", round(endemism_summ$r.squared, 4), "\n",
      "  Adj. R²:     ", round(endemism_summ$adj.r.squared, 4), "\n",
      "  F-statistic: ", round(endemism_summ$fstatistic[1], 2),
      " (df = ", endemism_summ$fstatistic[2], ", ", endemism_summ$fstatistic[3], ")\n",
      "  p-value:     ", format.pval(pf(endemism_summ$fstatistic[1],
                                         endemism_summ$fstatistic[2],
                                         endemism_summ$fstatistic[3],
                                         lower.tail = FALSE), digits = 4), "\n",
      "  Slope:       ", round(coef(endemism_model)[2], 4),
      " (SE = ", round(endemism_summ$coefficients[2, 2], 4), ")\n\n"
    )
  } else {
    "  ** Insufficient data for reliable endemism model (< 3 cells) **\n\n"
  },
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

