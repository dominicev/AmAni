{
  library(knitr)
  library(rstudioapi)
  
  if (exists(".amani_wd")) {
    setwd(.amani_wd)
  } else {
    # 1. Get the directory of the currently active script
    script_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
    
    # 2. Define the parent directory using a relative path ("..")
    parent_dir <- file.path(script_dir, "../..")
    
    # 3. Set the directory for knitr (Markdown chunks) and the current R session
    knitr::opts_knit$set(root.dir = parent_dir)
    setwd(parent_dir)
  }
  
  output_folder <- "./5-SensitivityTesting/output/"
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(sf)
  library(purrr)
  library(ggplot2)
  library(fitdistrplus)
  library(MASS)
  library(e1071)
  
}
getwd()


################STEP 1: Import animal datasets to use as "truth"###################
# Check if the file is there
# 

animalDataPath <- "./5-SensitivityTesting/data/animal_grid_Endemism.RData"
if (file.exists(animalDataPath)) {
  message("The animal data already exists. You can skip to STEP 2....")
} else {
  warning("Animal data not found. Continue running...")
}

# 1. Get the full paths so R knows exactly where to find them
file_paths <- list.files("./2-ModellingAnimals/data", 
                         pattern = "_cleaned.csv", 
                         recursive = TRUE, 
                         full.names = TRUE) # Crucial: gives full relative path

# 2. Use named paths so the .id column uses the actual file names
names(file_paths) <- basename(file_paths)

# 3. Read and combine into one master data frame
allAnimalDatasets <- map_df(file_paths, read.csv, .id = "source_file")

summary(allAnimalDatasets)
head(allAnimalDatasets)
length(allAnimalDatasets)

allAnimalDatasets<-allAnimalDatasets |> select(species,decimalLatitude,decimalLongitude,occurrenceStatus )
length(allAnimalDatasets)

#only present occurrences, no absences should be included
presentAnimals <- allAnimalDatasets %>% 
  filter(occurrenceStatus == "PRESENT")

#number of occurrences
dim(presentAnimals)

#number of species
length(
       union(unlist(presentAnimals["species"]),
       unlist(presentAnimals["species"]))
       )

################Plant import and setup###################

dataDir <- paste0("./1-FilteringPlantData/output")

occ_sf <- readRDS(paste0(dataDir, "/species_counts4_filtered.rds"))

# Extract the grid from the plant data to preserve Stage 0 grid_ids.
# Do NOT rebuild via st_make_grid() — different bounding box produces
# different grid_ids that break downstream joins with species_counts.
grid_size_m <- 22360
grid_sf <- occ_sf |>
  dplyr::group_by(grid_id) |>
  dplyr::slice(1) |>
  dplyr::ungroup() |>
  dplyr::select(grid_id, geometry)


amazon_shp <- st_read(
  paste0(dirname(getwd()), "/AmAni/0-CuratingPlantData/data/amazon_shapefile/amazon_sensulatissimo_gmm_v1.shp"),
  quiet = TRUE
) |> st_transform(st_crs(grid_sf))

################Convert animal data to spatial###################

animal_sf <- st_as_sf(
  presentAnimals,
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
summary(animal_grid)


# Define where you want to save it (adjust the folder path to match your directory)
save_path <- "./5-SensitivityTesting/data/animal_grid_data.RData"

# Save the object
save(animal_grid, file = save_path)

animal_endemism <- animal_grid |>
  group_by(grid_id) |>
  summarize(animal_richness = n_distinct(Morphotype), .groups = "drop")

# Save the object
 save_path_End <- "./5-SensitivityTesting/data/animal_grid_endemism.RData"

save(animal_endemism, file = save_path_End)



################STEP 2: Simulate completely sampled Amazon using randomized plant distribution data to create simulated truth###################

#Load animal data from above workflow

rm(animal_grid)
load("./5-SensitivityTesting/data/animal_grid_data.RData")
load("./5-SensitivityTesting/data/animal_grid_endemism.RData")


animal_grid_data <- grid_sf |>
  dplyr::select(grid_id, geometry) |>
  dplyr::left_join(
    animal_endemism |> sf::st_drop_geometry() |> dplyr::select(grid_id, animal_richness),
    by = "grid_id"
  ) |>
  mutate(animal_richness = replace_na(animal_richness, 0))

# Verify it loaded successfully by checking its summary again
summary(animal_grid)
summary(animal_grid_data)

#plot it to check manually later

timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
output_filename <- paste0("AnimalEmpir.RawEndem.Richn", timestamp, ".pdf")
output_path <- file.path(output_folder, output_filename)

# Create the PDF plot
pdf(output_path, width = 8, height = 6)
p<-ggplot() +
  geom_sf(data = animal_grid_data, aes(fill = animal_richness ), color = NA) +
  scale_fill_viridis_c(option = "plasma", na.value = "white") +
  theme_minimal() +
  labs(title = "Species-Level Richness (500 sq km grid)", fill = "Richness")
print(p)
dev.off()

###calculate raw endemism
rename_grid_id <- function(df) {
  if ("grid_id.x" %in% colnames(df)) {
    df <- dplyr::rename(df, grid_id = grid_id.x)
  }
  df
}

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

# Clip to Amazon boundary before plotting
animal_grid_data_clipped <- st_intersection(
  animal_endemism_grid %>% filter(animal_endemism_count > 0),
  amazon_shp
)

timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
output_filename <- paste0("AnimalEmpir.RawEndem.", timestamp, ".pdf")
output_path <- file.path(output_folder, output_filename)

pdf(output_path, width = 8, height = 6)

p<-ggplot(animal_grid_data_clipped) +
  geom_sf(aes(fill = animal_endemism_count), color = NA) +
  scale_fill_viridis_c(
    name = "Endemic\nspecies",
    option = "plasma"
  ) +
  theme_minimal() +
  labs(
    title = "Animal Endemism",
    x = NULL,
    y = NULL
  )
print(p)
dev.off()



###################STEP 3: Subsample complete true data to create testing datasets######################

  #Examine features of full dataset
  #Calculate species grid cell occupancy
species_ranges <- animal_grid %>%
  st_drop_geometry() %>%
  distinct(species, grid_id) %>%     # remove duplicate records within a grid
  count(species, name = "n_cells")

head(species_ranges)

ggplot(species_ranges, aes(x = n_cells)) +
  geom_histogram(bins = 50) +
  scale_x_log10() +
  theme_bw() +
  labs(
    x = "Occupied grid cells (log scale)",
    y = "Number of species"
  )



x <- species_ranges$n_cells

fit_gamma <- fitdist(x, "gamma")
fit_lnorm <- fitdist(x, "lnorm")
fit_exp <- fitdistr(x, "exponential")
lambda <- 1 / mean(x)

#lower is better 
AIC(fit_gamma)
AIC(fit_lnorm)
AIC(fit_exp)


data.frame(
  n_species = length(x),
  mean = mean(x),
  median = median(x),
  sd = sd(x),
  variance = var(x),
  min = min(x),
  max = max(x),
  skewness = skewness(x),
  kurtosis = kurtosis(x)
)


mu <- mean(log(species_ranges$n_cells))
sigma <- sd(log(species_ranges$n_cells))

mu
sigma


ggplot(species_ranges, aes(n_cells)) +
  geom_histogram(aes(y = after_stat(density)),
                 bins = 50,
                 fill = "grey80") +
  stat_function(
    fun = dlnorm,
    args = list(meanlog = mu,
                sdlog = sigma),
    linewidth = 1
  ) +
  scale_x_log10() +
  theme_bw()

quantile(species_ranges$n_cells,
         probs = c(0.025, 0.25, 0.5, 0.75, 0.975))

# a. Low to high endemism training datasets#######
#stratified (biased) resampling scheme over the species range-size distribution
#This will create true datasets with different endemism characteristics
quantile(species_ranges$n_cells,
         probs = c(0.05, 0.25, 0.5, 0.75, 0.95, 0.975, 0.99))

species_ranges <- species_ranges %>%
  mutate(
    quantile_class = ntile(n_cells, 5)  # 5 strata (quintiles)
  )


quantile_summary <- species_ranges %>%
  group_by(quantile_class) %>%
  summarise(
    n_species = n(),
    min_n_cells = min(n_cells),
    max_n_cells = max(n_cells),
    mean_n_cells = mean(n_cells),
    median_n_cells = median(n_cells),
    .groups = "drop"
  )

quantile_summary

#This function will allow me to get occurence data from the results of the biased sampling below.
rebuild_community <- function(animal_grid, sampled_df, scenario_name = NULL) {
  
  sp <- unique(sampled_df$species)
  
  out <- animal_grid %>%
    filter(species %in% sp) %>%
    
    # enforce 1 record per species per grid cell
    distinct(species, grid_id, .keep_all = TRUE)
  
  if (!is.null(scenario_name)) {
    out$scenario <- scenario_name
  }
  
  out
  
}
#the weights should add up to 1. They probably don't NEED to, but it makes sense to limit them to 1 so the weights across the different schemes are more easily compared

biased_species_sample <- function(data,
                                  n_quantiles = 5,
                                  weights,
                                  sample_frac = 0.8,
                                  seed = 1) {
  
  set.seed(seed)
  
  df <- data
  
  # 1. assign quantiles using ntile (robust to ties)
  df <- df %>%
    mutate(quantile_class = ntile(n_cells, n_quantiles))
  
  # 2. ensure weights are properly named
  if (is.numeric(weights)) {
    if (is.null(names(weights))) {
      names(weights) <- as.character(seq_along(weights))
    }
  }
  
  # 3. map weights to quantile classes safely
  df <- df %>%
    mutate(
      quantile_class_chr = as.character(quantile_class),
      weight = weights[quantile_class_chr]
    )
  
  # 4. safety check: missing weights
  if (any(is.na(df$weight))) {
    stop("Some quantile classes do not have matching weights.")
  }
  
  # 5. normalize probabilities
  df <- df %>%
    mutate(prob = weight / sum(weight))
  
  # 6. sample 80% (or user-defined fraction)
  n_target <- floor(sample_frac * nrow(df))
  
  sampled <- df %>%
    slice_sample(n = n_target, weight_by = prob, replace = FALSE)
  
  # 7. diagnostics: quantile composition shift
  quantile_shift <- rbind(
    full = prop.table(table(df$quantile_class)),
    sample = prop.table(table(sampled$quantile_class))
  )
  
  list(
    sampled_data = sampled,
    full_data = df,
    quantile_shift = quantile_shift,
    weights_used = weights
  )
}


#high endemism weights
baseWeight <- 0.2
extremeEndemicWeights <- baseWeight + 0.15
largeRangeWeight <- baseWeight - 0.099
extremeRangeWeight <- baseWeight - 0.199

weights <- c(
  "1" = extremeEndemicWeights,
  "2" = extremeEndemicWeights,
  "3" = baseWeight+0.1,
  "4" = largeRangeWeight,
  "5" = extremeRangeWeight
)
weights
sum(weights)

HighEndemismResult <- biased_species_sample(
  data = species_ranges,
  weights = weights,
  sample_frac = 0.8,
  seed = 1
)

#medium endemism weights
baseWeight <- 0.2
extremeEndemicWeights <- baseWeight + 0.075
largeRangeWeight <- baseWeight - 0.05
extremeRangeWeight <- baseWeight - 0.1

weights <- c(
  "1" = extremeEndemicWeights,
  "2" = extremeEndemicWeights,
  "3" = baseWeight,
  "4" = largeRangeWeight,
  "5" = extremeRangeWeight
)

mediumEndemismResult <- biased_species_sample(
  data = species_ranges,
  weights = weights,
  sample_frac = 0.8,
  seed = 1
)

#normal endemism weights
baseWeight <- 0.2
extremeEndemicWeights <- baseWeight + 0.1
largeRangeWeight <- baseWeight - 0.1999
extremeRangeWeight <- baseWeight - 0.1999

weights <- c(
  "1" = baseWeight,
  "2" = baseWeight,
  "3" = baseWeight,
  "4" = baseWeight,
  "5" = baseWeight
)
sum(weights)


normalEndemismResult <- biased_species_sample(
  data = species_ranges,
  weights = weights,
  sample_frac = 0.8,
  seed = 1
)


#low endemism weights
baseWeight <- 0.2
extremeEndemicWeights <- baseWeight - 0.075
largeRangeWeight <- baseWeight + 0.05
extremeRangeWeight <- baseWeight + 0.1

weights <- c(
  "1" = extremeEndemicWeights,
  "2" = extremeEndemicWeights,
  "3" = baseWeight,
  "4" = largeRangeWeight,
  "5" = extremeRangeWeight
)

sum(weights)


lowEndemismResult <- biased_species_sample(
  data = species_ranges,
  weights = weights,
  sample_frac = 0.8,
  seed = 1
)


#VERY low endemism weights
baseWeight <- 0.2
extremeEndemicWeights <- baseWeight - 0.1
largeRangeWeight <- baseWeight + 0.0999
extremeRangeWeight <- baseWeight + 0.1999

weights <- c(
  "1" = extremeEndemicWeights,
  "2" = extremeEndemicWeights,
  "3" = baseWeight-0.1,
  "4" = largeRangeWeight,
  "5" = extremeRangeWeight
)
weights
sum(weights)

veryLowEndemismResult <- biased_species_sample(
  data = species_ranges,
  weights = weights,
  sample_frac = 0.8,
  seed = 1
)



High <- HighEndemismResult$sampled_data %>% mutate(scenario = "high")
Medium <- mediumEndemismResult$sampled_data %>% mutate(scenario = "medium")
Normal <- normalEndemismResult$sampled_data %>% mutate(scenario = "normal")
Low <- lowEndemismResult$sampled_data %>% mutate(scenario = "low")
VeryLow <- veryLowEndemismResult$sampled_data %>% mutate(scenario = "very_low")

all_samples <- bind_rows(High, Medium, Normal, Low, VeryLow)

all_samples <- all_samples %>%
  mutate(scenario = factor(
    scenario,
    levels = c("high", "medium", "normal", "low", "very_low"),
    ordered = TRUE
  ))


scenario_summary1 <- all_samples %>%
  group_by(scenario) %>%
  summarise(
    
    #mean range size
    mean_range_size = mean(n_cells),
    
    #median range size
    median_range_size = median(n_cells),
    
    #std dev range size
    sd_range_size = sd(n_cells),
    
    #skewness of range size
    skew_range_size = skewness(n_cells),
    .groups = "drop"
  )

scenario_summary1

ggplot(all_samples, aes(x = log10(n_cells), fill = scenario)) +
  geom_density(alpha = 0.35) +
  theme_bw() +
  labs(
    x = "log10(range size in grid cells)",
    y = "Density"
  )


ggplot(all_samples, aes(x = scenario, y = n_cells)) +
  geom_boxplot() +
  scale_y_log10() +
  theme_bw() +
  labs(y = "Range size (grid cells, log scale)")


Low_comm <- rebuild_community(animal_grid, Low, "low")
High_comm <- rebuild_community(animal_grid, High, "high")
Medium_comm <- rebuild_community(animal_grid, Medium, "medium")
Normal_comm <- rebuild_community(animal_grid, Normal, "normal")
VeryLow_comm <- rebuild_community(animal_grid, VeryLow, "very_low")

all_comm <- bind_rows(
  High_comm, 
  Medium_comm,
  Normal_comm,
  Low_comm,
  VeryLow_comm
)

all_comm <- all_comm %>%
  mutate(scenario = factor(
    scenario,
    levels = c("high", "medium", "normal", "low", "very_low"),
    ordered = TRUE
  ))

#because the above are "true" datasets (i.e., each one has a different true number of total species in the simulation), I need to subsample these again randomly to simulate sampling. NOTE: That I sampled 80 % of species in the original sampling above, but that was how we modified the endemism, by sampling species with different endemism characteristics and omitting ones without the characteristics we wanted.
samplingfrac<-0.8


dataSet_list <- list(
  high = High_comm %>%
    slice_sample(n = floor(nrow(High_comm) * samplingfrac)),
  
  med = Medium_comm %>%
    slice_sample(n = floor(nrow(Medium_comm) * samplingfrac)),
  
  normal = Normal_comm %>%
    slice_sample(n = floor(nrow(Normal_comm) * samplingfrac)),
  
  low = Low_comm %>%
    slice_sample(n = floor(nrow(Low_comm) * samplingfrac)),
  
  veryLow = VeryLow_comm %>%
    slice_sample(n = floor(nrow(VeryLow_comm) * samplingfrac))
)

dataSet_list_DF <- bind_rows(dataSet_list, .id = "scenario")

# sapply(list(
#   high = High_comm,
#   medium = Medium_comm,
#   normal = Normal_comm,
#   low = Low_comm,
#   very_low = VeryLow_comm
# ), function(x) length(unique(x$species)))
# 
# # sapply(list(
# #   high = High_comm,
# #   medium = Medium_comm,
# #   normal = Normal_comm,
# #   low = Low_comm,
# #   very_low = VeryLow_comm
# # ), function(x) length(x$species))


scenario_summary2 <- all_comm %>%
  group_by(scenario) %>%
  summarise(
    
    # richness structure
    n_species = n_distinct(species),
    
    # spatial footprint
    n_grid_cells = n_distinct(grid_id),
    
    # total observations
    total_observations = n(),
    
    # per-species intensity
    mean_obs_per_species = n() / n_distinct(species),
    
    # A fixed metric in this simulation###NOTE: CHANGE IF USING IN ANOTHER SCENARIO WHERE YOU MODIFY THIS
    obs_per_speciesInCell = 1,
    
    # per-cell packing
    mean_obs_per_cell = n() / n_distinct(grid_id),
    
    # number of species divided by number of grid cells
    mean_rangeSize = n_distinct(species) / n_distinct(grid_id),
    
    .groups = "drop"
  )


scenario_summary3 <- dataSet_list_DF %>%
  group_by(scenario) %>%
  summarise(
    
    # richness structure
    n_species_sampled = n_distinct(species),
    
    # spatial footprint
    n_grid_cells_sampled = n_distinct(grid_id),
    
    # total observations
    total_observations_sampled = n(),
    
    # per-species intensity
    mean_obs_per_species_sampled = n() / n_distinct(species),
    
    # A fixed metric in this simulation###NOTE: CHANGE IF USING IN ANOTHER SCENARIO WHERE YOU MODIFY THIS
    obs_per_speciesInCell_sampled = 1,
    
    # per-cell packing
    mean_obs_per_cell_sampled = n() / n_distinct(grid_id),
    
    # number of species divided by number of grid cells
    mean_rangeSize_sampled = n_distinct(species) / n_distinct(grid_id),
    
    .groups = "drop"
  )


scenario_summary2_clean <- scenario_summary2 %>%
  st_drop_geometry()
scenario_summary3_clean <- scenario_summary3 %>%
  st_drop_geometry()


scenario_summary_full0 <- scenario_summary1 %>%
  left_join(scenario_summary2_clean, by = "scenario")

scenario_summary_full <- scenario_summary_full0 %>%
  left_join(scenario_summary3_clean, by = "scenario")



timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
outfile <- paste0(output_folder,paste0("endemismSamples.summary_", timestamp, ".csv"))
write.csv(scenario_summary_full, outfile, row.names = FALSE)

scenario_list <- list(
  high = HighEndemismResult$sampled_data,
  medium = mediumEndemismResult$sampled_data,
  normal = normalEndemismResult$sampled_data,
  low = lowEndemismResult$sampled_data,
  very_low = veryLowEndemismResult$sampled_data
)

for (nm in names(scenario_list)) {
  
  df <- scenario_list[[nm]]
  
  outfile <- paste0(output_folder,paste0(nm, "Endemism.80pec.Info_", timestamp, ".csv"))
  
  write.csv(df, outfile, row.names = FALSE)
}



for (nm in names(dataSet_list)) {
  
  df <- scenario_list[[nm]]
  
  outfile <- paste0(output_folder,paste0(nm, ".Endemism.sample_", timestamp, ".csv"))
  
  write.csv(df, outfile, row.names = FALSE)
}


# b. Low to high total sampling (% of species) training datasets#######
#The above scenarios modified endemism of the constituent species and created true datasets. in this one I am modifying the % of species sampled overall. This will create a sampled datasets

scenario_levels <- c("min", "quarter", "half", "high", "higher", "highest")

weights <- c(
  "1" = baseWeight,
  "2" = baseWeight,
  "3" = baseWeight,
  "4" = baseWeight,
  "5" = baseWeight
)

samplingFracs <- c(
 0.1,0.25, 0.5, 0.75, 0.85, 0.95
)

minSamplingData <- biased_species_sample(
  data = species_ranges,
  weights = weights,
  sample_frac = samplingFracs[[1]],
  seed = 1
)

quarterSamplingData <- biased_species_sample(
  data = species_ranges,
  weights = weights,
  sample_frac = samplingFracs[[2]],
  seed = 1
)

halfSamplingData <- biased_species_sample(
  data = species_ranges,
  weights = weights,
  sample_frac = samplingFracs[[3]],
  seed = 1
)

highSamplingData <- biased_species_sample(
  data = species_ranges,
  weights = weights,
  sample_frac = samplingFracs[[4]],
  seed = 1
)

higherSamplingData <- biased_species_sample(
  data = species_ranges,
  weights = weights,
  sample_frac = samplingFracs[[5]],
  seed = 1
)

highestSamplingData <- biased_species_sample(
  data = species_ranges,
  weights = weights,
  sample_frac = samplingFracs[[6]],
  seed = 1
)

####Combine results

minSampling <- minSamplingData$sampled_data %>% mutate(scenario = "min")
quarterSampling <- quarterSamplingData$sampled_data %>% mutate(scenario = "quarter")
halfSampling <- halfSamplingData$sampled_data %>% mutate(scenario = "half")
highSampling <- highSamplingData$sampled_data %>% mutate(scenario = "high")
higherSampling <- higherSamplingData$sampled_data %>% mutate(scenario = "higher")
highestSampling <- highestSamplingData$sampled_data %>% mutate(scenario = "highest")

all_samples <- bind_rows(minSampling, quarterSampling, halfSampling, highSampling, higherSampling, highestSampling)

all_samples <- all_samples %>%
  mutate(scenario = factor(
    scenario,
    levels = scenario_levels,
    ordered = TRUE
  ))


scenario_summary1 <- all_samples %>%
  group_by(scenario) %>%
  summarise(
    
    #mean range size
    mean_range_size = mean(n_cells),
    
    #median range size
    median_range_size = median(n_cells),
    
    #std dev range size
    sd_range_size = sd(n_cells),
    
    #skewness of range size
    skew_range_size = skewness(n_cells),
    .groups = "drop"
  )

scenario_summary1

ggplot(all_samples, aes(x = log10(n_cells), fill = scenario)) +
  geom_density(alpha = 0.35) +
  theme_bw() +
  labs(
    x = "log10(range size in grid cells)",
    y = "Density"
  )


ggplot(all_samples, aes(x = scenario, y = n_cells)) +
  geom_boxplot() +
  scale_y_log10() +
  theme_bw() +
  labs(y = "Range size (grid cells, log scale)")



minSampling_comm <- rebuild_community(animal_grid, minSampling,"min")
quarterSampling_comm <- rebuild_community(animal_grid, quarterSampling, "quarter")
halfSampling_comm <-rebuild_community(animal_grid,  halfSampling, "half")
highSampling_comm <- rebuild_community(animal_grid, highSampling, "high")
higherSampling_comm <- rebuild_community(animal_grid, higherSampling, "high")
highestSamplingData_comm <- rebuild_community(animal_grid, highestSampling, "highest")

all_comm <- bind_rows(
  minSampling_comm, 
  quarterSampling_comm,
  halfSampling_comm,
  highSampling_comm,
  higherSampling_comm,
  highestSamplingData_comm
)

all_comm <- all_comm %>%
  mutate(scenario = factor(
    scenario,
    levels = c("min", "quarter", "half", "high","higher", "highest"),
    ordered = TRUE
  ))


scenario_summary2 <- all_comm %>%
  group_by(scenario) %>%
  summarise(
    
    # richness structure
    n_species = n_distinct(species),
    
    # spatial footprint
    n_grid_cells = n_distinct(grid_id),
    
    # total observations
    total_observations = n(),
    
    # per-species intensity
    mean_obs_per_species = n() / n_distinct(species),
    
    # A fixed metric in this simulation###NOTE: CHANGE IF USING IN ANOTHER SCENARIO WHERE YOU MODIFY THIS
    obs_per_speciesInCell = 1,
    
    # per-cell packing
    mean_obs_per_cell = n() / n_distinct(grid_id),
    
    # number of species divided by number of grid cells
    mean_rangeSize = n_distinct(species) / n_distinct(grid_id),
    
    .groups = "drop"
  )

scenario_summary2_clean <- scenario_summary2 %>%
  st_drop_geometry()
scenario_summary_full <- scenario_summary1 %>%
  left_join(scenario_summary2_clean, by = "scenario")

scenario_summary_full <- scenario_summary_full %>%
  mutate(
    scenario = factor(scenario, levels = scenario_levels, ordered = TRUE)
  ) %>%
  arrange(scenario) %>%
  mutate(
    sampling_fraction = samplingFracs[as.integer(scenario)]
  )


timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
outfile <- paste0(output_folder,paste0("lowtoHighTotalSampling.summary_", timestamp, ".csv"))
write.csv(scenario_summary_full, outfile, row.names = FALSE)

scenario_list <- list(
  min_0.1 = minSampling, 
  med_0.25 = quarterSampling,
  half_0.5 = halfSampling,
  high_0.75 = highSampling,
  higher_0.85 = higherSampling,
  highest_0.95 = highestSampling
)

for (nm in names(scenario_list)) {
  
  df <- scenario_list[[nm]]
  
  outfile <- paste0(output_folder,paste0(nm, ".Sampling.Info_", timestamp, ".csv"))
  
  write.csv(df, outfile, row.names = FALSE)
}


dataSet_list <- list(
  min_0.1 = minSampling_comm, 
  med_0.25 = quarterSampling_comm,
  half_0.5 = halfSampling_comm,
  high_0.75 = highSampling_comm,
  higher_0.85 = higherSampling_comm,
  highest_0.95 = highestSamplingData_comm
)

for (nm in names(dataSet_list)) {
  
  df <- scenario_list[[nm]]
  
  outfile <- paste0(output_folder,paste0(nm, ".Sampling.sample_", timestamp, ".csv"))
  
  write.csv(df, outfile, row.names = FALSE)
}

# c. Low to high site sampling (N sites) of training datasets#######
#In this sampling regime, I am going to thin the data by site, instead of by taxon. This will create sampled datasets.


allGridIDs<-unique(animal_grid$grid_id)
samplingFrac <- 0.80

#Sample the grid IDs based on their sampling frequencing in the original dataset
grid_freq <- animal_grid %>%
  st_drop_geometry() %>%
  count(grid_id) %>%
  pull(n)

head(grid_freq)


filtered_animal_gridThreeA <- animal_grid %>%
  filter(grid_id %in% sample(allGridIDs,
                             3,
                             replace = FALSE,
                             prob = grid_freq)) %>%
  slice_sample(prop = samplingFrac)

filtered_animal_gridThreeB <- animal_grid %>%
  filter(grid_id %in% sample(allGridIDs,
                             3,
                             replace = FALSE,
                             prob = grid_freq)) %>%
  slice_sample(prop = samplingFrac)

filtered_animal_gridThreeC <- animal_grid %>%
  filter(grid_id %in% sample(allGridIDs,
                             3,
                             replace = FALSE,
                             prob = grid_freq)) %>%
  slice_sample(prop = samplingFrac)

filtered_animal_gridThreeD <- animal_grid %>%
  filter(grid_id %in% sample(allGridIDs,
                             3,
                             replace = FALSE,
                             prob = grid_freq)) %>%
  slice_sample(prop = samplingFrac)


filtered_animal_grid10A <- animal_grid %>%
  filter(grid_id %in% sample(allGridIDs,
                             10,
                             replace = FALSE,
                             prob = grid_freq)) %>%
  slice_sample(prop = samplingFrac)

filtered_animal_grid10B <- animal_grid %>%
  filter(grid_id %in% sample(allGridIDs,
                             10,
                             replace = FALSE,
                             prob = grid_freq)) %>%
  slice_sample(prop = samplingFrac)

filtered_animal_grid10C <- animal_grid %>%
  filter(grid_id %in% sample(allGridIDs,
                             10,
                             replace = FALSE,
                             prob = grid_freq)) %>%
  slice_sample(prop = samplingFrac)

filtered_animal_gridlow <- animal_grid %>%
  filter(grid_id %in% sample(allGridIDs,
                             floor(length(allGridIDs)*0.1),
                             replace = FALSE,
                             prob = grid_freq)) %>%
  slice_sample(prop = samplingFrac)

filtered_animal_gridMedium <- animal_grid %>%
  filter(grid_id %in% sample(allGridIDs,
                             floor(length(allGridIDs)*0.25),
                             replace = FALSE,
                             prob = grid_freq)) %>%
  slice_sample(prop = samplingFrac)

filtered_animal_gridHigh <- animal_grid %>%
  filter(grid_id %in% sample(allGridIDs,
                             floor(length(allGridIDs)*0.5),
                             replace = FALSE,
                             prob = grid_freq)) %>%
  slice_sample(prop = samplingFrac)

filtered_animal_gridHigher <- animal_grid %>%
  filter(grid_id %in% sample(allGridIDs,
                             floor(length(allGridIDs)*0.75),
                             replace = FALSE,
                             prob = grid_freq)) %>%
  slice_sample(prop = samplingFrac)

scenario_levels <- c("three", "ten", "low", "medium", "high", "higher")



dataSet_list <- list(
  threeA = filtered_animal_gridThreeA,
  threeB = filtered_animal_gridThreeB,
  threeC = filtered_animal_gridThreeC,
  threeD = filtered_animal_gridThreeD,
  
  tenA = filtered_animal_grid10A,
  tenB = filtered_animal_grid10B,
  tenC = filtered_animal_grid10C,
  
  low = filtered_animal_gridlow,
  medium = filtered_animal_gridMedium,
  high = filtered_animal_gridHigh,
  higher = filtered_animal_gridHigher
)


summary_table <- lapply(names(dataSet_list), function(dataset_name) {
  
  df <- dataSet_list[[dataset_name]] %>%
    st_drop_geometry()
  
  data.frame(
    dataset = dataset_name,
    
    # number of unique grid cells
    n_grid_cells = n_distinct(df$grid_id),
    
    # proport of unique grid cells sampled
    n_grid_cells = n_distinct(df$grid_id)/length(allGridIDs),
    
    # total observations (rows)
    total_observations = nrow(df),
    
    # average observations per species
    mean_obs_per_species =
      nrow(df) / n_distinct(df$species),
    
    # average observations per species-grid combination
    obs_per_speciesInCell =
      nrow(df) / n_distinct(df[c("species", "grid_id")]),
    
    # average observations per grid cell
    mean_obs_per_cell =
      nrow(df) / n_distinct(df$grid_id)
  )
  
}) %>%
  bind_rows()

summary_table

#####STILL NEED TO EXPORT SUMMARY TABLE AND DATASETS FROM THIS SECTION



####d. Clustering of samples in training dataset######
#When collecting actual data, the sites will probably be clustered. The result, is that the model is now extrapolating rather than interpolating. This could be less accurate.






















