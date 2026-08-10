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


animalDataPath<-"./5-SensitivityTesting/data/allGBIFCombined.19.6.2026.tsv"
if (file.exists(animalDataPath)) {
  message("The animal data already exists. You can skip to STEP 2....")
} else {
  warning("Animal data not found. Continue running...")
  source("./5-SensitivityTesting/scripts/0-CompileCleanedData.R")
}

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
###This can be skipped if the output files already exist
animalDataPath <- "./5-SensitivityTesting/data/simFromPlant/simulatedOccurencesForAnimal.rds"

presentAnimals<-readRDS(animalDataPath)
 
head(presentAnimals)

animal_sf <- presentAnimals
animal_grid<-animal_sf

# Map "scientificName" to "Morphotype" so downstream scripts work unchanged
animal_grid$Morphotype <- animal_grid$scientificName 

message("  Records within Amazon grid: ", nrow(animal_grid),
        " (", n_distinct(animal_grid$Morphotype), " species)")
summary(animal_grid)


save_path <- "./5-SensitivityTesting/data/animal_grid_data.RData"

# Save the object
save(animal_grid, file = save_path)

animal_endemism <- animal_grid |>
  group_by(grid_id) |>
  summarize(animal_richness = n_distinct(Morphotype), .groups = "drop")

# Save the object
 save_path_End <- "./5-SensitivityTesting/data/animal_grid_endemism.RData"

save(animal_endemism, file = save_path_End)


################TRUTH 1: BASELINE: Completely sampled Amazon using randomized plant distribution data to create simulated truth###################

#Load animal data from above workflow

rm(animal_grid)
load("./5-SensitivityTesting/data/animal_grid_data.RData")
load("./5-SensitivityTesting/data/animal_grid_endemism.RData")



#calculate Richness of grid IDs
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
output_filename <- paste0("An.Emp.Pla.Extr.Richn", timestamp, ".pdf")
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
output_filename <- paste0("An.Emp.Pla.Extr.End.", timestamp, ".pdf")
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

head(animal_grid_data_clipped)
true_animal_gridData<-animal_grid_data_clipped
true_animal_endemism<-animal_endemism
true_animal_occ_data<-animal_sf

################TRUTH 2: ENDEMISM RESAMPLES from complete true data to create true organisms with varying endemism######################

#We have two datasets so far:
#Sample - this is a combination of all the raw occurrence data from GBIF. It contains a variety of organisms and endemism types (high to low). We're calling this a "sample" because it is a sampling of the true-in-life unknown biodiversity. As such, it is our basis for parameterizing sample characteristics 
#True - this are the simulated occurence data based on the raw data (sample). We are going to consider this to be "true-in-simulation". It can serve as a base-line truth dataset and it can also serve as the dataset we draw species from to create other datasets.

load("./5-SensitivityTesting/data/rawData/animal_grid_data.RData")
load("./5-SensitivityTesting/data/rawData/animal_grid_endemism.RData")

  #Examine features of sampled dataset
  #Calculate species grid cell occupancy
head(animal_grid)
head(true_animal_occ_data)
true_animal_occ_data$species<-true_animal_occ_data$scientificName

true_species_ranges <- true_animal_occ_data %>%
  st_drop_geometry() %>%
  distinct(species  , grid_id) %>%     # remove duplicate records within a grid
  count(species , name = "n_cells")


sampled_species_ranges <- animal_grid %>%
  st_drop_geometry() %>%
  distinct(species  , grid_id) %>%     # remove duplicate records within a grid
  count(species , name = "n_cells")

head(sampled_species_ranges)
head(true_species_ranges)


ggplot(sampled_species_ranges, aes(x = n_cells)) +
  geom_histogram(bins = 50) +
  scale_x_log10() +
  theme_bw() +
  labs(
    x = "Occupied grid cells (log scale)",
    y = "Number of species"
  )

ggplot(true_species_ranges, aes(x = n_cells)) +
  geom_histogram(bins = 50) +
  scale_x_log10() +
  theme_bw() +
  labs(
    x = "Occupied grid cells (log scale)",
    y = "Number of species"
  )

# #Get a best fit function for this distribution ###Outdated
# x <- sampled_species_ranges$n_cells
# 
# fit_gamma <- fitdist(x, "gamma")
# fit_lnorm <- fitdist(x, "lnorm")
# fit_exp <- fitdistr(x, "exponential")
# lambda <- 1 / mean(x)
# 
# #lower is better 
# AIC(fit_gamma)
# AIC(fit_lnorm)
# AIC(fit_exp)
# 
# 
# data.frame(
#   n_species = length(x),
#   mean = mean(x),
#   median = median(x),
#   sd = sd(x),
#   variance = var(x),
#   min = min(x),
#   max = max(x),
#   skewness = skewness(x),
#   kurtosis = kurtosis(x)
# )
# 
# 
# mu <- mean(log(sampled_species_ranges$n_cells))
# sigma <- sd(log(sampled_species_ranges$n_cells))
# 
# mu
# sigma
# 
# 
# ggplot(sampled_species_ranges, aes(n_cells)) +
#   geom_histogram(aes(y = after_stat(density)),
#                  bins = 50,
#                  fill = "grey80") +
#   stat_function(
#     fun = dlnorm,
#     args = list(meanlog = mu,
#                 sdlog = sigma),
#     linewidth = 1
#   ) +
#   scale_x_log10() +
#   theme_bw()

quantile(sampled_species_ranges$n_cells,
         probs = c(0.025, 0.25, 0.5, 0.75, 0.975))

quantile(true_species_ranges$n_cells,
         probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
##From running the above code we can say that the sampled datasets have smaller ranges, which makes sense since the point of simulating the true ones based on the plant model, is that it would increase their ranges


################# a. Low to high endemism truth datasets#######

#stratified (biased) resampling scheme over the species range-size distribution
#This will create true datasets with different endemism characteristics
quantile(sampled_species_ranges$n_cells,
         probs = c(0.05, 0.25, 0.5, 0.75, 0.95, 0.975, 0.99))

sampled_species_ranges <- sampled_species_ranges %>%
  mutate(
    quantile_class = ntile(n_cells, 5)  # 5 strata (quintiles)
  )

true_species_ranges <- true_species_ranges %>%
  mutate(
    quantile_class = ntile(n_cells, 5)  # 5 strata (quintiles)
  )


sampled_quantile_summary <- sampled_species_ranges %>%
  group_by(quantile_class) %>%
  summarise(
    n_species = n(),
    min_n_cells = min(n_cells),
    max_n_cells = max(n_cells),
    mean_n_cells = mean(n_cells),
    median_n_cells = median(n_cells),
    .groups = "drop"
  )

true_quantile_summary <- true_species_ranges %>%
  group_by(quantile_class) %>%
  summarise(
    n_species = n(),
    min_n_cells = min(n_cells),
    max_n_cells = max(n_cells),
    mean_n_cells = mean(n_cells),
    median_n_cells = median(n_cells),
    .groups = "drop"
  )

#the difference between the ranges of the two datasets can inform the base weights I should use. To get the baseline subsample, I want to sample true_species_ranges by weighting these different quantile classes according to the quantile classes in sampled_quantile_summary. 
true_quantile_summary
sampled_quantile_summary

#in order to do this, we need to use the quantile classes of the sampled data to make the weights for the true data
  
# 1. Define the target proportions based on the sampled_quantile_summary.
# Since sampled classes 1 & 2 are both exactly 1 cell, that bin represents 40%.
# The others represent 20% each.
target_proportions <- tibble(
  sampled_bin = c("class1and2", "class3", "class4", "class5"),
  target_prop = c(0.40, 0.20, 0.20, 0.20)
)

# 2. Categorize the TRUE dataset using the SAMPLED boundaries
true_species_ranges <- true_species_ranges %>%
  mutate(
    sampled_bin = case_when(
      n_cells == 1 ~ "class1and2",
      n_cells > 1 & n_cells <= 3 ~ "class3",
      n_cells > 3 & n_cells <= 8 ~ "class4",
      n_cells > 8 ~ "class5"
    )
  )

# 3. Calculate the actual proportions of these bins in the TRUE dataset
true_proportions <- true_species_ranges %>%
  count(sampled_bin) %>%
  mutate(actual_prop = n / sum(n))

# 4. Calculate the required sampling weights
weights_df <- true_proportions %>%
  left_join(target_proportions, by = "sampled_bin") %>%
  mutate(sampling_weight = target_prop / actual_prop)

# View your weights
sum(weights_df$sampling_weight)
weights_df$sampling_weight<-weights_df$sampling_weight/sum(weights_df$sampling_weight)

head(animal_grid)

head(sampled_species_ranges)

sampled_quantile_summary
true_quantile_summary


weightKernel<-c(
  0.2,#high endemism weights
  0, #mid endemism weights
  - 0.15, #low endemism weights
  - 0.25 #very low endemism weights
)
sum(weightKernel)+weightKernel[[1]]<0
weights<-make_weights(weightKernel)
sum(weights)

HighEndemismResult <- biased_species_sample(
  data = true_species_ranges,
  weights = weights,
  sample_frac = 0.8,
  seed = 1
)
str(HighEndemismResult)

#medium endemism weights
weightKernel<-c(
  0.1,#high endemism weights
  0, #mid endemism weights
  - 0.05, #low endemism weights
  - 0.15 #very low endemism weights
)
sum(weightKernel)+weightKernel[[1]]
sum(weightKernel)+weightKernel[[1]]<0.001 ##should be TRUE
weights<-make_weights(weightKernel)
sum(weights)==1 ##should be TRUE

mediumEndemismResult <- biased_species_sample(
  data = true_species_ranges,
  weights = weights,
  sample_frac = 0.8,
  seed = 1
)

#sample endemism weights
weights<-make_weights()
sum(weights)==1 ##should be TRUE

sampleEndemismResult <- biased_species_sample(
  data = true_species_ranges,
  weights = weights,
  sample_frac = 0.8,
  seed = 1
)

#True endemism weights
trueEndemismResult <- biased_species_sample(
  data = true_species_ranges,
  weights = unlist(c(weights_df$actual_prop[[1]], weights_df$actual_prop)),
  sample_frac = 0.8,
  seed = 1
)


#low endemism weights
weightKernel<-c(
  - 0.075,#high endemism weights
  0, #mid endemism weights
  +0.05, #low endemism weights
  + 0.1 #very low endemism weights
)
sum(weightKernel)+weightKernel[[1]]
sum(weightKernel)+weightKernel[[1]]<0.001 ##should be TRUE
weights<-make_weights(weightKernel)

lowEndemismResult <- biased_species_sample(
  data = true_species_ranges,
  weights = weights,
  sample_frac = 0.8,
  seed = 1
)

#VERY low endemism weights
weightKernel<-c(
  - 0.2,#high endemism weights
  -0.3, #mid endemism weights
  0.3, #low endemism weights
  0.4 #very low endemism weights
)
sum(weightKernel)+weightKernel[[1]]
sum(weightKernel)+weightKernel[[1]]<0.001 ##should be TRUE
weights<-make_weights(weightKernel)


veryLowEndemismResult <- biased_species_sample(
  data = true_species_ranges,
  weights = weights,
  sample_frac = 0.8,
  seed = 1
)

High <- HighEndemismResult$sampled_data %>% mutate(scenario = "high")
Medium <- mediumEndemismResult$sampled_data %>% mutate(scenario = "medium")
Sample_w <- sampleEndemismResult$sampled_data %>% mutate(scenario = "sample")
True_w <- trueEndemismResult$sampled_data %>% mutate(scenario = "true")
Low <- lowEndemismResult$sampled_data %>% mutate(scenario = "low")
VeryLow <- veryLowEndemismResult$sampled_data %>% mutate(scenario = "very_low")

all_samples <- bind_rows(High, Medium, Low, VeryLow, Sample_w, True_w)

all_samples <- all_samples %>%
  mutate(scenario = factor(
    scenario,
    levels = c("high", "medium", "low", "very_low", "sample", "true"),
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


#this will recreate occurrence data based on the sampling scenarios
VeryLow_comm <- rebuild_community(true_animal_occ_data, VeryLow, "very_low")
Low_comm <- rebuild_community(true_animal_occ_data, Low, "low")
Medium_comm <- rebuild_community(true_animal_occ_data, Medium, "medium")
High_comm <- rebuild_community(true_animal_occ_data, High, "high")
sample_comm<-rebuild_community(true_animal_occ_data, Sample_w, "sample")
true_comm<-rebuild_community(true_animal_occ_data, True_w, "true")

#we're going to make a summary table of all of these organisms, so we have to combine them into a single list
all_comm <- list(
  High_comm, 
  Medium_comm,
  Low_comm,
  VeryLow_comm,
  sample_comm,
  true_comm
)

summary(all_comm)

#Combine the list of sf objects into a single sf data frame
all_comm_flat <- bind_rows(all_comm)

#this is a summary of all the datasets we can combine with the stats from the previous datasets

scenario_summary2 <- all_comm_flat %>%
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


scenario_summary_full0 <- scenario_summary1 %>%
  left_join(scenario_summary2_clean, by = "scenario")

#output the summary info for these datasets
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
outfile <- paste0(output_folder,paste0("endemismSamples.summary_", timestamp, ".csv"))
write.csv(scenario_summary_full0, outfile, row.names = FALSE)

#object with the endemism class weights for each scenario
scenario_list <- list(
  high = HighEndemismResult$weights_used    ,
  medium = mediumEndemismResult$weights_used  ,
  normal = normalEndemismResult$weights_used  ,
  low = lowEndemismResult$weights_used  ,
  very_low = veryLowEndemismResult$weights_used  
)

#the occurrence datasets
dataSet_list <- list(
  high = HighEndemismResult$sampled_data,
  medium = mediumEndemismResult$sampled_data,
  normal = normalEndemismResult$sampled_data,
  low = lowEndemismResult$sampled_data,
  very_low = veryLowEndemismResult$sampled_data
)

scenario_df <- bind_rows(scenario_list, .id = "scenario")

#export the weights
  {outfile <- paste0(output_folder,paste0("EndemismWeights.80pec._", timestamp, ".csv"))

  write.csv(scenario_df, outfile, row.names = FALSE)}


#export each of the datasets
for (nm in names(dataSet_list)) {
  
  df <- scenario_list[[nm]]
  
  outfile <- paste0(output_folder,paste0(nm, ".Endemism.sample_", timestamp, ".csv"))
  
  write.csv(df, outfile, row.names = FALSE)
}


################TRUTH 3: RICHNESS RESAMPLES from complete true data to create true organisms with varying richness######################

#From above, we have the following addition truth scenarios
#scenario_list$high - a high endemism organism
#scenario_list$medium - a medium endemism organism
#scenario_list$normal - an organism with typical endemism
#scenario_list$low - a low endemism organism
#scenario_list$very_low - a very low endemism organism

#Now we're going to make true-in-simulation organisms with the same proportional endemism, but different numbers of species.
#This is also going to be different from the above simulations because I am going to sample the original data with replacement.

#################b. Low to high richness truth datasets#######




HighestRichnessResult <- species_richness_sample(data = true_species_ranges,
                                      richness = 1000000,
                                      seed = 1)
HighRichnessResult <- species_richness_sample(data = true_species_ranges,
                                              richness = 100000,
                                              seed = 1)
MediumRichnessResult <- species_richness_sample(data = true_species_ranges,
                                              richness = 10000,
                                              seed = 1)
LowRichnessResult <- species_richness_sample(data = true_species_ranges,
                                              richness = 1000,
                                              seed = 1)
VerylowRichnessResult <- species_richness_sample(data = true_species_ranges,
                                              richness = 100,
                                              seed = 1)


highest_R <- HighestRichnessResult$sampled_data %>% mutate(scenario = "highest")
high_R <- HighRichnessResult$sampled_data %>% mutate(scenario = "high")
medium_R <- MediumRichnessResult$sampled_data %>% mutate(scenario = "medium")
low_R <- LowRichnessResult$sampled_data %>% mutate(scenario = "low")
veryLow_R <- VerylowRichnessResult$sampled_data %>% mutate(scenario = "very_low")

all_samples_R <- bind_rows(highest_R, high_R, medium_R, low_R, veryLow_R)

all_samples_R <- all_samples_R %>%
  mutate(scenario = factor(
    scenario,
    levels = c("highest", "high", "medium", "low", "very_low"),
    ordered = TRUE
  ))


scenario_summary_Richness <- all_samples_R %>%
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

scenario_summary_Richness

ggplot(all_samples_R, aes(x = log10(n_cells), fill = scenario)) +
  geom_density(alpha = 0.35) +
  theme_bw() +
  labs(
    x = "log10(range size in grid cells)",
    y = "Density"
  )


ggplot(all_samples_R, aes(x = scenario, y = n_cells)) +
  geom_boxplot() +
  scale_y_log10() +
  theme_bw() +
  labs(y = "Range size (grid cells, log scale)")


highestR_comm <- rebuild_community(true_animal_occ_data, highest_R, "highest")
highR_comm <- rebuild_community(true_animal_occ_data, high_R, "high")
mediumR_comm <- rebuild_community(true_animal_occ_data, medium_R, "medium")
lowR_comm <- rebuild_community(true_animal_occ_data, low_R, "low")
veryLowR_comm<-rebuild_community(true_animal_occ_data, veryLow_R, "veryLow")

all_comm_R <- bind_rows(
  highestR_comm, 
  highR_comm,
  mediumR_comm,
  lowR_comm,
  veryLowR_comm,
)

all_comm_R <- all_comm_R %>%
  mutate(scenario = factor(
    scenario,
    levels = c("highest", "high", "medium", "low", "very_low"),
    ordered = TRUE
  ))

scenario_summary_Richness <- all_comm_R %>%
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


timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
outfile <- paste0(output_folder,paste0("richnessSamples.summary_", timestamp, ".csv"))
write.csv(scenario_summary_Richness, outfile, row.names = FALSE)

#####NOTE: We still need to add something to export the richness datasets.




#object with the endemism class weights for each scenario
 dataSet_list_R<- list(
  highest =  HighestRichnessResult$sampled_data,
  high = HighRichnessResult$sampled_data,
  medium = MediumRichnessResult$sampled_data,
  low = LowRichnessResult$sampled_data,
  very_low = VerylowRichnessResult$sampled_data
)

#the occurrence datasets
 scenario_list_R <- list(
  highest =  HighestRichnessResult$sampled_data%>%count(),
  high = HighRichnessResult$sampled_data%>%count(),
  medium = MediumRichnessResult$sampled_data%>%count(),
  low = LowRichnessResult$sampled_data%>%count(),
  very_low = VerylowRichnessResult$sampled_data%>%count()
)

scenario_df_R <- bind_rows(dataSet_list_R, .id = "scenario")

#export the weights
{outfile <- paste0(output_folder,paste0("richnessScenarios.", timestamp, ".csv"))
  
  write.csv(scenario_df_R, outfile, row.names = FALSE)}


#export each of the datasets
for (nm in names(dataSet_list)) {
  
  df <- scenario_list[[nm]]
  
  outfile <- paste0(output_folder,paste0(nm, ".Endemism.sample_", timestamp, ".csv"))
  
  write.csv(df, outfile, row.names = FALSE)
}

#CONTINUE HERE
# c. IN PROGRESS: Low to high total sampling (% of species) training datasets#######
#The above scenarios modified endemism of the constituent species and created true datasets. in this one I am modifying the % of species sampled overall. T

scenario_levels <- c("min", "quarter", "half", "high", "higher", "highest")

weights <- make_weights()

samplingFracs <- c(
 0.1,0.25, 0.5, 0.75, 0.85, 0.95
)

minSamplingData <- biased_species_sample(
  data = true_species_ranges,
  weights = weights,
  sample_frac = samplingFracs[[1]],
  seed = 1
)

quarterSamplingData <- biased_species_sample(
  data = true_species_ranges,
  weights = weights,
  sample_frac = samplingFracs[[2]],
  seed = 1
)

halfSamplingData <- biased_species_sample(
  data = true_species_ranges,
  weights = weights,
  sample_frac = samplingFracs[[3]],
  seed = 1
)

highSamplingData <- biased_species_sample(
  data = true_species_ranges,
  weights = weights,
  sample_frac = samplingFracs[[4]],
  seed = 1
)

higherSamplingData <- biased_species_sample(
  data = true_species_ranges,
  weights = weights,
  sample_frac = samplingFracs[[5]],
  seed = 1
)

highestSamplingData <- biased_species_sample(
  data = true_species_ranges,
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

# d. IN PROGRESSLow to high site sampling (N sites) of training datasets#######
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






















