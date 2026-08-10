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
  
  output_folder <- "./5-SensitivityTesting/output/trueDatasets/"
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(sf)
  library(purrr)
  library(ggplot2)
  library(fitdistrplus)
  library(MASS)
  library(e1071)
  source("./5-SensitivityTesting/scripts/utils/helpers.R")
  
}

animalDataPath<-"./5-SensitivityTesting/data/rawData/allGBIFCombined.19.6.2026.tsv"
if (file.exists(animalDataPath)) {
  message("The animal data already exists. You can skip to making the true datasets....")
} else {
  warning("Animal data not found. Continue running...")
  source("./5-SensitivityTesting/scripts/0-CompileCleanedData.R")
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
}
  ################Convert animal data to spatial###################
  ###This can be skipped if the output files already exist
  animalDataPath2<-"./5-SensitivityTesting/data/simFromPlant/animal_grid_data.RData"
  if (file.exists(animalDataPath2)) {
    message("The animal data already exists. You can skip to making the true datasets....")
  } else {
    warning("Animal data not found. Continue running...")
  animalDataPath <- "./5-SensitivityTesting/data/simFromPlant/simulatedOccurencesForAnimal.rds"
  
  presentAnimals<-readRDS(animalDataPath)
  
  head(presentAnimals)
  
  animal_sf <- presentAnimals
  animal_grid<-animal_sf
  
  # Map "scientificName" to "Morphotype" so downstream scripts work unchanged
  animal_grid$Morphotype <- animal_grid$scientificName 
  animal_grid$species <- animal_grid$scientificName 
  
  message("  Records within Amazon grid: ", nrow(animal_grid),
          " (", n_distinct(animal_grid$Morphotype), " species)")
  summary(animal_grid)
  
  save_path <- "./5-SensitivityTesting/data/simFromPlant/animal_grid_data.RData"
  
  # Save the object
  save(animal_grid, file = save_path)
  
  animal_endemism <- animal_grid |>
    group_by(grid_id) |>
    summarize(animal_richness = n_distinct(Morphotype), .groups = "drop")
  
  # Save the object
  save_path_End <- "./5-SensitivityTesting/data/simFromPlant/animal_grid_endemism.RData"
  
  save(animal_endemism, file = save_path_End)
}



################TRUTH 2: ENDEMISM RESAMPLES from complete true data to create true organisms with varying endemism######################

#We have two datasets so far:
#Sample - this is a combination of all the raw occurrence data from GBIF. It contains a variety of organisms and endemism types (high to low). We're calling this a "sample" because it is a sampling of the true-in-life unknown biodiversity. As such, it is our basis for parameterizing sample characteristics 
#True - this are the simulated occurence data based on the raw data (sample). We are going to consider this to be "true-in-simulation". It can serve as a base-line truth dataset and it can also serve as the dataset we draw species from to create other datasets.

load("./5-SensitivityTesting/data/rawData/animal_grid_data.RData")
sampleGBIF_animal_occ_data<-animal_grid
load("./5-SensitivityTesting/data/rawData/animal_grid_endemism.RData")
sampleGBIF_gridEndemism<-animal_endemism

rm(animal_grid)
rm(animal_endemism)

load("./5-SensitivityTesting/data/simFromPlant/animal_grid_data.RData")
load("./5-SensitivityTesting/data/simFromPlant/animal_grid_endemism.RData")
true_animal_occ_data<-animal_grid
true_gridEndemism<-animal_endemism



#Examine features of sampled dataset
#Calculate species grid cell occupancy
head(animal_grid)
head(true_animal_occ_data)
true_animal_occ_data$species<-true_animal_occ_data$scientificName

true_species_ranges <- true_animal_occ_data %>%
  st_drop_geometry() %>%
  distinct(species  , grid_id) %>%     # remove duplicate records within a grid
  count(species , name = "n_cells")


sampled_species_ranges <- sampleGBIF_animal_occ_data %>%
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
true_species_ranges
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
  sample_frac = 0.5,
  seed = 1
)   
plotData <- HighEndemismResult$sampled_data[, .(species, n_cells, quantile_class, sampled_bin)]

ggplot(plotData, aes(x = n_cells)) +
  geom_histogram(bins = 50) +
  scale_x_log10() +
  theme_bw() +
  labs(
    x = "Occupied grid cells (log scale)",
    y = "Number of species"
  )



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
  sample_frac = 0.5,
  seed = 1
)

weightKernel<-c(
  0,#high endemism weights
  0, #mid endemism weights
  0.0, #low endemism weights
  0.0 #very low endemism weights
)
#sample endemism weights
weights<-make_weights(weightKernel)
sum(weights)==1 ##should be TRUE

sampleEndemismResult <- biased_species_sample(
  data = true_species_ranges,
  weights = weights,
  sample_frac = 0.5,
  seed = 1
)

#True endemism weights
trueEndemismResult <- biased_species_sample(
  data = true_species_ranges,
  weights = unlist(c(weights_df$actual_prop[[1]], weights_df$actual_prop)),
  sample_frac = 0.5,
  seed = 1
)
unlist(c(weights_df$actual_prop[[1]], weights_df$actual_prop))/make_weights()


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
  sample_frac = 0.5,
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
  sample_frac = 0.5,
  seed = 1,
  smooth = FALSE
)
plotData <- veryLowEndemismResult$sampled_data[, .(species, n_cells, quantile_class, sampled_bin)]

ggplot(plotData, aes(x = n_cells)) +
  geom_histogram(bins = 50) +
  scale_x_log10() +
  theme_bw() +
  labs(
    x = "Occupied grid cells (log scale)",
    y = "Number of species"
  )



str(veryLowEndemismResult)
 
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

str(VeryLow)
true_animal_occ_data["species"]




#this will recreate occurrence data based on the sampling scenarios
VeryLow_comm <- rebuild_community(setDT(true_animal_occ_data), VeryLow, "very_low")
Low_comm <- rebuild_community(setDT(true_animal_occ_data), Low, "low")
Medium_comm <- rebuild_community(setDT(true_animal_occ_data), Medium, "medium")
High_comm <- rebuild_community(setDT(true_animal_occ_data), High, "high")
sample_comm<-rebuild_community(setDT(true_animal_occ_data), Sample_w, "sample")
true_comm<-rebuild_community(setDT(true_animal_occ_data), True_w, "true")

#we're going to make a summary table of all of these organisms, so we have to combine them into a single list
all_comm <- list(
  high = High_comm, 
  medium = Medium_comm,
  low = Low_comm,
  veryLow = VeryLow_comm,
  sample = sample_comm,
  true = true_comm
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
str(scenario_summary2)

scenario_summary2_clean <- scenario_summary2 %>%
  st_drop_geometry()


scenario_summary_full0 <- scenario_summary1 %>%
  left_join(scenario_summary2_clean, by = "scenario")

#output the summary info for these datasets
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
outfile <- paste0(output_folder,paste0("endemismSamples.summary_", timestamp, ".csv"))
write.csv(scenario_summary_full0, outfile, row.names = FALSE)


# Extract weights from all endemism classes
weights_list <- list(
  high = HighEndemismResult$weights_used,
  medium = mediumEndemismResult$weights_used,
  sample = sampleEndemismResult$weights_used,
  true = trueEndemismResult$weights_used,
  low = lowEndemismResult$weights_used,
  very_low = veryLowEndemismResult$weights_used
)

# Convert to data frame for easier viewing/export
weights_df2 <- map_df(weights_list, ~as.data.frame(as.list(.)), .id = "scenario") |>
  pivot_longer(cols = -scenario, names_to = "endemism_class", values_to = "weight")

# Export weights
outfile <- paste0(output_folder, paste0("EndemismWeights_", timestamp, ".csv"))
write.csv(weights_df2, outfile, row.names = FALSE)

# The occurrence datasets


#export each of the datasets
for (nm in names(all_comm)) {
  
  df <- all_comm[[nm]]
  
  outfile <- paste0(output_folder,paste0(nm, ".Endemism.sample_", timestamp, ".RData"))
  
  save(df, file = outfile)
  
  #write.csv(df, outfile, row.names = FALSE)
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


highestR_comm <- rebuild_community(setDT(true_animal_occ_data), highest_R, "highest")
highR_comm <- rebuild_community(true_animal_occ_data, high_R, "high")
mediumR_comm <- rebuild_community(true_animal_occ_data, medium_R, "medium")
lowR_comm <- rebuild_community(true_animal_occ_data, low_R, "low")
veryLowR_comm<-rebuild_community(true_animal_occ_data, veryLow_R, "veryLow")



all_comm_R <- list(
  highest = highestR_comm, 
  high = highR_comm,
  medium = mediumR_comm,
  low = lowR_comm,
  veryLow = veryLowR_comm
)

all_comm_R_flat <- bind_rows(all_comm_R)

names(all_comm_R)
all_comm_R_flat <- all_comm_R_flat %>%
  mutate(scenario = factor(
    scenario,
    levels = c("highest", "high", "medium", "low", "very_low"),
    ordered = TRUE
  ))

scenario_summary_Richness <- all_comm_R_flat %>%
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

scenario_df_R <- bind_rows(scenario_list_R, .id = "scenario")

#export the richness information
{outfile <- paste0(output_folder,paste0("richnessScenarios.", timestamp, ".csv"))
  
  write.csv(scenario_df_R, outfile, row.names = FALSE)}


#export each of the datasets
for (nm in names(all_comm_R)) {
  
  df <- all_comm_R[[nm]]
  
  outfile <- paste0(output_folder,paste0(nm, ".Richness.sample_", timestamp, ".RData"))
  save(df, file = outfile)
  #write.csv(df, outfile, row.names = FALSE)
}

