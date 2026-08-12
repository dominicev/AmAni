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
  
  output_folder <- "./5-SensitivityTesting/output/sampleDatasets/"
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
  source("./5-SensitivityTesting/scripts/utils/samplingFunctions.R")
  source("./5-SensitivityTesting/scripts/utils/outputSample.R")
  
  
}

#Here is a list of different sampling conditions we can assess

#Sampling conditions:
 #DONE #Number of species - sampling species based on occurrence, we are unlikely to sampling 100% of all taxa. Locally common species will be sampled first, and locally rare species will be sampled last.
 
 #DONE #Proportion of sites - even GBIF data are unlikely to have all sites sampled. Accessible sites, or ones with infrastructure are more likely to be sampled. Parameter site sampling probability based on GBIF data
#DONE #Endemicity of sites - certain sites are likely to have more endemic species than others.
#DONE #Richness of sites - certain sites will be more diverse than others
#DONE  #Number of sites - Curated datasets are unlikely to have many sites sampled. This is different from "Proportion of sites" because it will rely on 3-20 sites only .
 #DONE #Clustering of sites - again, curated datasets are likely to have spatially clustered sites.
  #DONE #Sampling evenness - any datasets are unlikely to be sampled totally evenly (some sites a higher % of species are sampled, others, a lower %)

######Import datasets######

#Import some testing datasets (will eventually import other datasets)
load("./5-SensitivityTesting/data/rawData/animal_grid_data.RData")
sampleGBIF_animal_occ_data<-animal_grid
load("./5-SensitivityTesting/data/rawData/animal_grid_endemism.RData")
sampleGBIF_gridRichness<-animal_endemism
animal_occ_Samples<-sampleGBIF_animal_occ_data
rm(animal_grid)
rm(animal_endemism)

load("./5-SensitivityTesting/data/simFromPlant/animal_grid_data.RData")
load("./5-SensitivityTesting/data/simFromPlant/animal_grid_endemism.RData")
true_animal_occ_data<-animal_grid
true_gridRichness<-animal_endemism
true_animal_occ_data$species<-true_animal_occ_data$scientificName


rm(animal_grid)
rm(animal_endemism)

#testing sample function
#richness_sites: Options: low, medium, high; Default: NULL
#endemicity_sites: Options: low, medium, high; Default: NULL
##clustering: Options: low, medium, high, GBIF; Default: NULL
##sampling_evenness: Options: even (Default),veryLow, low
##
{
test<-sample_Sites(true_animal_occ_data, prop_species = 1,
             prop_sites = NULL, 
             n_sites = 20,
             endemicity_sites = NULL , 
             richness_sites = NULL, 
             clustering = "medium", 
             animal_occ_Samples = sampleGBIF_animal_occ_data, 
             sampling_evenness = "low")
  summary(test)
  str(test)
  attr(test, "min_weight")
  testData<-test$occurenceData
  # Use data.table for speed: drop geometry column, compute in two passes
  occ_dt <- as.data.table(st_drop_geometry(testData))
  occ_dt[, range_size   := .N, by = species]
  occ_dt[, species_weight := 1 / range_size]
  animal_weighted_endemism <- occ_dt[, .(
    weighted_endemism = sum(species_weight),
    species_richness  = .N
  ), by = grid_id]
  
  testPlotData<-left_join(animal_weighted_endemism, testData, by = "grid_id")
  
  testPlotData
length(unique(testPlotData$grid_id))

all_coords <- st_coordinates(st_centroid(testPlotData$geometry))

plot_df <- as.data.frame(all_coords) |>
  mutate(
    grid_id = testPlotData$grid_id,
    richness = testPlotData$species_richness,
    endemism = testPlotData$weighted_endemism  
    
    ) |>
  distinct(grid_id, .keep_all = TRUE) |>
  st_as_sf(coords = c("X", "Y"), crs = st_crs(testPlotData)) |>
  st_buffer(dist = 21180, nQuadSegs = 1)

#plot 1
ggplot() +
  geom_sf(data = plot_df, aes(fill = richness), color = NA) +
  scale_fill_viridis_c(option = "magma", name = "weight") +
  labs(
    title = paste0("sampled sites"),
    x = NULL, y = NULL
  ) +
  theme_minimal()
}

{
  #plot 2
ggplot() +
  geom_sf(data = plot_df, aes(fill = endemism), color = NA) +
  scale_fill_viridis_c(option = "magma", name = "weight") +
  labs(
    title = paste0("sampled sites"),
    x = NULL, y = NULL
  ) +
  theme_minimal()
}


#1. Assessing sampling of baseline data############################

#First, we want to assess how the different sampling parameters affect the accuracy of the model. After that, we want to assess how a realistic but effective combination of sampling parameters perform on datasets with different biological features (endemism, richness, etc)
#
#In this section, we will create sample replicates of a baseline dataset: the full true occurence data
#
##Sampling scenarios 1: Low to high total sampling (% of species) training datasets#######


scenario_levels <- c("min", "quarter", "half", "high", "higher", "highest")

samplingFracs <- c(
  0.1,0.25, 0.5, 0.75, 0.85, 0.95
)

# Loop over sampling fractions and create labeled results
sampling_results <- list()

for (i in seq_along(samplingFracs)) {
  result <- sample_Sites(
    true_animal_occ_data,
    prop_species = samplingFracs[i],
    prop_sites = 1, 
    n_sites = NULL,
    endemicity_sites = NULL, 
    richness_sites = NULL, 
    clustering = NULL, 
    animal_occ_Samples = sampleGBIF_animal_occ_data, 
    sampling_evenness = "even"
  )
  
  sampling_results[[scenario_levels[i]]] <- result
}

sampling_results

####Combine results, summarize, and export
#debug(output_sample_results)
results <- output_sample_results(
  sampling_results = sampling_results,
  scenario_levels  = scenario_levels,
  sampling_levels   = samplingFracs,
  true_occ_data    = true_animal_occ_data,
  output_folder    = output_folder,
  output_prefix    = "lowtoHighSpeciesSampling"
)

##Sampling scenarios 2: Low to high site sampling (N sites) of training datasets#######
#In this sampling regime, I am going to thin the data by site, instead of by taxon. This will create sampled datasets.


samplingLevels <- c(
   10,10,10,10,10, 20,20,20,20,20,20, 40,40,40, 80,80,80, 160,160, 320, 640, 1280)
scenario_levels <- lapply(as.character(  samplingLevels), function(x) paste0(x,"_", seq_along(samplingLevels)))
# Convert list to a matrix
mat <- do.call(rbind, scenario_levels)
scenario_Names <- diag(mat)

length(unique(sampleGBIF_animal_occ_data$grid_id))

# Loop over sampling fractions and create labeled results
sampling_results <- list()

for (i in seq_along(samplingLevels)) {
  result <- sample_Sites(
    true_animal_occ_data,
    prop_species = 1,
    prop_sites = 1, 
    n_sites = samplingLevels[i],
    endemicity_sites = NULL, 
    richness_sites = NULL, 
    clustering = NULL, 
    animal_occ_Samples = sampleGBIF_animal_occ_data, 
    sampling_evenness = "even"
  )
  
  sampling_results[[scenario_Names[i]]] <- result
}

length(sampling_results)
length(samplingLevels)

####Combine results, summarize, and export

results <- output_sample_results(
  sampling_results = sampling_results,
  scenario_levels  = scenario_Names,
  sampling_levels   = samplingLevels,
  true_occ_data    = true_animal_occ_data,
  output_folder    = output_folder,
  output_prefix    = "lowtoHighSiteSampling"
)

##Sampling scenarios 3: Sampling from a range of site endemicities in training datasets#######
#What if the sites we sample at don't represent the full range of endemicity? What if they are all very high, or very low endemicity?

#scenario_levels <- c("low", "medium", "high")
scenario_levels <- c("low", "low", "medium", "medium",  "high", "high")
samplingNames <- c("low1", "low2", "medium1", "medium2", "high1", "high2")
#endemicity_sites: Options: low, medium, high; Default: NULL

samplingLevels <- scenario_levels

# Loop over sampling fractions and create labeled results
sampling_results <- list()

for (i in seq_along(samplingLevels)) {
  result <- sample_Sites(
    true_animal_occ_data,
    prop_species = 1,
    prop_sites = 1, 
    n_sites = 200,
    endemicity_sites = samplingLevels[i], 
    richness_sites = NULL, 
    clustering = NULL, 
    animal_occ_Samples = sampleGBIF_animal_occ_data, 
    sampling_evenness = "even"
  )
  
  sampling_results[[samplingNames[i]]] <- result
}

sampling_results

####Combine results, summarize, and export

results <- output_sample_results(
  sampling_results = sampling_results,
  scenario_levels  = scenario_levels,
  sampling_levels   = samplingLevels,
  true_occ_data    = true_animal_occ_data,
  output_folder    = output_folder,
  output_prefix    = "lowtoHighSiteEnd"
)

##Sampling scenarios 4:  Sampling from a range of site richness in training datasets#######
#What if the sites we sample at don't represent the full range of richness? What if they are all very high, or very low richness?

scenario_levels <- c("low", "low", "medium", "medium",  "high", "high")
samplingNames <- c("low1", "low2", "medium1", "medium2", "high1", "high2")
#endemicity_sites: Options: low, medium, high; Default: NULL

samplingLevels <- scenario_levels

# Loop over sampling fractions and create labeled results
sampling_results <- list()

for (i in seq_along(samplingLevels)) {
  result <- sample_Sites(
    true_animal_occ_data,
    prop_species = 1,
    prop_sites = 1, 
    n_sites = 200,
    endemicity_sites = NULL, 
    richness_sites = samplingLevels[i], 
    clustering = NULL, 
    animal_occ_Samples = sampleGBIF_animal_occ_data, 
    sampling_evenness = "even"
  )
  
  sampling_results[[samplingNames[i]]] <- result
}

sampling_results

####Combine results, summarize, and export

results <- output_sample_results(
  sampling_results = sampling_results,
  scenario_levels  = scenario_levels,
  sampling_levels   = samplingLevels,
  true_occ_data    = true_animal_occ_data,
  output_folder    = output_folder,
  output_prefix    = "lowtoHighSiteRichness"
)



##Sampling scenarios 5:  Clustering of sites sampled in training datasets#######
#What if the sites we sample are clustered together, rather than randomly distributed across the amazon?
#I am going to run this multiple times at different levels to see how much they vary
scenario_levels <- c("low", "low", "low","low", "medium", "medium", "medium","medium", "high","high","high","high","high", "GBIF")
samplingNames <- c("low1", "low2", "low3","low4", "medium1", "medium2", "medium3","medium4", "high1","high2", "high3", "high4", "high5","GBIF")

samplingLevels <- scenario_levels

# Loop over sampling fractions and create labeled results
sampling_results <- list()

  for (i in seq_along(samplingLevels)) {
  result <- sample_Sites(
    true_animal_occ_data,
    prop_species = 1,
    prop_sites = 1, 
    n_sites = 200,
    endemicity_sites = NULL, 
    richness_sites = NULL, 
    clustering = samplingLevels[i], 
    animal_occ_Samples = sampleGBIF_animal_occ_data, 
    sampling_evenness = "even"
  )
  
  sampling_results[[samplingNames[i]]] <- result
  }

length(sampling_results)

####Combine results, summarize, and export
scenario_names <- c("low", "medium", "high", "GBIF")

samplingLevels <- c("low1", "low2", "low3","low4", "medium1", "medium2", "medium3","medium4", "high1","high2", "high3", "high4", "high5","GBIF")

results <- output_sample_results(
  sampling_results = sampling_results,
  scenario_levels  = scenario_levels,
  sampling_levels   = samplingLevels,
  true_occ_data    = true_animal_occ_data,
  output_folder    = output_folder,
  output_prefix    = "variedSiteClustering"
)


##Sampling scenarios 6: Sampling evenness by site in training datasets#######
#What if the sites we sample aren't sampled evenly? Some sites have a high proportion of species sampled, and other sites have a low propotion of species sampled.
#This one I am going to this twice, once with 10 sites sampled, and again with 1000 sites sampled. This is because I expect the sample size of sites to affect the impact of evenness

##sampling_evenness: Options: even (Default),veryLow, low

scenario_levels <- c("even", "veryLow", "low")
scenario_names <- c("even1", "veryLow1", "low1", "even2", "veryLow2", "low2")

samplingLevels <- scenario_levels

# Loop over sampling fractions and create labeled results
sampling_results <- list()

for (i in seq_along(samplingLevels)) {
  result <- sample_Sites(
    true_animal_occ_data,
    prop_species = 1,
    prop_sites = 1, 
    n_sites = 50,
    endemicity_sites = NULL, 
    richness_sites = NULL, 
    clustering = NULL, 
    animal_occ_Samples = sampleGBIF_animal_occ_data, 
    sampling_evenness = samplingLevels[i]
  )
  
  sampling_results[[scenario_names[i]]] <- result
}


for (i in seq_along(samplingLevels)) {
  result <- sample_Sites(
    true_animal_occ_data,
    prop_species = 1,
    prop_sites = 1, 
    n_sites = 1000,
    endemicity_sites = NULL, 
    richness_sites = NULL, 
    clustering = NULL, 
    animal_occ_Samples = sampleGBIF_animal_occ_data, 
    sampling_evenness = samplingLevels[i]
  )
  
  sampling_results[[scenario_names[i+length(samplingLevels)]]] <- result
}


length(sampling_results)

####Combine results, summarize, and export

scenario_levels <- c("even", "veryLow", "low", "even", "veryLow", "low")



results <- output_sample_results(
  sampling_results = sampling_results,
  scenario_levels  = scenario_levels,
  sampling_levels   = scenario_names,
  true_occ_data    = true_animal_occ_data,
  output_folder    = output_folder,
  output_prefix    = "variedSiteEvenness"
)



##Sampling scenarios 7: cockroach MDD  training datasets#######
#Here I am going to replicate something similar to the sampling qualities of the MDD cockroach data


scenario_levels<-c("blatt")
samplingLevels<-scenario_levels
sampling_results <- list()
  result <- sample_Sites(
    true_animal_occ_data,
    prop_species = 0.7,
    prop_sites = 1, 
    n_sites = 20,
    endemicity_sites = NULL, 
    richness_sites = NULL, 
    clustering = "high", 
    animal_occ_Samples = sampleGBIF_animal_occ_data, 
    sampling_evenness = "even"
  )
  
  sampling_results[["blatt"]] <- result

####Combine results, summarize, and export

results <- output_sample_results(
  sampling_results = sampling_results,
  scenario_levels  = scenario_levels,
  sampling_levels   = samplingLevels,
  true_occ_data    = true_animal_occ_data,
  output_folder    = output_folder,
  output_prefix    = "roachLike"
)

scenario_levels<-c("blatt")
samplingLevels<-scenario_levels
sampling_results <- list()
result <- sample_Sites(
  true_animal_occ_data,
  prop_species = 0.8,
  prop_sites = 1, 
  n_sites = 20,
  endemicity_sites = NULL, 
  richness_sites = NULL, 
  clustering = "high", 
  animal_occ_Samples = sampleGBIF_animal_occ_data, 
  sampling_evenness = "low"
)

sampling_results[["blatt"]] <- result

####Combine results, summarize, and export

results <- output_sample_results(
  sampling_results = sampling_results,
  scenario_levels  = scenario_levels,
  sampling_levels   = samplingLevels,
  true_occ_data    = true_animal_occ_data,
  output_folder    = output_folder,
  output_prefix    = "roachLike2"
)


scenario_levels<-c("blatt")
samplingLevels<-scenario_levels
sampling_results <- list()
result <- sample_Sites(
  true_animal_occ_data,
  prop_species = 0.6,
  prop_sites = 1, 
  n_sites = 20,
  endemicity_sites = NULL, 
  richness_sites = NULL, 
  clustering = "medium", 
  animal_occ_Samples = sampleGBIF_animal_occ_data, 
  sampling_evenness = "low"
)

sampling_results[["blatt"]] <- result

####Combine results, summarize, and export

results <- output_sample_results(
  sampling_results = sampling_results,
  scenario_levels  = scenario_levels,
  sampling_levels   = samplingLevels,
  true_occ_data    = true_animal_occ_data,
  output_folder    = output_folder,
  output_prefix    = "roachLike3"
)

##Sampling scenarios 8: BEETLE GBIF  training datasets#######
#Here I am going to replicate something similar to the sampling qualities of the GBIF beetle data

scenario_levels<-c("coleo")
samplingLevels<-scenario_levels
sampling_results <- list()

result <- sample_Sites(
  true_animal_occ_data,
  prop_species = 0.6,
  prop_sites = 1, 
  n_sites = 1400,
  endemicity_sites = NULL, 
  richness_sites = NULL, 
  clustering = "low", 
  animal_occ_Samples = sampleGBIF_animal_occ_data, 
  sampling_evenness = "veryLow"
)

sampling_results[[scenario_levels[1]]] <- result

####Combine results, summarize, and export

results <- output_sample_results(
  sampling_results = sampling_results,
  scenario_levels  = scenario_levels,
  sampling_levels   = samplingLevels,
  true_occ_data    = true_animal_occ_data,
  output_folder    = output_folder,
  output_prefix    = "beetleLike"
)


scenario_levels<-c("coleo")
samplingLevels<-scenario_levels
sampling_results <- list()

result <- sample_Sites(
  true_animal_occ_data,
  prop_species = 0.4,
  prop_sites = 1, 
  n_sites = 1000,
  endemicity_sites = NULL, 
  richness_sites = NULL, 
  clustering = "low", 
  animal_occ_Samples = sampleGBIF_animal_occ_data, 
  sampling_evenness = "veryLow"
)

sampling_results[[scenario_levels[1]]] <- result

####Combine results, summarize, and export

results <- output_sample_results(
  sampling_results = sampling_results,
  scenario_levels  = scenario_levels,
  sampling_levels   = samplingLevels,
  true_occ_data    = true_animal_occ_data,
  output_folder    = output_folder,
  output_prefix    = "beetleLike2"
)

scenario_levels<-c("coleo")
samplingLevels<-scenario_levels
sampling_results <- list()

result <- sample_Sites(
  true_animal_occ_data,
  prop_species = 0.3,
  prop_sites = 1, 
  n_sites = 2000,
  endemicity_sites = NULL, 
  richness_sites = NULL, 
  clustering = "low", 
  animal_occ_Samples = sampleGBIF_animal_occ_data, 
  sampling_evenness = "veryLow"
)

sampling_results[[scenario_levels[1]]] <- result

####Combine results, summarize, and export

results <- output_sample_results(
  sampling_results = sampling_results,
  scenario_levels  = scenario_levels,
  sampling_levels   = samplingLevels,
  true_occ_data    = true_animal_occ_data,
  output_folder    = output_folder,
  output_prefix    = "beetleLike3"
)

##Sampling scenarios 8: Butterfly GBIF  training datasets#######
#Here I am going to replicate something similar to the sampling qualities of the GBIF beetle data

scenario_levels<-c("butter")
samplingLevels<-scenario_levels
sampling_results <- list()

result <- sample_Sites(
  true_animal_occ_data,
  prop_species = 0.80,
  prop_sites = 1, 
  n_sites = 2000,
  endemicity_sites = NULL, 
  richness_sites = NULL, 
  clustering = "low", 
  animal_occ_Samples = sampleGBIF_animal_occ_data, 
  sampling_evenness = "veryLow"
)

sampling_results[[scenario_levels[1]]] <- result

####Combine results, summarize, and export

results <- output_sample_results(
  sampling_results = sampling_results,
  scenario_levels  = scenario_levels,
  sampling_levels   = samplingLevels,
  true_occ_data    = true_animal_occ_data,
  output_folder    = output_folder,
  output_prefix    = "butterflyLike"
)

scenario_levels<-c("butter")
samplingLevels<-scenario_levels
sampling_results <- list()

result <- sample_Sites(
  true_animal_occ_data,
  prop_species = 0.60,
  prop_sites = 1, 
  n_sites = 30,
  endemicity_sites = NULL, 
  richness_sites = NULL, 
  clustering = "medium", 
  animal_occ_Samples = sampleGBIF_animal_occ_data, 
  sampling_evenness = "low"
)

sampling_results[[scenario_levels[1]]] <- result

####Combine results, summarize, and export

results <- output_sample_results(
  sampling_results = sampling_results,
  scenario_levels  = scenario_levels,
  sampling_levels   = samplingLevels,
  true_occ_data    = true_animal_occ_data,
  output_folder    = output_folder,
  output_prefix    = "butterflyLike2"
)

scenario_levels<-c("butter")
samplingLevels<-scenario_levels
sampling_results <- list()

result <- sample_Sites(
  true_animal_occ_data,
  prop_species = 0.50,
  prop_sites = 1, 
  n_sites = 300,
  endemicity_sites = NULL, 
  richness_sites = NULL, 
  clustering = "medium", 
  animal_occ_Samples = sampleGBIF_animal_occ_data, 
  sampling_evenness = "low"
)

sampling_results[[scenario_levels[1]]] <- result

####Combine results, summarize, and export

results <- output_sample_results(
  sampling_results = sampling_results,
  scenario_levels  = scenario_levels,
  sampling_levels   = samplingLevels,
  true_occ_data    = true_animal_occ_data,
  output_folder    = output_folder,
  output_prefix    = "butterflyLike3"
)


scenario_levels<-c("butter")
samplingLevels<-scenario_levels
sampling_results <- list()

result <- sample_Sites(
  true_animal_occ_data,
  prop_species = 0.80,
  prop_sites = 1, 
  n_sites = 20,
  endemicity_sites = NULL, 
  richness_sites = NULL, 
  clustering = "high", 
  animal_occ_Samples = sampleGBIF_animal_occ_data, 
  sampling_evenness = "even"
)

sampling_results[[scenario_levels[1]]] <- result

####Combine results, summarize, and export

results <- output_sample_results(
  sampling_results = sampling_results,
  scenario_levels  = scenario_levels,
  sampling_levels   = samplingLevels,
  true_occ_data    = true_animal_occ_data,
  output_folder    = output_folder,
  output_prefix    = "butterflyLike4"
)

















