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
  
  data_folder <- "./5-SensitivityTesting/output/trueDatasets/"
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


######Import datasets######

#Import datasets 
list.files(data_folder, pattern = "\\.RData$")


# 1. Get the full file paths and the short file names
files_full <- list.files(data_folder, pattern = "\\.RData$", full.names = TRUE)
files_short <- list.files(data_folder, pattern = "\\.RData$", full.names = FALSE)

# 2. Loop through each file
for (i in seq_along(files_full)) {
  
  # Extract the first two terms (splitting by the period '.')
  # This grabs index 1 and 2, then pastes them back together with a dot
  name_parts <- unlist(strsplit(files_short[i], "\\."))
  custom_name <- paste(name_parts[2], name_parts[1], sep = ".")
  
  # Create a temporary sandboxed environment
  temp_env <- new.env()
  
  # Load the RData file into that temporary environment
  loaded_objects <- load(files_full[i], envir = temp_env)
  
  # Assuming there is only 1 main object saved inside each RData file:
  actual_object_name <- loaded_objects[1]
  
  # Assign the object from the temp environment to your Global Environment 
  # using your brand new custom name
  assign(custom_name, get(actual_object_name, envir = temp_env), envir = .GlobalEnv)
  
  # Clean up messaging (optional)
  message(paste("Loaded:", custom_name))
}

###The purpose of this document is to make sampled datasets based on simulated truth datasets that vary
###in their biological properties.
###We can think of this as simulating sampling of different types of organisms.
###Some organisms have very low spatial endemism (like many birds, especially larger ones)
###Some organisms have very high spatial endemism (like hemiptera, frogs, beetles, etc)
###Same goes for richness (high richness organisms like insects, low richness organisms like mammals)

###Does the model work better for one type of organism or another?


#2. Assessing sampling of endemism data############################

##Sampling scenarios 1: Organisms varying in their endemism level#######

scenario_levels <- c("Endemism.veryLow", "Endemism.low", "Endemism.medium", "Endemism.high", "Endemism.sample", "Endemism.true")
rm(animal_weighted_endemism)
scenario_dataSets <-list(Endemism.veryLow, Endemism.low, Endemism.medium, Endemism.high, Endemism.sample, Endemism.true)
# Loop over sampling fractions and create labeled results
sampling_results <- list()
for (i in seq_along(scenario_dataSets)) {
  result <- sample_Sites(
    true_animal_occ_data=scenario_dataSets[[i]],
    prop_species = 0.8,
    prop_sites = 1, 
    n_sites = 100,
    endemicity_sites = NULL, 
    richness_sites = NULL, 
    clustering = NULL, 
    animal_occ_Samples = sampleGBIF_animal_occ_data, 
    sampling_evenness = "even"
  )
  
  sampling_results[[scenario_levels[i]]] <- result
}

length(sampling_results)

####Combine results, summarize, and export
samplingLevels <- c("veryLow", "low", "medium", "high", "sample", "true")

results <- output_sample_results(
  sampling_results = sampling_results,
  scenario_levels  = scenario_levels,
  sampling_levels   = samplingLevels,
  true_occ_data    = true_animal_occ_data,
  output_folder    = output_folder,
  output_prefix    = "variEndemSampling"
)

#3. Assessing sampling of richness data############################

scenario_levels <- c("Richness.veryLow", "Richness.low", "Richness.medium", "Richness.high", "Richness.highest")
scenario_dataSets <-list(Richness.veryLow, Richness.low, Richness.medium, Richness.high, Richness.highest)
# Loop over sampling fractions and create labeled results
sampling_results <- list()
for (i in seq_along(scenario_dataSets)) {
  result <- sample_Sites(
    true_animal_occ_data=scenario_dataSets[[i]],
    prop_species = 0.8,
    prop_sites = 1, 
    n_sites = 100,
    endemicity_sites = NULL, 
    richness_sites = NULL, 
    clustering = NULL, 
    animal_occ_Samples = sampleGBIF_animal_occ_data, 
    sampling_evenness = "even"
  )
  
  sampling_results[[scenario_levels[i]]] <- result
}

length(sampling_results)

####Combine results, summarize, and export
samplingLevels <- c("veryLow", "low", "medium", "high", "highest")


results <- output_sample_results(
  sampling_results = sampling_results,
  scenario_levels  = scenario_levels,
  sampling_levels   = samplingLevels,
  true_occ_data    = true_animal_occ_data,
  output_folder    = output_folder,
  output_prefix    = "variRichSampling"
)


