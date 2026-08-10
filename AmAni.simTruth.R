#This notebook will take the combined raw data in "C:\Users\dominice\OneDrive - University of Illinois - Urbana\Projects\Peruvian Blattodea\Cockroaches in a coal mine\AmAni\5-SensitivityTesting\data\rawData"
#and use the animal model to generate complete distributions for these species.
#This is an improvement from the raw data because the species in this dataset will 
#have more complete range data. Whereas the raw data are likely to have very spotty range data. 
#Since we're going to subsample this to simulate sampling from the true population, it makes
#sense to start out with data that are more complete.

source("AmAni.manager.R")
source("AmAni.summarize.R")

# --- All combined GBIF data ---
#data.table::fread("2-ModellingAnimals/data/allCombined/allGBIFCombined.19.6.2026.csv", nrows = 3)

  run_animal_model("allCombined", clean = FALSE, redo = FALSE,
                   min_occurrences = 1, min_completeness = 0,
                   min_richness = 1,
                   filter_plant_completeness = 0,
                   #fit_params = c("amplitude", "meanlog","sdlog" ),
                   data_file = "data/allCombined/allGBIFCombined.19.6.2026.tsv")
  
  
  
  