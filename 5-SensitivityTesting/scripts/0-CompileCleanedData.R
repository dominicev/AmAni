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

write.table(allAnimalDatasets[-1], file = "./5-SensitivityTesting/data/allGBIFCombined_full.19.6.2026.tsv", sep = "\t", row.names = FALSE)


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


head(presentAnimals[-1])

write.table(presentAnimals[-1], file = "./5-SensitivityTesting/data/allGBIFCombined.19.6.2026.tsv", sep = "\t", row.names = FALSE)
write.table(presentAnimals[-1], file = "./5-SensitivityTesting/data/allGBIFCombined.19.6.2026.csv", sep = ",", row.names = FALSE)
