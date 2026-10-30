#run_animal_model notes based on simulation results
# run_animal_model(folder_name,
#                  data_file = sub("^.*?data/", "data/", full_path),
#                  useSamplingExtension = TRUE, #adds the linear model extension based on sampling parameter
#                  clean = FALSE,
#                  redo = FALSE,
#                  # Completeness floor applied upstream via cell subsetting,
#                  # so the model's own filters stay disabled.
#                  min_occurrences = 1,
#                  min_completeness = 0,
#                  min_richness = 1,
#                  estimateBounds=FALSE,
#                  simulate_ranges = FALSE,
#                  skipRawEndemism = TRUE,
#                  customOutputFolderPrefix = "p_",
#                  gammaErrorType = "narrow", #The alternative "full" propogates the error from the log-normal model downstream into two more simulations. "narrow" just takes the error ranges of linear model from the best estimate of the surrogate model. Although "full" is niavely more appropriate, simulations were done on a situation akin to "narrow" and it was accurate >80% of the time when fill_plants = TRUE.
#                  disable_outlier_filter = TRUE,
#                  customOutputFolderSuffix = paste0("_",toString(set),"_",floor_tag),
#                  fill_plants=FALSE, #1. fill_plants does seem to improve accuracy in simulation while fill_animals does not.
#                                     #This might be because the animals are much more poorly sampled than the plants (in the simulation, this is by design), and so the filled animal ranges are unreliable. 
#                                     #2. If you set fill_plants = TRUE and fill_animals = TRUE at the same, set fit_params = c("amplitude"). You're artificially forcing spatial extents to correlate between the datasets.
#                                     #3. If you are highly undersampling endemic sites and highly oversampling low endemism sites, then set fill_plants = TRUE. If you over sample high endemism sites (e.g., cloud forests), then fill_plants = TRUE will increase error.
#                  fill_animals=FALSE, #You might want to turn this on if you think you have very good sampling completeness in sites. Important: See above comments for fill_plants.
#                  fit_params = c("amplitude", "meanlog", "sdlog"),
#                  filter_plant_completeness = NULL
# )
# 
###Notes:
###
###Things to do soon:
### 1. Make final figures MOSTLY DONE
### 3. Incorporate ant data DONE
### 4. Implement plant-animal correlation maximizer
### 5. Include retained data plots in organism summary<<<DONE
### 6. Fix bug in max min confidence intervals
### 7. Sync with one drive
### 8. Update on github
### 9. Add option to remove morphotypes
### 10. Check why no vulnerable ot endangered species are being found
### 11. Fix filter maps
### 
### 
###Things to do later:
### 1. New simulations for linear model
### 2. Reconcile probable overestimation of gamma in empirical data and consistent underestimation of gamma in simulation
### 3. Differentiate Chao2 completeness from Chao1 completeness in simulations (add abundances?)
### 4. Incorporate new Peru Blattodea data
### 5. Use iChao1 instead of Chao1?
### 6. Consider options to speed up large range simulation (i.e., parallelizing the small species ranges individually first, and then breaking up inndividual large ranges across multple cores )....i noticed that moths parallelizes but flies didn't...so maybe this can be addressed



# AmAni.run.R — Run the full pipeline for each organism
source("AmAni.manager.R")
coresNow<-55


#----Initialize some toy data (based on cockroach occurences) for troubleshooting (if needed)
#
# data_file <- "data/DarwinCoreAnimalOccurences.8-1-2026.xlsx"
# fill_plants <- FALSE
# fill_animals <- FALSE
# redo <- FALSE
# min_occurrences <- 50
# min_completeness <- 0.2
# min_richness <- 50
# max_richness <- 200
# fit_params <- c("amplitude", "meanlog","sdlog" )
#
#
#
#


# --- Cockroaches (Madre de Dios) ---
# Curated dataset — no cell filtering needed (min_occurrences = 0 disables it)
# fit_params controls which log-normal parameters are estimated from the animal
# MDD data when scaling the plant spatial extent distribution:
#   c("amplitude")                      — default (current behaviour, amplitude only)
#   c("amplitude", "meanlog")           — also shift location (range sizes)
#   c("amplitude", "meanlog", "sdlog")  — fit all three (location + spread)
check_status("cockroaches")
run_animal_model("cockroaches", clean = FALSE, 
                 min_occurrences = 1, 
                 min_completeness = 0, 
                 min_richness = 1,
                 fill_plants = TRUE,
                 fill_animals = FALSE,
                 useSamplingExtension=TRUE,
                 gammaErrorType = "narrow",
                 redo = FALSE,
                 data_file = "data/cockroaches_curated/cockroaches_curated.csv",
                 fix_taxon_names = FALSE,
                 simulateParallel=TRUE,
                 simulateNCores=coresNow,
                 simulate_ranges = "SkipUpper",
                 skipRawEndemism = FALSE,
                 disable_outlier_filter = TRUE,
                 fit_params = c("amplitude", "meanlog","sdlog" ))
print("=== Completed run_animal_model: cockroaches ===")


# --- Cockroaches, all South America (GBIF) ---
# Very sparse GBIF data (267 spp, ~5k records). Auto filter likely triggers
# the safety guard (< 30 cells), so filtering is effectively skipped.
# check_status("allRoaches")
# run_animal_model("allRoaches", clean = TRUE, redo = FALSE,
#                  min_occurrences = "auto", min_completeness = "auto",
#                  min_richness = 15,
#                  redo = FALSE,
#                  data_file = "data/cockroaches all s america/0011428-250402121839773.csv")
# run_extinction_model("allRoaches")

# --- Butterflies (Curated and GBIF combined) ---
# Curated dataset — no cell filtering needed (min_occurrences = 0 disables it)
# fit_params controls which log-normal parameters are estimated from the animal
#  c("amplitude")                      — default (current behaviour, amplitude only)
#   c("amplitude", "meanlog")           — also shift location (range sizes)
#   c("amplitude", "meanlog", "sdlog")  — fit all three (location + spread)





# --- Odonata (GBIF) ---
source("AmAni.manager.R")
check_status("odonata")
run_animal_model("odonata", clean = FALSE, redo = FALSE,
                 animalFilterType = "extraStrict",
                 completenessMethod = "Chao1",
                 min_completeness = 0.6,
                 min_occurrences=4,
                 min_richness=20,
                 max_richness=135,
                 fill_plants = TRUE,
                 useSamplingExtension=TRUE,
                 fix_taxon_names = FALSE,
                 simulate_ranges = "SkipUpper",
                 simulateNCores=coresNow,
                 simulateParallel=TRUE,
                 gammaErrorType = "narrow",

                 fit_params = c("amplitude", "meanlog","sdlog" ),
                 data_file = "data/odonata/0059398-260226173443078_cleaned_taxaFixed.csv")
print("=== Completed run_animal_model: odonata ===")

# --- Odonata (Curated) ---
check_status("odonata-curated")
run_animal_model("odonata-curated", clean = FALSE, redo = FALSE,
                 animalFilterType = "extraStrict",
                  completenessMethod ="Chao1",
                 min_completeness = 0.70,
                 min_occurrences=4,
                 min_richness=20,
                 max_richness=135,
                 #disable_outlier_filter=TRUE,
                 fill_plants = TRUE,
                 fill_animals = FALSE,
                 useSamplingExtension=TRUE,
                 gammaErrorType = "narrow",
                 fix_taxon_names = FALSE,
                 simulateParallel=TRUE,
                 simulateNCores=coresNow,
                 simulate_ranges = "SkipUpper",
                 skipRawEndemism = FALSE,
                 fit_params = c("amplitude", "meanlog","sdlog" ),
                 data_file = "data/odonata_curated/allCombined_CuratedOccurences_taxaFixed.csv")
print("=== Completed run_animal_model: odonata-curated ===")

# --- Odonata (COMBINED) ---
check_status("odonata-combined")
run_animal_model("odonata-combined", clean = FALSE, redo = FALSE,
                 animalFilterType = "extraStrict",
                 completenessMethod ="Chao1",
                 min_completeness = 0.75,
                 min_occurrences=4,
                 min_richness=20,
                 max_richness=135,
                 fill_plants = TRUE,
                 fill_animals = FALSE,
                 useSamplingExtension=TRUE,
                 gammaErrorType = "narrow",

                 simulateParallel=TRUE,
                 simulateNCores=coresNow,
                 simulate_ranges = "SkipUpper",
                 skipRawEndemism = FALSE,
                 fit_params = c("amplitude", "meanlog","sdlog" ),
                 
                 fix_taxon_names = FALSE,
                 data_file = "data/odonata-combined/combined_odonata-combined_taxaFixed.csv")
print("=== Completed run_animal_model: odonata-combined ===")
                 
                 #fix_taxon_names = TRUE,
                 #secondary_data_file = "data/odonata/0059398-260226173443078_cleaned_taxaFixed.csv",
                 #data_file = "data/odonata_curated/allCombined_CuratedOccurences_taxaFixed.csv")
# could replace the above with the single input file data/odonata-combined/combined_odonata-combined_taxaFixed.csv
# but if any of the datasets get updated it is necessary to recombine them

# # --- Free memory between organisms ---
rm(list = ls())          # wipes the global env, incl. run_animal_model and coresNow
gc(full = TRUE)          # return the memory to the OS
source("AmAni.manager.R")  # reload the pipeline functions (rm removed them)
coresNow <- 55           # re-set the core count (rm removed it too)

# --- Hemiptera (GBIF) ---
check_status("hemiptera")
run_animal_model("hemiptera", clean = FALSE, redo = FALSE,
                 animalFilterType = "extraStrict",
                 completenessMethod ="Chao1",
                 min_completeness = 0.5,
                 min_occurrences=5,
                 min_richness=30,
                 
                 disable_outlier_filter = TRUE, #I am going to shut this off for the diverse and under observed groups
                 fill_plants = TRUE,
                 useSamplingExtension=TRUE,
                 gammaErrorType = "narrow",
                 fix_taxon_names = FALSE,
                 simulateNCores=coresNow,
                 simulateParallel=TRUE,
                 fill_animals = FALSE,
                 simulate_ranges = "SkipUpper",
                 fit_params = c("amplitude", "meanlog","sdlog" ),
                 data_file = "data/hemiptera/0059425-260226173443078_cleaned_taxaFixed.csv")
print("=== Completed run_animal_model: hemiptera ===")


# --- Free memory between organisms ---
rm(list = ls())          # wipes the global env, incl. run_animal_model and coresNow
gc(full = TRUE)          # return the memory to the OS
source("AmAni.manager.R")  # reload the pipeline functions (rm removed them)
coresNow <- 55           # re-set the core count (rm removed it too)


# --- Lepidoptera (GBIF) ---
check_status("lepidoptera")
run_animal_model("lepidoptera", clean = FALSE, redo = FALSE,
                 animalFilterType = "extraStrict",
                 completenessMethod ="Chao1",
                 min_completeness = 0.6,
                 min_occurrences=5,
                 min_richness=50,
                 
                 disable_outlier_filter = TRUE, #I am going to shut this off for the diverse and under observed groups
                 fill_plants = TRUE,
                 useSamplingExtension=TRUE,
                 gammaErrorType = "narrow",
                 fill_animals = FALSE,
                 fix_taxon_names = FALSE,
                 simulateNCores=coresNow,
                 simulateParallel=TRUE,
                 simulate_ranges = "SkipUpper",
                 
                 fit_params = c("amplitude", "meanlog","sdlog" ),
                 data_file = "data/lepidoptera/0060839-260226173443078_cleaned_taxaFixed.csv")
print("=== Completed run_animal_model: lepidoptera ===")

# --- Free memory between organisms ---
rm(list = ls())          # wipes the global env, incl. run_animal_model and coresNow
gc(full = TRUE)          # return the memory to the OS
source("AmAni.manager.R")  # reload the pipeline functions (rm removed them)
coresNow <- 55           # re-set the core count (rm removed it too)


# --- Lepidoptera combined no butterflies (GBIF+curated) ---
check_status("moths")
run_animal_model("moths", clean = FALSE, redo = FALSE,
                 animalFilterType = "extraStrict",
                 completenessMethod ="Chao1",
                 min_completeness = 0.4,
                 min_occurrences=30,
                 min_richness=50,
                 disable_outlier_filter = TRUE, #I am going to shut this off for the diverse and under observed groups
                 fill_plants = TRUE,
                 useSamplingExtension=TRUE,
                 gammaErrorType = "narrow",
                 fill_animals = FALSE,
                 fix_taxon_names = FALSE,
                 simulateNCores=coresNow,
                 simulateParallel=TRUE,
                 simulate_ranges = "SkipUpper",
                 
                 fit_params = c("amplitude", "meanlog","sdlog" ),
                 data_file = "data/Leps_combo_final/moths_combined_dedup.csv")
print("=== Completed run_animal_model: lepidoptera ===")

# --- Free memory between organisms ---
rm(list = ls())          # wipes the global env, incl. run_animal_model and coresNow
gc(full = TRUE)          # return the memory to the OS
source("AmAni.manager.R")  # reload the pipeline functions (rm removed them)
coresNow <- 55           # re-set the core count (rm removed it too)


# --- Butterflies (GBIF) ---
check_status("butterflies")
run_animal_model("butterflies", clean = FALSE, redo = FALSE,
                 
                 animalFilterType = "extraStrict",
                 completenessMethod ="Chao1",
                 min_completeness = 0.55,
                 min_occurrences=5,
                 min_richness=100,
                 
                 fill_plants = TRUE,
                 useSamplingExtension=TRUE,
                 gammaErrorType = "narrow",
                 fill_animals = FALSE,
                 fix_taxon_names = FALSE,
                 simulateNCores=coresNow,
                 simulate_ranges = "SkipUpper",
                 data_file = "data/butterflies/0059395-260226173443078_cleaned_taxaFixed.csv")
print("=== Completed run_animal_model: butterflies ===")

# --- Free memory between organisms ---
rm(list = ls())          # wipes the global env, incl. run_animal_model and coresNow
gc(full = TRUE)          # return the memory to the OS
source("AmAni.manager.R")  # reload the pipeline functions (rm removed them)
coresNow <- 55           # re-set the core count (rm removed it too)

check_status("butterfl-curated")
run_animal_model("butterfl-curated", clean = FALSE,
                 data_file = "data/butterfl-curated/combined_butterfl-curated_taxaFixed.csv",
                 # data_file = "data/butterfl-curated/combinedCuratedOccur_cleaned.csv",
                 # secondary_data_file = "data/butterflies/0078048-260226173443078_cleaned.csv",
                 
                 animalFilterType = "extraStrict",
                 completenessMethod ="Chao1",
                 min_completeness = 0.55,
                 min_occurrences=5,
                 min_richness=200,
                 max_richness=2000,
                 
                 fill_plants = TRUE,
                 useSamplingExtension=TRUE,
                 gammaErrorType = "narrow",
                 fix_taxon_names = FALSE,
                 simulateParallel=TRUE,
                 simulateNCores=coresNow,
                 simulate_ranges = "SkipUpper",
                 fit_params = c("amplitude", "meanlog","sdlog" )
)
print("=== Completed run_animal_model: butterfl-curated ===")


# --- Free memory between organisms ---
rm(list = ls())          # wipes the global env, incl. run_animal_model and coresNow
gc(full = TRUE)          # return the memory to the OS
source("AmAni.manager.R")  # reload the pipeline functions (rm removed them)
coresNow <- 55           # re-set the core count (rm removed it too)


# --- Beetles (GBIF) ---
check_status("beetles")
run_animal_model("beetles", clean = FALSE, redo = FALSE,
                 animalFilterType = "extraStrict",
                 completenessMethod ="Chao1",
                 min_completeness = 0.35,
                 min_occurrences=5,
                 min_richness=100,
                 
                 fill_plants = TRUE,
                 useSamplingExtension=TRUE,
                 gammaErrorType = "narrow",
                 fill_animals = FALSE,
                 fix_taxon_names = FALSE,
                 simulateNCores=coresNow,
                 simulate_ranges = "SkipUpper",
                 disable_outlier_filter = TRUE, #I am going to shut this off for the diverse and under observed groups
                 fit_params = c("amplitude", "meanlog","sdlog" ),
                 
                 data_file = "data/beetle data/0011381-250402121839773_cleaned_taxaFixed.csv")
print("=== Completed run_animal_model: beetles ===")

# --- Free memory between organisms ---
rm(list = ls())          # wipes the global env, incl. run_animal_model and coresNow
gc(full = TRUE)          # return the memory to the OS
source("AmAni.manager.R")  # reload the pipeline functions (rm removed them)
coresNow <- 55           # re-set the core count (rm removed it too)

# --- DungBeetles (GBIF) ---
check_status("DungBeetles")
run_animal_model("DungBeetles", clean = FALSE, redo = FALSE,
                 animalFilterType = "extraStrict",
                 completenessMethod ="Chao1",
                 min_completeness = 0.65,
                 min_occurrences=5,
                 min_richness=30,
                 fill_plants = TRUE,
                 useSamplingExtension=TRUE,
                 gammaErrorType = "narrow",
                 fill_animals = FALSE,
                 fix_taxon_names = FALSE,
                 simulateNCores=coresNow,
                 simulate_ranges = "SkipUpper",
                 fit_params = c("amplitude", "meanlog","sdlog" ),
                 data_file = "data/dungb/0065336-260226173443078_cleaned_taxaFixed.csv")
print("=== Completed run_animal_model: DungBeetles ===")


# --- Hymenoptera (GBIF) ---
check_status("hymenoptera")
run_animal_model("hymenoptera", clean = FALSE, redo = FALSE,
                 animalFilterType = "extraStrict",
                 completenessMethod ="Chao1, ACE",
                 min_completeness = 0.3,
                 min_occurrences=100,
                 min_richness=100,
                 fill_plants = TRUE,
                 useSamplingExtension=TRUE,
                 gammaErrorType = "narrow",
                 fill_animals = FALSE,
                 fix_taxon_names = FALSE,
                 simulateNCores=coresNow,
                 simulate_ranges = "SkipUpper",
                 disable_outlier_filter = TRUE, #I am going to shut this off for the diverse and under observed groups
                 fit_params = c("amplitude", "meanlog","sdlog" ),
                 data_file = "data/hymenoptera/0060837-260226173443078_cleaned_taxaFixed.csv")
print("=== Completed run_animal_model: hymenoptera ===")

# --- Free memory between organisms ---
rm(list = ls())          # wipes the global env, incl. run_animal_model and coresNow
gc(full = TRUE)          # return the memory to the OS
source("AmAni.manager.R")  # reload the pipeline functions (rm removed them)
coresNow <- 55           # re-set the core count (rm removed it too)

# --- Hymenoptera (GBIF-noBees) ---
check_status("other_hymenoptera")
run_animal_model("other_hymenoptera", clean = FALSE, redo = FALSE, 
                 animalFilterType = "extraStrict",
                 completenessMethod ="Chao1",
                 min_completeness = 0.15,
                 min_occurrences=30,
                 min_richness=60,
                 fill_plants = TRUE,
                 useSamplingExtension=TRUE,
                 gammaErrorType = "narrow",
                 fill_animals = FALSE,
                 fix_taxon_names = FALSE,
                 simulateNCores=coresNow,
                 simulate_ranges = "SkipUpper",
                 disable_outlier_filter = TRUE, #I am going to shut this off for the diverse and under observed groups
                 fit_params = c("amplitude", "meanlog","sdlog" ),
                 data_file = "data/hymenoptera/0060837-260226173443078_noAntsNoBees_cleaned_taxaFixed.csv")
print("=== Completed run_animal_model: hymenoptera ===")

# --- Free memory between organisms ---
rm(list = ls())          # wipes the global env, incl. run_animal_model and coresNow
gc(full = TRUE)          # return the memory to the OS
source("AmAni.manager.R")  # reload the pipeline functions (rm removed them)
coresNow <- 55           # re-set the core count (rm removed it too)


# --- Bees (GBIF) ---
check_status("bees")
run_animal_model("bees", clean = FALSE, redo = FALSE,
                 fill_plants = TRUE,
                 animalFilterType = "extraStrict",
                 completenessMethod ="Chao1",
                 min_completeness = 0.65,
                 min_occurrences=5,
                 min_richness=30,
                 
                 useSamplingExtension=TRUE,
                 gammaErrorType = "narrow",
                 fix_taxon_names = FALSE,
                 simulateNCores=coresNow,
                 simulateParallel=TRUE,
                 disable_outlier_filter = FALSE,
                 simulate_ranges = "SkipUpper",
                 fit_params = c("amplitude", "meanlog","sdlog" ),
                 data_file = "data/bees/0059409-260226173443078_cleaned_taxaFixed.csv")
print("=== Completed run_animal_model: bees ===")



# --- Flies (GBIF, zip) ---
check_status("flies")
run_animal_model("flies", clean = FALSE, redo = FALSE,
                 animalFilterType = "balanced",
                 completenessMethod ="Chao1",
                 min_completeness = 0.5,
                 min_richness = 100,
                 
                 fill_plants = TRUE,
                 useSamplingExtension=TRUE,
                 gammaErrorType = "narrow",
                 fill_animals = FALSE,
                 fix_taxon_names = FALSE,
                 simulateNCores=coresNow,
                 disable_outlier_filter = TRUE, #I am going to shut this off for the diverse and under observed groups
                 fit_params = c("amplitude", "meanlog","sdlog" ),
                 data_file = "data/flies/0060830-260226173443078_cleaned_taxaFixed.csv")
print("=== Completed run_animal_model: flies ===")

# --- Free memory between organisms ---
rm(list = ls())          # wipes the global env, incl. run_animal_model and coresNow
gc(full = TRUE)          # return the memory to the OS
source("AmAni.manager.R")  # reload the pipeline functions (rm removed them)
coresNow <- 55           # re-set the core count (rm removed it too)

# --- termites (GBIF, zip) ---
check_status("termites")
run_animal_model("termites", clean = FALSE, redo = FALSE,
                 animalFilterType = "extraStrict",
                 completenessMethod ="Chao1",
                 min_completeness = 0.5,
                 min_occurrences = 20,
                 min_richness = 10,
                 fill_plants = TRUE,
                 useSamplingExtension=TRUE,
                 gammaErrorType = "narrow",
                 fill_animals = FALSE,
                 fix_taxon_names = FALSE,
                 simulateNCores=coresNow,
                 disable_outlier_filter = TRUE, #I am going to shut this off for the diverse and under observed groups
                 fit_params = c("amplitude", "meanlog","sdlog" ),
                 data_file = "data/termites/0008738-260921141020460_cleaned_taxaFixed.csv")
print("=== Completed run_animal_model: termites ===")

# --- MANTODEA (GBIF, zip) ---
check_status("Mantodea")
run_animal_model("Mantodea", clean = FALSE, redo = FALSE,
                 animalFilterType = "strict",
                 completenessMethod ="Chao1",
                 min_completeness = 0.25,
                 min_occurrences = 20,
                 min_richness = 5,
                 fill_plants = TRUE,
                 useSamplingExtension=TRUE,
                 gammaErrorType = "narrow",
                 fill_animals = FALSE,
                 fix_taxon_names = TRUE,
                 simulateNCores=coresNow,
                 disable_outlier_filter = TRUE, #I am going to shut this off for the diverse and under observed groups
                 fit_params = c("amplitude", "meanlog","sdlog" ),
                 data_file = "data/Mantodea/0006411-260928105237408_cleaned_taxaFixed.csv")
print("=== Completed run_animal_model: mantodea ===")


# --- MIRIDAE (GBIF, zip) --- <<This data is too sparse to be meaningful. I am skipping it.
# check_status("Miridae")
# run_animal_model("Miridae", clean = TRUE, redo = FALSE,
#                  animalFilterType = "balanced",
#                  completenessMethod ="Chao1",
#                  min_completeness = 0.4,
#                  fill_plants = TRUE,
#                  useSamplingExtension=TRUE,
#                  gammaErrorType = "narrow",
#                  fill_animals = FALSE,
#                  fix_taxon_names = TRUE,
#                  simulateNCores=coresNow,
#                  disable_outlier_filter = TRUE, #I am going to shut this off for the diverse and under observed groups
#                  fit_params = c("amplitude", "meanlog","sdlog" ),
#                  data_file = "data/Miridae/0008207-260928105237408.zip")
# print("=== Completed run_animal_model: miridae ===")


# # --- Ants curated (GBIF, zip) ---
# check_status("ants_curated")
# run_animal_model("ants_curated", clean = FALSE, redo = FALSE,
#                  animalFilterType = "extraStrict",
#                  completenessMethod ="Chao1",
#                  min_completeness = 0.6,
#                  min_occurrences = 50,
#                  min_richness = 80,
#                  fill_plants = TRUE,
#                  useSamplingExtension=TRUE,
#                  gammaErrorType = "narrow",
#                  fill_animals = FALSE,
#                  fix_taxon_names = FALSE,
#                  simulateNCores=coresNow,
#                  disable_outlier_filter = TRUE, #I am going to shut this off for the diverse and under observed groups
#                  fit_params = c("amplitude", "meanlog","sdlog" ),
#                  data_file = "data/ants_curated/dbFG_JeromeOrivel.csv")
# print("=== Completed run_animal_model: ants ===")

# --- Ants combined (GBIF, zip) ---
check_status("ants_combined")
run_animal_model("ants_combined", clean = FALSE, redo = FALSE,
                 animalFilterType = "strict",
                 completenessMethod ="Chao1",
                 min_completeness = 0.7,
                 min_richness = 80,
                 useSamplingExtension=TRUE,
                 gammaErrorType = "narrow",
                 fill_animals = FALSE,
                 fix_taxon_names = FALSE,
                 simulateNCores=coresNow,
                 disable_outlier_filter = TRUE, #I am going to shut this off for the diverse and under observed groups
                 fit_params = c("amplitude", "meanlog","sdlog" ),
                 data_file = "data/ants_curated/ants_combined.csv")
print("=== Completed run_animal_model: ants ===")


# --- bees combined (GBIF, zip) ---
check_status("bees_combined")
run_animal_model("bees_combined", clean = FALSE, redo = FALSE,
                 animalFilterType = "extraStrict",
                 completenessMethod ="Chao1",
                 min_completeness = 0.7,
                 min_occurrences = 50,
                 min_richness = 30,
                 fill_plants = TRUE,
                 useSamplingExtension=TRUE,
                 gammaErrorType = "narrow",
                 fill_animals = FALSE,
                 fix_taxon_names = FALSE,
                 simulateNCores=coresNow,
                 disable_outlier_filter = TRUE, #I am going to shut this off for the diverse and under observed groups
                 fit_params = c("amplitude", "meanlog","sdlog" ),
                 data_file = "data/bees_curated/bees_combined.csv")
print("=== Completed run_animal_model: ants ===")

# --- Free memory between organisms ---
rm(list = ls())          # wipes the global env, incl. run_animal_model and coresNow
gc(full = TRUE)          # return the memory to the OS
source("AmAni.manager.R")  # reload the pipeline functions (rm removed them)
coresNow <- 55           # re-set the core count (rm removed it too)

# --- Mammals (GBIF, zip) ---
check_status("mammals")
run_animal_model("mammals", clean = FALSE, redo = FALSE,
                 animalFilterType = "extraStrict",
                 completenessMethod ="Chao1",
                 min_completeness = 0.8,
                 min_occurrences = 100,
                 
                 fill_plants = TRUE,
                 useSamplingExtension=TRUE,
                 gammaErrorType = "narrow",
                 fill_animals = FALSE,
                 fix_taxon_names = FALSE,
                 simulateNCores=coresNow,
                 disable_outlier_filter=FALSE,
                 simulate_ranges = "SkipUpper",
                 fit_params = c("amplitude", "meanlog","sdlog" ),
                 data_file = "data/mammals/0063424-260226173443078_cleaned_taxaFixed.csv")
print("=== Completed run_animal_model: mammals ===")

# --- Frogs (GBIF) ---
check_status("frogs")
run_animal_model("frogs", clean = FALSE, redo = FALSE,
                 animalFilterType = "extraStrict",
                 completenessMethod ="Chao1",
                 min_completeness = 0.85,
                 min_occurrences = 100,
                 
                 fill_plants = TRUE,
                 useSamplingExtension=TRUE,
                 gammaErrorType = "narrow",
                 fix_taxon_names = FALSE,
                 simulateNCores=coresNow,
                 simulateParallel=TRUE,
                 disable_outlier_filter=FALSE,
                 simulate_ranges = "SkipUpper",
                 fit_params = c("amplitude", "meanlog","sdlog" ),
                 data_file = "data/frogs/0059434-260226173443078_cleaned_taxaFixed.csv")

print("=== Completed run_animal_model: frogs ===")

# --- Tyrannidae (GBIF, zip) --- ##flycatchers
check_status("tyrannidae")
run_animal_model("tyrannidae", clean = FALSE, redo = FALSE,
                 animalFilterType = "extraStrict",
                 completenessMethod ="Chao1",
                 min_completeness = 0.85,
                 min_occurrences = 50,
                 
                 fill_plants = TRUE,
                 useSamplingExtension=TRUE,
                 gammaErrorType = "narrow",
                 fix_taxon_names = FALSE,
                 simulateNCores=coresNow,
                 simulateParallel=TRUE,
                 simulate_ranges = "SkipUpper",
                 disable_outlier_filter = FALSE,
                 fit_params = c("amplitude", "meanlog","sdlog" ),
                 data_file = "data/tyrannidae/0063441-260226173443078_cleaned_taxaFixed.csv")


print("=== Completed run_animal_model: tyrannidae ===")

# --- Antbirds (GBIF) ---
check_status("antbirds") 
###NOTE: The below parameters result in a 0.98 gamma expansion ratio, so the
###model almost perfectly predicts the actual ant bird richness. however, there
###are 199 cells in the model so i want to see what happens when we increase the
###completeness of the cells that contribute to the model.

# run_animal_model("antbirds", clean = TRUE, redo = FALSE,
#                  min_occurrences = "auto", min_completeness = "auto",
#                  min_richness = 1,
#                  filter_plant_completeness = 0.60,
#                  data_file = "data/antbirds/0000507-250424135625364.csv")

###NOTE: These parameters resulted in much higher estimates of ant-bird diversity
# run_animal_model("antbirds", clean = TRUE, redo = FALSE,
#                  min_occurrences = 300, min_completeness = 0.9,
#                  min_richness = 1,
#                  filter_plant_completeness = 0.80,
#                  data_file = "data/antbirds/0000507-250424135625364.csv")


run_animal_model("antbirds", clean = FALSE, redo = FALSE,
                 animalFilterType = "extraStrict",
                 completenessMethod ="Chao1",
                 min_completeness = 0.85,
                 min_occurrences = 50,
                 
                 fill_plants = TRUE,
                 useSamplingExtension=TRUE,
                 gammaErrorType = "narrow",
                 fill_animals = FALSE,
                 fix_taxon_names = FALSE,
                 simulateNCores=coresNow,
                 simulate_ranges = "SkipUpper",
                 fit_params = c("amplitude", "meanlog","sdlog" ),
                 data_file = "data/antbirds/0000507-250424135625364_cleaned_taxaFixed.csv")
print("=== Completed run_animal_model: antbirds ===")



# --- Finches (GBIF) ---
check_status("finches")
run_animal_model("finches", clean = FALSE, redo = FALSE,
                 animalFilterType = "extraStrict",
                 completenessMethod ="Chao1",
                 min_completeness = 0.6,
                 min_occurrences = 15,
                 min_richness = 5,
                 fill_plants = TRUE,
                 useSamplingExtension=TRUE,
                 gammaErrorType = "narrow",
                 fill_animals = FALSE,
                 fix_taxon_names = FALSE,
                 disable_outlier_filter=FALSE,
                 simulate_ranges = "SkipUpper",
                 fit_params = c("amplitude", "meanlog","sdlog" ),
                 data_file = "data/finches/0011407-250402121839773_cleaned_taxaFixed.csv")
print("=== Completed run_animal_model: finches ===")



# --- Squamates (squamates) ---
check_status("squamata")
run_animal_model("squamata", clean = FALSE, redo = FALSE,
                 animalFilterType = "extraStrict",
                 completenessMethod ="Chao1",
                 min_completeness = 0.85,
                 min_occurrences = 100,
                 
                 fill_plants = TRUE,
                 useSamplingExtension=TRUE,
                 gammaErrorType = "narrow",
                 fill_animals = FALSE,
                 fix_taxon_names = FALSE,
                 simulate_ranges = "SkipUpper",
                 fit_params = c("amplitude", "meanlog","sdlog" ),
                 data_file = "data/squamates/0063428-260226173443078_cleaned_taxaFixed.csv")
print("=== Completed run_animal_model: squamata ===")



# --- Free memory between organisms ---
rm(list = ls())          # wipes the global env, incl. run_animal_model and coresNow
gc(full = TRUE)          # return the memory to the OS
source("AmAni.manager.R")  # reload the pipeline functions (rm removed them)
coresNow <- 55           # re-set the core count (rm removed it too)


###run_extinction_model for each organism####

run_extinction_model("cockroaches")
run_extinction_model("butterfl-curated")
run_extinction_model("odonata")
run_extinction_model("odonata-combined")
run_extinction_model("odonata-curated")
run_extinction_model("hemiptera")

# --- Free memory between organisms ---
rm(list = ls())          # wipes the global env, incl. run_animal_model and coresNow
gc(full = TRUE)          # return the memory to the OS
source("AmAni.manager.R")  # reload the pipeline functions (rm removed them)
coresNow <- 55           # re-set the core count (rm removed it too)

#run_extinction_model("lepidoptera")
#run_extinction_model("butterflies")
run_extinction_model("beetles")
run_extinction_model("DungBeetles")
run_extinction_model("other_hymenoptera")
run_extinction_model("bees")
run_extinction_model("moths")

run_extinction_model("termites")
run_extinction_model("Mantodea")
#run_extinction_model("Miridae")
run_extinction_model("ants_curated")
run_extinction_model("ants_combined")
run_extinction_model("bees_combined")

# --- Free memory between organisms ---
rm(list = ls())          # wipes the global env, incl. run_animal_model and coresNow
gc(full = TRUE)          # return the memory to the OS
source("AmAni.manager.R")  # reload the pipeline functions (rm removed them)
coresNow <- 55           # re-set the core count (rm removed it too)

run_extinction_model("frogs")
run_extinction_model("tyrannidae")
run_extinction_model("antbirds")
run_extinction_model("finches")
run_extinction_model("mammals")

# --- Free memory between organisms ---
rm(list = ls())          # wipes the global env, incl. run_animal_model and coresNow
gc(full = TRUE)          # return the memory to the OS
source("AmAni.manager.R")  # reload the pipeline functions (rm removed them)
coresNow <- 55           # re-set the core count (rm removed it too)

run_extinction_model("squamata")
run_extinction_model("flies")

gg wp# --- Free memory between organisms ---
rm(list = ls())          # wipes the global env, incl. run_animal_model and coresNow
gc(full = TRUE)          # return the memory to the OS
source("AmAni.manager.R")  # reload the pipeline functions (rm removed them)
coresNow <- 55           # re-set the core count (rm removed it too)

# =============================================================================
# IUCN threat classification (based on range reduction from extinction model)
# =============================================================================
source("4-ExtModeller/scripts/utils/iucn_threat_classification.R")

extinction_model_organisms <- c(
  "butterfl-curated", "odonata", "odonata-combined", "odonata-curated",
  "hemiptera", "moths", "butterflies", "beetles", "DungBeetles",
  "hymenoptera", "bees", "flies", "frogs", "tyrannidae", "antbirds",
  "finches", "mammals", "cockroaches", "squamata"
)

for (organism in extinction_model_organisms) {
  dirs <- organism_dirs(organism)
  thresholds_path <- file.path(dirs$extinction_model, "summary_range_reduction_thresholds.csv")

  if (!file.exists(thresholds_path)) {
    message("Skipping IUCN classification for '", organism, "' — file not found: ", thresholds_path)
    next
  }

  message("Classifying IUCN threat categories for: ", organism)
  iucn_result <- compute_iucn_category_counts(dirs$extinction_model)
  plot_iucn_threat(
    iucn_result,
    output_dir = dirs$extinction_model,
    filename_prefix = paste0(organism, "_iucn_threat_classification"),
    organism_label = organism
  )
}



# --- Plant verification ---
# Runs the full pipeline on plant data as a sanity check.
# Gamma richness should approximate the known plant richness.
# run_plant_model()
# 
#run_extinction_model("plants")

# =============================================================================
# Shape-shifted re-runs (fit_params examples)
# =============================================================================
# Re-run specific organisms with shape-shifting to compare against the
# default amplitude-only scaling. These estimate meanlog and/or sdlog from
# the animal MDD data and warp the plant range-size distribution accordingly.
#
# Amplitude + location shift (range sizes shifted to match animal):
# run_animal_model("butterflies", clean = TRUE, redo = FALSE,
#                  min_occurrences = "auto", min_completeness = "auto",
#                  min_richness = 800, disable_outlier_filter = TRUE,
#                  filter_plant_completeness = 0.60,
#                  data_file = "data/butterflies/0059395-260226173443078.csv",
#                  fit_params = c("amplitude", "meanlog"))
#
# All three parameters (location + spread + amplitude):
# run_animal_model("beetles", clean = TRUE, redo = FALSE,
#                  min_occurrences = "auto", min_completeness = "auto",
#                  min_richness = 25, filter_plant_completeness = 0.60,
#                  data_file = "data/beetle data/0011381-250402121839773.csv",
#                  fit_params = c("amplitude", "meanlog", "sdlog"))


# =============================================================================
# Summaries
# =============================================================================
source("AmAni.summarize.R")

# Full cross-organism summary
summarize_all()


# Group summaries (insects only, vertebrates only, etc.)
summarize_all(results_name = "insects", group = "insects")
summarize_all(results_name = "vertebrates",group = "vertebrates")
summarize_all(results_name = "comparisons",group = "comparisons")
summarize_all(include = c("butterflies", "butterfl-curated", "odonata", "odonata-curated", "odonata-combined"))

# 
# hyperdiverse group summaries
summarize_all(include = c("butterflies", "butterfl-curated", "lepidoptera", "beetles", "DungBeetles", "hymenoptera", "flies", "hemiptera"))

# mesodiverse insect  summaries
summarize_all(include = c("odonata","odonata-curated", "odonata-combined",  "cockroaches", "allRoaches", "bees"))

# verebrate  summaries
summarize_all(include = c("squamata","frogs", "mammals", "tyrannidae", "finches", "antbirds"))

# Individual organisms
summarize_organism("cockroaches")

summarize_organism("butterfl-curated")
summarize_organism("odonata")
summarize_organism("odonata-combined")
summarize_organism("odonata-curated")

summarize_organism("termites")
summarize_organism("Mantodea")
#
#summarize_organism("Miridae")
summarize_organism("ants_curated")
summarize_organism("ants_combined")
summarize_organism("bees_combined")

summarize_organism("hemiptera")
summarize_organism("lepidoptera")
summarize_organism("moths")


summarize_organism("butterflies")
summarize_organism("beetles")
summarize_organism("DungBeetles")
#summarize_organism("hymenoptera")
summarize_organism("other_hymenoptera")

#summarize_organism("bees")
summarize_organism("flies")
summarize_organism("frogs")
summarize_organism("tyrannidae")
summarize_organism("antbirds")
summarize_organism("finches")
summarize_organism("mammals")
summarize_organism("squamata")
# Publication figures (Figures 1-3): PDF + PNG in output/publication_figures/
make_publication_figures()
# make_publication_figures(figures = "fig3")   # just one figure

# Spatial extent overlap comparison (5 smallest animal extents)
 plot_spatial_extent_overlap()
# plot_spatial_extent_overlap(group = "insects")

# Workflow compilation diagram
# compile_workflow_diagram("butterflies")
# compile_workflow_diagram("cockroaches")


