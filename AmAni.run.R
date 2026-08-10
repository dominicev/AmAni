#Next steps:
#0. make a simulated dataset of fake animals where you know the true gamma diversity
#this will allow us to know what parameters of the animal model result in the most accurate estimates.
#in particular, i think it might not just be important that the overlapping cells have
#high completeness in both the plant and animals, but maybe it matters which cells are overlapping 
#(i.e., if the datasets overlap in areas where the plant diversity is very low, or unusually high)
#1. for those results that have a lot of overlapping grid cells with plants (we 
#know this N from looking at the richness and endemism correlation plots), increase
#the minimum completeness threshold for the plant and/or the animal until the N
#is no lower than ~10
#2. add curated datasets for as many insects as possible
#
#
#
#


# AmAni.run.R — Run the full pipeline for each organism
source("AmAni.manager.R")
source("AmAni.summarize.R")


# --- Cockroaches (Madre de Dios) ---
# Curated dataset — no cell filtering needed (min_occurrences = 0 disables it)
# fit_params controls which log-normal parameters are estimated from the animal
# MDD data when scaling the plant spatial extent distribution:
#   c("amplitude")                      — default (current behaviour, amplitude only)
#   c("amplitude", "meanlog")           — also shift location (range sizes)
#   c("amplitude", "meanlog", "sdlog")  — fit all three (location + spread)
check_status("cockroaches")
run_animal_model("cockroaches", clean = FALSE, min_occurrences = 0,
                 #fill_plants = TRUE,
                 fill_animals = TRUE,
                 #alphaRichnessMethod = "extend",
                 redo = FALSE,
                 data_file = "data/DarwinCoreAnimalOccurences.8-1-2026.xlsx")
run_extinction_model("cockroaches")


# --- Cockroaches, all South America (GBIF) ---
# Very sparse GBIF data (267 spp, ~5k records). Auto filter likely triggers
# the safety guard (< 30 cells), so filtering is effectively skipped.
check_status("allRoaches")
run_animal_model("allRoaches", clean = TRUE, redo = FALSE,
                 min_occurrences = "auto", min_completeness = "auto",
                 min_richness = 15,
                 filter_plant_completeness = 0.60,
                 redo = FALSE,
                 data_file = "data/cockroaches all s america/0011428-250402121839773.csv")
run_extinction_model("allRoaches")

# --- Butterflies (Curated and GBIF combined) ---
# Curated dataset — no cell filtering needed (min_occurrences = 0 disables it)
# fit_params controls which log-normal parameters are estimated from the animal
#  c("amplitude")                      — default (current behaviour, amplitude only)
#   c("amplitude", "meanlog")           — also shift location (range sizes)
#   c("amplitude", "meanlog", "sdlog")  — fit all three (location + spread)
check_status("butterfl-curated")
run_animal_model("butterfl-curated", clean = FALSE,
                 data_file = "data/butterfl-curated/combinedCuratedOccur_cleaned.csv",
                 secondary_data_file = "data/butterflies/0078048-260226173443078_cleaned.csv",
                 fill_plants = TRUE,
                 fill_animals = TRUE,
                  min_occurrences = 600,
                  min_completeness = 0.2,
                  min_richness = 300,
                  max_richness = 3100,
                 redo = FALSE,
                 fit_params = c("amplitude", "meanlog","sdlog" ),
                 )
run_extinction_model("butterfl-curated")




# --- Odonata (GBIF) ---
check_status("odonata")
run_animal_model("odonata", clean = TRUE, redo = FALSE,
                 min_occurrences = "auto", min_completeness = "auto",
                 min_richness = 15,
                 filter_plant_completeness = 0.60,

                 fit_params = c("amplitude", "meanlog","sdlog" ),
                 data_file = "data/odonata/0059398-260226173443078.csv")
run_extinction_model("odonata")




# --- Hemiptera (GBIF) ---
check_status("hemiptera")
run_animal_model("hemiptera", clean = TRUE, redo = FALSE,
                 min_occurrences = "auto", min_completeness = "auto",
                 min_richness = 20,
                 filter_plant_completeness = 0.60,
                 fit_params = c("amplitude", "meanlog","sdlog" ),
                 data_file = "data/hemiptera/0059425-260226173443078.csv")
run_extinction_model("hemiptera")





# --- Lepidoptera (GBIF) ---
check_status("lepidoptera")
run_animal_model("lepidoptera", clean = TRUE, redo = FALSE,
                 min_occurrences = "auto", min_completeness = "auto",
                 min_richness = 800,
                 filter_plant_completeness = 0.70,
                 fit_params = c("amplitude", "meanlog","sdlog" ),
                 data_file = "data/lepidoptera/0060839-260226173443078.csv")
run_extinction_model("lepidoptera")

# --- Butterflies (GBIF) ---
check_status("butterflies")
run_animal_model("butterflies", clean = TRUE, redo = FALSE,
                 min_occurrences = "auto", min_completeness = 0.7,
                 min_richness = 800,
                 disable_outlier_filter = TRUE,
                 filter_plant_completeness = 0.7,
                 fit_params = c("amplitude", "meanlog","sdlog" ),
                 data_file = "data/butterflies/0078048-260226173443078.zip")
run_extinction_model("butterflies")

# --- Beetles (GBIF) ---
check_status("beetles")
run_animal_model("beetles", clean = TRUE, redo = FALSE,
                 min_occurrences = "auto", min_completeness = "auto",
                 min_richness = 25,
                 filter_plant_completeness = 0.70,
                 fit_params = c("amplitude", "meanlog","sdlog" ),
                 data_file = "data/beetle data/0011381-250402121839773.csv")
run_extinction_model("beetles")


# --- DungBeetles (GBIF) ---
check_status("DungBeetles")
run_animal_model("DungBeetles", clean = TRUE, redo = FALSE,
                 min_occurrences = "auto", min_completeness = "auto",
                 min_richness = 25,
                 filter_plant_completeness = 0.70,
                 fit_params = c("amplitude", "meanlog","sdlog" ),
                 data_file = "data/dungb/0065336-260226173443078_cleaned.csv")
run_extinction_model("DungBeetles")


# --- Hymenoptera (GBIF) ---
check_status("hymenoptera")
run_animal_model("hymenoptera", clean = TRUE, redo = FALSE,
                 min_occurrences = "auto", min_completeness = "auto",
                 min_richness = 20,
                 filter_plant_completeness = 0.70,
                 fit_params = c("amplitude", "meanlog","sdlog" ),
                 data_file = "data/hymenoptera/0060837-260226173443078.csv")
run_extinction_model("hymenoptera")

# --- Bees (GBIF) ---
check_status("bees")
run_animal_model("bees", clean = TRUE, redo = FALSE,
                 min_occurrences = "auto", min_completeness = "auto",
                 min_richness = 15,
                 filter_plant_completeness = 0.60,
                 data_file = "data/bees/0059409-260226173443078.csv")
run_extinction_model("bees")

# --- Flies (GBIF, zip) ---
check_status("flies")
run_animal_model("flies", clean = TRUE, redo = FALSE,
                 min_occurrences = "auto", min_completeness = "auto",
                 min_richness = 20,
                 filter_plant_completeness = 0.60,
                 data_file = "data/flies/0060830-260226173443078.zip")
run_extinction_model("flies")

# --- Mammals (GBIF, zip) ---
check_status("mammals")
run_animal_model("mammals", clean = TRUE, redo = FALSE,
                 min_occurrences = "auto", min_completeness = "auto",
                 min_richness = 10,
                 filter_plant_completeness = 0.60,
                 fit_params = c("amplitude", "meanlog","sdlog" ),
                 data_file = "data/mammals/0063424-260226173443078.zip")
run_extinction_model("mammals")

# --- Frogs (GBIF) ---
check_status("frogs")
run_animal_model("frogs", clean = TRUE, redo = FALSE,
                 min_occurrences = "auto", min_completeness = "auto",
                 min_richness = 5,
                 filter_plant_completeness = 0.60,
                 fit_params = c("amplitude", "meanlog","sdlog" ),
                 data_file = "data/frogs/0059434-260226173443078.csv")
run_extinction_model("frogs")


# --- Tyrannidae (GBIF, zip) --- ##flycatchers
check_status("tyrannidae")
run_animal_model("tyrannidae", clean = TRUE, redo = FALSE,
                 min_occurrences = "auto", min_completeness = "auto",
                 min_richness = 5,
                 filter_plant_completeness = 0.60,
                 fit_params = c("amplitude", "meanlog","sdlog" ),
                 data_file = "data/tyrannidae/0063441-260226173443078.zip")
run_extinction_model("tyrannidae")

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


run_animal_model("antbirds", clean = TRUE, redo = FALSE,
                 min_occurrences = 300, min_completeness = 0.9,
                 min_richness = 1,
                 filter_plant_completeness = 0.80,
                 fit_params = c("amplitude", "meanlog","sdlog" ),
                 data_file = "data/antbirds/0000507-250424135625364.csv")

run_extinction_model("antbirds")

# --- Finches (GBIF) ---
check_status("finches")
run_animal_model("finches", clean = TRUE, redo = FALSE,
                 min_occurrences = 100, min_completeness = "auto",
                 min_richness = 1,
                 filter_plant_completeness = 0.50,
                 data_file = "data/finches/0011407-250402121839773.csv")
run_extinction_model("finches")


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

# Full cross-organism summary
summarize_all()


# Group summaries (insects only, vertebrates only, etc.)
# summarize_all(group = "insects")
# summarize_all(group = "vertebrates")
# 
# hyperdiverse group summaries
summarize_all(include = c("butterflies", "lepidoptera", "beetles", "hymenoptera", "flies", "hemiptera"))

# mesodiverse insect  summaries
summarize_all(include = c("odonata", "cockroaches", "allRoaches", "bees"))

# verebrate  summaries
summarize_all(include = c("frogs", "mammals", "tyrannidae", "finches", "antbirds"))

# Individual organisms
summarize_organism("butterfl-curated")
summarize_organism("butterflies")
summarize_organism("cockroaches")
summarize_organism("odonata")
summarize_organism("flies")
summarize_organism("mammals")
summarize_organism("frogs")
summarize_organism("tyrannidae")
summarize_organism("beetles")

# Spatial extent overlap comparison (5 smallest animal extents)
# plot_spatial_extent_overlap()
# plot_spatial_extent_overlap(group = "insects")

# Workflow compilation diagram
# compile_workflow_diagram("butterflies")
# compile_workflow_diagram("cockroaches")


