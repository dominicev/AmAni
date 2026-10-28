#This notebook will take the simulated sampling from the simulated truth datasets and run the animal model for them using different run_animal_model parameters
#
#
#####NOTE:::: I think I need to do additional replicates and vary min_occurrences, while keeping min_completeness at 0 (the CHAO1 completeness estimation will fail on the simulated data because I did not simulate abundances)
library(sf)

source("AmAni.summarize.R")
library(dplyr)
library(purrr)


getwd()
# --- import 1 dataset to confirm it can run on the format ---
animal.dir<-"./2-ModellingAnimals/data"
data.dir<-"./5-SensitivityTesting/output/sampleDatasets/"

files_full <- list.files(data.dir, pattern = "\\.RData$", full.names = TRUE)
files_short <- list.files(data.dir, pattern = "\\.RData$", full.names = FALSE)

files_full_sub <- files_full[grepl("like|variedSiteEvenness", files_full, ignore.case = TRUE)]
files_short_sub <- files_short[grepl("like|variedSiteEvenness", files_short, ignore.case = TRUE)]
length(files_full_sub)


fn<-1

settingCombinations1 <- tidyr::expand_grid(
  # Fixed parameters (same across all runs)
  # clean                  = FALSE,
  # force_clean          = FALSE,
  # redo                 = TRUE,
  # min_occurrences      = "auto",
  # min_completeness     = "auto",
  # max_richness         = NULL,
  # disable_outlier_filter = TRUE,
  # estimateBounds       = FALSE,
  # skipRawEndemism      = TRUE,

  # Varying parameters (one row per unique combination)
  fill_plants               = c(FALSE, TRUE),
  fill_animals              = c(FALSE, TRUE),
  filter_plant_completeness = NULL,
  fit_params                = list(
    # c("amplitude"),
    # c("amplitude", "meanlog"),
    # c("amplitude", "sdlog"),
    c("amplitude", "meanlog", "sdlog")
  )
)


settingCombinations2 <- tidyr::expand_grid(
  # Fixed parameters (same across all runs)
  # clean                  = FALSE,
  # force_clean          = FALSE,
  # redo                 = TRUE,
  # min_occurrences      = "auto",
  # min_completeness     = "auto",
  # max_richness         = NULL,
  # disable_outlier_filter = TRUE,
  # estimateBounds       = FALSE,
  # skipRawEndemism      = TRUE,
  
  # Varying parameters (one row per unique combination)
  fill_plants               = FALSE,
  fill_animals              = FALSE,
  filter_plant_completeness = list(NULL, 0.5, 0.75, 0.9),
  fit_params                = list(
    # c("amplitude"),
    # c("amplitude", "meanlog"),
    # c("amplitude", "sdlog"),
    c("amplitude", "meanlog", "sdlog")
  )
)

settingCombinations3 <- tidyr::expand_grid(
  # Fixed parameters (same across all runs)
  # clean                  = FALSE,
  # force_clean          = FALSE,
  # redo                 = TRUE,
  # min_occurrences      = "auto",
  # min_completeness     = "auto",
  # max_richness         = NULL,
  # disable_outlier_filter = TRUE,
  # estimateBounds       = FALSE,
  # skipRawEndemism      = TRUE,
  
  # Varying parameters (one row per unique combination)
  fill_plants               = FALSE,
  fill_animals              = FALSE,
  filter_plant_completeness = NULL,
  fit_params                = list(
     c("amplitude"),
     c("amplitude", "meanlog"),
     c("amplitude", "sdlog"),
    c("amplitude", "meanlog", "sdlog")
  )
)

settingCombinations<-bind_rows(list(settingCombinations3, settingCombinations2, settingCombinations1))

fn=5
set=3
source("AmAni.manager.R")
files_full_sub[12:15*2]


# --- Re-run ONLY the 10 incomplete base x set combinations ------------------
# Each row is one (input file, parameter set) pair that never produced a
# gamma_richnessfromSPRICH*.txt output. `fn` is resolved to the most recent
# sample file matching each base simulation; `set` indexes settingCombinations.
# (Change tail(hits, 1L) to hits[1L] below to use the earliest file instead.)
target_specs <- data.frame(
  pattern = c(
    "blatt_blatt_roachLike\\.sample",             # sim_blatt_roachLike
    "butter_butter_butterflyLike\\.sample",       # sim_butter_butterflyLike
    "butter_butter_butterflyLike2\\.sample",      # sim_butter_butterflyLike2
    "butter_butter_butterflyLike3\\.sample",      # sim_butter_butterflyLike3
    "butter_butter_butterflyLike4\\.sample",      # sim_butter_butterflyLike4
    "coleo_coleo_beetleLike\\.sample",            # sim_coleo_beetleLike
    "coleo_coleo_beetleLike2\\.sample",           # sim_coleo_beetleLike2
    "coleo_coleo_beetleLike3\\.sample",           # sim_coleo_beetleLike3
    "low_run1_low2_variedSiteEvenness\\.sample",  # sim_low_run1_variedSiteEvenness
    "veryLow_veryLow1_variedSiteEvenness\\.sample" # sim_veryLow_variedSiteEvenness
  ),
  set = c(11L, 12L, 12L, 12L, 12L, 12L, 12L, 12L, 11L, 5L),
  stringsAsFactors = FALSE
)
# Resolve each pattern to the most recent matching file index in files_full_sub
target_specs$fn <- vapply(target_specs$pattern, function(p) {
  hits <- grep(p, basename(files_full_sub))
  if (length(hits) == 0L) NA_integer_ else tail(hits, 1L)
}, integer(1))
# Fail loudly if any target could not be matched to an input file
stopifnot(!anyNA(target_specs$fn))

for (t in seq_len(nrow(target_specs))) {
  fn  <- target_specs$fn[t]
  set <- target_specs$set[t]
  tryCatch({
  
    loaded_objects <- load(files_full_sub[fn])
    fnElem<-unlist(lapply(  unlist(strsplit(basename(files_full_sub[fn]), split="\\.")), function(x){strsplit(x, "_")}))
    
    # alternatve code for fnElem 
    # m <- regmatches(basename(files_full_sub[fn]),
    #             regexpr("^(?P<term>[^.]+)\\.(?P<scen>[^_]+)_(?P<level>[^_]+)_(?P<rep>\\d+)\\.RData$",
    #                     basename(files_full_sub[fn]), perl = TRUE))
    
    
    
    # names(df)
    # unique(df$scenario )
    # head(df)
    # 
    occurrenceStatusVector<-rep("PRESENT",length(df$scientificName))
    df[["occurrenceStatus"]]<-occurrenceStatusVector
    
    # 1. Transform JUST the geometry column (sf knows how to handle an sfc column directly)
    geom_wgs84 <- st_transform(df$geometry, 4326)
    
    # 2. Turn off s2 temporarily for lightning-fast flat math
    sf_use_s2(FALSE)
    centroids <- st_coordinates(st_centroid(geom_wgs84))
    sf_use_s2(TRUE) # turn it back on
    
    # 3. Use data.table syntax (:=) to add the columns instantly in-place
    df[, `:=`(
      decimalLongitude = centroids[, "X"],
      decimalLatitude  = centroids[, "Y"]
    )]
    
    # 4. Subset for DarwinCore Export
    # We drop 'geometry' because CSVs cannot store spatial polygons
    final_export <- df[, c("scientificName","species", "decimalLatitude", "decimalLongitude","occurrenceStatus")]
    
    # 1. Build the folder name and target directory path
    folder_name <- paste0("sim_", unique(df$scenario), "_", 
                          if(length(fnElem)==7){ fnElem[3]}
                             else {fnElem[4]}
                          , "_",toString(set))
    target_dir  <- file.path(animal.dir, folder_name)
    
    # 2. Build the full file path where the CSV will live
    file_name   <- paste0(folder_name, "_Data.csv")
    full_path   <- file.path(target_dir, file_name)
    
    # 3. Create the folder if it doesn't exist yet
    # (recursive = TRUE ensures all parent folders are built if missing)
    if (!dir.exists(target_dir)) {
      dir.create(target_dir, recursive = TRUE)
    }
    
    
    
    # 4. Now export the CSV safely
    write.csv(final_export, full_path, row.names = FALSE)
    
#undebug(run_animal_model) ###NOTE: If there are only a handful of sites and the overlap with plant data happens to be very small, then the model will fail with the error:
#
#Error in numericDeriv(form[[3L]], names(ind), env, central = nDcentral) : 
#    Missing value or an infinity produced when evaluating the model
outdir<-paste0("./p_output/",folder_name, "/animal_model/")
existingOut<-list.files(outdir, pattern = "gamma_richnessfromSPRICH.*\\.txt$")
if(length(existingOut)==1){
  message(paste0(
    folder_name, " simulation completed already. Continuing..."
      ))

}else{
    #undebug(simulateScaledSpeciesRichness)
    #undebug(run_animal_model)
    run_animal_model(folder_name,
                     data_file = sub("^.*?data/", "data/", full_path),
                     clean = FALSE,
                     redo = FALSE,
                     min_occurrences = "auto",
                     min_completeness = "auto",
                     min_richness = 1,
                     estimateBounds=FALSE,
                     simulate_ranges = FALSE,
                     skipRawEndemism = TRUE,
                     customOutputFolderPrefix = "p_",
                     disable_outlier_filter = TRUE,
                     customOutputFolderSuffix = paste0("_",toString(set)),
                     fill_plants=settingCombinations$fill_plants[[set]],
                     fill_animals=settingCombinations$fill_animals[[set]],
                     fit_params = settingCombinations$fit_params[[set]],
                     filter_plant_completeness = settingCombinations$filter_plant_completeness[[set]]
                       )
}
  }, error = function(e) {
    message("Error at iteration ", fn, ": ", conditionMessage(e))
  })
}

#these functions will let me parameterize min_completeness in the loop below
compute_site_completeness <- function(occ_data,
                                      true_richness,
                                      site_col    = "grid_id",
                                      species_col = "species") {
  
  # Observed richness per site in the sample
  obs_richness <- aggregate(
    list(n_obs = occ_data[[species_col]]),
    by  = setNames(list(occ_data[[site_col]]), site_col),
    FUN = function(x) length(unique(x))
  )
  
  # Join true richness and take the ratio
  site_completeness <- merge(obs_richness, true_richness,
                             by = site_col, all.x = TRUE)
  site_completeness$completeness <- site_completeness$n_obs /
    site_completeness$n_true
  
  return(site_completeness)
}

compute_completeness_stats <- function(completeness_computed) {
  list(
    comp_mean   = mean(completeness_computed$completeness),
    comp_sd     = sd(completeness_computed$completeness),
    comp_spread = max(completeness_computed$completeness) - min(completeness_computed$completeness),
    comp_max = max(completeness_computed$completeness)
  )
}


# Minimum per-cell completeness floors to test.
# A cell is retained only if its TRUE completeness (observed distinct species /
# true species richness in that cell) meets or exceeds the floor. Because the
# simulated data have no meaningful abundances (every species has a single
# occurrence per cell), CHAO1-based completeness is unusable. Since true richness
# varies enormously across cells, a single min_occurrences cutoff cannot target a
# completeness level; instead we threshold each cell's completeness directly using
# ground-truth richness. Filtering happens upstream (cell subsetting), so the
# model's own min_occurrences / min_completeness filters stay disabled.
completeness_floors <- c(0.3, 0.5, 0.7, 0.9)

#load the true richness
true_env <- new.env()
true_names <- load(
  file.path("C:/Users/dominice/OneDrive - University of Illinois - Urbana",
            "Projects/Peruvian Blattodea/Cockroaches in a coal mine",
            "AmAni/5-SensitivityTesting/data/simFromPlant/animal_grid_data.RData"),
  envir = true_env
)
true_data <- as.data.frame(get(true_names[1], envir = true_env))

# True species richness per grid cell (distinct species present in truth)
true_richness <- aggregate(
  list(n_true = true_data$species),
  by  = list(grid_id = true_data$grid_id),
  FUN = function(x) length(unique(x))
)


# --- Select datasets for the completeness sweep --------------------------
# All sampled datasets EXCEPT the roachLike simulations and the "even"
# variedSiteEvenness simulations. Where a base simulation was sampled more than
# once (multiple timestamps in the filename), keep only the most recent file.
#
# NOTE (2026-08-20): Each input trial has TWO sampled versions (two timestamps),
# so 22 files qualify after exclusions. To keep runtime manageable we run only
# the 11 most-recent versions here. The other 11 (older timestamps) can be run
# later -- e.g. swap tail(...) for head(...) in the dedup below, or drop the
# dedup entirely to process all 22.
sens_base <- sub("\\.sample.*$", "", basename(files_full_sub))
sens_idx  <- which(!grepl("roachLike", basename(files_full_sub)) &
                     !grepl("^even_",   basename(files_full_sub)))
sens_idx  <- vapply(unique(sens_base[sens_idx]), function(b) {
  hits <- sens_idx[sens_base[sens_idx] == b]
  # Filenames sort chronologically, so the last after ordering is the newest
  tail(hits[order(basename(files_full_sub)[hits])], 1L)
}, integer(1))

message("Completeness sweep: ", length(sens_idx), " datasets")
message(paste(" -", basename(files_full_sub)[sens_idx], collapse = "\n"))

# Fixed label for this experiment's output folders (run parameters are hardcoded
# below; the per-cell completeness floor is what varies).
set <- "comp"

for (fn in sens_idx) {
  tryCatch({
    
    loaded_objects <- load(files_full_sub[fn])
    
    # Per-cell true completeness = observed distinct species / true richness.
    # FIX: use `<-` (was `<`, a comparison) and pass the loaded data frame `df`,
    # not `loaded_objects` (the character vector of object names from load()).
    sc_cells <- compute_site_completeness(df, true_richness)
    comp_stats <- compute_completeness_stats(sc_cells)
    meanSamplingCompletionOfAllCells <- comp_stats$comp_mean
    
    fnElem<-unlist(lapply(  unlist(strsplit(basename(files_full_sub[fn]), split="\\.")), function(x){strsplit(x, "_")}))
    
    # alternatve code for fnElem 
    # m <- regmatches(basename(files_full_sub[fn]),
    #             regexpr("^(?P<term>[^.]+)\\.(?P<scen>[^_]+)_(?P<level>[^_]+)_(?P<rep>\\d+)\\.RData$",
    #                     basename(files_full_sub[fn]), perl = TRUE))
    
    occurrenceStatusVector<-rep("PRESENT",length(df$scientificName))
    df[["occurrenceStatus"]]<-occurrenceStatusVector
    
    # 1. Transform JUST the geometry column (sf knows how to handle an sfc column directly)
    geom_wgs84 <- st_transform(df$geometry, 4326)
    
    # 2. Turn off s2 temporarily for lightning-fast flat math
    sf_use_s2(FALSE)
    centroids <- st_coordinates(st_centroid(geom_wgs84))
    sf_use_s2(TRUE) # turn it back on
    
    # 3. Use data.table syntax (:=) to add the columns instantly in-place
    df[, `:=`(
      decimalLongitude = centroids[, "X"],
      decimalLatitude  = centroids[, "Y"]
    )]
    
    # 4. Build a DarwinCore export table, keeping grid_id so we can subset cells
    #    by completeness floor below. grid_id is dropped before each CSV is written.
    export_all <- as.data.frame(
      df[, c("grid_id", "scientificName", "species",
             "decimalLatitude", "decimalLongitude", "occurrenceStatus")]
    )
    
    # Base folder identity for this dataset (scenario + level + set); the
    # completeness floor is appended per run below so runs never collide.
    base_folder <- paste0("sim_", unique(df$scenario), "_", 
                          if(length(fnElem)==7){ fnElem[3]}
                          else {fnElem[4]}
                          , "_",toString(set))
    
    # --- Sweep each minimum per-cell completeness floor ---------------------
    for (floor_c in completeness_floors) {
      
      # Cells whose TRUE completeness meets or exceeds the floor
      keep_cells <- sc_cells$grid_id[!is.na(sc_cells$completeness) &
                                       sc_cells$completeness >= floor_c]
      
      # A high floor can leave no qualifying cells; skip rather than run an
      # empty dataset (which would fail the model fit).
      if (length(keep_cells) == 0L) {
        message(base_folder, " floor=", floor_c,
                ": no cells meet the completeness floor. Skipping.")
        next
      }
      
      export_filt <- export_all[export_all$grid_id %in% keep_cells,
                                setdiff(names(export_all), "grid_id")]
      
      # Tag every path with the completeness floor (e.g. 0.7 -> "c07")
      floor_tag   <- paste0("c", sub("\\.", "", format(floor_c, nsmall = 1)))
      folder_name <- paste0(base_folder, "_", floor_tag)
      target_dir  <- file.path(animal.dir, folder_name)
      file_name   <- paste0(folder_name, "_Data.csv")
      full_path   <- file.path(target_dir, file_name)
      
      # Create the folder if it doesn't exist yet
      if (!dir.exists(target_dir)) {
        dir.create(target_dir, recursive = TRUE)
      }
      
      # Write the filtered intermediate file, but NEVER overwrite an existing one
      if (file.exists(full_path)) {
        message(folder_name, ": intermediate file already exists. Not overwriting.")
      } else {
        write.csv(export_filt, full_path, row.names = FALSE)
        message(folder_name, ": wrote ", nrow(export_filt), " rows from ",
                length(keep_cells), " of ", nrow(sc_cells),
                " cells (completeness floor ", floor_c, ").")
      }
      
      #undebug(run_animal_model) ###NOTE: If there are only a handful of sites and the overlap with plant data happens to be very small, then the model will fail with the error:
      #
      #Error in numericDeriv(form[[3L]], names(ind), env, central = nDcentral) : 
      #    Missing value or an infinity produced when evaluating the model
      outdir<-paste0("./p_output/",folder_name, "/animal_model/")
      existingOut<-list.files(outdir, pattern = "gamma_richnessfromSPRICH.*\\.txt$")
      if(length(existingOut)==1){
        message(paste0(
          folder_name, " simulation completed already. Continuing..."
        ))
        
      }else{
        #undebug(simulateScaledSpeciesRichness)
        #undebug(run_animal_model)
        run_animal_model(folder_name,
                         data_file = sub("^.*?data/", "data/", full_path),
                         clean = FALSE,
                         redo = FALSE,
                         # Completeness floor applied upstream via cell subsetting,
                         # so the model's own filters stay disabled.
                         min_occurrences = 1,
                         min_completeness = 0,
                         min_richness = 1,
                         estimateBounds=FALSE,
                         simulate_ranges = FALSE,
                         skipRawEndemism = TRUE,
                         customOutputFolderPrefix = "p_",
                         disable_outlier_filter = TRUE,
                         customOutputFolderSuffix = paste0("_",toString(set),"_",floor_tag),
                         fill_plants=FALSE,
                         fill_animals=FALSE,
                         fit_params = c("amplitude", "meanlog", "sdlog"),
                         filter_plant_completeness = NULL
        )
      }
    }
  }, error = function(e) {
    message("Error at iteration ", fn, ": ", conditionMessage(e))
  })
}


###COMBINE ALL RESULTS TOGETHER#####
fn<-3
source("AmAni.manager.R")
list_idx <- 0
for (fn in seq_along(files_full_sub)) {
  print(fn)
  loaded_objects <- load(files_full_sub[fn])
  
  fnElem<-unlist(lapply(  unlist(strsplit(basename(files_full_sub[fn]), split="\\.")), function(x){strsplit(x, "_")}))
  
  # alternatve code for fnElem 
  # m <- regmatches(basename(files_full_sub[fn]),
  #             regexpr("^(?P<term>[^.]+)\\.(?P<scen>[^_]+)_(?P<level>[^_]+)_(?P<rep>\\d+)\\.RData$",
  #                     basename(files_full_sub[fn]), perl = TRUE))
  
  # The test name always sits immediately before "sample" in fnElem,
  # regardless of how many run/replicate tokens precede it.
  test_name_fromFile   <- fnElem[which(fnElem == "sample") - 1]
  search_term_fromFile <- test_name_fromFile
  
  # If df failed to load (NULL), build an all-NA row for every set and move on
  if (is.null(df)) {
    message("df is NULL for file: ", basename(files_full_sub[fn]), " — recording NA rows.")
    prev_valid <- Filter(Negate(is.null), comparisonList)
    for (set in seq_len(nrow(settingCombinations))) {
      list_idx <- list_idx + 1
      if (length(prev_valid) > 0) {
        na_row <- prev_valid[[length(prev_valid)]][1, , drop = FALSE]
        na_row[1, ] <- NA
        na_row[[1]] <- search_term_fromFile
      } else {
        na_row <- tibble(search_term = search_term_fromFile)
      }
      comparisonList[[list_idx]] <- bind_cols(na_row, settingCombinations[set, ])
    }
    next
  }
  
  # names(df)
  # unique(df$scenario )
  # head(df)
  # 
  occurrenceStatusVector<-rep("PRESENT",length(df$scientificName))
  df[["occurrenceStatus"]]<-occurrenceStatusVector
  
  # 1. Transform JUST the geometry column (sf knows how to handle an sfc column directly)
  geom_wgs84 <- st_transform(df$geometry, 4326)
  
  # 2. Turn off s2 temporarily for lightning-fast flat math
  sf_use_s2(FALSE)
  centroids <- st_coordinates(st_centroid(geom_wgs84))
  sf_use_s2(TRUE) # turn it back on
  
  # 3. Use data.table syntax (:=) to add the columns instantly in-place
  df[, `:=`(
    decimalLongitude = centroids[, "X"],
    decimalLatitude  = centroids[, "Y"]
  )]
  
  # 4. Subset for DarwinCore Export
  # We drop 'geometry' because CSVs cannot store spatial polygons
  final_export <- df[, c("scientificName","species", "decimalLatitude", "decimalLongitude","occurrenceStatus")]
  
  for (set in seq_len(nrow(settingCombinations))) {
    list_idx <- list_idx + 1
  
  # 1. Build the folder name and target directory path
  folder_name <- paste0("sim_", unique(df$scenario), "_", test_name_fromFile, "_", toString(set))
  target_dir  <- file.path(animal.dir, folder_name)
  
  # 2. Build the full file path where the CSV will live
  file_name   <- paste0(folder_name, "_Data.csv")
  full_path   <- file.path(target_dir, file_name)
  
  
    # Try single-suffix first, then double-suffix (from customOutputFolderSuffix)
    outdir <- paste0("./p_output/", folder_name, "/animal_model/")
    if (!dir.exists(outdir)) {
      folder_name_alt <- paste0(folder_name, "_", toString(set))
      outdir <- paste0("./p_output/", folder_name_alt, "/animal_model/")
    }
    scenarioFromOut <- unique(df$scenario)  # e.g. "10_1", "10_2", "high" — use df directly to avoid mis-parsing
    
    #get results from the predictive model
    gammaFileInfo <- extract_gamma_file(outdir)
    
    #get the summary of the input dataset from the predictive model
    datasetSummary<-extract_dataset_summary_file(outdir)
    
    #Get the summary of the input dataset from the sampling procedure
      # Extract the search term ("variEndemSampling")
      search_term <- test_name_fromFile
    
      # Find the file using correct Regex
      summaryInfoFile<-list.files(
      path = data.dir, 
      pattern = paste0("^", search_term, "\\.summary_", 
                       strsplit(strsplit( basename(files_full_sub[fn]), split ="sample_" )[[1]][2], split =".RData" )
                       ,"\\.csv$"), 
      full.names = TRUE
      )
      if (length(summaryInfoFile) == 0) {
        message("No dataset_summary files found in: ", data.dir)
        result_row<-NA
        truth_row<-NA
        trueSpeciesRichness<- NA
        }else{
      summaryInfo<-read.csv(summaryInfoFile)
      result_row <- summaryInfo[summaryInfo$scenario == scenarioFromOut, , drop = FALSE]
    
    #Get the summary of the true dataset
      # 1. Define the list of terms used for full dataset
      true_full_terms <- c(
        "lowtoHighSpeciesSampling",
        "lowtoHighSiteSampling",
        "lowtoHighSiteEnd",
        "lowtoHighSiteRichness",
        "variedSiteClustering",
        "variedSiteEvenness",
        "roachLike",
        "beetleLike",
        "butterflyLike"
      )
      
      # 2. Check and set the variable
      if (search_term %in% true_full_terms) {
        datasetType <- "trueFull"
        datasetlevel<-NULL
        
      } else {
      datasetlevel<-unlist(strsplit(scenarioFromOut, "\\."))[2]
      datasetType<-unlist(strsplit(scenarioFromOut, "\\."))[1]
      }
      
        #create blank truth row
      truthfileEx  <- "./5-SensitivityTesting/output/trueDatasets/endemismSamples.summary_20260625_082631.csv"
      truthInfo  <- read.csv(truthfileEx)
        truth_row  <- truthInfo[1, , drop = FALSE]
        truth_row[] <- NA
      
      switch(datasetType,
             "Endemism" = {
               truthfile  <- "./5-SensitivityTesting/output/trueDatasets/endemismSamples.summary_20260625_082631.csv"
               truthInfo  <- read.csv(truthfile)
               truth_row  <- truthInfo[truthInfo$scenario == datasetlevel, , drop = FALSE]
               trueSpeciesRichness<-truth_row["n_species"]
               
             },
             
             "Richness" = {
               truthfile  <- "./5-SensitivityTesting/output/trueDatasets/richnessSamples.summary_20260625_075324.csv"
               truthInfo  <- read.csv(truthfile)
               truth_row  <- truthInfo[truthInfo$scenario == datasetlevel, , drop = FALSE]
               trueSpeciesRichness<-truth_row["n_species"]
             },
             
             "trueFull" = {
               truth_row[] <- NA
               trueSpeciesRichness <- extract_gamma_file("./5-SensitivityTesting/data/simFromPlant")$best_estimate 
             }
      )
    
    #Combine all metadata for analysis
      }
      # Prefix overlapping columns to avoid name collisions across iterations
      if (!all(is.na(result_row))) {
        names(result_row) <- paste0("sample_", names(result_row))
      }
      if (!all(is.na(truth_row))) {
        names(truth_row) <- paste0("true_", names(truth_row))
      }
      
      comparisonList[[list_idx]]<-bind_cols(
        search_term,
          gammaFileInfo,
          datasetSummary,
          result_row,
          truth_row,
          tibble(trueSpeciesRichness = trueSpeciesRichness), # Converts the single value to a column
          settingCombinations[set, ],
          tibble(input_file = basename(files_full_sub[fn]))
        )
      
      # Also import results from archived runs in this output folder
      archive_parent <- file.path(outdir, "archive")
      if (dir.exists(archive_parent)) {
        archive_dirs <- list.dirs(archive_parent, recursive = FALSE)
        for (arch_dir in archive_dirs) {
          arch_gamma   <- extract_gamma_file(arch_dir)
          arch_summary <- extract_dataset_summary_file(arch_dir)
          # Skip if the archive produced no gamma file
          if (all(is.na(arch_gamma))) next
          list_idx <- list_idx + 1
          comparisonList[[list_idx]] <- bind_cols(
            search_term,
            arch_gamma,
            arch_summary,
            result_row,
            truth_row,
            tibble(trueSpeciesRichness = trueSpeciesRichness),
            settingCombinations[set, ],
            tibble(input_file = basename(files_full_sub[fn]))
          )
        }
      }

  } # end set loop

  # --- Scrape the completeness-sweep runs for this dataset -------------------
  # The sweep above filtered grid cells by TRUE completeness before modeling and
  # wrote each run to a folder tagged with its floor (e.g. "..._comp_c07").
  # Reconstruct that folder identity exactly as the sweep did, pull the model
  # results, and tag each row with the completeness floor plus the number of grid
  # cells before and after filtering.
  comp_base_folder <- paste0("sim_", unique(df$scenario), "_",
                             if (length(fnElem) == 7) fnElem[3] else fnElem[4],
                             "_comp")

  # TRUE per-cell completeness for this dataset; the full cell count is the
  # pre-filter total, and each floor determines how many cells survive.
  sc_cells_scrape <- compute_site_completeness(df, true_richness)
  n_cells_before  <- nrow(sc_cells_scrape)

  for (floor_c in completeness_floors) {
    # Same floor -> tag mapping as the sweep (0.7 -> "c07")
    floor_tag   <- paste0("c", sub("\\.", "", format(floor_c, nsmall = 1)))
    comp_folder <- paste0(comp_base_folder, "_", floor_tag)
    # run_animal_model appended customOutputFolderSuffix ("_comp_{floor_tag}") to
    # the folder name, producing a double-tagged path: comp_c03_comp_c03.
    comp_outdir <- paste0("./p_output/", comp_folder, "_comp_", floor_tag, "/animal_model/")

    # This dataset/floor combination may not have been run; skip quietly
    if (!dir.exists(comp_outdir)) next

    comp_gamma <- extract_gamma_file(comp_outdir)
    if (all(is.na(comp_gamma))) next
    comp_summary <- extract_dataset_summary_file(comp_outdir)

    # Cells retained at this floor (matches the sweep's keep_cells rule)
    n_cells_after <- sum(!is.na(sc_cells_scrape$completeness) &
                           sc_cells_scrape$completeness >= floor_c)

    list_idx <- list_idx + 1
    comparisonList[[list_idx]] <- bind_cols(
      search_term_fromFile,
      comp_gamma,
      comp_summary,
      result_row,
      truth_row,
      tibble(trueSpeciesRichness = trueSpeciesRichness),
      tibble(
        completeness_floor  = floor_c,
        n_grid_cells_before = n_cells_before,
        n_grid_cells_after  = n_cells_after
      ),
      tibble(input_file = basename(files_full_sub[fn]))
    )
  }
}

saveRDS(comparisonList, "./5-SensitivityTesting//results//p_FULLRAW_comparisonList.rds")
comparisonList[40]$fill_plants

names(comparisonList[[39]])
names(comparisonList[[40]])
summary(comparisonList)
names(comparisonList2)
comparisonList2 <- comparisonList[sapply(comparisonList, function(x) {
  # Keep any entry that has a best_estimate, even if sample_scenario is missing
  # (replicates whose summary CSV timestamp didn't match are still valid model outputs)
  is.data.frame(x) && "best_estimate" %in% names(x)
})]
#comparisonList3<-comparisonList2[!is.na(lapply(comparisonList2, function(x) x$best_estimate))] 
#summary(comparisonList2)
#summary(comparisonList3)


# More robust approach that avoids the warning
comparison_list_clean <- comparisonList2[!sapply(comparisonList2, is.null)]
summary(comparison_list_clean)
#comparison_list_clean
# Manually bind with more control
# 
#x$<-comparison_list_clean
names(comparison_list_clean[[1]])
comparison_listFull <- do.call(rbind, lapply(comparison_list_clean, function(x) {
  if (nrow(x) == 0) {
    return(NULL)
  }
  
  # Handle trueSpeciesRichness which may be nested data frame or numeric
  true_rich <- x$trueSpeciesRichness
  if (is.data.frame(true_rich)) {
    true_rich <- true_rich[[1, 1]]
  }
  
  get_col <- function(df, col, default = NA) {
    if (col %in% names(df)) df[[col]] else default
  }
  sc <- coalesce(get_col(x, "sample_scenario"), get_col(x, "sample_original_scenario_name"))
  data.frame(
    scenario = paste0(sc, "-", x$...1),
    test = x$...1,
    best_estimate = x$best_estimate,
    lower_bound = x$lower_bound,
    upper_bound = x$upper_bound,
    param_names_fit = get_col(x, "param_names_fit"),
    delta_meanlog = get_col(x, "delta_meanlog"),
    delta_meanlog_status = get_col(x, "delta_meanlog_status"),
    ratio_sdlog = get_col(x, "ratio_sdlog"),
    sample_endemism_min = get_col(x, "sample_endemism_min"),
    sample_endemism_max = get_col(x, "sample_endemism_max"),
    sample_median_range_size = get_col(x, "sample_median_range_size"),
    sample_mean_range_size = get_col(x, "sample_mean_range_size"),
    sample_sd_range_size = get_col(x, "sample_sd_range_size"),
    n_species = coalesce(get_col(x, "sample_n_species"), NA_integer_),
    n_grid_cells = coalesce(get_col(x, "sample_n_grid_cells"), NA_integer_),
    total_observations = coalesce(get_col(x, "sample_total_observations"), NA_integer_),
    prop_species = get_col(x, "sample_prop_species"),
    prop_sites = get_col(x, "sample_prop_sites"),
    clustering_label = get_col(x, "sample_clustering_label"),
    sample_numberStartSites = get_col(x, "sample_numberStartSites"),
    sample_diffusion = get_col(x, "sample_diffusion"),
    sampling_fraction = get_col(x, "sample_sampling_fraction"),
    sampling_evenness = get_col(x, "sample_sampling_evenness"),
    min_completeness = get_col(x, "min_completeness"),
    completeness_floor = get_col(x, "completeness_floor"),
    n_grid_cells_before = get_col(x, "n_grid_cells_before"),
    n_grid_cells_after = get_col(x, "n_grid_cells_after"),
    
    mean_range_size = coalesce(get_col(x, "mean_range_size"), NA_integer_),
    mean_rangeSize = coalesce(get_col(x, "mean_rangeSize"), NA_integer_),
    median_range_size = coalesce(get_col(x, "median_range_size"), NA_integer_),
    sd_range_size = coalesce(get_col(x, "sd_range_size"), NA_integer_),
    trueSpeciesRichness = true_rich,
    fill_plants = get_col(x, "fill_plants"),
    fill_animals = get_col(x, "fill_animals"),
    filter_plant_completeness = {
      val <- get_col(x, "filter_plant_completeness")
      if (is.list(val)) val <- val[[1]]
      if (is.null(val)) NA else val
    },
    fit_params = {
      val <- get_col(x, "fit_params")
      if (is.list(val)) val <- val[[1]]
      if (is.null(val)) NA_character_ else toString(val)
    },
    input_file = get_col(x, "input_file"),
    stringsAsFactors = FALSE
  )
}))

comparison_listFull
str(comparison_listFull)
write.csv(comparison_listFull, "./5-SensitivityTesting//results//p_comparison_results.csv", row.names = FALSE)



comparison_df <- do.call(rbind, lapply(comparison_list_clean, function(x) {
  if (nrow(x) == 0) {
    return(NULL)
  }
  
  # Handle trueSpeciesRichness which may be nested data frame or numeric
  true_rich <- x$trueSpeciesRichness
  if (is.data.frame(true_rich)) {
    true_rich <- true_rich[[1, 1]]
  }
  
  get_col <- function(df, col, default = NA) {
    if (col %in% names(df)) df[[col]] else default
  }
  data.frame(
    test = x$...1,
    scenario = coalesce(get_col(x, "sample_scenario"), get_col(x, "sample_original_scenario_name")),
    best_estimate = x$best_estimate,
    trueSpeciesRichness = true_rich,
    stringsAsFactors = FALSE
  )
}))

# Calculate difference
comparison_df$difference <- comparison_df$best_estimate - comparison_df$trueSpeciesRichness
comparison_df$diff_prop <- comparison_df$difference/comparison_df$trueSpeciesRichness


# Reset row names for clarity
rownames(comparison_df) <- NULL

print(comparison_df)



# Keep only the columns we need
plot_df <- comparison_df[, c("test","scenario", "best_estimate", "trueSpeciesRichness", "difference", "diff_prop")]

print(plot_df)

# Helper: assign a sort rank to a scenario label.
# "true", "sample", "GBIF" (case-insensitive) go to the end;
# common ordinal words get numeric ranks; others sort alphabetically in the middle.
scenario_sort_key <- function(scen) {
  # Strip any "Xxx." prefix (e.g. "Endemism.") to get the core label
  core <- sub("^[A-Za-z]+\\.", "", scen)
  ordinal_rank <- c(
    lowest = 1, min = 1, minimum = 1, 
    verylow  = 2, vlow = 2, quarter = 1,
    low      = 3, 
    medium   = 4, med = 4, even = 4, half = 3,
    high     = 5,
    higher = 6, veryhigh = 6,
    highest  = 7
  )
  end_labels <- c("true", "sample", "gbif")
  rank <- vapply(tolower(core), function(s) {
    if (s %in% end_labels) return(1000L)
    r <- ordinal_rank[s]
    if (!is.na(r)) return(as.integer(r * 10))
    # Plain numeric labels (e.g. "10", "100")
    num <- suppressWarnings(as.numeric(s))
    if (!is.na(num)) return(as.integer(num * 10))
    # Replicate-numbered numeric labels (e.g. "10_1", "20_6"): sort by base level, then replicate
    num_base <- suppressWarnings(as.numeric(sub("_[0-9]+$", "", s)))
    rep_idx  <- suppressWarnings(as.integer(sub("^[0-9]+_", "", s)))
    if (!is.na(num_base) && !is.na(rep_idx)) return(as.integer(num_base * 10000L + rep_idx))
    # fallback: alphabetical via utf8 value of first char
    500L + utf8ToInt(substr(s, 1, 1))
  }, integer(1))
  rank
}

# Single-scenario organism-type tests are merged into one combined facet
single_scen_tests <- c("beetleLike", "butterflyLike", "roachLike")

# Build ordered scenario factor within each facet group.
# Replicates (same test + scenario) get a numeric suffix so they appear as
# separate bars rather than being aggregated at the same x position.
plot_df_ordered <- plot_df |>
  filter(!is.na(scenario)) |>
  mutate(
    facet_group = ifelse(test %in% single_scen_tests, "organismType", test),
    scen_rank   = scenario_sort_key(scenario)
  ) |>
  group_by(facet_group, scenario) |>
  mutate(
    rep_num = row_number(),
    n_reps  = n()
  ) |>
  ungroup() |>
  mutate(
    # Reformat "10_1" → "10run1" for display; leave word-based scenarios unchanged
    scenario_display = if_else(
      grepl("^[0-9]+_[0-9]+$", scenario),
      sub("^([0-9]+)_([0-9]+)$", "\\1run\\2", scenario),
      scenario
    ),
    # For word-based scenarios with true duplicates, append [rep] suffix
    scenario_label = if_else(
      n_reps > 1 & !grepl("^[0-9]+_[0-9]+$", scenario),
      paste0(scenario_display, " [", rep_num, "]"),
      scenario_display
    ),
    # Sort key: scenario level first, then replicate index
    sort_key       = scen_rank * 1000L + rep_num,
    scenario_label = reorder(scenario_label, sort_key)
  )

# ── Plot 1: Best Estimate vs Actual, faceted by test ──────────────────────────
plot_data1 <- plot_df_ordered |>
  pivot_longer(
    cols      = c(best_estimate, trueSpeciesRichness, difference),
    names_to  = "metric",
    values_to = "value"
  ) |>
  mutate(
    metric = factor(recode_values(metric,
      "best_estimate"       ~ "Best Estimate",
      "trueSpeciesRichness" ~ "True Species Richness",
      "difference"          ~ "Difference (Est - Actual)"
    ), levels = c("Best Estimate", "True Species Richness", "Difference (Est - Actual)"))
  )

p1 <- ggplot(plot_data1, aes(x = scenario_label, y = value, fill = metric)) +
  geom_col(position = "dodge") +
  facet_wrap(~ facet_group, scales = "free_x") +
  labs(
    title = "Comparison of Best Estimate vs Actual Species Count",
    x     = "Scenario",
    y     = "Value",
    fill  = "Metric"
  ) +
  theme_minimal() +
  theme(
    axis.text.x      = element_text(angle = 45, hjust = 1),
    strip.text       = element_text(face = "bold"),
    panel.spacing    = unit(1, "lines")
  )

# ── Plot 2: Proportional difference, faceted by test ─────────────────────────
plot_data2 <- plot_df_ordered |>
  pivot_longer(
    cols      = diff_prop,
    names_to  = "metric",
    values_to = "value"
  ) |>
  mutate(metric = "Proportional difference (difference / actual)")

p2 <- ggplot(plot_data2, aes(x = scenario_label, y = value, fill = metric)) +
  geom_col(position = "dodge") +
  facet_wrap(~ facet_group, scales = "free_x") +
  labs(
    title = "Relative difference (0 is best, −1 or 1 is worst)",
    x     = "Scenario",
    y     = "Value",
    fill  = "Metric"
  ) +
  theme_minimal() +
  theme(
    axis.text.x   = element_text(angle = 45, hjust = 1),
    strip.text    = element_text(face = "bold"),
    panel.spacing = unit(1, "lines"),
    legend.position = "none"
  )

# ── PDF output ────────────────────────────────────────────────────────────────
pdf_dir  <- "C:/Users/dominice/OneDrive - University of Illinois - Urbana/Projects/Peruvian Blattodea/Cockroaches in a coal mine/AmAni/5-SensitivityTesting/results"
pdf_path <- file.path(pdf_dir, paste0("p_sensitivity_plots_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".pdf"))

pdf(pdf_path, width = 20, height = 8)
print(p1)
print(p2)
dev.off()

message("Plots saved to: ", pdf_path)
