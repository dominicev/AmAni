#This notebook will take the simulated sampling from the simulated truth datasets and run the animal model for them. I won't run the extinction model (yet?) because I haven't simulated the true extinction. I currently don't think it fully makes sense to simulate extinction. (a) the reasons the extinction simulator would be inaccurate are not things the simulation would currently account for (species ranges changing, unknown pre-deforestation diversity levels).
library(sf)

source("AmAni.summarize.R")

#SOME NOTES AFTER FIRST ROUND OF TESTS
#1. REPLICATES WITH LOW NUMBER OF SITES TEND TO FAIL BECAUSE OF LOW OVERLAP WITH THE PLANT DATA #MADE MORE DATASETS
# I CAN FIX THIS BY INCREASING THE MINIMUM NUMBER OF SITES
#2. I NEED MORE REPLICATES OF THE FOLLOWING TESTS
#   A. SITE ENDEMISM #MADE DATASETS
#   B. SITE RICHNESS #MADE DATASETS
#   C. MEDIUM AND HIGH SITE CLUSTERING #MADE DATASETS



getwd()
# --- import 1 dataset to confirm it can run on the format ---
animal.dir<-"./2-ModellingAnimals/data"
data.dir<-"./5-SensitivityTesting/output/sampleDatasets/"

files_full <- list.files(data.dir, pattern = "\\.RData$", full.names = TRUE)
files_short <- list.files(data.dir, pattern = "\\.RData$", full.names = FALSE)
length(files_full)
files_full[41]
comparisonList<-list()



fn<-1

source("AmAni.manager.R")
for (fn in 1:length(files_full)) {
  tryCatch({
  
    loaded_objects <- load(files_full[fn])
    fnElem<-unlist(lapply(  unlist(strsplit(basename(files_full[fn]), split="\\.")), function(x){strsplit(x, "_")}))
    
    # alternatve code for fnElem 
    # m <- regmatches(basename(files_full[fn]),
    #             regexpr("^(?P<term>[^.]+)\\.(?P<scen>[^_]+)_(?P<level>[^_]+)_(?P<rep>\\d+)\\.RData$",
    #                     basename(files_full[fn]), perl = TRUE))
    
    
    
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
                          )
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
outdir<-paste0("./output/",folder_name, "/animal_model/")
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
                     min_occurrences = 1,
                     min_completeness = 0,
                       min_richness = 1,
                       filter_plant_completeness = 0,
                     disable_outlier_filter = TRUE,
                       fit_params = c("amplitude", "meanlog","sdlog"),
                     estimateBounds=FALSE,
                     skipRawEndemism = FALSE
                       )
}
  }, error = function(e) {
    message("Error at iteration ", fn, ": ", conditionMessage(e))
  })
}



fn<-3
source("AmAni.manager.R")
for (fn in 1:length(files_full)) {
  print(fn)
  loaded_objects <- load(files_full[fn])
  
  fnElem<-unlist(lapply(  unlist(strsplit(basename(files_full[fn]), split="\\.")), function(x){strsplit(x, "_")}))
  
  # alternatve code for fnElem 
  # m <- regmatches(basename(files_full[fn]),
  #             regexpr("^(?P<term>[^.]+)\\.(?P<scen>[^_]+)_(?P<level>[^_]+)_(?P<rep>\\d+)\\.RData$",
  #                     basename(files_full[fn]), perl = TRUE))
  
  # The test name always sits immediately before "sample" in fnElem,
  # regardless of how many run/replicate tokens precede it.
  test_name_fromFile   <- fnElem[which(fnElem == "sample") - 1]
  search_term_fromFile <- test_name_fromFile
  
  # If df failed to load (NULL), build an all-NA row and move on
  if (is.null(df)) {
    message("df is NULL for file: ", basename(files_full[fn]), " — recording NA row.")
    prev_valid <- Filter(Negate(is.null), comparisonList)
    if (length(prev_valid) > 0) {
      na_row <- prev_valid[[length(prev_valid)]][1, , drop = FALSE]
      na_row[1, ] <- NA
      na_row[[1]] <- search_term_fromFile  # restore file-name-derived value
    } else {
      na_row <- tibble(search_term = search_term_fromFile)
    }
    comparisonList[[fn]] <- na_row
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
  
  # 1. Build the folder name and target directory path
  folder_name <- paste0("sim_", unique(df$scenario), "_", test_name_fromFile)
  target_dir  <- file.path(animal.dir, folder_name)
  
  # 2. Build the full file path where the CSV will live
  file_name   <- paste0(folder_name, "_Data.csv")
  full_path   <- file.path(target_dir, file_name)
  
  
    outdir<-paste0("./output/", folder_name,"/animal_model/")
    scenarioFromOut <- unique(df$scenario)  # e.g. "10_1", "10_2", "high" — use df directly to avoid mis-parsing
    
    #get results from the predictive model
    gammaFileInfo <- extract_gamma_file(outdir)
    
    #get the summary of the input dataset from the predictive model
    datasetSummary<-extract_dataset_summary_file(outdir)
    
    #Get the summary of the input dataset from the sampling procedure
      # Extract the search term ("variEndemSampling")
      search_term <- tail(unlist(strsplit(folder_name, "_")), 1)
    
      # Find the file using correct Regex
      summaryInfoFile<-list.files(
      path = data.dir, 
      pattern = paste0("^", search_term, "\\.summary_", 
                       strsplit(strsplit( basename(files_full[fn]), split ="sample_" )[[1]][2], split =".RData" )
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
      
      comparisonList[[fn]]<-bind_cols(
        search_term,
          gammaFileInfo,
          datasetSummary,
          result_row,
          truth_row,
          tibble(trueSpeciesRichness = trueSpeciesRichness) # Converts the single value to a column
        )
      
}

saveRDS(comparisonList, "./5-SensitivityTesting//results//FULLRAW_comparisonList.rds")


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
summary(comparisonList3)


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
    
    mean_range_size = coalesce(get_col(x, "mean_range_size"), NA_integer_),
    mean_rangeSize = coalesce(get_col(x, "mean_rangeSize"), NA_integer_),
    median_range_size = coalesce(get_col(x, "median_range_size"), NA_integer_),
    sd_range_size = coalesce(get_col(x, "sd_range_size"), NA_integer_),
    trueSpeciesRichness = true_rich,
    stringsAsFactors = FALSE
  )
}))

str(comparison_listFull)
write.csv(comparison_listFull, "./5-SensitivityTesting//results//comparison_results.csv", row.names = FALSE)



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
    metric = recode(factor(metric,
                           levels = c("best_estimate", "trueSpeciesRichness", "difference")),
      best_estimate        = "Best Estimate",
      trueSpeciesRichness  = "True Species Richness",
      difference           = "Difference (Est - Actual)"
    )
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
pdf_path <- file.path(pdf_dir, paste0("sensitivity_plots_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".pdf"))

pdf(pdf_path, width = 20, height = 8)
print(p1)
print(p2)
dev.off()

message("Plots saved to: ", pdf_path)
