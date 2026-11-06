# ============================================================
# AmAni.manager.R — Master workflow coordinator
# ============================================================
#
# Pipeline overview
# -----------------
#
# Stage 0: Curating Plant Data (0-CuratingPlantData/)
#   - Download GBIF plant occurrence records for the Amazon basin
#   - Filter occurrences spatially to Amazon boundary shapefile
#   - Standardize taxonomy against World Flora Online backbone
#   - Aggregate species counts onto a 22.36 km grid (EPSG:3857)
#   - Calculate raw and weighted endemism per grid cell
#   - Key outputs: species_counts3.rds, species_endemism2.rds
#
# Stage 1: Filtering Plant Data (1-FilteringPlantData/)
#   - Apply outlier filters (95th percentile + 3x IQR) to endemism values
#   - Remove grid cells affected by large cities (pop >= 10,000)
#   - Produce filtered species counts, raw endemism, and weighted endemism
#   - Key output: FINAL_allEndemism_FullyFiltered.rds (one row per grid cell)
#
# Stage 2: Modelling Animals (2-ModellingAnimals/)
#   - Import animal occurrence data (cockroaches, etc.)
#   - Map occurrences onto the same plant-derived grid
#   - Calculate animal species richness and endemism per grid cell
#   - Fit spatial extent models and simulate occurrence counts via gamma scaling
#   - Key output: simulated occurrence tables for extinction modelling
#
# Stage 3: Extinction Predictor (3-ExtPredictor/)
#   - Extract deforestation patch geometry from RAISG 2001-2020 raster
#   - Tile raster into 3,630 overlapping windows; extract patches via
#     morphological operations and connected-component labeling
#   - Fit lognormal and power-law distributions to patch area
#   - Key outputs: patch area frequency tables and model parameters
#
# Stage 4: Extinction Modeller (4-ExtModeller/)
#   - Simulate species extinction under deforestation scenarios using
#     patch size distributions (Stage 3) and species grids (Stages 1-2)
#   - Three spatial models: neutral (random), clustered (observed-weighted),
#     and shifted (toroidal null)
#   - 100 replicates per scenario x multiple frequency table variants
#   - Key outputs: extinction counts, per-species status, summary comparisons
#
# Data flow:
#   GBIF records -> [0] filter/grid -> [1] quality filter ->
#     [2] animal overlay + simulation -> [4] extinction simulation
#   RAISG raster -> [3] patch extraction -> [4] extinction simulation
#
# Stages 0, 1, and 3 produce stable upstream data that rarely changes.
# The functions below re-run the animal (Stage 2) and extinction (Stage 4)
# steps while drawing on those pre-existing outputs. Results are stored
# per-organism under output/{organism}/.
#
# ============================================================
# Usage:
#   source("AmAni.manager.R")
#   check_status()                          # overview of all organisms
#   check_status("cockroaches")             # detailed status for one organism
#   run_animal_model("cockroaches")         # re-run Stage 2
#   run_extinction_model("cockroaches")     # re-run Stage 4 (all models)
#   run_extinction_model("cockroaches",
#                        models = "neutral") # run only neutral model
# ============================================================

# --- Project paths -----------------------------------------------------------

project_root <- dirname(rstudioapi::getSourceEditorContext()$path)

stage_dirs <- list(
  curate   = file.path(project_root, "0-CuratingPlantData"),
  filter   = file.path(project_root, "1-FilteringPlantData"),
  animals  = file.path(project_root, "2-ModellingAnimals"),
  predict  = file.path(project_root, "3-ExtPredictor"),
  model    = file.path(project_root, "4-ExtModeller")
)

output_root <- file.path(project_root, "output")

# --- Load GBIF cleaning utility ---------------------------------------------
source(file.path(stage_dirs$animals, "scripts", "00_cleanGBIF.R"))

# --- Load per-organism IUCN threat classification ----------------------------
source(file.path(stage_dirs$model, "scripts", "utils", "classify_iucn_organism.R"))


# --- Helper: verify upstream outputs ----------------------------------------

verify_upstream <- function(stage_name, required_files) {
  missing <- required_files[!file.exists(required_files)]
  if (length(missing) > 0) {
    stop(
      "Missing upstream outputs for ", stage_name, ":\n",
      paste(" -", missing, collapse = "\n"),
      "\nRe-run the upstream stage(s) before proceeding.",
      call. = FALSE
    )
  }
  message(stage_name, ": all ", length(required_files), " upstream files verified.")
  invisible(TRUE)
}


# --- Helper: set stage working directory for sourced scripts ------------------

set_stage_wd <- function(stage_dir) {
  .amani_wd <<- stage_dir
  setwd(stage_dir)
}

restore_wd <- function() {
  for (v in c(".amani_wd", ".amani_output_dir", ".amani_animal_model_dir", ".amani_data_file",".amani_estimateBounds", ".amani_gammaErrorType", ".amani_redo", ".amani_min_occurrences", ".amani_min_completeness", ".amani_min_richness", ".amani_max_richness", ".amani_disable_outlier_filter", ".amani_filter_plant_completeness", ".amani_filter_summary", ".amani_fit_params", ".amani_alphaRichnessMethod", ".amani_animalFilterType", ".amani_delta_meanlog", ".amani_ratio_sdlog", ".amani_plant_mdd_meanlog", ".amani_fill_plants")) {
    if (exists(v, envir = .GlobalEnv)) rm(list = v, envir = .GlobalEnv)
  }
  setwd(project_root)
}


# --- Helper: archive old outputs before a fresh run --------------------------

archive_outputs <- function(output_dir) {
  # Move all existing files in output_dir into archive/<timestamp>/
  # so fresh runs start with a clean folder.
  
  if (!dir.exists(output_dir)) return(invisible(NULL))
  
  existing <- list.files(output_dir, all.files = FALSE, full.names = TRUE)
  # Exclude the archive directory itself
  existing <- existing[basename(existing) != "archive"]
  
  if (length(existing) == 0) {
    message("  No existing outputs to archive in: ", output_dir)
    return(invisible(NULL))
  }
  
  stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  archive_dir <- file.path(output_dir, "archive", stamp)
  dir.create(archive_dir, recursive = TRUE)
  
  moved <- file.rename(existing, file.path(archive_dir, basename(existing)))
  n_ok <- sum(moved)
  message("  Archived ", n_ok, "/", length(existing),
          " items to: ", archive_dir)
  invisible(archive_dir)
}


# --- Helper: resolve organism output paths -----------------------------------

organism_dirs <- function(organism) {
  base <- file.path(output_root, organism)
  list(
    animal_model     = file.path(base, "animal_model"),
    extinction_model = file.path(base, "extinction_model")
  )
}

# --- Helper: extract gamma richness parameters from output -----------------------------------
# Function to extract gamma richness parameters from SPRICH output file
extract_gamma_file <- function(directory, verbose = TRUE) {
  
  # Find all matching files
  files <- list.files(directory, 
                      pattern = "gamma_richnessfromSPRICH.*\\.txt$",
                      full.names = TRUE)
  
  if (length(files) == 0) {
    message("No gamma_richnessfromSPRICH files found in: ", directory)
    result <- data.frame(
      param_names_fit = NA,
      extent_scaling_factor = NA,
      delta_meanlog = NA,
      delta_meanlog_status = NA,
      ratio_sdlog = NA,
      ratio_sdlog_status = NA,
      best_estimate = NA,
      upper_bound = NA,
      lower_bound = NA,
      stringsAsFactors = FALSE
    )
  }else{
  
  # Use most recent if multiple files exist
  if (length(files) > 1) {
    file_info <- file.info(files)
    target_file <- rownames(file_info)[which.max(file_info$mtime)]
  } else {
    target_file <- files[1]
  }
  
  if (verbose) {
    message("Reading: ", basename(target_file))
  }
  
  # Read file
  lines <- readLines(target_file)
  
  # Extract fit parameters line
  param_line <- grep("^Fit parameters:", lines, value = TRUE)
  if (length(param_line) == 0) {
    stop("Could not find 'Fit parameters:' line in file")
  }
  
  # Parse parameters (e.g., "[amplitude, meanlog, sdlog]")
  param_str <- gsub(".*\\[(.*)\\].*", "\\1", param_line)
  param_names <- trimws(strsplit(param_str, ",")[[1]])
  
  # Extract extent scaling factor
  extent_line <- grep("^Extent scaling factor:", lines, value = TRUE)
  extent_value <- as.numeric(gsub(".*: ([0-9.]+).*", "\\1", extent_line))
  
  # Helper function to extract numeric value after colon
  extract_numeric <- function(pattern, lines) {
    line <- grep(pattern, lines, value = TRUE)
    numeric<-as.numeric(strsplit(trimws(strsplit(line, ":")[[1]][[2]]), " ")[[1]][[1]])
    return(numeric)
  }
  
  # Helper function to extract status (active/inactive) after numeric value
  extract_status <- function(pattern, lines) {
    line <- grep(pattern, lines, value = TRUE)
    if (length(line) == 0) return(NA)
    # Extract text in parentheses
    status <- gsub(".*\\(([^)]+)\\).*", "\\1", line[1])
    return(if (grepl("\\(", line[1])) status else NA_character_)
  }
  
  # Extract delta_meanlog and status
  delta_meanlog <- extract_numeric("^delta_meanlog:", lines)
  delta_meanlog_status <- extract_status("^delta_meanlog:", lines)
  
  # Extract ratio_sdlog and status
  ratio_sdlog <- extract_numeric("^ratio_sdlog:", lines)
  ratio_sdlog_status <- extract_status("^ratio_sdlog:", lines)
  
  
  # Extract best estimate (number of species)
  best_line <- grep("^Best estimate:", lines, value = TRUE)
  best_value <- as.numeric(strsplit(trimws(strsplit(best_line, ":")[[1]][[2]]), " ")[[1]][[1]])
  
  # Extract bounds
  upper_line <- grep("^Upper bound:", lines, value = TRUE)
  upper_value <-as.numeric(strsplit(trimws(strsplit(upper_line, ":")[[1]][[2]]), " ")[[1]][[1]])
  
  lower_line <- grep("^Lower bound:", lines, value = TRUE)
  lower_value <- as.numeric(strsplit(trimws(strsplit(lower_line, ":")[[1]][[2]]), " ")[[1]][[1]])
  
  # Build output data frame
  result <- data.frame(
    param_names_fit = paste(param_names, collapse = ", "),
    extent_scaling_factor = extent_value,
    delta_meanlog = delta_meanlog,
    delta_meanlog_status = delta_meanlog_status,
    ratio_sdlog = ratio_sdlog,
    ratio_sdlog_status = ratio_sdlog_status,
    best_estimate = best_value,
    upper_bound = upper_value,
    lower_bound = lower_value,
    stringsAsFactors = FALSE
  )
  
  if (verbose) {
    message("Extracted: best estimate = ", best_value, " species")
  }}
  
  return(result)
}

# --- Helper: extract datasetsummary parameters from output -----------------------------------
extract_dataset_summary_file <- function(directory, verbose = TRUE) {
  #' Extract Dataset Summary Statistics
  #'
  #' Parses dataset_summary_*.txt files to extract filtering parameters
  #' and before/after data quality metrics.
  #'
  #' @param directory Character. Path to directory containing dataset_summary_*.txt files.
  #' @param verbose Logical. If TRUE, print status messages.
  #'
  #' @return Data frame (1 row) with filter thresholds (min_occurrences, min_completeness,
  #'         min_richness, outlier_cutoff, track_b_entry, filter_applied) and paired
  #'         before/after statistics for grid cells, total records, total species,
  #'         and richness distribution (min, Q1, median, mean, Q3, max).
  #'         Also includes species loss counts (species_retained, species_lost).
  #'
  #' @details
  #' Searches for the most recently modified file matching "dataset_summary_*.txt".
  #' Extracts both threshold parameters and summary statistics in structured format
  #' for downstream validation and reporting.
  
  files <- list.files(directory, 
                      pattern = "dataset_summary.*\\.txt$",
                      full.names = TRUE)
  
  if (length(files) == 0) {
    message("No dataset_summary files found in: ", directory)
    result <- data.frame(
      min_occurrences = NA,
      min_completeness = NA,
      min_richness = NA,
      outlier_cutoff = NA,
      track_b_entry = NA,
      filter_applied = NA,
      before_grid_cells = NA,
      before_total_records = NA,
      before_total_species = NA,
      before_richness_min = NA,
      before_richness_Q1 = NA,
      before_richness_median = NA,
      before_richness_mean = NA,
      before_richness_Q3 = NA,
      before_richness_max = NA,
      after_grid_cells = NA,
      after_total_records = NA,
      after_total_species = NA,
      after_richness_min = NA,
      after_richness_Q1 = NA,
      after_richness_median = NA,
      after_richness_mean = NA,
      after_richness_Q3 = NA,
      after_richness_max = NA,
      species_retained = NA,
      species_retained_total = NA,
      species_lost = NA,
      stringsAsFactors = FALSE
    )
  }else{
  
  if (length(files) > 1) {
    file_info <- file.info(files)
    target_file <- rownames(file_info)[which.max(file_info$mtime)]
  } else {
    target_file <- files[1]
  }
  
  if (verbose) {
    message("Reading: ", basename(target_file))
  }
  
  lines <- readLines(target_file)
  
  # Extract numeric value immediately after colon
  extract_value <- function(pattern, lines) {
    line <- grep(pattern, lines, value = TRUE)
    if (length(line) == 0) return(NA_real_)
    value_str <- trimws(strsplit(line[1], ":")[[1]][2])
    return(as.numeric(strsplit(value_str, " ")[[1]][1]))
  }
  
  # Extract text value immediately after colon
  extract_quoted <- function(pattern, lines) {
    line <- grep(pattern, lines, value = TRUE)
    if (length(line) == 0) return(NA_character_)
    return(trimws(strsplit(line[1], ":")[[1]][2]))
  }
  
  # Extract distribution statistics from comma-separated key=value pairs
  extract_stats <- function(line) {
    if (is.na(line) || length(line) == 0) {
      return(list(min = NA, Q1 = NA, med = NA, mean = NA, Q3 = NA, max = NA))
    }
    
    parts_str <- strsplit(line, ":")[[1]][2]
    parts <- strsplit(parts_str, ",")[[1]]
    
    stats <- list()
    for (part in parts) {
      part_clean <- trimws(part)
      if (grepl("=", part_clean)) {
        split_part <- strsplit(part_clean, "=")[[1]]
        stat_name <- trimws(split_part[1])
        stat_value <- as.numeric(trimws(split_part[2]))
        stats[[stat_name]] <- stat_value
      }
    }
    
    return(list(
      min = stats$min,
      Q1 = stats$Q1,
      med = stats$med,
      mean = stats$mean,
      Q3 = stats$Q3,
      max = stats$max
    ))
  }
  
  # Filter thresholds
  min_occurrences <- extract_value("^  min_occurrences:", lines)
  min_completeness <- extract_value("^  min_completeness:", lines)
  min_richness <- extract_value("^  min_richness:", lines)
  outlier_cutoff <- extract_quoted("^  outlier cutoff:", lines)
  track_b_entry <- extract_quoted("^  Track B entry:", lines)
  filter_applied <- extract_quoted("^  Filter applied:", lines)
  
  # Before filtering stats
  before_grid_cells <- extract_value("^  Grid cells:.*$", lines)
  before_records <- extract_value("^  Total records:.*$", lines)
  before_species <- extract_value("^  Total species:.*$", lines)
  
  richness_lines <- grep("^  Richness/cell:", lines, value = TRUE)
  before_richness <- extract_stats(richness_lines[1])
  after_richness <- extract_stats(richness_lines[length(richness_lines)])
  
  # After filtering stats
  after_cells_line <- lines[grep("^  Grid cells:.*well-sampled", lines)]
  after_grid_cells <- if (length(after_cells_line) > 0) {
    as.numeric(strsplit(trimws(strsplit(after_cells_line[1], ":")[[1]][2]), " ")[[1]][1])
  } else {
    NA_real_
  }
  
  after_records_line <- lines[grep("^  Total records:.*\\(.*%", lines)]
  after_records <- if (length(after_records_line) > 0) {
    as.numeric(strsplit(trimws(strsplit(after_records_line[1], ":")[[1]][2]), " ")[[1]][1])
  } else {
    NA_real_
  }
  
  after_species_line <- lines[grep("^  Total species:.*\\(.*%", lines)]
  after_species <- if (length(after_species_line) > 0) {
    as.numeric(strsplit(trimws(strsplit(after_species_line[1], ":")[[1]][2]), " ")[[1]][1])
  } else {
    NA_real_
  }
  
  # Species loss
  species_retained_line <- grep("^  Species retained:", lines, value = TRUE)
  if (length(species_retained_line) > 0) {
    species_retained_str <- trimws(strsplit(species_retained_line[1], ":")[[1]][2])
    species_retained <- as.numeric(strsplit(species_retained_str, "/")[[1]][1])
    species_retained_total <- as.numeric(strsplit(species_retained_str, "/")[[1]][2])
  } else {
    species_retained <- NA_real_
    species_retained_total <- NA_real_
  }
  
  species_lost_line <- grep("^  Species lost:", lines, value = TRUE)
  if (length(species_lost_line) > 0) {
    species_lost_str <- trimws(strsplit(species_lost_line[1], ":")[[1]][2])
    species_lost <- as.numeric(strsplit(species_lost_str, " ")[[1]][1])
  } else {
    species_lost <- NA_real_
  }
  
  result <- data.frame(
    min_occurrences = min_occurrences,
    min_completeness = min_completeness,
    min_richness = min_richness,
    outlier_cutoff = outlier_cutoff,
    track_b_entry = track_b_entry,
    filter_applied = filter_applied,
    before_grid_cells = before_grid_cells,
    before_total_records = before_records,
    before_total_species = before_species,
    before_richness_min = before_richness$min,
    before_richness_Q1 = before_richness$Q1,
    before_richness_median = before_richness$med,
    before_richness_mean = before_richness$mean,
    before_richness_Q3 = before_richness$Q3,
    before_richness_max = before_richness$max,
    after_grid_cells = after_grid_cells,
    after_total_records = after_records,
    after_total_species = after_species,
    after_richness_min = after_richness$min,
    after_richness_Q1 = after_richness$Q1,
    after_richness_median = after_richness$med,
    after_richness_mean = after_richness$mean,
    after_richness_Q3 = after_richness$Q3,
    after_richness_max = after_richness$max,
    species_retained = species_retained,
    species_retained_total = species_retained_total,
    species_lost = species_lost,
    stringsAsFactors = FALSE
  )
  
  if (verbose) {
    message("Extracted: ", after_species, " species retained from ", 
            before_species, " (", species_lost, " lost)")
  }}
  
  return(result)
}

message("Function 'extract_dataset_summary_file' is ready to use.")




# --- check_status() ----------------------------------------------------------

check_status <- function(organism = NULL) {
  
  # Shared upstream stages (always the same regardless of organism)
  shared <- list(
    "Stage 0: Plant curation" = c(
      "species_counts3.rds"   = file.path(stage_dirs$curate, "output", "species_counts3.rds"),
      "species_endemism2.rds" = file.path(stage_dirs$curate, "output", "species_endemism2.rds")
    ),
    "Stage 1: Plant filtering" = c(
      "species_counts4_filtered.rds"       = file.path(stage_dirs$filter, "output", "species_counts4_filtered.rds"),
      "FINAL_allEndemism_FullyFiltered.rds" = file.path(stage_dirs$filter, "output", "FINAL_allEndemism_FullyFiltered.rds")
    ),
    "Stage 3: Deforestation patches" = c(
      "freq: boundary_filtered" = file.path(stage_dirs$predict, "output", "patch_area_frequency_logbins_boundary_filtered.csv"),
      "freq: overfiltered"      = file.path(stage_dirs$predict, "output", "patch_area_frequency_logbins_overfiltered.csv"),
      "freq: underfiltered"     = file.path(stage_dirs$predict, "output", "patch_area_frequency_logbins_underfiltered.csv"),
      "model: boundary_filtered" = file.path(stage_dirs$predict, "output", "lognormal_model_parameters_boundary_filtered.csv"),
      "model: overfiltered"     = file.path(stage_dirs$predict, "output", "lognormal_model_parameters_overfiltered.csv"),
      "model: underfiltered"    = file.path(stage_dirs$predict, "output", "lognormal_model_parameters_underfiltered.csv"),
      "patch metrics directory" = file.path(stage_dirs$predict, "output", "deforestTileOutStats")
    )
  )
  
  for (stage_name in names(shared)) {
    paths <- shared[[stage_name]]
    status <- ifelse(file.exists(paths), "\u2713", "\u2717")
    message("\n", stage_name)
    for (i in seq_along(paths)) message("  ", status[i], " ", names(paths)[i])
  }
  
  # Per-organism stages
  if (is.null(organism)) {
    # List all organisms found
    organisms <- list.dirs(output_root, recursive = FALSE, full.names = FALSE)
    if (length(organisms) == 0) {
      message("\nNo organism output directories found in: ", output_root)
    } else {
      message("\nOrganism runs found: ", paste(organisms, collapse = ", "))
      for (org in organisms) {
        dirs <- organism_dirs(org)
        has_animal <- file.exists(file.path(dirs$animal_model, "simulatedOccurencesForAnimal.rds"))
        has_ext    <- length(list.files(dirs$extinction_model, pattern = "_all_replicates\\.csv$")) > 0
        message("  ", org, ": animal_model=",
                ifelse(has_animal, "\u2713", "\u2717"),
                "  extinction_model=",
                ifelse(has_ext, "\u2713", "\u2717"))
      }
    }
  } else {
    dirs <- organism_dirs(organism)
    message("\n", organism, " — Stage 2: Animal modelling")
    am_path <- dirs$animal_model
    if (dir.exists(am_path)) {
      rds <- list.files(am_path, pattern = "\\.rds$")
      pdfs <- list.files(am_path, pattern = "\\.pdf$")
      csvs <- list.files(am_path, pattern = "\\.csv$")
      has_sim <- "simulatedOccurencesForAnimal.rds" %in% rds
      message("  ", ifelse(has_sim, "\u2713", "\u2717"),
              " simulatedOccurencesForAnimal.rds")
      message("  ", length(pdfs), " PDF plots, ", length(csvs), " CSVs, ",
              length(rds), " RDS files")
    } else {
      message("  \u2717 directory does not exist")
    }
    
    message("\n", organism, " — Stage 4: Extinction modelling")
    em_path <- dirs$extinction_model
    if (dir.exists(em_path)) {
      for (model in c("neutral", "clustered", "shifted")) {
        f <- file.path(em_path, paste0(model, "_all_replicates.csv"))
        message("  ", ifelse(file.exists(f), "\u2713", "\u2717"),
                " ", model, " replicates")
      }
      summ <- file.path(em_path, "summary_all_models.csv")
      message("  ", ifelse(file.exists(summ), "\u2713", "\u2717"),
              " summary table")
    } else {
      message("  \u2717 directory does not exist")
    }
  }
  
  invisible(NULL)
}


# --- run_animal_model() ------------------------------------------------------
# Runs the Stage 2 animal modelling pipeline for a given organism.
#
# Parameters:
#   organism             - name for the output directory under output/.
#   data_file            - path to input data (.csv, .zip, .tsv, or .xlsx).
#   secondary_data_file  - optional second input file to combine with data_file.
#   simulateParallel     - logical; run simulations in parallel.
#   simulateNCores       - number of cores to use when simulateParallel = TRUE.
#   fill_plants          - logical; use interpolated (TRUE) or empirical (FALSE) plant occurrences.
#   fill_animals         - logical; buffer/interpolate animal occurrences (TRUE) or use empirical only (FALSE; recommended). 
#   useSamplingExtension - logical; whether to apply a the linear model extension based on sampling parameters (Default = TRUE;recommended)
#   simulate_ranges      - TRUE, "BestOnly", or "SkipUpper"; controls how many animal geographic range simulations are run. Set to "BestOnly" if you don't want to simulate three different sets of animal ranges or "SkipUpper" which just makes a best simulation and a conservative simulation (low estimate). TRUE (default) simulates ranges for the best estimate and the two 95% confidence estimates.
#   clean                - logical; run 00_cleanGBIF.R before modelling (.csv/.zip only) to clean GBIF data that is potentially messy
#   force_clean          - logical; re-run cleaning even if a cleaned file already exists.
#   fix_taxon_names      - logical; resolve taxon name spelling/synonyms via GBIF backbone (requires clean = TRUE).
#   taxon_kingdom        - kingdom string used to restrict taxon name matching, e.g. "Animalia".
#   taxon_review_file    - optional path to write a review CSV of unresolved/spurious name matches.
#   redo                 - logical; recompute raw endemism and plant richness from scratch. 
#   alphaRichnessMethod  - method for alpha richness estimation; default "empirical". NOT CURRENTLY IN USE
#   
#   ANIMAL-CELL FILTERING PARAMETERS
#   animalFilterType     - "balanced", "strict", or "extraStrict"; controls how strictly animal occurrences/cells are filtered. "balanced (Default) will not fully respect your minimum parameters but will make sure to include sites above those minima if another filter would otherwise remove it. "strict" will only filter based on completeness (but will automatically turn off if resulting in <30 input cells). "extraStrict" does the same as "strict" but may over-filter sites resulting in little data to analyze (ignore 30 cell rule)
#   completenessMethod   - Method for calculating completeness per grid-cell: "Chao1", "ACE", "Coverage", c("Chao1", "ACE"), or "all". The multi parameter input and "all" will keep cells passed by all the filters.
#   min_occurrences      - minimum occurrences per cell for "well-sampled" cells; "auto" or a number.
#   min_completeness     - minimum Chao1 completeness for well-sampled cells; "auto" or a number. Good values tend to be between 0.6 and 0.8.
#   min_richness         - minimum observed richness per cell to keep; default 20.
#   max_richness         - maximum plausible richness per cell; NULL uses the 99.9th percentile.
#   disable_outlier_filter - logical; if TRUE, skip the high-richness outlier cutoff.
#   
#   filter_plant_completeness - numeric threshold to restrict plant cells by Chao1 completeness, or NULL.
#   fit_params           - character vector of log-normal parameters to fit, subset of c("amplitude", "meanlog", "sdlog").
#   estimateBounds       - logical; whether to estimate upper/lower bounds. Set to FALSE for a faster run (without bounds based on simulated species at least)
#   gammaErrorType       - "full" or "narrow"; controls gamma-diversity error/bounds reporting. "full" reports the accumulated error from both models. Simulations were tested on "narrow" errors and they prove to be relatively accurate.
#   skipRawEndemism      - logical; skip raw endemism computation, which isn't used in the full model. Recommend to skip for time.
#   customOutputFolderPrefix - optional string prefixed to the "output" segment of the output path.
#   customOutputFolderSuffix - optional string appended to the organism output folder name.


# Case- and whitespace-insensitive matching of character arguments to their
# canonical spellings. `x` may be a vector or a comma-separated string.
# Returns canonical values; errors on anything unrecognised.
.amani_canonical <- function(x, choices, arg) {
  vals <- trimws(unlist(strsplit(as.character(x), ",")))
  idx <- match(tolower(vals), tolower(choices))
  if (anyNA(idx)) {
    stop(sprintf("Invalid value(s) for '%s': %s. Valid options: %s",
                 arg, paste(shQuote(vals[is.na(idx)]), collapse = ", "),
                 paste(shQuote(choices), collapse = ", ")), call. = FALSE)
  }
  unique(choices[idx])
}

run_animal_model <- function(organism,
                             data_file,
                             secondary_data_file = NULL,
                             fill_plants = TRUE,
                             fill_animals = FALSE,
                             alphaRichnessMethod = "empirical",
                             useSamplingExtension = TRUE,
                             simulate_ranges = TRUE, 
                             clean = FALSE,
                             force_clean = FALSE,
                             fix_taxon_names = FALSE,
                             taxon_kingdom = "Animalia",
                             taxon_review_file = NULL,
                             redo = TRUE,
                             min_occurrences = "auto",
                             min_completeness = "auto",
                             min_richness = 20,
                             max_richness = NULL,
                             disable_outlier_filter = FALSE,
                             filter_plant_completeness = NULL,
                             fit_params = c("amplitude"),
                             estimateBounds = TRUE,
                             gammaErrorType = c("full", "narrow"),
                             skipRawEndemism = TRUE,
                             animalFilterType = c("balanced", "strict", "extraStrict"),
                             completenessMethod = "Chao1",
                             customOutputFolderPrefix = NULL,
                             customOutputFolderSuffix = NULL,
                             simulateParallel=FALSE,
                             simulateNCores=1
                             ) {
  # Run the animal modelling pipeline (Stage 2).
  #
  # The entry-point script is chosen automatically based on data_file extension:
  #   .xlsx      -> scripts/01_CockroachAnalysis.R  (Darwin Core Excel / original format)
  #   .csv/.zip  -> scripts/01_GBIFAnalysis.R       (standard GBIF CSV/zip download)
  #
  # Args:
  #   organism:   name used for the output directory under output/
  #   data_file:  path to the input data file (relative to 2-ModellingAnimals/).
  #               Accepts .csv, .zip (GBIF zip download), or .xlsx.
  #   secondary_data_file: optional second input file (relative path, same as
  #               data_file). Default NULL preserves current behaviour. When set,
  #               data_file and secondary_data_file are read and row-bound by
  #               matching column names (missing columns filled with NA) into a
  #               single combined CSV, and the rest of the pipeline runs on that
  #               combined file. Not supported for .xlsx inputs. Either file's
  #               binomial name column may be labelled 'species' or
  #               'scientificName'; both are folded into a single 'species'
  #               column and 'scientificName' is dropped from the output.
  #   fill_plants:if TRUE, uses the plant occurrence data with interpolated range
  #               occurrences. If FALSE (default) only the empirical occurrences
  #               are used
  #   fill_animals:if TRUE, buffers and interpolates the animal occurrence data.
  #               If FALSE (default) only the empirical occurrences are used.            
  #   clean:      if TRUE, run 00_cleanGBIF.R on the file before modelling.
  #               Applies to .csv and .zip files; ignored for .xlsx.
  #   force_clean: if TRUE, re-run cleaning even if a *_cleaned.csv already
  #               exists. Useful when cleaning code has changed.
  #   fix_taxon_names: if TRUE, run CleanTaxonNames.R after the GBIF cleaning
  #               step to resolve spelling variants and synonyms via the GBIF
  #               backbone. Requires clean = TRUE (or a pre-cleaned input). The
  #               corrected data are saved as *_taxaFixed.csv and used for all
  #               downstream modelling. Default FALSE.
  #   taxon_kingdom: kingdom string passed to clean_taxon_names() to restrict
  #               fuzzy matching (e.g. "Animalia", "Plantae"). Default "Animalia".
  #   taxon_review_file: optional path for the review CSV of unresolved or
  #               spurious name matches (relative to 2-ModellingAnimals/). NULL
  #               (default) writes nothing.
  #   redo:       if TRUE (default), recompute raw endemism and plant richness
  #               from scratch. Set to FALSE to reuse cached intermediate objects.
  #   min_occurrences: minimum number of occurrences per cell to be considered
  #               "well-sampled" (Track A). Set to "auto" (default) to estimate
  #               from the species discovery rate curve (the occurrence count
  #               where the rate drops below 0.25 new species per occurrence).
  #               Set to a positive number to override. Set to 0 to disable
  #               cell filtering entirely.
  #   min_completeness: minimum Chao1 sampling completeness (S_obs / S_est)
  #               for Track A cells. Set to "auto" (default) to estimate from
  #               the data (the completeness at which the median estimated
  #               missing fraction drops below 50%). Set to a numeric value to
  #               override.
  #   min_richness: minimum observed species richness for Track A cells.
  #               Default 20. This should be set manually as a very conservative
  #               estimate of how many species of this organism group one would
  #               expect in the least diverse region of the Amazon, including
  #               deforested areas.
  #   max_richness: maximum plausible species richness per grid cell. If NULL
  #               (default), the 99.9th percentile of observed richness is used
  #               as the outlier cutoff. Set to a positive number to use a
  #               domain-informed upper bound instead (e.g., the maximum number
  #               of species of this group that could realistically co-occur in
  #               a ~22km grid cell).
  #   disable_outlier_filter: if TRUE, skip the outlier cutoff entirely — no cells
  #               are dropped for having too-high richness. Useful when the high-
  #               richness tail reflects genuine diversity rather than data artefacts.
  #               Default FALSE.
  #   filter_plant_completeness: if set to a numeric value (e.g. 0.60),
  #               restrict the plant-animal correlation analyses to plant
  #               grid cells with Chao1 completeness >= this value. NULL
  #               (default) uses all plant cells.
  #   fit_params: character vector — which log-normal parameters to estimate
  #               from the animal MDD data when scaling the plant spatial
  #               extent distribution. Subset of c("amplitude", "meanlog", "sdlog").
  #                 "amplitude" — scale total species count (default; existing behaviour)
  #                 "meanlog"   — shift the distribution's location (range sizes)
  #                 "sdlog"     — rescale the distribution's spread (endemic vs widespread)
  #               Parameters not listed are inherited from the plant fit.
  #               Default: c("amplitude").
  
  if (missing(organism)) stop("Specify an organism name, e.g. run_animal_model('cockroaches')", call. = FALSE)
  if (missing(data_file)) stop("Specify a data_file path.", call. = FALSE)

  # Which gamma-diversity error mode to report: "full" (upper + lower bound
  # re-simulation) or "narrow" (conservative bounds only). Passed to
  # 02_animalAnalysis.R via .amani_gammaErrorType.
  if (length(gammaErrorType) > 1) gammaErrorType <- gammaErrorType[1]
  if (length(animalFilterType) > 1) animalFilterType <- animalFilterType[1]
  gammaErrorType <- .amani_canonical(gammaErrorType, c("full", "narrow"), "gammaErrorType")
  animalFilterType <- .amani_canonical(animalFilterType, c("balanced", "strict", "extraStrict"), "animalFilterType")

  # simulate_ranges: logical, or "BestOnly"/"SkipUpper" in any capitalisation
  if (!is.logical(simulate_ranges)) {
    simulate_ranges <- .amani_canonical(simulate_ranges, c("BestOnly", "SkipUpper", "TRUE", "FALSE"), "simulate_ranges")
    if (length(simulate_ranges) != 1) stop("simulate_ranges must be a single value.", call. = FALSE)
    if (simulate_ranges %in% c("TRUE", "FALSE")) simulate_ranges <- as.logical(simulate_ranges)
  }

  # Completeness methods: "Chao1", "ACE", "Coverage", "all", or combinations
  completenessMethod <- .amani_canonical(completenessMethod, c("Chao1", "ACE", "Coverage", "all"), "completenessMethod")
  if ("all" %in% completenessMethod) completenessMethod <- "all"
  if (length(completenessMethod) == 2 && all(c("Chao1", "ACE") %in% completenessMethod)) completenessMethod <- "Chao1, ACE"

  fit_params <- .amani_canonical(fit_params, c("amplitude", "meanlog", "sdlog"), "fit_params")

  # Strict filter modes rely on explicit thresholds; fill any unset ones with defaults.
  if (animalFilterType %in% c("strict", "extraStrict")) {
    is_unset <- function(arg_missing, val) {
      arg_missing || is.null(val) || (length(val) == 1 && (is.na(val) || identical(val, "auto")))
    }
    strict_defaults <- list(min_occurrences = 5, min_richness = 25, min_completeness = 0.7)
    unset <- c(
      min_occurrences  = is_unset(missing(min_occurrences),  min_occurrences),
      min_richness     = is_unset(missing(min_richness),     min_richness),
      min_completeness = is_unset(missing(min_completeness), min_completeness)
    )
    if (any(unset)) {
      to_fill <- names(unset)[unset]
      warning(sprintf(
        "animalFilterType = '%s' requires explicit minimum thresholds, but these were not set: %s. Using defaults: %s. These defaults may not be appropriate for '%s'; refine them manually for this organism group.",
        animalFilterType,
        paste(to_fill, collapse = ", "),
        paste(sprintf("%s = %s", to_fill, unlist(strict_defaults[to_fill])), collapse = ", "),
        organism
      ), call. = FALSE)
      for (nm in to_fill) assign(nm, strict_defaults[[nm]])
    }
  }
  
  # Choose entry script based on file extension
  ext <- tolower(tools::file_ext(data_file))
  script <- switch(ext,
                   xlsx = "scripts/01_CockroachAnalysis.R",
                   csv  = "scripts/01_GBIFAnalysis.R",
                   tsv  = "scripts/01_GBIFAnalysis.R",
                   zip  = "scripts/01_GBIFAnalysis.R",
                   stop("Unsupported data file format: .", ext,
                        "\n  Supported: .xlsx (Darwin Core Excel), .tsv, .csv or .zip (GBIF download)",
                        call. = FALSE)
  )
  
  # Optionally combine a secondary data file before anything else.
  #   Combining is only supported for the GBIF CSV/TSV/zip path. Files are
  #   read and row-bound by matching column names (fill = TRUE), then written
  #   to a single combined CSV that the rest of the pipeline processes.
  if (!is.null(secondary_data_file)) {
    sec_ext <- tolower(tools::file_ext(secondary_data_file))
    if (ext == "xlsx" || sec_ext == "xlsx") {
      stop("secondary_data_file is not supported for XLSX (Darwin Core Excel) input.",
           call. = FALSE)
    }

    # Reader that mirrors the script's fread_auto (handles zip + tsv).
    #   Also normalises the species column at read time: files may label
    #   the binomial column 'species' (standard GBIF downloads) or
    #   'scientificName' (e.g. curated files). Whichever header is present
    #   is treated as the species column, so rbindlist() below aligns the
    #   two files on a single 'species' column instead of leaving one file's
    #   values stranded under a differently-named column. 'scientificName' is
    #   dropped afterward so the combined output carries only 'species'.
    read_any <- function(rel_path) {
      p <- file.path(stage_dirs$animals, rel_path)
      e <- tolower(tools::file_ext(rel_path))
      dt <- if (e == "zip") {
        inner <- grep("\\.(csv|tsv)$", unzip(p, list = TRUE)$Name,
                      ignore.case = TRUE, value = TRUE)[1]
        if (is.na(inner)) stop("No CSV/TSV found inside: ", rel_path, call. = FALSE)
        inner_sep <- if (tolower(tools::file_ext(inner)) == "tsv") "\t" else ","
        cmd <- sprintf('unzip -p "%s" "%s"', normalizePath(p, winslash = "/"), inner)
        data.table::fread(cmd = cmd, sep = inner_sep)
      } else {
        data.table::fread(file = p, sep = if (e == "tsv") "\t" else ",")
      }

      # Flexible species-column matching: accept either 'species' or
      # 'scientificName' as the header, and standardise to 'species'.
      has_species <- "species" %in% names(dt)
      has_sciname <- "scientificName" %in% names(dt)
      if (!has_species && has_sciname) {
        dt[, species := scientificName]
        message("  ", basename(rel_path), ": using 'scientificName' as the species column.")
      } else if (has_species && has_sciname) {
        n_fixed <- dt[is.na(species) | species == "", .N]
        if (n_fixed > 0) {
          dt[is.na(species) | species == "", species := scientificName]
          message("  ", basename(rel_path), ": filled ", n_fixed,
                  " missing 'species' values from 'scientificName'.")
        }
      }

      # Drop 'scientificName' now that its values have been folded into
      # 'species' -- the combined file should carry a single species column.
      if (has_sciname) dt[, scientificName := NULL]

      dt
    }

    message("  Combining primary + secondary data files ...")
    combined <- data.table::rbindlist(
      list(read_any(data_file), read_any(secondary_data_file)),
      use.names = TRUE, fill = TRUE
    )

    combined_rel <- file.path("data",organism, sprintf("combined_%s.csv", organism))
    combined_abs <- file.path(stage_dirs$animals, combined_rel)
    if (!dir.exists(dirname(combined_abs))) dir.create(dirname(combined_abs), recursive = TRUE)
    data.table::fwrite(combined, combined_abs)

    message("  Wrote combined file (", nrow(combined), " rows): ", combined_rel)

    # Downstream reads a single combined CSV.
    data_file <- combined_rel
    ext <- "csv"
    script <- "scripts/01_GBIFAnalysis.R"
  }

  # Optionally clean GBIF CSV/zip before modelling
  #   Skip if the file is already a cleaned output (avoids re-cleaning when
  #   the tool sandbox or a retry re-enters this block).
  if (clean && ext %in% c("csv", "zip") &&
      (!grepl("_cleaned\\.csv$", data_file) || force_clean)) {
    if (force_clean) message("  force_clean = TRUE — re-cleaning regardless of existing output.")
    raw_path <- file.path(stage_dirs$animals, data_file)
    cleaned_path <- clean_gbif(raw_path)
    # Update data_file and ext to point to the cleaned CSV
    data_file <- sub(paste0("^", stage_dirs$animals, "/?"), "", cleaned_path)
    ext <- "csv"
    message("  Using cleaned file: ", data_file)
  } else if (grepl("_cleaned\\.csv$", data_file) && !force_clean) {
    message("  Input already cleaned — skipping clean_gbif.")
  }

  # Optionally standardise taxon names via the GBIF backbone
  #   Runs after clean_gbif so it operates on standardised, filtered records.
  #   Saves a *_taxaFixed.csv alongside the cleaned file and updates data_file.
  #   Requires a cleaned input: either clean = TRUE was set (so clean_gbif() has
  #   already run), or the data_file already carries a _cleaned.csv suffix.
  #   Running on raw combined data is unsafe because curated rows often lack a
  #   populated `species` column (names stored in `scientificName` instead).
  if (fix_taxon_names && ext == "csv") {

    source(file.path(stage_dirs$animals, "scripts/utils/CleanTaxonNames.R"),
           local = TRUE)

    abs_cleaned <- file.path(stage_dirs$animals, data_file)
    taxa_fixed_path <- sub("(\\.csv)?$", "_taxaFixed.csv",
                           abs_cleaned, ignore.case = TRUE)

    review_abs <- if (!is.null(taxon_review_file))
      file.path(stage_dirs$animals, taxon_review_file)
    else
      NULL

    message("  Running clean_taxon_names() — querying GBIF backbone ...")
    clean_taxon_names(
      occ          = abs_cleaned,
      kingdom      = taxon_kingdom,
      output_file  = taxa_fixed_path,
      review_file  = review_abs
    )

    data_file <- sub(paste0("^", stage_dirs$animals, "/?"), "", taxa_fixed_path)
    message("  Using taxa-fixed file: ", data_file)
  }
  
  # Verify upstream plant data exists
  upstream_files <- c(
    file.path(stage_dirs$filter, "output", "species_counts4_filtered.rds"),
    file.path(stage_dirs$filter, "output", "FINAL_allEndemism_FullyFiltered.rds"),
    file.path(stage_dirs$filter, "output", "species_c4_filt_interpolated.rds"),
    file.path(stage_dirs$filter, "output", "INTERP_allEndemism_FullyFiltered.rds")
    
  )
  verify_upstream("Stage 2 (animal modelling)", upstream_files)
  
  # Set up organism-specific output directory
  dirs <- organism_dirs(organism)
  if (!is.null(customOutputFolderPrefix)) {
    dirs$animal_model <- sub(
      "(^|[/\\\\])output([/\\\\])",
      paste0("\\1", customOutputFolderPrefix, "output\\2"),
      dirs$animal_model
    )
  }
  if (!is.null(customOutputFolderSuffix)) {
    # Append the suffix to the organism directory (parent of animal_model),
    # e.g. p_output/ORGANISM -> p_output/ORGANISM_suffix
    org_dir <- dirname(dirs$animal_model)
    dirs$animal_model <- file.path(
      file.path(dirname(org_dir), paste0(basename(org_dir), customOutputFolderSuffix)),
      basename(dirs$animal_model)
    )
  }
  
  
  
  
  if (!dir.exists(dirs$animal_model)) dir.create(dirs$animal_model, recursive = TRUE)
  archive_outputs(dirs$animal_model)
  
  # Set global overrides for sourced scripts
  set_stage_wd(stage_dirs$animals)
  .amani_output_dir <<- dirs$animal_model
  .amani_simulateParallel<<-simulateParallel
  .amani_simulateNCores<<-simulateNCores
  .amani_redo <<- redo
  .amani_completenessMethod <<- completenessMethod
  .amani_min_occurrences <<- min_occurrences
  .amani_min_completeness <<- min_completeness
  .amani_min_richness <<- min_richness
  .amani_max_richness <<- max_richness
  .amani_disable_outlier_filter <<- disable_outlier_filter
  .amani_filter_plant_completeness <<- filter_plant_completeness
  .amani_fit_params <<- fit_params
  .amani_alphaRichnessMethod <<- alphaRichnessMethod
  .amani_data_file <<- data_file
  .amani_estimateBounds <<- estimateBounds
  .amani_gammaErrorType <<- gammaErrorType
  .amani_skipRawEndemism<<- skipRawEndemism
  .amani_animalFilterType<<-animalFilterType
  .amani_fillAnimals<<-fill_animals
  .amani_simulateRanges<<-simulate_ranges
  .amani_useSamplingExtension<<-useSamplingExtension
  
  # Expose fill_plants so the sourced scripts can key the plant-richness cache
  # by plant-data type (filled vs unfilled).
  .amani_fill_plants<<-fill_plants
  #choose the plant occurrence data
  dataDir <- file.path(stage_dirs$filter, "output")
  if (fill_plants){
    occ_sf<<-readRDS(paste0(dataDir,"/species_c4_filt_interpolated.rds"))
    endemismData<<-readRDS(paste0(dataDir,"/INTERP_allEndemism_FullyFiltered.rds"))
    
    }else{
    occ_sf<<-readRDS(paste0(dataDir,"/species_counts4_filtered.rds"))
    endemismData<<-readRDS(paste0(dataDir,"/FINAL_allEndemism_FullyFiltered.rds"))
    
  }
  
  
  on.exit(restore_wd())
  
  message("Running animal model for: ", organism)
  message("  Data:   ", data_file, " (format: ", ext, ")")
  message("  Script: ", script)
  message("  Output: ", dirs$animal_model)
  message("  Fit params: [", paste(fit_params, collapse = ", "), "]")
  message("Using plant data type: ", if(fill_plants){"range interpolated (filled)"}else{"raw empirical"})
  
  source(script, local = FALSE)
  
  message("\nAnimal modelling complete.")
}


# --- run_extinction_model() ---------------------------------------------------

run_extinction_model <- function(organism,
                                 models = c("neutral", "clustered", "shifted"),
                                 summarize = TRUE) {
  # Run the extinction simulation pipeline (Stage 4).
  #
  # Reads pre-existing data from Stages 1, 2, and 3 and writes results to
  # output/{organism}/extinction_model/.
  #
  # Args:
  #   organism:  name of the organism (must match a directory under output/)
  #   models:    character vector — subset of c("neutral", "clustered", "shifted")
  #   summarize: if TRUE (default), run 03_summarize_results.R afterward

  if (missing(organism)) stop("Specify an organism name, e.g. run_extinction_model('cockroaches')", call. = FALSE)

  models <- match.arg(models, c("neutral", "clustered", "shifted"), several.ok = TRUE)
  dirs <- organism_dirs(organism)

  # Verify upstream data from Stages 1, 2, and 3
  # Stage 4 uses three filtering regimes from Stage 3 (overfiltered, boundary_filtered, underfiltered)
  upstream_files <- c(
    file.path(stage_dirs$filter, "output", "species_counts4_filtered.rds"),
    file.path(dirs$animal_model, "simulatedOccurencesForAnimal.rds"),
    file.path(stage_dirs$predict, "output", "patch_area_frequency_logbins_boundary_filtered.csv"),
    file.path(stage_dirs$predict, "output", "patch_area_frequency_logbins_overfiltered.csv"),
    file.path(stage_dirs$predict, "output", "patch_area_frequency_logbins_underfiltered.csv"),
    file.path(stage_dirs$predict, "output", "lognormal_model_parameters_boundary_filtered.csv"),
    file.path(stage_dirs$predict, "output", "lognormal_model_parameters_overfiltered.csv"),
    file.path(stage_dirs$predict, "output", "lognormal_model_parameters_underfiltered.csv")
  )
  verify_upstream("Stage 4 (extinction modelling)", upstream_files)

  # Set up organism-specific output directory
  if (!dir.exists(dirs$extinction_model)) dir.create(dirs$extinction_model, recursive = TRUE)
  archive_outputs(dirs$extinction_model)

  # Set global overrides for sourced scripts
  set_stage_wd(stage_dirs$model)
  .amani_output_dir <<- dirs$extinction_model
  .amani_animal_model_dir <<- dirs$animal_model
  on.exit(restore_wd())

  message("Running extinction model for: ", organism)
  message("  Animal model input: ", dirs$animal_model)
  message("  Output: ", dirs$extinction_model)

  # 02_simulate_clustered.R runs both clustered and shifted models
  scripts_to_run <- character(0)
  if ("neutral" %in% models) {
    scripts_to_run <- c(scripts_to_run, "scripts/01_simulate_neutral.R")
  }
  if ("clustered" %in% models || "shifted" %in% models) {
    scripts_to_run <- c(scripts_to_run, "scripts/02_simulate_clustered.R")
  }

  for (script in scripts_to_run) {
    message("\nRunning: ", script)
    source(script, local = FALSE)
  }

  if (summarize) {
    message("\nSummarizing results...")
    source("scripts/03_summarize_results.R", local = FALSE)
  }

  # IUCN threat classification (reads the summaries written above)
  classify_iucn_organism(organism, extinction_dir = dirs$extinction_model)

  message("\nExtinction modelling complete.")
}


# --- run_plant_model() --------------------------------------------------------
#
# Run the animal model pipeline using plant data as the "animal" input.
# This serves as a verification step: because the spatial extent model is
# built from plant data, applying the pipeline to plants should produce a
# gamma richness estimate close to the known plant richness.
#
# The function creates a pseudo-animal grid from the filtered plant species
# counts, then feeds it through the same Stage 2 + Stage 4 pipeline.

run_plant_model <- function(redo = TRUE,
                            min_occurrences = 0,
                            min_completeness = 0,
                            min_richness = 1,
                            filter_plant_completeness = NULL,
                            run_extinction = TRUE) {

  organism <- "plants_verification"

  # Verify upstream plant data exists
  plant_counts_file <- file.path(stage_dirs$filter, "output",
                                  "species_counts4_filtered.rds")
  plant_endemism_file <- file.path(stage_dirs$filter, "output",
                                    "FINAL_allEndemism_FullyFiltered.rds")
  verify_upstream("Plant verification", c(plant_counts_file, plant_endemism_file))

  message("=== Running pipeline on PLANTS (verification) ===")
  message("  This uses plant species as the 'animal' input to verify the pipeline.")

  # Load plant data and reshape it to look like GBIF animal data
  occ_sf <- readRDS(plant_counts_file)

  # Build a pseudo animal_grid: one row per (species, grid_cell)
  plant_as_animal <- occ_sf |>
    sf::st_drop_geometry() |>
    dplyr::distinct(scientificName, grid_id)

  # Build the grid reference from the same plant data
  grid_sf <- occ_sf |>
    dplyr::group_by(grid_id) |>
    dplyr::slice(1) |>
    dplyr::ungroup() |>
    dplyr::select(grid_id, geometry)

  # Compute richness per cell (matching the animal pipeline's expected format)
  animal_grid_data <- grid_sf |>
    dplyr::left_join(
      plant_as_animal |>
        dplyr::count(grid_id, name = "animal_richness"),
      by = "grid_id"
    ) |>
    dplyr::mutate(animal_richness = tidyr::replace_na(animal_richness, 0L))

  # Also create 'animal_grid' with Morphotype column (needed by ModelAnimal)
  animal_grid_global <- occ_sf |>
    dplyr::select(scientificName, grid_id, geometry) |>
    dplyr::rename(Morphotype = scientificName)

  # Set up output directory
  dirs <- organism_dirs(organism)
  if (!dir.exists(dirs$animal_model)) dir.create(dirs$animal_model, recursive = TRUE)
  archive_outputs(dirs$animal_model)

  # Push objects and parameters to global env for sourced scripts
  set_stage_wd(stage_dirs$animals)
  .amani_output_dir <<- dirs$animal_model
  .amani_redo <<- redo
  .amani_completenessMethod <<- completenessMethod
  .amani_min_occurrences <<- min_occurrences
  .amani_min_completeness <<- min_completeness
  .amani_min_richness <<- min_richness
  .amani_disable_outlier_filter <<- FALSE
  .amani_filter_plant_completeness <<- filter_plant_completeness
  .amani_simulateParallel<<-simulateParallel
  .amani_simulateNCores<<-simulateNCores
  # Plant verification always uses amplitude-only (self-consistency check)
  .amani_fit_params <<- c("amplitude")
  .amani_alphaRichnessMethod<<-alphaRichnessMethod
  .amani_fillAnimals<<-fill_animals

  # Push the pre-built data frames so the pipeline uses them directly
  animal_grid_data <<- animal_grid_data
  animal_grid <<- animal_grid_global

  on.exit(restore_wd())

  message("  Plant species as pseudo-animals: ", dplyr::n_distinct(plant_as_animal$scientificName))
  message("  Grid cells with data: ", sum(animal_grid_data$animal_richness > 0))
  message("  Output: ", dirs$animal_model)

  # Run the downstream analysis (skip the GBIF import step)
  source("scripts/utils/ModelAnimal.R")
  modelAnimal(animal_grid_data,
              redoeverything = redo,
              fill_animals = .amani_fillAnimals,
              alphaRichnessMethod = .amani_alphaRichnessMethod,
              fit_params = c("amplitude"))

  message("\nPlant verification animal model complete.")

  if (run_extinction) {
    run_extinction_model(organism)
  }

  invisible(dirs)
}
