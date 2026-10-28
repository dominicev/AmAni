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
  
  data_folder <- "./5-SensitivityTesting/results/"
  rawResultsFile <- "./5-SensitivityTesting/results/p_comparison_results.csv"
  output_folder <- "./5-SensitivityTesting/results/"
  source("./2-ModellingAnimals/scripts/utils/computeCellStats.R")
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


# =============================================================================
# STEP 1: DATA PREP
# =============================================================================

df <- fread(rawResultsFile)
names(df)
# --- Deduplicate pseudoreplicates ---
# Same input_file + parameter combo run multiple times → keep most complete row
# completeness_floor distinguishes the site-completeness sweep runs; without it
# the four floor variants collapse into one row (the 0.3 floor wins the sort).
dedup_keys <- c("input_file", "fit_params", "fill_plants", "fill_animals", "filter_plant_completeness", "completeness_floor")
n_before <- nrow(df)
df[, has_result := !is.na(best_estimate)]
df[, ci_width := ifelse(has_result, upper_bound - lower_bound, Inf)]
setorderv(df, c(dedup_keys, "has_result", "ci_width"), order = c(rep(1, length(dedup_keys)), -1, 1))
df <- df[, .SD[1], by = dedup_keys]
df[, c("has_result", "ci_width") := NULL]
message("Deduplication: ", n_before, " -> ", nrow(df), " rows (", n_before - nrow(df), " duplicates removed)")

unique(df$test)
# --- Create DVs ---
# Signed bias (for correction): positive = overestimate, negative = underestimate
df$bias_magn <- df$best_estimate - df$trueSpeciesRichness
# Ratio bias (for proportional correction): >1 = overestimate, <1 = underestimate
df$bias_prop <- df$best_estimate / df$trueSpeciesRichness
# Unsigned accuracy (for ranking "how accurate"): closer to 0 = better
df$abs_error <- abs(df$bias_magn)
df$abs_prop_error <- abs(df$bias_prop - 1)

# Direction flag for point styling in plots
df$estimate_direction <- ifelse(df$bias_magn > 0, "overestimate", "underestimate")

# Fix string factors to numeric
# sampling_evenness values represent the lower bound on proportion of species sampled per site
df$sampling_evenness <- c(1, 0.5, 0.1)[match(df$sampling_evenness, c("even", "low", "veryLow"))]
df$cellsPerSpp <- df$n_grid_cells / df$n_species
df$SppPerStartCell <- df$n_species/ df$n_grid_cells

#Log transform some variables
df$n_grid_cells <- log(df$n_grid_cells)
df$sample_numberStartSites <- log(df$sample_numberStartSites)
df$trueSpeciesRichness <- log(df$trueSpeciesRichness)


# Recode filter_plant_completeness: NA means "not applied"
df$filter_plant_completeness <- ifelse(is.na(df$filter_plant_completeness), "none",
                                        as.character(df$filter_plant_completeness))
df$n_cells_filtered<-log(df$n_grid_cells_before-df$n_grid_cells_after+.01)
df$completeness_floor

# Convert logical columns to factors
df$fill_plants <- as.character(df$fill_plants)
df$fill_animals <- as.character(df$fill_animals)


# --- Drop zero/near-zero variance columns ---
# prop_sites has no variance (all 1 or NA)
zero_var_cols <- names(df)[sapply(df, function(col) {
  vals <- na.omit(col)
  length(unique(vals)) < 2
})]
message("Dropping zero-variance columns: ", paste(zero_var_cols, collapse = ", "))

# --- Flag near-collinear IVs (r > 0.95) ---
# sample_sd_range_size is nearly perfectly correlated with sample_mean_range_size
# Keep sample_mean_range_size (more interpretable)
# Also going to keep  n_grid_cells because it is critical to the filtering steps
collinear_drops <- c( "sample_sd_range_size",
                     "sample_median_range_size",
                     "mean_rangeSize", "median_range_size", "sd_range_size")
message("Dropping near-collinear columns: ", paste(collinear_drops, collapse = ", "))

names(df)

# Also drop metadata columns that are not IVs
metadata_drops <- c("scenario", "test", "param_names_fit", 
                     "delta_meanlog_status",  # status flag, not a parameter
                     "sampling_fraction",     # derived, potentially collinear
                     "min_completeness",      # derived
                     "lower_bound", "upper_bound",
                     "input_file")            # identifier, not a treatment

all_drops <- unique(c(zero_var_cols, collinear_drops, metadata_drops))

# --- Subset to analysis-ready data ---
# Keep only rows where the model produced a result
df_analysis <- df[!is.na(df$bias_magn), ]
message("Rows with results: ", nrow(df_analysis), " / ", nrow(df), 
        " (", nrow(df) - nrow(df_analysis), " model failures dropped)")

# --- Define DVs and IVs ---
dv_names <- c("bias_magn", "bias_prop", "abs_error", "abs_prop_error")

# IVs = all treatment columns not in drops or DVs
my_treatments <- c("n_species", "n_grid_cells","sample_mean_range_size", "total_observations",
                   "prop_species",
                   "clustering_label", "sample_numberStartSites", "sample_diffusion",
                   "sampling_evenness", "trueSpeciesRichness",
                   "delta_meanlog", "ratio_sdlog", 
                   "sample_endemism_min", "sample_endemism_max",
                   "cellsPerSpp", "SppPerStartCell",
                   "fill_plants", "fill_animals", "filter_plant_completeness", "fit_params",
                   "n_grid_cells_before", "n_grid_cells_after", "n_grid_cells_before", "completeness_floor", "n_cells_filtered"
                   )

# Final check: which treatments actually exist and have variance
my_treatments <- intersect(my_treatments, names(df_analysis))
my_treatments <- setdiff(my_treatments, all_drops)
message("Final IV set (", length(my_treatments), "): ", paste(my_treatments, collapse = ", "))

# =============================================================================
# STEP 1b: BALANCED SUBSETTING FOR NEW-PARAMETER COMPARISONS
# =============================================================================
# The experiment uses a one-factor-at-a-time (OFAT) design around a single
# baseline cell. Each of the four "new parameter" axes below was varied one at
# a time while the other three were held at their baseline values. Consequently
# the baseline level of each axis (e.g. fit_params = "amplitude, meanlog, sdlog")
# accumulates rows from every OTHER axis's variation and is massively
# overrepresented relative to the non-baseline levels.
#
# To compare levels of a focal axis fairly, we hold the OTHER new-param axes at
# baseline so every focal level is evaluated on the same matched set of
# replicates. Simulation IVs (n_species, sample_mean_range_size, ...) are left
# pooled across all new-param combinations.
# Baseline (default) value of each one-at-a-time new-param axis.
new_param_baseline <- list(
  fill_plants               = "FALSE",
  fill_animals              = "FALSE",
  filter_plant_completeness = "none",
  fit_params                = "amplitude, meanlog, sdlog",
  completeness_floor = NA
)

#' Restrict data to a balanced subset for comparing a focal axis
#'
#' For a focal variable that is one of the four new-param axes (or the derived
#' "fill_interaction", which varies fill_plants and fill_animals jointly), returns
#' the rows holding every OTHER new-param axis at its baseline value. For any
#' other (simulation) IV the data are returned unchanged.
#'
#' @param data data.frame or data.table.
#' @param focal Character name of the focal comparison variable.
#' @param baseline Named list of baseline values (default new_param_baseline).
#' @param verbose Logical. Report the number of rows retained.
#' @return Subset of data (same class as input).
balance_focal <- function(data, focal, baseline = new_param_baseline,
                          verbose = TRUE) {

  # fill_interaction is a joint recoding of fill_plants and fill_animals
  focal_axes <- if (identical(focal, "fill_interaction")) {
    c("fill_plants", "fill_animals")
  } else {
    focal
  }

  # Axes to pin at baseline = new-param axes that are NOT the focal axis
  hold_axes <- setdiff(names(baseline), focal_axes)

  # If focal is not a new-param axis, hold_axes == all axes -> leave pooled
  if (length(hold_axes) == length(baseline)) {
    return(data)
  }

  keep <- rep(TRUE, nrow(data))
  for (ax in hold_axes) {
    keep <- keep & data[[ax]] %in% baseline[[ax]]
  }

  if (verbose) {
    message("  [balance] ", focal, ": ", sum(keep), "/", length(keep),
            " rows (holding ", paste(hold_axes, collapse = ", "),
            " at baseline)")
  }

  data[keep, ]
}

# =============================================================================
# STEP 2a: IMPORT SAMPLE DATASETS AND COMPUTE SITE COMPLETENESS
# =============================================================================
# The `input_file` column of df_analysis names the .RData sample dataset that
# produced each row. These files live in the sampleDatasets output folder.
# Each .RData contains a single data.table named `df` with columns including
# `grid_id` (site) and `species` (taxon).
#
# Because these are simulated data, the TRUE community at every site is known
# (identical across replicates). We therefore compute *actual* sampling
# completeness per site = observed richness / true richness, rather than
# estimating it with Chao1.

# --- True (simulated) community: one row per species x grid_id presence ---
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

#' Compute actual per-site sampling completeness for one sample dataset
#'
#' @param occ_data data.frame of sampled occurrence records.
#' @param true_richness data.frame with columns `grid_id` and `n_true`
#'   (true species richness per site).
#' @param site_col Character. Site / grid-cell ID column. Default "grid_id".
#' @param species_col Character. Species column. Default "species".
#' @return data.frame: one row per sampled site with `grid_id`, `n_obs`,
#'   `n_true`, and `completeness` (= n_obs / n_true).
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

# Folder holding the .RData sample datasets referenced by input_file
sample_dataset_dir <- file.path(
  "C:/Users/dominice/OneDrive - University of Illinois - Urbana",
  "Projects/Peruvian Blattodea/Cockroaches in a coal mine",
  "AmAni/5-SensitivityTesting/output/sampleDatasets"
)

#' Load the occurrence data.table from one sample-dataset .RData file
#'
#' Loads into a scratch environment so the global `df` is never clobbered.
#'
#' @param input_file Character. Basename of the .RData file (as stored in
#'   df_analysis$input_file).
#' @param dataset_dir Character. Directory containing the .RData files.
#' @param verbose Logical. Print status messages. Default TRUE.
#' @return The occurrence data.table stored inside the .RData file.
load_sample_dataset <- function(input_file,
                                dataset_dir = sample_dataset_dir,
                                verbose     = TRUE) {

  input_path <- file.path(dataset_dir, input_file)
  if (!file.exists(input_path)) {
    stop("Sample dataset not found: ", input_path)
  }

  if (verbose) message("[load_sample_dataset] Loading: ", input_file)

  # Load into an isolated environment to protect the global `df`
  load_env <- new.env()
  loaded_names <- load(input_path, envir = load_env)

  # The file stores a single object (a data.table named `df`). Coerce to a plain
  # data.frame: compute_cell_stats() relies on data.frame `[, cols]` semantics,
  # which data.table overrides.
  occ_data <- as.data.frame(get(loaded_names[1], envir = load_env))

  if (verbose) {
    message("[load_sample_dataset] Loaded object '", loaded_names[1],
            "' with ", nrow(occ_data), " rows.")
  }

  return(occ_data)
}

#' Load compute the desired stats from the cell-by-cell completeness values
#'
#' Loads into a scratch environment so the global `df` is never clobbered.
#'
#' @param completeness_computed Output of compute_site_completeness.
#' @return A list of three statistics
compute_completeness_stats <- function(completeness_computed) {
  list(
    comp_mean   = mean(completeness_computed$completeness),
    comp_sd     = sd(completeness_computed$completeness),
    comp_spread = max(completeness_computed$completeness) - min(completeness_computed$completeness)
  )
}
# --- Import one roachLike dataset and compute site completeness ---
roach_file <- "blatt_blatt_roachLike.sample_20260626_105416.RData"
roach_occ  <- load_sample_dataset(roach_file)
roach_completeness <- compute_site_completeness(roach_occ, true_richness)
compute_completeness_stats(roach_completeness)

names(df_analysis)
min(df_analysis$total_observations/exp(df_analysis$n_grid_cells))
max(df_analysis$total_observations/exp(df_analysis$n_grid_cells))
mean(df_analysis$total_observations/exp(df_analysis$n_grid_cells))


# --- Loop over all unique input_files and compute completeness stats ---
# Build a lookup table (one row per unique input_file) then join to df_analysis.
unique_files <- unique(df_analysis$input_file)
message("Computing completeness for ", length(unique_files), " unique input files...")

completeness_lookup <- do.call(rbind, lapply(unique_files, function(fname) {
  occ_data      <- load_sample_dataset(fname, verbose = FALSE)
  site_comp     <- compute_site_completeness(occ_data, true_richness)
  comp_stats    <- compute_completeness_stats(site_comp)
  data.frame(
    input_file  = fname,
    comp_mean   = comp_stats$comp_mean,
    comp_sd     = comp_stats$comp_sd,
    comp_spread = comp_stats$comp_spread
  )
}))

message("Completeness lookup complete. Files processed: ", nrow(completeness_lookup))

# Add columns to df_analysis by merging on input_file.
# Drop any pre-existing completeness columns to avoid duplicates on re-runs.
df_analysis <- df_analysis[, setdiff(names(df_analysis),
                                     c("comp_mean", "comp_sd", "comp_spread")),
                           with = FALSE]
df_analysis <- merge(df_analysis, completeness_lookup,
                     by = "input_file", all.x = TRUE)

message("Columns added. NAs in comp_mean: ",
        sum(is.na(df_analysis$comp_mean)))

# Register completeness stats as IVs for the sensitivity analysis. Done here
# (rather than in the STEP 1 my_treatments definition) because these columns
# only exist after the merge above; the earlier intersect() would drop them.
my_treatments <- unique(c(my_treatments, "comp_mean", "comp_sd", "comp_spread"))
message("Final IV set (", length(my_treatments), "): ",
        paste(my_treatments, collapse = ", "))

# =============================================================================
# STEP 2b: ACCURACY AND BIAS SUMMARY BY LEVEL
# =============================================================================

#' Summarize accuracy and bias per level of each IV
#'
#' For continuous IVs, bins into quantile-based groups so levels are balanced.
#' Reports n, mean bias (= correction factor), mean absolute error (= accuracy),
#' SD, and 95% CI per level. Identifies the "most accurate" level.
#'
#' @param data data.frame of analysis-ready data
#' @param iv_names character vector of IV column names
#' @param dv_bias character name of the signed bias DV (e.g., "bias_magn")
#' @param dv_accuracy character name of the unsigned accuracy DV (e.g., "abs_error")
#' @param n_bins integer number of quantile bins for continuous IVs (default 4)
#' @param verbose logical print status messages
#' @return list with two elements: $by_level (full breakdown) and $best_levels (most accurate level per IV)
summarize_accuracy_by_level <- function(data, iv_names, dv_bias, dv_accuracy,
                                        n_bins = 4, verbose = TRUE) {
  
  results_list <- list()
  
  for (iv in iv_names) {
    
    # Skip if column doesn't exist or has no variance
    if (!iv %in% names(data)) next
    # Balance OFAT comparisons of the four new-param axes: hold the other
    # new-param axes at baseline so each focal level uses matched replicates.
    # Simulation IVs are returned unchanged (pooled).
    iv_data <- balance_focal(data, iv, verbose = FALSE)
    vals <- iv_data[[iv]]
    if (all(is.na(vals)) || length(unique(na.omit(vals))) < 2) next
    
    # Determine grouping variable
    if (is.numeric(vals) && length(unique(na.omit(vals))) > n_bins * 2) {
      # Continuous IV: bin into quantiles
      breaks <- quantile(vals, probs = seq(0, 1, length.out = n_bins + 1), na.rm = TRUE)
      breaks <- unique(breaks) # collapse ties
      if (length(breaks) < 3) next # not enough unique quantiles
      group_var <- cut(vals, breaks = breaks, include.lowest = TRUE, dig.lab = 5)
      iv_type <- "continuous (binned)"
    } else {
      # Categorical or discrete IV: use as-is
      group_var <- as.factor(vals)
      iv_type <- "categorical"
    }
    
    # Compute per-level summaries (from the balanced subset for this IV)
    bias_vals <- iv_data[[dv_bias]]
    acc_vals  <- iv_data[[dv_accuracy]]
    
    levels_in_data <- levels(group_var)
    for (lev in levels_in_data) {
      
      mask <- !is.na(group_var) & group_var == lev & !is.na(bias_vals)
      n_lev <- sum(mask)
      if (n_lev < 2) next
      
      bias_lev <- bias_vals[mask]
      acc_lev  <- acc_vals[mask]
      
      # 95% CI for mean bias (the correction factor)
      se_bias <- sd(bias_lev, na.rm = TRUE) / sqrt(n_lev)
      ci_lower <- mean(bias_lev, na.rm = TRUE) - 1.96 * se_bias
      ci_upper <- mean(bias_lev, na.rm = TRUE) + 1.96 * se_bias
      
      results_list[[length(results_list) + 1]] <- data.frame(
        IV          = iv,
        IV_Type     = iv_type,
        Level       = as.character(lev),
        N           = n_lev,
        Mean_Bias   = mean(bias_lev, na.rm = TRUE),   # correction factor
        SD_Bias     = sd(bias_lev, na.rm = TRUE),
        CI_Lower    = ci_lower,
        CI_Upper    = ci_upper,
        Mean_AbsErr = mean(acc_lev, na.rm = TRUE),     # accuracy (lower = better)
        Min_AbsErr  = min(acc_lev, na.rm = TRUE),
        Max_AbsErr  = max(acc_lev, na.rm = TRUE),
        stringsAsFactors = FALSE
      )
    }
    
    if (verbose) message("  Summarized: ", iv, " (", iv_type, ", ", 
                          length(levels_in_data), " levels)")
  }
  
  if (length(results_list) == 0) {
    warning("No valid IV levels found.")
    return(NULL)
  }
  
  by_level <- do.call(rbind, results_list)
  rownames(by_level) <- NULL
  
  # Identify best level per IV (lowest mean absolute error)
  best_levels <- do.call(rbind, lapply(split(by_level, by_level$IV), function(sub) {
    sub[which.min(sub$Mean_AbsErr), ]
  }))
  rownames(best_levels) <- NULL
  
  return(list(by_level = by_level, best_levels = best_levels))
}


# =============================================================================
# STEP 2a: VARIABLE IMPORTANCE (multivariate)
# =============================================================================

#' Rank variable importance using Type III partial eta-squared
#'
#' Fits a single linear model with all IVs to control for confounding.
#' Uses car::Anova Type III SS so results are valid with unbalanced data.
#' Categorical IVs with too few levels or too few observations are skipped.
#'
#' @param data data.frame (or data.table)
#' @param dv character name of the dependent variable
#' @param iv_names character vector of IV names
#' @param min_obs_per_level integer minimum observations per factor level (default 3)
#' @param max_na_frac numeric max fraction of NAs allowed per IV (default 0.5);
#'   IVs exceeding this are skipped to preserve complete-case sample size
#' @param verbose logical
#' @return data.frame with IV, Partial_Eta_Sq, F_Value, P_Value, DF, N_Complete
rank_variable_importance <- function(data, dv, iv_names, min_obs_per_level = 3,
                                     max_na_frac = 0.5, verbose = TRUE) {
  
  library(car)
  
  # Filter IVs: must exist, have variance, not too many NAs,
  # and categorical IVs need sufficient obs per level
  valid_ivs <- character(0)
  for (iv in iv_names) {
    if (!iv %in% names(data)) next
    vals <- data[[iv]]
    non_na <- na.omit(vals)
    if (length(unique(non_na)) < 2) next
    
    # Skip IVs with too many NAs (they destroy complete cases)
    na_frac <- sum(is.na(vals)) / length(vals)
    if (na_frac > max_na_frac) {
      if (verbose) message("  Skipping ", iv, " (", round(na_frac * 100, 1), "% NA)")
      next
    }
    
    if (!is.numeric(vals)) {
      tab <- table(vals)
      if (any(tab < min_obs_per_level)) {
        if (verbose) message("  Skipping ", iv, " (factor level with <", 
                              min_obs_per_level, " obs)")
        next
      }
    }
    valid_ivs <- c(valid_ivs, iv)
  }
  
  if (length(valid_ivs) == 0) {
    warning("No valid IVs remaining.")
    return(NULL)
  }
  
  if (verbose) message("Fitting model with ", length(valid_ivs), " IVs for DV: ", dv)
  
  # Build formula
  formula_str <- paste0("`", dv, "` ~ ", paste0("`", valid_ivs, "`", collapse = " + "))
  
  # Subset to complete cases for these columns
  model_data <- data[, c(dv, valid_ivs), with = FALSE]
  model_data <- na.omit(model_data)
  
  if (verbose) message("  Complete cases: ", nrow(model_data))
  
  if (nrow(model_data) < length(valid_ivs) + 5) {
    warning("Too few complete cases (", nrow(model_data), 
            ") relative to IVs (", length(valid_ivs), "). Reduce IV set.")
    return(NULL)
  }
  
  # Convert character columns to factors with sum contrasts for Type III SS
  char_cols <- names(model_data)[sapply(model_data, is.character)]
  for (cc in char_cols) {
    model_data[[cc]] <- as.factor(model_data[[cc]])
    contrasts(model_data[[cc]]) <- contr.sum(nlevels(model_data[[cc]]))
  }
  
  model <- lm(as.formula(formula_str), data = model_data)
  
  anova_result <- tryCatch(
    car::Anova(model, type = "III"),
    error = function(e) {
      warning("Type III ANOVA failed: ", e$message, ". Falling back to Type II.")
      car::Anova(model, type = "II")
    }
  )
  
  # Calculate partial eta-squared: SS_effect / (SS_effect + SS_residual)
  ss <- anova_result[["Sum Sq"]]
  f_vals <- anova_result[["F value"]]
  p_vals <- anova_result[["Pr(>F)"]]
  df_vals <- anova_result[["Df"]]
  term_names <- rownames(anova_result)
  
  # Find residual SS
  resid_idx <- which(term_names == "Residuals")
  ss_resid <- ss[resid_idx]
  
  # Remove Intercept and Residuals rows
  keep <- !term_names %in% c("(Intercept)", "Residuals")
  
  result <- data.frame(
    DV              = dv,
    IV              = gsub("`", "", term_names[keep]),
    Partial_Eta_Sq  = ss[keep] / (ss[keep] + ss_resid),
    F_Value         = f_vals[keep],
    P_Value         = p_vals[keep],
    DF              = df_vals[keep],
    N_Complete      = nrow(model_data),
    stringsAsFactors = FALSE
  )
  
  # Sort by partial eta-squared descending
  result <- result[order(-result$Partial_Eta_Sq), ]
  rownames(result) <- NULL
  
  if (verbose) {
    message("  Top 10 IVs by partial eta-squared:")
    top5 <- head(result, 10)
    for (i in seq_len(nrow(top5))) {
      message("    ", top5$IV[i], ": eta2p = ", round(top5$Partial_Eta_Sq[i], 3),
              ", p = ", signif(top5$P_Value[i], 3))
    }
  }
  
  return(result)
}


# =============================================================================
# RUN ANALYSES
# =============================================================================

df_analysis

# --- Step 2b: Per-level accuracy summaries ---
message("\n===== MAGNITUDE BIAS SUMMARY =====")
magn_summary <- summarize_accuracy_by_level(
  data = df_analysis, iv_names = my_treatments,
  dv_bias = "bias_magn", dv_accuracy = "abs_error",
  n_bins = 4, verbose = TRUE
)

message("\n===== PROPORTIONAL BIAS SUMMARY =====")
prop_summary <- summarize_accuracy_by_level(
  data = df_analysis, iv_names = my_treatments,
  dv_bias = "bias_prop", dv_accuracy = "abs_prop_error",
  n_bins = 4, verbose = TRUE
)
prop_summary$by_level

# Most accurate parameter levels
message("\n--- Best levels (magnitude) ---")
magn_summary$best_levels

message("\n--- Best levels (proportional) ---")
prop_summary$best_levels

varianceinBias<-lapply(
  quantile( df_analysis$bias_prop, probs = c(0.1, 0.9), na.rm = TRUE), function(x)round(x, digits = 3))

# Overall mean bias (the "global correction factor")
message("\n--- Overall correction factors ---")
message("  Magnitude bias:     mean = ", round(mean(df_analysis$bias_magn, na.rm = TRUE), 1),
        ", median = ", round(median(df_analysis$bias_magn, na.rm = TRUE), 1))
message("  Proportional bias:  mean = ", round(mean(df_analysis$bias_prop, na.rm = TRUE), 3),
        ", median = ", round(median(df_analysis$bias_prop, na.rm = TRUE), 3))
message("  Variance in bias:  95% range = ", varianceinBias[[1]], " to ", varianceinBias[[2]])


# --- Step 2a: Variable importance ---
message("\n===== VARIABLE IMPORTANCE (Magnitude) =====")
vi_magn <- rank_variable_importance(
  data = df_analysis, dv = "bias_magn",
  iv_names = my_treatments, verbose = TRUE
)

message("\n===== VARIABLE IMPORTANCE (Proportional) =====")
vi_prop <- rank_variable_importance(
  data = df_analysis, dv = "bias_prop",
  iv_names = my_treatments, verbose = TRUE
)


# =============================================================================
# STEP 3: VISUALIZATIONS
# =============================================================================

# --- 3a. Variable importance comparison plot ---
if (!is.null(vi_magn) && !is.null(vi_prop)) {
  vi_combined <- rbind(vi_magn, vi_prop)
  vi_combined$Significant <- ifelse(vi_combined$P_Value < 0.05, "p < 0.05", "n.s.")
  
  ggplot(vi_combined, aes(x = Partial_Eta_Sq, 
                           y = reorder(IV, Partial_Eta_Sq),
                           fill = Significant)) +
    geom_col() +
    facet_wrap(~ DV, scales = "free_x") +
    scale_fill_manual(values = c("p < 0.05" = "steelblue", "n.s." = "grey70")) +
    labs(x = "Partial Eta-Squared (Type III)", y = NULL,
         title = "Variable Importance: Controlling for Confounding") +
    theme_minimal(base_size = 12) +
    theme(legend.position = "bottom")
}


# --- 3b. Bias correction plot (top IVs only) ---
# Show mean bias ± CI for the top 6 most important IVs (magnitude)
if (!is.null(vi_magn) && !is.null(magn_summary)) {
  top_ivs <- head(vi_magn$IV, 6)
  bias_top <- magn_summary$by_level[magn_summary$by_level$IV %in% top_ivs, ]
  
  if (nrow(bias_top) > 0) {
    # Order IV factor by importance
    bias_top$IV <- factor(bias_top$IV, levels = rev(top_ivs))
    
    # Determine direction for each level based on sign of Mean_Bias
    bias_top$estimate_direction <- ifelse(bias_top$Mean_Bias > 0, "overestimate", "underestimate")
    
    ggplot(bias_top, aes(x = Mean_Bias, y = Level)) +
      geom_point(aes(size = N, fill = estimate_direction), shape = 21, color = "black") +
      scale_fill_manual(values = c("overestimate" = "white", "underestimate" = "grey60"),
                        name = "Direction") +
      geom_errorbarh(aes(xmin = CI_Lower, xmax = CI_Upper), height = 0.2) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
      facet_wrap(~ IV, scales = "free_y", ncol = 2) +
      labs(x = "Mean Bias (estimate - truth)", y = NULL, size = "n",
           title = "Bias Correction by Parameter Level",
           subtitle = "Top 6 IVs by importance. Red line = no bias.") +
      theme_minimal(base_size = 11) +
      theme(strip.text = element_text(face = "bold"))
  }
}


# --- 3c. Calibration plot ---
ggplot(df_analysis, aes(x = trueSpeciesRichness, y = best_estimate)) +
  geom_point(aes(fill = estimate_direction), shape = 21, alpha = 0.5,
             color = "black", size = 2) +
  scale_fill_manual(values = c("overestimate" = "white", "underestimate" = "grey60"),
                    name = "Direction") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(x = "True Species Richness", y = "Best Estimate",
       title = "Calibration: Predicted vs. True Species Richness",
       subtitle = "Dashed red = perfect prediction") +
  theme_minimal(base_size = 12)


# =============================================================================
# STEP 4: INTERACTION TESTING (scaffold)
# =============================================================================

#' Test pairwise interactions between top IVs
#'
#' Two-stage approach to avoid aliasing in interaction models:
#'   1. Pre-screen: drop interaction pairs where both IVs are numeric and |r| > cor_threshold.
#'      This catches the most common multicollinearity case (e.g., n_species * total_obs).
#'   2. Post-fit: use alias() to detect any remaining aliased terms (e.g., from categorical
#'      IVs with sparse level combinations), drop them, and refit (up to 3 passes).
#'
#' Continuous IVs are centered and scaled before fitting so interaction terms are
#' not aliased with main effects.
#'
#' @param data data.frame or data.table
#' @param dv character DV name
#' @param top_ivs character vector of IV names to test interactions for
#' @param cor_threshold numeric drop numeric-numeric pairs with |r| above this (default 0.7)
#' @param verbose logical
#' @return data.frame of interaction terms with Partial_Eta_Sq, F_Value, P_Value
test_interactions <- function(data, dv, top_ivs, cor_threshold = 0.7, verbose = TRUE) {
  
  library(car)
  
  if (length(top_ivs) < 2) {
    message("Need at least 2 IVs to test interactions.")
    return(NULL)
  }
  
  # Subset to complete cases
  model_data <- data[, c(dv, top_ivs), with = FALSE]
  model_data <- as.data.frame(na.omit(model_data))
  
  # --- Center and scale continuous IVs ---
  # Required so interaction terms are not aliased with main effects
  for (iv in top_ivs) {
    if (is.numeric(model_data[[iv]])) {
      model_data[[iv]] <- scale(model_data[[iv]], center = TRUE, scale = TRUE)[, 1]
    } else {
      model_data[[iv]] <- as.factor(model_data[[iv]])
      contrasts(model_data[[iv]]) <- contr.sum(nlevels(model_data[[iv]]))
    }
  }
  
  # --- Stage 1: pre-screen numeric-numeric pairs by correlation ---
  pairs_all  <- combn(top_ivs, 2, simplify = FALSE)
  pairs_keep <- list()
  for (pair in pairs_all) {
    iv1 <- model_data[[pair[1]]]; iv2 <- model_data[[pair[2]]]
    if (is.numeric(iv1) && is.numeric(iv2)) {
      r <- abs(cor(iv1, iv2, use = "complete.obs"))
      if (r > cor_threshold) {
        if (verbose) message("  Pre-screen dropping ", pair[1], ":", pair[2],
                              " (|r|=", round(r, 2), ")")
        next
      }
    }
    pairs_keep <- c(pairs_keep, list(pair))
  }
  
  if (length(pairs_keep) == 0) {
    warning("All pairs dropped due to collinearity. Lower cor_threshold.")
    return(NULL)
  }
  
  n_params <- length(top_ivs) + length(pairs_keep)
  if (nrow(model_data) < n_params + 10) {
    warning("Too few cases (", nrow(model_data), ") for ", n_params, " terms.")
    return(NULL)
  }
  
  if (verbose) message("  Complete cases: ", nrow(model_data),
                        ", candidate interaction terms: ", length(pairs_keep))
  
  # --- Stage 2: fit, detect aliased terms, drop and refit (up to 3 passes) ---
  build_and_clean_model <- function(pairs_candidate) {
    main_effects <- paste0(top_ivs, collapse = " + ")
    int_terms    <- sapply(pairs_candidate, function(p) paste0(p[1], ":", p[2]))
    fml <- paste0(dv, " ~ ", main_effects, " + ", paste(int_terms, collapse = " + "))
    m   <- lm(as.formula(fml), data = model_data)
    
    for (pass in 1:3) {
      aliased <- alias(m)$Complete
      if (is.null(aliased) || nrow(aliased) == 0) break
      aliased_terms <- rownames(aliased)
      if (verbose) message("  Pass ", pass, ": removing aliased terms: ",
                            paste(aliased_terms, collapse = ", "))
      pairs_candidate <- Filter(function(p) {
        !paste0(p[1], ":", p[2]) %in% aliased_terms &&
        !paste0(p[2], ":", p[1]) %in% aliased_terms
      }, pairs_candidate)
      if (length(pairs_candidate) == 0) return(NULL)
      int_terms <- sapply(pairs_candidate, function(p) paste0(p[1], ":", p[2]))
      fml <- paste0(dv, " ~ ", main_effects, " + ", paste(int_terms, collapse = " + "))
      m   <- lm(as.formula(fml), data = model_data)
    }
    list(model = m, pairs = pairs_candidate)
  }
  
  result_obj <- build_and_clean_model(pairs_keep)
  if (is.null(result_obj)) {
    warning("No non-aliased interaction terms remain.")
    return(NULL)
  }
  model <- result_obj$model
  
  anova_result <- tryCatch(
    car::Anova(model, type = "III"),
    error = function(e) {
      message("  Type III failed: ", e$message, ". Trying Type II.")
      tryCatch(car::Anova(model, type = "II"),
               error = function(e2) { warning("Type II also failed: ", e2$message); NULL })
    }
  )
  if (is.null(anova_result)) return(NULL)
  
  ss <- anova_result[["Sum Sq"]]
  term_names <- rownames(anova_result)
  ss_resid <- ss[term_names == "Residuals"]
  interaction_mask <- grepl(":", term_names)
  if (!any(interaction_mask)) { message("No interaction terms in output."); return(NULL) }
  
  result <- data.frame(
    DV             = dv,
    Interaction    = term_names[interaction_mask],
    Partial_Eta_Sq = ss[interaction_mask] / (ss[interaction_mask] + ss_resid),
    F_Value        = anova_result[["F value"]][interaction_mask],
    P_Value        = anova_result[["Pr(>F)"]][interaction_mask],
    stringsAsFactors = FALSE
  )
  result <- result[order(-result$Partial_Eta_Sq), ]
  rownames(result) <- NULL
  
  if (verbose) {
    sig <- result[result$P_Value < 0.05, ]
    if (nrow(sig) > 0) {
      message("  Significant interactions (p < 0.05):")
      for (i in seq_len(nrow(sig)))
        message("    ", sig$Interaction[i], ": eta2p=", round(sig$Partial_Eta_Sq[i], 3),
                ", p=", signif(sig$P_Value[i], 3))
    } else {
      message("  No significant interactions at p < 0.05")
    }
  }
  return(result)
}

# Example usage (uncomment when ready):
 top_6 <- head(vi_magn$IV, 6)
 interaction_results <- test_interactions(df_analysis, "bias_magn", top_6)

 interaction_results
# =============================================================================
# STEP 5: LINEAR MODELS OF ACCURACY BY EXPERIMENT
# =============================================================================
# Fit simple linear models within each experiment (test), using Type III SS
# to handle unequal sample sizes across factor levels.

library(car)
library(broom)
 
 df
 # ===== FOR DEBUGGING: Initialize with variedSiteEvenness          =====
 data <- df_analysis
 dv <- "abs_error"
 iv_names <- my_treatments
 verbose <- TRUE
 exp <- "variedSiteEvenness"
 # ================================================================
 # 
lm_by_experiment <- function(data, dv, iv_names, focal_iv_map = NULL,
                              plot_dir = NULL, verbose = TRUE) {
  # focal_iv_map: optional named character vector mapping experiment name -> focal IV column.
  # If provided, only that IV is used for the given experiment instead of all iv_names.
  
  
  
  experiments <- unique(data$test)
  results_list <- list()
  
  for (exp in experiments) {
    exp_data <- data[data$test == exp, ]
    
    # Use focal IV only if specified for this experiment, otherwise use all iv_names
    candidate_ivs <- if (!is.null(focal_iv_map) && exp %in% names(focal_iv_map)) {
      focal_iv_map[[exp]]
    } else {
      iv_names
    }
    
    # Identify usable IVs for this experiment
    valid_ivs <- character(0)
    for (iv in candidate_ivs) {
      if (!iv %in% names(exp_data)) next
      vals <- na.omit(exp_data[[iv]])
      if (length(unique(vals)) < 2) next
      valid_ivs <- c(valid_ivs, iv)
    }
    if (length(valid_ivs) == 0) {
      if (verbose) message("  ", exp, ": no valid IVs, skipping")
      next
    }
    
    # # Drop IVs that have NAs within this experiment so na.omit doesn't
    # # unnecessarily discard rows due to columns irrelevant to that experiment
    # valid_ivs <- valid_ivs[sapply(valid_ivs, function(iv) !any(is.na(exp_data[[iv]])))]
    # if (length(valid_ivs) == 0) {
    #   if (verbose) message("  ", exp, ": all IVs have NAs in this experiment, skipping")
    #   next
    # }
    
    # Include estimate_direction for plot styling (not used in model)
    model_cols <- c(dv, valid_ivs)
    if ("estimate_direction" %in% names(exp_data)) model_cols <- c(model_cols, "estimate_direction")
    model_data <- as.data.frame(na.omit(exp_data[, model_cols, with = FALSE]))
    # Need strictly more observations than parameters to fit
    if (nrow(model_data) <= length(valid_ivs)) {
      if (verbose) message("  ", exp, ": n=", nrow(model_data),
                            " <= ", length(valid_ivs), " IVs, skipping")
      next
    }
    
    # Set sum contrasts for categorical IVs (required for Type III SS)
    char_cols <- names(model_data)[sapply(model_data, function(x) is.character(x) || is.factor(x))]
    for (cc in char_cols) {
      model_data[[cc]] <- as.factor(model_data[[cc]])
      if (nlevels(model_data[[cc]]) >= 2)
        contrasts(model_data[[cc]]) <- contr.sum(nlevels(model_data[[cc]]))
    }
    
    fml <- paste0("`", dv, "` ~ ", paste0("`", valid_ivs, "`", collapse = " + "))
    fit <- tryCatch(lm(as.formula(fml), data = model_data), error = function(e) NULL)
    if (is.null(fit)) next
    
    # Type III ANOVA for unequal sample sizes
    aov_res <- tryCatch(car::Anova(fit, type = "III"), error = function(e) NULL)
    if (is.null(aov_res)) next
    
    ss <- aov_res[["Sum Sq"]]
    term_names <- rownames(aov_res)
    ss_resid <- ss[term_names == "Residuals"]
    keep <- !term_names %in% c("(Intercept)", "Residuals")
    
    exp_result <- data.frame(
      Experiment     = exp,
      DV             = dv,
      IV             = gsub("`", "", term_names[keep]),
      Partial_Eta_Sq = ss[keep] / (ss[keep] + ss_resid),
      F_Value        = aov_res[["F value"]][keep],
      P_Value        = aov_res[["Pr(>F)"]][keep],
      N              = nrow(model_data),
      R2_adj         = summary(fit)$adj.r.squared,
      stringsAsFactors = FALSE
    )
    results_list[[length(results_list) + 1]] <- exp_result
    
    # --- Regression plots ---
    if (!is.null(plot_dir)) {
      library(patchwork)
      dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)
      
      iv_plots <- lapply(valid_ivs, function(iv) {
        is_cat <- !is.numeric(model_data[[iv]])
        p <- ggplot(model_data, aes(x = .data[[iv]], y = .data[[dv]]))
        if (is_cat) {
          p <- p +
            geom_jitter(aes(fill = estimate_direction), shape = 21, width = 0.2,
                        alpha = 0.6, color = "black", size = 2) +
            stat_summary(fun = mean, geom = "crossbar", width = 0.4,
                         color = "steelblue", linewidth = 0.8) +
            scale_fill_manual(values = c("overestimate" = "white", "underestimate" = "grey60"),
                              name = "Direction")
        } else {
          p <- p +
            geom_point(aes(fill = estimate_direction), shape = 21,
                       alpha = 0.6, color = "black", size = 2) +
            scale_fill_manual(values = c("overestimate" = "white", "underestimate" = "grey60"),
                              name = "Direction") +
            geom_smooth(method = "lm", se = TRUE, color = "steelblue", formula = y ~ x)
        }
        # Annotate with partial eta-sq and p from exp_result
        iv_row <- exp_result[exp_result$IV == iv, ]
        subtitle <- if (nrow(iv_row) == 1) {
          paste0("eta2p = ", round(iv_row$Partial_Eta_Sq, 3),
                 ",  p = ", signif(iv_row$P_Value, 2))
        } else NULL
        p + labs(x = iv, y = dv, subtitle = subtitle) +
          theme_minimal(base_size = 11)
      })
      
      n_plots <- length(iv_plots)
      ncols   <- min(n_plots, 3)
      nrows   <- ceiling(n_plots / ncols)
      combined <- wrap_plots(iv_plots, ncol = ncols) +
        plot_annotation(title   = paste0(exp, "  |  DV: ", dv),
                        caption = paste0("adj. R² = ", round(unique(exp_result$R2_adj), 3),
                                         "  |  n = ", nrow(model_data)),
                        theme   = theme_minimal(base_size = 13))
      
      fname <- file.path(plot_dir, paste0(exp, "_", dv, ".svg"))#change to PNG?
      ggsave(fname, combined,
             width  = 4 * ncols,
             height = 4 * nrows,
             limitsize = FALSE)
      if (verbose) message("    Plot saved: ", fname)
    }
    
    if (verbose) {
      sig <- sum(exp_result$P_Value < 0.05, na.rm = TRUE)
      message("  ", exp, ": ", length(valid_ivs), " IVs, n=", nrow(model_data),
              ", adj.R2=", round(summary(fit)$adj.r.squared, 3),
              ", ", sig, " significant terms")
    }
  }
  
  if (length(results_list) == 0) return(NULL)
  do.call(rbind, results_list) |> (\(d) { rownames(d) <- NULL; d })()
}
df_analysis[df_analysis$test %in% "variedSiteEvenness"]
names(df_analysis)
# --- Focal IV per experiment ---
# For each experiment, name the single column that was the primary treatment axis.
# Experiments not listed here will fall back to using all of my_treatments.
# Use NA (or omit) for experiments with no single focal IV.
focal_iv_map <- c(
  # Original simulation experiments
  lowtoHighSiteSampling    = "n_grid_cells",
  variEndemSampling        = "sample_mean_range_size",
  variedSiteEvenness       = "sampling_evenness",
  variedSiteClustering     = "sample_numberStartSites",
  lowtoHighSpeciesSampling = "prop_species",
  lowtoHighSiteEnd         = "cellsPerSpp",
  lowtoHighSiteRichness    = "SppPerStartCell",
  variRichSampling         = "trueSpeciesRichness",
  animal_filtering_mag = "n_cells_filtered",
  animal_filtering_compl = "completeness_floor"
  
  
  # Organism-like tests not listed → fall back to all my_treatments
)
# Remove any NA entries so unspecified experiments fall back to all IVs
focal_iv_map <- focal_iv_map[!is.na(focal_iv_map)]

lm_plot_dir <- file.path(output_folder, "p_new_param_plots")

message("\n===== LINEAR MODELS BY EXPERIMENT (Magnitude Accuracy) =====")
lm_magn_by_exp <- lm_by_experiment(df_analysis, "abs_error", my_treatments,
                                    focal_iv_map = focal_iv_map,
                                    plot_dir = lm_plot_dir)
summary(df_analysis)
df_analysis$trueSpeciesRichness


message("\n===== LINEAR MODELS BY EXPERIMENT (Proportional Accuracy) =====")
lm_prop_by_exp <- lm_by_experiment(df_analysis, "abs_prop_error", my_treatments,
                                    focal_iv_map = focal_iv_map,
                                    plot_dir = lm_plot_dir)

# --- Summary: adj R2 per experiment ---
message("\n--- Model fit summary (adj. R-squared) ---")
if (!is.null(lm_magn_by_exp)) {
  r2_magn <- lm_magn_by_exp |>
    dplyr::distinct(Experiment, R2_adj) |>
    dplyr::mutate(DV = "abs_error")
}
if (!is.null(lm_prop_by_exp)) {
  r2_prop <- lm_prop_by_exp |>
    dplyr::distinct(Experiment, R2_adj) |>
    dplyr::mutate(DV = "abs_prop_error")
}
r2_summary <- rbind(r2_magn, r2_prop)
print(r2_summary)

# --- Visualization: top significant predictors across experiments ---
lm_all <- rbind(lm_magn_by_exp, lm_prop_by_exp)
lm_sig <- lm_all[lm_all$P_Value < 0.05, ]

if (nrow(lm_sig) > 0) {
  # Count how many experiments each IV is significant in, per DV
  sig_counts <- lm_sig |>
    dplyr::count(DV, IV, name = "N_Sig_Experiments") |>
    dplyr::arrange(DV, dplyr::desc(N_Sig_Experiments))
  
  ggplot(sig_counts, aes(x = N_Sig_Experiments,
                          y = reorder(IV, N_Sig_Experiments))) +
    geom_col() +
    facet_wrap(~ DV, scales = "free_x") +
    labs(x = "Number of experiments where p < 0.05",
         y = NULL,
         title = "Consistently significant predictors of accuracy across experiments") +
    theme_minimal(base_size = 12)
}

# =============================================================================
# STEP 5b: PER-VARIABLE PLOTS (all IVs, pooled across tests)
# =============================================================================
# Generates the same style of regression/boxplot as the per-experiment LM plots,
# but pooling all data to show overall effect of each IV.

p_iv_plot_dir <- file.path(output_folder, "p_iv_plots")
dir.create(p_iv_plot_dir, showWarnings = FALSE, recursive = TRUE)

for (dv in c("abs_error", "abs_prop_error")) {
  
  iv_plots <- list()
  
  for (iv in my_treatments) {
    if (!iv %in% names(df_analysis)) next
    # Balance new-param axis comparisons; simulation IVs stay pooled
    iv_data <- balance_focal(df_analysis, iv, verbose = FALSE)
    vals <- na.omit(iv_data[[iv]])
    if (length(unique(vals)) < 2) next
    
    is_cat <- !is.numeric(iv_data[[iv]])
    
    p <- ggplot(iv_data, aes(x = .data[[iv]], y = .data[[dv]]))
    
    if (is_cat) {
      p <- p +
        geom_boxplot(outlier.shape = NA) +
        geom_jitter(aes(fill = estimate_direction), shape = 21, width = 0.2,
                    alpha = 0.6, color = "black", size = 2) +
        stat_summary(fun = mean, geom = "crossbar", width = 0.4,
                     color = "steelblue", linewidth = 0.8) +
        scale_fill_manual(values = c("overestimate" = "white", "underestimate" = "grey60"),
                          name = "Direction")
    } else {
      p <- p +
        geom_point(aes(fill = estimate_direction), shape = 21,
                   alpha = 0.6, color = "black", size = 2) +
        scale_fill_manual(values = c("overestimate" = "white", "underestimate" = "grey60"),
                          name = "Direction") +
        geom_smooth(method = "lm", se = TRUE, color = "steelblue", formula = y ~ x)
    }
    
    p <- p +
      labs(x = iv, y = dv) +
      theme_minimal(base_size = 11) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    iv_plots[[iv]] <- p
  }
  
  if (length(iv_plots) > 0) {
    library(patchwork)
    n_plots <- length(iv_plots)
    ncols <- min(n_plots, 3)
    nrows <- ceiling(n_plots / ncols)
    combined <- wrap_plots(iv_plots, ncol = ncols) +
      plot_annotation(
        title = paste0("All tests pooled  |  DV: ", dv),
        caption = paste0("n = ", nrow(df_analysis)),
        theme = theme_minimal(base_size = 13)
      )
    fname <- file.path(p_iv_plot_dir, paste0("pooled_all_IVs_", dv, ".svg"))
    ggsave(fname, combined,
           width = 4 * ncols, height = 4 * nrows, limitsize = FALSE)
    message("Saved: ", fname)
    print(combined)
  }
}

# =============================================================================
# STEP 6: ANALYSIS OF NEW PARAMETERS (fill, filter, fit_params)
# =============================================================================



# Categorical new parameters -> boxplot/jitter treatment (loop 6c).
new_param_vars <- c("fill_plants", "fill_animals", "filter_plant_completeness", "fit_params")

# Continuous new parameters -> linear-model treatment (loop 6c-cont).
continuous_param_vars <- c("comp_mean", "comp_sd", "comp_spread", "completeness_floor", "n_cells_filtered" )

# --- 6a. Success/failure rates by new parameter combinations ---
fill_filter_vars <- c("fill_plants", "fill_animals", "filter_plant_completeness")
failure_summary <- df[, .(
  total = .N,
  succeeded = sum(!is.na(best_estimate)),
  failed = sum(is.na(best_estimate))
), by = fill_filter_vars]
message("\n===== FILL / FILTER SUCCESS RATES =====")
print(failure_summary)

# --- 6b. Create organism group variable ---
# Organism is encoded in the scenario prefix (e.g. "blatt-roachLike"):
#   blatt -> roach, coleo -> beetle, butter -> butterfly.
# variedSiteEvenness runs (even/low/veryLow prefixes) are an experiment type,
# not an organism, so they are kept as their own group.
df_analysis[, organism_group := data.table::fcase(
  grepl("variedSiteEvenness", scenario), "variedSiteEvenness",
  grepl("^blatt",  scenario),            "roach",
  grepl("^coleo",  scenario),            "beetle",
  grepl("^butter", scenario),            "butterfly",
  default = NA_character_
)]
message("Organism groups: ", paste(unique(df_analysis$organism_group), collapse = ", "))

# --- 6c. Plots for each new parameter (saved to SVG) ---
p_plot_dir <- file.path(output_folder, "p_new_param_plots")
dir.create(p_plot_dir, showWarnings = FALSE, recursive = TRUE)

# Shared aesthetics matching lm_by_experiment plots
point_shape_aes <- aes(fill = estimate_direction)
direction_scale <- scale_fill_manual(
  values = c("overestimate" = "white", "underestimate" = "grey60"),
  name = "Direction"
)

# --- 6c-pre. Create fill_plants x fill_animals interaction variable ---
# Encoding: first letter = fill_plants (P), second letter = fill_animals (A)
# e.g. "TF" means fill_plants=TRUE, fill_animals=FALSE
df_analysis[, fill_interaction := factor(
  paste0(
    ifelse(fill_plants  == "TRUE", "T", "F"),
    ifelse(fill_animals == "TRUE", "T", "F")
  ),
  levels = c("FF", "FT", "TF", "TT"),
  labels = c(
    "FF\n(Plants=F, Animals=F)",
    "FT\n(Plants=F, Animals=T)",
    "TF\n(Plants=T, Animals=F)",
    "TT\n(Plants=T, Animals=T)"
  )
)]
message("fill_interaction distribution:\n",
        paste(capture.output(table(df_analysis$fill_interaction)), collapse = "\n"))

# Append to param loop so it gets the same accuracy + bias plots
new_param_vars <- c(new_param_vars, "fill_interaction")

for (param in new_param_vars) {
  if (!param %in% names(df_analysis) || length(unique(na.omit(df_analysis[[param]]))) < 2) {
    message("  Skipping ", param, " (no variance in analysis data)")
    next
  }

  # Balance: hold the other new-param axes at baseline for this focal axis
  param_data <- balance_focal(df_analysis, param, verbose = TRUE)

  # Accuracy
  p_acc <- ggplot(param_data, aes(x = .data[[param]], y = abs_prop_error)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(point_shape_aes, shape = 21, width = 0.2, alpha = 0.6,
                color = "black", size = 2) +
    direction_scale +
    stat_summary(fun = mean, geom = "crossbar", width = 0.4,
                 color = "steelblue", linewidth = 0.8) +
    labs(x = param, y = "Absolute Proportional Error",
         title = paste("Accuracy by", param)) +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(file.path(p_plot_dir, paste0(param, "_accuracy.svg")), p_acc,
         width = 6, height = 5)
  print(p_acc)
  
  # Bias
  p_bias <- ggplot(param_data, aes(x = .data[[param]], y = bias_prop)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(point_shape_aes, shape = 21, width = 0.2, alpha = 0.6,
                color = "black", size = 2) +
    direction_scale +
    stat_summary(fun = mean, geom = "crossbar", width = 0.4,
                 color = "steelblue", linewidth = 0.8) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
    labs(x = param, y = "Proportional Bias (estimate / truth)",
         title = paste("Bias by", param),
         subtitle = "Red line = no bias") +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(file.path(p_plot_dir, paste0(param, "_bias.svg")), p_bias,
         width = 6, height = 5)
  print(p_bias)
  
  # Faceted by organism group
  n_groups <- length(unique(param_data$organism_group))
  p_facet <- ggplot(param_data, aes(x = .data[[param]], y = bias_prop)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(point_shape_aes, shape = 21, width = 0.2, alpha = 0.6,
                color = "black", size = 2) +
    direction_scale +
    stat_summary(fun = mean, geom = "crossbar", width = 0.4,
                 color = "steelblue", linewidth = 0.8) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
    facet_wrap(~ organism_group, scales = "free_y") +
    labs(x = param, y = "Proportional Bias",
         title = paste("Bias by", param, "across organism groups")) +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(file.path(p_plot_dir, paste0(param, "_bias_by_organism.svg")), p_facet,
         width = 4 * min(n_groups, 3), height = 4 * ceiling(n_groups / 3),
         limitsize = FALSE)
  print(p_facet)
  
  message("  Saved plots for: ", param)
}

# --- 6c-cont. Plots for continuous new parameters (linear models) ---
# Same panels as 6c (accuracy, bias, bias-by-organism) but for continuous IVs:
# scatter + linear fit, styled to match the continuous IV plots from Step 5.
# Each plot's subtitle reports the bivariate lm slope, R^2, and slope p-value.
# 
# 

#NOTE: The beetle data look like they have high completeness

for (param in continuous_param_vars) {
  if (!param %in% names(df_analysis) ||
      length(unique(na.omit(df_analysis[[param]]))) < 2) {
    message("  Skipping ", param, " (no variance in analysis data)")
    next
  }

  #' Bivariate lm stats for a subtitle string (slope, R^2, slope p-value)
  lm_subtitle <- function(x, y) {
    fit <- lm(y ~ x)
    sm  <- summary(fit)
    paste0("slope = ", signif(coef(fit)[2], 3),
           ",  R2 = ", round(sm$r.squared, 3),
           ",  p = ", signif(coef(sm)[2, 4], 2))
  }

  # Accuracy
  p_acc <- ggplot(df_analysis, aes(x = .data[[param]], y = abs_prop_error)) +
    geom_point(point_shape_aes, shape = 21, alpha = 0.6,
               color = "black", size = 2) +
    direction_scale +
    geom_smooth(method = "lm", se = TRUE, color = "steelblue", formula = y ~ x) +
    labs(x = param, y = "Absolute Proportional Error",
         title = paste("Accuracy by", param),
         subtitle = lm_subtitle(df_analysis[[param]], df_analysis$abs_prop_error)) +
    theme_minimal(base_size = 11)
  ggsave(file.path(p_plot_dir, paste0(param, "_accuracy.svg")), p_acc,
         width = 6, height = 5)
  print(p_acc)

  # Bias
  p_bias <- ggplot(df_analysis, aes(x = .data[[param]], y = bias_prop)) +
    geom_point(point_shape_aes, shape = 21, alpha = 0.6,
               color = "black", size = 2) +
    direction_scale +
    geom_smooth(method = "lm", se = TRUE, color = "steelblue", formula = y ~ x) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
    labs(x = param, y = "Proportional Bias (estimate / truth)",
         title = paste("Bias by", param),
         subtitle = paste0(lm_subtitle(df_analysis[[param]], df_analysis$bias_prop),
                           "   |   Red line = no bias")) +
    theme_minimal(base_size = 11)
  ggsave(file.path(p_plot_dir, paste0(param, "_bias.svg")), p_bias,
         width = 6, height = 5)
  print(p_bias)

  # Faceted by organism group (linear fit per group)
  n_groups <- length(unique(df_analysis$organism_group))

  # Per-facet lm stats: fit bias_prop ~ param within each organism group.
  facet_stats <- do.call(rbind, lapply(
    split(df_analysis, df_analysis$organism_group),
    function(d) {
      d <- d[!is.na(d[[param]]) & !is.na(d$bias_prop), ]
      # Need >= 3 points and variation in x to fit a meaningful line
      if (nrow(d) < 3 || length(unique(d[[param]])) < 2) return(NULL)
      fit <- lm(d$bias_prop ~ d[[param]])
      sm  <- summary(fit)
      data.frame(
        organism_group = d$organism_group[1],
        label = paste0("slope = ", signif(coef(fit)[2], 3),
                       "\nR2 = ", round(sm$r.squared, 3),
                       "\np = ", signif(coef(sm)[2, 4], 2))
      )
    }
  ))

  p_facet <- ggplot(df_analysis, aes(x = .data[[param]], y = bias_prop)) +
    geom_point(point_shape_aes, shape = 21, alpha = 0.6,
               color = "black", size = 2) +
    direction_scale +
    geom_smooth(method = "lm", se = TRUE, color = "steelblue", formula = y ~ x) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
    geom_text(data = facet_stats,
              aes(x = -Inf, y = Inf, label = label),
              hjust = -0.05, vjust = 1.1, size = 3, color = "grey20",
              inherit.aes = FALSE) +
    facet_wrap(~ organism_group, scales = "free_y") +
    labs(x = param, y = "Proportional Bias",
         title = paste("Bias by", param, "across organism groups")) +
    theme_minimal(base_size = 11)
  ggsave(file.path(p_plot_dir, paste0(param, "_bias_by_organism.svg")), p_facet,
         width = 4 * min(n_groups, 3), height = 4 * ceiling(n_groups / 3),
         limitsize = FALSE)
  print(p_facet)

  message("  Saved continuous-parameter plots for: ", param)
}

# --- 6d. Kruskal-Wallis tests for each new parameter ---
message("\n===== STATISTICAL TESTS FOR NEW PARAMETERS =====")
for (param in new_param_vars) {
  if (!param %in% names(df_analysis) || length(unique(na.omit(df_analysis[[param]]))) < 2) next
  # Balance the focal axis (hold other new-param axes at baseline) before testing
  kw_data <- balance_focal(df_analysis, param, verbose = FALSE)
  kw <- kruskal.test(as.formula(paste("abs_prop_error ~", param)), data = kw_data)
  message(param, ": chi-sq = ", round(kw$statistic, 2),
          ", df = ", kw$parameter, ", p = ", signif(kw$p.value, 3),
          ", n = ", nrow(kw_data))
}

# --- 6e. Organism-by-organism best parameter combination ---
message("\n===== BEST PARAMETER COMBO BY ORGANISM GROUP =====")

# Summarize accuracy for every combination of new params within each organism group
combo_summary <- df_analysis[, .(
  N = .N,
  mean_abs_prop_error = mean(abs_prop_error, na.rm = TRUE),
  median_abs_prop_error = median(abs_prop_error, na.rm = TRUE),
  mean_bias_prop = mean(bias_prop, na.rm = TRUE),
  median_bias_prop = median(bias_prop, na.rm = TRUE),
  mean_abs_error = mean(abs_error, na.rm = TRUE)
), by = c("organism_group", new_param_vars)]
setorder(combo_summary, organism_group, mean_abs_prop_error)

# Best combo per organism group (lowest mean absolute proportional error)
best_combos <- combo_summary[, .SD[1], by = organism_group]
message("Best parameter combinations (by mean absolute proportional error):")
print(best_combos)

# Full table for inspection
message("\nFull combo summary (sorted by organism group, then accuracy):")
print(combo_summary)

# Save combo summary
fwrite(combo_summary, file.path(output_folder, "p_organism_param_combo_summary.csv"))
fwrite(best_combos, file.path(output_folder, "p_best_combos_by_organism.csv"))
message("Saved: p_organism_param_combo_summary.csv, p_best_combos_by_organism.csv")

# --- 6f. Heatmap of accuracy by parameter combo per organism group ---
# Create a label for each parameter combination
combo_summary[, combo_label := paste(
  ifelse(fill_plants == "TRUE", "FillPl", ""),
  ifelse(fill_animals == "TRUE", "FillAn", ""),
  ifelse(filter_plant_completeness != "none", paste0("Filt", filter_plant_completeness), ""),
  fit_params,
  sep = " | "
)]
combo_summary[, combo_label := gsub("^(\\s*\\|\\s*)+", "", combo_label)]
combo_summary[, combo_label := gsub("(\\s*\\|\\s*)+$", "", combo_label)]

p_heatmap <- ggplot(combo_summary, aes(x = organism_group, y = combo_label,
                                       fill = mean_abs_prop_error)) +
  geom_tile() +
  geom_text(aes(label = round(mean_abs_prop_error, 2)), size = 3) +
  scale_fill_viridis_c(direction = -1, name = "Mean |Prop Error|") +
  labs(x = "Organism Group", y = "Parameter Combination",
       title = "Accuracy by Parameter Combination and Organism Group") +
  theme_minimal(base_size = 11) +
  theme(axis.text.y = element_text(size = 8),
        axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(p_plot_dir, "param_combo_heatmap.svg"), p_heatmap,
       width = 8, height = 6 + 0.25 * length(unique(combo_summary$combo_label)),
       limitsize = FALSE)
print(p_heatmap)


# =============================================================================
# STEP 7: EXPORT OBJECTS FOR DOWNSTREAM NOTEBOOKS ----
# =============================================================================
# Serialises all key analysis objects produced above so that Notebooks 6 and 7
# can load them directly without re-running the full pipeline.
#
# Exported files (all land in output_folder = 5-SensitivityTesting/results/):
#   p_df_analysis.csv            — full analysis dataset (all IVs + DVs)
#   p_vi_magn.csv / p_vi_prop.csv — variable importance tables (partial eta-sq)
#   p_lm_all.csv / p_lm_sig.csv  — per-experiment LM results (all / p<0.05 only)
#   p_sig_counts.csv             — count of experiments each IV is significant in
#   p_interaction_results.csv    — top-6-IV pairwise interaction test results
#   p_magn_summary_bylevel.csv   — per-level magnitude bias summary
#   p_prop_summary_bylevel.csv   — per-level proportional bias summary
# =============================================================================
{
  fwrite(df_analysis, file.path(output_folder, "p_df_analysis.csv"))

  if (!is.null(vi_magn))
    write.csv(vi_magn,  file.path(output_folder, "p_vi_magn.csv"),  row.names = FALSE)
  if (!is.null(vi_prop))
    write.csv(vi_prop,  file.path(output_folder, "p_vi_prop.csv"),  row.names = FALSE)

  if (!is.null(lm_all))
    write.csv(lm_all,  file.path(output_folder, "p_lm_all.csv"),  row.names = FALSE)
  if (!is.null(lm_sig))
    write.csv(lm_sig,  file.path(output_folder, "p_lm_sig.csv"),  row.names = FALSE)
  if (exists("sig_counts"))
    write.csv(sig_counts, file.path(output_folder, "p_sig_counts.csv"), row.names = FALSE)

  if (exists("interaction_results") && !is.null(interaction_results))
    write.csv(interaction_results,
              file.path(output_folder, "p_interaction_results.csv"), row.names = FALSE)

  if (!is.null(magn_summary) && !is.null(magn_summary$by_level))
    write.csv(magn_summary$by_level,
              file.path(output_folder, "p_magn_summary_bylevel.csv"), row.names = FALSE)
  if (!is.null(prop_summary) && !is.null(prop_summary$by_level))
    write.csv(prop_summary$by_level,
              file.path(output_folder, "p_prop_summary_bylevel.csv"), row.names = FALSE)

  message("Step 7 export complete → ", output_folder)
}








