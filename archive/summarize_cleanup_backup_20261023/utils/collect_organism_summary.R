# utils/collect_organism_summary.R
# Sourced by AmAni.summarize.R — relies on the packages (tidyverse, sf,
# patchwork, ggrepel) loaded there, on organism_label() from
# utils/organism_label.R, and on the private helpers in
# utils/pub_data_helpers.R.
#
# Combines, per organism: the before/after cell-filtering stats from
# dataset_summary_*.txt, gamma richness estimates (total and surrogate-only
# models) from gamma_richnessfromSPRICH*.txt, and the sampling-extension
# diagnostics from sampling_extension_report_*.txt. Organisms without a given
# report simply get NA in the corresponding columns. Parsing logic mirrors
# summarize_all() so the two stay consistent.

collect_organism_summary <- function(organisms = NULL, output_dir = "output") {
  if (is.null(organisms)) organisms <- .organism_dirs(output_dir)

  purrr::map_dfr(organisms, function(org) {
    am_dir <- file.path(output_dir, org, "animal_model")
    if (!dir.exists(am_dir)) return(NULL)

    # -- Dataset summary (filter stats) --
    ds_file <- .latest_file(am_dir, "^dataset_summary_.*\\.txt$")
    if (!is.null(ds_file)) {
      ds <- readLines(ds_file)

      pre_cells   <- .extract_val(ds, "^\\s*Grid cells:")
      pre_records <- .extract_val(ds, "^\\s*Total records:")
      pre_species <- .extract_val(ds, "^\\s*Total species:")

      # After-filtering section, if a filter was applied
      after_idx <- grep("After filtering", ds)
      if (length(after_idx) > 0) {
        ds_after <- ds[after_idx[1]:length(ds)]
        post_cells   <- .extract_val(ds_after, "^\\s*Grid cells:")
        post_records <- .extract_val(ds_after, "^\\s*Total records:")
        post_species <- .extract_val(ds_after, "^\\s*Total species:")
      } else {
        post_cells <- post_records <- post_species <- NA_real_
      }

      filter_applied <- any(grepl("Filter applied:\\s*YES", ds))
      min_occ  <- .extract_val(ds, "min_occurrences:")
      min_comp <- .extract_val(ds, "min_completeness:")
      min_rich <- .extract_val(ds, "min_richness:")
      outlier_cutoff <- .extract_val(ds, "outlier cutoff:", as_type = "character")
      track_b_entry  <- .extract_val(ds, "Track B entry:", as_type = "character")

      # "  Grid cells:  117 (117 well-sampled + 55 super-diverse; 2 outliers dropped)"
      gc_line <- grep("^\\s*Grid cells:", ds_after, value = TRUE)
      if (length(gc_line) > 0) {
        grab <- function(re) {
          m <- regmatches(gc_line[1], regexpr(re, gc_line[1]))
          if (length(m)) as.integer(sub("\\s.*", "", m)) else NA_integer_
        }
        well_sampled    <- grab("([0-9]+)\\s+well-sampled")
        super_diverse   <- grab("([0-9]+)\\s+super-diverse")
        outliers_dropped <- grab("([0-9]+)\\s+outliers? dropped")
      } else {
        well_sampled <- super_diverse <- outliers_dropped <- NA_integer_
      }

      # "  Species lost:  25 (14.1%)" — keep both count and percentage
      sl <- .extract_val(ds, "^\\s*Species lost:", as_type = "character")
      filter_species_lost <- .extract_val(ds, "^\\s*Species lost:")
      filter_species_lost_pct <- {
        m <- regmatches(sl, regexpr("[0-9]+\\.?[0-9]*(?=%)", sl, perl = TRUE))
        if (length(m)) as.numeric(m) else NA_real_
      }
    } else {
      pre_cells <- pre_records <- pre_species <- NA_real_
      post_cells <- post_records <- post_species <- NA_real_
      filter_applied <- FALSE
      min_occ <- min_comp <- min_rich <- NA_real_
      outlier_cutoff <- track_b_entry <- NA_character_
      well_sampled <- super_diverse <- outliers_dropped <- NA_integer_
      filter_species_lost <- filter_species_lost_pct <- NA_real_
    }

    # -- Gamma richness --
    gm_file <- .latest_file(am_dir, "^gamma_richnessfromSPRICH.*\\.txt$")
    if (!is.null(gm_file)) {
      gm <- readLines(gm_file)
      scaling_factor <- .extract_val(gm, "Scaling factor:")
      fit_params <- .extract_val(gm, "Fit parameters:", as_type = "character")

      delta_str  <- .extract_val(gm, "^delta_meanlog:", as_type = "character")
      delta_meanlog <- .extract_val(gm, "^delta_meanlog:")
      delta_active <- grepl("active", delta_str, ignore.case = TRUE)
      ratio_str  <- .extract_val(gm, "^ratio_sdlog:", as_type = "character")
      ratio_sdlog <- .extract_val(gm, "^ratio_sdlog:")
      ratio_active <- grepl("active", ratio_str, ignore.case = TRUE)

      # New format: total model (surrogate + extension) is primary
      gamma_best  <- .extract_val(gm, "^Best estimate from total model")
      gamma_upper <- .extract_val(gm, "^Upper bound:\\s")
      gamma_lower <- .extract_val(gm, "^Lower bound:\\s")

      # Fall back to old format if total-model fields are missing
      if (is.na(gamma_best)) {
        gamma_best  <- .extract_val(gm, "Best estimate:")
        gamma_upper <- .extract_val(gm, "Upper bound:")
        gamma_lower <- .extract_val(gm, "Lower bound:")
      }

      surrogate_best  <- .extract_val(gm, "^Best estimate from plant surrogate")
      surrogate_upper <- .extract_val(gm, "^Upper bound \\(narrow")
      surrogate_lower <- .extract_val(gm, "^Lower bound \\(narrow")
      richness_shortfall <- .extract_val(gm, "^Expected proportional richness shortfall")
    } else {
      scaling_factor <- gamma_best <- gamma_upper <- gamma_lower <- NA_real_
      surrogate_best <- surrogate_upper <- surrogate_lower <- NA_real_
      richness_shortfall <- NA_real_
      fit_params <- NA_character_
      delta_meanlog <- ratio_sdlog <- NA_real_
      delta_active <- ratio_active <- FALSE
    }

    # -- Sampling extension diagnostics --
    se_file <- .latest_file(am_dir, "^sampling_extension_report_.*\\.txt$")
    if (!is.null(se_file)) {
      se <- readLines(se_file)
      se_chao2_s_obs   <- .extract_val(se, "S_obs =")
      se_chao2_s_chao2 <- .extract_val(se, "S_Chao2 =")
      se_mean_range    <- .extract_val(se, "sample_mean_range_size:")
      se_partial_gamma <- .extract_val(se, "partialGamma:")
      se_abs_prop_error <- .extract_val(se, "abs_prop_error:")
      se_correction_factor <- .extract_val(se, "correction_factor:")
      se_full_gamma_best  <- .extract_val(se, "fullGamma \\(best\\):")
      se_full_gamma_lower <- .extract_val(se, "fullGamma \\(lower")
      se_full_gamma_upper <- .extract_val(se, "fullGamma \\(upper")
    } else {
      se_chao2_s_obs <- se_chao2_s_chao2 <- se_mean_range <- NA_real_
      se_partial_gamma <- se_abs_prop_error <- se_correction_factor <- NA_real_
      se_full_gamma_best <- se_full_gamma_lower <- se_full_gamma_upper <- NA_real_
    }

    # Observed species: post-filter if a filter was applied, else pre-filter
    observed_spp <- if (!is.na(post_species) && filter_applied) post_species else pre_species

    tibble::tibble(
      organism_dir          = org,
      organism              = organism_label(org),
      # --- Pre-filter ---
      pre_cells             = pre_cells,
      pre_records           = pre_records,
      pre_species           = pre_species,
      # --- Filter settings ---
      filter_applied        = filter_applied,
      min_occurrences       = min_occ,
      min_completeness      = min_comp,
      min_richness          = min_rich,
      outlier_cutoff        = outlier_cutoff,
      track_b_entry         = track_b_entry,
      # --- Post-filter ---
      post_cells            = post_cells,
      post_records          = post_records,
      post_species          = post_species,
      well_sampled          = well_sampled,
      super_diverse         = super_diverse,
      outliers_dropped      = outliers_dropped,
      filter_species_lost   = filter_species_lost,
      filter_species_lost_pct = filter_species_lost_pct,
      # --- Gamma richness ---
      observed_spp          = observed_spp,
      gamma_best            = gamma_best,
      gamma_lower           = gamma_lower,
      gamma_upper           = gamma_upper,
      surrogate_best        = surrogate_best,
      surrogate_lower       = surrogate_lower,
      surrogate_upper       = surrogate_upper,
      richness_shortfall    = richness_shortfall,
      gamma_difference      = gamma_best - pre_species,
      gamma_ratio           = round(gamma_best / pre_species, 2),
      # --- Model internals ---
      scaling_factor        = scaling_factor,
      fit_params            = fit_params,
      delta_meanlog         = delta_meanlog,
      delta_meanlog_active  = delta_active,
      ratio_sdlog           = ratio_sdlog,
      ratio_sdlog_active    = ratio_active,
      # --- Sampling extension diagnostics ---
      se_chao2_s_obs        = se_chao2_s_obs,
      se_chao2_s_chao2      = se_chao2_s_chao2,
      se_mean_range_size    = se_mean_range,
      se_partial_gamma      = se_partial_gamma,
      se_abs_prop_error     = se_abs_prop_error,
      se_correction_factor  = se_correction_factor,
      se_full_gamma_best    = se_full_gamma_best,
      se_full_gamma_lower   = se_full_gamma_lower,
      se_full_gamma_upper   = se_full_gamma_upper
    )
  })
}
