
source("scripts/utils/helpers.R")

fillingScriptDir <- paste0(dirname(getwd()), "/1-FilteringPlantData/scripts")
source(paste0(fillingScriptDir, "/computeCompleteness.R"))
source(paste0(fillingScriptDir, "/fillingRanges.R"))

modelAnimal<-function(
    animal_occurrenceData,
    fill_animals = FALSE,
    alphaRichnessMethod = "empirical",
  diagnostic =FALSE, #Set to TRUE if you want animalAnalysis to check that it's working
  spatial = "explicit", #If set to NULL, species will be simulated in any grid cell. The default option prefers grid cells where there are a lot of plant data
  redoeverything = FALSE, #if TRUE it will recalculate raw endemism and plant_richness. The latter takes a while
  customScale = NULL, #Alternatively, set to a numeric value between 0 and 1. This will skip calculation of upper and lower bounds and should only be used for troubleshooting. The default value will calculate the scaling factor for you (which is the crux of the analysis)
  fit_params = c("amplitude"), #Which log-normal parameters to estimate from the animal data.
  # Subset of c("amplitude", "meanlog", "sdlog"):
  #   "amplitude" — scale total species count (always included; existing behaviour)
  #   "meanlog"   — shift the distribution's location (range sizes)
  #   "sdlog"     — rescale the distribution's spread (endemic-vs-widespread ratio)
  # Parameters not listed are inherited from the plant fit.
  # Examples:
  #   c("amplitude")                   — current default (amplitude only)
  #   c("amplitude", "meanlog")        — amplitude + location shift
  #   c("amplitude", "meanlog", "sdlog") — fit all three parameters
  estimateBounds = .amani_estimateBounds #TRUE or FALSE...FALSE speeds up everything 3-fold but then no error ranges will be estimated
  
){
  # Clear stale objects from previous organisms so nothing carries over
  # plant_richness is NOT cleared here — it's organism-independent and expensive
  # to compute. It IS plant-data-type dependent (filled vs unfilled), so
  # compute_plant_richness() guards reuse with .plant_richness_fill_state and
  # reloads/recomputes from the correct cache when fill_plants changes.
  stale <- c("scaled_results", "animal_endemism_filtered", "overlap_bbox",
             "animal_endemism_grid", "animal_endemism",
             "combined_endemism", "scaling_factors", "ratio_median",
             "fit_pred_shifted", "lognorm_fit_animal",
             ".amani_delta_meanlog", ".amani_ratio_sdlog", ".amani_plant_mdd_meanlog",
             ".amani_filled_species")
  for (v in stale) if (exists(v, envir = .GlobalEnv)) rm(list = v, envir = .GlobalEnv)

  # Push parameters to global env so sourced scripts can see them
  animal_grid_data <<- animal_occurrenceData
  diagnostic <<- diagnostic
  spatial <<- spatial
  redoeverything <<- redoeverything
  customScale <<- customScale
  .amani_fit_params <<- fit_params
  .amani_estimateBounds <<- estimateBounds

  # ── Combined cell filter ─────────────────────────────────────────────
  #
  # Two tracks of cells are kept:
  #
  # Track A — "Well-sampled": cells with adequate sampling completeness.
  #   Requires: Chao1 completeness >= min_completeness,
  #             occurrences      >= min_occurrences,
  #             observed species >= min_richness,
  #             AND richness     <= outlier cutoff (see max_richness below).
  #
  # Track B — "Super-diverse": cells with very high richness that may be
  #   under-sampled but are too important to discard.
  #   Requires: richness > 97th percentile AND <= outlier cutoff.
  #   Cells above the outlier cutoff are treated as true outliers
  #   (likely data artefacts) and dropped.
  #
  # Filter parameters:
  #   min_completeness  Chao1 S_obs / S_est. Set to "auto" to estimate from
  #                     the data, or a numeric value to override. Default "auto".
  #   min_occurrences   Minimum occurrences per cell. Set to "auto" to
  #                     estimate from the discovery rate curve, or numeric to
  #                     override. Default "auto".
  #   min_richness      Minimum observed species. This should be set manually
  #                     as a very conservative estimate of the minimum species
  #                     count one would expect in the least diverse region of
  #                     the Amazon, including deforested areas. Default 20.
  #   max_richness      Maximum plausible species richness per grid cell.
  #                     If NULL (default), the 99.9th percentile of observed
  #                     richness is used. If a positive number, that value is
  #                     used directly as the outlier cutoff — any cell with
  #                     more species is dropped. This is useful when you have
  #                     a domain-informed upper bound on how many species
  #                     could realistically occur in a single grid cell.
  # ─────────────────────────────────────────────────────────────────────
  run_filter <- if (exists(".amani_min_occurrences", envir = .GlobalEnv)) {
    # "auto" or a positive number both trigger filtering
    val <- .amani_min_occurrences
    !is.null(val) && (identical(val, "auto") || (is.numeric(val) && val > 0))
  } else {
    FALSE
  }

  if (run_filter && exists("animal_grid", envir = .GlobalEnv)) {

    # --- Per-cell statistics (computed first; needed by auto-estimation) ---
    spp_abund <- animal_grid |>
      sf::st_drop_geometry() |>
      dplyr::count(grid_id, Morphotype, name = "abundance")

    cell_stats <- spp_abund |>
      dplyr::group_by(grid_id) |>
      dplyr::summarise(
        n_occurrences = sum(abundance),
        n_species     = dplyr::n(),
        f1 = sum(abundance == 1),
        f2 = sum(abundance == 2),
        .groups = "drop"
      ) |>
      dplyr::mutate(
        chao1 = n_species + dplyr::if_else(
          f2 > 0, f1^2 / (2 * f2), f1 * (f1 - 1) / 2
        ),
        completeness = n_species / chao1
      )

    n_before <- nrow(cell_stats)

    # --- Resolve thresholds (auto or manual) ---
    raw_occ  <- if (exists(".amani_min_occurrences",  envir = .GlobalEnv)) .amani_min_occurrences  else "auto"
    raw_comp <- if (exists(".amani_min_completeness", envir = .GlobalEnv)) .amani_min_completeness else "auto"

    if (identical(raw_occ, "auto") || identical(raw_comp, "auto")) {
      auto <- estimate_filter_thresholds(cell_stats)
    }

    min_occurrences  <- if (identical(raw_occ, "auto"))  auto$min_occurrences  else as.numeric(raw_occ)
    min_completeness <- if (identical(raw_comp, "auto")) auto$min_completeness else as.numeric(raw_comp)
    min_richness     <- if (exists(".amani_min_richness", envir = .GlobalEnv)) .amani_min_richness else 20

    # --- Richness fences ---
    credible_for_fence <- cell_stats |> dplyr::filter(n_species >= 3)
    upper_fence_97 <- quantile(credible_for_fence$n_species, 0.97)

    # Outlier cutoff: disabled, user-supplied max_richness, or 99.9th percentile
    disable_outlier <- isTRUE(
      if (exists(".amani_disable_outlier_filter", envir = .GlobalEnv))
        .amani_disable_outlier_filter else FALSE
    )
    raw_max_rich <- if (exists(".amani_max_richness", envir = .GlobalEnv)) .amani_max_richness else NULL
    if (disable_outlier) {
      outlier_cutoff <- Inf
      outlier_source <- "DISABLED"
      message("Outlier site filter disabled")
    } else if (!is.null(raw_max_rich) && is.numeric(raw_max_rich) && raw_max_rich > 0) {
      outlier_cutoff <- raw_max_rich
      outlier_source <- paste0("user-specified (", outlier_cutoff, " spp)")
    } else {
      outlier_cutoff <- quantile(credible_for_fence$n_species, 0.999)
      outlier_source <- paste0("99.9th percentile (", round(outlier_cutoff), " spp)")
    }

    # --- Track A: well-sampled cells ---
    well_sampled_ids <- cell_stats$grid_id[
      cell_stats$completeness   >= min_completeness &
      cell_stats$n_occurrences  >= min_occurrences  &
      cell_stats$n_species      >= min_richness     &
      cell_stats$n_species      <= outlier_cutoff
    ]

    # --- Track B: super-diverse cells (above 97th, below outlier cutoff) ---
    super_diverse_ids <- cell_stats$grid_id[
      cell_stats$n_species > upper_fence_97 &
      cell_stats$n_species <= outlier_cutoff
    ]

    # --- True outliers (above outlier cutoff) — dropped ---
    true_outlier_ids <- cell_stats$grid_id[
      cell_stats$n_species > outlier_cutoff
    ]

    keep_cells <- union(well_sampled_ids, super_diverse_ids)

    # --- Save pre-filter stats for the downstream report ---
    kept_stats <- cell_stats |> dplyr::filter(grid_id %in% keep_cells)
    .amani_filter_summary <<- list(
      pre_filter = list(
        n_cells       = n_before,
        n_records     = sum(cell_stats$n_occurrences),
        n_species     = length(unique(spp_abund$Morphotype)),
        richness_summary  = summary(cell_stats$n_species),
        occ_summary       = summary(cell_stats$n_occurrences),
        comp_summary      = summary(cell_stats$completeness)
      ),
      post_filter = list(
        n_cells       = length(keep_cells),
        n_well_sampled = length(well_sampled_ids),
        n_super_diverse = length(super_diverse_ids),
        n_outliers    = length(true_outlier_ids),
        n_records     = sum(kept_stats$n_occurrences),
        n_species     = length(unique(spp_abund$Morphotype[spp_abund$grid_id %in% keep_cells])),
        richness_summary  = summary(kept_stats$n_species),
        occ_summary       = summary(kept_stats$n_occurrences),
        comp_summary      = summary(kept_stats$completeness)
      ),
      thresholds = list(
        min_occurrences  = min_occurrences,
        min_completeness = min_completeness,
        min_richness     = min_richness,
        outlier_cutoff   = outlier_cutoff,
        outlier_source   = outlier_source,
        upper_fence_97   = upper_fence_97
      ),
      filter_applied = TRUE
    )

    # Safety guard: if fewer than 30 cells would survive, skip the animal
    # filter entirely — not enough data to justify it.
    if (length(keep_cells) < 30) {
      warning("  Animal cell filter would retain only ", length(keep_cells),
              " cells (< 30). Skipping filter to preserve data.",
              call. = FALSE)
      .amani_filter_summary$filter_applied <<- FALSE
    } else {
      animal_grid <<- animal_grid[animal_grid$grid_id %in% keep_cells, ]
      animal_grid_data <<- animal_grid_data[animal_grid_data$grid_id %in% keep_cells, ]
    }

    message("  Combined cell filter:")
    message("    Outlier cutoff: ", outlier_source)
    message("    Well-sampled (completeness >= ", round(min_completeness, 2),
            ", occ >= ", min_occurrences, ", spp >= ", min_richness,
            ", spp <= ", round(outlier_cutoff), "): ",
            length(well_sampled_ids), " cells")
    message("    Super-diverse (spp ", round(upper_fence_97), "–",
            round(outlier_cutoff), "): ", length(super_diverse_ids), " cells")
    message("    True outliers dropped (spp > ", round(outlier_cutoff),
            "): ", length(true_outlier_ids), " cells")
    if (length(keep_cells) < 30) {
      message("    *** FILTER SKIPPED — only ", length(keep_cells),
              " cells would remain (minimum 30 required) ***")
      message("    Grid cells: ", n_before, " (unchanged)")
    } else {
      message("    Grid cells: ", n_before, " -> ", length(keep_cells))
    }
    message("    Records: ", nrow(animal_grid))
  }

  if(exists("scaled_results", envir = .GlobalEnv) && !redoeverything){
    NULL}else{
      if(.amani_skipRawEndemism) NULL else source("scripts/utils/rawEndemismCalculator.R")
      
    }
  
  
  # if(alphaRichnessMethod == "extend"){
  # 
  #   # ── Species-list padding based on Chao1 sampling completeness ──────────
  #   #
  #   # For each grid cell, estimate the number of undetected species as
  #   #   n_pad = round(chao1) - Sobs
  #   # then generate that many synthetic names by sampling a base name from
  #   # the global species pool (with replacement) and appending a globally-
  #   # unique incrementing suffix ("_padN").  Padded rows inherit the cell
  #   # geometry and non-species metadata from a representative observed row
  #   # and are marked source_pad = "padded" so downstream code can flag them.
  #   # ─────────────────────────────────────────────────────────────────────────
  #   message("  Species-list padding (alphaRichnessMethod = 'extend'):")
  #   
  #   sp_col <- detect_species_col(sf::st_drop_geometry(animal_grid))
  # 
  #   # Per-cell abundance counts → Sobs, f1, f2, Chao1
  #   spp_abund_ext <- animal_grid |>
  #     sf::st_drop_geometry() |>
  #     dplyr::count(grid_id, .data[[sp_col]], name = "abundance")
  # 
  #   cell_est <- spp_abund_ext |>
  #     dplyr::group_by(grid_id) |>
  #     dplyr::summarise(
  #       Sobs = dplyr::n(),
  #       f1   = sum(abundance == 1L),
  #       f2   = sum(abundance == 2L),
  #       .groups = "drop"
  #     ) |>
  #     dplyr::mutate(
  #       chao1 = Sobs + dplyr::if_else(
  #         f2 > 0,
  #         f1^2 / (2 * f2),
  #         f1 * (f1 - 1) / 2
  #       ),
  #       n_pad = pmax(0L, round(chao1 - Sobs))
  #     )
  # 
  #   total_pads <- sum(cell_est$n_pad)
  # 
  #   if (total_pads > 0L) {
  #     # Global pool of base names to draw from (all observed morphotypes)
  #     sp_pool <- unique(sf::st_drop_geometry(animal_grid)[[sp_col]])
  # 
  #     # Sample one base name per padding slot (with replacement), then
  #     # attach a globally-unique suffix so no two padded names collide
  #     set.seed(SEEDnumber)
  #     base_names <- sample(sp_pool, total_pads, replace = TRUE)
  #     pad_names  <- paste0(base_names, "_pad", seq_len(total_pads))
  # 
  #     # namePadList: named vector mapping padded name -> base species,
  #     # useful for tracing the origin of padded rows later
  #     namePadList <<- stats::setNames(base_names, pad_names)
  # 
  #     # One template row per grid cell (carries geometry + metadata)
  #     cell_template <- animal_grid |>
  #       dplyr::group_by(grid_id) |>
  #       dplyr::slice(1L) |>
  #       dplyr::ungroup()
  # 
  #     # Repeat template rows to match padding slots, then overwrite species col
  #     pad_grid_ids <- rep(cell_est$grid_id, times = cell_est$n_pad)
  #     template_idx <- match(pad_grid_ids, cell_template$grid_id)
  # 
  #     pad_rows <- cell_template[template_idx, ] |>
  #       dplyr::mutate(
  #         !!sp_col   := pad_names,
  #         source_pad  = "padded"
  #       )
  # 
  #     animal_grid <<- dplyr::bind_rows(animal_grid, pad_rows)
  # 
  #     message("    Cells with padding:   ", sum(cell_est$n_pad > 0),
  #             " of ", nrow(cell_est))
  #     message("    Padded rows added:    ", total_pads)
  #     message("    Global pool size:     ", length(sp_pool))
  #   } else {
  #     namePadList <<- character(0)
  #     message("    No padding required — all cells estimated as fully sampled.")
  #   }
  # }
  
  if(fill_animals){
    # Species already range-filled by the pre-extension block above (if it ran).
    # Anything in this set is skipped here so no species is filled twice.
    already_filled <- if (exists(".amani_filled_species", envir = .GlobalEnv)) {
      .amani_filled_species
    } else {
      character(0)
    }

    # animal_grid still has its pre-fill columns (species/Morphotype) here,
    # not `scientificName` — detect the actual species column.
    pre_fill_species_col <- detect_species_col(sf::st_drop_geometry(animal_grid))
    all_species    <- unique(animal_grid[[pre_fill_species_col]])
    species_to_fill <- setdiff(all_species, already_filled)

    if (length(species_to_fill) == 0) {
      message("  Range filling (post-extension): all ", length(all_species),
              " species already filled; skipping.")
    } else {
      message("  Range filling (post-extension): ", length(species_to_fill),
              " of ", length(all_species), " species need filling (",
              length(intersect(all_species, already_filled)),
              " already filled, skipped).")

      # Rows for already-filled species pass through unchanged; only the
      # unfilled (e.g. newly padded) species go through the filling algorithm.
      keep_rows   <- animal_grid[animal_grid[[pre_fill_species_col]] %in% already_filled, ]
      tofill_rows <- animal_grid[animal_grid[[pre_fill_species_col]] %in% species_to_fill, ]

      # Build covariates from the full grid so every cell is represented,
      # regardless of which species are being filled.
      grid_covariates <- build_grid_cell_covariates(animal_grid, completeness = "compute")
      n_species_to_fill <- length(species_to_fill)
      newly_filled <- fill_all_species_ranges(
        tofill_rows, grid_covariates$grid_unique_cells,
        parallel = if(length(unique(tofill_rows$grid_id))>10 && n_species_to_fill>200){ TRUE }else{ FALSE },
        n_cores = max(1L, parallel::detectCores() - 1L),
        verbose = FALSE)

      # fill_all_species_ranges() always names its output column `scientificName`
      # (see fillingRanges.R); restore `Morphotype` so downstream code in
      # 02_animalAnalysis.R that still expects it keeps working.
      newly_filled <- newly_filled |> dplyr::mutate(Morphotype = scientificName)

      # Combine untouched (already-filled) species with the freshly filled ones.
      animal_grid <<- dplyr::bind_rows(keep_rows, newly_filled)
    }

    # Update the tracker so it reflects every species now filled.
    .amani_filled_species <<- unique(animal_grid[[detect_species_col(sf::st_drop_geometry(animal_grid))]])
  }
  
#if fill_animals = TRUE, then run the range filling algorithm.

 
  
  source("scripts/02_animalAnalysis.R")
  
}


