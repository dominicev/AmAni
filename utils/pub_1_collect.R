# utils/pub_1_collect.R — Publication figures, step 1: data collection
#
# Collectors return tidy tibbles tagged with the raw directory name
# (`organism_dir`) and/or the display label (`organism`). Parsing logic mirrors
# summarize_all() so the two stay consistent; once stable, summarize_all()
# should call these instead of its inline code. Uses organism_label(),
# ORGANISM_GROUPS, .latest_file(), .extract_val(), .organism_dirs() and
# load_dataset_summary_figure_inputs() from utils/func_core.R.
#
#   collect_pub_data()               bundles the four collectors below:
#     collect_extinction_summary()   summary_all_models.csv, all models x scenarios
#     collect_organism_summary()     dataset_summary + gamma + sampling_extension
#     collect_dataset_distributions() per-cell five-number summaries (before/after)
#     collect_iucn()                 iucn_threat_classification_data.rds $counts
#
#   Per-cell data from plotdata_cell_filter_report.rds (only organisms with
#   that file are included):
#     collect_retained_cells()       retained cells ($delta) + plant/animal stats
#     collect_cell_distributions()   all animal cells ($animal) + retained flag
#     collect_plant_completeness_raw() Chao1 completeness from RAW plant counts
#     retained_cell_stats()          one row per retained cell, across organisms

# --- Bundle ---------------------------------------------------------------------

# Bundle everything the publication figures need into one named list:
#   $organisms  organism_dir + display label + coarse taxon group
#   $extinction species/population loss summaries (all models x scenarios)
#   $summary    per-organism filtering, gamma and sampling-extension stats
#   $dists      per-cell five-number summaries before/after filtering
#   $iucn       IUCN category counts (all models x scenarios)
#
# `include` filters to specific organism directory names; `group` adds a
# named group from ORGANISM_GROUPS (e.g. group = "vertebrates"). Like the
# rest of the project, the working directory is assumed to be the project
# root, so output_dir defaults to "output".

collect_pub_data <- function(include = NULL, group = NULL,
                             output_dir = "output") {
  if (!is.null(group)) {
    if (!group %in% names(ORGANISM_GROUPS))
      stop("Unknown group '", group, "'. Available: ",
           paste(names(ORGANISM_GROUPS), collapse = ", "), call. = FALSE)
    include <- union(include, ORGANISM_GROUPS[[group]])
  }

  dirs <- .organism_dirs(output_dir)
  if (!is.null(include)) {
    missing <- setdiff(include, dirs)
    if (length(missing) > 0)
      warning("No output directory for: ", paste(missing, collapse = ", "),
              call. = FALSE)
    dirs <- intersect(dirs, include)
  }
  if (length(dirs) == 0)
    stop("No organism output directories found in '", output_dir, "'.",
         call. = FALSE)

  # Coarse taxonomic grouping for colouring (Vertebrates vs everything else)
  is_vert <- dirs %in% ORGANISM_GROUPS$vertebrates

  list(
    organisms = tibble::tibble(
      organism_dir = dirs,
      organism     = organism_label(dirs),
      taxon_group  = ifelse(is_vert, "Vertebrates", "Invertebrates")
    ),
    extinction = collect_extinction_summary(dirs, output_dir),
    summary    = collect_organism_summary(dirs, output_dir),
    dists      = collect_dataset_distributions(dirs, output_dir),
    iucn       = collect_iucn(dirs, output_dir)
  )
}

# --- Per-organism report collectors ------------------------------------------------

# One row per organism x spatial model x deforestation frequency scenario,
# with replicate summaries of species and population losses, read from each
# organism's extinction_model/summary_all_models.csv (latest run only).
collect_extinction_summary <- function(organisms = NULL, output_dir = "output") {
  if (is.null(organisms)) organisms <- .organism_dirs(output_dir)
  purrr::map_dfr(organisms, function(org) {
    f <- .latest_file(file.path(output_dir, org, "extinction_model"),
                      "^summary_all_models\\.csv$")
    if (is.null(f)) return(NULL)
    read.csv(f) |>
      dplyr::mutate(organism_dir = org,
                    organism = organism_label(org), .before = 1)
  })
}

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

# Per-cell five-number summaries (richness / occurrences / completeness)
# before vs after filtering, bound across organisms. One row per organism
# x variable x stage. Reuses scrape_dataset_summary_stats() so this stays
# consistent with the per-organism dataset summary figures.

collect_dataset_distributions <- function(organisms = NULL, output_dir = "output") {
  if (is.null(organisms)) organisms <- .organism_dirs(output_dir)
  purrr::map_dfr(organisms, function(org) {
    inputs <- load_dataset_summary_figure_inputs(
      file.path(output_dir, org, "animal_model"))
    if (is.null(inputs)) return(NULL)
    inputs$dists |>
      dplyr::mutate(organism_dir = org,
                    organism = organism_label(org), .before = 1)
  })
}

# IUCN threat classification counts, one row per organism x spatial model
# x scenario x IUCN category, read from each organism's
# extinction_model/iucn_threat_classification_data.rds (latest run only).
# All models/scenarios are kept; figures subset as needed. Organisms
# without the rds (no IUCN-matched species) are absent.

collect_iucn <- function(organisms = NULL, output_dir = "output") {
  if (is.null(organisms)) organisms <- .organism_dirs(output_dir)
  purrr::map_dfr(organisms, function(org) {
    f <- .latest_file(file.path(output_dir, org, "extinction_model"),
                      "iucn_threat_classification_data\\.rds$")
    if (is.null(f)) return(NULL)
    readRDS(f)$counts |>
      dplyr::select(-dplyr::any_of("x")) |>  # drop rds row-id artifact
      dplyr::mutate(organism_dir = org,
                    organism = organism_label(org), .before = 1)
  })
}

# --- Per-cell collectors (cell-filter reports) -------------------------------------

# Paths to every organism's plotdata_cell_filter_report.rds, named by
# organism directory.
.cell_filter_reports <- function(output_dir = PUB_SETTINGS$output_dir) {
  files <- Sys.glob(file.path(output_dir, "*", "animal_model", "plot_data",
                              "plotdata_cell_filter_report.rds"))
  names(files) <- basename(dirname(dirname(dirname(files))))
  files
}

# Retained cells per organism, from the $delta table: carries per-cell plant
# and animal stats so figures can encode them. `organism` is the directory name.
collect_retained_cells <- function(output_dir = PUB_SETTINGS$output_dir) {
  purrr::imap(.cell_filter_reports(output_dir), \(f, org) {
    readRDS(f)$delta |>
      filter(retained) |>
      mutate(organism = org) |>
      select(organism, grid_id, plant_completeness,
             animal_richness, animal_completeness)
  }) |>
    bind_rows()
}

# Every animal cell per organism, from the $animal table (spans every grid
# cell and carries the `retained` flag), so before/after both come from one file.
collect_cell_distributions <- function(output_dir = PUB_SETTINGS$output_dir) {
  purrr::imap(.cell_filter_reports(output_dir), \(f, org) {
    readRDS(f)$animal |>
      st_drop_geometry() |>
      mutate(organism_dir = org,
             organism = organism_label(org), .before = 1) |>
      select(organism_dir, organism, grid_id, retained,
             animal_completeness, animal_richness)
  }) |>
    bind_rows()
}

# Per-cell Chao1 completeness (and completeness-weighted richness) from the
# RAW filtered plant counts. The plant_completeness stored in the cell-filter
# reports is computed from the FILLED plant data (when fill_plants = TRUE),
# whose inflated record counts leave no singletons, so it is exactly 1 in
# every cell (see the note at the plant_completeness calculation in
# 02_animalAnalysis.R). Figures that encode plant completeness must use this
# raw version instead.
collect_plant_completeness_raw <- function(
    counts_path = file.path("1-FilteringPlantData", "output",
                            "species_counts4_filtered.rds")) {
  readRDS(counts_path) |>
    st_drop_geometry() |>
    group_by(grid_id) |>
    summarise(
      plant_richness = dplyr::n_distinct(scientificName),
      f1 = sum(count == 1),
      f2 = sum(count == 2),
      .groups = "drop"
    ) |>
    mutate(
      plant_chao1 = plant_richness + if_else(
        f2 > 0, f1^2 / (2 * f2), f1 * (f1 - 1) / 2),
      plant_completeness_raw = plant_richness / plant_chao1,
      # completeness x richness: observed richness discounted by how
      # completely it was sampled
      plant_complete_richness = plant_richness * plant_completeness_raw
    ) |>
    select(grid_id, plant_completeness_raw, plant_complete_richness)
}

# One row per grid cell: geometry from any retaining organism plus
# cross-organism summaries (organism count, first plant stats, summed
# animal richness across retaining organisms) and raw plant completeness.
retained_cell_stats <- function(cells,
                                raw_completeness = collect_plant_completeness_raw()) {
  stats <- cells |>
    st_drop_geometry() |>
    group_by(grid_id) |>
    summarise(n_retained        = n_distinct(organism),
              plant_completeness = dplyr::first(plant_completeness),
              animal_richness    = sum(animal_richness),
              animal_completeness = dplyr::first(animal_completeness),
              .groups = "drop")
  cells |>
    group_by(grid_id) |>
    slice(1) |>
    ungroup() |>
    select(grid_id) |>
    left_join(stats, by = "grid_id") |>
    left_join(raw_completeness, by = "grid_id")
}
