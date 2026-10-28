# utils/summarize_organism.R
# Extracted from AmAni.summarize.R. Sourced by AmAni.summarize.R — relies on
# the packages (tidyverse, sf, patchwork, ggrepel) loaded there.

# =============================================================================
# --- 6. Single-organism summary function ------------------------------------
# =============================================================================
#
# summarize_organism("cockroaches") produces:
#   - A multi-panel compilation figure (PDF) combining all key plots
#   - A summary statistics table (CSV) with the same metrics as the
#     cross-organism comparison
#
# Designed to be called interactively after sourcing this file.

summarize_organism <- function(organism, output_dir = NULL) {

  if (!exists("project_root")) stop("Source AmAni.manager.R first.", call. = FALSE)

  if (is.null(output_dir)) output_dir <- file.path(project_root, "output")

  org_am  <- file.path(output_dir, organism, "animal_model")
  org_ext <- file.path(output_dir, organism, "extinction_model")
  org_out <- file.path(output_dir, organism, "organism_summary")
  dir.create(org_out, showWarnings = FALSE, recursive = TRUE)

  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")

  load_pd <- function(filename) {
    path <- file.path(org_am, "plot_data", filename)
    if (file.exists(path)) readRDS(path) else NULL
  }

  if (!dir.exists(file.path(org_am, "plot_data"))) {
    stop("No plot_data/ directory found for '", organism,
         "'. Run the animal model pipeline first.", call. = FALSE)
  }

  message("Summarizing single organism: ", organism)

  extract_val_local <- function(lines, pattern) {
    line <- grep(pattern, lines, value = TRUE, ignore.case = TRUE)
    if (length(line) == 0) return(NA_real_)
    raw <- sub(".*:\\s*", "", line[1])
    m <- regmatches(raw, regexpr("-?[0-9]+\\.?[0-9]*", raw))
    if (length(m) == 0) return(NA_real_)
    as.numeric(m)
  }

  extract_val_char <- function(lines, pattern) {
    line <- grep(pattern, lines, value = TRUE, ignore.case = TRUE)
    if (length(line) == 0) return(NA_character_)
    trimws(sub(".*:\\s*", "", line[1]))
  }

  ds_file <- list.files(org_am, pattern = "^dataset_summary_.*\\.txt$", full.names = TRUE)
  if (length(ds_file) > 0) {
    ds_file <- ds_file[which.max(file.mtime(ds_file))]
    ds <- readLines(ds_file)
    pre_cells   <- extract_val_local(ds, "^\\s*Grid cells:")
    pre_records <- extract_val_local(ds, "^\\s*Total records:")
    pre_species <- extract_val_local(ds, "^\\s*Total species:")
    after_idx <- grep("After filtering", ds)
    if (length(after_idx) > 0) {
      ds_after <- ds[after_idx[1]:length(ds)]
      post_cells   <- extract_val_local(ds_after, "^\\s*Grid cells:")
      post_records <- extract_val_local(ds_after, "^\\s*Total records:")
      post_species <- extract_val_local(ds_after, "^\\s*Total species:")
    } else {
      ds_after <- character(0)
      post_cells <- post_records <- post_species <- NA_real_
    }
    filter_applied <- any(grepl("Filter applied:\\s*YES", ds))
    min_occ  <- extract_val_local(ds, "min_occurrences:")
    min_comp <- extract_val_local(ds, "min_completeness:")
    min_rich <- extract_val_local(ds, "min_richness:")
    outlier_cutoff <- extract_val_char(ds, "outlier cutoff:")
    track_b_entry  <- extract_val_char(ds, "Track B entry:")

    gc_line <- grep("^\\s*Grid cells:", ds_after, value = TRUE)
    if (length(gc_line) > 0) {
      well_sampled_m   <- regmatches(gc_line[1], regexpr("([0-9]+)\\s+well-sampled", gc_line[1]))
      super_diverse_m  <- regmatches(gc_line[1], regexpr("([0-9]+)\\s+super-diverse", gc_line[1]))
      outliers_drop_m  <- regmatches(gc_line[1], regexpr("([0-9]+)\\s+outliers? dropped", gc_line[1]))
      well_sampled   <- if (length(well_sampled_m))  as.integer(sub("\\s.*", "", well_sampled_m))  else NA_integer_
      super_diverse  <- if (length(super_diverse_m)) as.integer(sub("\\s.*", "", super_diverse_m)) else NA_integer_
      outliers_dropped <- if (length(outliers_drop_m)) as.integer(sub("\\s.*", "", outliers_drop_m)) else NA_integer_
    } else {
      well_sampled <- super_diverse <- outliers_dropped <- NA_integer_
    }
    spp_loss_pct <- extract_val_local(ds, "Species lost:")
  } else {
    pre_cells <- pre_records <- pre_species <- NA_real_
    post_cells <- post_records <- post_species <- NA_real_
    filter_applied <- FALSE
    min_occ <- min_comp <- min_rich <- NA_real_
    outlier_cutoff <- track_b_entry <- NA_character_
    well_sampled <- super_diverse <- outliers_dropped <- NA_integer_
    spp_loss_pct <- NA_real_
  }

  gm_file <- list.files(org_am, pattern = "^gamma_richnessfromSPRICH.*\\.txt$", full.names = TRUE)
  if (length(gm_file) > 0) {
    gm_file <- gm_file[which.max(file.mtime(gm_file))]
    gm <- readLines(gm_file)
    extent_scaling_factor <- extract_val_local(gm, "Extent scaling factor:|Scaling factor:")

    fit_params <- extract_val_char(gm, "Fit parameters:")
    gm_delta_meanlog     <- extract_val_local(gm, "^delta_meanlog:")
    gm_delta_meanlog_str <- extract_val_char(gm, "^delta_meanlog:")
    gm_delta_active      <- grepl("active", gm_delta_meanlog_str, ignore.case = TRUE)
    gm_ratio_sdlog       <- extract_val_local(gm, "^ratio_sdlog:")
    gm_ratio_sdlog_str   <- extract_val_char(gm, "^ratio_sdlog:")
    gm_ratio_active      <- grepl("active", gm_ratio_sdlog_str, ignore.case = TRUE)
    fill_plants_gm <- grepl("Plant range filling.*was used", paste(gm, collapse = " "))

    # New format: total model (surrogate + extension) as primary
    gamma_best  <- extract_val_local(gm, "^Best estimate from total model")
    gamma_upper <- extract_val_local(gm, "^Upper bound:\\s")
    gamma_lower <- extract_val_local(gm, "^Lower bound:\\s")

    # Fall back to old format
    if (is.na(gamma_best)) {
      gamma_best  <- extract_val_local(gm, "Best estimate:")
      gamma_upper <- extract_val_local(gm, "Upper bound:")
      gamma_lower <- extract_val_local(gm, "Lower bound:")
    }

    # Surrogate-only estimates (new format)
    surrogate_best  <- extract_val_local(gm, "^Best estimate from plant surrogate")
    surrogate_upper <- extract_val_local(gm, "^Upper bound \\(narrow")
    surrogate_lower <- extract_val_local(gm, "^Lower bound \\(narrow")
    richness_shortfall <- extract_val_local(gm, "^Expected proportional richness shortfall")
  } else {
    extent_scaling_factor <- gamma_best <- gamma_upper <- gamma_lower <- NA_real_
    surrogate_best <- surrogate_upper <- surrogate_lower <- richness_shortfall <- NA_real_
    fit_params <- NA_character_
    gm_delta_meanlog <- gm_ratio_sdlog <- NA_real_
    gm_delta_active <- gm_ratio_active <- FALSE
    fill_plants_gm <- NA
  }

  # --- Sampling extension report ---
  se_file <- list.files(org_am, pattern = "^sampling_extension_report_.*\\.txt$", full.names = TRUE)
  if (length(se_file) > 0) {
    se_file <- se_file[which.max(file.mtime(se_file))]
    se <- readLines(se_file)

    se_fill_plants       <- extract_val_char(se, "fill_Plants:")
    se_comp_spread       <- extract_val_local(se, "comp_spread:")
    se_comp_mean         <- extract_val_local(se, "comp_mean:")
    se_comp_sd           <- extract_val_local(se, "comp_sd:")
    se_prop_species      <- extract_val_local(se, "prop_species \\(proxy\\):")
    se_chao2_s_obs       <- extract_val_local(se, "S_obs =")
    se_chao2_s_chao2     <- extract_val_local(se, "S_Chao2 =")
    se_n_grid_cells      <- extract_val_local(se, "n_grid_cells \\(raw\\):")
    se_mean_range        <- extract_val_local(se, "sample_mean_range_size:")
    se_delta_meanlog     <- extract_val_local(se, "^\\s*delta_meanlog:")
    se_ratio_sdlog       <- extract_val_local(se, "^\\s*ratio_sdlog:")
    se_partial_gamma     <- extract_val_local(se, "partialGamma:")
    se_abs_prop_error    <- extract_val_local(se, "abs_prop_error:")
    se_correction_factor <- extract_val_local(se, "correction_factor:")
    se_full_gamma_best   <- extract_val_local(se, "fullGamma \\(best\\):")
    se_full_gamma_lower  <- extract_val_local(se, "fullGamma \\(lower")
    se_full_gamma_upper  <- extract_val_local(se, "fullGamma \\(upper")
  } else {
    se_fill_plants <- NA_character_
    se_comp_spread <- se_comp_mean <- se_comp_sd <- se_prop_species <- NA_real_
    se_chao2_s_obs <- se_chao2_s_chao2 <- NA_real_
    se_n_grid_cells <- se_mean_range <- NA_real_
    se_delta_meanlog <- se_ratio_sdlog <- NA_real_
    se_partial_gamma <- se_abs_prop_error <- se_correction_factor <- NA_real_
    se_full_gamma_best <- se_full_gamma_lower <- se_full_gamma_upper <- NA_real_
  }

  cs_file <- list.files(org_am, pattern = "^correlation_summary_.*\\.txt$", full.names = TRUE)
  if (length(cs_file) > 0) {
    cs_file <- cs_file[which.max(file.mtime(cs_file))]
    cs <- readLines(cs_file)
    animal_cells <- extract_val_local(cs, "Animal cells used:")
    plant_comp_filter <- extract_val_char(cs, "Plant completeness filter:")

    rich_idx <- grep("Richness correlation", cs)
    if (length(rich_idx) > 0) {
      cs_rich <- cs[rich_idx[1]:min(rich_idx[1] + 10, length(cs))]
      rich_cells  <- extract_val_local(cs_rich, "Cells:")
      rich_r2     <- extract_val_local(cs_rich, "^\\s*R.:")
      rich_adj_r2 <- extract_val_local(cs_rich, "Adj\\. R.:")
      rich_p      <- extract_val_local(cs_rich, "p-value:")
      rich_slope  <- extract_val_local(cs_rich, "Slope:")
    } else {
      rich_cells <- rich_r2 <- rich_adj_r2 <- rich_p <- rich_slope <- NA_real_
    }
    end_idx <- grep("Endemism correlation", cs)
    if (length(end_idx) > 0) {
      cs_end <- cs[end_idx[1]:min(end_idx[1] + 10, length(cs))]
      end_r2     <- extract_val_local(cs_end, "^\\s*R.:")
      end_adj_r2 <- extract_val_local(cs_end, "Adj\\. R.:")
      end_p      <- extract_val_local(cs_end, "p-value:")
      end_slope  <- extract_val_local(cs_end, "Slope:")
    } else {
      end_r2 <- end_adj_r2 <- end_p <- end_slope <- NA_real_
    }
  } else {
    animal_cells <- rich_cells <- rich_r2 <- rich_adj_r2 <- rich_p <- rich_slope <- NA_real_
    end_r2 <- end_adj_r2 <- end_p <- end_slope <- NA_real_
    plant_comp_filter <- NA_character_
  }

  observed_spp <- if (!is.na(post_species) && filter_applied) post_species else pre_species

  richness_map   <- load_pd("plotdata_animal_richness_map.rds")
  endemism_map   <- load_pd("plotdata_animal_endemism_map.rds")
  endem_corr     <- load_pd("plotdata_endemism_correlation.rds")
  rich_corr      <- load_pd("plotdata_richness_correlation.rds")
  endemism_ratio_data <- load_pd("plotdata_endemism_ratios.rds")
  if (is.null(endemism_ratio_data)) endemism_ratio_data <- load_pd("plotdata_scaling_factors.rds")
  # endem_sim not loaded -- endemism-ratio simulation kept in individual output only
  extent_model   <- load_pd("plotdata_plant_extent_model.rds")
  animal_extent  <- load_pd("plotdata_animal_extent_distribution.rds")
  scaled_map     <- load_pd("plotdata_scaled_endemism_map.rds")
  sim_rich_map   <- load_pd("plotdata_simulated_richness_map.rds")

  panels <- list()

  # Panels are only added when their underlying plot_data exists and is
  # non-empty; empty plots are skipped so they don't appear in the compilation.
  has_rows <- function(d) !is.null(d) && nrow(d) > 0

  # --- Dataset before/after cell-filtering summary (first panel) ---
  ds_summary_inputs <- load_dataset_summary_figure_inputs(org_am)
  if (!is.null(ds_summary_inputs)) {
    panels$dataset_summary <- build_dataset_summary_figure(
      ds_summary_inputs$thresholds, ds_summary_inputs$counts, ds_summary_inputs$dists,
      title = "Dataset before vs after cell filtering",
      subtitle = NULL)
  }

  # --- Richness & endemism in the area of overlap ---
  # Panels are only added when their underlying data exist; empty plots are
  # skipped entirely rather than filled with a placeholder.
  if (has_rows(richness_map)) {
    panels$richness_map <- ggplot(richness_map) +
      geom_sf(aes(fill = animal_richness), color = NA) +
      scale_fill_viridis_c(option = "magma", direction = -1, name = "Species\nrichness") +
      theme_minimal() + labs(title = "Animal Richness")
  }
  if (has_rows(endemism_map)) {
    panels$endemism_map <- ggplot(endemism_map) +
      geom_sf(aes(fill = animal_endemism_count), color = NA) +
      scale_fill_viridis_c(option = "plasma", direction = -1, name = "Endemic\nspecies") +
      theme_minimal() + labs(title = "Animal Endemism")
  }

  # --- Plant vs animal richness ---
  if (has_rows(rich_corr)) {
    r2_val <- tryCatch(
      summary(lm(morphotype_richness ~ plant_richness, data = rich_corr))$r.squared,
      error = function(e) NA_real_)
    panels$rich_corr <- ggplot(rich_corr, aes(x = plant_richness, y = morphotype_richness)) +
      geom_point(alpha = 0.6, color = "darkgreen") +
      geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.8) +
      annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.5,
               label = sprintf("R2 = %.3f  (n = %d)", r2_val, nrow(rich_corr)), size = 3) +
      theme_minimal() +
      labs(title = "Plant vs Animal Richness", x = "Plant richness", y = "Animal richness")
  }

  # --- Plant vs animal endemism ---
  endem_corr_pos <- if (has_rows(endem_corr)) {
    endem_corr |> filter(raw_endemism > 0, animal_endemism_count > 0)
  } else {
    NULL
  }
  if (has_rows(endem_corr_pos)) {
    r2_val <- tryCatch(
      summary(lm(log(animal_endemism_count) ~ log(raw_endemism), data = endem_corr_pos))$r.squared,
      error = function(e) NA_real_)
    panels$endem_corr <- ggplot(endem_corr_pos, aes(x = raw_endemism, y = animal_endemism_count)) +
      geom_point(alpha = 0.6, color = "darkgreen") +
      geom_smooth(method = "lm", se = TRUE, color = "red", linewidth = 0.8) +
      scale_x_log10() + scale_y_log10() +
      annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.5,
               label = sprintf("R2 = %.3f  (n = %d)", r2_val, nrow(endem_corr_pos)), size = 3) +
      theme_minimal() +
      labs(title = "Plant vs Animal Endemism (log-log)",
           x = "Plant raw endemism (log)", y = "Animal endemism count (log)")
  }
  if (!is.null(endemism_ratio_data)) {
    if ("scaling_factor" %in% names(endemism_ratio_data) && !"endemism_ratio" %in% names(endemism_ratio_data)) {
      endemism_ratio_data <- dplyr::rename(endemism_ratio_data, endemism_ratio = scaling_factor)
    }
    panels$endemism_ratios <- ggplot(endemism_ratio_data, aes(x = endemism_ratio)) +
      geom_histogram(bins = 30, fill = "darkgreen", color = "black") +
      theme_minimal() +
      labs(title = "Endemism Ratios (Animal / Plant, < 2)", x = "Endemism ratio", y = "Frequency")
  }

  # Panel F (Simulated endemic species): SKIPPED -- endemism-ratio simulation
  # kept in individual output only. See endemism_ratio_simulation_results_*.txt.

  # --- Animal spatial spread distribution (spat.Spread) ---
  # Placed before the plant spatial extent model. Shows how many morphotypes
  # occupy each number of grid cells (the animal-side spatial spread).
  if (has_rows(animal_extent)) {
    panels$animal_extent <- ggplot(animal_extent, aes(x = num_grid_cells, y = n)) +
      geom_col(fill = "lightblue") +
      theme_minimal() +
      labs(title = "Animal Spatial Spread",
           x = "Number of grid cells", y = "Number of species")
  }

  # --- Plant spatial extent model ---
  if (!is.null(extent_model)) {
    plant_df <- extent_model$plant_df
    ci_df    <- extent_model$ci_df
    ci_df$y_pred      <- extent_model$y_pred
    ci_df$y_predUpper <- extent_model$y_predUpper
    panels$extent_model <- ggplot() +
      geom_ribbon(data = ci_df, aes(x = num_grid_cells, ymin = lower, ymax = upper),
                  fill = "darkgreen", alpha = 0.2) +
      geom_line(data = ci_df, aes(x = num_grid_cells, y = median), color = "darkgreen", linewidth = 0.8) +
      geom_line(data = ci_df, aes(x = num_grid_cells, y = y_pred), color = "steelblue", linetype = "dashed", linewidth = 0.6) +
      geom_line(data = ci_df, aes(x = num_grid_cells, y = y_predUpper), color = "steelblue", linetype = "dashed", linewidth = 0.6) +
      geom_point(data = plant_df, aes(x = num_grid_cells, y = n), color = "darkgreen", size = 1.5) +
      scale_x_log10() + scale_y_log10() + theme_minimal() +
      labs(title = "Estimated spatial extent model\n(surrogate)", x = "Grid cells (log)", y = "Species (log)")
  }

  # --- Scaled plant endemism ---
  if (has_rows(scaled_map) && any(!is.na(scaled_map$scaled_plant_endemism))) {
    panels$scaled_map <- ggplot(scaled_map) +
      geom_sf(aes(fill = scaled_plant_endemism), color = NA) +
      scale_fill_viridis_c(option = "plasma", na.value = "white", name = "Scaled\nendemism") +
      theme_minimal() + labs(title = "Scaled Plant Endemism")
  }

  # --- Simulated species richness ---
  if (has_rows(sim_rich_map)) {
    panels$sim_richness_map <- ggplot(sim_rich_map) +
      geom_sf(aes(fill = sim_richness), color = NA) +
      scale_fill_viridis_c(option = "magma", direction = -1, na.value = "white", name = "Simulated\nrichness") +
      theme_minimal() + labs(title = "Simulated animal species")
  }

  # --- Spatial extent overlap (single organism) ---
  spatial_extent_plot <- plot_spatial_extent_overlap(
    include = organism, output_dir = output_dir, save_pdf = FALSE)
  if (!is.null(spatial_extent_plot)) {
    panels$spatial_extent <- spatial_extent_plot
  }

  # --- Observed vs estimated gamma richness ---
  # Observed baseline uses the pre-filter total species; the analyzed portion
  # is noted as text.
  observed_baseline <- if (!is.na(pre_species)) pre_species else observed_spp
  if (!is.na(observed_baseline) && !is.na(gamma_best)) {
    expansion_ratioHere <- round((gamma_best/observed_baseline) , digits = 1) #NOTE this is not actually the percentage undiscovered. it's the proportion of species predictyed based on the osberved data
    # Three columns: observed, surrogate-only gamma, then the full/total-model
    # gamma. The surrogate estimate sits between observed and the total, and
    # both estimate columns carry their lower/upper bounds as error bars.
    gamma_levels <- c("Observed\n(total)", "Surrogate\ngamma", "Estimated\ngamma")
    gamma_df <- tibble(
      metric  = factor(gamma_levels, levels = gamma_levels),
      species = c(observed_baseline, surrogate_best, gamma_best),
      lower   = c(NA_real_, surrogate_lower, gamma_lower),
      upper   = c(NA_real_, surrogate_upper, gamma_upper)) |>
      dplyr::filter(!is.na(species))
    panels$gamma_richness <- ggplot(gamma_df, aes(x = metric, y = species, fill = metric)) +
      geom_col(width = 0.6) +
      geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, na.rm = TRUE) +
      scale_fill_manual(values = c("Observed\n(total)" = "grey50",
                                   "Surrogate\nmodel\ngamma"  = "#74a892",
                                   "Total\nmodel\nestimated\ngamma"  = "darkgreen"),
                        guide = "none") +
      annotate("text", x = "Estimated\ngamma",
               y = max(gamma_upper, gamma_best, na.rm = TRUE) * 0.5,
               label = paste0(expansion_ratioHere, "X\nrichness\nexpansion"),
               size = 3.2, fontface = "bold") +
      annotate("text", x = "Observed\n(total)", y = observed_baseline,
               label = paste0("n analyzed = ", scales::comma(observed_spp)),
               vjust = -0.6, size = 2.8) +
      theme_minimal() + labs(title = "Total Amazonian richness", x = NULL, y = "Species count")
  }

  # --- Extinct species estimation ---
  ext_file <- file.path(org_ext, "summary_all_models.csv")
  extinct_clustered <- NA_real_
  if (file.exists(ext_file)) {
    ext_data <- read.csv(ext_file)
    # Extinct count for the composition column: FreqBoundary / Clustered median.
    ext_row <- ext_data[ext_data$scenario == "FreqBoundary" &
                          ext_data$model == "Clustered", ]
    if (nrow(ext_row) > 0) extinct_clustered <- ext_row$median_spp_lost[1]
  }
  # NOTE: the standalone "Species Lost Under Deforestation" panel
  # (panels$extinction) built directly from ext_data remains disabled below;
  # only the extinct_clustered scalar used by the composition panel is
  # restored here.
  # if (file.exists(ext_file)) {
  #   panels$extinction <- ggplot(ext_data, aes(x = scenario, y = median_spp_lost, fill = model)) +
  #     geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  #     geom_errorbar(aes(ymin = q25_spp_lost, ymax = q75_spp_lost),
  #                   position = position_dodge(width = 0.8), width = 0.25) +
  #     theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  #     labs(title = "Species Lost Under Deforestation", x = NULL, y = "Species lost (median, IQR)", fill = "Spatial model")
  # }

  # --- IUCN threat classification panel ---
  # Reuses the same functions that generate each organism's standalone
  # <organism>_iucn_threat_classification.pdf, so the panel is identical to
  # that plot. Output is written to a temp dir to avoid re-touching the
  # organism's extinction_model/ files.
  thresh_path <- file.path(org_ext, "summary_range_reduction_thresholds.csv")
  if (file.exists(thresh_path)) {
    source(file.path(project_root, "4-ExtModeller", "scripts", "utils",
                     "iucn_threat_classification.R"))
    iucn_classification <- compute_iucn_category_counts(org_ext)
    panels$iucn_threat <- plot_iucn_threat(
      iucn_classification,
      output_dir = tempdir(),
      filename_prefix = paste0("iucn_threat_", organism),
      organism_label = organism_label(organism))
  }

  # --- Composition + IUCN threat stacked columns ---
  # Column 1 ("Composition"): proportion of species observed, predicted-but-
  # unsampled, and projected extinct. Predicted uses the midpoint of the gamma
  # error bars (NOT the total best estimate). Column 2 ("IUCN threat"): the
  # IUCN Red List category breakdown (Clustered / FreqBoundary), colored with
  # the IUCN colour scheme.
  iucn_file <- list.files(org_ext,
    pattern = "iucn_threat_classification_data\\.rds$", full.names = TRUE)

  comp_df <- NULL
  if (!is.na(observed_baseline) && !is.na(gamma_best)) {
    gamma_mid <- if (!is.na(gamma_lower) && !is.na(gamma_upper)) {
      (gamma_lower + gamma_upper) / 2
    } else {
      gamma_best
    }
    comp_df <- tibble::tibble(
      column   = "Composition",
      category = c("Observed", "Predicted (unsampled)", "Extinct"),
      count    = c(observed_baseline,
                   pmax(gamma_mid - observed_baseline, 0),
                   dplyr::coalesce(extinct_clustered, 0)))
  }

  iucn_df <- NULL
  if (length(iucn_file) > 0) {
    iucn_data <- readRDS(iucn_file[which.max(file.mtime(iucn_file))])
    iucn_counts <- iucn_data$counts |>
      dplyr::filter(model == "Clustered", scenario == "FreqBoundary")
    if (nrow(iucn_counts) > 0) {
      iucn_df <- tibble::tibble(
        column   = "IUCN threat",
        category = as.character(iucn_counts$iucn_category),
        count    = iucn_counts$count)
    }
  }

  stacked_df <- dplyr::bind_rows(comp_df, iucn_df)
  if (!is.null(stacked_df) && nrow(stacked_df) > 0) {
    cat_levels <- c("Observed", "Predicted (unsampled)", "Extinct",
                    "NT/LC", "VU", "EN", "CR/EX")
    stacked_df <- stacked_df |>
      dplyr::mutate(
        column   = factor(column, levels = c("Composition", "IUCN threat")),
        category = factor(category, levels = cat_levels))
    stacked_cols <- c(
      "Observed"              = "grey50",
      "Predicted (unsampled)" = "#2d6a4f",
      "Extinct"               = "#c1121f",
      "NT/LC"                 = "#60C659",
      "VU"                    = "#F9E814",
      "EN"                    = "#FC7F3F",
      "CR/EX"                 = "#D81E05")
    panels$composition_iucn <- ggplot(stacked_df,
        aes(x = column, y = count, fill = category)) +
      geom_col(position = "fill", width = 0.7) +
      scale_y_continuous(labels = scales::percent) +
      scale_fill_manual(values = stacked_cols, drop = FALSE, name = NULL) +
      theme_minimal() +
      labs(title = "Species composition & IUCN threat\nbased on past deforestation\nonly",
           x = NULL, y = "Percent of total")
  }

  if (length(panels) == 0) {
    message("  No plot data found for '", organism, "'. Returning stats only.")
  }

  # The dataset before/after summary (if present) is pulled out and placed
  # on its own full-width row at the top of the PDF, since it is itself a
  # multi-panel compilation and would be cramped inside the uniform grid.
  # The remaining panels keep the fixed 3-column layout so the PDF
  # dimensions below (which assume 3 columns) match what patchwork renders.
  has_ds_panel <- "dataset_summary" %in% names(panels)
  ds_panel <- if (has_ds_panel) panels$dataset_summary else NULL
  rest_panels <- panels[setdiff(names(panels), "dataset_summary")]

  n_rest <- length(rest_panels)
  ncol_rest <- min(max(n_rest, 1), 3)
  rest_grid <- if (n_rest > 0) wrap_plots(rest_panels, ncol = ncol_rest) else NULL

  if (!is.null(ds_panel) && !is.null(rest_grid)) {
    # wrap_elements() treats each nested patchwork as a single graphical
    # element -- combining two already-composed patchworks directly with
    # "/" flattens their sub-panels and breaks the nrow/ncol accounting.
    combined <- (wrap_elements(full = ds_panel) / wrap_elements(full = rest_grid)) +
      plot_layout(heights = c(1.4, ceiling(n_rest / ncol_rest)))
  } else if (!is.null(ds_panel)) {
    combined <- ds_panel
  } else {
    combined <- rest_grid
  }
  combined <- combined +
    plot_annotation(title = organism_label(organism),
                    theme = theme(plot.title = element_text(size = 16, face = "bold")))

  fig_file <- file.path(org_out, paste0("compilation_", organism, "_", timestamp, ".pdf"))
  fig_h <- ceiling(n_rest / 3) * 5 + (if (!is.null(ds_panel)) 7 else 0)
  fig_w <- max(ncol_rest, 1) * 5
  pdf(fig_file, width = fig_w, height = fig_h)
  print(combined)
  dev.off()
  message("  Wrote: ", fig_file)

  # Two-page cell-filtering report (filter stats + maps; completeness correlation)
  cfr_file <- write_cell_filter_report_pdf(organism, org_am, org_out, timestamp)

  stats_tbl <- tibble(
    organism = organism,
    pre_cells = pre_cells, pre_records = pre_records, pre_species = pre_species,
    filter_applied = filter_applied,
    min_occurrences = min_occ, min_completeness = min_comp, min_richness = min_rich,
    outlier_cutoff = outlier_cutoff, track_b_entry = track_b_entry,
    post_cells = post_cells, post_records = post_records, post_species = post_species,
    well_sampled = well_sampled, super_diverse = super_diverse,
    outliers_dropped = outliers_dropped, spp_loss_pct = spp_loss_pct,
    model_cells = animal_cells, extent_scaling_factor = extent_scaling_factor,
    fit_params = fit_params, fill_plants = fill_plants_gm,
    delta_meanlog = gm_delta_meanlog, delta_meanlog_active = gm_delta_active,
    ratio_sdlog = gm_ratio_sdlog, ratio_sdlog_active = gm_ratio_active,
    observed_spp = observed_spp,
    gamma_best = gamma_best, gamma_lower = gamma_lower, gamma_upper = gamma_upper,
    surrogate_best = surrogate_best, surrogate_lower = surrogate_lower,
    surrogate_upper = surrogate_upper, richness_shortfall = richness_shortfall,
    gamma_difference = gamma_best - observed_spp,
    gamma_ratio = round(gamma_best / observed_spp, 2),
    se_fill_plants = se_fill_plants,
    se_comp_spread = se_comp_spread, se_comp_mean = se_comp_mean, se_comp_sd = se_comp_sd,
    se_prop_species = se_prop_species,
    se_chao2_s_obs = se_chao2_s_obs, se_chao2_s_chao2 = se_chao2_s_chao2,
    se_n_grid_cells = se_n_grid_cells, se_mean_range_size = se_mean_range,
    se_delta_meanlog = se_delta_meanlog, se_ratio_sdlog = se_ratio_sdlog,
    se_partial_gamma = se_partial_gamma, se_abs_prop_error = se_abs_prop_error,
    se_correction_factor = se_correction_factor,
    se_full_gamma_best = se_full_gamma_best,
    se_full_gamma_lower = se_full_gamma_lower,
    se_full_gamma_upper = se_full_gamma_upper,
    plant_comp_filter = plant_comp_filter,
    rich_corr_cells = rich_cells, rich_r2 = rich_r2, rich_adj_r2 = rich_adj_r2,
    rich_p = rich_p, rich_slope = rich_slope,
    end_r2 = end_r2, end_adj_r2 = end_adj_r2, end_p = end_p, end_slope = end_slope)

  if (file.exists(ext_file)) {
    ext_tbl <- read.csv(ext_file) |> mutate(organism = organism, .before = 1)
    csv_file <- file.path(org_out, paste0("summary_", organism, "_", timestamp, ".csv"))
    write.csv(stats_tbl, csv_file, row.names = FALSE)
    message("  Wrote: ", csv_file)
    ext_csv <- file.path(org_out, paste0("extinction_", organism, "_", timestamp, ".csv"))
    write.csv(ext_tbl, ext_csv, row.names = FALSE)
    message("  Wrote: ", ext_csv)
  } else {
    csv_file <- file.path(org_out, paste0("summary_", organism, "_", timestamp, ".csv"))
    write.csv(stats_tbl, csv_file, row.names = FALSE)
    message("  Wrote: ", csv_file)
  }

  cat("\n=== ", organism_label(organism), " -- Dataset & Model Summary ===\n\n")
  stats_tbl |>
    dplyr::select(organism, pre_species, post_species, model_cells,
                  extent_scaling_factor, fit_params, observed_spp,
                  gamma_best, gamma_difference, gamma_ratio,
                  se_correction_factor, richness_shortfall,
                  rich_r2, rich_p) |>
    as.data.frame() |>
    print(row.names = FALSE)

  invisible(list(stats = stats_tbl, figure_file = fig_file))
}
