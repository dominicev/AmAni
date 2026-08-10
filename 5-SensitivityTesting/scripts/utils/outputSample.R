#' Process and Export Sampling Scenario Results
#'
#' Combines occurrence data across scenarios, computes range-size and community
#' summaries, extracts sampling attributes, and exports everything to CSV/RData.
#'
#' Requires: data.table, dplyr, sf, purrr, ggplot2, e1071
#'
#' @param sampling_results Named list of sampling results (each element has
#'   `$occurenceData` and sampling attributes).
#' @param scenario_levels Character vector of scenario names in order.
#' @param sampling_levels Numeric vector of sampling fractions (same length as
#'   `scenario_levels`).
#' @param true_occ_data Full true occurrence data frame (passed to
#'   `rebuild_community2()`).
#' @param output_folder Path to output directory. Created if it doesn't exist.
#' @param output_prefix String used as the base name for all output files.
#' @param obs_per_speciesInCell Numeric scalar recorded in the community summary
#'   (default 1).
#' @param plot Logical; whether to print diagnostic ggplots (default TRUE).
#' @return Invisibly returns a list with `scenario_summary_full`, `all_samples`,
#'   `all_comm`, `comm_list`, and `sample_species_ranges`.
output_sample_results <- function(sampling_results,
                                  scenario_levels,
                                  sampling_levels,
                                  true_occ_data,
                                  output_folder,
                                  output_prefix = "lowtoHighSpeciesSampling",
                                  obs_per_speciesInCell = 1,
                                  plot = TRUE) {

  # Ensure output folder exists

  if (!dir.exists(output_folder)) {
    dir.create(output_folder, recursive = TRUE)
    message("Created output folder: ", output_folder)
  }


  # --- 1. Combine occurrence data across scenarios ---
  message("Combining occurrence data across scenarios...")

  all_samples <- data.table::rbindlist(
    lapply(sampling_results, `[[`, "occurenceData"),
    idcol = "scenario"
  )

  all_samples <- all_samples %>%
    dplyr::mutate(scenario = factor(
      scenario,
      levels = scenario_levels,
      ordered = TRUE
    ))


  # --- 2. Compute species range summaries ---
  message("Computing species range summaries...")

  sample_species_ranges <- all_samples %>%
    sf::st_drop_geometry() %>%
    dplyr::distinct(species, grid_id, scenario) %>%
    dplyr::count(species, scenario, name = "n_cells")

  scenario_summary1 <- sample_species_ranges %>%
    dplyr::group_by(scenario) %>%
    dplyr::summarise(
      mean_range_size   = mean(n_cells),
      median_range_size = median(n_cells),
      sd_range_size     = sd(n_cells),
      skew_range_size   = e1071::skewness(n_cells),
      .groups = "drop"
    )

  message("Range size summary:")
  print(scenario_summary1)

  # Diagnostic plots
  if (plot) {
    message("Generating diagnostic plots...")

    p1 <- ggplot2::ggplot(sample_species_ranges,
                          ggplot2::aes(x = log10(n_cells), fill = scenario)) +
      ggplot2::geom_density(alpha = 0.35) +
      ggplot2::theme_bw() +
      ggplot2::labs(x = "log10(range size in grid cells)", y = "Density")
    print(p1)

    p2 <- ggplot2::ggplot(sample_species_ranges,
                          ggplot2::aes(x = scenario, y = n_cells)) +
      ggplot2::geom_boxplot() +
      ggplot2::scale_y_log10() +
      ggplot2::theme_bw() +
      ggplot2::labs(y = "Range size (grid cells, log scale)")
    print(p2)
  }


  # --- 3. Rebuild community matrices ---
  message("Rebuilding community matrices...")

  comm_list <- purrr::imap(sampling_results, function(result, scenario_name) {
    rebuild_community2(
      data.table::setDT(true_occ_data),
      data.table::setDT(result$occurenceData),
      scenario_name
    )
  })

  all_comm <- data.table::rbindlist(comm_list, idcol = "temp_scenario_id")
  all_comm[, temp_scenario_id := NULL]
  all_comm[, scenario := factor(scenario, levels = scenario_levels, ordered = TRUE)]


  # --- 4. Compute community summaries ---
  message("Computing community summaries...")

  scenario_summary2 <- all_comm %>%
    dplyr::group_by(scenario) %>%
    dplyr::summarise(
      n_species             = dplyr::n_distinct(species),
      n_grid_cells          = dplyr::n_distinct(grid_id),
      total_observations    = dplyr::n(),
      mean_obs_per_species  = dplyr::n() / dplyr::n_distinct(species),
      obs_per_speciesInCell = obs_per_speciesInCell,
      mean_obs_per_cell     = dplyr::n() / dplyr::n_distinct(grid_id),
      mean_rangeSize        = dplyr::n_distinct(species) / dplyr::n_distinct(grid_id),
      .groups = "drop"
    )


  # --- 5. Extract sampling attributes ---
  message("Extracting sampling attributes...")

  sampling_attrs <- purrr::imap(sampling_results, function(result, scenario_name) {
    attrs <- attributes(result)
    get_attr <- function(x) if (is.null(x) || length(x) == 0) NA else x[1]

    data.frame(
      scenario          = scenario_name,
      prop_species      = get_attr(attrs$prop_species),
      prop_sites        = get_attr(attrs$prop_sites),
      n_sites_requested = get_attr(attrs$n_sites_requested),
      n_sites_actual    = get_attr(attrs$n_sites_actual),
      endemicity_label  = get_attr(attrs$endemicity_label),
      endemism_min      = get_attr(attrs$endemism_min),
      endemism_max      = get_attr(attrs$endemism_max),
      richness_label    = get_attr(attrs$richness_label),
      richness_min      = get_attr(attrs$richness_min),
      richness_max      = get_attr(attrs$richness_max),
      clustering_label  = get_attr(attrs$clustering_label),
      diffusion         = get_attr(attrs$diffusion),
      numberStartSites  = get_attr(attrs$numberStartSites),
      min_weight        = get_attr(attrs$min_weight),
      sampling_evenness = get_attr(attrs$sampling_evenness),
      min_prop_species  = get_attr(attrs$min_prop_species),
      stringsAsFactors  = FALSE
    )
  }) %>%
    dplyr::bind_rows()


  # --- 6. Join summaries and export ---
  message("Joining summaries and exporting...")

  scenario_summary_full <- scenario_summary1 %>%
    dplyr::left_join(sf::st_drop_geometry(scenario_summary2), by = "scenario") %>%
    dplyr::left_join(sampling_attrs, by = "scenario") %>%
    dplyr::mutate(scenario = factor(scenario, levels = scenario_levels, ordered = TRUE)) %>%
    dplyr::arrange(scenario) %>%
    dplyr::mutate(sampling_fraction = sampling_levels[as.integer(scenario)])

  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")

  # Summary CSV
  summary_file <- paste0(output_folder, output_prefix, ".summary_", timestamp, ".csv")
  write.csv(scenario_summary_full, summary_file, row.names = FALSE)
  message("Writing: ", summary_file)

  # Export names combine scenario label and sampling fraction
  export_names <- setNames(
    paste0(scenario_levels, "_", sampling_levels),
    scenario_levels
  )

  # Per-scenario occurrence CSVs
  for (scenario_name in scenario_levels) {
    df <- sampling_results[[scenario_name]]$occurenceData
    outfile <- paste0(
      output_folder, export_names[scenario_name],
      output_prefix, ".summary_.Info_", timestamp, ".csv"
    )
    write.csv(df, outfile, row.names = FALSE)
    message("Writing: ", outfile)
  }

  # Per-scenario community RData files
  for (scenario_name in scenario_levels) {
    df <- comm_list[[scenario_name]]
    outfile <- paste0(
      output_folder, export_names[scenario_name],
      output_prefix, ".sample_", timestamp, ".RData"
    )
    save(df, file = outfile)
    message("Writing: ", outfile)
  }

  message("Done! All files written to: ", output_folder)

  invisible(list(
    scenario_summary_full = scenario_summary_full,
    all_samples           = all_samples,
    all_comm              = all_comm,
    comm_list             = comm_list,
    sample_species_ranges = sample_species_ranges
  ))
}


#' Process and Export Sampling Scenario Results
#'
#' Combines occurrence data across scenarios, computes range-size and community
#' summaries, extracts sampling attributes, and exports everything to CSV/RData.
#'
#' Requires: data.table, dplyr, sf, purrr, ggplot2, e1071
#'
#' @param sampling_results Named list of sampling results (each element has
#'   `$occurenceData` and sampling attributes).
#' @param scenario_levels Character vector of scenario names in order.
#' @param sampling_levels Numeric vector of sampling fractions (same length as
#'   `scenario_levels`).
#' @param true_occ_data Full true occurrence data frame (passed to
#'   `rebuild_community()`).
#' @param output_folder Path to output directory. Created if it doesn't exist.
#' @param output_prefix String used as the base name for all output files.
#' @param obs_per_speciesInCell Numeric scalar recorded in the community summary
#'   (default 1).
#' @param plot Logical; whether to print diagnostic ggplots (default TRUE).
#' @return Invisibly returns a list with `scenario_summary_full`, `all_samples`,
#'   `all_comm`, `comm_list`, and `sample_species_ranges`.
output_sample_results <- function(sampling_results,
                                  scenario_levels,
                                  sampling_levels,
                                  true_occ_data,
                                  output_folder,
                                  output_prefix = "lowtoHighSpeciesSampling",
                                  obs_per_speciesInCell = 1,
                                  plot = TRUE) {
  
  # Ensure output folder exists
  if (!dir.exists(output_folder)) {
    dir.create(output_folder, recursive = TRUE)
    message("Created output folder: ", output_folder)
  }
  
  # --- 0. INTERNAL SANITIZATION ---
  # Safely handle duplicated scenario names (e.g. "low", "low") without breaking
  # old scripts. This forces list names and factor levels to be strictly unique.
  internal_scenario_levels <- make.unique(as.character(scenario_levels), sep = "_run")
  
  if (length(sampling_results) == length(internal_scenario_levels)) {
    names(sampling_results) <- internal_scenario_levels
  } else {
    warning("Length of sampling_results differs from scenario_levels. Extraction may fail.")
  }
  
  # --- 1. Combine occurrence data across scenarios ---
  message("Combining occurrence data across scenarios...")
  
  all_samples <- data.table::rbindlist(
    lapply(sampling_results, `[[`, "occurenceData"),
    idcol = "scenario"
  )
  
  all_samples <- all_samples %>%
    dplyr::mutate(scenario = factor(
      scenario,
      levels = internal_scenario_levels,
      ordered = TRUE
    ))
  
  
  # --- 2. Compute species range summaries ---
  message("Computing species range summaries...")
  
  sample_species_ranges <- all_samples %>%
    sf::st_drop_geometry() %>%
    dplyr::distinct(species, grid_id, scenario) %>%
    dplyr::count(species, scenario, name = "n_cells")
  
  scenario_summary1 <- sample_species_ranges %>%
    dplyr::group_by(scenario) %>%
    dplyr::summarise(
      mean_range_size   = mean(n_cells),
      median_range_size = median(n_cells),
      sd_range_size     = sd(n_cells),
      skew_range_size   = e1071::skewness(n_cells),
      .groups = "drop"
    )
  
  message("Range size summary:")
  print(scenario_summary1)
  
  # Diagnostic plots
  if (plot) {
    message("Generating diagnostic plots...")
    
    p1 <- ggplot2::ggplot(sample_species_ranges,
                          ggplot2::aes(x = log10(n_cells), fill = scenario)) +
      ggplot2::geom_density(alpha = 0.35) +
      ggplot2::theme_bw() +
      ggplot2::labs(x = "log10(range size in grid cells)", y = "Density")
    print(p1)
    
    p2 <- ggplot2::ggplot(sample_species_ranges,
                          ggplot2::aes(x = scenario, y = n_cells)) +
      ggplot2::geom_boxplot() +
      ggplot2::scale_y_log10() +
      ggplot2::theme_bw() +
      ggplot2::labs(y = "Range size (grid cells, log scale)")
    print(p2)
  }
  
  
  # --- 3. Rebuild community matrices ---
  message("Rebuilding community matrices...")
  
  comm_list <- purrr::imap(sampling_results, function(result, scenario_name) {
    rebuild_community2(
      data.table::setDT(true_occ_data),
      data.table::setDT(result$occurenceData),
      scenario_name
    )
  })
  
  all_comm <- data.table::rbindlist(comm_list, idcol = "temp_scenario_id")
  all_comm[, temp_scenario_id := NULL]
  all_comm[, scenario := factor(scenario, levels = internal_scenario_levels, ordered = TRUE)]
  
  
  # --- 4. Compute community summaries ---
  message("Computing community summaries...")
  
  scenario_summary2 <- all_comm %>%
    dplyr::group_by(scenario) %>%
    dplyr::summarise(
      n_species             = dplyr::n_distinct(species),
      n_grid_cells          = dplyr::n_distinct(grid_id),
      total_observations    = dplyr::n(),
      mean_obs_per_species  = dplyr::n() / dplyr::n_distinct(species),
      obs_per_speciesInCell = obs_per_speciesInCell,
      mean_obs_per_cell     = dplyr::n() / dplyr::n_distinct(grid_id),
      mean_rangeSize        = dplyr::n_distinct(species) / dplyr::n_distinct(grid_id),
      .groups = "drop"
    )
  
  
  # --- 5. Extract sampling attributes ---
  message("Extracting sampling attributes...")
  
  sampling_attrs <- purrr::imap(sampling_results, function(result, scenario_name) {
    attrs <- attributes(result)
    get_attr <- function(x) if (is.null(x) || length(x) == 0) NA else x[1]
    
    data.frame(
      scenario          = scenario_name,
      prop_species      = get_attr(attrs$prop_species),
      prop_sites        = get_attr(attrs$prop_sites),
      n_sites_requested = get_attr(attrs$n_sites_requested),
      n_sites_actual    = get_attr(attrs$n_sites_actual),
      endemicity_label  = get_attr(attrs$endemicity_label),
      endemism_min      = get_attr(attrs$endemism_min),
      endemism_max      = get_attr(attrs$endemism_max),
      richness_label    = get_attr(attrs$richness_label),
      richness_min      = get_attr(attrs$richness_min),
      richness_max      = get_attr(attrs$richness_max),
      clustering_label  = get_attr(attrs$clustering_label),
      diffusion         = get_attr(attrs$diffusion),
      numberStartSites  = get_attr(attrs$numberStartSites),
      min_weight        = get_attr(attrs$min_weight),
      sampling_evenness = get_attr(attrs$sampling_evenness),
      min_prop_species  = get_attr(attrs$min_prop_species),
      stringsAsFactors  = FALSE
    )
  }) %>%
    dplyr::bind_rows()
  
  
  # --- 6. Join summaries and export ---
  message("Joining summaries and exporting...")
  
  scenario_summary_full <- scenario_summary1 %>%
    dplyr::left_join(sf::st_drop_geometry(scenario_summary2), by = "scenario") %>%
    dplyr::left_join(sampling_attrs, by = "scenario") %>%
    dplyr::mutate(scenario = factor(scenario, levels = internal_scenario_levels, ordered = TRUE)) %>%
    dplyr::arrange(scenario) %>%
    dplyr::mutate(
      sampling_fraction      = sampling_levels[as.integer(scenario)],
      original_scenario_name = scenario_levels[as.integer(scenario)] # Retain base names for grouping
    ) %>%
    dplyr::relocate(original_scenario_name, .after = scenario)
  
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  
  # Summary CSV
  summary_file <- paste0(output_folder, "/", output_prefix, ".summary_", timestamp, ".csv")
  write.csv(scenario_summary_full, summary_file, row.names = FALSE)
  message("Writing: ", summary_file)
  
  # Export names combine internal scenario label and sampling fraction
  export_names <- setNames(
    paste0(internal_scenario_levels, "_", sampling_levels),
    internal_scenario_levels
  )
  
  # Per-scenario occurrence CSVs
  for (scenario_name in internal_scenario_levels) {
    df <- sampling_results[[scenario_name]]$occurenceData
    outfile <- paste0(
      output_folder, "/", export_names[scenario_name], "_",
      output_prefix, ".summary_.Info_", timestamp, ".csv"
    )
    write.csv(df, outfile, row.names = FALSE)
    message("Writing: ", outfile)
  }
  
  # Per-scenario community RData files
  for (scenario_name in internal_scenario_levels) {
    df <- comm_list[[scenario_name]]
    outfile <- paste0(
      output_folder, "/", export_names[scenario_name], "_",
      output_prefix, ".sample_", timestamp, ".RData"
    )
    save(df, file = outfile)
    message("Writing: ", outfile)
  }
  
  message("Done! All files written to: ", output_folder)
  
  invisible(list(
    scenario_summary_full = scenario_summary_full,
    all_samples           = all_samples,
    all_comm              = all_comm,
    comm_list             = comm_list,
    sample_species_ranges = sample_species_ranges
  ))
}
