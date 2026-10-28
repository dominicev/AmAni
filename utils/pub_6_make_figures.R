# utils/pub_6_make_figures.R — Publication figures, step 6: build and save all
#
#   make_publication_figures()   collect data once, draw every publication
#                                figure and write PDF (+ PNG preview) to
#                                <output_dir>/<results_dir>/
#
# Each figure is drawn inside tryCatch so one failure (e.g. no network for a
# PhyloPic, rgl unavailable for the 3D map) doesn't stop the rest; failures
# are reported in the returned table.
#
# Requires the pub_0 ... pub_5 files to be sourced (done by AmAni.summarize.R).

make_publication_figures <- function(include = PUB_SETTINGS$include,
                                     group = PUB_SETTINGS$group,
                                     output_dir = PUB_SETTINGS$output_dir,
                                     dir = file.path(output_dir, PUB_SETTINGS$results_dir),
                                     figures = c("fig1", "fig2", "fig3"),
                                     value_col = "gamma_lower",
                                     include_3d = TRUE) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  w <- PUB_SETTINGS$fig_width_double
  message("Collecting publication data ...")
  pub <- collect_pub_data(include = include, group = group, output_dir = output_dir)

  results <- list()
  build <- function(name, fig, expr, height = NULL, width = w) {
    if (!fig %in% figures) return(invisible())
    message("Drawing ", name, " ...")
    res <- tryCatch({
      save_pub_fig(force(expr), name, dir = dir, width = width, height = height)
      "ok"
    }, error = function(e) paste("FAILED:", conditionMessage(e)))
    if (res != "ok") warning(name, ": ", res, call. = FALSE)
    results[[name]] <<- res
  }

  # --- Figure 1: sampling and filtering ---------------------------------------
  if ("fig1" %in% figures) {
    retained <- collect_retained_cells(output_dir)
    raw_comp <- collect_plant_completeness_raw()
    stats <- retained_cell_stats(retained, raw_comp)
    build("fig1_retained_cells_by_organism", "fig1",
          fig_retained_cells(retained, raw_comp), height = 120)
    build("fig1_retained_cells_combined", "fig1",
          fig_retained_cells_combined(retained, raw_comp), height = 110)
    build("fig1_retained_cells_combined_labelled", "fig1",
          label_retained_cells(fig_retained_cells_combined(retained, raw_comp), stats),
          height = 110)
    build("fig1_sampling_retention", "fig1", fig_sampling_retention(pub), height = 100)
    build("fig1_cell_distributions", "fig1", fig_cell_distributions(pub), height = 110)
    if (include_3d) res <- tryCatch({
      fig_retained_cells_3d(retained, raw_comp,
                            file = file.path(dir, "fig1_retained_cells_3d.png"),
                            labels = TRUE)
      "ok"
    }, error = function(e) paste("FAILED:", conditionMessage(e)))
    if (include_3d) {
      if (res != "ok") warning("fig1_retained_cells_3d: ", res, call. = FALSE)
      results[["fig1_retained_cells_3d"]] <- res
    }
  }

  # --- Figure 2: gamma richness -------------------------------------------------
  build("fig2_gamma_treemap", "fig2", fig_gamma_treemap(pub, value_col = value_col),
        height = w)
  build("fig2_gamma_uncertainty_treemap", "fig2",
        fig_gamma_uncertainty_treemap(pub, organism_col = value_col), height = w)

  # --- Figure 3: extinction and IUCN ----------------------------------------------
  build("fig3_extinction", "fig3", fig_extinction(pub, value_col), height = 175)
  build("fig3a_species_lost", "fig3", fig_species_lost(pub, value_col), height = 90)
  build("fig3b_extinction_sensitivity", "fig3",
        fig_extinction_sensitivity(pub, value_col), height = 90)
  build("fig3c_iucn_ribbons", "fig3", fig_iucn_stack(pub, value_col), height = 80)

  out <- tibble::tibble(figure = names(results), status = unlist(results))
  message("Publication figures written to ", normalizePath(dir))
  print(out)
  invisible(list(pub = pub, status = out, dir = dir))
}
