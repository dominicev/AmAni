# utils/plot_spatial_extent_overlap.R
# Extracted from AmAni.summarize.R. Sourced by AmAni.summarize.R — relies on
# the packages (tidyverse, sf, patchwork, ggrepel) loaded there.

# =============================================================================
# --- 7. Spatial extent overlap comparison ------------------------------------
# =============================================================================
#
# plot_spatial_extent_overlap() loads per-organism spatial extent distributions
# (animal + plant) and produces a faceted comparison plot for the organisms
# with the smallest animal spatial extents.
#
# Requires plotdata_spatial_extent_distributions.rds (saved by the pipeline
# after the next re-run) or falls back to plotdata_plant_extent_model.rds
# (plant data only — available now for most organisms).

plot_spatial_extent_overlap <- function(include = NULL, group = NULL,
                                        max_cells = 10,
                                        output_dir = NULL, save_pdf = TRUE) {

  if (!exists("project_root")) source("AmAni.manager.R")
  if (is.null(output_dir)) output_dir <- file.path(project_root, "output")

  all_dirs <- list.dirs(output_dir, full.names = FALSE, recursive = FALSE)
  all_dirs <- all_dirs[!all_dirs %in% c("cross_organism_summary", "shared")]

  if (!is.null(group) && group %in% names(ORGANISM_GROUPS)) {
    include <- union(include, ORGANISM_GROUPS[[group]])
  }
  if (!is.null(include)) {
    all_dirs <- all_dirs[all_dirs %in% include]
  }

  # Load combined distributions (preferred) or plant-only fallback
  combined_data <- purrr::map_dfr(all_dirs, function(org) {
    pd_dir <- file.path(output_dir, org, "animal_model", "plot_data")

    # Preferred: saved combined distribution (animal + plant)
    combined_path <- file.path(pd_dir, "plotdata_spatial_extent_distributions.rds")
    if (file.exists(combined_path)) {
      d <- readRDS(combined_path)
      return(dplyr::mutate(d, organism = org))
    }

    # Fallback: plant-only from extent model
    ext_path <- file.path(pd_dir, "plotdata_plant_extent_model.rds")
    if (file.exists(ext_path)) {
      d <- readRDS(ext_path)
      return(dplyr::mutate(d$plant_df, organism = org))
    }
    NULL
  })

  if (nrow(combined_data) == 0) {
    message("No spatial extent data found.")
    return(invisible(NULL))
  }

  # Show all matched organisms, filtered to 1..max_cells on x-axis
  focus_orgs <- unique(combined_data$organism)

  plot_data <- combined_data |>
    dplyr::filter(num_grid_cells <= max_cells) |>
    dplyr::mutate(
      organism = organism_label(organism),
      num_grid_cells = factor(num_grid_cells, levels = seq_len(max_cells))
    )

  n_orgs <- dplyr::n_distinct(plot_data$organism)

  p <- ggplot(plot_data, aes(x = num_grid_cells, y = n, fill = group)) +
    geom_col(position = "dodge", color = "grey30", linewidth = 0.2) +
    scale_fill_manual(
      values = c("Plant species" = "#2d6a4f", "Animal species" = "#457b9d"),
      na.value = "grey60", name = NULL
    ) +
    facet_wrap(~organism, scales = "free_y",
               ncol = min(n_orgs, 4)) +
    labs(
      title = paste0("Spatial Extent Distribution (1\u2013", max_cells, " grid cells)"),
      subtitle = if (!is.null(group)) paste0("Group: ", group) else NULL,
      x = "Number of grid cells occupied",
      y = "Number of species"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(size = 7),
      legend.position = "top"
    )

  if (save_pdf) {
    summary_dir <- file.path(output_dir, "cross_organism_summary")
    dir.create(summary_dir, showWarnings = FALSE, recursive = TRUE)
    ts <- format(Sys.time(), "%Y%m%d_%H%M%S")
    out_file <- file.path(summary_dir,
                          paste0("spatial_extent_overlap_", ts, ".pdf"))
    pdf(out_file, width = 12, height = 8)
    print(p)
    dev.off()
    message("Wrote: ", out_file)
  }

  p
}
