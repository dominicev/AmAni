# utils/write_cell_filter_report_pdf.R
# Extracted from AmAni.summarize.R. Sourced by AmAni.summarize.R — relies on
# the packages (tidyverse, sf, patchwork, ggrepel) loaded there.

# Writes <org_out>/cell_filter_report_<organism>_<timestamp>.pdf (2 pages).
# Returns the file path, or NULL if the inputs are unavailable.
write_cell_filter_report_pdf <- function(organism, org_am, org_out, timestamp,
                                         page_w = 10, page_h = 16) {

  cfr_path <- file.path(org_am, "plot_data", "plotdata_cell_filter_report.rds")
  if (!file.exists(cfr_path)) {
    message("  No plotdata_cell_filter_report.rds for '", organism,
            "'; skipping cell-filter report PDF.")
    return(NULL)
  }
  cfr <- readRDS(cfr_path)

  ds_inputs <- load_dataset_summary_figure_inputs(org_am)
  stats_fig <- if (!is.null(ds_inputs)) {
    build_dataset_summary_figure(
      ds_inputs$thresholds, ds_inputs$counts, ds_inputs$dists,
      title = "Dataset before vs after cell filtering", subtitle = NULL)
  }

  maps <- build_cell_filter_map_panels(cfr)
  # 7 maps in a 2-column grid; the delta map is last, leaving one empty slot
  maps_grid <- wrap_plots(maps, ncol = 2) +
    plot_annotation(tag_levels = "A")
  n_map_rows <- ceiling(length(maps) / 2)

  page1 <- if (!is.null(stats_fig)) {
    (wrap_elements(full = stats_fig) / wrap_elements(full = maps_grid)) +
      plot_layout(heights = c(1.5, n_map_rows))
  } else {
    maps_grid
  }
  page1 <- page1 +
    plot_annotation(
      title = paste0(organism_label(organism), ": cell filtering report"),
      theme = theme(plot.title = element_text(size = 16, face = "bold")))

  p_corr <- build_completeness_correlation_plot(cfr)
  page2 <- if (!is.null(p_corr)) {
    (p_corr / plot_spacer()) + plot_layout(heights = c(1, 0.8)) +
      plot_annotation(
        title = paste0(organism_label(organism), ": plant vs animal completeness"),
        theme = theme(plot.title = element_text(size = 16, face = "bold")))
  }

  out_file <- file.path(org_out,
                        paste0("cell_filter_report_", organism, "_", timestamp, ".pdf"))
  pdf(out_file, width = page_w, height = page_h)
  on.exit(dev.off(), add = TRUE)
  print(page1)
  if (!is.null(page2)) {
    print(page2)
  } else {
    message("  Too few overlapping cells for a completeness correlation; page 2 omitted.")
  }
  message("  Wrote: ", out_file)
  out_file
}
