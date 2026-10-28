# utils/collect_dataset_distributions.R
# Sourced by AmAni.summarize.R — relies on the packages (tidyverse, sf,
# patchwork, ggrepel) loaded there, on organism_label() from
# utils/organism_label.R, on load_dataset_summary_figure_inputs() from
# utils/load_dataset_summary_figure_inputs.R, and on the private helpers in
# utils/pub_data_helpers.R.

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
