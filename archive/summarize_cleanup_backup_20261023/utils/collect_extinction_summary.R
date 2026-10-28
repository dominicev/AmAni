# utils/collect_extinction_summary.R
# Sourced by AmAni.summarize.R — relies on the packages (tidyverse, sf,
# patchwork, ggrepel) loaded there, on organism_label() from
# utils/organism_label.R, and on the private helpers in
# utils/pub_data_helpers.R.

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
