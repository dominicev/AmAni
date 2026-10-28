# utils/collect_iucn.R
# Sourced by AmAni.summarize.R — relies on the packages (tidyverse, sf,
# patchwork, ggrepel) loaded there, on organism_label() from
# utils/organism_label.R, and on the private helpers in
# utils/pub_data_helpers.R.

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
