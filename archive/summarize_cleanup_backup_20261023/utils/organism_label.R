# utils/organism_label.R
# Extracted from AmAni.summarize.R. Sourced by AmAni.summarize.R — relies on
# the packages (tidyverse, sf, patchwork, ggrepel) loaded there.

# --- Organism display labels --------------------------------------------------
# Maps internal directory names to human-readable labels for plots.
organism_label <- function(x) {
  dplyr::case_when(
    x == "cockroaches" ~ "Cockroaches.Curated",
    x == "allRoaches"  ~ "cockroaches",
    TRUE ~ paste0(toupper(substr(x, 1, 1)), substr(x, 2, nchar(x)))
  )
}
