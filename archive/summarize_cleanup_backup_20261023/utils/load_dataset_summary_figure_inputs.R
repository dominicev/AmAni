# utils/load_dataset_summary_figure_inputs.R
# Extracted from AmAni.summarize.R. Sourced by AmAni.summarize.R — relies on
# the packages (tidyverse, sf, patchwork, ggrepel) loaded there.

load_dataset_summary_figure_inputs <- function(am_dir) {
  ds_file <- list.files(am_dir, pattern = "^dataset_summary_.*\\.txt$", full.names = TRUE)
  if (length(ds_file) == 0) return(NULL)
  ds_file <- ds_file[which.max(file.mtime(ds_file))]
  scrape_dataset_summary_stats(readLines(ds_file))
}
