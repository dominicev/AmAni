# utils/pub_data_helpers.R
# Sourced by AmAni.summarize.R — relies on the packages (tidyverse, sf,
# patchwork, ggrepel) loaded there.
#
# Private helpers shared by the publication-figure data collectors
# (collect_extinction_summary.R, collect_organism_summary.R,
# collect_dataset_distributions.R, collect_iucn.R, collect_pub_data.R).

# Most recent file matching a pattern in a directory, or NULL if none exists.
# Files are timestamped, so "most recent" picks the latest run.
.latest_file <- function(dir, pattern) {
  f <- list.files(dir, pattern = pattern, full.names = TRUE)
  if (length(f) == 0) return(NULL)
  f[which.max(file.mtime(f))]
}

# Extract a value from text-report lines by regex pattern (same semantics as
# the extract_val() helper inside summarize_all()).
.extract_val <- function(lines, pattern, as_type = "numeric") {
  line <- grep(pattern, lines, value = TRUE, ignore.case = TRUE)
  if (length(line) == 0) return(NA)
  raw <- trimws(sub(".*:\\s*", "", line[1]))
  if (as_type == "numeric") {
    m <- regmatches(raw, regexpr("[0-9]+\\.?[0-9]*", raw))
    if (length(m) == 0) return(NA_real_)
    as.numeric(m)
  } else {
    raw
  }
}

# All organism output directories (dirs with an animal_model/ subdirectory).
# This skips non-organism folders like "shared", "comparisons", "insects",
# "vertebrates" and "cross_organism_summary" automatically.
.organism_dirs <- function(output_dir = "output") {
  dirs <- list.dirs(output_dir, full.names = FALSE, recursive = FALSE)
  dirs[dir.exists(file.path(output_dir, dirs, "animal_model"))]
}
