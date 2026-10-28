# utils/scrape_dataset_summary_stats.R
# Extracted from AmAni.summarize.R. Sourced by AmAni.summarize.R — relies on
# the packages (tidyverse, sf, patchwork, ggrepel) loaded there.

# --- Dataset before/after cell-filtering summary figure ---------------------
#
# scrape_dataset_summary_stats() parses a dataset_summary_*.txt report (see
# 2-ModellingAnimals scripts) into the three inputs required by
# build_dataset_summary_figure(): the filter-threshold text card, the
# grid-cell/record/species retention counts, and the per-cell five-number
# summaries (richness, occurrences, completeness) before vs after filtering.
#
# load_dataset_summary_figure_inputs() is a convenience wrapper that finds
# the most recent dataset_summary_*.txt in an organism's animal_model/
# directory and scrapes it directly.

scrape_dataset_summary_stats <- function(ds_lines) {

  # --- Filter thresholds block (between the two section headers) ---
  start_idx <- grep("Filter thresholds", ds_lines)
  before_idx <- grep("Before filtering", ds_lines)
  thresh_lines <- character(0)
  if (length(start_idx) > 0 && length(before_idx) > 0) {
    block <- ds_lines[(start_idx[1] + 1):(before_idx[1] - 1)]
    thresh_lines <- block[grepl(":", block) & !grepl("^\\s*Filter applied:", block)]
  }
  thresholds <- tibble::tibble(
    label = trimws(sub(":.*", "", thresh_lines)),
    value = trimws(sub("^[^:]*:\\s*", "", thresh_lines))
  )

  # --- Split remaining lines into before/after blocks ---
  after_idx <- grep("After filtering", ds_lines)
  before_block <- if (length(after_idx) > 0) ds_lines[1:(after_idx[1] - 1)] else ds_lines
  after_block  <- if (length(after_idx) > 0) ds_lines[after_idx[1]:length(ds_lines)] else character(0)

  grab_num <- function(block, pattern) {
    l <- grep(pattern, block, value = TRUE)
    if (length(l) == 0) return(NA_real_)
    raw <- sub(".*:\\s*", "", l[1])
    m <- regmatches(raw, regexpr("[0-9]+\\.?[0-9]*", raw))
    if (length(m) == 0) NA_real_ else as.numeric(m)
  }

  counts <- tibble::tribble(
    ~metric,      ~before,                                        ~after,
    "Grid cells", grab_num(before_block, "^\\s*Grid cells:"),      grab_num(after_block, "^\\s*Grid cells:"),
    "Records",    grab_num(before_block, "^\\s*Total records:"),  grab_num(after_block, "^\\s*Total records:"),
    "Species",    grab_num(before_block, "^\\s*Total species:"),  grab_num(after_block, "^\\s*Total species:")
  )

  # --- Per-cell five-number-summary lines, e.g.:
  # "  Richness/cell:     min=1, Q1=1, med=2, mean=6.3, Q3=7, max=132" ---
  parse_stat_line <- function(block, label) {
    l <- grep(paste0("^\\s*", label), block, value = TRUE)
    if (length(l) == 0) return(NULL)
    raw <- sub(".*:\\s*", "", l[1])
    parts <- strsplit(raw, ",\\s*")[[1]]
    keys <- trimws(sub("=.*", "", parts))
    vals <- as.numeric(sub(".*=", "", parts))
    setNames(vals, keys)
  }

  build_dist_rows <- function(variable, label) {
    b <- parse_stat_line(before_block, label)
    a <- parse_stat_line(after_block, label)
    dplyr::bind_rows(
      if (!is.null(b)) tibble::tibble(
        variable = variable, stage = "Before",
        min = b[["min"]], q1 = b[["Q1"]], med = b[["med"]],
        q3 = b[["Q3"]], max = b[["max"]], mean = b[["mean"]]),
      if (!is.null(a)) tibble::tibble(
        variable = variable, stage = "After",
        min = a[["min"]], q1 = a[["Q1"]], med = a[["med"]],
        q3 = a[["Q3"]], max = a[["max"]], mean = a[["mean"]])
    )
  }

  empty_dists <- tibble::tibble(
    variable = character(), stage = character(),
    min = double(), q1 = double(), med = double(),
    q3 = double(), max = double(), mean = double()
  )

  dists <- dplyr::bind_rows(
    build_dist_rows("Richness", "Richness/cell:"),
    build_dist_rows("Occurrences", "Occurrences/cell:"),
    build_dist_rows("Completeness", "Completeness/cell:")
  )
  # bind_rows() of all-NULL inputs drops columns entirely; keep the schema
  # explicit so downstream mutate()/factor() calls don't silently resolve
  # `variable`/`stage` against some unrelated object in the environment.
  if (nrow(dists) == 0) dists <- empty_dists

  list(thresholds = thresholds, counts = counts, dists = dists)
}
