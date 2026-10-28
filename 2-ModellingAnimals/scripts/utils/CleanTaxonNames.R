# CleanTaxonNames.R — Taxonomic name standardisation for DarwinCore occurrence data
#
# Queries the GBIF backbone taxonomy to:
#   - Correct orthographic / spelling variants   (matchType == "VARIANT")
#   - Replace synonyms with their accepted names (status == "SYNONYM", any matchType)
#
# Spurious cross-taxon swaps (different genus AND different epithet) are NOT
# applied automatically; they are written to a review file for manual inspection.
# Duplicate occurrence records are intentionally left untouched.
#
# Usage:
#   source("2-ModellingAnimals/scripts/utils/CleanTaxonNames.R")
#
#   # From a data.frame (e.g. output of clean_gbif()):
#   cleaned <- clean_taxon_names(occ_df, kingdom = "Animalia")
#
#   # From a CSV path:
#   cleaned <- clean_taxon_names("path/to/cleaned.csv", kingdom = "Animalia",
#                                output_file = "path/to/taxa_fixed.csv")
#
# Required packages: tidyverse, rgbif, stringdist

library(tidyverse)
library(rgbif)
library(stringdist)

#' Standardise taxon names in a DarwinCore occurrence data frame
#'
#' @param occ         A data.frame with DarwinCore columns **or** a path to a
#'                    CSV file that will be read with [readr::read_csv()].
#' @param species_col Name of the column holding species binomials.
#'                    Defaults to \code{"species"}.
#' @param kingdom     Optional kingdom string passed to the GBIF backbone query
#'                    (e.g. \code{"Animalia"}, \code{"Plantae"}). Restricts
#'                    fuzzy matching so names aren't swapped across kingdoms.
#'                    \code{NULL} (default) applies no kingdom filter.
#' @param output_file Path for the cleaned CSV. Only written when non-NULL.
#'                    If \code{occ} was supplied as a file path and
#'                    \code{output_file} is NULL, a \code{_taxaFixed} suffix is
#'                    added automatically.
#' @param review_file Path for the review CSV (spurious matches + unresolved
#'                    names). Defaults to \code{NULL} (not written).
#' @param bucket_size Number of names sent per GBIF API request. Default 50.
#' @param sleep       Seconds to pause between API batches. Default 1.
#' @param verbose     Print progress messages? Default \code{TRUE}.
#'
#' @return The corrected data.frame (invisibly when \code{output_file} is set).
clean_taxon_names <- function(
    occ,
    species_col  = "species",
    kingdom      = NULL,
    output_file  = NULL,
    review_file  = NULL,
    bucket_size  = "auto",
    sleep        = "auto",
    max_retries  = 3,
    verbose      = TRUE
) {

  # ── 1. Load data ──────────────────────────────────────────────────────────
  input_was_path <- is.character(occ) && length(occ) == 1
  if (input_was_path) {
    if (verbose) message("Reading: ", occ)
    occ <- readr::read_csv(occ, show_col_types = FALSE)

    # Derive a sensible output path when not provided
    if (is.null(output_file)) {
      output_file <- sub("(\\.csv)?$", "_taxaFixed.csv", occ, ignore.case = TRUE)
    }
  }

  if (!is.data.frame(occ)) stop("`occ` must be a data.frame or a path to a CSV file.")
  if (!species_col %in% names(occ)) {
    stop("Column '", species_col, "' not found in `occ`. ",
         "Adjust `species_col` to match your data.")
  }

  n_records <- nrow(occ)
  if (verbose) message("Records in: ", n_records)

  # ── 1b. Sanity-check the species column ───────────────────────────────────
  na_frac <- mean(is.na(occ[[species_col]]) | occ[[species_col]] == "")
  if (na_frac > 0.1) {
    warning(
      sprintf(
        "%.0f%% of '%s' values are NA or empty. ",
        na_frac * 100, species_col
      ),
      "This usually means either:\n",
      "  (a) the wrong column name was passed via `species_col`, or\n",
      "  (b) the data has not been run through clean_gbif() yet and some rows\n",
      "      store names in a different column (e.g. 'scientificName').\n",
      "Those rows will pass through unchanged with NA species.",
      call. = FALSE
    )
  }

  # ── 2. Extract unique, non-empty species names ────────────────────────────
  unique_names <- occ |>
    dplyr::distinct(.data[[species_col]]) |>
    dplyr::filter(!is.na(.data[[species_col]]),
                  .data[[species_col]] != "") |>
    dplyr::pull(species_col)

  if (verbose) message("Unique names to query: ", length(unique_names))

  # ── 3. Query GBIF backbone taxonomy ──────────────────────────────────────
  name_data <- tibble::tibble(name = unique_names)
  if (!is.null(kingdom)) name_data$kingdom <- kingdom

  # ── 3a. Dynamic bucket sizing based on pool size ─────────────────────────
  n_names <- length(unique_names)

  if (identical(bucket_size, "auto")) {
    bucket_size <- if (n_names <= 50)   max(5L, n_names)
                   else if (n_names <= 200)  25L
                   else if (n_names <= 500)  40L
                   else if (n_names <= 1000) 50L
                   else                      75L
  }

  if (identical(sleep, "auto")) {
    sleep <- if (n_names <= 50)   0.5
             else if (n_names <= 200)  1.0
             else if (n_names <= 500)  2.0
             else if (n_names <= 1000) 3.0
             else                      4.0
  }

  n_batches <- ceiling(n_names / bucket_size)
  est_secs  <- n_batches * sleep

  if (verbose) {
    message(sprintf(
      "  GBIF query: %d names \u2192 %d batches (bucket_size=%d, sleep=%.1fs, est. ~%.0fs)",
      n_names, n_batches, bucket_size, sleep, est_secs
    ))
  }

  # ── 3b. Query with automatic retry + backoff ────────────────────────────
  gbif_match <- NULL
  for (attempt in seq_len(max_retries)) {
    gbif_match <- tryCatch(
      rgbif::name_backbone_checklist(
        name_data   = name_data,
        bucket_size = bucket_size,
        sleep       = sleep
      ),
      error = function(e) e
    )

    if (!inherits(gbif_match, "error")) {
      if (verbose && attempt > 1) message("  \u2713 Succeeded on attempt ", attempt)
      break
    }

    if (attempt == max_retries) {
      stop("GBIF backbone query failed after ", max_retries, " attempts. ",
           "Last error: ", conditionMessage(gbif_match), call. = FALSE)
    }

    # Back off: halve bucket size (min 5), double sleep
    bucket_size <- max(5L, bucket_size %/% 2L)
    sleep       <- sleep * 2
    n_batches   <- ceiling(n_names / bucket_size)
    est_secs    <- n_batches * sleep

    if (verbose) {
      message(sprintf(
        "  \u2717 Attempt %d failed (%s). Retrying with bucket_size=%d, sleep=%.1fs (est. ~%.0fs)...",
        attempt, conditionMessage(gbif_match), bucket_size, sleep, est_secs
      ))
    }
  }

  # ── 4. Classify matches ───────────────────────────────────────────────────
  # Four outcomes per verbatim name:
  #   (a) EXACT + ACCEPTED  → already correct, no change.
  #   (b) EXACT + SYNONYM   → replace with the accepted species name.
  #   (c) VARIANT + any     → potential spelling fix or synonym; evaluate
  #                            for spurious cross-taxon swaps before applying.
  #   (d) HIGHERRANK / NONE → cannot resolve to species level; flag for review.

  match_summary <- gbif_match |>
    dplyr::transmute(
      verbatim_name,
      matchType,
      confidence,
      rank,
      gbif_species = species,   # accepted species name from GBIF
      status
    )

  # (a) Exact + accepted: no action needed
  exact_accepted <- match_summary |>
    dplyr::filter(matchType == "EXACT",
                  status == "ACCEPTED",
                  rank == "SPECIES")

  # (b) Exact + synonym: straightforward accepted-name swap
  exact_synonym <- match_summary |>
    dplyr::filter(matchType == "EXACT",
                  status == "SYNONYM",
                  rank == "SPECIES",
                  !is.na(gbif_species),
                  gbif_species != verbatim_name)

  # (c) Variant matches: need spurious-match filtering
  variant_candidates <- match_summary |>
    dplyr::filter(matchType == "VARIANT",
                  rank == "SPECIES",
                  !is.na(gbif_species),
                  gbif_species != verbatim_name) |>
    dplyr::mutate(
      genus_in      = word(verbatim_name, 1),
      genus_out     = word(gbif_species, 1),
      epithet_in    = word(verbatim_name, 2),
      epithet_out   = word(gbif_species, 2),
      genus_changed = genus_in != genus_out,
      epithet_dist  = stringdist::stringdist(epithet_in, epithet_out, method = "lv"),
      # Spurious = genus changes AND epithet is different — likely a cross-taxon match.
      # A pure spelling fix (epithet_dist <= 1) or a genus respelling while keeping
      # the same epithet (synonym move) are both trustworthy.
      spurious = genus_changed & epithet_dist > 1
    )

  variant_ok       <- dplyr::filter(variant_candidates, !spurious)
  variant_spurious <- dplyr::filter(variant_candidates,  spurious)

  # (d) No species-level resolution
  no_species_match <- match_summary |>
    dplyr::filter(
      matchType %in% c("NONE", "HIGHERRANK") |
        rank != "SPECIES" |
        is.na(gbif_species)
    ) |>
    dplyr::filter(!verbatim_name %in% exact_accepted$verbatim_name,
                  !verbatim_name %in% exact_synonym$verbatim_name,
                  !verbatim_name %in% variant_candidates$verbatim_name)

  # ── 5. Report ─────────────────────────────────────────────────────────────
  if (verbose) {
    message("\n── GBIF match summary ──────────────────────────────────")
    message("  Exact + accepted (no change):      ", nrow(exact_accepted))
    message("  Exact + synonym  (name replaced):  ", nrow(exact_synonym))
    message("  Variant fixes applied:             ", nrow(variant_ok))
    message("  Variant matches flagged (spurious):", nrow(variant_spurious))
    message("  No species-level match (kept as-is):", nrow(no_species_match))
    message("────────────────────────────────────────────────────────\n")
  }

  # ── 6. Build replacement lookup (exact synonyms + safe variants) ──────────
  replacements <- dplyr::bind_rows(
    exact_synonym |> dplyr::select(verbatim_name, gbif_species),
    variant_ok    |> dplyr::select(verbatim_name, gbif_species)
  )

  # ── 7. Apply replacements ─────────────────────────────────────────────────
  occ_clean <- occ |>
    dplyr::left_join(
      replacements |> dplyr::rename(!!species_col := verbatim_name,
                                    .gbif_accepted = gbif_species),
      by = species_col
    ) |>
    dplyr::mutate(
      !!species_col := dplyr::coalesce(.gbif_accepted, .data[[species_col]])
    ) |>
    dplyr::select(-.gbif_accepted)

  n_records_out <- nrow(occ_clean)
  if (verbose) {
    message("Records out: ", n_records_out,
            " (unchanged — no duplicates removed)")
    n_replaced <- occ_clean |>
      dplyr::count(.data[[species_col]]) |>
      nrow()
    message("Unique accepted names after cleaning: ", n_replaced)
  }

  # ── 8. Write outputs ───────────────────────────────────────────────────────
  if (!is.null(output_file)) {
    if (verbose) message("Writing cleaned data: ", output_file)
    readr::write_csv(occ_clean, output_file)
  }

  if (!is.null(review_file)) {
    review <- dplyr::bind_rows(
      variant_spurious |>
        dplyr::transmute(
          verbatim_name,
          gbif_suggestion = gbif_species,
          reason = "spurious_cross_taxon_variant",
          confidence,
          genus_changed,
          epithet_dist
        ),
      no_species_match |>
        dplyr::transmute(
          verbatim_name,
          gbif_suggestion = NA_character_,
          reason = paste0("no_species_match_", matchType),
          confidence,
          genus_changed = NA,
          epithet_dist  = NA_real_
        )
    ) |>
      dplyr::arrange(reason, verbatim_name)

    if (verbose) message("Writing review file:  ", review_file,
                         "  (", nrow(review), " names to inspect)")
    readr::write_csv(review, review_file)
  }

  if (!is.null(output_file)) invisible(occ_clean) else occ_clean
}
