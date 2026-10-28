# 00_cleanGBIF.R — Pre-processing filters for standard GBIF occurrence CSVs
#
# Reads a raw GBIF download CSV file, applies quality filters, and writes
# a cleaned CSV to the same directory. The cleaned file can then be passed
# to run_animal_model() as the data_file.
#
# Usage (from AmAni.run.R or interactively):
#   source("2-ModellingAnimals/scripts/00_cleanGBIF.R")
#   clean_gbif("2-ModellingAnimals/data/beetle data/0011381-250402121839773.csv")
#
# Filters applied:
#   1.  Remove records without species-level identification
#   2.  Remove records with blank or NA species names
#   3.  Remove species names that look like undetermined IDs (e.g. "sp.", "cf.", "aff.",
#         "agg.", "sensu lato", "s.l.", "group", "complex", "?", "×").
#         Trailing numeric morphotype suffixes (e.g. "Aus bus 1") are intentionally kept.
#   4.  Remove records missing coordinates
#   5.  Remove records with coordinate uncertainty > max_uncertainty_m (default 10 km)
#   6.  Remove ABSENT records
#   7.  Remove records flagged with high-impact GBIF issues
#   8.  Remove cultivated / introduced / managed records (establishmentMeans)
#   9.  Remove records with implausible collection years (< min_year or > current year)
#   10. Standardize species names (trim whitespace, collapse internal spaces)
#   11. Remove exact spatial duplicates (same species + identical lat/lon)

library(data.table)

clean_gbif <- function(
    input_file,
    output_file = NULL,
    max_uncertainty_m = 10000,
    min_year = 1800,
    verbose = TRUE
) {
  
  if (is.null(output_file)) {
    # For zip inputs, derive the output name from the CSV inside the archive
    if (tolower(tools::file_ext(input_file)) == "zip") {
      #csv_name <- grep("\\.csv$", unzip(input_file, list = TRUE)$Name, value = TRUE)[1] ##original, works
      csv_name <- grep("\\.(csv|tsv)$", unzip(input_file, list = TRUE)$Name, ignore.case = TRUE, value = TRUE) #fix to add TSV
      output_file <- file.path(dirname(input_file),
                               sub("\\.csv$", "_cleaned.csv", basename(csv_name)))
    } else {
      output_file <- sub("\\.csv$", "_cleaned.csv", input_file)
    }
  }
  
  if (verbose) message("Reading: ", input_file)
  
  # Support reading directly from zip archives
  fread_auto <- function(path, ...) {
    ext <- tolower(tools::file_ext(path))
    # GBIF .tsv downloads are unquoted; CSVs quote fields that contain commas
    sep_for <- function(e) if (e == "tsv") "\t" else ","
    quote_for <- function(e) if (e == "tsv") "" else "\""
    
    if (ext == "zip") {
      csv_name <- grep("\\.(csv|tsv)$", unzip(path, list = TRUE)$Name, value = TRUE)
      if (length(csv_name) == 0) stop("No CSV found inside: ", path, call. = FALSE)
      if (length(csv_name) > 1) csv_name <- csv_name[1]
      
      # Check inner extension for zipped files
      inner_ext <- tolower(tools::file_ext(csv_name))
      
      cmd <- sprintf('unzip -p "%s" "%s"', normalizePath(path, winslash = "/"), csv_name)
      data.table::fread(cmd = cmd, sep = sep_for(inner_ext), quote = quote_for(inner_ext), ...)
    } else {
      data.table::fread(file = path, sep = sep_for(ext), quote = quote_for(ext), ...)
    }
  }
  
  #Peek at the first row of the file to see what headers actually exist
  file_headers <- names(fread_auto(input_file, nrows = 0))
  print(paste0("There are ", toString( length(file_headers)), " header columns in the raw data file."))
  
  
  # Define the columns your model ideally wants
  desired_columns <- c("gbifID", "species", "taxonRank", "scientificName",
                       "decimalLatitude", "decimalLongitude",
                       "coordinateUncertaintyInMeters", "occurrenceStatus", "issue",
                       "establishmentMeans", "year")
  
  #Only select the columns that are actually present in the file
  columns_to_keep <- intersect(desired_columns, file_headers)
  
  #Read the file safely
  raw <- fread_auto(input_file, select = columns_to_keep)
  
  # Fall back to scientificName when no species column is available
  if (!"species" %in% names(raw) && "scientificName" %in% names(raw)) {
    if (verbose) message("  'species' column not found; using scientificName as species")
    data.table::setnames(raw, "scientificName", "species")
  } else if ("species" %in% names(raw) && "scientificName" %in% names(raw)) {
    # Fill blank/NA species values from scientificName
    raw[is.na(species) | species == "", species := scientificName]
  }
  
  n_start <- nrow(raw)
  if (verbose) message("  Raw records: ", n_start)
  
  # --- Filter log helper ---
  filter_log <- function(dt, label) {
    n <- nrow(dt)
    if (verbose) message("  After ", label, ": ", n,
                         " (-", n_start - n, " cumulative)")
    n_start <<- n  # update for next delta (cumulative)
    dt
  }
  
  # 1. Species-level identification only
  dt <- raw[taxonRank == "SPECIES"]
  dt <- filter_log(dt, "taxonRank == SPECIES")
  
  # 2. Non-empty species name
  dt <- dt[!is.na(species) & species != ""]
  dt <- filter_log(dt, "non-empty species name")
  
  # 3. Remove undetermined / provisional IDs
  #    Covers: sp./spp./cf./aff./nr., bare trailing "sp", aggregates (agg.),
  #    sensu lato (sensu lato / s.l. / sl), informal group / complex designators,
  #    uncertain-ID markers (?) and hybrid markers (×).
  #    Trailing numeric morphotype suffixes like "Aus bus 1" are intentionally kept —
  #    none of these patterns match a trailing number.
  bad_pattern <- paste0(
    "\\b(sp\\.|spp\\.|cf\\.|aff\\.|nr\\.)\\b",  # original undetermined epithets
    "|\\bsp$",                                    # bare "sp" at end of name
    "|\\bagg\\.?\\b",                             # aggregation (agg / agg.)
    "|\\bsensu\\s+lato\\b",                       # sensu lato (full form)
    "|\\bs\\.\\s*l\\.?\\b",                       # s.l. / s. l.
    "|\\bsl\\b",                                  # sl (abbreviated sensu lato)
    "|\\bgroup\\b",                               # informal species group
    "|\\bcomplex\\b",                             # informal species complex
    "|\\?",                                       # uncertain identification marker
    "|\u00d7"                                     # × hybrid marker (U+00D7)
  )
  dt <- dt[!grepl(bad_pattern, species, ignore.case = TRUE)]
  dt <- filter_log(dt, "remove undetermined / uncertain IDs")
  
  # 4. Valid coordinates
  dt <- dt[!is.na(decimalLatitude) & !is.na(decimalLongitude)]
  dt <- filter_log(dt, "valid coordinates")
  
  # 5. Coordinate uncertainty threshold
  #    Records with NA uncertainty are kept (many museum records lack this field)
  dt <- dt[is.na(coordinateUncertaintyInMeters) |
             coordinateUncertaintyInMeters <= max_uncertainty_m]
  dt <- filter_log(dt, paste0("uncertainty <= ", max_uncertainty_m, " m or NA"))
  
  # 6. Presence records only
  dt <- dt[occurrenceStatus == "PRESENT" | occurrenceStatus == ""]
  dt <- filter_log(dt, "PRESENT records")
  
  # 7. Remove records with high-impact GBIF issues
  #    These flags indicate the coordinates are likely wrong
  bad_issues <- c("ZERO_COORDINATE",
                  "COORDINATE_INVALID",
                  "COUNTRY_COORDINATE_MISMATCH",
                  "COORDINATE_OUT_OF_RANGE")
  has_bad <- Reduce(`|`, lapply(bad_issues, function(flag) grepl(flag, dt$issue)))
  dt <- dt[!has_bad]
  dt <- filter_log(dt, "remove bad GBIF issue flags")
  
  # 8. Remove cultivated / introduced / managed records
  #    Only applied when the field is present and non-empty; NA / blank records
  #    (the vast majority for invertebrates) are kept.
  if ("establishmentMeans" %in% names(dt)) {
    remove_means <- c("MANAGED", "CULTIVATED", "INTRODUCED", "INVASIVE")
    dt <- dt[!(establishmentMeans %in% remove_means)]
    dt <- filter_log(dt, "remove cultivated/introduced (establishmentMeans)")
  }
  
  # 9. Remove records with implausible collection years
  #    Keeps NA years (many older records lack this field).
  if ("year" %in% names(dt)) {
    current_year <- as.integer(format(Sys.Date(), "%Y"))
    dt <- dt[is.na(year) | year == 0L |
               (year >= min_year & year <= current_year)]
    dt <- filter_log(dt, paste0("plausible year (", min_year, "\u2013", current_year, ") or NA"))
  }
  
  # 10. Standardize species names (before duplicate removal so comparisons are clean)
  dt[, species := trimws(species)]
  dt[, species := gsub("\\s+", " ", species)]
  
  # 11. Remove exact spatial duplicates (same species + identical lat/lon)
  #     Retains the first occurrence of each species-location combination.
  dt <- unique(dt, by = c("species", "decimalLatitude", "decimalLongitude"))
  dt <- filter_log(dt, "remove exact spatial duplicates")
  
  # --- Summary ---
  if (verbose) {
    message("\n  Final: ", nrow(dt), " records, ",
            uniqueN(dt$species), " species")
    message("  Writing: ", output_file)
  }
  
  # Write cleaned file (keep only the columns 01_GBIFAnalysis.R needs)
  out <- dt[, .(species, decimalLatitude, decimalLongitude, occurrenceStatus)]
  fwrite(out, output_file)
  
  invisible(output_file)
}
