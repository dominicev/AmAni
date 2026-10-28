## AmAni.cleanUp.R
## Purpose: minimize disk space used by larger outputs and inputs
##          of the AmAni pipeline.

## ---- Libraries -------------------------------------------------------------

## fs: file/folder size inspection and deletion
library(fs)

## (Add others here as needed, e.g.)
## library(archive)   # compress/decompress files

## ---- Working directory -----------------------------------------------------

setwd("~/Desktop/AmAni")
getwd()

## ---- Helpers ---------------------------------------------------------------

# Summarize space usage of the main output/input folders (recursive contents,
# via `du`; falls back to NA if a dir is missing or du fails)
size_report <- function(dirs = c("output", "p_output", "sim_output",
                                 "archive", "0-CuratingPlantData",
                                 "1-FilteringPlantData", "2-ModellingAnimals",
                                 "3-ExtPredictor", "4-ExtModeller",
                                 "5-SensitivityTesting")) {
  dirs <- dirs[fs::dir_exists(dirs)]
  bytes <- vapply(dirs, function(d) {
    out <- suppressWarnings(system2("du", c("-sb", "--", shQuote(d)),
                                    stdout = TRUE, stderr = FALSE))
    if (length(out)) as.numeric(sub("\t.*", "", out[1])) else NA_real_
  }, numeric(1))
  out <- data.frame(dir = dirs, size_GB = round(bytes / 1e9, 2))
  out[order(-out$size_GB), ]
}

## ---- Workflow -------------------------------------------------------------
## Run each step manually. Destructive steps default to dry_run = TRUE and
## only print what they would do. Current (non-archived) outputs are never touched.

## Step 1: Stale root files (.deb installer, .RDataTmp*) ----------------------
## Deletes only files not modified in the last `min_age_days` days.
clean_root_files <- function(min_age_days = 7, dry_run = TRUE) {
  f <- fs::dir_info(".", all = TRUE, regexp = "(\\.RDataTmp[0-9]*|\\.deb)$")
  f$age_days <- as.numeric(difftime(Sys.time(), f$modification_time, units = "days"))
  f$action <- ifelse(f$age_days >= min_age_days, "delete", "keep (recent)")
  print(data.frame(file = f$path, size = f$size, age_days = round(f$age_days),
                   action = f$action))
  if (!dry_run) fs::file_delete(f$path[f$action == "delete"])
  invisible(f)
}

## Step 2: Thin archived large files ------------------------------------------
## Within each taxon's animal_model/archive/, group large files by their name
## minus the timestamp (e.g. "simulatedOccurencesForAnimal.csv"), order by the
## timestamp in the filename, keep oldest / middle / newest, delete the rest.
plan_archive_thinning <- function(root = "output", min_size = "100MB") {
  arch <- fs::dir_ls(root, recurse = 2, type = "directory", regexp = "animal_model/archive$")
  f <- fs::dir_info(arch, recurse = TRUE, type = "file")
  f <- f[f$size >= fs::fs_bytes(min_size), ]
  f$taxon <- sub(paste0("^", root, "/([^/]+)/.*"), "\\1", f$path)
  f$stamp <- sub(".*([0-9]{8}_[0-9]{6}).*", "\\1", fs::path_file(f$path))
  f$group <- gsub("[._]?[0-9]{8}_[0-9]{6}", "", fs::path_file(f$path))
  f <- f[order(f$taxon, f$group, f$stamp), ]
  f$action <- "delete"
  for (g in split(seq_len(nrow(f)), paste(f$taxon, f$group))) {
    keep <- unique(g[c(1, ceiling((length(g) + 1) / 2), length(g))])
    f$action[keep] <- "keep"
  }
  data.frame(taxon = f$taxon, group = f$group, stamp = f$stamp,
             size = f$size, action = f$action, path = f$path)
}

apply_archive_thinning <- function(plan, dry_run = TRUE) {
  del <- plan$path[plan$action == "delete"]
  message(length(del), " files, ", format(sum(plan$size[plan$action == "delete"])),
          " to delete", if (dry_run) " (dry run)")
  if (!dry_run) fs::file_delete(del)
  invisible(del)
}

## Step 3: Compress archived runs ---------------------------------------------
## Each <root>/<taxon or sim>/animal_model/archive/<run>/ folder -> <run>.tar.gz,
## original folder removed only after the tarball is verified readable.
## Covers the real-taxon runs (output/) and the simulation runs (p_output/,
## sim_output/), which share the same layout.
compress_archive_runs <- function(root = c("output", "p_output", "sim_output"),
                                  dry_run = TRUE) {
  root <- root[fs::dir_exists(root)]
  arch <- fs::dir_ls(root, recurse = 2, type = "directory", regexp = "animal_model/archive$")
  if (!length(arch)) return(invisible(character()))
  runs <- fs::dir_ls(arch, type = "directory")
  for (r in runs) {
    tgz <- paste0(r, ".tar.gz")
    message(r, " -> ", tgz, if (dry_run) " (dry run)")
    if (dry_run) next
    status <- system2("tar", c("-czf", shQuote(tgz), "-C", shQuote(fs::path_dir(r)),
                               shQuote(fs::path_file(r))))
    ok <- status == 0 && system2("tar", c("-tzf", shQuote(tgz)), stdout = FALSE) == 0
    if (ok) fs::dir_delete(r) else warning("Compression failed, kept: ", r)
  }
  invisible(runs)
}

## Step 4: Compress superseded input files ------------------------------------
## For every <base>_cleaned.csv / <base>_cleaned_taxaFixed.csv, the other
## versions of <base> in the same folder (e.g. the raw GBIF <base>.csv) are
## added to a zip. If a zip containing <base> already exists in that folder
## (e.g. the original GBIF download), files go into it; otherwise <base>.zip
## is created. Originals are removed only after they're verified in the zip.
## compress_cleaned = TRUE also zips <base>_cleaned.csv when a
## _cleaned_taxaFixed.csv version exists.
plan_input_compression <- function(root = "2-ModellingAnimals/data",
                                   compress_cleaned = FALSE) {
  cleaned <- fs::dir_ls(root, recurse = TRUE, type = "file",
                        regexp = "_cleaned(_taxafixed)?\\.csv$", ignore.case = TRUE)
  bases <- unique(data.frame(
    dir  = fs::path_dir(cleaned),
    base = sub("_cleaned(_taxafixed)?\\.csv$", "", fs::path_file(cleaned), ignore.case = TRUE)
  ))
  rows <- lapply(seq_len(nrow(bases)), function(i) {
    d <- bases$dir[i]; b <- bases$base[i]
    sib <- fs::dir_ls(d, type = "file")
    nm <- fs::path_file(sib)
    is_ver <- startsWith(nm, b) & grepl("\\.(csv|tsv|txt)$", nm, ignore.case = TRUE)
    is_clean <- grepl("_cleaned", nm, ignore.case = TRUE)
    has_fixed <- any(grepl(paste0("^", b, "_cleaned_taxafixed\\.csv$"), nm, ignore.case = TRUE))
    is_plain_clean <- grepl(paste0("^", b, "_cleaned\\.csv$"), nm, ignore.case = TRUE)
    targ <- sib[is_ver & (!is_clean | (compress_cleaned & has_fixed & is_plain_clean))]
    if (!length(targ)) return(NULL)
    # Existing zip: name contains base, or its contents do
    zips <- sib[grepl("\\.zip$", nm, ignore.case = TRUE)]
    hit <- vapply(zips, function(z) {
      grepl(b, fs::path_file(z), fixed = TRUE) ||
        any(startsWith(utils::unzip(z, list = TRUE)$Name, b))
    }, logical(1))
    zipfile <- if (any(hit)) zips[hit][1] else fs::path(d, paste0(b, ".zip"))
    data.frame(file = targ, size = fs::file_size(targ),
               zip = as.character(zipfile),   # scalar: recycle across rows
               zip_exists = fs::file_exists(zipfile))
  })
  plan <- do.call(rbind, rows)
  # A file can be a "version" of two bases (e.g. <base>_noAntsNoBees.csv is
  # also the raw file of base <base>_noAntsNoBees); zip each file only once.
  plan[!duplicated(plan$file), ]
}

apply_input_compression <- function(plan, dry_run = TRUE) {
  message(nrow(plan), " files, ", format(sum(plan$size)), " to compress",
          if (dry_run) " (dry run)")
  if (dry_run) return(invisible(plan))
  for (z in unique(plan$zip)) {
    files <- plan$file[plan$zip == z]
    in_zip <- if (fs::file_exists(z)) utils::unzip(z, list = TRUE) else
      data.frame(Name = character(), Length = numeric())
    # Skip re-adding files already in the zip with identical size
    already <- fs::path_file(files) %in% in_zip$Name &
      as.numeric(fs::file_size(files)) ==
        in_zip$Length[match(fs::path_file(files), in_zip$Name)]
    to_add <- files[!already]
    if (length(to_add)) {
      status <- system2("zip", c("-j", "-9", shQuote(z), shQuote(to_add)))
      if (status != 0) { warning("zip failed for ", z); next }
    }
    final <- utils::unzip(z, list = TRUE)
    ok <- fs::path_file(files) %in% final$Name &
      as.numeric(fs::file_size(files)) == final$Length[match(fs::path_file(files), final$Name)]
    fs::file_delete(files[ok])
    if (any(!ok)) warning("Not verified in zip, kept: ", paste(files[!ok], collapse = ", "))
  }
  invisible(plan)
}

## Step 5: Exact duplicate files ----------------------------------------------
## Hash files >= min_size (only those sharing a byte size with another file),
## group identical ones. Review the plan and edit `action` before applying;
## by default the first path alphabetically in each group is kept. Deletions
## are replaced by hard links to the kept copy, so every path that scripts
## read (e.g. 1-SimTruthScenarios.R, 04_area_accounting.R) keeps working.
plan_duplicates <- function(dirs = c("0-CuratingPlantData", "1-FilteringPlantData",
                                     "2-ModellingAnimals", "3-ExtPredictor",
                                     "4-ExtModeller", "5-SensitivityTesting"),
                            min_size = "50MB") {
  f <- fs::dir_info(dirs, recurse = TRUE, type = "file")
  f <- f[as.numeric(f$size) >= as.numeric(fs::fs_bytes(min_size)), c("path", "size")]
  f$size <- as.numeric(f$size)  # drop fs_bytes class: %in% fails on it
  f <- f[f$size %in% f$size[duplicated(f$size)], ]
  if (!nrow(f)) return(NULL)
  f$md5 <- unname(tools::md5sum(f$path))
  f <- f[f$md5 %in% f$md5[duplicated(f$md5)], ]
  f <- f[order(f$md5, f$path), ]
  f$action <- ifelse(duplicated(f$md5), "delete", "keep")
  as.data.frame(f)
}

apply_duplicates <- function(plan, dry_run = TRUE) {
  # Never delete a file whose group has no kept copy
  safe <- plan$md5 %in% plan$md5[plan$action == "keep"]
  del <- plan$path[plan$action == "delete" & safe]
  message(length(del), " files, ", format(sum(plan$size[plan$path %in% del])),
          " to replace with hard links", if (dry_run) " (dry run)")
  if (dry_run) return(invisible(del))
  keep <- plan$path[plan$action == "keep"]
  for (p in del) {
    src <- keep[plan$md5[match(p, plan$path)] == plan$md5[plan$action == "keep"]][1]
    tmplink <- paste0(p, ".duplink")
    tryCatch(fs::link_create(src, tmplink, symbolic = FALSE),
             error = function(e) stop("Hard link failed for ", p, "; original kept."))
    fs::file_delete(p)
    fs::file_move(tmplink, p)
  }
  invisible(del)
}

## Step 6: gzip large raw inputs ----------------------------------------------
## gzip -k, verify with gzip -t, then remove the original.
## NOTE: classification.csv is read in
##   0-CuratingPlantData/scripts/curatingRawPlantData.R (line ~249, WFO.file =)
## Update that path to ".csv.gz" after compressing (fread reads .gz natively,
## requires R.utils).
gzip_inputs <- function(files = "0-CuratingPlantData/data/_DwC_backbone_R/classification.csv",
                        dry_run = TRUE) {
  files <- files[fs::file_exists(files)]
  print(data.frame(file = files, size = fs::file_size(files)))
  if (dry_run) { message("(dry run)"); return(invisible(files)) }
  for (f in files) {
    gz <- paste0(f, ".gz")
    ok <- system2("gzip", c("-k", "-9", shQuote(f))) == 0 &&
      system2("gzip", c("-t", shQuote(gz))) == 0
    if (ok) fs::file_delete(f) else warning("gzip failed, kept: ", f)
  }
  invisible(files)
}

## Step 7: Unzipped copies of files that are already in a zip -----------------
## For each zip, look for its entries extracted next to it, either loose in
## the zip's folder or in a subfolder named after the zip. An extracted file
## is deleted only if its name and byte size match the zip entry. Emptied
## folders are then removed.
plan_unzipped_copies <- function(dirs = c("0-CuratingPlantData", "2-ModellingAnimals",
                                          "3-ExtPredictor", "5-SensitivityTesting")) {
  zips <- fs::dir_ls(dirs, recurse = TRUE, type = "file", regexp = "\\.zip$",
                     ignore.case = TRUE)
  rows <- lapply(zips, function(z) {
    ent <- tryCatch(utils::unzip(z, list = TRUE),
                    error = function(e) {
                      message("Skipping unreadable zip: ", z)
                      NULL
                    })
    if (is.null(ent)) return(NULL)
    ent <- ent[ent$Length > 0, ]
    zd <- fs::path_dir(z)
    cand <- unique(c(fs::path(zd, ent$Name),
                     fs::path(zd, fs::path_ext_remove(fs::path_file(z)), ent$Name),
                     fs::path(zd, fs::path_ext_remove(fs::path_file(z)), fs::path_file(ent$Name))))
    cand <- cand[fs::file_exists(cand)]
    if (!length(cand)) return(NULL)
    zsize <- ent$Length[match(fs::path_file(cand), fs::path_file(ent$Name))]
    data.frame(zip = z, file = cand, size = as.numeric(fs::file_size(cand)),
               zip_size = zsize, match = as.numeric(fs::file_size(cand)) == zsize)
  })
  do.call(rbind, rows)
}

apply_unzipped_copies <- function(plan, dry_run = TRUE) {
  del <- plan$file[plan$match]
  message(length(del), " files, ", format(fs::fs_bytes(sum(plan$size[plan$match]))),
          " to delete", if (dry_run) " (dry run)")
  if (any(!plan$match)) message(sum(!plan$match), " name matches with differing size kept")
  if (dry_run) return(invisible(del))
  fs::file_delete(del)
  for (d in unique(fs::path_dir(del))) {
    if (length(fs::dir_ls(d, all = TRUE)) == 0) fs::dir_delete(d)
  }
  invisible(del)
}

## Future: convert large simulated range CSVs to a more efficient format -------
## The simulatedOccurencesForAnimal.*.csv files (~166 GB across output/, plus
## ~43 GB of sim CSVs in p_output/ and sim_output/) should eventually be
## converted to Parquet (arrow) or compressed CSV (.csv.gz / .csv.zst).
## Not done yet: the current files stay CSV until downstream readers are updated.

## ---- Run steps one at a time ----------------------------------------------
## size_report()
## clean_root_files()                    # review, then: clean_root_files(dry_run = FALSE)
## thin_plan <- plan_archive_thinning()  # review thin_plan
## apply_archive_thinning(thin_plan)     # then: apply_archive_thinning(thin_plan, dry_run = FALSE)
## compress_archive_runs()               # then: compress_archive_runs(dry_run = FALSE)
##                                       # (output/, p_output/ and sim_output/; or one
##                                       #  at a time, e.g. compress_archive_runs("p_output"))
## unz_plan <- plan_unzipped_copies()    # review unz_plan (run before step 4: removes
##                                       # files already in zips from the input plan)
## apply_unzipped_copies(unz_plan)       # then: apply_unzipped_copies(unz_plan, dry_run = FALSE)
## input_plan <- plan_input_compression() # review input_plan (regenerate after step 7)
## apply_input_compression(input_plan)   # then: apply_input_compression(input_plan, dry_run = FALSE)
## dup_plan <- plan_duplicates()         # review dup_plan; edit dup_plan$action to choose which copy stays
## apply_duplicates(dup_plan)            # then: apply_duplicates(dup_plan, dry_run = FALSE)
## gzip_inputs()                         # then: gzip_inputs(dry_run = FALSE); update curatingRawPlantData.R path
