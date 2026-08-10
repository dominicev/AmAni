# define a function to extract a random window, label patches, filter by minimum area,
# compute patch metrics, and write outputs to disk

# ============================================================================
# TEST VARIABLES - for line-by-line debugging (comment out when running function)
# ============================================================================
 # raster_path <- "def_00_20.tif" # update with actual path
 # seed <- 123
 # min_area_m2 <- 5e6
 # win_nrow <- 4000
 # win_ncol <- 4000
 # out_dir <- 
 # returnQ <- TRUE
 # apply_morph <- TRUE
 # morph_size <- 5 #5 is recommended
 # out_dir = getwd()
# ============================================================================

extract_patch_metrics <- function(
    raster_path,                 # SpatRaster: binary raster (1 = deforested, NA = background)
    seed = NULL,     # numeric or NULL: random seed for reproducibility (random method only)
    output_id = NULL,    # character or numeric: file naming prefix; defaults to seed if NULL
    min_area_m2 = 5e6,         # numeric: minimum patch area to retain (5e6 is 50 km^2)
    win_nrow = 4000,           # numeric: the height of the window
    win_ncol = 4000,           # numeric: the width of the window
    out_dir = getwd(),          # character: directory where outputs will be written
    returnQ = TRUE,
    apply_morph = FALSE,        # logical: apply morphological operations (dilation then erosion)
    morph_size = 5,             # numeric: size of the erosion kernel (dilation will be larger)
    small_patch_threshold = 50, # numeric: remove patches smaller than this many cells BEFORE morphology
    method = "random",          # character: "random" or "deterministic" sampling method
    start_row = NULL,           # numeric: starting row for deterministic method (required if method = "deterministic")
    start_col = NULL            # numeric: starting column for deterministic method (required if method = "deterministic")
) {
  # validate method parameter
  if (!method %in% c("random", "deterministic")) {
    stop("method must be either 'random' or 'deterministic'")
  }
  
  # validate deterministic method parameters
  if (method == "deterministic") {
    if (is.null(start_row) || is.null(start_col)) {
      stop("start_row and start_col must be provided when method = 'deterministic'")
    }
  }
  
  # Resolve output identifier for file naming (backward compatible)
  if (is.null(output_id)) output_id <- seed
  if (is.null(output_id)) stop("Either seed or output_id must be provided for file naming")
  
  # ensure output directory exists before attempting to write files
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  # optionally set random seed for reproducibility (only relevant for random method)
  # using NULL explicitly unsets any previously fixed RNG state
  if (method == "random") {
    set.seed(seed)
    print(seed)
  }

  
  ###importRasterObject
  # Only load raster metadata initially to minimize memory footprint
  raster_in<-rast(raster_path)
  
    # extract full raster dimensions
  nr <- nrow(raster_in)       # total number of rows in the raster
  nc <- ncol(raster_in)       # total number of columns in the raster
  
  {# sample a valid top-left corner for the window
  if (method == "random") {
    # randomly sample starting position
    start_row <- sample.int(nr - win_nrow + 1, 1)  # starting row index
    start_col <- sample.int(nc - win_ncol + 1, 1)  # starting column index
  } else {
    # use provided deterministic starting position
    # validate that the window fits within the raster bounds
    if (start_row < 1 || start_row > (nr - win_nrow + 1)) {
      stop(paste0("start_row must be between 1 and ", nr - win_nrow + 1))
    }
    if (start_col < 1 || start_col > (nc - win_ncol + 1)) {
      stop(paste0("start_col must be between 1 and ", nc - win_ncol + 1))
    }
  }
  
  # define row and column indices that make up the window
  rows <- start_row:(start_row + win_nrow - 1)   # row indices for subsetting
  cols <- start_col:(start_col + win_ncol - 1)   # column indices for subsetting
  
  # extract the raster window without dropping spatial metadata
  # Use crop with an extent to read only the needed portion from disk
  ext_to_read <- ext(raster_in)
  ext_to_read[1] <- xFromCol(raster_in, cols[1])      # xmin
  ext_to_read[2] <- xFromCol(raster_in, cols[length(cols)])  # xmax
  ext_to_read[3] <- yFromRow(raster_in, rows[length(rows)])  # ymin
  ext_to_read[4] <- yFromRow(raster_in, rows[1])      # ymax
  
  def_bin_sub <- crop(raster_in, ext_to_read)
  
  # Clear the full raster from memory immediately
  rm(raster_in)
  gc(verbose = FALSE)
  
  # attempt to trim NA-only outer margins to reduce computation
  # trim() fails if no NA border exists, so we guard with try()
  def_bin_sub_trimmed <- try(
    trim(def_bin_sub),        # remove only outer NA borders
    silent = TRUE
  )
  
  # if trimming succeeded, use the trimmed version; otherwise keep the cropped version
  if (!inherits(def_bin_sub_trimmed, "try-error")) {
    def_bin_sub <- def_bin_sub_trimmed
  }
  rm(def_bin_sub_trimmed)
  
  # coerce raster values to integer to reduce memory footprint
  def_bin_sub <- as.int(def_bin_sub)
  
  def_bin_original <- def_bin_sub
  
  # # store original for comparison plotting
  #   plot_title <- ifelse(method == "random", 
  #                        "randomly selected plot", 
  #                        "deterministically selected plot")
  #   plot(def_bin_original, main = plot_title, 
  #      col = c("white", "blue"), legend = FALSE)
    }
  
  #apply morphological operations to remove minuscule patches that will not factor into analysis
  #this will also hopefully speed up later calculations
  {
   
    # count_before <- sum(values(def_bin_sub) == 1, na.rm = TRUE)
    # 
    # cat("Cells before processing:", count_before, "\n")
    # 
    cat("Step 1: Severe thinning (erosion)...\n")
    # Create a large kernel for severe erosion
    kernel_erode <- matrix(1, nrow = morph_size+4, ncol = morph_size+4)
    
    # Create temporary file for intermediate result to avoid keeping in memory
    temp_eroded <- tempfile(fileext = ".tif")
    
    # Severe erosion: aggressively shrink patches
    # Write directly to disk to avoid memory accumulation
    def_bin_thinned <- focal(
      def_bin_sub,
      w = kernel_erode,
      fun = min,              # contract to minimum value in neighborhood
      na.rm = FALSE,
      NAonly = FALSE,
      filename = temp_eroded,
      overwrite = TRUE,
      wopt = list(datatype = "INT1U", gdal = c("COMPRESS=DEFLATE"))
    )
    
    # Force garbage collection and clear def_bin_sub from memory
    rm(def_bin_sub)
    gc(verbose = FALSE)
    
    cat("Step 2: Severe dilation...\n")
    # Create an even larger kernel for severe dilation
    kernel_dilate <- matrix(1, nrow = morph_size+8, ncol = morph_size+8)
    
    # Create temp file for dilation result
    temp_dilated <- tempfile(fileext = ".tif")
    
    # Severe dilation: expand patches back out
    # Write directly to disk
    def_bin_mask <- focal(
      def_bin_thinned,
      w = kernel_dilate,
      fun = max,              # expand to maximum value in neighborhood
      na.rm = TRUE,
      NAonly = FALSE,
      filename = temp_dilated,
      overwrite = TRUE,
      wopt = list(datatype = "INT1U", gdal = c("COMPRESS=DEFLATE"))
    )
    
    # Clear eroded raster from memory
    rm(def_bin_thinned)
    gc(verbose = FALSE)
    # Delete temp file
    unlink(temp_eroded)
    
    cat("Step 3: Applying mask to original...\n")
    # Use the mask to filter the original raster
    # Where mask = 1, keep original values; where mask = NA, set to NA
    def_bin_masked <- def_bin_original * def_bin_mask
    
    # Clear mask from memory
    rm(def_bin_mask)
    gc(verbose = FALSE)
    # Delete temp file
    unlink(temp_dilated)

    # # plot comparison: original, result, and difference
    # par(mfrow = c(1, 3), mar = c(2, 2, 3, 2))
    # plot(def_bin_original, main = "Original", 
    #      col = c("white", "darkred"), legend = FALSE)
    # plot(def_bin_masked, main = "Masked Original", 
    #      col = c("white", "purple"), legend = FALSE)
    # plot(def_bin_diff, main = "Removed by Mask", 
    #      col = c("white", "orange"), legend = FALSE)
    # par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1)  # reset to default
    # 
    # # Count cells after masking
    # count_after <- sum(values(def_bin_masked) == 1, na.rm = TRUE)
    # count_removed <- sum(values(def_bin_diff) == 1, na.rm = TRUE)
    # cat("Cells after masking:", count_after, "\n")
    # cat("Cells removed:", count_removed, "\n")
    def_bin_sub<-def_bin_masked
    }
  
  #optionally: use morphological dilation to expand remaining patches...
  # this will...
  #   a. connect adjacent patches that are effectively touching
  #   b. slightly exaggerate the size of patches (probably more prominent for smaller patches; yet, this may somewhat compensate for the patches we lost above with the filtering) 
  #   c. hopefully increase the computational efficiency of the code below because the total number of patches is lower.
  if (apply_morph){
    
    kernel_dilate <- matrix(1, nrow = morph_size, ncol = morph_size)
    
    # Create temp file for inflation result
    temp_inflated <- tempfile(fileext = ".tif")
    
    # Dilate to connect adjacent patches
    # Write directly to disk to save memory
    def_bin_inflated<- focal(
       def_bin_masked,
      w = kernel_dilate,
      fun = max,              # expand to max value in neighborhood
      na.rm = TRUE,
      NAonly = TRUE,
      filename = temp_inflated,
      overwrite = TRUE,
      wopt = list(datatype = "INT1U", gdal = c("COMPRESS=DEFLATE"))
    )
    
    # Clear masked raster from memory
    rm(def_bin_masked)
    gc(verbose = FALSE)
    
    # plot comparison: original, result, and difference
    # par(mfrow = c(2, 2), mar = c(2, 2, 3, 2))
    # plot(def_bin_original, main = "Original",
    #      col = c("white", "darkred"), legend = FALSE)
    # plot(def_bin_masked, main = "Masked Original",
    #      col = c("white", "purple"), legend = FALSE)
    # plot(def_bin_inflated, main = "inflated",
    #      col = c("white", "black"), legend = FALSE)
    # plot(def_bin_thinned, main = "thinned",
    #      col = c("white", "black"), legend = FALSE)
    # par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1)  # reset to default

    def_bin_sub<-def_bin_inflated
    
    # Clean up temp file (will be cleaned automatically on exit, but being explicit)
    # Note: we keep the raster object in memory, just not the temp file
  }
  
  
  # Check if there's any deforestation left after morphological operations
  n_deforest_cells <- sum(values(def_bin_sub) == 1, na.rm = TRUE)
  
  if (n_deforest_cells == 0) {
    cat("No deforestation cells remaining after morphological operations\n")
    done_flag_file <- file.path(out_dir, paste0("w", as.character(output_id), ".done"))
    file.create(done_flag_file)
    return(invisible(NULL))
  }
  
  cat("Found", n_deforest_cells, "deforestation cells\n")
  
  # label connected deforestation patches using 8-direction (queen) connectivity
  cat("Labeling final patches...\n")
  patches_sub <- patches(
    def_bin_sub,              # binary raster subset
    directions = 8            # queen adjacency
  )
  
  # Clear def_bin_sub from memory after patches are created
  rm(def_bin_sub)
  gc(verbose = FALSE)
  
  # Check if patches() returned a valid result
  if (is.null(patches_sub) || ncell(patches_sub) == 0) {
    cat("patches() returned NULL or empty result\n")
    done_flag_file <- file.path(out_dir, paste0("w", as.character(output_id), ".done"))
    file.create(done_flag_file)
    return(invisible(NULL))
  }
  
  # reproject patches to an equal-area CRS in meters for correct area calculations
  patches_m <- project(
    patches_sub,              # labeled patch raster in lon/lat
    "EPSG:5880",              # equal-area CRS for South America
    method = "near"           # nearest-neighbor to preserve patch IDs
  )
  
  # compute per-cell area (in square meters) in the projected raster
  cell_area_rast <- cellSize(
    patches_m,                # projected raster
    unit = "m"                # square meters
  )
  
  # collapse per-cell areas to a single mean value (cells are uniform after projection)
  cell_area_m2 <- mean(
    values(cell_area_rast),   # numeric vector of cell areas
    na.rm = TRUE              # ignore NA cells
  )
  
  # extract patch IDs for all raster cells as a vector
  vals <- values(patches_m)   # patch ID per cell, NA for background
  
  # remove NA values so we only count real patch cells
  vals_no_na <- vals[!is.na(vals)]
  
  # count how many cells belong to each patch ID
  cell_counts <- table(vals_no_na)  # named vector: names = patch IDs
  
  # convert cell counts to patch areas in square meters
  patch_area_m2 <- as.numeric(cell_counts) * cell_area_m2
  
  # identify patch IDs that meet the minimum area threshold
  keep_ids <- as.integer(
    names(cell_counts)[patch_area_m2 >= min_area_m2]
  )
  
  
  # handle the edge case where no patches meet the size threshold
  if (length(keep_ids) == 0) {
    
    # emit a warning so this condition is explicit
    warning("No patches meet the minimum area threshold")
    done_flag_file <- file.path(out_dir, paste0("w", as.character(output_id), ".done"))
    file.create(done_flag_file)
    
    # return an empty result invisibly
    return(invisible(NULL))
  }
  
  # create a filtered raster that removes all small patches
  patches_filt <- patches_sub                  # copy labeled raster
  patches_filt[!(patches_filt %in% keep_ids)] <- NA  # drop small patches
  
  # convert filtered patches to polygons, dissolving by patch ID
  patch_poly <- as.polygons(
    patches_filt,              # raster of retained patch IDs
    dissolve = TRUE,           # one polygon per patch
    values = TRUE              # keep patch ID attribute
  )
  
  # compute patch areas from polygons (in km² for reporting)
  patch_area_km2 <- expanse(
    patch_poly,                # polygon patches
    unit = "km"                # square kilometers
  )
  
  # compute patch perimeters from vector geometry (in km)
  patch_perim_km <- perim(
    patch_poly                # polygon patches
  )
  
  # compute patch centroids
  patch_centroids <- centroids(
    patch_poly                # polygon patches
  )
  
  # extract centroid coordinates into a matrix
  centroid_coords <- crds(
    patch_centroids            # centroid geometries
  )
  
  # assemble final patch metrics table
  patch_metrics <- data.frame(
    patch_id     = keep_ids,                                  # patch ID
    area_km2     = patch_area_km2,                            # area in km²
    perimeter_km = patch_perim_km,                            # perimeter in km
    centroid_x   = centroid_coords[, 1],                      # X coordinate
    centroid_y   = centroid_coords[, 2],                      # Y coordinate
    compactness  = (4 * pi * patch_area_km2) /                # shape metric
      (patch_perim_km^2)
  )
  
  # drop incomplete rows defensively
  patch_metrics <- patch_metrics[complete.cases(patch_metrics), ]
  
  # write patch metrics table to CSV
  write.csv(
    patch_metrics,
    file = file.path(out_dir, paste0("w", as.character(output_id),".patch_metrics.csv")),
    row.names = FALSE
  )

  
  # write filtered patch raster to disk for reproducibility
  writeRaster(
    patches_filt,
    filename = file.path(out_dir, paste0("w", as.character(output_id),".patches_filtered.tif")),
    overwrite = TRUE
  )
  
  # write patch polygons to disk as a GeoPackage
  writeVector(
    patch_poly,
    filename = file.path(out_dir, paste0("w", as.character(output_id),".patch_polygons.gpkg")),
    overwrite = TRUE
  )
  
  # write a small completion flag to signal a fully successful run
  # this prevents partially written outputs from being treated as valid
  done_flag_file <- file.path(out_dir, paste0("w", as.character(output_id), ".done"))
  file.create(done_flag_file)
  
  # return results invisibly for interactive use
  if (returnQ) {
    invisible(
      list(
        metrics  = patch_metrics,
        raster   = patches_filt,
        polygons = patch_poly
      )
    )
  } else {
    invisible(NULL)
  }
}
