######Make function for sampling

#function for sampling a proportion of species in an occurrence dataset
#prop_species : a number between 0 and 1. the proportion of species to sample

sample_Species <-function(true_animal_occ_data, prop_species = 0.8){
  dt <- as.data.table(true_animal_occ_data)
  # Get unique species
  all_species <- unique(dt$species)
  # Sample species
  sampled_species_names <- sample(
    all_species, 
    size = round(length(all_species) *prop_species )
  )
  
  setkey(dt, species)
  speciesSample <- dt[.(sampled_species_names)]
  
  # # 1. Get a random sample of unique species names
  # sampled_species_names <- true_animal_occ_data %>%
  #   distinct(species) %>%
  #   slice_sample(prop = prop_species) %>%
  #   pull(species) # converts the column into a vector
  # 
  # #Filter the original dataset to only include those species
  # speciesSample <- true_animal_occ_data %>%
  #   filter(species %in% sampled_species_names)
  # 
  # # Self-check: This will now equal exactly 0.8
  # length(unique(speciesSample$species)) / length(unique(true_animal_occ_data$species))
  
  return(speciesSample)
}

sample_Species_perSite <- function(siteSample, prop_species = 0.8, min_prop_species = prop_species) {
  dt <- as.data.table(siteSample)
  
  # Pre-generate one random proportion per grid cell so it can be captured
  grid_ids <- unique(dt$grid_id)
  prop_draws <- data.table(
    grid_id  = grid_ids,
    prop_drawn = runif(length(grid_ids), min = min_prop_species, max = prop_species)
  )
  
  # Sample a randomized proportion of species per grid cell (each cell gets a different proportion)
  siteSpeciesSample <- dt[prop_draws, on = "grid_id"][
    , .SD[species %in% sample(unique(species), max(1, round(length(unique(species)) * prop_drawn[1])))],
    by = grid_id
  ]
  
  return(list(
    sample     = siteSpeciesSample,
    prop_draws = prop_draws
  ))
}


#function for sampling a proportion of sites in an occurence dataset
#
#
#
#
#richness_sites: Options: low, medium, high; Default: NULL
#endemicity_sites: Options: low, medium, high; Default: NULL
##clustering: Options: low, medium, high, GBIF; Default: NULL
##sampling_evenness: Options: even (Default),veryLow, low
sample_Sites <-function(true_animal_occ_data, prop_species = 1,  prop_sites = 0.8, n_sites = NULL, endemicity_sites = NULL, richness_sites = NULL, clustering = NULL, animal_occ_Samples, sampling_evenness = "even"){
  
  message("sample_Sites: starting")
  message("  params: prop_species=", prop_species, ", prop_sites=", prop_sites,
          ", n_sites=", ifelse(is.null(n_sites), "NULL", n_sites),
          ", endemicity=", ifelse(is.null(endemicity_sites), "NULL", endemicity_sites),
          ", richness=", ifelse(is.null(richness_sites), "NULL", richness_sites),
          ", clustering=", ifelse(is.null(clustering), "NULL", clustering),
          ", evenness=", sampling_evenness)
  
  occ_data<-true_animal_occ_data
  message("  input data: ", nrow(occ_data), " rows, ",
          length(unique(occ_data$grid_id)), " unique grid cells, ",
          length(unique(occ_data$species)), " unique species")
  
  # --- parameter metadata (populated throughout the function) ---
  # param_prop_species       <- prop_species
  # param_prop_sites         <- prop_sites
  param_n_sites_requested  <- if (is.null(n_sites)) NA_integer_ else as.integer(n_sites)
  param_n_sites_actual     <- NA_integer_
  param_endemicity_label   <- if (is.null(endemicity_sites)) NA_character_ else endemicity_sites
  param_endemism_min       <- NA_real_
  param_endemism_max       <- NA_real_
  param_richness_label     <- if (is.null(richness_sites)) NA_character_ else richness_sites
  param_richness_min       <- NA_real_
  param_richness_max       <- NA_real_
  param_clustering_label   <- if (is.null(clustering)) NA_character_ else clustering
  param_diffusion          <- NA_real_
  param_numberStartSites   <- NA_integer_
  param_min_weight         <- NA_real_
  param_sampling_evenness  <- sampling_evenness
  param_min_prop_species   <- NA_real_
  
  #Calculate endemicity and richness of grid cells if specified

    message("  ... calculating grid cell diversity stats")
    
    # Use data.table for speed: drop geometry column, compute in two passes
    occ_dt <- as.data.table(st_drop_geometry(true_animal_occ_data))
    occ_dt[, range_size   := .N, by = species]
    occ_dt[, species_weight := 1 / range_size]
    animal_weighted_endemism <- occ_dt[, .(
      weighted_endemism = sum(species_weight),
      species_richness  = .N
    ), by = grid_id]

  #sample sites based on their weighted endemism
  
  if(is.null(endemicity_sites)){    
    message("  ... skipping endemicity-based sampling")
    param_endemism_max <- max(animal_weighted_endemism$weighted_endemism)
    param_endemism_min <- min(animal_weighted_endemism$weighted_endemism)
    
  }else{
    endemicity_cells <- switch(endemicity_sites,
      "low" = {
        param_endemism_max <- quantile(animal_weighted_endemism$weighted_endemism, 0.25)
        animal_weighted_endemism[weighted_endemism < param_endemism_max, grid_id]
      },
      "medium" = {
        param_endemism_min <- quantile(animal_weighted_endemism$weighted_endemism, 0.375)
        param_endemism_max <- quantile(animal_weighted_endemism$weighted_endemism, 0.625)
        animal_weighted_endemism[weighted_endemism >= param_endemism_min & weighted_endemism <= param_endemism_max, grid_id]
      },
      "high" = {
        param_endemism_min <- quantile(animal_weighted_endemism$weighted_endemism, 0.75)
        animal_weighted_endemism[weighted_endemism > param_endemism_min, grid_id]
      }
    )
    # %chin% is data.table's fast character %in%
    occ_data <- occ_data[occ_data$grid_id %in% endemicity_cells, ]
    message("  endemicity filter ('", endemicity_sites, "'): ",
            length(endemicity_cells), " cells retained, ",
            nrow(occ_data), " rows remain")
  }
  
  
  #sample sites based on their richness
  if(is.null(richness_sites)){    
    message("  ... skipping richness-based sampling")
    param_richness_max <- max(animal_weighted_endemism$species_richness)
    param_richness_min <- min(animal_weighted_endemism$species_richness)
  }else{
    richness_cells <- switch(richness_sites,
      "low" = {
        param_richness_max <- quantile(animal_weighted_endemism$species_richness, 0.25)
        animal_weighted_endemism[species_richness < param_richness_max, grid_id]
      },
      "medium" = {
        param_richness_min <- quantile(animal_weighted_endemism$species_richness, 0.375)
        param_richness_max <- quantile(animal_weighted_endemism$species_richness, 0.625)
        animal_weighted_endemism[species_richness >= param_richness_min & species_richness <= param_richness_max, grid_id]
      },
      "high" = {
        param_richness_min <- quantile(animal_weighted_endemism$species_richness, 0.75)
        animal_weighted_endemism[species_richness > param_richness_min, grid_id]
      }
    )
    
    occ_data <- occ_data[occ_data$grid_id %in% richness_cells, ]
    message("  richness filter ('", richness_sites, "'): ",
            length(richness_cells), " cells retained, ",
            nrow(occ_data), " rows remain")
    
  }
  
  #if sampling is spatially clustered, make a weights vector
  if(is.null(clustering)){    
    message("  sites will be sampled randomly with respect to space")
    weights_full <- rep(1, length(unique(occ_data$grid_id)))
  }else{
    message("  applying spatial clustering weights (mode: '", clustering, "')")
    #1. use GBIF data to get base weights for the current version of occ_data
    valid_grid_ids <- occ_data$grid_id
    filtered_occ_Samples <- animal_occ_Samples %>% 
      filter(grid_id %in% valid_grid_ids)
    
    
    baseProbabilities <- as.data.table(filtered_occ_Samples)[, .(n = .N), by = grid_id]
    message("  GBIF base probabilities cover ", nrow(baseProbabilities), " of ",
            length(unique(occ_data$grid_id)), " grid cells")
    
    # Compute centroids once on unique grid cells only (occ_data has many rows per grid)
    unique_grids_sf <- occ_data[!duplicated(occ_data$grid_id), ]
    all_coords      <- st_coordinates(st_centroid(unique_grids_sf$geometry))
    unique_grid_ids <- unique_grids_sf$grid_id
    

    #diffusion value notes:
    # 500000 covers about 1/4 the amazon
    # 250000 covers about 1/20th the amazon
    # 100000 covers about 1/100th the amazon
    
    switch(clustering, 
           "low" = {
             cluster_weights_dt <- clusterWeights(baseProbabilities,diffusion=250000, numberStartSites=50, all_coords = all_coords, unique_grid_ids=unique_grid_ids)
             cluster_weights_dt <- cluster_weights_dt[match(baseProbabilities$grid_id, cluster_weights_dt$grid_id)]
             weights <- (cluster_weights_dt$cluster_weight * baseProbabilities$n)+rep(0.1, length(baseProbabilities$n))
             param_diffusion        <- 250000
             param_numberStartSites <- 50
             param_min_weight       <- 0.1
           },
           "medium" = {
             cluster_weights_dt <- clusterWeights(baseProbabilities, diffusion=200000, numberStartSites=2, all_coords = all_coords, unique_grid_ids=unique_grid_ids)
             cluster_weights_dt <- cluster_weights_dt[match(baseProbabilities$grid_id, cluster_weights_dt$grid_id)]
             weights <- (cluster_weights_dt$cluster_weight * baseProbabilities$n)+rep(0.05, length(baseProbabilities$n))
             param_diffusion        <- 200000
             param_numberStartSites <- 2
             param_min_weight       <- 0.05
           },
           "high" = {
             cluster_weights_dt <- clusterWeights(baseProbabilities, diffusion=100000, numberStartSites=1, all_coords = all_coords, unique_grid_ids=unique_grid_ids)
             cluster_weights_dt <- cluster_weights_dt[match(baseProbabilities$grid_id, cluster_weights_dt$grid_id)]
             weights <- (cluster_weights_dt$cluster_weight * baseProbabilities$n)+rep(0.01, length(baseProbabilities$n))
             param_diffusion        <- 100000
             param_numberStartSites <- 1
             param_min_weight       <- 0.01
             # 
             # grid_geoms <- occ_data %>% 
             #   select(grid_id) %>% 
             #   distinct(grid_id)
             # plot_df_flat <- occ_data %>%
             #   st_drop_geometry() %>%                  # Drops the heavy spatial attributes
             #   distinct(grid_id, .keep_all = TRUE) %>% # Lightning fast on a regular dataframe
             #   left_join(cluster_weights_dt, by = "grid_id")
             # plot_df <- grid_geoms %>% 
             #   left_join(plot_df_flat, by = "grid_id")
             # 
             # ggplot() +
             #   geom_sf(data = plot_df, aes(fill = cluster_weight), color = NA) +
             #   #geom_sf(data = filter(plot_df, is_start), color = "red", linewidth = 1, fill = NA) +
             #   scale_fill_viridis_c(option = "magma", name = "weight", na.value = "transparent") +
             #   labs(
             #     title = paste0("Cluster sampling weights | diffusion = ", cluster_weights_dt$diffusion[1], " m, min_weight = ", cluster_weights_dt$min_weight[1]),
             #     x = NULL, y = NULL
             #   ) +
             #   theme_minimal()
             # 
             
             
             
             
           },
           "GBIF" = {
             weights <- baseProbabilities$n
             }
    )
    
    message("  clustering weights computed (range: ",
            round(min(weights), 4), " - ", round(max(weights), 4), ")")
    #clustering = "low"
    # Extend weights to cover all unique grid_ids in occ_data using data.table (fast)
    min_wt <- if(clustering == "GBIF") 1 else cluster_weights_dt$min_weight[1]
    
    # Create a data.table with unique grid_ids from occ_data
    occ_distinct <- data.table(grid_id = unique(occ_data$grid_id))
    
    # Create a data.table with the weights we calculated
    weights_dt <- data.table(grid_id = baseProbabilities$grid_id, weight = weights)
    
    # Right join on occ_distinct so every grid in occ_data gets a weight
    weights_full <- weights_dt[occ_distinct, on = "grid_id"][
      is.na(weight), weight := min_wt
    ][, weight]
    
    # n_filled is the count of grids in occ_data absent from baseProbabilities — no second join needed
    n_filled <- nrow(occ_distinct) - nrow(baseProbabilities)
    if(n_filled > 0) message("  ", n_filled, " grid cells had no GBIF data; assigned min_weight=", min_wt)
    
    # # # Plot cluster weights for diagnostics
    # # # Note: all_coords and unique_grid_ids are now per unique grid cell (not per row)
    # startingCells <- attr(cluster_weights_dt, "startingCells")
    # plot_df <- as.data.frame(all_coords) |>
    #   mutate(
    #     grid_id = unique_grid_ids,
    #     cluster_weight = weights_full,
    #     is_start = grid_id %in% startingCells
    #   ) |>
    #   distinct(grid_id, .keep_all = TRUE) |>
    #   st_as_sf(coords = c("X", "Y"), crs = st_crs(occ_data)) |>
    #   st_buffer(dist = 31180, nQuadSegs = 1)
    # 
    # ggplot() +
    #   geom_sf(data = plot_df, aes(fill = cluster_weight), color = NA) +
    #   geom_sf(data = filter(plot_df, is_start), color = "red", linewidth = 1, fill = NA) +
    #   scale_fill_viridis_c(option = "magma", name = "weight") +
    #   labs(
    #     title = paste0("Cluster sampling weights  |  diffusion = ", cluster_weights_dt$diffusion[1], " m,  min_weight = ", cluster_weights_dt$min_weight[1]),
    #     x = NULL, y = NULL
    #   ) +
    #   theme_minimal()
    
    
  }
  
  # Create a DT with unique grid_ids and their corresponding weights
  grid_weights <- data.table(
    grid_id = unique(occ_data$grid_id),
    weight = weights_full
  )
  
  
  if(is.null(n_sites)){
    n_sample <- round(nrow(grid_weights) * prop_sites)
    message("  sampling ", n_sample, " of ", nrow(grid_weights), " grid cells (prop_sites=", prop_sites, ")")
    sampled_grid_ids <- grid_weights[, sample(grid_id, size = n_sample, prob = weight)]
  } else {
    message("  sampling ", n_sites, " of ", nrow(grid_weights), " grid cells (n_sites)")
    sampled_grid_ids <- grid_weights[, sample(grid_id, size = n_sites, prob = weight)]
  }
  param_n_sites_actual <- length(sampled_grid_ids)
  # if(param_n_sites_actual==n_sites){ message("n_sites is equal to the number of grid_ids being sampling....it's working")} else{ message(paste0("PROBLEM: n_sites = ", toString(n_sites), " but length(sampled_grid_ids) =",toString(param_n_sites_actual), "...this isn't right...fix it" ))
  # }
  
  
  str(sampled_grid_ids)

  
  # %chin% is data.table's fast character %in% for sf row subsetting
  siteSample <- occ_data[occ_data$grid_id %in% sampled_grid_ids, ]
  message("  site sample: ", nrow(siteSample), " rows across ",
          length(sampled_grid_ids), " grid cells")
  
  message("  applying species sampling (evenness='", sampling_evenness, "')")
  prop_draws_df <- NULL
  switch(sampling_evenness,
         
         "even" = {
           siteSpeciesSample <- sample_Species(siteSample, prop_species = prop_species)
           param_min_prop_species <- prop_species
           
         },
         "low" = {
           res <- sample_Species_perSite(siteSample, prop_species = prop_species, min_prop_species = prop_species / 2)
           siteSpeciesSample      <- res$sample
           prop_draws_df          <- res$prop_draws
           param_min_prop_species <- prop_species / 2
         },
         "veryLow" = {
           res <- sample_Species_perSite(siteSample, prop_species = prop_species, min_prop_species = prop_species / 10)
           siteSpeciesSample      <- res$sample
           prop_draws_df          <- res$prop_draws
           param_min_prop_species <- prop_species / 10
         }
         
  )
  
  
  message("  final sample: ", nrow(siteSpeciesSample), " rows, ",
          length(unique(siteSpeciesSample$species)), " species, ",
          length(unique(siteSpeciesSample$grid_id)), " grid cells")
  message("sample_Sites: done")
  
  # Create output as a list with occurrence data nested and prop_draws if applicable
  result <- list(occurenceData = siteSpeciesSample[, c("grid_id", "geometry", "species")])
  
  # Add prop_draws only if it exists (from "low" or "veryLow" evenness)
  if (!is.null(prop_draws_df)) {
    result$prop_draws <- prop_draws_df
  }

  # Attach all run parameters as attributes on the list
  attr(result, "prop_species")       <- prop_species
  attr(result, "prop_sites")         <- prop_sites
  attr(result, "n_sites_requested")  <- n_sites
  attr(result, "n_sites_actual")     <- param_n_sites_actual
  attr(result, "endemicity_label")   <- param_endemicity_label
  attr(result, "endemism_min")       <- param_endemism_min
  attr(result, "endemism_max")       <- param_endemism_max
  attr(result, "richness_label")     <- param_richness_label
  attr(result, "richness_min")       <- param_richness_min
  attr(result, "richness_max")       <- param_richness_max
  attr(result, "clustering_label")   <- param_clustering_label
  attr(result, "diffusion")          <- param_diffusion
  attr(result, "numberStartSites")   <- param_numberStartSites
  attr(result, "min_weight")         <- param_min_weight
  attr(result, "sampling_evenness")  <- param_sampling_evenness
  attr(result, "min_prop_species")   <- param_min_prop_species
  
  return(result)
}

clusterWeights<-function(baseProbabilities, diffusion=50000, numberStartSites=1, all_coords = c(), unique_grid_ids=c()){
  min_weight=1/max(baseProbabilities)
  #2. choose some starting cells based on the base weights (which are calculated based on the number of empirical observations per cell)
  startingCells <- sample(baseProbabilities$grid_id, prob = baseProbabilities$n, size = numberStartSites)
  
  # Starting cell coordinates (unique_grid_ids already has one row per grid)
  start_coords <- all_coords[unique_grid_ids %in% startingCells, , drop = FALSE]
  
  # Vectorized Euclidean distance: n_unique_grids x n_start_sites
  dx <- outer(all_coords[, 1], start_coords[, 1], "-")
  dy <- outer(all_coords[, 2], start_coords[, 2], "-")
  dist_matrix <- sqrt(dx^2 + dy^2)
  
  # pmin across columns is faster than apply(..., 1, min) for multiple start sites;
  # for a single start site it reduces to a plain vector
  dist_to_nearest <- if (ncol(dist_matrix) == 1L) {
    dist_matrix[, 1]
  } else {
    do.call(pmin, as.data.frame(dist_matrix))
  }
  
  # Exponential decay: weight = min_weight + (1 - min_weight) * exp(-d / diffusion)
  cluster_weights <- min_weight + (1 - min_weight) * exp(-dist_to_nearest / diffusion)
  
  # One row per unique grid cell
  result <- data.table(
    grid_id          = unique_grid_ids,
    cluster_weight   = cluster_weights,
    min_weight       = min_weight,
    diffusion        = diffusion,
    numberStartSites = numberStartSites
  )
  
  # Store startingCells as an attribute instead of using <<-
  attr(result, "startingCells") <- startingCells
  
  return(result)
}


