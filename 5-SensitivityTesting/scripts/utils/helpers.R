###for calculating raw endemism
rename_grid_id <- function(df) {
  if ("grid_id.x" %in% colnames(df)) {
    df <- dplyr::rename(df, grid_id = grid_id.x)
  }
  df
}



#the below function seems to be faster than the two above
#

rebuild_community <- function(animal_grid, sampled_df, scenario_name = NULL) {
  
  # Note: animal_grid must already be a data.table object! 
  # If not, convert it outside this function using setDT(animal_grid)
  
  sp <- unique(sampled_df$species)
  # %chin% is a fast version of %in% for character vectors
  # unique(..., by = ...) is highly optimized in data.table
  out <- unique(animal_grid[species %chin% sp], by = c("species", "grid_id"))

  if (!is.null(scenario_name)) {
    # Modify by reference (zero memory copy)
    out[, scenario := scenario_name]
  }
  
  return(out)
}

rebuild_community2 <- function(animal_grid, sampled_df, scenario_name = NULL) {
  
  # Note: animal_grid must already be a data.table object! 
  # If not, convert it outside this function using setDT(animal_grid)
  
  sp <- unique(sampled_df$species)
  ids <- unique(sampled_df$grid_id)
  # %chin% is a fast version of %in% for character vectors
  # unique(..., by = ...) is highly optimized in data.table
  out <- unique(animal_grid[species %chin% sp], by = c("species", "grid_id"))
  out2 <- unique(out[grid_id %in% ids], by = c("species", "grid_id"))
  
  if (!is.null(scenario_name)) {
    # Modify by reference (zero memory copy)
    out2[, scenario := scenario_name]
  }
  
  return(out2)
}



#the weights should add up to 1. They probably don't NEED to, but it makes sense to limit them to 1 so the weights across the different schemes are more easily compared

biased_species_sample <- function(data,
                                  n_quantiles = 5,
                                  weights,
                                  sample_frac = 0.8,
                                  seed = 1,
                                  smooth = TRUE) {
  set.seed(seed)
  
  df <- data
  
  # 1. assign quantile classes (kept for diagnostics regardless of smooth mode)
  df <- df %>%
    mutate(quantile_class = ntile(n_cells, n_quantiles))
  
  # 2. ensure weights are properly named
  if (is.numeric(weights)) {
    if (is.null(names(weights))) {
      names(weights) <- as.character(seq_along(weights))
    }
  }
  
  if (smooth) {
    # --- Smooth weighting: interpolate continuously across percentile rank ---
    
    # Percentile rank of each species by n_cells (0–1 scale)
    df <- df %>%
      mutate(pct_rank = rank(n_cells, ties.method = "average") / n())
    
    # Control points at the midpoint of each quantile class
    midpoints <- (seq_len(n_quantiles) - 0.5) / n_quantiles
    
    # Linearly interpolate weights; rule = 2 clamps at boundaries
    df <- df %>%
      mutate(weight = approx(midpoints, as.numeric(weights), 
                             xout = pct_rank, rule = 2)$y)
    
  } else {
    # --- Original discrete weighting ---
    df <- df %>%
      mutate(
        quantile_class_chr = as.character(quantile_class),
        weight = weights[quantile_class_chr]
      )
  }
  
  # 3. safety check: missing or non-positive weights
  if (any(is.na(df$weight))) {
    stop("Some species received NA weights.")
  }
  
  # 4. normalize probabilities
  df <- df %>%
    mutate(prob = weight / sum(weight))
  
  # 5. sample at a user-defined fraction (Default: 80% of species)
  n_target <- floor(sample_frac * nrow(df))
  
  sampled <- df %>%
    slice_sample(n = n_target, weight_by = prob, replace = FALSE)
  
  # 6. diagnostics: quantile composition shift (discrete bins, for interpretability)
  quantile_shift <- rbind(
    full = prop.table(table(df$quantile_class)),
    sample = prop.table(table(sampled$quantile_class))
  )
  
  list(
    sampled_data = sampled,
    full_data = df,
    quantile_shift = quantile_shift,
    weights_used = weights,
    smooth = smooth
  )
}


species_richness_sample <- function(data,
                                    richness = 10000,
                                    seed = 1) {
  set.seed(seed)
  
  df <- data
  weights<-make_weights()
  # 1. assign quantiles using ntile (robust to ties)
  df <- df %>%
    mutate(quantile_class = ntile(n_cells, 5))
  
  # 2. ensure weights are properly named
  if (is.numeric(weights)) {
    if (is.null(names(weights))) {
      names(weights) <- as.character(seq_along(weights))
    }
  }
  
  # 3. map weights to quantile classes safely
  df <- df %>%
    mutate(
      quantile_class_chr = as.character(quantile_class),
      weight = weights[quantile_class_chr]
    )
  
  # 4. safety check: missing weights
  if (any(is.na(df$weight))) {
    stop("Some quantile classes do not have matching weights.")
  }
  
  # 5. normalize probabilities
  df <- df %>%
    mutate(prob = weight / sum(weight))
  
  # 6. sample at a user-defined number of species
  n_target <- richness
  
  sampled <- df %>%
    slice_sample(n = n_target, weight_by = prob, replace = TRUE)#sampling with replacement so we can get hyperdiverse taxa with more than the source population
  
  # 7. diagnostics: quantile composition shift
  quantile_shift <- rbind(
    full = prop.table(table(df$quantile_class)),
    sample = prop.table(table(sampled$quantile_class))
  )
  
  list(
    sampled_data = sampled,
    full_data = df,
    quantile_shift = quantile_shift,
    weights_used = weights,
    richness = richness
  )
}


#with this function, you input relative weights and it rescales the base weights appropriately. the relative weights should be addative differences between -1 and +1. 0 would be no change, positive values increase the weights and negative values decrease them. The input vector should sum to 0.
#NOTE: The input vector has 4 items but the output returns 5. This is because there are 5 quantile classes but the first two have the same endemism

make_weights<-function(weightKernel=c(
  0,#high endemism weights
  0, #mid endemism weights
  0, #low endemism weights
  0 #very low endemism weights
)){
  
  #baseWeights are the weights that would recreate the sampled data (raw animal occurences)
  #make
  baseWeights<-weights_df$sampling_weight
  
  #weightKernel is the multiplier to create an experimental rep that will increase or decrease certain endemism classes
  mult_weightKernel<-weightKernel+c(1, 1, 1, 1)
  
  #If the below is TRUE, the weights add up to zero
  (sum(weightKernel)+weightKernel[[1]])<0.001
  weights<-mult_weightKernel*baseWeights
  weights<-unlist(c(weights[[1]], weights))
  weights <- weights / sum(weights)
  return(weights)
}