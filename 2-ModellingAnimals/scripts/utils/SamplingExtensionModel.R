warning("WARNING: SamplingExtensionModel.R and extend_richness_from_sampling() will fail unless fit_params == c(\"amplitude\", \"meanlog\",\"sdlog\" )
")

p_model<-readRDS(file.path(project_root, "5-SensitivityTesting", "results", "p_model_final.rds"))

testCoefficients<-readRDS(file.path(project_root, "5-SensitivityTesting", "results", "p_model_final.rds"))

source(file.path(project_root, "2-ModellingAnimals", "scripts", "utils", "computeCellStats.R"))
# =============================================================================
# extend_richness_from_sampling_predict()
# -----------------------------------------------------------------------------
# Ported from 5-SensitivityTesting/scripts/7-PredictiveModel.R (Section 4e).
# Predicts the SIGNED proportional error of a richness estimate from observable
# sampling characteristics, then converts it to a completed ("full") gamma
# richness (predictors + partialGamma in, completed richness out). The signed
# error, signed_prop_error = (estimate - true) / true, distinguishes over- from
# under-estimates and inverts as fullGamma = partialGamma / (1 + signed error).
#
#   * Predictors are supplied in NATURAL units (not pre-scaled). The function
#     standardizes them internally using the training mean/sd from
#     p_model_scale_params.csv, so callers never z-score by hand.
#   * Returns completed richness (fullGamma) plus a lower/upper richness
#     interval derived from the model's prediction (or confidence) interval.
#   * Any predictor left NULL defaults to that predictor's training median for
#     the chosen fill_Plants level (mirrors the Section 4c worked examples). But
#     still, this shouldn't be run with anything NULL.
#
# Requires (loaded once below; paths are relative to the PROJECT ROOT, matching
# the p_model load at the top of this file):
#   5-SensitivityTesting/results/p_model_signed.rds
#   5-SensitivityTesting/results/p_model_scale_params.csv
#   5-SensitivityTesting/results/p_df_analysis_predictions.csv  (extrapolation ranges;
#     used because pred_ranges is not exported by 7-PredictiveModel.R)
#
# NOTE: predictors are supplied in their raw, natural units. n_grid_cells is the
# only model predictor that was log-transformed during training
# (5-IDParamSens.2.R: df$n_grid_cells <- log(df$n_grid_cells)), so it is given
# here as a raw cell count and log()'d internally before scaling. All other
# predictors (comp_*, prop_species, ratio_sdlog, delta_meanlog,
# sample_mean_range_size) enter the model untransformed. The extrapolation
# warning flags inputs that fall outside the training range.
# =============================================================================

# Model predictors that were log-transformed in training (natural log).
.p_log_ivs <- c("n_grid_cells")

.p_results_dir <- file.path(project_root, "5-SensitivityTesting", "results")
p_model_final  <- readRDS(file.path(.p_results_dir, "p_model_final.rds"))
# Signed-error models (Section 4e): used for the completed-richness inversion.
p_model_signed <- readRDS(file.path(.p_results_dir, "p_model_signed.rds"))
p_scale_params <- read.csv(file.path(.p_results_dir, "p_model_scale_params.csv"),
                           stringsAsFactors = FALSE)
p_train_data   <- as.data.frame(
  data.table::fread(file.path(.p_results_dir, "p_df_analysis_predictions.csv"))
)

extend_richness_from_sampling <- function(
    partialGamma,
    comp_spread            = NULL,
    comp_mean              = NULL,
    prop_species           = NULL,
    ratio_sdlog            = .amani_ratio_sdlog,
    n_grid_cells           = NULL,
    comp_sd                = NULL,
    sample_mean_range_size = .amani_delta_meanlog,
    delta_meanlog          = NULL,
    fill_Plants            = .amani_fill_plants,
    interval               = c("prediction", "confidence"),
    level                  = 0.95) {

  interval <- match.arg(interval)
  fp       <- if (isTRUE(fill_Plants)) "TRUE" else "FALSE"
  fit      <- p_model_signed[[fp]]
  sp       <- p_scale_params[p_scale_params$fill_plants == fp, ]
  ext_ivs  <- sp$iv
  train_fp <- p_train_data[p_train_data$model_fill_plants == fp, ]

  # Supplied predictors in raw natural units; NULL -> training median for this level.
  supplied <- list(
    comp_spread = comp_spread, comp_mean = comp_mean, prop_species = prop_species,
    ratio_sdlog = ratio_sdlog, n_grid_cells = n_grid_cells, comp_sd = comp_sd,
    sample_mean_range_size = sample_mean_range_size, delta_meanlog = delta_meanlog
  )

  # Apply the training-time log transform to any supplied log-scale predictor so
  # it matches the model scale. Defaults pulled from train_fp below are already
  # stored on that scale, so they are left as-is.
  for (v in .p_log_ivs) {
    if (!is.null(supplied[[v]])) supplied[[v]] <- log(supplied[[v]])
  }

  # Allow vector inputs (e.g. several sites); recycle to a common length.
  n_row <- max(length(partialGamma),
               vapply(supplied, length, integer(1)))

  newdata <- data.frame(row.names = seq_len(n_row))
  for (v in ext_ivs) {
    val <- supplied[[v]]
    if (is.null(val)) val <- stats::median(train_fp[[v]], na.rm = TRUE)
    newdata[[v]] <- rep_len(val, n_row)
  }

  # Standardize with training mean/sd; flag extrapolation beyond training range.
  scaled  <- data.frame(row.names = seq_len(n_row))
  warn_iv <- character(0)
  for (v in ext_ivs) {
    m  <- sp$mean[sp$iv == v]
    s  <- sp$sd[sp$iv == v]
    lo <- min(train_fp[[v]], na.rm = TRUE)
    hi <- max(train_fp[[v]], na.rm = TRUE)
    if (any(newdata[[v]] < lo | newdata[[v]] > hi, na.rm = TRUE))
      warn_iv <- c(warn_iv, v)
    scaled[[paste0(v, "_s")]] <- (newdata[[v]] - m) / s
  }
  if (length(warn_iv))
    warning("Extrapolation: input(s) outside training range for ",
            paste(warn_iv, collapse = ", "),
            " (fill_Plants = ", fp, "). Prediction is unreliable.")

  pr <- predict(fit, newdata = scaled, interval = interval, level = level)

  # Convert the SIGNED proportional error to a completed richness. The signed
  # model targets signed_prop_error = (estimate - true) / true (Section 4e of
  # 7-PredictiveModel.R), which inverts to:
  #   fullGamma = partialGamma / (1 + signed_prop_error).
  # Negative error -> underestimate (fullGamma > partialGamma); positive ->
  # overestimate. 1 + error <= 0 implies non-positive true richness -> NA.
  pg      <- rep_len(partialGamma, n_row)
  to_full <- function(e) {
    out <- pg / (1 + e)
    out[!is.finite(out) | (1 + e) <= 0] <- NA_real_
    out
  }
  full_fit <- to_full(pr[, "fit"])
  # fullGamma decreases in signed error, so sort the transformed bounds.
  full_lwr <- pmin(to_full(pr[, "lwr"]), to_full(pr[, "upr"]))
  full_upr <- pmax(to_full(pr[, "lwr"]), to_full(pr[, "upr"]))

  # --- Constrain interval bounds (mirrors 7-PredictiveModel.R Section 4e) -----
  # (1) The lower bound cannot fall below the observed partialGamma estimate.
  full_lwr <- pmax(full_lwr, pg, na.rm = TRUE)
  # (2) For fill_Plants = TRUE, cap the upper bound at 5x partialGamma (liberal
  #     ceiling); divergent / NA upper bounds are treated as infinite first.
  if (identical(fp, "TRUE")) {
    upr_capped <- full_upr
    upr_capped[is.na(upr_capped)] <- Inf
    full_upr <- pmin(upr_capped, pg * 5)
  }

  data.frame(
    fill_Plants       = fp,
    partialGamma      = pg,
    signed_prop_error = round(pr[, "fit"], 4),
    fullGamma         = full_fit,
    fullGamma_lwr     = full_lwr,   # lower richness bound (>= partialGamma)
    fullGamma_upr     = full_upr    # upper richness bound (TRUE: <= 5x pg)
  )
}


