# =============================================================================
# 7-PredictiveModel.R
# -----------------------------------------------------------------------------
# PURPOSE
#   Build a single multivariate linear model that predicts the *proportional*
#   number of species missing from a species-richness estimate, given observable
#   sampling characteristics.  The fitted model produces a practical prediction
#   equation that can be applied to real occurrence datasets without knowing the true community.
#
#   Phase 1: fit a parsimonious 3-variable model — easy to visualise in 3D
#             and interpret for a biological audience.
#   Phase 2: expand to all variables that significantly improve fit in combination.
#
# INPUTS  (exported by Scripts 5 & 6)
#   results/p_df_analysis.csv          — full analysis dataset (all IVs + DVs)
#   results/p_vi_magn.csv              — variable importance (magnitude bias)
#   results/p_vi_prop.csv              — variable importance (proportional bias)
#   results/p_sig_counts.csv           — # experiments each IV is significant in
#   results/p_lm_all.csv               — per-experiment LM table
#   results/p_interaction_summary.csv  — significant interactions from Notebook 6
#
# OUTPUTS
#   results/p_model_3var.rds             — fitted 3-variable model object
#   results/p_model_final.rds            — fitted final (extended) model object
#   results/p_model_coefficients.csv     — coefficient table (original + scaled units)
#   results/p_model_scale_params.csv     — centering/scaling parameters per predictor
#   results/p_df_analysis_predictions.csv — training data + model predictions + residuals
#   results/p_model_plots/               — diagnostic and visualisation plots
#
# FLOW
#   1. Setup & load
#   2. Variable selection (candidate pool + collinearity check)
#   3. Simple 3-variable model + diagnostics + 3D visualisation
#   4. Extended model (all significant predictors + key interactions)
#   5. Prediction function and worked examples
#   6. Model export
#   7. Notes and caveats
# =============================================================================


# SETUP ----
{
  library(knitr)
  library(rstudioapi)

  if (exists(".amani_wd")) {
    setwd(.amani_wd)
  } else {
    script_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
    parent_dir <- file.path(script_dir, "../..")
    knitr::opts_knit$set(root.dir = parent_dir)
    setwd(parent_dir)
  }

  # Packages
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(car)      # Type III ANOVA + VIF
  # library(plotly) or library(rgl)  — for 3D visualisation in Section 3c
  # library(MASS)    — for stepAIC in Section 4b
  # library(MuMIn)   — optional: dredge() for all-subsets selection

  # Paths
  results_folder <- "./5-SensitivityTesting/results/"
  model_plot_dir <- file.path(results_folder, "p_model_plots")
  dir.create(model_plot_dir, showWarnings = FALSE, recursive = TRUE)

  # Load exports
  df_analysis   <- fread(file.path(results_folder, "p_df_analysis.csv"))
  vi_prop       <- read.csv(file.path(results_folder, "p_vi_prop.csv"))
  lm_all        <- read.csv(file.path(results_folder, "p_lm_all.csv"))
  int_summary   <- read.csv(file.path(results_folder, "p_interaction_summary.csv"))
  # int_summary is output of Notebook 6; load after running that script.

  # Two-version analysis: fit all models separately for fill_plants = TRUE / FALSE.
  # The 84 floor-experiment rows where fill_plants is NA are dropped in both subsets.
  # NOTE: vi_prop was computed on pooled data; it is used only for predictor RANKING
  # (Section 2a) and is therefore acceptable as a common reference for both subsets.
  fp_levels <- c("FALSE", "TRUE")

  # Primary DV for the predictive model
  # abs_prop_error: unsigned proportional error — tells "how wrong" the estimate is
  # bias_prop: signed proportional bias — tells "which direction and how far"
  # Decision: use abs_prop_error as the primary prediction target.
  # NOTE: confirm this choice matches the downstream use case
  #   (e.g., if the goal is correction, bias_prop is more appropriate).
  primary_dv <- "abs_prop_error"
}


# =============================================================================
# SECTION 1: VARIABLE SELECTION ----


## 1a. Candidate predictor pool ----
# NOTE on ranking: the outline proposed scoring candidates by
#   partial_eta_sq (vi_prop) x N_Sig_Experiments (sig_counts).
# In practice sig_counts is degenerate here — Script 5 fit most per-experiment
# models with a single FOCAL IV each, so only `sampling_evenness` ever lands in
# sig_counts. We therefore rank candidates by vi_prop partial eta-squared alone
# (computed on bias_prop; used as a proxy ordering for abs_prop_error).
{
  # --- Exclusion sets -------------------------------------------------------
  # Simulation-only: require the TRUE community, unknowable for real datasets.
  sim_only <- c("trueSpeciesRichness", "n_species", "cellsPerSpp", "SppPerStartCell")

  # Entirely-NA placeholder columns (unused legacy range-size fields).
  allNA <- names(df_analysis)[sapply(df_analysis, function(x) all(is.na(x)))]

  # Populated only in the completeness-floor experiment (~80% NA when pooled);
  # keeping them would collapse the complete-case sample. Revisit in a
  # floor-specific model if needed.
  floor_only <- c("completeness_floor", "n_cells_filtered",
                  "n_grid_cells_before", "n_grid_cells_after", "min_completeness")

  # Metadata, identifiers, DVs/DV-derived, constants, and categoricals reserved
  # for Phase 2 (this Phase-1 pool is numeric only).
  non_ivs <- c("input_file", "fit_params", "scenario", "test", "best_estimate",
               "lower_bound", "upper_bound", "param_names_fit",
               "delta_meanlog_status", "estimate_direction", "organism_group",
               "prop_sites", "sampling_fraction",
               "bias_magn", "bias_prop", "abs_error", "abs_prop_error",
               "fill_plants", "fill_animals", "filter_plant_completeness",
               "fill_interaction", "clustering_label",
               "sample_endemism_min", "sample_endemism_max")

  numeric_cols  <- names(df_analysis)[sapply(df_analysis, is.numeric)]
  candidate_ivs <- setdiff(numeric_cols,
                           unique(c(sim_only, allNA, floor_only, non_ivs)))

  # --- Annotate each candidate ---------------------------------------------
  # source: main_effect for all Phase-1 candidates (interactions added in 3a).
  # availability: whether the IV is directly OBSERVABLE from occurrence records
  #   or a PROXY that must be ESTIMATED for real data.
  #   • comp_mean/comp_sd/comp_spread and prop_species are computed here from the
  #     known simulated community. For a real dataset they are NOT directly
  #     observed and would need estimation (e.g. Chao-style completeness,
  #     estimated species pool). Flagged as "proxy (estimable)".
  proxy_ivs <- c("comp_mean", "comp_sd", "comp_spread", "prop_species")
  candidate_meta <- data.frame(
    iv           = candidate_ivs,
    source       = "main_effect",
    data_type    = "continuous",
    n_complete   = sapply(candidate_ivs, function(v) sum(!is.na(df_analysis[[v]]))),
    vi_prop_eta  = vi_prop$Partial_Eta_Sq[match(candidate_ivs, vi_prop$IV)],
    availability = ifelse(candidate_ivs %in% proxy_ivs,
                          "proxy (estimable)", "observable"),
    stringsAsFactors = FALSE
  )
  candidate_meta <- candidate_meta[order(-candidate_meta$vi_prop_eta), ]
  rownames(candidate_meta) <- NULL
  message("Candidate pool (", length(candidate_ivs), " numeric IVs):")
  print(candidate_meta)
}

## 1b. Collinearity check ----
# Pairwise correlations on the candidate pool; flag |r| > 0.80 and keep one
# representative per collinear cluster. Representative chosen for interpretability
# / availability and, where tied, higher vi_prop importance.
{
  cor_ivs <- candidate_meta$iv
  cormat  <- cor(as.matrix(df_analysis[, ..cor_ivs]), use = "pairwise.complete.obs")

  # Heatmap of the correlation matrix.
  cor_long <- as.data.frame(as.table(cormat))
  names(cor_long) <- c("v1", "v2", "r")
  print(
    ggplot(cor_long, aes(v1, v2, fill = r)) +
      geom_tile(color = "grey90") +
      geom_text(aes(label = round(r, 2)), size = 2.6) +
      scale_fill_gradient2(limits = c(-1, 1), low = "#b2182b",
                           mid = "white", high = "#2166ac") +
      labs(x = NULL, y = NULL,
           title = "Candidate IV correlation matrix",
           subtitle = "Clusters with |r| > 0.80 collapsed to one representative") +
      theme_minimal(base_size = 10) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  )

  hi <- which(abs(cormat) > 0.80 & upper.tri(cormat), arr.ind = TRUE)
  flagged <- data.frame(v1 = rownames(cormat)[hi[, 1]],
                        v2 = colnames(cormat)[hi[, 2]],
                        r  = round(cormat[hi], 2))
  message("Collinear pairs (|r| > 0.80):")
  print(flagged)

  # --- Resolution decisions (documented) ------------------------------------
  # Range-size cluster: mean/sd/median range size are r 0.73-0.97.
  #   Keep sample_mean_range_size (most interpretable proxy for endemism).
  # Sampling-extent cluster: n_grid_cells / total_observations /
  #   sample_numberStartSites / sample_diffusion are r 0.80-0.95.
  #   Keep n_grid_cells (clearest "number of sampled sites"; fully observed;
  #   central to the filtering steps).
  # Evenness/spread cluster: sampling_evenness ~ comp_spread = -0.91.
  #   Keep comp_spread (higher vi_prop importance; part of the completeness set).
  collinear_drops <- c("sample_sd_range_size", "sample_median_range_size",
                       "total_observations", "sample_numberStartSites",
                       "sample_diffusion", "sampling_evenness")

  final_candidates <- setdiff(candidate_meta$iv, collinear_drops)

  # Confirm the survivors are mutually below threshold.
  cm_final <- cor(as.matrix(df_analysis[, ..final_candidates]),
                  use = "pairwise.complete.obs")
  message("Final candidates (", length(final_candidates), "): ",
          paste(final_candidates, collapse = ", "))
  message("Max |r| among final candidates: ",
          round(max(abs(cm_final[upper.tri(cm_final)])), 2),
          "  (VIF re-checked after fitting in Section 2c).")
}


# =============================================================================
# SECTION 2: SIMPLE 3-VARIABLE MODEL ----

# Select the 3 predictors with the largest consistent effect (highest partial
# eta-sq in vi_prop AND significant in the most experiments in sig_counts) that
# are also available from real occurrence records.
# 3 variables → easy 3D surface plot and a clean introductory result.


## 2a. Select the 3 predictors ----
# Rank final_candidates by vi_prop partial eta-squared (pooled-data proxy for
# abs_prop_error importance). The same top-3 set is used as a starting point
# for both fill_plants subsets; within-subset importance is re-examined in 2b.
{
  rank_df <- candidate_meta[candidate_meta$iv %in% final_candidates,
                             c("iv", "vi_prop_eta", "availability")]
  rank_df <- rank_df[order(-rank_df$vi_prop_eta), ]
  rownames(rank_df) <- NULL

  message("Candidates ranked by vi_prop partial eta-squared:")
  print(rank_df)

  # Top 3 predictors
  IV1 <- rank_df$iv[1]   # comp_spread  — completeness heterogeneity
  IV2 <- rank_df$iv[2]   # comp_mean    — mean site-level completeness
  IV3 <- rank_df$iv[3]   # prop_species — community representation breadth
  message("\nSelected: IV1 = ", IV1, ", IV2 = ", IV2, ", IV3 = ", IV3)
}


## 2b. Fit the 3-variable model (both fill_plants levels) ----
# All IVs are continuous and fully observed (n = 404 minus floor-experiment NAs).
# Models fitted separately for fill_plants = FALSE and TRUE.
# Centering/scaling is done within each subset (separate scale_params).
{
  model_3var      <- list()   # fitted lm objects, keyed by fill_plants level
  scale_params_3var <- list() # per-predictor mean + sd used in scaling
  df_3var         <- list()   # complete-case modelling data per level

  for (fp in fp_levels) {
    sub <- as.data.frame(df_analysis[df_analysis$fill_plants == fp, ])
    sub <- sub[complete.cases(sub[, c(primary_dv, IV1, IV2, IV3)]), ]

    sp <- data.frame(
      iv   = c(IV1, IV2, IV3),
      mean = c(mean(sub[[IV1]]), mean(sub[[IV2]]), mean(sub[[IV3]])),
      sd   = c(sd(sub[[IV1]]),   sd(sub[[IV2]]),   sd(sub[[IV3]]))
    )

    iv1s <- paste0(IV1, "_s"); sub[[iv1s]] <- as.numeric(scale(sub[[IV1]]))
    iv2s <- paste0(IV2, "_s"); sub[[iv2s]] <- as.numeric(scale(sub[[IV2]]))
    iv3s <- paste0(IV3, "_s"); sub[[iv3s]] <- as.numeric(scale(sub[[IV3]]))

    fml <- as.formula(paste0(primary_dv, " ~ ", iv1s, " + ", iv2s, " + ", iv3s))
    fit <- lm(fml, data = sub)

    model_3var[[fp]]        <- fit
    scale_params_3var[[fp]] <- sp
    df_3var[[fp]]           <- sub

    s <- summary(fit)
    message("\n========== fill_plants = ", fp, "  (n = ", nrow(sub), ") ==========")
    message("Formula (scaled IVs): ", deparse(formula(fit)))
    message("R\u00b2 = ", round(s$r.squared, 3),
            "   adj. R\u00b2 = ", round(s$adj.r.squared, 3))
    coef_tbl <- cbind(round(coef(s), 4),
                      round(confint(fit), 3))
    colnames(coef_tbl)[5:6] <- c("2.5%", "97.5%")
    print(coef_tbl)

    if (s$adj.r.squared < 0.15)
      warning("  adj. R\u00b2 < 0.15 for fill_plants = ", fp,
              " — 3D surface will be uninformative. Consider adding IVs.")
  }
}


## 2c. Model diagnostics (both fill_plants levels) ----
{
  library(patchwork)

  for (fp in fp_levels) {
    fit <- model_3var[[fp]]
    sub <- df_3var[[fp]]
    n   <- nrow(sub)
    cook_cut <- 4 / n

    diag_df <- data.frame(
      fitted    = fitted(fit),
      residuals = residuals(fit),
      std_resid = rstandard(fit),
      sqrt_abs  = sqrt(abs(rstandard(fit))),
      cooks     = cooks.distance(fit)
    )

    p_rf <- ggplot(diag_df, aes(fitted, residuals)) +
      geom_point(alpha = 0.5) + geom_hline(yintercept = 0, linetype = "dashed") +
      geom_smooth(method = "loess", se = FALSE, color = "steelblue", formula = y~x) +
      labs(x = "Fitted", y = "Residuals", title = "Residuals vs Fitted") +
      theme_minimal(base_size = 10)

    p_qq <- ggplot(diag_df, aes(sample = std_resid)) +
      stat_qq(alpha = 0.5) + stat_qq_line(color = "steelblue") +
      labs(x = "Theoretical quantiles", y = "Standardised residuals",
           title = "Normal Q-Q") +
      theme_minimal(base_size = 10)

    p_sl <- ggplot(diag_df, aes(fitted, sqrt_abs)) +
      geom_point(alpha = 0.5) +
      geom_smooth(method = "loess", se = FALSE, color = "steelblue", formula = y~x) +
      labs(x = "Fitted", y = "\u221a|Std. residuals|", title = "Scale-Location") +
      theme_minimal(base_size = 10)

    p_ck <- ggplot(diag_df, aes(seq_len(n), cooks)) +
      geom_col(alpha = 0.6) +
      geom_hline(yintercept = cook_cut, linetype = "dashed", color = "red") +
      labs(x = "Observation index", y = "Cook\u2019s D",
           title = paste0("Cook's Distance  (cut = ", round(cook_cut, 4), ")")) +
      theme_minimal(base_size = 10)

    print(
      wrap_plots(p_rf, p_qq, p_sl, p_ck, ncol = 2) +
        plot_annotation(
          title    = paste0("Diagnostics — 3-variable model  |  fill_plants = ", fp),
          subtitle = paste0("n = ", n, ",  adj. R\u00b2 = ",
                            round(summary(fit)$adj.r.squared, 3)),
          theme    = theme_minimal(base_size = 12)
        )
    )

    # VIF (car::vif on the fitted model)
    message("\nVIF — fill_plants = ", fp, ":")
    print(round(car::vif(fit), 2))

    # Save diagnostic plot
    fname <- file.path(model_plot_dir,
                       paste0("diag_3var_fill_plants_", fp, ".png"))
    ggsave(fname, width = 10, height = 8, dpi = 150)
    message("Saved: ", fname)
  }

  # NOTE: abs_prop_error is bounded (0, 1) and right-skewed by construction.
  # If Q-Q plots show heavy tails or Scale-Location shows a rising trend,
  # consider re-fitting with a Gamma GLM (log link) in Section 3:
  #   glm(abs_prop_error ~ ..., family = Gamma(link = "log"))
}


## 2d. 3D visualisation (placeholder — requires plotly or rgl) ----
# Continuous × continuous × continuous → use plotly response surface.
# IV1 (comp_spread) on x, IV2 (comp_mean) on y, z = predicted abs_prop_error.
# IV3 (prop_species) shown via colour of observed points and surface sliders.
# Run this block interactively; static PNG fallback saved for manuscript.
{
  library(plotly)

  for (fp in fp_levels) {
    fit <- model_3var[[fp]]
    sub <- df_3var[[fp]]
    sp  <- scale_params_3var[[fp]]

    iv1_seq <- seq(min(sub[[IV1]]), max(sub[[IV1]]), length.out = 40)
    iv2_seq <- seq(min(sub[[IV2]]), max(sub[[IV2]]), length.out = 40)
    grid    <- expand.grid(IV1 = iv1_seq, IV2 = iv2_seq)
    grid[[paste0(IV1,"_s")]] <- (grid$IV1 - sp$mean[1]) / sp$sd[1]
    grid[[paste0(IV2,"_s")]] <- (grid$IV2 - sp$mean[2]) / sp$sd[2]
    grid[[paste0(IV3,"_s")]] <- 0   # IV3 held at its mean (0 after scaling)
    grid$z <- predict(fit, newdata = grid)
    z_mat  <- matrix(grid$z, 40, 40)

    fig <- plot_ly(x = iv1_seq, y = iv2_seq, z = z_mat, type = "surface") |>
      layout(title = paste0("fill_plants = ", fp),
             scene = list(xaxis = list(title = IV1),
                          yaxis = list(title = IV2),
                          zaxis = list(title = primary_dv)))
    htmlwidgets::saveWidget(fig,
      file.path(model_plot_dir, paste0("surface_3var_fp", fp, ".html")))
  }
}


# =============================================================================
# SECTION 3: EXTENDED MODEL ----

# Expand beyond the 3-variable model to include all predictors that significantly
# improve fit. Interaction terms from Notebook 6 are added here.


## 3a. Full candidate models (both fill_plants levels) + pooled interaction ----
# Two parallel extended models — one per fill_plants level — each using all 8
# final_candidates. Because fill_plants is constant within a subset, the H1
# interaction (fill_plants x prop_species) cannot be fit there; it is fit
# separately on the POOLED data (fill_plants NA rows dropped) so the interaction
# coefficient can be reported and formally tested.
{
  ext_ivs <- final_candidates   # 8 IVs surviving Section 1b

  # Helper: scale a set of IVs in-place, appending "_s" columns.
  scale_ivs <- function(d, ivs) {
    for (v in ivs) d[[paste0(v, "_s")]] <- as.numeric(scale(d[[v]]))
    d
  }

  model_ext        <- list()
  scale_params_ext <- list()
  df_ext           <- list()

  for (fp in fp_levels) {
    sub <- as.data.frame(df_analysis[df_analysis$fill_plants == fp, ])
    sub <- sub[complete.cases(sub[, c(primary_dv, ext_ivs)]), ]
    scale_params_ext[[fp]] <- data.frame(
      iv   = ext_ivs,
      mean = sapply(ext_ivs, function(v) mean(sub[[v]])),
      sd   = sapply(ext_ivs, function(v) sd(sub[[v]]))
    )
    sub <- scale_ivs(sub, ext_ivs)
    fml <- as.formula(paste0(primary_dv, " ~ ",
                             paste0(ext_ivs, "_s", collapse = " + ")))
    fit <- lm(fml, data = sub)
    model_ext[[fp]] <- fit
    df_ext[[fp]]    <- sub

    message("\n===== Extended model | fill_plants = ", fp,
            "  (n = ", nrow(sub), ") =====")
    message("adj. R\u00b2 = ", round(summary(fit)$adj.r.squared, 3),
            "   AIC = ", round(AIC(fit), 1))
    print(car::Anova(fit, type = "II"))
  }

  # --- Pooled model with the explicit H1 interaction term ---
  pooled <- as.data.frame(df_analysis[df_analysis$fill_plants %in% fp_levels, ])
  pooled <- pooled[complete.cases(pooled[, c(primary_dv, ext_ivs)]), ]
  scale_params_pooled <- data.frame(
    iv   = ext_ivs,
    mean = sapply(ext_ivs, function(v) mean(pooled[[v]])),
    sd   = sapply(ext_ivs, function(v) sd(pooled[[v]]))
  )
  pooled <- scale_ivs(pooled, ext_ivs)
  pooled$fill_plants <- factor(pooled$fill_plants, levels = c("FALSE", "TRUE"))
  # Treatment contrasts so fill_plantsTRUE is interpreted vs the FALSE reference.
  contrasts(pooled$fill_plants) <- contr.treatment(2)

  fml_pool <- as.formula(paste0(
    primary_dv, " ~ ", paste0(ext_ivs, "_s", collapse = " + "),
    " + fill_plants + fill_plants:prop_species_s"
  ))
  model_pooled_int <- lm(fml_pool, data = pooled)

  message("\n===== Pooled model + H1 interaction  (n = ", nrow(pooled), ") =====")
  message("adj. R\u00b2 = ", round(summary(model_pooled_int)$adj.r.squared, 3))
  print(car::Anova(model_pooled_int, type = "III")["prop_species_s:fill_plants", ])
  message("Interaction beta (prop_species_s:fill_plantsTRUE) = ",
          round(coef(model_pooled_int)["prop_species_s:fill_plants2"], 4))
}


## 3b. Model selection (AIC) per fill_plants level ----
# Compare, within each level and on the SAME complete-case rows:
#   Model 0 (intercept), Model 3 (Section 2 three-var), Model F (all 8),
#   Model S (backward stepAIC from F). Select final by AIC + interpret-ability.
{
  model_final <- list()
  aic_tables  <- list()

  for (fp in fp_levels) {
    sub <- df_ext[[fp]]
    m0 <- lm(as.formula(paste0(primary_dv, " ~ 1")), data = sub)
    m3 <- lm(as.formula(paste0(primary_dv, " ~ ",
             paste0(c(IV1, IV2, IV3), "_s", collapse = " + "))), data = sub)
    mF <- model_ext[[fp]]
    mS <- MASS::stepAIC(mF, direction = "backward", trace = FALSE)

    # Biological override: keep comp_mean even if stepAIC drops it.
    if (!("comp_mean_s" %in% attr(terms(mS), "term.labels"))) {
      message("  [override] stepAIC dropped comp_mean for fill_plants = ", fp,
              "; adding it back (biologically important).")
      mS <- update(mS, . ~ . + comp_mean_s)
    }

    aic_tables[[fp]] <- data.frame(
      model = c("M0_null", "M3_threevar", "MF_full", "MS_step"),
      k     = c(length(coef(m0)), length(coef(m3)),
                length(coef(mF)), length(coef(mS))),
      AIC   = c(AIC(m0), AIC(m3), AIC(mF), AIC(mS)),
      adjR2 = c(summary(m0)$adj.r.squared, summary(m3)$adj.r.squared,
                summary(mF)$adj.r.squared, summary(mS)$adj.r.squared)
    )
    aic_tables[[fp]]$dAIC <- aic_tables[[fp]]$AIC - min(aic_tables[[fp]]$AIC)

    # Final = stepAIC-selected model (with any override applied).
    model_final[[fp]] <- mS

    message("\n----- AIC comparison | fill_plants = ", fp, " -----")
    print(aic_tables[[fp]])
    message("Final model terms: ",
            paste(attr(terms(mS), "term.labels"), collapse = " + "))
  }
}


## 3c. Final model diagnostics ----
{
  library(patchwork)

  diag_panel <- function(fit, label) {
    n <- length(residuals(fit)); cut <- 4 / n
    dd <- data.frame(fitted = fitted(fit), residuals = residuals(fit),
                     std = rstandard(fit), sqrt_abs = sqrt(abs(rstandard(fit))),
                     cooks = cooks.distance(fit))
    p1 <- ggplot(dd, aes(fitted, residuals)) + geom_point(alpha = 0.5) +
      geom_hline(yintercept = 0, linetype = "dashed") +
      geom_smooth(method = "loess", se = FALSE, color = "steelblue", formula = y~x) +
      labs(title = "Residuals vs Fitted", x = "Fitted", y = "Residuals") +
      theme_minimal(base_size = 10)
    p2 <- ggplot(dd, aes(sample = std)) + stat_qq(alpha = 0.5) +
      stat_qq_line(color = "steelblue") +
      labs(title = "Normal Q-Q", x = "Theoretical", y = "Std. resid.") +
      theme_minimal(base_size = 10)
    p3 <- ggplot(dd, aes(fitted, sqrt_abs)) + geom_point(alpha = 0.5) +
      geom_smooth(method = "loess", se = FALSE, color = "steelblue", formula = y~x) +
      labs(title = "Scale-Location", x = "Fitted", y = "\u221a|Std. resid.|") +
      theme_minimal(base_size = 10)
    p4 <- ggplot(dd, aes(seq_len(n), cooks)) + geom_col(alpha = 0.6) +
      geom_hline(yintercept = cut, linetype = "dashed", color = "red") +
      labs(title = paste0("Cook's D (cut=", round(cut, 4), ")"),
           x = "Obs", y = "Cook's D") + theme_minimal(base_size = 10)
    wrap_plots(p1, p2, p3, p4, ncol = 2) +
      plot_annotation(title = label, theme = theme_minimal(base_size = 12))
  }

  for (fp in fp_levels) {
    print(diag_panel(model_final[[fp]],
                     paste0("Final extended model | fill_plants = ", fp)))
    message("VIF — final model, fill_plants = ", fp, ":")
    print(round(car::vif(model_final[[fp]]), 2))
    ggsave(file.path(model_plot_dir, paste0("diag_final_fill_plants_", fp, ".png")),
           width = 10, height = 8, dpi = 150)
  }

  # NOTE: abs_prop_error is bounded in (0,1). If diagnostics show strong
  # heteroscedasticity, refit as Gamma(link="log") and compare.
}


## 3d. Final coefficient table ----
# One tidy table per fill_plants level: scaled beta (per 1 SD), original-unit
# beta (= beta_scaled / sd), 95% CI, SE, t, p, plus a plain-language reading.
{
  coef_table_list <- list()

  for (fp in fp_levels) {
    fit <- model_final[[fp]]
    sp  <- scale_params_ext[[fp]]
    sm  <- summary(fit)$coefficients
    ci  <- confint(fit)
    terms_s <- rownames(sm)

    # Map "<iv>_s" back to raw IV + its SD for original-unit conversion.
    raw_iv <- sub("_s$", "", terms_s)
    sd_iv  <- sp$sd[match(raw_iv, sp$iv)]   # NA for intercept

    tbl <- data.frame(
      predictor     = terms_s,
      beta_scaled   = round(sm[, "Estimate"], 4),
      beta_original = round(ifelse(is.na(sd_iv), NA, sm[, "Estimate"] / sd_iv), 5),
      ci_lo         = round(ci[, 1], 4),
      ci_hi         = round(ci[, 2], 4),
      SE            = round(sm[, "Std. Error"], 4),
      t             = round(sm[, "t value"], 2),
      p             = signif(sm[, "Pr(>|t|)"], 3),
      stringsAsFactors = FALSE
    )
    tbl$sig <- ifelse(tbl$p < 0.05, "*", "")
    tbl$plain <- ifelse(
      is.na(sd_iv), "",
      paste0("+1 SD in ", raw_iv, " -> ",
             round(sm[, "Estimate"] * 100, 2),
             " pp change in expected abs_prop_error"))
    rownames(tbl) <- NULL
    coef_table_list[[fp]] <- tbl

    message("\n----- Coefficient table | fill_plants = ", fp, " -----")
    print(tbl)
  }

  # Combine for export (tag the fill_plants level).
  coef_table <- do.call(rbind, lapply(fp_levels, function(fp) {
    cbind(fill_plants = fp, coef_table_list[[fp]])
  }))
  rownames(coef_table) <- NULL
}


# =============================================================================
# SECTION 4: PREDICTION FUNCTION AND WORKED EXAMPLES ----

# Uses model_final (full 8-var model for fill_plants = FALSE; stepAIC+override
# 5-var model for TRUE). We keep the full FALSE model for PREDICTION ACCURACY.
# CAVEAT: the completeness-family predictors are jointly collinear (VIF > 10 in
# the FALSE model), so INDIVIDUAL coefficients are unstable and must NOT be read
# as marginal effects — the model is fit-for-purpose for prediction only.

## 4a. Training ranges + medians (for extrapolation checks & example defaults) ----
{
  pred_ranges <- setNames(lapply(fp_levels, function(fp) {
    d <- df_ext[[fp]]
    data.frame(iv = ext_ivs,
               min = sapply(ext_ivs, function(v) min(d[[v]])),
               max = sapply(ext_ivs, function(v) max(d[[v]])))
  }), fp_levels)

  train_median <- setNames(lapply(fp_levels, function(fp) {
    setNames(sapply(ext_ivs, function(v) median(df_ext[[fp]][[v]])), ext_ivs)
  }), fp_levels)
}

## 4b. Prediction function ----
# newdata: data.frame with RAW (unscaled) columns for all 8 ext_ivs.
# fill_plants: which model to apply ("FALSE" or "TRUE").
# Returns point estimate + interval; warns on extrapolation beyond training range.
predict_missing_prop <- function(newdata, fill_plants = c("FALSE", "TRUE"),
                                  interval = c("prediction", "confidence"),
                                  level = 0.95) {
  fill_plants <- match.arg(fill_plants)
  interval    <- match.arg(interval)
  fit <- model_final[[fill_plants]]
  sp  <- scale_params_ext[[fill_plants]]
  rng <- pred_ranges[[fill_plants]]

  newdata <- as.data.frame(newdata)
  missing_iv <- setdiff(ext_ivs, names(newdata))
  if (length(missing_iv))
    stop("newdata missing required IV column(s): ", paste(missing_iv, collapse = ", "))

  # Scale each raw IV with the training mean/sd; flag extrapolation.
  scaled  <- data.frame(row.names = seq_len(nrow(newdata)))
  warn_iv <- character(0)
  for (i in seq_len(nrow(sp))) {
    v   <- sp$iv[i]
    raw <- newdata[[v]]
    lo  <- rng$min[rng$iv == v]; hi <- rng$max[rng$iv == v]
    if (any(raw < lo | raw > hi, na.rm = TRUE)) warn_iv <- c(warn_iv, v)
    scaled[[paste0(v, "_s")]] <- (raw - sp$mean[i]) / sp$sd[i]
  }
  if (length(warn_iv))
    warning("Extrapolation: input(s) outside training range for ",
            paste(warn_iv, collapse = ", "),
            " (fill_plants = ", fill_plants, "). Prediction is unreliable.")

  pr <- predict(fit, newdata = scaled, interval = interval, level = level)
  data.frame(fill_plants = fill_plants,
             pred_abs_prop_error = round(pr[, "fit"], 4),
             lwr = round(pr[, "lwr"], 4),
             upr = round(pr[, "upr"], 4))
}

## 4c. Worked examples ----
# Vary the sampling-quality IVs (n_grid_cells, comp_mean, prop_species) across
# their 10th / 50th / 90th training percentiles; hold the rest at the median.
{
  qtl <- function(fp, v, p) as.numeric(quantile(df_ext[[fp]][[v]], p))
  make_example <- function(fp, p, label) {
    x <- as.list(train_median[[fp]])
    for (v in c("n_grid_cells", "comp_mean", "prop_species")) x[[v]] <- qtl(fp, v, p)
    cbind(scenario = label, as.data.frame(x, stringsAsFactors = FALSE))
  }

  # Build the poor/medium/well gradient on the FALSE training distribution and
  # score it under BOTH models, so the effect of plant filling is visible.
  examples <- rbind(
    make_example("FALSE", 0.10, "poor"),
    make_example("FALSE", 0.50, "medium"),
    make_example("FALSE", 0.90, "well")
  )

  worked <- do.call(rbind, lapply(fp_levels, function(fp) {
    cbind(examples[, "scenario", drop = FALSE],
          predict_missing_prop(examples[, ext_ivs], fill_plants = fp))
  }))
  message("\n----- Worked examples (predicted abs_prop_error + 95% PI) -----")
  print(worked)

  # Example 4 — real Peruvian Blattodea (roach) dataset.
  # Completeness IVs come from the observed occurrence data (roach_completeness,
  # computed in Script 5). The remaining IVs would likewise be computed from the
  # real dataset; here they default to training medians and are flagged.
  real_example <- tryCatch({
    rc <- compute_completeness_stats(
      compute_site_completeness(roach_occ, true_richness))
    x <- as.list(train_median[["FALSE"]])
    x$comp_mean <- rc$comp_mean; x$comp_sd <- rc$comp_sd
    x$comp_spread <- rc$comp_spread
    message("\n----- Example 4: real roach dataset (comp_* observed; other IVs = training median) -----")
    predict_missing_prop(as.data.frame(x)[, ext_ivs], fill_plants = "FALSE")
  }, error = function(e) {
    message("Example 4 skipped (need roach_occ + true_richness from Script 5): ",
            conditionMessage(e))
    NULL
  })
  if (!is.null(real_example)) print(real_example)
}

## 4d. Sensitivity sweep (marginal contribution of prop_species) ----
# Sweep prop_species across its observed range, holding other IVs at the median,
# under both fill_plants models. Shows the H1 flattening as a prediction curve.
{
  sweep_v <- "prop_species"
  sweep_df <- do.call(rbind, lapply(fp_levels, function(fp) {
    seqv <- seq(pred_ranges[[fp]]$min[pred_ranges[[fp]]$iv == sweep_v],
                pred_ranges[[fp]]$max[pred_ranges[[fp]]$iv == sweep_v],
                length.out = 60)
    nd <- as.data.frame(lapply(train_median[[fp]], rep, length(seqv)))
    nd[[sweep_v]] <- seqv
    cbind(prop_species = seqv,
          suppressWarnings(predict_missing_prop(nd[, ext_ivs], fill_plants = fp)))
  }))

  print(
    ggplot(sweep_df, aes(prop_species, pred_abs_prop_error, color = fill_plants)) +
      geom_ribbon(aes(ymin = lwr, ymax = upr, fill = fill_plants),
                  alpha = 0.15, color = NA) +
      geom_line(linewidth = 1) +
      scale_color_manual(values = c("FALSE" = "grey40", "TRUE" = "#1b9e77")) +
      scale_fill_manual(values = c("FALSE" = "grey40", "TRUE" = "#1b9e77")) +
      labs(x = "prop_species (community representation)",
           y = "Predicted abs_prop_error",
           title = "Sensitivity: predicted error vs prop_species",
           subtitle = "Other IVs held at training median; band = 95% prediction interval") +
      theme_minimal(base_size = 11)
  )
  ggsave(file.path(model_plot_dir, "sensitivity_prop_species.png"),
         width = 8, height = 5, dpi = 150)
}


# =============================================================================
# SECTION 4e: TOTAL EXPECTED SPECIES FROM SIMULATION RESULTS ----
# Apply the two final per-fill_plants models (model_final[["FALSE"]] and
# model_final[["TRUE"]]) to the sensitivity-test results to recover the total
# (true) expected richness implied by each sampled best_estimate.
#
# --- IV mapping (results column -> model IV) and transforms --------------------
#   prop_species           -> prop_species            (as-is)
#   ratio_sdlog            -> ratio_sdlog             (as-is)
#   n_grid_cells           -> n_grid_cells            (LOG: training used log)
#   sample_mean_range_size -> sample_mean_range_size  (as-is)
#   delta_meanlog          -> delta_meanlog           (as-is)
#   comp_spread/comp_mean/comp_sd  constructed from sampling_evenness,
#     prop_species, and completeness_floor (see build_comp_stats()).
# Only n_grid_cells is log-transformed: training-vs-results ranges match under
# log (log(6) = 1.792, log(2000) = 7.601); every other numeric IV is on native
# scale. The comp_* stats live on the [0, 1] completeness scale (never logged).
# predict_missing_prop() centres/scales internally with scale_params_ext, so IVs
# are passed here in raw (pre-scaling but post-log) units.
#
# --- Signed error DV and inverting to a species count -------------------------
# The original models were fit to abs_prop_error = |best_estimate - true| / true,
# which cannot separate over- from under-estimates. Here we REFIT each of the two
# final models to a SIGNED proportional error, reusing the same selected
# predictors per fill_plants level:
#     signed_prop_error = (best_estimate - true) / true   (= bias_magn / true).
# So abs_prop_error = |signed_prop_error|. This inverts to:
#     total_expected = best_estimate / (1 + predicted_signed_prop_error).
# Negative predicted error -> overestimate (total < best_estimate); positive ->
# underestimate (total > best_estimate). Because total decreases in the signed
# error, the LOWER total bound comes from the UPPER error bound and vice versa.
# Rows implying non-positive true richness (1 + error <= 0) are flagged NA.

#' Construct completeness summary stats from the sampling design
#'
#' Builds the three completeness predictors (comp_mean, comp_spread, comp_sd)
#' that the models expect but that are not stored directly in the results table.
#' Logic mirrors how the sampler generates per-site completeness:
#' \itemize{
#'   \item \code{even}: perfect evenness, so spread = 0, mean = prop_species,
#'     sd = 0.
#'   \item \code{low} / \code{veryLow}: completeness is spread between a floor and
#'     prop_species. When \code{completeness_floor} is supplied the spread is
#'     \code{prop_species - completeness_floor}; otherwise it defaults to
#'     \code{prop_species - prop_species / divisor} (divisor = 2 for "low",
#'     10 for "veryLow"). The mean sits at the midpoint and the sd assumes a
#'     uniform distribution over the spread (sd = spread / sqrt(12)).
#' }
#'
#' @param sampling_evenness character vector: "even", "low", or "veryLow".
#' @param prop_species numeric vector of community representation.
#' @param completeness_floor numeric vector; NA where undefined.
#' @return data.frame with columns comp_mean, comp_spread, comp_sd.
build_comp_stats <- function(sampling_evenness, prop_species, completeness_floor) {
  n           <- length(prop_species)
  comp_spread <- numeric(n)
  comp_mean   <- numeric(n)
  comp_sd     <- numeric(n)

  for (i in seq_len(n)) {
    sampling_evenness_i  <- sampling_evenness[i]
    prop_species_i       <- prop_species[i]
    completeness_floor_i <- completeness_floor[i]

    if (identical(sampling_evenness_i, "even")) {
      # Perfectly even sampling: no spread in per-site completeness.
      comp_spread[i] <- 0
      comp_mean[i]   <- prop_species_i
      comp_sd[i]     <- 0
    } else {
      # "low" halves prop_species; "veryLow" divides by 10 for the default spread.
      divisor <- if (identical(sampling_evenness_i, "veryLow")) 10 else 2
      if (!is.na(completeness_floor_i)) {
        comp_spread[i] <- prop_species_i - completeness_floor_i
      } else {
        comp_spread[i] <- prop_species_i - prop_species_i / divisor
      }
      comp_mean[i] <- prop_species_i - comp_spread[i] / 2
      # Uniform over an interval of width comp_spread -> sd = width / sqrt(12).
      comp_sd[i]   <- comp_spread[i] / sqrt(12)
    }
  }

  data.frame(comp_mean = comp_mean, comp_spread = comp_spread, comp_sd = comp_sd)
}

{
  verbose <- TRUE
  fp_results <- file.path(results_folder, "p_comparison_results.csv")
  if (verbose) message("Loading simulation results: ", fp_results)
  res_sim <- as.data.frame(fread(fp_results))

  # Keep only rows that map to one of the two fitted models and have an estimate.
  keep <- res_sim$fill_plants %in% fp_levels & !is.na(res_sim$best_estimate)
  n_drop <- sum(!keep)
  if (verbose && n_drop > 0)
    message("Dropping ", n_drop, " row(s) with NA fill_plants / best_estimate ",
            "(no matching model).")
  res_pred <- res_sim[keep, ]

  # Construct the three completeness predictors from the sampling design.
  comp <- build_comp_stats(res_pred$sampling_evenness,
                           res_pred$prop_species,
                           res_pred$completeness_floor)

  # Assemble raw (pre-scaling) IVs in the units the models were trained on.
  newdata <- data.frame(
    comp_spread            = comp$comp_spread,
    comp_mean              = comp$comp_mean,
    prop_species           = res_pred$prop_species,
    ratio_sdlog            = res_pred$ratio_sdlog,
    n_grid_cells           = log(res_pred$n_grid_cells),
    comp_sd                = comp$comp_sd,
    sample_mean_range_size = res_pred$sample_mean_range_size,
    delta_meanlog          = res_pred$delta_meanlog
  )

  # Refit the two final models on the SIGNED proportional error, reusing each
  # model's selected predictors (only the DV changes).
  model_signed <- list()
  for (fp in fp_levels) {
    train_fp <- df_ext[[fp]]
    train_fp$signed_prop_error <- train_fp$bias_magn / exp(train_fp$trueSpeciesRichness)
    rhs   <- attr(terms(model_final[[fp]]), "term.labels")
    fml_s <- reformulate(rhs, response = "signed_prop_error")
    model_signed[[fp]] <- lm(fml_s, data = train_fp)
    if (verbose) message("Refit signed model | fill_plants = ", fp,
                         "  (adj R\u00b2 = ",
                         round(summary(model_signed[[fp]])$adj.r.squared, 3), ")")
  }

  # Predict signed error + 95% prediction intervals with the matching model.
  # Scale each raw IV with the training mean/sd from scale_params_ext.
  pred <- data.frame(pred_signed = rep(NA_real_, nrow(newdata)),
                     pred_lwr = NA_real_, pred_upr = NA_real_)
  for (fp in fp_levels) {
    idx <- res_pred$fill_plants == fp
    if (!any(idx)) next
    if (verbose) message("Predicting ", sum(idx), " row(s) with model_signed[[\"",
                         fp, "\"]].")
    sp     <- scale_params_ext[[fp]]
    scaled <- data.frame(row.names = seq_len(sum(idx)))
    for (i in seq_len(nrow(sp)))
      scaled[[paste0(sp$iv[i], "_s")]] <-
        (newdata[idx, sp$iv[i]] - sp$mean[i]) / sp$sd[i]
    pr <- predict(model_signed[[fp]], newdata = scaled,
                  interval = "prediction", level = 0.95)
    pred$pred_signed[idx] <- pr[, "fit"]
    pred$pred_lwr[idx]    <- pr[, "lwr"]
    pred$pred_upr[idx]    <- pr[, "upr"]
  }

  # Invert to total expected richness: total = best_estimate / (1 + error).
  # Guard against 1 + error <= 0 (implies non-positive true richness).
  invert_total <- function(best_estimate, err) {
    out <- best_estimate / (1 + err)
    out[!is.finite(out) | (1 + err) <= 0] <- NA_real_
    out
  }
  res_pred$pred_signed_prop_error <- pred$pred_signed
  res_pred$total_expected_species <- invert_total(res_pred$best_estimate,
                                                  pred$pred_signed)
  # total decreases in signed error: lower total <- upper error, and vice versa.
  res_pred$total_expected_lwr <- invert_total(res_pred$best_estimate, pred$pred_upr)
  res_pred$total_expected_upr <- invert_total(res_pred$best_estimate, pred$pred_lwr)

  # --- Constrain interval bounds --------------------------------------------
  # (1) The lower bound cannot fall below the observed best_estimate: the true
  #     total is not expected to be smaller than the value already estimated.
  res_pred$total_expected_lwr <- pmax(res_pred$total_expected_lwr,
                                      res_pred$best_estimate, na.rm = TRUE)
  # (2) For fill_plants = TRUE, cap the upper bound at 5x best_estimate (a
  #     liberal ceiling given behavior seen in other simulations). Divergent /
  #     NA upper bounds are treated as infinite and pulled down to the cap.
  cap_idx <- res_pred$fill_plants == "TRUE"
  upr_true <- res_pred$total_expected_upr[cap_idx]
  upr_true[is.na(upr_true)] <- Inf
  res_pred$total_expected_upr[cap_idx] <-
    pmin(upr_true, res_pred$best_estimate[cap_idx] * 5)

  n_bad <- sum(is.na(res_pred$total_expected_species))
  if (verbose && n_bad > 0)
    message("Flagged ", n_bad, " row(s) with 1 + error <= 0 (total set NA).")

  # Persist and log input parameters alongside results.
  out_cols <- c("scenario", "test", "fill_plants", "sampling_evenness",
                "best_estimate", "trueSpeciesRichness", "prop_species",
                "completeness_floor", "pred_signed_prop_error",
                "total_expected_species", "total_expected_lwr",
                "total_expected_upr", "input_file")
  out_cols <- intersect(out_cols, names(res_pred))
  fp_out <- file.path(results_folder, "p_total_expected_species.csv")
  write.csv(res_pred[, out_cols], fp_out, row.names = FALSE)
  if (verbose) message("Wrote ", nrow(res_pred), " predictions to ", fp_out)

  print(utils::head(res_pred[, out_cols], 10))

  # --- Recovery plot: predicted total vs comp_mean, faceted by fill_plants ---
  # Attach the constructed completeness stats for plotting.
  res_pred$comp_mean   <- comp$comp_mean
  res_pred$comp_spread <- comp$comp_spread

  # true richness reference line (constant across the simulation results).
  true_ref <- stats::median(res_pred$trueSpeciesRichness, na.rm = TRUE)

  # Points at the same comp_mean overlap; jitter x slightly for visibility.
  pdat <- res_pred[is.finite(res_pred$total_expected_species) &
                     is.finite(res_pred$comp_mean), ]
  set.seed(4172)
  pdat$comp_mean_j <- pdat$comp_mean + runif(nrow(pdat), -0.02, 0.02)

  print(
    ggplot(pdat, aes(comp_mean_j, total_expected_species, color = comp_spread)) +
      geom_hline(yintercept = true_ref, linetype = "dashed", color = "grey30") +
      geom_linerange(aes(ymin = total_expected_lwr, ymax = total_expected_upr),
                     alpha = 0.3) +
      geom_point(size = 1.5, alpha = 0.5) +
      facet_wrap(~ fill_plants, labeller = label_both) +
      scale_y_log10(labels = scales::label_comma()) +
      scale_color_viridis_c(option = "cividis")+
      labs(x = "comp_mean (jittered)",
           y = "Predicted total expected species",
           color = "comp_spread",
           title = "Recovery of total expected species (signed-error model)",
           subtitle = paste0("Lines = 95% PI (lower \u2265 best_estimate; ",
                             "TRUE upper capped at 5\u00d7 best_estimate); ",
                             "dashed = true richness")) +
      theme_minimal(base_size = 11)
  )
  ggsave(file.path(model_plot_dir, "total_expected_species_recovery.svg"),
         width = 9, height = 5, dpi = 150)
}


# =============================================================================
# SECTION 5: MODEL EXPORT ----

{
  # Model objects are LISTS keyed by fill_plants level ("FALSE"/"TRUE"), plus
  # the pooled model that carries the explicit H1 interaction term.
  saveRDS(model_final,      file.path(results_folder, "p_model_final.rds"))
  saveRDS(model_3var,       file.path(results_folder, "p_model_3var.rds"))
  saveRDS(model_ext,        file.path(results_folder, "p_model_extended.rds"))
  saveRDS(model_pooled_int, file.path(results_folder, "p_model_pooled_interaction.rds"))
  # Signed-error refits (Section 4e): DV = signed_prop_error = (est - true)/true.
  saveRDS(model_signed,     file.path(results_folder, "p_model_signed.rds"))

  # Coefficient table (built in 3d; already tagged with fill_plants level)
  write.csv(coef_table,  file.path(results_folder, "p_model_coefficients.csv"),
            row.names = FALSE)

  # Scaling parameters (mean + sd per predictor used in predict_missing_prop),
  # one block per fill_plants level.
  scale_params <- do.call(rbind, lapply(fp_levels, function(fp)
    cbind(fill_plants = fp, scale_params_ext[[fp]])))
  write.csv(scale_params, file.path(results_folder, "p_model_scale_params.csv"),
            row.names = FALSE)

  # Training data with predictions and residuals from the FINAL model of each
  # level (floor-experiment rows with fill_plants = NA are not scored).
  df_analysis_pred <- do.call(rbind, lapply(fp_levels, function(fp) {
    d <- df_ext[[fp]]
    d$model_fill_plants     <- fp
    d$pred_abs_prop_error   <- as.numeric(fitted(model_final[[fp]]))
    d$resid_abs_prop_error  <- as.numeric(residuals(model_final[[fp]]))
    d
  }))
  fwrite(df_analysis_pred, file.path(results_folder, "p_df_analysis_predictions.csv"))

  message("Section 5 export complete → ", results_folder)
}


# =============================================================================
# SECTION 6: NOTES AND CAVEATS ----
# =============================================================================
{
  # --- On the training data ---
  # • Model trained on simulated data; extrapolation to real amazonian animal
  #   data assumes the simulation parameter space spans the real dataset's properties.
  #   Before applying, verify that real-data IV values fall within the training range.
  #
  #
  # --- On applying the model to counts ---
  # • The model predicts abs_prop_error (a proportion), not a species count.
  #   To estimate the number of missing species:
  #     n_missing_est ≈ round(predicted_abs_prop_error × observed_richness)
  #   This approximation assumes the error is symmetric around the estimate;
  #   check bias_prop to confirm direction before applying.
  #
}
