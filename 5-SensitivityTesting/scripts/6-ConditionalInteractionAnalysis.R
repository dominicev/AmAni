# =============================================================================
# 6-ConditionalInteractionAnalysis.R
# -----------------------------------------------------------------------------
# PURPOSE
#   Explore *why* some IVs show weak or inconsistent main effects on model
#   accuracy / bias (as catalogued in Script 5). Specifically, tests whether
#   the benefit of treatment parameters (range filling, plant filtering) is
#   conditional on the sampling context — e.g., does range filling help more
#   when sites are strongly sampled?  Does combined fill + filter help more
#   when there are more total grid cells?
#
# INPUTS  (exported by Script 5, Step 7)
#   results/p_df_analysis.csv
#   results/p_vi_magn.csv, p_vi_prop.csv
#   results/p_lm_all.csv, p_sig_counts.csv
#   results/p_interaction_results.csv
#
# OUTPUTS
#   results/p_conditional_screen.csv   — pairwise interaction screen table
#   results/p_conditional_summary.csv  — significant conditional effects (FDR-corrected)
#   results/p_conditional_plots/       — one plot per tested interaction
#
# FLOW
#   1. Setup & load
#   2. Overview: which IVs are weakly explained in Script 5?
#   3. Hypothesis inventory (commit before looking at results)
#   4. Pairwise interaction screen (data-driven)
#   5. Formal hypothesis tests (one subsection per hypothesis)
#   6. Multi-correction summary
#   7. Export findings for Notebook 7
# =============================================================================


# SETUP ----
{
  # Working-directory resolution (mirrors Script 5 convention)
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
  library(purrr)
  library(car)      # for Type III ANOVA in formal tests

  # Paths
  results_folder <- "./5-SensitivityTesting/results/"
  plot_folder    <- file.path(results_folder, "p_conditional_plots")
  dir.create(plot_folder, showWarnings = FALSE, recursive = TRUE)

  # Load exports from Script 5
  df_analysis        <- fread(file.path(results_folder, "p_df_analysis.csv"))
  vi_magn            <- read.csv(file.path(results_folder, "p_vi_magn.csv"))
  vi_prop            <- read.csv(file.path(results_folder, "p_vi_prop.csv"))
  lm_all             <- read.csv(file.path(results_folder, "p_lm_all.csv"))
  sig_counts         <- read.csv(file.path(results_folder, "p_sig_counts.csv"))
  interaction_res_5  <- read.csv(file.path(results_folder, "p_interaction_results.csv"))
  # ^ top-6-IV interaction results already computed in Script 5 Step 4;
  #   used here as a starting reference, not re-run.

  # IV role classification (mirrors Script 5 definitions)
  # treatment_ivs: parameters we directly control in the model pipeline
  treatment_ivs <- c("fill_plants", "fill_animals", "filter_plant_completeness",
                     "completeness_floor", "fit_params", "fill_interaction",
                     "n_cells_filtered")

  # sampling_ivs: properties of the sample / simulation that we cannot control.
  # Split by type because Section 3 handles categorical and continuous
  # moderators differently.
  sampling_ivs_cat <- c(
    "clustering_label"          # spatial clustering of start sites (4 levels)
  )
  sampling_ivs_cont <- c(
    "n_grid_cells", "sample_numberStartSites",
    "comp_mean", "comp_sd", "comp_spread", "sampling_evenness",
    "n_species", "total_observations",
    "cellsPerSpp", "SppPerStartCell", "prop_species", "sample_diffusion",
    # range-size structure of the sample
    "sample_mean_range_size", "sample_median_range_size", "sample_sd_range_size",
    # underlying range-size distribution shape (log-normal parameters)
    "delta_meanlog", "ratio_sdlog"
  )
  sampling_ivs <- c(sampling_ivs_cat, sampling_ivs_cont)

  # Deliberately EXCLUDED sampling candidates (verified against df_analysis):
  #   • prop_sites, sample_endemism_min, sample_endemism_max — constant (1 value),
  #     so there is no variance to relate to anything.
  #   • mean_range_size, mean_rangeSize, median_range_size, sd_range_size —
  #     entirely NA (unused placeholder columns; the real range-size signal lives
  #     in the sample_*_range_size columns above).
  #   • n_grid_cells_before/after, min_completeness — only populated for the
  #     filter experiments (>75% NA), so they enter conditional models unevenly;
  #     revisit per-experiment if needed.
  # NOTE: sample_diffusion is retained but is NA for ~37% of rows; models using
  #   it will drop those rows, so interpret its interactions on the reduced set.

  # Guard: keep only IVs that exist AND have real variance (>1 non-NA value).
  # This self-cleans the lists if the upstream export schema changes.
  has_variance <- function(x) {
    v <- df_analysis[[x]]
    !is.null(v) && sum(!is.na(v)) > 0 && dplyr::n_distinct(v, na.rm = TRUE) > 1
  }
  all_declared <- c(treatment_ivs, sampling_ivs)
  dropped_ivs  <- setdiff(all_declared, Filter(has_variance, all_declared))
  if (length(dropped_ivs)) {
    message("Dropped IVs (missing or no variance): ",
            paste(dropped_ivs, collapse = ", "))
  }
  treatment_ivs     <- Filter(has_variance, treatment_ivs)
  sampling_ivs_cat  <- Filter(has_variance, sampling_ivs_cat)
  sampling_ivs_cont <- Filter(has_variance, sampling_ivs_cont)
  sampling_ivs      <- c(sampling_ivs_cat, sampling_ivs_cont)

  # trueSpeciesRichness is NOT a moderator (it is constant in this design), but
  # it is the scaling constant that converts the proportional DVs back to real
  # units (species counts): e.g. abs_prop_error * true_richness ≈ absolute error
  # in number of species. Keep it available for reporting/back-transformation.
  if ("trueSpeciesRichness" %in% names(df_analysis)) {
    true_richness <- unique(df_analysis$trueSpeciesRichness)
    if (length(true_richness) > 1) {
      # Design changed to vary richness: keep the per-row column for scaling.
      message("trueSpeciesRichness varies (", length(true_richness),
              " values); scale per-row rather than by a single constant.")
    } else {
      message("trueSpeciesRichness (scaling constant) = ", round(true_richness, 3))
    }
  }

  # Primary DVs for this notebook
  dv_acc  <- "abs_prop_error"
  dv_bias <- "bias_prop"
}


# =============================================================================
# SETUP (cont.): SUPERVISED (PLS) COMBINATIONS FOR COLLINEAR IV GROUPS ----
#
# Many sampling IVs are strongly collinear (verified via the correlation matrix
# on the non-floor rows: e.g. n_grid_cells / n_species / total_observations all
# r > 0.8; the three range-size measures r 0.74-0.97; sample_diffusion /
# sample_numberStartSites / SppPerStartCell form a clustering cluster). Instead
# of unsupervised PCA (which captures the largest X-variance direction and can
# dilute an outcome-relevant variable), we summarize each cluster with PARTIAL
# LEAST SQUARES components oriented toward TOTAL SIGNED ACCURACY (bias_prop, the
# signed proportional error; negative = underestimate). PLS finds the linear
# combination of a group that best covaries with that outcome.
#   • pls_filter  : completeness_floor + n_cells_filtered     (TREATMENT side; floor experiment)
#   • pls_range   : sample_{mean,median,sd}_range_size + cellsPerSpp
#   • pls_cluster : sample_diffusion + sample_numberStartSites + SppPerStartCell
#   • pls_comp    : comp_mean + comp_sd + comp_spread + prop_species + sampling_evenness
#   • pls_extent  : n_grid_cells + n_species + total_observations
# delta_meanlog / ratio_sdlog (near-independent) and clustering_label
# (categorical) remain individual IVs.
#
# ** LEAKAGE WARNING ** Because these components are built to covary with
# bias_prop, naively screening them against bias_prop would be circular and
# optimistic. To keep the comparison fair we generate CROSS-FITTED scores: each
# row's score comes from a PLS model trained on OTHER rows (5-fold), with fold
# signs aligned to a full-data reference. Even so, the PLS axes retain a built-in
# advantage toward bias_prop relative to unsupervised axes; interpret their
# bias_prop results cautiously and lean on the abs_prop_error (non-target)
# column for a cleaner read.
{
  library(pls)
  pls_target <- dv_bias  # "bias_prop" — total signed accuracy

  pls_groups <- list(
    pls_filter  = list(vars = c("completeness_floor", "n_cells_filtered"),
                       keep = 1, role = "treatment"),
    pls_range   = list(vars = c("sample_mean_range_size", "sample_median_range_size",
                                "sample_sd_range_size", "cellsPerSpp"),
                       keep = 2, role = "sampling"),
    pls_cluster = list(vars = c("sample_diffusion", "sample_numberStartSites",
                                "SppPerStartCell"),
                       keep = 1, role = "sampling"),
    pls_comp    = list(vars = c("comp_mean", "comp_sd", "comp_spread",
                                "prop_species", "sampling_evenness"),
                       keep = 2, role = "sampling"),
    pls_extent  = list(vars = c("n_grid_cells", "n_species", "total_observations"),
                       keep = 1, role = "sampling")
  )

  # Cross-fitted PLS scores: k-fold OOF scores, sign-aligned to a full-data
  # reference projection, oriented positively toward the target for readability.
  # Returns NA-padded scores aligned to every row of `data`.
  build_pls <- function(data, vars, target, keep, prefix, k = 5, seed = 4127) {
    cc  <- stats::complete.cases(data[, c(vars, target), drop = FALSE])
    idx <- which(cc)
    Xy  <- data[idx, c(vars, target), drop = FALSE]
    names(Xy)[ncol(Xy)] <- "y"
    fml <- reformulate(vars, "y")

    ref      <- plsr(fml, data = Xy, ncomp = keep, scale = TRUE)
    ref_proj <- ref$projection

    set.seed(seed)
    folds <- sample(rep_len(seq_len(k), length(idx)))
    oof   <- matrix(NA_real_, length(idx), keep)
    for (f in seq_len(k)) {
      tr <- folds != f; te <- folds == f
      if (sum(tr) < length(vars) + 2 || !any(te)) next
      fit <- plsr(fml, data = Xy[tr, , drop = FALSE], ncomp = keep, scale = TRUE)
      sc  <- predict(fit, newdata = Xy[te, vars, drop = FALSE],
                     type = "scores")[, seq_len(keep), drop = FALSE]
      for (kk in seq_len(keep)) {
        if (sum(fit$projection[, kk] * ref_proj[, kk]) < 0) sc[, kk] <- -sc[, kk]
      }
      oof[te, ] <- sc
    }
    for (kk in seq_len(keep)) {
      if (stats::cor(oof[, kk], Xy$y, use = "complete.obs") < 0) oof[, kk] <- -oof[, kk]
    }
    # Honest (cross-fitted) fraction of target variance the component captures.
    r2_cv <- apply(oof, 2, function(s)
      stats::cor(s, Xy$y, use = "complete.obs")^2)
    message("  [", prefix, "] n=", length(idx), "  cross-fitted R\u00b2 toward ",
            target, ": ", paste0("C", seq_len(keep), "=",
                                 round(100 * r2_cv, 1), "%", collapse = ", "))
    out <- matrix(NA_real_, nrow(data), keep,
                  dimnames = list(NULL, paste0(prefix, seq_len(keep))))
    out[idx, ] <- oof
    as.data.frame(out)
  }

  dfp           <- as.data.frame(df_analysis)
  pls_treatment <- character(0)
  pls_sampling  <- character(0)
  for (g in names(pls_groups)) {
    grp    <- pls_groups[[g]]
    scores <- build_pls(dfp, grp$vars, pls_target, grp$keep, g)
    for (nm in names(scores)) {
      df_analysis[[nm]] <- scores[[nm]]
      if (grp$role == "treatment") pls_treatment <- c(pls_treatment, nm)
      else                         pls_sampling  <- c(pls_sampling, nm)
    }
  }

  # Register PLS axes so Section 3 screens them alongside the raw variables.
  treatment_ivs     <- c(treatment_ivs, pls_treatment)
  sampling_ivs_cont <- c(sampling_ivs_cont, pls_sampling)
  sampling_ivs      <- c(sampling_ivs_cat, sampling_ivs_cont)
  # Treatments whose values exist only in the floor experiment (used by
  # balance_for_focal to isolate that branch instead of pinning to baseline).
  floor_treatments  <- c("completeness_floor", "n_cells_filtered", pls_treatment)
  message("Added ", length(c(pls_treatment, pls_sampling)), " cross-fitted PLS axes (",
          length(pls_treatment), " treatment, ", length(pls_sampling), " sampling).")
}


# =============================================================================
# SETUP (cont.): FINAL MINIMAL VARIABLE SET FOR THE SCREEN ----
#
# RESULT OF THE DIMENSION-REDUCTION EXPERIMENTS. We tried two ways to fix 
# collinear sampling IVs: unsupervised PCA (earlier) and the supervised,
# cross-fitted PLS toward bias_prop retained above. NEITHER improved fit.
# The composite axes never beat the single best individual variable — e.g. the
# comp composite reached partial eta-sq ~0.11 for fill_plants on abs_prop_error,
# versus ~0.20 for prop_species alone (PCA's version was weaker still, ~0.07),
# and no composite survived FDR. The PLS/PCA code is kept above for the record,
# but the composites are EXCLUDED from the final screen.
#
# Instead we minimize the test count directly: keep the smallest defensible set
# of non-redundant variables — one representative per collinear cluster plus the
# (near-)independent informative ones. This shrinks the FDR family from 384
# pairs to 48, which controls false positives far better than screening every
# correlated variable.
#
# Dropped, and why:
#   • pls_* / PCA composites       — did not outperform individual variables.
#   • fill_interaction             — redundant recoding of fill_plants x fill_animals.
#   • filter_plant_completeness    — its experiment held the sampling context
#                                    fixed, so its interactions are structurally
#                                    inestimable (aliased; 6 NA pairs earlier).
#   • n_cells_filtered             — redundant with completeness_floor.
#   • n_species / total_observations / cellsPerSpp / SppPerStartCell /
#     sample_numberStartSites      — collinear (r > 0.8) with n_grid_cells (rep).
#   • sample_median/sd_range_size  — collinear with sample_mean_range_size (rep).
#   • comp_sd / sampling_evenness /
#     sample_diffusion             — collinear with comp_spread (rep).
#   • delta_meanlog / clustering_label — no screen signal; dropped to limit tests.
{
  treatment_ivs <- c(
    "fill_plants",         # H1 focal
    "fill_animals",        # H4 (null) focal
    "fit_params",          # fitting parameterisation
    "completeness_floor"   # floor-experiment knob (continuous)
  )
  sampling_ivs_cat  <- character(0)                # clustering_label dropped
  sampling_ivs_cont <- c(
    "prop_species",            # species representation (strong, distinct)
    "comp_mean",               # mean site completeness (H1/H4)
    "comp_spread",             # completeness heterogeneity / evenness cluster rep
    "n_grid_cells",            # sampling extent / intensity cluster rep (H2)
    "sample_mean_range_size",  # range-size cluster rep
    "ratio_sdlog"              # range-shape parameter (independent; earlier hit)
  )
  sampling_ivs <- c(sampling_ivs_cat, sampling_ivs_cont)

  message("Final screen set: ", length(treatment_ivs), " treatments x ",
          length(sampling_ivs), " sampling IVs x 2 DVs = ",
          length(treatment_ivs) * length(sampling_ivs) * 2, " tests.")
}


# =============================================================================
# SECTION 1: OVERVIEW — WEAK MAIN EFFECTS FROM SCRIPT 5 ----

# Before testing interactions, establish which IVs and experiments have large
# amounts of unexplained variance. High residual variance in Script 5 is the
# *motivation* for looking for conditional effects here.

{
  # Thresholds used to flag "weak" explanatory power. Kept as named objects so
  # the same cutoffs are reused in 1a/1b and are easy to revise in one place.
  weak_r2_cut  <- 0.20   # experiments below this have mostly unexplained variance
  weak_eta_cut <- 0.05   # IVs below this explain little variance
  sig_cut      <- 0.05   # p-value threshold for "non-significant"

  # --- 1a. R² distribution across experiments ---
  # adj. R² is a property of each (Experiment, DV) model, so it is duplicated
  # across the IV rows in lm_all; collapse to one row per model before plotting.
  # Focus on the proportional-error DVs, which are the outcomes carried forward
  # in this notebook (abs_prop_error for accuracy, bias_prop for bias).
  r2_by_exp <- lm_all |>
    filter(DV %in% c("abs_prop_error", "bias_prop")) |>
    distinct(Experiment, DV, R2_adj) |>
    arrange(DV, R2_adj)

  weak_experiments <- r2_by_exp |>
    filter(R2_adj < weak_r2_cut)

  print(
    ggplot(r2_by_exp,
           aes(x = R2_adj, y = reorder(Experiment, R2_adj))) +
      geom_col() +
      geom_vline(xintercept = weak_r2_cut, linetype = "dashed") +
      facet_wrap(~ DV) +
      labs(x = "Adjusted R\u00b2", y = "Experiment",
           title = "Explained variance per experiment (Script 5 main-effect models)",
           subtitle = paste0("Dashed line = weak-fit threshold (R\u00b2 < ",
                             weak_r2_cut, ")"))
  )

  # --- 1b. IVs with consistently LOW importance ---
  # An IV is a "null-main-effect" candidate if, for its proportional DV, it has
  # small partial eta-sq AND is non-significant. Check the magnitude table
  # (bias_magn / abs_error) and the proportional table (bias_prop / abs_prop_error)
  # and keep IVs that look weak in BOTH.
  weak_in <- function(vi_tbl) {
    vi_tbl |>
      filter(!is.na(Partial_Eta_Sq), !is.na(P_Value)) |>
      group_by(IV) |>
      summarise(
        max_eta = max(Partial_Eta_Sq),
        min_p   = min(P_Value),
        .groups = "drop"
      ) |>
      filter(max_eta < weak_eta_cut, min_p > sig_cut) |>
      pull(IV)
  }

  weak_ivs <- intersect(weak_in(vi_magn), weak_in(vi_prop))

  message("IVs with weak main effects in BOTH magnitude and proportional models:")
  print(weak_ivs)

  # --- 1c. Interaction results already in hand (Script 5 Step 4) ---
  # These top-6-IV interactions were computed in Script 5. They provide a
  # reference for which pairings already show signal (and which do not),
  # motivating the new hypotheses in Section 2.
  interaction_res_5 |>
    arrange(P_Value)
}


# =============================================================================
# SECTION 2: HYPOTHESIS INVENTORY ----

# Directional hypotheses for the treatment × sampling interactions.
# These were seeded by Script 5 (bias_magn interactions, Section 1c) and then
# refined by a PRELIMINARY within-group correlation look at abs_prop_error for
# the four focal pairs. Because that preliminary look has occurred, the formal
# tests in Section 4 are confirmatory-leaning-exploratory: report effect sizes
# and CIs, and do not treat any single p-value as pre-registered confirmation.
#
# Each hypothesis specifies the treatment IV, sampling IV, predicted direction,
# and primary DV.

{
  # Preliminary observations (within-level cor(sampling_iv, abs_prop_error)):
  #
  #   H1  fill_plants × prop_species:
  #       error falls with prop_species in BOTH groups, but more steeply when
  #       plants are NOT filled (r = -0.77 vs -0.49, n = 256 vs 64). Filling
  #       FLATTENS the prop_species dependence rather than helping most at high
  #       prop_species — so the interaction is real but its direction is the
  #       reverse of the naive "helps more when high" framing. RESTATED below.
  #
  #   H2  fill_interaction × n_grid_cells:
  #       slope of error on n_grid_cells is negative only when PLANTS are filled
  #       (TF r = -0.15, TT r = -0.16) and positive otherwise (FF/FT r ≈ +0.12).
  #       The moderator is plant filling, not the TT combination specifically —
  #       consistent with H1. NOTE: 84 rows have an empty fill_interaction value
  #       and must be resolved before formal testing (Section 4b).
  #
  #   H3  filter_plant_completeness × comp_spread:
  #       error RISES with comp_spread in all groups (r ≈ 0.36–0.48); filtering
  #       barely steepens it vs "none". Weak signal. DATA FLAG: the three
  #       filter levels return identical r (0.477), suggesting near-collinear /
  #       duplicated subsets — verify before interpreting.
  #
  #   H4  fill_animals × comp_mean:
  #       within-group correlations are essentially identical (r = -0.287 vs
  #       -0.284) → no interaction. Confirms the Script 5 null.
  #
  # Revised hypotheses:
  #
  # H1 (RESTATED): the negative prop_species→abs_prop_error slope is SHALLOWER
  #     when fill_plants = TRUE. Interpretation: plant filling substitutes for
  #     community completeness — it helps most when prop_species is LOW, and its
  #     marginal benefit shrinks when the community is already well represented.
  #     (Prior: strong; direction flipped by preliminary look.)
  #
  # H2 (REFRAMED): abs_prop_error DECREASES with n_grid_cells when plants are
  #     filled and INCREASES otherwise. Test on fill_plants (or fill_interaction
  #     collapsed to plant status) rather than the TT cell alone.
  #     (Prior: moderate; blocked by the empty-level data issue.)
  #
  # H3 (DOWNGRADED): filter_plant_completeness × comp_spread is weak and
  #     confounded by the identical-correlation data flag. Treat as exploratory;
  #     resolve the subset issue before formal testing.
  #
  # H4 (NULL, retained): fill_animals benefit does NOT depend on comp_mean.
  #     Documented as a pre-specified null.
  #
  # Register hypotheses in a machine-readable table so Sections 4-5 can iterate
  # over them and track which survive FDR correction.
  #   direction = expected sign of the sampling_iv→error slope WITHIN the
  #     treated group relative to the untreated group:
  #     "attenuate" = treatment makes the slope shallower (toward 0);
  #     "steepen"   = treatment makes the slope more negative;
  #     "none"      = no expected difference.
  #   status = pre-registration standing after the preliminary look.
  hypotheses <- tibble::tribble(
    ~id,   ~treatment_iv,               ~sampling_iv,   ~dv,               ~direction,   ~prior,      ~status,
    "H1",  "fill_plants",               "prop_species", "abs_prop_error",  "attenuate",  "strong",    "restated",
    "H2",  "fill_plants",               "n_grid_cells", "abs_prop_error",  "steepen",    "moderate",  "reframed",
    "H3",  "filter_plant_completeness", "comp_spread",  "abs_prop_error",  "steepen",    "weak",      "downgraded",
    "H4",  "fill_animals",              "comp_mean",    "abs_prop_error",  "none",       "weak",      "null"
  )

  # DATA-QUALITY items to resolve before Section 4:
  #   • fill_interaction has 84 rows with an empty ("") level.
  #   • filter_plant_completeness levels 0.5/0.75/0.9 show identical within-group
  #     correlations — check for duplicated or collinear rows.
  hypotheses
}


# =============================================================================
# SECTION 3: PAIRWISE INTERACTION SCREEN ----
#
# Data-driven sweep over all (treatment_iv × sampling_iv) combinations to
# DISCOVER candidate interactions worth formal testing — not to confirm them.
#
# DESIGN NOTE (correction to the original outline): the outline used two
# different effect sizes — |Δr| for categorical×continuous pairs and partial
# eta-sq for continuous×continuous pairs — but fed both into one heatmap. Those
# scales are not comparable. Here the PRIMARY, common effect size for every
# pair is the partial eta-squared of the interaction term from an OLS model
# (3c), which is defined identically regardless of variable type. The
# within-group correlation comparison (3b) is retained as an INTERPRETABILITY
# companion for categorical-treatment pairs: it supplies the direction/sign
# (does the treatment steepen or flatten the sampling→error slope?) that a
# positive-only eta-sq cannot.
#
# Every pair is screened on a BALANCED subset (holding the other treatment axes
# at baseline) so a treatment×sampling interaction is not confounded by other
# treatments. FDR-correct before interpreting any single result.


## 3a. IV roles, balancing, and helper functions ----
{
  scr <- as.data.frame(df_analysis)

  # --- Baseline cell + balancing (mirrors Script 5's balance_focal) ---
  # OFAT design: each new-param axis was varied while the others were held at
  # this baseline. To compare a focal treatment fairly we pin the OTHER axes.
  new_param_baseline <- list(
    fill_plants               = "FALSE",
    fill_animals              = "FALSE",
    filter_plant_completeness = "none",
    fit_params                = "amplitude, meanlog, sdlog",
    completeness_floor        = NA
  )

  # balance_for_focal(): like Script 5's balance_focal for the standard OFAT
  # axes, with a special case for the completeness-floor experiment. That
  # experiment (rows with a non-NA completeness_floor) is a SEPARATE branch:
  # its fill_* and fit_params are NA by design, so it cannot be pinned to the
  # baseline cell. For those focal treatments we simply isolate the floor rows.
  balance_for_focal <- function(data, focal) {
    floor_focal <- if (exists("floor_treatments")) floor_treatments else
      c("completeness_floor", "n_cells_filtered")
    if (focal %in% floor_focal) {
      return(data[!is.na(data[[focal]]), , drop = FALSE])
    }
    focal_axes <- if (identical(focal, "fill_interaction")) {
      c("fill_plants", "fill_animals")
    } else focal
    hold_axes <- setdiff(names(new_param_baseline), focal_axes)
    keep <- rep(TRUE, nrow(data))
    for (ax in hold_axes) keep <- keep & data[[ax]] %in% new_param_baseline[[ax]]
    data[keep, , drop = FALSE]
  }

  # --- IV role classification for the screen ---
  is_cat <- function(x) is.character(scr[[x]]) || is.logical(scr[[x]]) ||
    is.factor(scr[[x]])
  treat_cat  <- Filter(is_cat, treatment_ivs)   # fill_plants, fill_animals,
                                                #   filter_plant_completeness,
                                                #   fit_params, fill_interaction
  treat_cont <- setdiff(treatment_ivs, treat_cat)  # completeness_floor, n_cells_filtered
  screen_dvs <- c(dv_acc, dv_bias)

  # --- Fisher-z helpers for the correlation companion (3b) ---
  # 95% CI for one correlation.
  cor_ci <- function(r, n) {
    if (n < 4 || is.na(r)) return(c(NA_real_, NA_real_))
    se <- 1 / sqrt(n - 3)
    tanh(atanh(r) + c(-1, 1) * qnorm(0.975) * se)
  }
  # Two-sided test that two independent correlations are equal.
  cor_diff_test <- function(r1, n1, r2, n2) {
    if (any(c(n1, n2) < 4) || any(is.na(c(r1, r2)))) return(NA_real_)
    z <- (atanh(r1) - atanh(r2)) / sqrt(1 / (n1 - 3) + 1 / (n2 - 3))
    2 * pnorm(-abs(z))
  }

  # --- Core screen: partial eta-sq of the interaction term (3c) ---
  # Uses extra-sum-of-squares: SS_int = RSS(add) - RSS(full);
  # partial eta-sq = SS_int / (SS_int + RSS_full) = 1 - RSS_full / RSS_add.
  # Invariant to centring/scaling of the predictors.
  screen_pair <- function(data, treat, samp, dv, min_n = 20, min_cell = 5) {
    dd <- data[stats::complete.cases(data[, c(treat, samp, dv)]),
               c(treat, samp, dv), drop = FALSE]
    names(dd) <- c("T", "S", "Y")
    if (is.character(dd$T)) dd <- dd[dd$T != "", , drop = FALSE]
    if (is.character(dd$S)) dd <- dd[dd$S != "", , drop = FALSE]
    if (nrow(dd) < min_n) return(NULL)
    if (dplyr::n_distinct(dd$T) < 2 || dplyr::n_distinct(dd$S) < 2) return(NULL)
    # Require a minimum cell size for categorical levels so lm is estimable.
    if (is.character(dd$T) && min(table(dd$T)) < min_cell) return(NULL)
    if (is.character(dd$S) && min(table(dd$S)) < min_cell) return(NULL)
    if (is.numeric(dd$T)) dd$T <- as.numeric(scale(dd$T))
    if (is.numeric(dd$S)) dd$S <- as.numeric(scale(dd$S))

    out <- tryCatch({
      m_add  <- lm(Y ~ T + S, data = dd)
      m_full <- lm(Y ~ T * S, data = dd)
      cmp <- anova(m_add, m_full)
      rss_add  <- cmp$RSS[1]
      rss_full <- cmp$RSS[2]
      eta <- (rss_add - rss_full) / rss_add
      tibble::tibble(
        treatment_iv   = treat, sampling_iv = samp, dv = dv,
        n = nrow(dd), df_int = cmp$Df[2],
        partial_eta_sq = eta, F_value = cmp$F[2], p_raw = cmp$`Pr(>F)`[2]
      )
    }, error = function(e) NULL)
    out
  }
}

## 3b. Within-group correlation comparison (categorical treatment × continuous sampling) ----
# Interpretability companion: within each treatment level, how does the
# sampling IV correlate with the DV? A large difference between levels means the
# treatment moderates the sampling→error relationship, and its SIGN tells us the
# direction. Continuous sampling IVs only (correlation is undefined for the
# categorical clustering_label).
{
  cont_samp <- sampling_ivs_cont
  grid_b <- expand.grid(treatment_iv = treat_cat, sampling_iv = cont_samp,
                        dv = screen_dvs, stringsAsFactors = FALSE)

  cor_diff_results <- purrr::pmap_dfr(grid_b, function(treatment_iv, sampling_iv, dv) {
    bd <- balance_for_focal(scr, treatment_iv)
    dd <- bd[stats::complete.cases(bd[, c(treatment_iv, sampling_iv, dv)]),
             c(treatment_iv, sampling_iv, dv), drop = FALSE]
    names(dd) <- c("lev", "S", "Y")
    dd <- dd[!(is.character(dd$lev) & dd$lev == ""), , drop = FALSE]
    if (nrow(dd) < 20) return(NULL)

    lev_tbl <- dd |>
      dplyr::group_by(lev) |>
      dplyr::filter(dplyr::n() >= 10, dplyr::n_distinct(S) > 1) |>
      dplyr::summarise(n = dplyr::n(), r = cor(S, Y), .groups = "drop")
    if (nrow(lev_tbl) < 2) return(NULL)

    hi <- which.max(lev_tbl$r)
    lo <- which.min(lev_tbl$r)
    tibble::tibble(
      treatment_iv = treatment_iv, sampling_iv = sampling_iv, dv = dv,
      n_levels  = nrow(lev_tbl),
      level_hi  = as.character(lev_tbl$lev[hi]), r_hi = lev_tbl$r[hi],
      level_lo  = as.character(lev_tbl$lev[lo]), r_lo = lev_tbl$r[lo],
      delta_r   = lev_tbl$r[hi] - lev_tbl$r[lo],
      p_cor_diff = cor_diff_test(lev_tbl$r[hi], lev_tbl$n[hi],
                                 lev_tbl$r[lo], lev_tbl$n[lo])
    )
  })

  # NOTE: for treatments with >2 levels this compares the two EXTREME levels
  # chosen post-hoc, so p_cor_diff is anti-conservative and is descriptive only.
  # The formal ranking uses partial_eta_sq from 3c (all levels, one test).
  cor_diff_results |> dplyr::arrange(dplyr::desc(abs(delta_r)))
}

## 3c. Unified interaction screen: partial eta-squared for all treatment × sampling pairs ----
# Primary, type-agnostic effect size. Covers categorical and continuous
# treatments crossed with categorical (clustering_label) and continuous
# sampling IVs.
{
  grid_c <- expand.grid(treatment_iv = treatment_ivs, sampling_iv = sampling_ivs,
                        dv = screen_dvs, stringsAsFactors = FALSE)

  screen_results <- purrr::pmap_dfr(grid_c, function(treatment_iv, sampling_iv, dv) {
    bd <- balance_for_focal(scr, treatment_iv)
    screen_pair(bd, treatment_iv, sampling_iv, dv)
  })

  # Annotate variable types and attach the directional companion from 3b.
  screen_results <- screen_results |>
    dplyr::mutate(
      treat_type = ifelse(treatment_iv %in% treat_cat, "categorical", "continuous"),
      samp_type  = ifelse(sampling_iv %in% sampling_ivs_cat, "categorical", "continuous")
    ) |>
    dplyr::left_join(
      dplyr::select(cor_diff_results, treatment_iv, sampling_iv, dv,
                    delta_r, p_cor_diff),
      by = c("treatment_iv", "sampling_iv", "dv")
    )

  message("Screened ", nrow(screen_results), " estimable pairs of ",
          nrow(grid_c), " requested.")
  screen_results |> dplyr::arrange(dplyr::desc(partial_eta_sq)) |> head(10)
}

## 3d. Screen results: FDR correction, ranked table, and heatmap ----
{
  # BH correction across the whole screen (one family of tests).
  screen_results <- screen_results |>
    dplyr::mutate(
      p_fdr      = p.adjust(p_raw, method = "BH"),
      noteworthy = !is.na(p_fdr) & p_fdr < 0.10
    )

  # Tag pairs that were NOT pre-registered in Section 2 as exploratory. A pair
  # matches a hypothesis if its (treatment_iv, sampling_iv) appears in the table.
  hyp_pairs <- paste(hypotheses$treatment_iv, hypotheses$sampling_iv, sep = " x ")
  screen_results <- screen_results |>
    dplyr::mutate(
      pair          = paste(treatment_iv, sampling_iv, sep = " x "),
      preregistered = pair %in% hyp_pairs,
      finding = dplyr::case_when(
        !noteworthy                ~ "ns",
        noteworthy & preregistered ~ "confirmatory",
        TRUE                       ~ "exploratory"
      )
    )

  # Ranked table (top 20 by effect size).
  screen_ranked <- screen_results |>
    dplyr::arrange(dplyr::desc(partial_eta_sq)) |>
    dplyr::select(treatment_iv, sampling_iv, dv, n, partial_eta_sq,
                  delta_r, p_raw, p_fdr, finding) |>
    head(20)
  print(screen_ranked)

  message("Noteworthy pairs (p_fdr < 0.10): ", sum(screen_results$noteworthy),
          " (", sum(screen_results$finding == "confirmatory"), " confirmatory, ",
          sum(screen_results$finding == "exploratory"), " exploratory); ",
          sum(is.na(screen_results$p_raw)), " pairs inestimable (NA p).")

  # Heatmap: treatment (rows) × sampling (cols), fill = partial eta-sq,
  # star marks p_fdr < 0.10, faceted by DV.
  print(
    ggplot(screen_results,
           aes(x = sampling_iv, y = treatment_iv, fill = partial_eta_sq)) +
      geom_tile(color = "grey90") +
      geom_point(data = subset(screen_results, noteworthy),
                 shape = 8, size = 2.2, color = "cyan3", stroke = 0.9) +
      facet_wrap(~ dv, ncol = 1) +
      scale_fill_viridis_c(option = "magma", direction = -1,
                           name = expression(partial~eta^2)) +
      labs(
        x = "Sampling IV", y = "Treatment IV",
        title = "Conditional interaction screen (balanced subsets)",
        subtitle = "Star = BH-FDR < 0.10; effect size = partial eta-squared of the interaction"
      ) +
      theme_minimal(base_size = 10) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  )

  # Save outputs.
  write.csv(
    dplyr::select(screen_results, treatment_iv, sampling_iv, dv, treat_type,
                  samp_type, n, df_int, partial_eta_sq, F_value, delta_r,
                  p_raw, p_fdr, noteworthy, preregistered, finding),
    file.path(results_folder, "p_conditional_screen.csv"),
    row.names = FALSE
  )
  ggsave(
    file.path(plot_folder, "screen_heatmap.png"),
    width = 10, height = 8, dpi = 150
  )
  message("Section 3 screen saved → ",
          file.path(results_folder, "p_conditional_screen.csv"))
}


# =============================================================================
# SECTION 4: FORMAL HYPOTHESIS TESTS ----

# One subsection per interaction that the Section 3 screen actually flagged as
# worth a formal look. The screen (48-pair BH-FDR family) returned exactly ONE
# pair surviving correction — fill_plants × prop_species on abs_prop_error, the
# preregistered H1 — plus one non-surviving-but-suggestive candidate,
# fill_plants × sample_mean_range_size. Everything else in the family was null
# (p_fdr ≈ 1), including the other Section 2 hypotheses:
#   • H2 fill_plants × n_grid_cells   → partial eta-sq 0.015-0.037, p_fdr ≈ 0.80-1
#   • H3 filter_plant_completeness    → inestimable/aliased (dropped from screen)
#   • H4 fill_animals × comp_mean     → partial eta-sq 0.005, p_fdr = 1 (null holds)
# Those are therefore reported as NOT tested further here (their formal result is
# the screen row); we do not manufacture confirmatory tests for null pairs.
#
# Structure per interaction:
#   i.   Balance the dataset (balance_for_focal, Section 3a).
#   ii.  Fit the interaction model.
#   iii. Report: interaction coefficient + 95% CI, Type III partial eta-sq, p.
#   iv.  Visualize (saved to plot_folder).
#   v.   Plain-language interpretation.

## 4-setup. Shared helpers + results collector ----
{
  # Partial eta-squared for one term from a Type III ANOVA table (car::Anova).
  partial_eta_typeIII <- function(aov_tbl, term) {
    ss_term <- aov_tbl[term, "Sum Sq"]
    ss_res  <- aov_tbl["Residuals", "Sum Sq"]
    ss_term / (ss_term + ss_res)
  }

  # Pull the screen's FDR verdict for a pair/DV so each formal test can state
  # whether it survived multiplicity correction (single source of truth).
  screen_lookup <- function(treat, samp, dv) {
    screen_results |>
      dplyr::filter(treatment_iv == treat, sampling_iv == samp, dv == !!dv) |>
      dplyr::select(partial_eta_sq, p_raw, p_fdr, noteworthy, finding)
  }

  # Collector: one row per formal test, fed to Sections 5-6.
  formal_tests <- list()
}

## 4a. H1 (CONFIRMATORY) — fill_plants × prop_species → abs_prop_error ----
# The only pair surviving FDR in the screen (p_fdr = 0.015, partial eta-sq 0.20).
# Preliminary look (Section 2) predicted "attenuate": plant filling should make
# the negative prop_species→error slope SHALLOWER.
{
  bd <- balance_for_focal(scr, "fill_plants")
  dd <- bd[stats::complete.cases(bd[, c("fill_plants", "prop_species", dv_acc)]),
           c("fill_plants", "prop_species", dv_acc)]
  names(dd)[3] <- "Y"
  dd$fill_plants <- factor(dd$fill_plants, levels = c("FALSE", "TRUE"))

  m_h1  <- lm(Y ~ fill_plants * prop_species, data = dd)
  a_h1  <- car::Anova(m_h1, type = 3)
  int_term <- "fill_plantsTRUE:prop_species"

  # Report
  ci_h1  <- confint(m_h1)[int_term, ]
  eta_h1 <- partial_eta_typeIII(a_h1, "fill_plants:prop_species")
  p_h1   <- a_h1["fill_plants:prop_species", "Pr(>F)"]

  message("H1  fill_plants x prop_species -> abs_prop_error  (n = ", nrow(dd), ")")
  message("  interaction beta = ", round(coef(m_h1)[int_term], 3),
          "  95% CI [", round(ci_h1[1], 3), ", ", round(ci_h1[2], 3), "]")
  message("  Type III partial eta-sq = ", round(eta_h1, 3), "   p = ", signif(p_h1, 3))
  print(screen_lookup("fill_plants", "prop_species", dv_acc))

  # Within-group slopes (direction of the moderation).
  slopes_h1 <- dd |>
    dplyr::group_by(fill_plants) |>
    dplyr::summarise(n = dplyr::n(),
                     slope = coef(lm(Y ~ prop_species))[2],
                     r     = cor(prop_species, Y), .groups = "drop")
  print(as.data.frame(slopes_h1))

  # Visualise
  print(
    ggplot(dd, aes(prop_species, Y, color = fill_plants)) +
      geom_point(alpha = 0.7) +
      geom_smooth(method = "lm", se = TRUE) +
      scale_color_manual(values = c("FALSE" = "grey40", "TRUE" = "#1b9e77"),
                         name = "fill_plants") +
      labs(x = "prop_species (community representation)",
           y = "abs_prop_error",
           title = "H1: plant filling flattens the prop_species -> error slope",
           subtitle = paste0("Interaction beta = ", round(coef(m_h1)[int_term], 3),
                             " [", round(ci_h1[1], 2), ", ", round(ci_h1[2], 2),
                             "], partial eta-sq = ", round(eta_h1, 2),
                             ", p = ", signif(p_h1, 2), " (survives FDR)")) +
      theme_minimal(base_size = 11)
  )
  ggsave(file.path(plot_folder, "H1_fill_plants_x_prop_species.png"),
         width = 8, height = 5, dpi = 150)

  formal_tests[["H1"]] <- data.frame(
    hypothesis = "H1", treatment_iv = "fill_plants", sampling_iv = "prop_species",
    dv = dv_acc, term = int_term,
    beta = coef(m_h1)[int_term], ci_lo = ci_h1[1], ci_hi = ci_h1[2],
    partial_eta_sq = eta_h1, p_raw = p_h1,
    finding = "confirmatory", stringsAsFactors = FALSE
  )
}

## 4b. (EXPLORATORY) — fill_plants × sample_mean_range_size ----
# Next-strongest screen pair, but did NOT survive FDR (p_fdr = 0.41 for both
# DVs). Nominal Type III p ≈ 0.02-0.03; partial eta-sq ≈ 0.08. Reported for
# completeness and flagged EXPLORATORY — do not treat as a confirmed effect.
{
  bd <- balance_for_focal(scr, "fill_plants")
  for (dv in c(dv_acc, dv_bias)) {
    dd <- bd[stats::complete.cases(bd[, c("fill_plants", "sample_mean_range_size", dv)]),
             c("fill_plants", "sample_mean_range_size", dv)]
    names(dd)[3] <- "Y"
    dd$fill_plants <- factor(dd$fill_plants, levels = c("FALSE", "TRUE"))

    m  <- lm(Y ~ fill_plants * sample_mean_range_size, data = dd)
    a  <- car::Anova(m, type = 3)
    it <- "fill_plantsTRUE:sample_mean_range_size"
    ci  <- confint(m)[it, ]
    eta <- partial_eta_typeIII(a, "fill_plants:sample_mean_range_size")
    p   <- a["fill_plants:sample_mean_range_size", "Pr(>F)"]

    message("EXPLORATORY  fill_plants x sample_mean_range_size -> ", dv,
            "  (n = ", nrow(dd), ")")
    message("  interaction beta = ", round(coef(m)[it], 3),
            "  95% CI [", round(ci[1], 3), ", ", round(ci[2], 3), "]",
            "  partial eta-sq = ", round(eta, 3), "  p = ", signif(p, 3),
            "  (does NOT survive FDR)")

    formal_tests[[paste0("E_", dv)]] <- data.frame(
      hypothesis = "exploratory", treatment_iv = "fill_plants",
      sampling_iv = "sample_mean_range_size", dv = dv, term = it,
      beta = coef(m)[it], ci_lo = ci[1], ci_hi = ci[2],
      partial_eta_sq = eta, p_raw = p,
      finding = "exploratory", stringsAsFactors = FALSE
    )
  }

  # Visualise both DVs together for inspection.
  ex_long <- bd |>
    dplyr::select(fill_plants, sample_mean_range_size,
                  dplyr::all_of(c(dv_acc, dv_bias))) |>
    tidyr::pivot_longer(dplyr::all_of(c(dv_acc, dv_bias)),
                        names_to = "dv", values_to = "Y") |>
    dplyr::filter(!is.na(sample_mean_range_size), !is.na(Y))
  ex_long$fill_plants <- factor(ex_long$fill_plants, levels = c("FALSE", "TRUE"))

  print(
    ggplot(ex_long, aes(sample_mean_range_size, Y, color = fill_plants)) +
      geom_point(alpha = 0.6) +
      geom_smooth(method = "lm", se = TRUE) +
      facet_wrap(~ dv, scales = "free_y") +
      scale_color_manual(values = c("FALSE" = "grey40", "TRUE" = "#1b9e77"),
                         name = "fill_plants") +
      labs(x = "sample_mean_range_size", y = "error metric",
           title = "Exploratory: fill_plants x sample_mean_range_size",
           subtitle = "Nominally significant but does NOT survive FDR correction") +
      theme_minimal(base_size = 11)
  )
  ggsave(file.path(plot_folder, "EXP_fill_plants_x_range_size.png"),
         width = 9, height = 4.5, dpi = 150)
}

## 4c. Prespecified nulls / inestimable pairs — screen is the final word ----
{
  # H2, H3, H4 were not upheld by the screen; their screen rows ARE their formal
  # results. Print them together so the record is explicit.
  null_pairs <- dplyr::bind_rows(
    screen_lookup("fill_plants", "n_grid_cells", dv_acc) |>
      dplyr::mutate(hypothesis = "H2", pair = "fill_plants x n_grid_cells"),
    screen_lookup("fill_animals", "comp_mean", dv_acc) |>
      dplyr::mutate(hypothesis = "H4", pair = "fill_animals x comp_mean")
  )
  message("Prespecified pairs not upheld by the screen (reported as null):")
  print(as.data.frame(null_pairs))
  message("H3 (filter_plant_completeness x comp_spread) was aliased/inestimable ",
          "in this OFAT design and was excluded from the screen; not tested.")
}


# =============================================================================
# SECTION 5: MULTIPLE-COMPARISON SUMMARY ----
#
# Compile every formal Section 4 test into one table and BH-correct across it.
# CAVEAT: this is a small, secondary correction. H1's primary multiplicity
# control is the 48-pair screen family in Section 3 (where it already survived,
# p_fdr = 0.015); the two exploratory tests share H1's balanced subset (n = 64),
# so they are not independent looks. Read this table as bookkeeping, not as a
# second, stronger gate.

{
  formal_summary <- dplyr::bind_rows(formal_tests) |>
    dplyr::mutate(
      p_fdr    = p.adjust(p_raw, method = "BH"),
      survives = p_fdr < 0.10
    ) |>
    dplyr::arrange(p_raw) |>
    dplyr::select(hypothesis, treatment_iv, sampling_iv, dv, term,
                  beta, ci_lo, ci_hi, partial_eta_sq, p_raw, p_fdr,
                  finding, survives)
  rownames(formal_summary) <- NULL
  print(formal_summary)

  n_prereg   <- sum(formal_summary$finding == "confirmatory")
  n_explore  <- sum(formal_summary$finding == "exploratory")
  message("Formal tests: ", nrow(formal_summary), " total (",
          n_prereg, " confirmatory, ", n_explore, " exploratory).")
  message("Survive BH-FDR < 0.10 within Section 4: ",
          sum(formal_summary$survives),
          " (confirmatory: ",
          sum(formal_summary$survives & formal_summary$finding == "confirmatory"),
          ", exploratory: ",
          sum(formal_summary$survives & formal_summary$finding == "exploratory"), ").")
}


# =============================================================================
# SECTION 6: FINDINGS SUMMARY — EXPORT FOR NOTEBOOK 7 ----
#
# Turn the formal results into the interaction-term candidate pool for
# Notebook 7. Only H1 is carried as a confirmed interaction; the exploratory
# range-size pair is exported flagged so Notebook 7 can decide whether to
# include it.

{
  interaction_summary <- formal_summary |>
    dplyr::transmute(
      term          = gsub("TRUE", "", term),   # e.g. "fill_plants:prop_species"
      type          = "interaction",
      treatment_iv, sampling_iv, dv,
      partial_eta_sq, p_fdr,
      finding,
      direction = dplyr::case_when(
        hypothesis == "H1" ~
          "plant filling flattens the negative prop_species->error slope; benefit largest when prop_species is LOW",
        beta < 0 ~ paste0("fill_plants makes the ", sampling_iv,
                          "->", dv, " slope more negative"),
        TRUE     ~ paste0("fill_plants makes the ", sampling_iv,
                          "->", dv, " slope more positive")
      )
    ) |>
    dplyr::arrange(dplyr::desc(partial_eta_sq))
  print(interaction_summary)

  write.csv(interaction_summary,
            file.path(results_folder, "p_interaction_summary.csv"),
            row.names = FALSE)
  message("Section 6 export complete (", nrow(interaction_summary),
          " interaction terms) → ",
          file.path(results_folder, "p_interaction_summary.csv"))
}
