# Plan: Improve 4-IDParamSens.R Analysis Pipeline

## Diagnosis

The current file has several issues that will mislead conclusions:

### 1. Univariate effect sizes are confounded
`calculate_effect_sizes` tests each IV against each DV independently. With strong collinearity (e.g., `n_grid_cells` ↔ `sample_mean_range_size` at r=1.0, `n_species` ↔ `total_observations` at r=0.80), a univariate effect size for one IV absorbs the effect of correlated IVs. You can't tell which IV actually matters.

### 2. "Direction" doesn't mean "most accurate"
For `accuracyMagn` (= estimate − truth), the most accurate level is the one closest to **0**, not the highest. For `accuracyProp` (= estimate / truth), the most accurate level is closest to **1.0**. The current code reports the level with the highest group mean, which is wrong for both DVs.

### 3. Unbalanced design biases eta-squared
`sampling_evenness` has 112 rows at 1.0 vs. 6 at 0.1. Standard Type I SS eta-squared is biased with unbalanced data. Some IVs (like `prop_sites`) have zero variance and should be dropped automatically.

### 4. No bias correction model
The user's Goal B ("how much should I correct the estimate?") requires a calibration model, not just effect sizes. Mean bias per IV level tells you the expected correction factor.

### 5. Small n (107 usable rows) with many IVs (16 treatments)
Fitting 16+ univariate models doesn't help prioritize. A multivariate approach with variable importance is needed, but must be done carefully given the sample size.

---

## Proposed Changes

### Step 1: Data prep improvements
- **Drop zero-variance columns** automatically (e.g., `prop_sites`).
- **Flag and handle near-collinear IVs** — drop one of each near-collinear pair (r > 0.95), or combine them. Specifically: `sample_mean_range_size` and `sample_sd_range_size` are redundant with `n_grid_cells`.
- **Create proper accuracy DVs**:
  - `bias_magn` = `best_estimate - trueSpeciesRichness` (signed, for correction)
  - `abs_error` = `abs(bias_magn)` (unsigned, for "how accurate")
  - `bias_prop` = `best_estimate / trueSpeciesRichness` (signed ratio, for correction)
  - `abs_prop_error` = `abs(bias_prop - 1)` (unsigned, for "how accurate")
- **Investigate the 47 NAs** (~30% of data) — are they model failures? If systematic, they carry information.

### Step 2: Replace `calculate_effect_sizes` with two purpose-built analyses

#### 2a. `rank_variable_importance()` — answers Goal A
- Fit a single multiple regression of each DV on all (non-collinear) IVs.
- Report **partial eta-squared** (Type III SS via `car::Anova`) or standardized coefficients.
- This controls for confounding and handles imbalance properly.
- Also compute **permutation-based variable importance** from a random forest (`ranger`) as a complementary nonlinear check.
- Output: ranked table of IVs by importance, with p-values and partial effect sizes.

#### 2b. `summarize_accuracy_by_level()` — answers Goal A (detail) and Goal B
For each IV (categorical or binned-continuous):
- Report **n** per level (so the user sees imbalance).
- Report **mean bias** (`bias_magn` and `bias_prop`) per level — this IS the correction factor.
- Report **mean absolute error** per level — this IS accuracy.
- Report **SD** and **CI** of bias per level.
- For continuous IVs, bin into quantile-based groups (e.g., terciles) so levels are balanced.
- "Most accurate level" = level with lowest mean absolute error.

### Step 3: Visualization upgrades
- **Variable importance plot**: horizontal bar chart of partial eta-squared or VI scores, colored by whether it's a treatment.
- **Bias correction plot**: faceted dotplot showing mean bias ± CI per level of each important IV.
- **Calibration curve**: predicted vs. actual species richness, with a 1:1 reference line.

### Step 4: Interaction-ready structure (scaffold, not full implementation)
- Add a `test_interactions()` function that takes the top-k important IVs from Step 2a and fits interaction terms.
- Use `car::Anova(..., type = "III")` for proper SS with interactions and imbalance.
- Flag significant interactions for follow-up.

---

## Implementation Order

1. Data prep (Step 1) — clean up df, drop/flag problematic columns, create DVs
2. `summarize_accuracy_by_level()` (Step 2b) — most immediately useful, answers both goals
3. `rank_variable_importance()` (Step 2a) — multivariate importance ranking
4. Visualizations (Step 3)
5. Interaction scaffold (Step 4) — stub function for future use

## Packages needed
- `car` (Type III ANOVA)
- `ranger` (random forest variable importance, optional)
- Already loaded: `data.table`, `ggplot2`

## Notes
- Keep `calculate_effect_sizes` in the file for reference but mark it as superseded.
- All new functions will follow the AGENTS.md conventions (verbose, roxygen2, base R style, no pipes).
