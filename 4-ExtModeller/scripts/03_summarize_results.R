# =====================================================================
# 03_summarize_results.R — Compile and compare all simulation results
# =====================================================================

{
  library(knitr)
  library(rstudioapi)
  if (exists(".amani_wd")) {
    setwd(.amani_wd)
  } else {
    project_root <- dirname(dirname(rstudioapi::getActiveDocumentContext()$path))
    knitr::opts_knit$set(root.dir = project_root)
    setwd(project_root)
  }

  library(dplyr)
  library(ggplot2)
}

# Output directory (overridden by manager when running organism-specific pipelines)
if (exists(".amani_output_dir")) {
  outdir <- .amani_output_dir
} else {
  outdir <- "output"
}

############################
# LOAD ALL REPLICATE CSVs
############################

neutral_reps  <- read.csv(paste0(outdir, "/neutral_all_replicates.csv"))
clustered_reps <- read.csv(paste0(outdir, "/clustered_all_replicates.csv"))
shifted_reps  <- read.csv(paste0(outdir, "/shifted_all_replicates.csv"))

neutral_reps$model  <- "Neutral"
clustered_reps$model <- "Clustered"
shifted_reps$model  <- "Shifted"

all_reps <- rbind(neutral_reps, clustered_reps, shifted_reps)
all_reps$model <- factor(all_reps$model, levels = c("Neutral", "Shifted", "Clustered"))

cat("Loaded", nrow(all_reps), "total replicates across",
    length(unique(all_reps$model)), "models and",
    length(unique(all_reps$scenario)), "scenarios\n")


############################
# SUMMARY TABLE
############################

summary_tbl <- all_reps |>
  group_by(model, scenario) |>
  summarise(
    median_spp_lost = median(n_spp_lost_total),
    mean_spp_lost   = round(mean(n_spp_lost_total), 1),
    sd_spp_lost     = round(sd(n_spp_lost_total), 1),
    q25_spp_lost    = quantile(n_spp_lost_total, 0.25),
    q75_spp_lost    = quantile(n_spp_lost_total, 0.75),
    min_spp_lost    = min(n_spp_lost_total),
    max_spp_lost    = max(n_spp_lost_total),
    median_pop_lost = median(n_pop_lost_total),
    mean_pop_lost   = round(mean(n_pop_lost_total), 1),
    min_pop_lost    = min(n_pop_lost_total),
    max_pop_lost    = max(n_pop_lost_total),
    .groups = "drop"
  )

cat("\n=== FULL SUMMARY TABLE ===\n")
summary_tbl

write.csv(summary_tbl, paste0(outdir, "/summary_all_models.csv"), row.names = FALSE)
cat("\nSaved output/summary_all_models.csv\n")


############################
# SPECIES EXTINCTION PLOT
############################

timestamp_str <- format(Sys.time(), "%Y%m%d_%H%M%S")

p_spp <- ggplot(all_reps, aes(x = model, y = n_spp_lost_total, fill = model)) +
  geom_boxplot(outlier.size = 0.8) +
  stat_summary(fun.min = min, fun.max = max, geom = "errorbar",
               width = 0.3, linewidth = 0.4) +
  facet_wrap(~scenario, scales = "free_y") +
  scale_fill_manual(values = c(
    "Neutral"   = "#E69F00",
    "Shifted"   = "#56B4E9",
    "Clustered" = "#009E73"
  )) +
  labs(
    x = NULL, y = "Species lost",
    title = "Species extinctions across deforestation models"
  ) +
  theme_minimal() +
  theme(legend.position = "none")
p_spp

ggsave(paste0(outdir, "/species_extinctions_", timestamp_str, ".pdf"),
       p_spp, width = 8, height = 5)
cat("Saved species extinction plot to PDF\n")


############################
# POPULATION LOSS PLOT
############################

p_pop <- ggplot(all_reps, aes(x = model, y = n_pop_lost_total, fill = model)) +
  geom_boxplot(outlier.size = 0.8) +
  stat_summary(fun.min = min, fun.max = max, geom = "errorbar",
               width = 0.3, linewidth = 0.4) +
  facet_wrap(~scenario, scales = "free_y") +
  scale_fill_manual(values = c(
    "Neutral"   = "#E69F00",
    "Shifted"   = "#56B4E9",
    "Clustered" = "#009E73"
  )) +
  labs(
    x = NULL, y = "Populations lost",
    title = "Population losses across deforestation models"
  ) +
  theme_minimal() +
  theme(legend.position = "none")
p_pop

ggsave(paste0(outdir, "/population_losses_", timestamp_str, ".pdf"),
       p_pop, width = 8, height = 5)
cat("Saved population loss plot to PDF\n")


############################
# SPECIES EXTINCTION STATUS COMPARISON
############################

neutral_spp  <- read.csv(paste0(outdir, "/neutral_species_extinction_status.csv"))
clustered_spp <- read.csv(paste0(outdir, "/clustered_species_extinction_status.csv"))

# compare single-run extinction lists
ext_comparison <- data.frame(
  scenario = c("FreqOver", "FreqBoundary", "FreqUnder"),
  n_extinct_neutral = c(
    sum(neutral_spp$extinct_over),
    sum(neutral_spp$extinct_boundary),
    sum(neutral_spp$extinct_under)
  ),
  n_extinct_clustered = c(
    sum(clustered_spp$extinct_over),
    sum(clustered_spp$extinct_boundary),
    sum(clustered_spp$extinct_under)
  ),
  pct_single_grid_neutral = round(c(
    100 * mean(neutral_spp$n_grids[neutral_spp$extinct_over] == 1, na.rm = TRUE),
    100 * mean(neutral_spp$n_grids[neutral_spp$extinct_boundary] == 1, na.rm = TRUE),
    100 * mean(neutral_spp$n_grids[neutral_spp$extinct_under] == 1, na.rm = TRUE)
  ), 1),
  pct_single_grid_clustered = round(c(
    100 * mean(clustered_spp$n_grids[clustered_spp$extinct_over] == 1, na.rm = TRUE),
    100 * mean(clustered_spp$n_grids[clustered_spp$extinct_boundary] == 1, na.rm = TRUE),
    100 * mean(clustered_spp$n_grids[clustered_spp$extinct_under] == 1, na.rm = TRUE)
  ), 1)
)

cat("\n=== EXTINCTION STATUS COMPARISON (single-run) ===\n")
ext_comparison

write.csv(ext_comparison, paste0(outdir, "/summary_extinction_status_comparison.csv"), row.names = FALSE)
cat("Saved output/summary_extinction_status_comparison.csv\n")


############################
# RANGE SIZE DISTRIBUTION OF EXTINCT SPECIES
############################

# combine both status files with a model label
neutral_spp$model <- "Neutral"
clustered_spp$model <- "Clustered"

# melt to long form for plotting
neutral_long <- neutral_spp |>
  tidyr::pivot_longer(
    cols = c(extinct_over, extinct_boundary, extinct_under),
    names_to = "scenario", values_to = "extinct"
  ) |>
  mutate(scenario = sub("extinct_", "", scenario))

clustered_long <- clustered_spp |>
  tidyr::pivot_longer(
    cols = c(extinct_over, extinct_boundary, extinct_under),
    names_to = "scenario", values_to = "extinct"
  ) |>
  mutate(scenario = sub("extinct_", "", scenario))

both_long <- rbind(neutral_long, clustered_long)

extinct_only <- both_long |> filter(extinct)

p_range <- ggplot(extinct_only, aes(x = n_grids, fill = model)) +
  geom_histogram(binwidth = 1, position = "dodge") +
  facet_wrap(~scenario, scales = "free_y") +
  scale_fill_manual(values = c("Neutral" = "#E69F00", "Clustered" = "#009E73")) +
  labs(
    x = "Range size (number of grid cells)",
    y = "Number of extinct species",
    title = "Range sizes of extinct species: Neutral vs Clustered"
  ) +
  theme_minimal()
p_range

ggsave(paste0(outdir, "/extinct_range_sizes_", timestamp_str, ".pdf"),
       p_range, width = 8, height = 5)
cat("Saved range size plot to PDF\n")


############################
# RANGE REDUCTION THRESHOLDS (from replicates)
############################
# The replicate CSVs now include n_reduced_20..n_reduced_100 columns
# (species at each range-reduction threshold). Summarize across replicates.

threshold_cols <- c("n_reduced_20", "n_reduced_50", "n_reduced_75",
                    "n_reduced_90", "n_reduced_95", "n_reduced_99",
                    "n_reduced_100")

if (all(threshold_cols %in% names(all_reps))) {

  threshold_summary <- all_reps |>
    group_by(model, scenario) |>
    summarise(
      across(all_of(threshold_cols),
             list(median = median, q25 = ~quantile(.x, 0.25), q75 = ~quantile(.x, 0.75)),
             .names = "{.col}_{.fn}"),
      .groups = "drop"
    )

  write.csv(threshold_summary,
            paste0(outdir, "/summary_range_reduction_thresholds.csv"),
            row.names = FALSE)
  cat("\nSaved output/summary_range_reduction_thresholds.csv\n")

  cat("\n=== RANGE REDUCTION THRESHOLD SUMMARY (medians) ===\n")
  threshold_summary |>
    dplyr::select(model, scenario, ends_with("_median")) |>
    as.data.frame() |>
    print(row.names = FALSE)

} else {
  cat("\nRange reduction columns not found in replicate data.\n")
  cat("Re-run extinction models to generate threshold counts.\n")
}


############################
# PER-SPECIES RANGE REDUCTION (from single detailed runs)
############################

neutral_rr_file   <- paste0(outdir, "/neutral_species_range_reduction.csv")
clustered_rr_file <- paste0(outdir, "/clustered_species_range_reduction.csv")

if (file.exists(neutral_rr_file) && file.exists(clustered_rr_file)) {

  neutral_rr   <- read.csv(neutral_rr_file)
  clustered_rr <- read.csv(clustered_rr_file)

  neutral_rr$model   <- "Neutral"
  clustered_rr$model <- "Clustered"

  all_rr <- rbind(neutral_rr, clustered_rr)

  # Summary table: count species at each threshold per model × scenario
  rr_summary <- all_rr |>
    group_by(model, scenario) |>
    summarise(
      n_species       = n(),
      reduced_ge_20   = sum(pct_range_lost >= 0.20),
      reduced_ge_50   = sum(pct_range_lost >= 0.50),
      reduced_ge_75   = sum(pct_range_lost >= 0.75),
      reduced_ge_90   = sum(pct_range_lost >= 0.90),
      reduced_ge_95   = sum(pct_range_lost >= 0.95),
      reduced_ge_99   = sum(pct_range_lost >= 0.99),
      extinct_100     = sum(pct_range_lost >= 1.00),
      .groups = "drop"
    )

  write.csv(rr_summary,
            paste0(outdir, "/summary_range_reduction_single_run.csv"),
            row.names = FALSE)
  cat("\nSaved output/summary_range_reduction_single_run.csv\n")

  cat("\n=== RANGE REDUCTION SUMMARY (single detailed runs) ===\n")
  rr_summary |> as.data.frame() |> print(row.names = FALSE)

  # Save the combined per-species detail
  write.csv(all_rr, paste0(outdir, "/combined_species_range_reduction.csv"),
            row.names = FALSE)
  cat("Saved output/combined_species_range_reduction.csv\n")

} else {
  cat("\nPer-species range reduction CSVs not found. Run the detailed simulations first.\n")
}
