# AmAni.summarize.R — Cross-organism summary of extinction model results
#
# Provides two entry-point functions:
#   summarize_all()        — cross-organism summary (table + figures)
#   summarize_organism()   — single-organism compilation (table + figure)
#
# Reads summary CSVs, gamma richness estimates, and saved plot data from each
# organism's output directory and produces:
#   1. A combined summary table (CSV)
#   2. Comparison figures across organisms (PDF)
#   3. Compilation figures regenerated from saved plot data (PDF)
#
# Designed to be sourced at the end of AmAni.run.R or run standalone.

library(tidyverse)
library(sf)
library(patchwork)
library(ggrepel)

# --- Organism display labels --------------------------------------------------
# Maps internal directory names to human-readable labels for plots.
organism_label <- function(x) {
  dplyr::case_when(
    x == "cockroaches" ~ "Cockroaches.Curated",
    x == "allRoaches"  ~ "cockroaches",
    TRUE ~ paste0(toupper(substr(x, 1, 1)), substr(x, 2, nchar(x)))
  )
}

# Predefined organism groups for easy filtering in summarize_all().
# Usage: summarize_all(group = "insects")
ORGANISM_GROUPS <- list(
  insects     = c("cockroaches", "allRoaches", "beetles", "lepidoptera",
                  "hemiptera", "odonata", "butterflies", "hymenoptera",
                  "bees", "flies", "dungbeetles"),
  vertebrates = c("antbirds", "finches", "frogs", "mammals", "tyrannidae",
                  "squamates"),
  birds       = c("antbirds", "finches", "tyrannidae"),
  herps       = c("frogs", "squamates")
)

summarize_all <- function(include = NULL, group = NULL) {

if (!exists("project_root")) source("AmAni.manager.R")

# --- Discover which organisms have completed results ------------------------

output_dir <- file.path(project_root, "output")
organism_dirs_all <- list.dirs(output_dir, full.names = FALSE, recursive = FALSE)

# Only include organisms that have a summary_all_models.csv
organisms <- organism_dirs_all[sapply(organism_dirs_all, function(org) {
  file.exists(file.path(output_dir, org, "extinction_model", "summary_all_models.csv"))
})]

# --- Apply organism subset filter (include / group) ---
if (!is.null(group)) {
  if (!group %in% names(ORGANISM_GROUPS))
    stop("Unknown group '", group, "'. Available: ",
         paste(names(ORGANISM_GROUPS), collapse = ", "), call. = FALSE)
  include <- union(include, ORGANISM_GROUPS[[group]])
}
if (!is.null(include)) {
  organisms <- organisms[organisms %in% include]
}

if (length(organisms) == 0) {
  message("No completed organism results found. Skipping summary.")
} else {

message("Summarizing results for: ", paste(organisms, collapse = ", "))

summary_dir <- file.path(output_dir, "cross_organism_summary")
dir.create(summary_dir, showWarnings = FALSE, recursive = TRUE)

timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")


# --- 1. Collect extinction model summaries ----------------------------------

extinction_summary <- map_dfr(organisms, function(org) {
  read.csv(file.path(output_dir, org, "extinction_model", "summary_all_models.csv")) |>
    mutate(organism = org, .before = 1)
})

extinction_summary$organism <- organism_label(extinction_summary$organism)

write.csv(extinction_summary,
          file.path(summary_dir, "combined_extinction_summary.csv"),
          row.names = FALSE)
message("  Wrote: combined_extinction_summary.csv")


# --- 2. Collect gamma richness estimates ------------------------------------

gamma_summary <- map_dfr(organisms, function(org) {
  gamma_files <- list.files(
    file.path(output_dir, org, "animal_model"),
    pattern = "^gamma_richnessfromSPRICH.*\\.txt$",
    full.names = TRUE
  )
  if (length(gamma_files) == 0) {
    return(tibble(organism = org, gamma_richness = NA_real_,
                  upper_bound = NA_real_, lower_bound = NA_real_))
  }
  gamma_file <- gamma_files[which.max(file.mtime(gamma_files))]
  lines <- readLines(gamma_file)

  parse_val <- function(pattern) {
    line <- grep(pattern, lines, value = TRUE, ignore.case = TRUE)
    if (length(line) == 0) return(NA_real_)
    as.numeric(gsub("[^0-9.]", "", sub(".*:\\s*", "", line[1])))
  }

  tibble(
    organism = org,
    gamma_richness = parse_val("^Best estimate|^gammaRichness"),
    upper_bound = parse_val("^Upper bound"),
    lower_bound = parse_val("^Lower bound")
  )
})

gamma_summary$organism <- organism_label(gamma_summary$organism)

write.csv(gamma_summary,
          file.path(summary_dir, "combined_gamma_richness.csv"),
          row.names = FALSE)
message("  Wrote: combined_gamma_richness.csv")


# --- 3. Extinction comparison figures ---------------------------------------

# 3a. Median species lost by organism and spatial model (FreqBoundary only)
plot_data <- extinction_summary |>
  filter(scenario == "FreqBoundary")

if (nrow(plot_data) > 0) {
  pdf(file.path(summary_dir, paste0("species_lost_by_organism_", timestamp, ".pdf")),
      width = 10, height = 6)

  p <- ggplot(plot_data, aes(x = organism, y = median_spp_lost, fill = model)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7) +
    geom_errorbar(aes(ymin = q25_spp_lost, ymax = q75_spp_lost),
                  position = position_dodge(width = 0.8), width = 0.25) +
    labs(title = "Median Species Lost Under Deforestation (FreqBoundary)",
         x = NULL, y = "Species lost (median, IQR)",
         fill = "Spatial model") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
  print(p)
  dev.off()
  message("  Wrote: species_lost_by_organism_", timestamp, ".pdf")
}

# 3b. All scenarios faceted by organism
pdf(file.path(summary_dir, paste0("species_lost_all_scenarios_", timestamp, ".pdf")),
    width = 12, height = 8)

p <- ggplot(extinction_summary,
            aes(x = scenario, y = median_spp_lost, fill = model)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_errorbar(aes(ymin = q25_spp_lost, ymax = q75_spp_lost),
                position = position_dodge(width = 0.8), width = 0.25) +
  facet_wrap(~organism, scales = "free_y") +
  labs(title = "Species Lost Under Deforestation — All Scenarios",
       x = NULL, y = "Species lost (median, IQR)",
       fill = "Spatial model") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(p)
dev.off()
message("  Wrote: species_lost_all_scenarios_", timestamp, ".pdf")

# 3c. Gamma richness comparison
if (any(!is.na(gamma_summary$gamma_richness))) {
  pdf(file.path(summary_dir, paste0("gamma_richness_comparison_", timestamp, ".pdf")),
      width = 8, height = 5)

  gamma_plot <- gamma_summary |> filter(!is.na(gamma_richness))

  p <- ggplot(gamma_plot, aes(x = organism, y = gamma_richness)) +
    geom_col(width = 0.6) +
    geom_errorbar(aes(ymin = lower_bound, ymax = upper_bound), width = 0.25) +
    labs(title = "Estimated Amazonian Gamma Richness by Organism",
         x = NULL, y = "Estimated total species") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
  print(p)
  dev.off()
  message("  Wrote: gamma_richness_comparison_", timestamp, ".pdf")
}

# 3d. Clustered vs neutral ratio
ratio_data <- extinction_summary |>
  filter(model %in% c("Neutral", "Clustered")) |>
  dplyr::select(organism, model, scenario, median_spp_lost) |>
  pivot_wider(names_from = model, values_from = median_spp_lost) |>
  mutate(reduction_pct = (1 - Clustered / Neutral) * 100)

if (nrow(ratio_data) > 0) {
  write.csv(ratio_data,
            file.path(summary_dir, "clustered_vs_neutral_reduction.csv"),
            row.names = FALSE)
  message("  Wrote: clustered_vs_neutral_reduction.csv")

  pdf(file.path(summary_dir, paste0("clustered_reduction_", timestamp, ".pdf")),
      width = 10, height = 6)

  p <- ggplot(ratio_data, aes(x = organism, y = reduction_pct, fill = scenario)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7) +
    labs(title = "Reduction in Species Loss: Clustered vs. Neutral Model",
         subtitle = "Positive = clustered model predicts fewer extinctions",
         x = NULL, y = "% reduction in species lost",
         fill = "Frequency scenario") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
  print(p)
  dev.off()
  message("  Wrote: clustered_reduction_", timestamp, ".pdf")
}


# --- 4. Compilation figures from saved plot data ----------------------------
#
# Each organism's animal_model/plot_data/ directory contains RDS files with
# the data underlying its key plots. We load them all, tag with organism name,
# and produce multi-panel compilation figures.

message("\n--- Generating compilation figures from saved plot data ---")

# Helper: load an RDS from an organism's plot_data/ dir, or return NULL
load_plot_data <- function(org, filename) {
  path <- file.path(output_dir, org, "animal_model", "plot_data", filename)
  if (file.exists(path)) readRDS(path) else NULL
}

# Discover organisms with plot_data directories
orgs_with_plots <- organism_dirs_all[sapply(organism_dirs_all, function(org) {
  dir.exists(file.path(output_dir, org, "animal_model", "plot_data"))
})]

if (length(orgs_with_plots) == 0) {
  message("  No organisms have plot_data/ directories yet. Skipping compilation figures.")
  message("  (Re-run the animal model pipeline to generate plot data.)")
} else {

message("  Organisms with plot data: ", paste(orgs_with_plots, collapse = ", "))

# ---- 4a. Animal richness maps (faceted by organism) ----

richness_maps <- map_dfr(orgs_with_plots, function(org) {
  d <- load_plot_data(org, "plotdata_animal_richness_map.rds")
  if (!is.null(d)) mutate(d, organism = organism_label(org)) else NULL
})

if (nrow(richness_maps) > 0) {
  n_orgs <- n_distinct(richness_maps$organism)
  ncol_facet <- min(n_orgs, 3)
  fig_w <- 5 * ncol_facet
  fig_h <- 5 * ceiling(n_orgs / ncol_facet)

  pdf(file.path(summary_dir, paste0("compilation_animal_richness_maps_", timestamp, ".pdf")),
      width = fig_w, height = fig_h)
  p <- ggplot(richness_maps) +
    geom_sf(aes(fill = animal_richness), color = NA) +
    scale_fill_viridis_c(option = "magma", direction = -1, name = "Species\nrichness") +
    facet_wrap(~organism, ncol = ncol_facet) +
    theme_minimal() +
    labs(title = "Animal Richness per Grid Cell")
  print(p)
  dev.off()
  message("  Wrote: compilation_animal_richness_maps_", timestamp, ".pdf")
}


# ---- 4b. Animal endemism maps (faceted by organism) ----

endemism_maps <- map_dfr(orgs_with_plots, function(org) {
  d <- load_plot_data(org, "plotdata_animal_endemism_map.rds")
  if (!is.null(d)) mutate(d, organism = organism_label(org)) else NULL
})

if (nrow(endemism_maps) > 0) {
  n_orgs <- n_distinct(endemism_maps$organism)
  ncol_facet <- min(n_orgs, 3)
  fig_w <- 5 * ncol_facet
  fig_h <- 5 * ceiling(n_orgs / ncol_facet)

  pdf(file.path(summary_dir, paste0("compilation_animal_endemism_maps_", timestamp, ".pdf")),
      width = fig_w, height = fig_h)
  p <- ggplot(endemism_maps) +
    geom_sf(aes(fill = animal_endemism_count), color = NA) +
    scale_fill_viridis_c(option = "plasma", direction = -1, name = "Endemic\nspecies") +
    facet_wrap(~organism, ncol = ncol_facet) +
    theme_minimal() +
    labs(title = "Animal Endemism per Grid Cell (single-cell endemics)")
  print(p)
  dev.off()
  message("  Wrote: compilation_animal_endemism_maps_", timestamp, ".pdf")
}


# ---- 4c. Endemism correlation: animal ~ plant (faceted) ----

endem_corr <- map_dfr(orgs_with_plots, function(org) {
  d <- load_plot_data(org, "plotdata_endemism_correlation.rds")
  if (!is.null(d)) mutate(d, organism = organism_label(org)) else NULL
})

if (nrow(endem_corr) > 0) {
  # Compute R² per organism for annotation
  endem_r2 <- endem_corr |>
    group_by(organism) |>
    summarise(
      r2 = tryCatch(summary(lm(animal_endemism_count ~ raw_endemism))$r.squared, error = function(e) NA_real_),
      n = n(),
      .groups = "drop"
    ) |>
    mutate(label = sprintf("R² = %.3f  (n = %d)", r2, n))

  pdf(file.path(summary_dir, paste0("compilation_endemism_correlation_", timestamp, ".pdf")),
      width = 12, height = 8)
  p <- ggplot(endem_corr, aes(x = raw_endemism, y = animal_endemism_count)) +
    geom_point(alpha = 0.6, color = "darkgreen") +
    geom_smooth(method = "lm", se = TRUE, color = "red", linewidth = 0.8) +
    geom_text(data = endem_r2, aes(label = label),
              x = -Inf, y = Inf, hjust = -0.1, vjust = 1.5,
              size = 3, inherit.aes = FALSE) +
    facet_wrap(~organism, scales = "free") +
    theme_minimal() +
    labs(title = "Plant vs Animal Endemism per Grid Cell",
         x = "Plant raw endemism", y = "Animal endemism count")
  print(p)
  dev.off()
  message("  Wrote: compilation_endemism_correlation_", timestamp, ".pdf")
}


# ---- 4d. Richness correlation: animal ~ plant (faceted) ----

rich_corr <- map_dfr(orgs_with_plots, function(org) {
  d <- load_plot_data(org, "plotdata_richness_correlation.rds")
  if (!is.null(d)) mutate(d, organism = organism_label(org)) else NULL
})

if (nrow(rich_corr) > 0) {
  rich_r2 <- rich_corr |>
    group_by(organism) |>
    summarise(
      r2 = tryCatch(summary(lm(morphotype_richness ~ plant_richness))$r.squared, error = function(e) NA_real_),
      n = n(),
      .groups = "drop"
    ) |>
    mutate(label = sprintf("R² = %.3f  (n = %d)", r2, n))

  pdf(file.path(summary_dir, paste0("compilation_richness_correlation_", timestamp, ".pdf")),
      width = 12, height = 8)
  p <- ggplot(rich_corr, aes(x = plant_richness, y = morphotype_richness)) +
    geom_point(alpha = 0.6, color = "darkgreen") +
    geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.8) +
    geom_text(data = rich_r2, aes(label = label),
              x = -Inf, y = Inf, hjust = -0.1, vjust = 1.5,
              size = 3, inherit.aes = FALSE) +
    facet_wrap(~organism, scales = "free") +
    theme_minimal() +
    labs(title = "Plant vs Animal Richness per Grid Cell",
         x = "Plant richness", y = "Animal morphotype richness")
  print(p)
  dev.off()
  message("  Wrote: compilation_richness_correlation_", timestamp, ".pdf")
}


# ---- 4e. Scaling factor distributions (faceted) ----

scale_dists <- map_dfr(orgs_with_plots, function(org) {
  d <- load_plot_data(org, "plotdata_scaling_factors.rds")
  if (!is.null(d)) mutate(d, organism = organism_label(org)) else NULL
})

if (nrow(scale_dists) > 0) {
  pdf(file.path(summary_dir, paste0("compilation_scaling_factors_", timestamp, ".pdf")),
      width = 12, height = 8)
  p <- ggplot(scale_dists, aes(x = scaling_factor)) +
    geom_histogram(bins = 30, fill = "darkgreen", color = "black") +
    facet_wrap(~organism, scales = "free") +
    theme_minimal() +
    labs(title = "Distribution of Endemism Scaling Factors (< 2)",
         x = "Scaling factor", y = "Frequency")
  print(p)
  dev.off()
  message("  Wrote: compilation_scaling_factors_", timestamp, ".pdf")
}


# ---- 4f. Simulated endemic species distributions (faceted) ----
# SKIPPED — The endemism-ratio simulation is a rough extrapolation kept in
# each organism's individual output but excluded from cross-organism summaries.
# See endemism_ratio_simulation_results_*.txt in each animal_model/ folder.


# ---- 4g. Plant extent model with CIs (faceted) ----

extent_obs <- map_dfr(orgs_with_plots, function(org) {
  d <- load_plot_data(org, "plotdata_plant_extent_model.rds")
  if (!is.null(d)) mutate(d$plant_df, organism = organism_label(org)) else NULL
})

extent_ci <- map_dfr(orgs_with_plots, function(org) {
  d <- load_plot_data(org, "plotdata_plant_extent_model.rds")
  if (is.null(d)) return(NULL)
  ci <- d$ci_df
  ci$y_pred <- d$y_pred
  ci$y_predUpper <- d$y_predUpper
  ci$organism <- organism_label(org)
  ci
})

if (nrow(extent_obs) > 0 && nrow(extent_ci) > 0) {
  pdf(file.path(summary_dir, paste0("compilation_plant_extent_model_", timestamp, ".pdf")),
      width = 12, height = 8)
  p <- ggplot() +
    geom_ribbon(data = extent_ci, aes(x = num_grid_cells, ymin = lower, ymax = upper),
                fill = "darkgreen", alpha = 0.2) +
    geom_line(data = extent_ci, aes(x = num_grid_cells, y = median),
              color = "darkgreen", linewidth = 0.8) +
    geom_line(data = extent_ci, aes(x = num_grid_cells, y = y_pred),
              color = "steelblue", linetype = "dashed", linewidth = 0.6) +
    geom_line(data = extent_ci, aes(x = num_grid_cells, y = y_predUpper),
              color = "steelblue", linetype = "dashed", linewidth = 0.6) +
    geom_point(data = extent_obs, aes(x = num_grid_cells, y = n),
               color = "darkgreen", size = 1.5) +
    scale_x_log10() +
    scale_y_log10() +
    facet_wrap(~organism, scales = "free") +
    theme_minimal() +
    labs(title = "Plant Spatial Extent Model with Bootstrapped CIs",
         x = "Number of grid cells (log scale)", y = "Number of plant species (log scale)")
  print(p)
  dev.off()
  message("  Wrote: compilation_plant_extent_model_", timestamp, ".pdf")
}


# ---- 4h. Scaled endemism maps (faceted by organism) ----

scaled_maps <- map_dfr(orgs_with_plots, function(org) {
  d <- load_plot_data(org, "plotdata_scaled_endemism_map.rds")
  if (!is.null(d)) mutate(d, organism = organism_label(org)) else NULL
})

if (nrow(scaled_maps) > 0) {
  n_orgs <- n_distinct(scaled_maps$organism)
  ncol_facet <- min(n_orgs, 3)
  fig_w <- 5 * ncol_facet
  fig_h <- 5 * ceiling(n_orgs / ncol_facet)

  pdf(file.path(summary_dir, paste0("compilation_scaled_endemism_maps_", timestamp, ".pdf")),
      width = fig_w, height = fig_h)
  p <- ggplot(scaled_maps) +
    geom_sf(aes(fill = scaled_plant_endemism), color = NA) +
    scale_fill_viridis_c(option = "plasma", na.value = "white", name = "Scaled\nendemism") +
    facet_wrap(~organism, ncol = ncol_facet) +
    theme_minimal() +
    labs(title = "Scaled Plant Endemism (Adjusted by Animal Offset)")
  print(p)
  dev.off()
  message("  Wrote: compilation_scaled_endemism_maps_", timestamp, ".pdf")
}


# ---- 4i. Simulated richness maps (faceted by organism) ----

sim_rich_maps <- map_dfr(orgs_with_plots, function(org) {
  d <- load_plot_data(org, "plotdata_simulated_richness_map.rds")
  if (!is.null(d)) mutate(d, organism = organism_label(org)) else NULL
})

if (nrow(sim_rich_maps) > 0) {
  n_orgs <- n_distinct(sim_rich_maps$organism)
  ncol_facet <- min(n_orgs, 3)
  fig_w <- 5 * ncol_facet
  fig_h <- 5 * ceiling(n_orgs / ncol_facet)

  pdf(file.path(summary_dir, paste0("compilation_simulated_richness_maps_", timestamp, ".pdf")),
      width = fig_w, height = fig_h)
  p <- ggplot(sim_rich_maps) +
    geom_sf(aes(fill = sim_richness), color = NA) +
    scale_fill_viridis_c(option = "magma", direction = -1, na.value = "white",
                         name = "Simulated\nrichness") +
    facet_wrap(~organism, ncol = ncol_facet) +
    theme_minimal() +
    labs(title = "Simulated Species Richness per Grid Cell")
  print(p)
  dev.off()
  message("  Wrote: compilation_simulated_richness_maps_", timestamp, ".pdf")
}


# ---- 4j. Plant endemism maps (shared — not organism-specific) ----
# These are the same across all organisms, so just use the first one found.

plant_raw <- load_plot_data(orgs_with_plots[1], "plotdata_plant_raw_endemism.rds")
plant_wt  <- load_plot_data(orgs_with_plots[1], "plotdata_plant_weighted_endemism.rds")

if (!is.null(plant_raw) && !is.null(plant_wt)) {
  pdf(file.path(summary_dir, paste0("compilation_plant_endemism_maps_", timestamp, ".pdf")),
      width = 14, height = 6)

  p1 <- ggplot(plant_raw) +
    geom_sf(aes(fill = raw_endemism), color = NA) +
    scale_fill_viridis_c(option = "plasma", na.value = "white", name = "Raw\nendemism") +
    theme_minimal() +
    labs(title = "Plant Raw Endemism")

  p2 <- ggplot(plant_wt) +
    geom_sf(aes(fill = weighted_endemism), color = NA) +
    scale_fill_viridis_c(option = "plasma", na.value = "white", name = "Weighted\nendemism") +
    theme_minimal() +
    labs(title = "Plant Weighted Endemism")

  print(p1 + p2)
  dev.off()
  message("  Wrote: compilation_plant_endemism_maps_", timestamp, ".pdf")
}

message("  Compilation figures complete.")

} # end if orgs_with_plots


# --- 5. Cross-organism dataset & model summary table -------------------------
#
# Parses the text report files from each organism's animal_model/ directory
# to compile: pre/post filter stats, scaling factor, gamma richness,
# observed richness, and correlation statistics into one table + slope graph.

message("\n--- Compiling cross-organism dataset summary ---")

# Helper: extract a value from a text file by regex pattern
extract_val <- function(lines, pattern, as_type = "numeric") {
  line <- grep(pattern, lines, value = TRUE, ignore.case = TRUE)
  if (length(line) == 0) return(NA)
  raw <- sub(".*:\\s*", "", line[1])
  raw <- trimws(raw)
  if (as_type == "numeric") {
    # Pull first number (possibly with decimal)
    m <- regmatches(raw, regexpr("[0-9]+\\.?[0-9]*", raw))
    if (length(m) == 0) return(NA_real_)
    return(as.numeric(m))
  }
  raw
}

# Helper: extract "N (M well-sampled + K super-diverse; J outliers dropped)"
extract_post_cells <- function(lines) {
  line <- grep("Grid cells:", lines, value = TRUE)
  # Pick the line from the "After filtering" section (second occurrence)
  if (length(line) < 2) return(NA_integer_)
  m <- regmatches(line[2], regexpr("^\\s*Grid cells:\\s*([0-9]+)", line[2]))
  if (length(m) == 0) return(NA_integer_)
  as.integer(sub(".*:\\s*", "", m))
}

organism_summary <- map_dfr(organism_dirs_all, function(org) {
  am_dir <- file.path(output_dir, org, "animal_model")
  if (!dir.exists(am_dir)) return(NULL)

  # --- Dataset summary (filter stats) ---
  ds_file <- list.files(am_dir, pattern = "^dataset_summary_.*\\.txt$", full.names = TRUE)
  if (length(ds_file) > 0) {
    ds_file <- ds_file[which.max(file.mtime(ds_file))]
    ds <- readLines(ds_file)

    pre_cells   <- extract_val(ds, "^\\s*Grid cells:")
    pre_records <- extract_val(ds, "^\\s*Total records:")
    pre_species <- extract_val(ds, "^\\s*Total species:")

    # After filtering section starts after "--- After filtering ---"
    after_idx <- grep("After filtering", ds)
    if (length(after_idx) > 0) {
      ds_after <- ds[after_idx[1]:length(ds)]
      post_cells   <- extract_val(ds_after, "^\\s*Grid cells:")
      post_records <- extract_val(ds_after, "^\\s*Total records:")
      post_species <- extract_val(ds_after, "^\\s*Total species:")
    } else {
      post_cells <- post_records <- post_species <- NA_real_
    }

    filter_applied <- any(grepl("Filter applied:\\s*YES", ds))
    min_occ  <- extract_val(ds, "min_occurrences:")
    min_comp <- extract_val(ds, "min_completeness:")
    min_rich <- extract_val(ds, "min_richness:")
  } else {
    # No dataset summary (e.g. cockroaches with no filter)
    pre_cells <- pre_records <- pre_species <- NA_real_
    post_cells <- post_records <- post_species <- NA_real_
    filter_applied <- FALSE
    min_occ <- min_comp <- min_rich <- NA_real_
  }

  # --- Gamma richness + scaling factor ---
  gm_file <- list.files(am_dir, pattern = "^gamma_richnessfromSPRICH.*\\.txt$", full.names = TRUE)
  if (length(gm_file) > 0) {
    gm_file <- gm_file[which.max(file.mtime(gm_file))]
    gm <- readLines(gm_file)
    scaling_factor <- extract_val(gm, "Scaling factor:")
    gamma_best     <- extract_val(gm, "Best estimate:")
    gamma_upper    <- extract_val(gm, "Upper bound:")
    gamma_lower    <- extract_val(gm, "Lower bound:")
  } else {
    scaling_factor <- gamma_best <- gamma_upper <- gamma_lower <- NA_real_
  }

  # --- Correlation summary ---
  cs_file <- list.files(am_dir, pattern = "^correlation_summary_.*\\.txt$", full.names = TRUE)
  if (length(cs_file) > 0) {
    cs_file <- cs_file[which.max(file.mtime(cs_file))]
    cs <- readLines(cs_file)

    # Animal cells used
    animal_cells <- extract_val(cs, "Animal cells used:")

    # Richness correlation
    rich_idx <- grep("Richness correlation", cs)
    if (length(rich_idx) > 0) {
      cs_rich <- cs[rich_idx[1]:min(rich_idx[1] + 10, length(cs))]
      rich_cells <- extract_val(cs_rich, "Cells:")
      rich_r2    <- extract_val(cs_rich, "R..:")
      rich_p     <- extract_val(cs_rich, "p-value:")
    } else {
      rich_cells <- rich_r2 <- rich_p <- NA_real_
    }

    # Endemism correlation
    end_idx <- grep("Endemism correlation", cs)
    if (length(end_idx) > 0) {
      cs_end <- cs[end_idx[1]:min(end_idx[1] + 10, length(cs))]
      end_r2 <- extract_val(cs_end, "R..:")
      end_p  <- extract_val(cs_end, "p-value:")
    } else {
      end_r2 <- end_p <- NA_real_
    }
  } else {
    animal_cells <- rich_cells <- rich_r2 <- rich_p <- end_r2 <- end_p <- NA_real_
  }

  # Use observed species: post-filter if filter was applied, otherwise pre-filter
  observed_spp <- if (!is.na(post_species) && filter_applied) post_species else pre_species

  tibble(
    organism         = org,
    pre_cells        = pre_cells,
    pre_records      = pre_records,
    pre_species      = pre_species,
    filter_applied   = filter_applied,
    min_occurrences  = min_occ,
    min_completeness = min_comp,
    min_richness     = min_rich,
    post_cells       = post_cells,
    post_records     = post_records,
    post_species     = post_species,
    model_cells      = animal_cells,
    scaling_factor   = scaling_factor,
    observed_spp     = observed_spp,
    gamma_best       = gamma_best,
    gamma_lower      = gamma_lower,
    gamma_upper      = gamma_upper,
    gamma_difference = gamma_best - observed_spp,
    gamma_ratio      = round(gamma_best / observed_spp, 2),
    rich_corr_cells  = rich_cells,
    rich_r2          = rich_r2,
    rich_p           = rich_p,
    end_r2           = end_r2,
    end_p            = end_p
  )
})

if (nrow(organism_summary) > 0) {

  organism_summary$organism <- organism_label(organism_summary$organism)
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  write.csv(organism_summary,
            file.path(summary_dir, paste0(paste0("combined_organism_summary.", timestamp),".csv" )),
            row.names = FALSE)
  message("  Wrote: combined_organism_summary csv file")

  # Print a compact version to console
  cat("\n=== Cross-Organism Dataset & Model Summary ===\n\n")
  organism_summary |>
    dplyr::select(organism, pre_species, post_species, model_cells,
                  scaling_factor, observed_spp, gamma_best, gamma_difference, gamma_ratio,
                  rich_r2, rich_p) |>
    as.data.frame() |>
    print(row.names = FALSE)

  # --- Slope graph: observed richness → estimated gamma richness ---
  slope_data <- organism_summary |>
    dplyr::filter(!is.na(observed_spp) & !is.na(gamma_best)) |>
    dplyr::select(organism, observed_spp, gamma_best) |>
    tidyr::pivot_longer(cols = c(observed_spp, gamma_best),
                        names_to = "metric", values_to = "species") |>
    dplyr::mutate(
      metric = factor(metric,
                      levels = c("observed_spp", "gamma_best"),
                      labels = c("Observed\n(in grid)", "Estimated\ngamma"))
    )

  if (nrow(slope_data) > 0) {
    pdf(file.path(summary_dir, paste0("slope_observed_vs_gamma_", timestamp, ".pdf")),
        width = 8, height = 6)

    p <- ggplot(slope_data, aes(x = metric, y = species, group = organism)) +
      geom_line(aes(color = organism), linewidth = 0.8, alpha = 0.7) +
      geom_point(aes(color = organism), size = 3) +
      ggrepel::geom_text_repel(
                data = slope_data |> dplyr::filter(metric == "Estimated\ngamma"),
                aes(label = organism, color = organism),
                hjust = 0, nudge_x = 0.08, direction = "y",
                size = 3, show.legend = FALSE, segment.size = 0.3,
                max.overlaps = Inf) +
      scale_y_log10(labels = scales::comma) +
      labs(title = "Observed vs Estimated Gamma Richness by Organism",
           subtitle = "Log scale; each line connects observed species in grid to estimated Amazon total",
           x = NULL, y = "Species (log scale)",
           color = "Organism") +
      theme_minimal() +
      theme(legend.position = "none",
            panel.grid.minor = element_blank())
    print(p)
    dev.off()
    message("  Wrote: slope_observed_vs_gamma_", timestamp, ".pdf")
  }

  # --- Bar chart: gamma expansion ratio (gamma / observed) ---
  ratio_plot_data <- organism_summary |>
    dplyr::filter(!is.na(gamma_ratio)) |>
    dplyr::arrange(gamma_ratio)

  if (nrow(ratio_plot_data) > 0) {
    pdf(file.path(summary_dir, paste0("gamma_expansion_ratio_", timestamp, ".pdf")),
        width = 10, height = 6)

    p <- ggplot(ratio_plot_data,
                aes(y = reorder(organism, gamma_ratio), x = gamma_ratio)) +
      geom_col(width = 0.6) +
      geom_vline(xintercept = 1, linetype = "dashed", color = "red") +
      geom_text(aes(label = paste0(gamma_ratio, "x")),
                hjust = -0.2, size = 3) +
      labs(title = "Gamma Expansion Ratio (Estimated / Observed)",
           subtitle = "Red dashed line = 1x (no change). Higher = more unsampled species predicted.",
           y = NULL, x = "Gamma / Observed species") +
      theme_minimal()
    print(p)
    dev.off()
    message("  Wrote: gamma_expansion_ratio_", timestamp, ".pdf")
  }
}


# --- Print summary to console -----------------------------------------------

message("\n=== Cross-organism Summary ===")
message("Organisms with results: ", length(organisms))
message("Output directory: ", summary_dir)

cat("\nGamma Richness Estimates:\n")
print(gamma_summary, n = Inf)

cat("\nExtinction Summary (FreqBoundary, median species lost):\n")
extinction_summary |>
  filter(scenario == "FreqBoundary") |>
  dplyr::select(organism, model, median_spp_lost, q25_spp_lost, q75_spp_lost) |>
  as.data.frame() |>
  print()

} # end if organisms found

} # end summarize_all()


# =============================================================================
# --- 6. Single-organism summary function ------------------------------------
# =============================================================================
#
# summarize_organism("cockroaches") produces:
#   - A multi-panel compilation figure (PDF) combining all key plots
#   - A summary statistics table (CSV) with the same metrics as the
#     cross-organism comparison
#
# Designed to be called interactively after sourcing this file.

summarize_organism <- function(organism, output_dir = NULL) {

  if (!exists("project_root")) stop("Source AmAni.manager.R first.", call. = FALSE)

  if (is.null(output_dir)) output_dir <- file.path(project_root, "output")

  org_am  <- file.path(output_dir, organism, "animal_model")
  org_ext <- file.path(output_dir, organism, "extinction_model")
  org_out <- file.path(output_dir, organism, "organism_summary")
  dir.create(org_out, showWarnings = FALSE, recursive = TRUE)

  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")

  load_pd <- function(filename) {
    path <- file.path(org_am, "plot_data", filename)
    if (file.exists(path)) readRDS(path) else NULL
  }

  if (!dir.exists(file.path(org_am, "plot_data"))) {
    stop("No plot_data/ directory found for '", organism,
         "'. Run the animal model pipeline first.", call. = FALSE)
  }

  message("Summarizing single organism: ", organism)

  extract_val_local <- function(lines, pattern) {
    line <- grep(pattern, lines, value = TRUE, ignore.case = TRUE)
    if (length(line) == 0) return(NA_real_)
    raw <- sub(".*:\\s*", "", line[1])
    m <- regmatches(raw, regexpr("[0-9]+\\.?[0-9]*", raw))
    if (length(m) == 0) return(NA_real_)
    as.numeric(m)
  }

  ds_file <- list.files(org_am, pattern = "^dataset_summary_.*\\.txt$", full.names = TRUE)
  if (length(ds_file) > 0) {
    ds_file <- ds_file[which.max(file.mtime(ds_file))]
    ds <- readLines(ds_file)
    pre_cells   <- extract_val_local(ds, "^\\s*Grid cells:")
    pre_records <- extract_val_local(ds, "^\\s*Total records:")
    pre_species <- extract_val_local(ds, "^\\s*Total species:")
    after_idx <- grep("After filtering", ds)
    if (length(after_idx) > 0) {
      ds_after <- ds[after_idx[1]:length(ds)]
      post_cells   <- extract_val_local(ds_after, "^\\s*Grid cells:")
      post_records <- extract_val_local(ds_after, "^\\s*Total records:")
      post_species <- extract_val_local(ds_after, "^\\s*Total species:")
    } else {
      post_cells <- post_records <- post_species <- NA_real_
    }
    filter_applied <- any(grepl("Filter applied:\\s*YES", ds))
  } else {
    pre_cells <- pre_records <- pre_species <- NA_real_
    post_cells <- post_records <- post_species <- NA_real_
    filter_applied <- FALSE
  }

  gm_file <- list.files(org_am, pattern = "^gamma_richnessfromSPRICH.*\\.txt$", full.names = TRUE)
  if (length(gm_file) > 0) {
    gm_file <- gm_file[which.max(file.mtime(gm_file))]
    gm <- readLines(gm_file)
    extent_scaling_factor <- extract_val_local(gm, "Extent scaling factor:|Scaling factor:")
    gamma_best     <- extract_val_local(gm, "Best estimate:")
    gamma_upper    <- extract_val_local(gm, "Upper bound:")
    gamma_lower    <- extract_val_local(gm, "Lower bound:")
  } else {
    extent_scaling_factor <- gamma_best <- gamma_upper <- gamma_lower <- NA_real_
  }

  cs_file <- list.files(org_am, pattern = "^correlation_summary_.*\\.txt$", full.names = TRUE)
  if (length(cs_file) > 0) {
    cs_file <- cs_file[which.max(file.mtime(cs_file))]
    cs <- readLines(cs_file)
    animal_cells <- extract_val_local(cs, "Animal cells used:")
    rich_idx <- grep("Richness correlation", cs)
    if (length(rich_idx) > 0) {
      cs_rich <- cs[rich_idx[1]:min(rich_idx[1] + 10, length(cs))]
      rich_cells <- extract_val_local(cs_rich, "Cells:")
      rich_r2    <- extract_val_local(cs_rich, "R..:")
      rich_p     <- extract_val_local(cs_rich, "p-value:")
    } else {
      rich_cells <- rich_r2 <- rich_p <- NA_real_
    }
    end_idx <- grep("Endemism correlation", cs)
    if (length(end_idx) > 0) {
      cs_end <- cs[end_idx[1]:min(end_idx[1] + 10, length(cs))]
      end_r2 <- extract_val_local(cs_end, "R..:")
      end_p  <- extract_val_local(cs_end, "p-value:")
    } else {
      end_r2 <- end_p <- NA_real_
    }
  } else {
    animal_cells <- rich_cells <- rich_r2 <- rich_p <- end_r2 <- end_p <- NA_real_
  }

  observed_spp <- if (!is.na(post_species) && filter_applied) post_species else pre_species

  richness_map   <- load_pd("plotdata_animal_richness_map.rds")
  endemism_map   <- load_pd("plotdata_animal_endemism_map.rds")
  endem_corr     <- load_pd("plotdata_endemism_correlation.rds")
  rich_corr      <- load_pd("plotdata_richness_correlation.rds")
  endemism_ratio_data <- load_pd("plotdata_endemism_ratios.rds")
  if (is.null(endemism_ratio_data)) endemism_ratio_data <- load_pd("plotdata_scaling_factors.rds")
  # endem_sim not loaded -- endemism-ratio simulation kept in individual output only
  extent_model   <- load_pd("plotdata_plant_extent_model.rds")
  scaled_map     <- load_pd("plotdata_scaled_endemism_map.rds")
  sim_rich_map   <- load_pd("plotdata_simulated_richness_map.rds")

  panels <- list()

  if (!is.null(richness_map)) {
    panels$richness_map <- ggplot(richness_map) +
      geom_sf(aes(fill = animal_richness), color = NA) +
      scale_fill_viridis_c(option = "magma", direction = -1, name = "Species\nrichness") +
      theme_minimal() + labs(title = "Animal Richness")
  }
  if (!is.null(endemism_map)) {
    panels$endemism_map <- ggplot(endemism_map) +
      geom_sf(aes(fill = animal_endemism_count), color = NA) +
      scale_fill_viridis_c(option = "plasma", direction = -1, name = "Endemic\nspecies") +
      theme_minimal() + labs(title = "Animal Endemism")
  }
  if (!is.null(rich_corr)) {
    r2_val <- tryCatch(
      summary(lm(morphotype_richness ~ plant_richness, data = rich_corr))$r.squared,
      error = function(e) NA_real_)
    panels$rich_corr <- ggplot(rich_corr, aes(x = plant_richness, y = morphotype_richness)) +
      geom_point(alpha = 0.6, color = "darkgreen") +
      geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.8) +
      annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.5,
               label = sprintf("R2 = %.3f  (n = %d)", r2_val, nrow(rich_corr)), size = 3) +
      theme_minimal() +
      labs(title = "Plant vs Animal Richness", x = "Plant richness", y = "Animal richness")
  }
  if (!is.null(endem_corr)) {
    endem_corr_pos <- endem_corr |> filter(raw_endemism > 0, animal_endemism_count > 0)
    r2_val <- tryCatch(
      summary(lm(log(animal_endemism_count) ~ log(raw_endemism), data = endem_corr_pos))$r.squared,
      error = function(e) NA_real_)
    panels$endem_corr <- ggplot(endem_corr_pos, aes(x = raw_endemism, y = animal_endemism_count)) +
      geom_point(alpha = 0.6, color = "darkgreen") +
      geom_smooth(method = "lm", se = TRUE, color = "red", linewidth = 0.8) +
      scale_x_log10() + scale_y_log10() +
      annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.5,
               label = sprintf("R2 = %.3f  (n = %d)", r2_val, nrow(endem_corr_pos)), size = 3) +
      theme_minimal() +
      labs(title = "Plant vs Animal Endemism (log-log)",
           x = "Plant raw endemism (log)", y = "Animal endemism count (log)")
  }
  if (!is.null(endemism_ratio_data)) {
    if ("scaling_factor" %in% names(endemism_ratio_data) && !"endemism_ratio" %in% names(endemism_ratio_data)) {
      endemism_ratio_data <- dplyr::rename(endemism_ratio_data, endemism_ratio = scaling_factor)
    }
    panels$endemism_ratios <- ggplot(endemism_ratio_data, aes(x = endemism_ratio)) +
      geom_histogram(bins = 30, fill = "darkgreen", color = "black") +
      theme_minimal() +
      labs(title = "Endemism Ratios (Animal / Plant, < 2)", x = "Endemism ratio", y = "Frequency")
  }

  # Panel F (Simulated endemic species): SKIPPED -- endemism-ratio simulation
  # kept in individual output only. See endemism_ratio_simulation_results_*.txt.

  if (!is.null(extent_model)) {
    plant_df <- extent_model$plant_df
    ci_df    <- extent_model$ci_df
    ci_df$y_pred      <- extent_model$y_pred
    ci_df$y_predUpper <- extent_model$y_predUpper
    panels$extent_model <- ggplot() +
      geom_ribbon(data = ci_df, aes(x = num_grid_cells, ymin = lower, ymax = upper),
                  fill = "darkgreen", alpha = 0.2) +
      geom_line(data = ci_df, aes(x = num_grid_cells, y = median), color = "darkgreen", linewidth = 0.8) +
      geom_line(data = ci_df, aes(x = num_grid_cells, y = y_pred), color = "steelblue", linetype = "dashed", linewidth = 0.6) +
      geom_line(data = ci_df, aes(x = num_grid_cells, y = y_predUpper), color = "steelblue", linetype = "dashed", linewidth = 0.6) +
      geom_point(data = plant_df, aes(x = num_grid_cells, y = n), color = "darkgreen", size = 1.5) +
      scale_x_log10() + scale_y_log10() + theme_minimal() +
      labs(title = "Plant Spatial Extent Model", x = "Grid cells (log)", y = "Plant species (log)")
  }
  if (!is.null(scaled_map)) {
    panels$scaled_map <- ggplot(scaled_map) +
      geom_sf(aes(fill = scaled_plant_endemism), color = NA) +
      scale_fill_viridis_c(option = "plasma", na.value = "white", name = "Scaled\nendemism") +
      theme_minimal() + labs(title = "Scaled Plant Endemism")
  }
  if (!is.null(sim_rich_map)) {
    panels$sim_richness_map <- ggplot(sim_rich_map) +
      geom_sf(aes(fill = sim_richness), color = NA) +
      scale_fill_viridis_c(option = "magma", direction = -1, na.value = "white", name = "Simulated\nrichness") +
      theme_minimal() + labs(title = "Simulated Species Richness")
  }
  if (!is.na(observed_spp) && !is.na(gamma_best)) {
    pct_undiscovered <- round((gamma_best - observed_spp) / gamma_best * 100, 1)
    gamma_df <- tibble(
      metric  = factor(c("Observed", "Estimated\ngamma"), levels = c("Observed", "Estimated\ngamma")),
      species = c(observed_spp, gamma_best),
      lower   = c(NA_real_, gamma_lower),
      upper   = c(NA_real_, gamma_upper))
    panels$gamma_richness <- ggplot(gamma_df, aes(x = metric, y = species)) +
      geom_col(width = 0.5, fill = c("grey50", "darkgreen")) +
      geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, na.rm = TRUE) +
      annotate("text", x = 1.5, y = max(gamma_upper, gamma_best, na.rm = TRUE) * 0.5,
               label = paste0(pct_undiscovered, "%% est.\nundiscovered"), size = 3.5, fontface = "bold") +
      theme_minimal() + labs(title = "Observed vs Estimated Gamma Richness", x = NULL, y = "Species count")
  }
  ext_file <- file.path(org_ext, "summary_all_models.csv")
  if (file.exists(ext_file)) {
    ext_data <- read.csv(ext_file)
    panels$extinction <- ggplot(ext_data, aes(x = scenario, y = median_spp_lost, fill = model)) +
      geom_col(position = position_dodge(width = 0.8), width = 0.7) +
      geom_errorbar(aes(ymin = q25_spp_lost, ymax = q75_spp_lost),
                    position = position_dodge(width = 0.8), width = 0.25) +
      theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(title = "Species Lost Under Deforestation", x = NULL, y = "Species lost (median, IQR)", fill = "Spatial model")
  }

  if (length(panels) == 0) {
    message("  No plot data available for ", organism, ". Skipping figure.")
    return(invisible(NULL))
  }

  combined <- patchwork::wrap_plots(panels, ncol = 3) +
    patchwork::plot_annotation(
      title = paste0("Summary: ", organism_label(organism)),
      subtitle = paste0("Generated ", Sys.time()),
      theme = theme(plot.title = element_text(size = 16, face = "bold")))

  fig_file <- file.path(org_out, paste0("compilation_", organism, "_", timestamp, ".pdf"))
  n_panels <- length(panels)
  fig_h <- ceiling(n_panels / 3) * 5
  fig_w <- min(n_panels, 3) * 5
  pdf(fig_file, width = fig_w, height = fig_h)
  print(combined)
  dev.off()
  message("  Wrote: ", fig_file)

  stats_tbl <- tibble(
    organism = organism, pre_cells = pre_cells, pre_records = pre_records,
    pre_species = pre_species, filter_applied = filter_applied,
    post_cells = post_cells, post_records = post_records, post_species = post_species,
    model_cells = animal_cells, extent_scaling_factor = extent_scaling_factor,
    observed_spp = observed_spp, gamma_best = gamma_best, gamma_lower = gamma_lower,
    gamma_upper = gamma_upper, gamma_difference = gamma_best - observed_spp,
    gamma_ratio = round(gamma_best / observed_spp, 2),
    rich_corr_cells = rich_cells, rich_r2 = rich_r2, rich_p = rich_p,
    end_r2 = end_r2, end_p = end_p)

  if (file.exists(ext_file)) {
    ext_tbl <- read.csv(ext_file) |> mutate(organism = organism, .before = 1)
    csv_file <- file.path(org_out, paste0("summary_", organism, "_", timestamp, ".csv"))
    write.csv(stats_tbl, csv_file, row.names = FALSE)
    message("  Wrote: ", csv_file)
    ext_csv <- file.path(org_out, paste0("extinction_", organism, "_", timestamp, ".csv"))
    write.csv(ext_tbl, ext_csv, row.names = FALSE)
    message("  Wrote: ", ext_csv)
  } else {
    csv_file <- file.path(org_out, paste0("summary_", organism, "_", timestamp, ".csv"))
    write.csv(stats_tbl, csv_file, row.names = FALSE)
    message("  Wrote: ", csv_file)
  }

  cat("\n=== ", organism_label(organism), " -- Dataset & Model Summary ===\n\n")
  stats_tbl |>
    dplyr::select(organism, pre_species, post_species, model_cells,
                  extent_scaling_factor, observed_spp, gamma_best, gamma_difference,
                  gamma_ratio, rich_r2, rich_p) |>
    as.data.frame() |>
    print(row.names = FALSE)

  invisible(list(stats = stats_tbl, figure_file = fig_file))
}


# =============================================================================
# --- 7. Spatial extent overlap comparison ------------------------------------
# =============================================================================
#
# plot_spatial_extent_overlap() loads per-organism spatial extent distributions
# (animal + plant) and produces a faceted comparison plot for the organisms
# with the smallest animal spatial extents.
#
# Requires plotdata_spatial_extent_distributions.rds (saved by the pipeline
# after the next re-run) or falls back to plotdata_plant_extent_model.rds
# (plant data only — available now for most organisms).

plot_spatial_extent_overlap <- function(include = NULL, group = NULL,
                                        max_cells = 10,
                                        output_dir = NULL, save_pdf = TRUE) {

  if (!exists("project_root")) source("AmAni.manager.R")
  if (is.null(output_dir)) output_dir <- file.path(project_root, "output")

  all_dirs <- list.dirs(output_dir, full.names = FALSE, recursive = FALSE)
  all_dirs <- all_dirs[!all_dirs %in% c("cross_organism_summary", "shared")]

  if (!is.null(group) && group %in% names(ORGANISM_GROUPS)) {
    include <- union(include, ORGANISM_GROUPS[[group]])
  }
  if (!is.null(include)) {
    all_dirs <- all_dirs[all_dirs %in% include]
  }

  # Load combined distributions (preferred) or plant-only fallback
  combined_data <- purrr::map_dfr(all_dirs, function(org) {
    pd_dir <- file.path(output_dir, org, "animal_model", "plot_data")

    # Preferred: saved combined distribution (animal + plant)
    combined_path <- file.path(pd_dir, "plotdata_spatial_extent_distributions.rds")
    if (file.exists(combined_path)) {
      d <- readRDS(combined_path)
      return(dplyr::mutate(d, organism = org))
    }

    # Fallback: plant-only from extent model
    ext_path <- file.path(pd_dir, "plotdata_plant_extent_model.rds")
    if (file.exists(ext_path)) {
      d <- readRDS(ext_path)
      return(dplyr::mutate(d$plant_df, organism = org))
    }
    NULL
  })

  if (nrow(combined_data) == 0) {
    message("No spatial extent data found.")
    return(invisible(NULL))
  }

  # Show all matched organisms, filtered to 1..max_cells on x-axis
  focus_orgs <- unique(combined_data$organism)

  plot_data <- combined_data |>
    dplyr::filter(num_grid_cells <= max_cells) |>
    dplyr::mutate(
      organism = organism_label(organism),
      num_grid_cells = factor(num_grid_cells, levels = seq_len(max_cells))
    )

  n_orgs <- dplyr::n_distinct(plot_data$organism)

  p <- ggplot(plot_data, aes(x = num_grid_cells, y = n, fill = group)) +
    geom_col(position = "dodge", color = "grey30", linewidth = 0.2) +
    scale_fill_manual(
      values = c("Plant species" = "#2d6a4f", "Animal morphotypes" = "#457b9d"),
      na.value = "grey60", name = NULL
    ) +
    facet_wrap(~organism, scales = "free_y",
               ncol = min(n_orgs, 4)) +
    labs(
      title = paste0("Spatial Extent Distribution (1\u2013", max_cells, " grid cells)"),
      subtitle = if (!is.null(group)) paste0("Group: ", group) else NULL,
      x = "Number of grid cells occupied",
      y = "Number of species"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(size = 7),
      legend.position = "top"
    )

  if (save_pdf) {
    summary_dir <- file.path(output_dir, "cross_organism_summary")
    dir.create(summary_dir, showWarnings = FALSE, recursive = TRUE)
    ts <- format(Sys.time(), "%Y%m%d_%H%M%S")
    out_file <- file.path(summary_dir,
                          paste0("spatial_extent_overlap_", ts, ".pdf"))
    pdf(out_file, width = 12, height = 8)
    print(p)
    dev.off()
    message("Wrote: ", out_file)
  }

  p
}


# =============================================================================
# --- 8. Workflow compilation diagram -----------------------------------------
# =============================================================================
#
# compile_workflow_diagram() produces a ggplot-based flowchart illustrating
# the major input/output stages of the pipeline, annotated with example
# values from a chosen organism.

compile_workflow_diagram <- function(organism = "butterflies",
                                     output_dir = NULL,
                                     save_pdf = TRUE) {

  if (!exists("project_root")) source("AmAni.manager.R")
  if (is.null(output_dir)) output_dir <- file.path(project_root, "output")

  # --- Helper: extract first number after a colon from text lines ---
  .ev <- function(lines, pattern) {
    l <- grep(pattern, lines, value = TRUE, ignore.case = TRUE)
    if (!length(l)) return("?")
    raw <- sub(".*:\\s*", "", l[1])
    m <- regmatches(raw, regexpr("[0-9][0-9,]*\\.?[0-9]*", raw))
    if (!length(m)) "?" else m
  }

  org_am  <- file.path(output_dir, organism, "animal_model")
  org_ext <- file.path(output_dir, organism, "extinction_model")

  # --- Load example values ---
  pre_cells <- pre_spp <- post_cells <- post_spp <- "?"
  ds_files <- list.files(org_am, "^dataset_summary.*\\.txt$", full.names = TRUE)
  if (length(ds_files)) {
    ds <- readLines(ds_files[which.max(file.mtime(ds_files))])
    pre_cells <- .ev(ds, "^\\s*Grid cells:")
    pre_spp   <- .ev(ds, "^\\s*Total species:")
    ai <- grep("After filtering", ds)
    if (length(ai)) {
      da <- ds[ai[1]:length(ds)]
      post_cells <- .ev(da, "^\\s*Grid cells:")
      post_spp   <- .ev(da, "^\\s*Total species:")
    }
  }

  gamma_best <- "?"
  gm_files <- list.files(org_am, "^gamma_richness.*\\.txt$", full.names = TRUE)
  if (length(gm_files)) {
    gm <- readLines(gm_files[which.max(file.mtime(gm_files))])
    gamma_best <- .ev(gm, "Best estimate:")
  }

  med_lost <- "?"
  ext_file <- file.path(org_ext, "summary_all_models.csv")
  if (file.exists(ext_file)) {
    ext <- read.csv(ext_file)
    r <- ext[ext$scenario == "FreqBoundary" & ext$model == "Clustered", ]
    if (nrow(r)) med_lost <- as.character(round(r$median_spp_lost[1]))
  }

  org_display <- organism_label(organism)

  # --- Box definitions: cx, cy, label, fill ---
  bw <- 2.0   # box half-width
  bh <- 0.65  # box half-height

  boxes <- dplyr::tribble(
    ~cx,   ~cy,  ~label,                                                       ~fill,
    # Plant lane (top)
    1.5,    3.0,  "GBIF Plants\n(Amazon, ~9k spp)",                            "#2d6a4f",
    5.0,    3.0,  "Stage 0\u20131\nPlant grid + endemism\n(22 km cells)",      "#52b788",

    # Animal lane (middle)
    1.5,    1.5,  paste0(org_display, "\n",
                         pre_spp, " spp / ", pre_cells, " cells"),             "#1d3557",
    5.0,    1.5,  paste0("Stage 2\n",
                         post_spp, " spp / ", post_cells, " cells\n",
                         "Extent model + scaling"),                            "#457b9d",
    8.5,    1.5,  paste0("Simulated Occ.\n\u2248",
                         gamma_best, " spp\n(gamma est.)"),                    "#264653",

    # Deforestation lane (bottom)
    1.5,    0.0,  "RAISG Raster\n(2001\u20132020)",                            "#6d4c41",
    5.0,    0.0,  "Stage 3\nPatch extraction\nLognormal + power-law",          "#a1887f",

    # Extinction + results (right)
    11.5,   0.75, paste0("Stage 4: Extinction\nNeutral / Clustered\n",
                         "/ Shifted (100 reps)"),                              "#c1121f",
    14.5,   0.75, paste0("Results\n\u2248", med_lost,
                         " spp lost\n(FreqBound., Clustered)"),                "#e76f51"
  )

  boxes <- boxes |>
    dplyr::mutate(
      xmin = cx - bw / 2, xmax = cx + bw / 2,
      ymin = cy - bh,      ymax = cy + bh
    )

  # --- Arrow definitions ---
  arrows <- dplyr::tribble(
    ~x1,    ~y1,   ~x2,    ~y2,
    # Plant GBIF → Stage 0-1
    2.5,    3.0,   3.0,    3.0,
    # Stage 0-1 → Stage 2 (plant endemism reference)
    5.0,    2.35,  5.0,    2.15,
    # Animal data → Stage 2
    2.5,    1.5,   3.0,    1.5,
    # Stage 2 → Simulated occ.
    6.0,    1.5,   6.5,    1.5,
    # Simulated occ → Stage 4
    9.5,    1.5,   10.5,   1.15,
    # Deforestation → Stage 3
    2.5,    0.0,   3.0,    0.0,
    # Stage 3 → Stage 4
    6.0,    0.0,   10.5,   0.4,
    # Stage 4 → Results
    12.5,   0.75,  13.5,   0.75
  )

  # --- Build plot ---
  p <- ggplot() +
    geom_rect(data = boxes,
              aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              fill = boxes$fill, color = "white", linewidth = 0.9,
              alpha = 0.92) +
    geom_text(data = boxes,
              aes(x = cx, y = cy, label = label),
              size = 2.5, color = "white", fontface = "bold",
              lineheight = 1.0) +
    geom_segment(data = arrows,
                 aes(x = x1, y = y1, xend = x2, yend = y2),
                 arrow = arrow(length = unit(0.18, "cm"), type = "closed"),
                 color = "grey40", linewidth = 0.55) +
    # Lane labels
    annotate("text", x = 0.0, y = 3.0, label = "Plants",
             angle = 90, size = 3.2, color = "#2d6a4f", fontface = "bold") +
    annotate("text", x = 0.0, y = 1.5, label = "Animals",
             angle = 90, size = 3.2, color = "#1d3557", fontface = "bold") +
    annotate("text", x = 0.0, y = 0.0, label = "Deforestation",
             angle = 90, size = 3.2, color = "#6d4c41", fontface = "bold") +
    coord_cartesian(xlim = c(-0.5, 16.5), ylim = c(-1.0, 3.9),
                    expand = FALSE) +
    labs(title = paste0("AmAni Pipeline Overview: ", org_display),
         subtitle = "Major input/output stages (minor filtering steps omitted)") +
    theme_void() +
    theme(
      plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 9, hjust = 0.5, color = "grey40"),
      plot.background = element_rect(fill = "white", color = NA),
      plot.margin = margin(10, 10, 10, 10)
    )

  if (save_pdf) {
    summary_dir <- file.path(output_dir, "cross_organism_summary")
    dir.create(summary_dir, showWarnings = FALSE, recursive = TRUE)
    ts <- format(Sys.time(), "%Y%m%d_%H%M%S")
    out_file <- file.path(summary_dir,
                          paste0("workflow_diagram_", organism, "_", ts, ".pdf"))
    pdf(out_file, width = 18, height = 6)
    print(p)
    dev.off()
    message("Wrote: ", out_file)
  }

  p
}

