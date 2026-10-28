# utils/summarize_all.R
# Extracted from AmAni.summarize.R. Sourced by AmAni.summarize.R — relies on
# the packages (tidyverse, sf, patchwork, ggrepel) loaded there.

summarize_all <- function(include = NULL, group = NULL,
                          results_name = "cross_organism_summary",
                          archive = TRUE) {

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
  # Also restrict the full directory list so downstream sections (compilation
  # figures in section 4, dataset summary in section 5) respect the filter.
  organism_dirs_all <- organism_dirs_all[organism_dirs_all %in% include]
}

if (length(organisms) == 0) {
  message("No completed organism results found. Skipping summary.")
} else {

message("Summarizing results for: ", paste(organisms, collapse = ", "))

timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")

summary_dir <- file.path(output_dir, results_name)
dir.create(summary_dir, showWarnings = FALSE, recursive = TRUE)

# --- Archive existing results before writing new files ----------------------
# Move any pre-existing output into _archive/<timestamp>/ so each run starts
# clean without destroying prior results. The _archive dir itself is skipped.
if (archive) {
  archive_root <- file.path(summary_dir, "_archive")
  existing <- list.files(summary_dir, full.names = TRUE, all.files = FALSE)
  existing <- existing[basename(existing) != "_archive"]
  if (length(existing) > 0) {
    archive_dir <- file.path(archive_root, timestamp)
    dir.create(archive_dir, showWarnings = FALSE, recursive = TRUE)
    ok <- file.rename(existing, file.path(archive_dir, basename(existing)))
    # file.rename can fail across filesystems; fall back to copy + delete.
    if (any(!ok)) {
      failed <- existing[!ok]
      file.copy(failed, archive_dir, recursive = TRUE)
      unlink(failed, recursive = TRUE)
    }
    message("  Archived ", length(existing), " existing item(s) to ",
            file.path(results_name, "_archive", timestamp))
  }
}


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

# Include any organism with completed gamma richness output, even if its
# extinction model hasn't finished yet — extinction status only affects the
# "extinct" component of the stacked composition plot further below, not
# whether the organism appears in the gamma richness plots at all.
gamma_organisms <- organism_dirs_all[sapply(organism_dirs_all, function(org) {
  length(list.files(file.path(output_dir, org, "animal_model"),
                    pattern = "^gamma_richnessfromSPRICH.*\\.txt$")) > 0
})]

gamma_summary <- map_dfr(gamma_organisms, function(org) {
  gamma_files <- list.files(
    file.path(output_dir, org, "animal_model"),
    pattern = "^gamma_richnessfromSPRICH.*\\.txt$",
    full.names = TRUE
  )
  if (length(gamma_files) == 0) {
    return(tibble(organism = org, gamma_richness = NA_real_,
                  upper_bound = NA_real_, lower_bound = NA_real_,
                  surrogate_richness = NA_real_,
                  surrogate_upper = NA_real_, surrogate_lower = NA_real_,
                  richness_shortfall = NA_real_))
  }
  gamma_file <- gamma_files[which.max(file.mtime(gamma_files))]
  lines <- readLines(gamma_file)

  parse_val <- function(pattern) {
    line <- grep(pattern, lines, value = TRUE, ignore.case = TRUE)
    if (length(line) == 0) return(NA_real_)
    raw <- sub(".*:\\s*", "", line[1])
    m <- regmatches(raw, regexpr("-?[0-9]+\\.?[0-9]*", raw))
    if (length(m) == 0) return(NA_real_)
    as.numeric(m)
  }

  # New format: total model (surrogate + extension) is the primary estimate
  total_best  <- parse_val("^Best estimate from total model")
  total_upper <- parse_val("^Upper bound:\\s")
  total_lower <- parse_val("^Lower bound:\\s")

  # New format: surrogate-only estimates
  surr_best  <- parse_val("^Best estimate from plant surrogate")
  surr_upper <- parse_val("^Upper bound \\(narrow")
  surr_lower <- parse_val("^Lower bound \\(narrow")

  shortfall <- parse_val("^Expected proportional richness shortfall")

  # Fall back to old format if total model fields are missing
  if (is.na(total_best)) {
    total_best  <- parse_val("^Best estimate|^gammaRichness")
    total_upper <- parse_val("^Upper bound")
    total_lower <- parse_val("^Lower bound")
  }

  # Observed species counts from the dataset summary: the pre-filter total
  # (used as the observed baseline in gamma-change plots) and the analyzed
  # portion (post-filter if a filter was applied, else pre-filter).
  observed_total <- NA_real_
  observed_analyzed <- NA_real_
  ds_files <- list.files(
    file.path(output_dir, org, "animal_model"),
    pattern = "^dataset_summary_.*\\.txt$", full.names = TRUE)
  if (length(ds_files) > 0) {
    ds <- readLines(ds_files[which.max(file.mtime(ds_files))])
    grab_species <- function(txt) {
      l <- grep("^\\s*Total species:", txt, value = TRUE)
      if (!length(l)) return(NA_real_)
      raw <- sub(".*:\\s*", "", l[1])
      m <- regmatches(raw, regexpr("[0-9]+\\.?[0-9]*", raw))
      if (!length(m)) NA_real_ else as.numeric(m)
    }
    observed_total <- grab_species(ds)
    filt <- any(grepl("Filter applied:\\s*YES", ds))
    ai <- grep("After filtering", ds)
    if (length(ai) > 0) {
      post_sp <- grab_species(ds[ai[1]:length(ds)])
      observed_analyzed <- if (!is.na(post_sp) && filt) post_sp else observed_total
    } else {
      observed_analyzed <- observed_total
    }
  }

  tibble(
    organism = org,
    gamma_richness = total_best,
    upper_bound = total_upper,
    lower_bound = total_lower,
    surrogate_richness = surr_best,
    surrogate_upper = surr_upper,
    surrogate_lower = surr_lower,
    richness_shortfall = shortfall,
    observed_total = observed_total,
    observed_analyzed = observed_analyzed
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

# 3c. Gamma richness comparison — observed (pre-filter total) vs estimated gamma
if (any(!is.na(gamma_summary$gamma_richness))) {
  pdf(file.path(summary_dir, paste0("gamma_richness_comparison_", timestamp, ".pdf")),
      width = 9, height = 5)

  gamma_plot <- gamma_summary |> filter(!is.na(gamma_richness))

  # Long form: observed species column beside the estimated gamma column.
  gamma_long <- gamma_plot |>
    dplyr::select(organism, observed_total, gamma_richness, lower_bound, upper_bound) |>
    tidyr::pivot_longer(cols = c(observed_total, gamma_richness),
                        names_to = "metric", values_to = "species") |>
    dplyr::mutate(
      # Error bars only apply to the estimated gamma bar.
      lower_bound = ifelse(metric == "gamma_richness", lower_bound, NA_real_),
      upper_bound = ifelse(metric == "gamma_richness", upper_bound, NA_real_),
      metric = factor(metric,
                      levels = c("observed_total", "gamma_richness"),
                      labels = c("Observed (total)", "Estimated gamma"))
    )

  p <- ggplot(gamma_long, aes(x = organism, y = species, fill = metric)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7) +
    geom_errorbar(aes(ymin = lower_bound, ymax = upper_bound),
                  position = position_dodge(width = 0.8), width = 0.25, na.rm = TRUE) +
    labs(title = "Estimated Amazonian Gamma Richness by Organism",
         subtitle = "Observed = total species recorded (pre-filter)",
         x = NULL, y = "Species", fill = NULL) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
  print(p)
  dev.off()
  message("  Wrote: gamma_richness_comparison_", timestamp, ".pdf")

  # 3c-ii. Single stacked column per organism, partitioned into observed,
  # additional predicted (unsampled), and projected extinct species.
  extinct_by_org <- extinction_summary |>
    dplyr::filter(scenario == "FreqBoundary", model == "Clustered") |>
    dplyr::select(organism, extinct = median_spp_lost)

  stacked <- gamma_plot |>
    dplyr::select(organism, observed_total, gamma_richness) |>
    dplyr::left_join(extinct_by_org, by = "organism") |>
    dplyr::transmute(
      organism,
      `Observed`             = observed_total,
      `Predicted (unsampled)` = pmax(gamma_richness - observed_total, 0),
      `Extinct (projected)`  = dplyr::coalesce(extinct, 0)
    ) |>
    tidyr::pivot_longer(cols = -organism,
                        names_to = "category", values_to = "species") |>
    dplyr::mutate(category = factor(category,
      levels = c("Observed", "Predicted (unsampled)", "Extinct (projected)")))

  if (any(!is.na(stacked$species))) {
    # Total species per organism (incl. projected extinctions), labeled in the
    # middle of each proportional column.
    stacked_totals <- stacked |>
      dplyr::group_by(organism) |>
      dplyr::summarise(total = sum(species, na.rm = TRUE), .groups = "drop")

    pdf(file.path(summary_dir, paste0("gamma_richness_stacked_", timestamp, ".pdf")),
        width = 9, height = 5)
    p <- ggplot(stacked, aes(x = organism, y = species, fill = category)) +
      geom_col(position = "fill", width = 0.7) +
      geom_text(data = stacked_totals,
                aes(x = organism, y = 0.5, label = scales::comma(round(total))),
                inherit.aes = FALSE, size = 3, fontface = "bold", color = "white") +
      scale_y_continuous(labels = scales::percent) +
      scale_fill_manual(values = c(
        "Observed" = "grey50",
        "Predicted (unsampled)" = "#2d6a4f",
        "Extinct (projected)" = "#c1121f")) +
      labs(title = "Species Composition by Organism (proportional)",
           subtitle = "Columns scaled to 100%; label = total species incl. projected extinctions. Extinct = FreqBoundary/Clustered median",
           x = NULL, y = "Percent of total", fill = NULL) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 30, hjust = 1))
    print(p)
    dev.off()
    message("  Wrote: gamma_richness_stacked_", timestamp, ".pdf")
  }
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


# 3e. IUCN threat classification stacked bars across organisms
# Collects per-organism IUCN threat data (Clustered / FreqBoundary) and
# produces two summary figures:
#   i.  Proportional stacked bar
#   ii. Absolute stacked bar on log scale with species counts labeled

iucn_all <- map_dfr(organisms, function(org) {
  iucn_file <- list.files(
    file.path(output_dir, org, "extinction_model"),
    pattern = "iucn_threat_classification_data\\.rds$", full.names = TRUE)
  if (length(iucn_file) == 0) return(NULL)
  iucn_data <- readRDS(iucn_file[which.max(file.mtime(iucn_file))])
  iucn_data$counts |>
    dplyr::filter(model == "Clustered", scenario == "FreqBoundary") |>
    dplyr::mutate(organism = organism_label(org))
})

if (nrow(iucn_all) > 0) {
  iucn_cat_levels <- c("NT/LC", "VU", "EN", "CR/EX")
  iucn_cols <- c("NT/LC" = "#60C659", "VU" = "#F9E814",
                 "EN" = "#FC7F3F", "CR/EX" = "#D81E05")
  iucn_all <- iucn_all |>
    dplyr::mutate(
      count = round(count, 0),
      iucn_category = factor(iucn_category, levels = iucn_cat_levels))

  # 3e-i. Proportional axes
  pdf(file.path(summary_dir, paste0("iucn_threat_stack_proportional_", timestamp, ".pdf")),
      width = 10, height = 6)
  p <- ggplot(iucn_all, aes(x = organism, y = count, fill = iucn_category)) +
    geom_col(position = "fill", width = 0.7) +
    scale_y_continuous(labels = scales::percent, limits = c(0, NA)) +
    scale_fill_manual(values = iucn_cols, drop = FALSE, name = "IUCN Category") +
    labs(title = "IUCN Threat Classification by Organism (Proportional)/nBased on past deforestation only",
         subtitle = "Clustered model, FreqBoundary scenario",
         x = NULL, y = "Proportion of species") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
  print(p)
  dev.off()
  message("  Wrote: iucn_threat_stack_proportional_", timestamp, ".pdf")

  # 3e-ii. Absolute axes, y-axis starting at 0, with species counts as dodged
  # bars. IUCN categories ordered from least concern to most threatened.
  iucn_dodge <- iucn_all |>
    dplyr::filter(count > 0)

  pdf(file.path(summary_dir, paste0("iucn_threat_absolute_", timestamp, ".pdf")),
      width = 10, height = 6)
  p <- ggplot(iucn_dodge, aes(x = organism, y = count, fill = iucn_category)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7) +
    geom_text(aes(label = scales::comma(count)),
              position = position_dodge(width = 0.8), vjust = -0.3,
              size = 2.5) +
    scale_y_continuous(labels = scales::comma, limits = c(0, NA),
                       expand = expansion(mult = c(0, 0.1))) +
    scale_fill_manual(values = iucn_cols, drop = FALSE, name = "IUCN Category") +
    labs(title = "IUCN Threat Classification by Organism\nnBased on past deforestation only",
         subtitle = "Clustered model, FreqBoundary scenario",
         x = NULL, y = "Number of species") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
  print(p)
  dev.off()
  message("  Wrote: iucn_threat_absolute_", timestamp, ".pdf")
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
    labs(title = "Spatial Extent Model with Bootstrapped CIs\n(surrogate)",
         x = "Number of grid cells (log scale)", y = "Number of species (log scale)")
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


# --- 4.5 Compilation of dataset before/after filtering summaries -----------
#
# One multi-page PDF containing every organism's dataset-filtering summary
# figure (thresholds + retention + per-cell distributions), one page per
# organism -- a "compilation of compilations".

message("\n--- Compiling dataset-filtering summaries across organisms ---")

ds_summary_figs <- purrr::compact(purrr::map(organism_dirs_all, function(org) {
  am_dir <- file.path(output_dir, org, "animal_model")
  inputs <- load_dataset_summary_figure_inputs(am_dir)
  if (is.null(inputs)) return(NULL)
  fig <- build_dataset_summary_figure(
    inputs$thresholds, inputs$counts, inputs$dists,
    title = paste("Dataset before vs after cell filtering \u2014", organism_label(org)))
  if (is.null(fig)) return(NULL)
  fig
}))

if (length(ds_summary_figs) > 0) {
  ds_compilation_file <- file.path(summary_dir,
    paste0("compilation_dataset_summaries_", timestamp, ".pdf"))
  pdf(ds_compilation_file, width = 10, height = 10)
  for (fig in ds_summary_figs) print(fig)
  dev.off()
  message("  Wrote: ", ds_compilation_file, " (", length(ds_summary_figs), " organism pages)")
} else {
  message("  No dataset_summary_*.txt files found. Skipping.")
}


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
      ds_after <- character(0)
      post_cells <- post_records <- post_species <- NA_real_
    }

    filter_applied <- any(grepl("Filter applied:\\s*YES", ds))
    min_occ  <- extract_val(ds, "min_occurrences:")
    min_comp <- extract_val(ds, "min_completeness:")
    min_rich <- extract_val(ds, "min_richness:")

    # Outlier cutoff and Track B entry
    outlier_cutoff <- extract_val(ds, "outlier cutoff:", as_type = "character")
    track_b_entry  <- extract_val(ds, "Track B entry:", as_type = "character")

    # Well-sampled / super-diverse / outliers dropped from after-filter Grid cells line
    gc_line <- grep("^\\s*Grid cells:", ds_after, value = TRUE)
    if (length(gc_line) > 0) {
      well_sampled_m   <- regmatches(gc_line[1], regexpr("([0-9]+)\\s+well-sampled", gc_line[1]))
      super_diverse_m  <- regmatches(gc_line[1], regexpr("([0-9]+)\\s+super-diverse", gc_line[1]))
      outliers_drop_m  <- regmatches(gc_line[1], regexpr("([0-9]+)\\s+outliers? dropped", gc_line[1]))
      well_sampled   <- if (length(well_sampled_m))  as.integer(sub("\\s.*", "", well_sampled_m))  else NA_integer_
      super_diverse  <- if (length(super_diverse_m)) as.integer(sub("\\s.*", "", super_diverse_m)) else NA_integer_
      outliers_dropped <- if (length(outliers_drop_m)) as.integer(sub("\\s.*", "", outliers_drop_m)) else NA_integer_
    } else {
      well_sampled <- super_diverse <- outliers_dropped <- NA_integer_
    }

    # Species loss
    spp_loss_pct <- extract_val(ds, "Species lost:")
  } else {
    pre_cells <- pre_records <- pre_species <- NA_real_
    post_cells <- post_records <- post_species <- NA_real_
    filter_applied <- FALSE
    min_occ <- min_comp <- min_rich <- NA_real_
    outlier_cutoff <- track_b_entry <- NA_character_
    well_sampled <- super_diverse <- outliers_dropped <- NA_integer_
    spp_loss_pct <- NA_real_
  }

  # --- Gamma richness + scaling factor ---
  gm_file <- list.files(am_dir, pattern = "^gamma_richnessfromSPRICH.*\\.txt$", full.names = TRUE)
  if (length(gm_file) > 0) {
    gm_file <- gm_file[which.max(file.mtime(gm_file))]
    gm <- readLines(gm_file)
    scaling_factor <- extract_val(gm, "Scaling factor:")

    # Which parameters were fit (e.g. "[amplitude]" or "[amplitude, meanlog, sdlog]")
    fit_params <- extract_val(gm, "Fit parameters:", as_type = "character")

    # delta_meanlog / ratio_sdlog and their active/off status
    gm_delta_meanlog     <- extract_val(gm, "^delta_meanlog:")
    gm_delta_meanlog_str <- extract_val(gm, "^delta_meanlog:", as_type = "character")
    gm_delta_active      <- grepl("active", gm_delta_meanlog_str, ignore.case = TRUE)
    gm_ratio_sdlog       <- extract_val(gm, "^ratio_sdlog:")
    gm_ratio_sdlog_str   <- extract_val(gm, "^ratio_sdlog:", as_type = "character")
    gm_ratio_active      <- grepl("active", gm_ratio_sdlog_str, ignore.case = TRUE)

    # Plant range filling
    fill_plants_gm <- grepl("Plant range filling.*was used", paste(gm, collapse = " "))

    # New format: total model (surrogate + extension) as primary
    gamma_best  <- extract_val(gm, "^Best estimate from total model")
    gamma_upper <- extract_val(gm, "^Upper bound:\\s")
    gamma_lower <- extract_val(gm, "^Lower bound:\\s")

    # Fall back to old format
    if (is.na(gamma_best)) {
      gamma_best  <- extract_val(gm, "Best estimate:")
      gamma_upper <- extract_val(gm, "Upper bound:")
      gamma_lower <- extract_val(gm, "Lower bound:")
    }

    # Surrogate-only estimates (new format)
    surrogate_best  <- extract_val(gm, "^Best estimate from plant surrogate")
    surrogate_upper <- extract_val(gm, "^Upper bound \\(narrow")
    surrogate_lower <- extract_val(gm, "^Lower bound \\(narrow")
    richness_shortfall <- extract_val(gm, "^Expected proportional richness shortfall")
  } else {
    scaling_factor <- gamma_best <- gamma_upper <- gamma_lower <- NA_real_
    surrogate_best <- surrogate_upper <- surrogate_lower <- richness_shortfall <- NA_real_
    fit_params <- NA_character_
    gm_delta_meanlog <- gm_ratio_sdlog <- NA_real_
    gm_delta_active <- gm_ratio_active <- FALSE
    fill_plants_gm <- NA
  }

  # --- Sampling extension report ---
  se_file <- list.files(am_dir, pattern = "^sampling_extension_report_.*\\.txt$", full.names = TRUE)
  if (length(se_file) > 0) {
    se_file <- se_file[which.max(file.mtime(se_file))]
    se <- readLines(se_file)

    fill_plants       <- extract_val(se, "fill_Plants:", as_type = "character")
    comp_spread       <- extract_val(se, "comp_spread:")
    comp_mean         <- extract_val(se, "comp_mean:")
    comp_sd           <- extract_val(se, "comp_sd:")
    prop_species      <- extract_val(se, "prop_species \\(proxy\\):")
    chao2_s_obs       <- extract_val(se, "S_obs =")
    chao2_s_chao2     <- extract_val(se, "S_Chao2 =")
    se_n_grid_cells   <- extract_val(se, "n_grid_cells \\(raw\\):")
    se_mean_range     <- extract_val(se, "sample_mean_range_size:")
    se_delta_meanlog  <- extract_val(se, "^\\s*delta_meanlog:")
    se_ratio_sdlog    <- extract_val(se, "^\\s*ratio_sdlog:")
    partial_gamma     <- extract_val(se, "partialGamma:")
    abs_prop_error    <- extract_val(se, "abs_prop_error:")
    correction_factor <- extract_val(se, "correction_factor:")
    full_gamma_best   <- extract_val(se, "fullGamma \\(best\\):")
    full_gamma_lower  <- extract_val(se, "fullGamma \\(lower")
    full_gamma_upper  <- extract_val(se, "fullGamma \\(upper")
  } else {
    fill_plants <- NA_character_
    comp_spread <- comp_mean <- comp_sd <- prop_species <- NA_real_
    chao2_s_obs <- chao2_s_chao2 <- NA_real_
    se_n_grid_cells <- se_mean_range <- NA_real_
    se_delta_meanlog <- se_ratio_sdlog <- NA_real_
    partial_gamma <- abs_prop_error <- correction_factor <- NA_real_
    full_gamma_best <- full_gamma_lower <- full_gamma_upper <- NA_real_
  }

  # --- Correlation summary ---
  cs_file <- list.files(am_dir, pattern = "^correlation_summary_.*\\.txt$", full.names = TRUE)
  if (length(cs_file) > 0) {
    cs_file <- cs_file[which.max(file.mtime(cs_file))]
    cs <- readLines(cs_file)

    # Animal cells used
    animal_cells <- extract_val(cs, "Animal cells used:")

    # Filter settings
    plant_comp_filter <- extract_val(cs, "Plant completeness filter:", as_type = "character")

    # Richness correlation
    rich_idx <- grep("Richness correlation", cs)
    if (length(rich_idx) > 0) {
      cs_rich <- cs[rich_idx[1]:min(rich_idx[1] + 10, length(cs))]
      rich_cells  <- extract_val(cs_rich, "Cells:")
      rich_r2     <- extract_val(cs_rich, "^\\s*R.:")
      rich_adj_r2 <- extract_val(cs_rich, "Adj\\. R.:")
      rich_p      <- extract_val(cs_rich, "p-value:")
      rich_slope  <- extract_val(cs_rich, "Slope:")
    } else {
      rich_cells <- rich_r2 <- rich_adj_r2 <- rich_p <- rich_slope <- NA_real_
    }

    # Endemism correlation
    end_idx <- grep("Endemism correlation", cs)
    if (length(end_idx) > 0) {
      cs_end <- cs[end_idx[1]:min(end_idx[1] + 10, length(cs))]
      end_r2     <- extract_val(cs_end, "^\\s*R.:")
      end_adj_r2 <- extract_val(cs_end, "Adj\\. R.:")
      end_p      <- extract_val(cs_end, "p-value:")
      end_slope  <- extract_val(cs_end, "Slope:")
    } else {
      end_r2 <- end_adj_r2 <- end_p <- end_slope <- NA_real_
    }
  } else {
    animal_cells <- rich_cells <- rich_r2 <- rich_adj_r2 <- rich_p <- rich_slope <- NA_real_
    end_r2 <- end_adj_r2 <- end_p <- end_slope <- NA_real_
    plant_comp_filter <- NA_character_
  }

  # Use observed species: post-filter if filter was applied, otherwise pre-filter
  observed_spp <- if (!is.na(post_species) && filter_applied) post_species else pre_species

  tibble(
    organism         = org,
    # --- Pre-filter ---
    pre_cells        = pre_cells,
    pre_records      = pre_records,
    pre_species      = pre_species,
    # --- Filter settings ---
    filter_applied   = filter_applied,
    min_occurrences  = min_occ,
    min_completeness = min_comp,
    min_richness     = min_rich,
    outlier_cutoff   = outlier_cutoff,
    track_b_entry    = track_b_entry,
    # --- Post-filter ---
    post_cells       = post_cells,
    post_records     = post_records,
    post_species     = post_species,
    well_sampled     = well_sampled,
    super_diverse    = super_diverse,
    outliers_dropped = outliers_dropped,
    spp_loss_pct     = spp_loss_pct,
    # --- Model cells & scaling ---
    model_cells      = animal_cells,
    scaling_factor   = scaling_factor,
    fit_params       = fit_params,
    fill_plants      = fill_plants_gm,
    delta_meanlog    = gm_delta_meanlog,
    delta_meanlog_active = gm_delta_active,
    ratio_sdlog      = gm_ratio_sdlog,
    ratio_sdlog_active = gm_ratio_active,
    # --- Gamma richness ---
    observed_spp     = observed_spp,
    gamma_best       = gamma_best,
    gamma_lower      = gamma_lower,
    gamma_upper      = gamma_upper,
    surrogate_best   = surrogate_best,
    surrogate_lower  = surrogate_lower,
    surrogate_upper  = surrogate_upper,
    richness_shortfall = richness_shortfall,
    gamma_difference = gamma_best - observed_spp,
    gamma_ratio      = round(gamma_best / observed_spp, 2),
    # --- Sampling extension model ---
    se_fill_plants       = fill_plants,
    se_comp_spread       = comp_spread,
    se_comp_mean         = comp_mean,
    se_comp_sd           = comp_sd,
    se_prop_species      = prop_species,
    se_chao2_s_obs       = chao2_s_obs,
    se_chao2_s_chao2     = chao2_s_chao2,
    se_n_grid_cells      = se_n_grid_cells,
    se_mean_range_size   = se_mean_range,
    se_delta_meanlog     = se_delta_meanlog,
    se_ratio_sdlog       = se_ratio_sdlog,
    se_partial_gamma     = partial_gamma,
    se_abs_prop_error    = abs_prop_error,
    se_correction_factor = correction_factor,
    se_full_gamma_best   = full_gamma_best,
    se_full_gamma_lower  = full_gamma_lower,
    se_full_gamma_upper  = full_gamma_upper,
    # --- Correlation summary ---
    plant_comp_filter = plant_comp_filter,
    rich_corr_cells  = rich_cells,
    rich_r2          = rich_r2,
    rich_adj_r2      = rich_adj_r2,
    rich_p           = rich_p,
    rich_slope       = rich_slope,
    end_r2           = end_r2,
    end_adj_r2       = end_adj_r2,
    end_p            = end_p,
    end_slope        = end_slope
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
                  scaling_factor, fit_params, observed_spp,
                  gamma_best, gamma_difference, gamma_ratio,
                  se_correction_factor, richness_shortfall,
                  rich_r2, rich_p) |>
    as.data.frame() |>
    print(row.names = FALSE)

  # --- Slope graph: observed richness → surrogate gamma → full-model gamma ---
  # Observed uses the pre-filter total species; the analyzed portion is noted
  # as text next to the observed column. Uncertainty is shown as a translucent
  # ribbon (lower/upper bound at each step) instead of error-bar whiskers,
  # colored to match each organism's line; at Observed (no uncertainty) the
  # ribbon collapses to a single point and fans out from there.
  slope_base <- organism_summary |>
    dplyr::filter(!is.na(pre_species) & !is.na(gamma_best))

  slope_data <- dplyr::bind_rows(
    slope_base |> dplyr::transmute(
      organism, metric = "Observed\n(total)", species = pre_species,
      lower = pre_species, upper = pre_species, analyzed = observed_spp),
    slope_base |> dplyr::transmute(
      organism, metric = "Surrogate\ngamma", species = surrogate_best,
      lower = dplyr::coalesce(surrogate_lower, surrogate_best),
      upper = dplyr::coalesce(surrogate_upper, surrogate_best),
      analyzed = NA_real_),
    slope_base |> dplyr::transmute(
      organism, metric = "Full-model\ngamma", species = gamma_best,
      lower = dplyr::coalesce(gamma_lower, gamma_best),
      upper = dplyr::coalesce(gamma_upper, gamma_best),
      analyzed = NA_real_)
  ) |>
    dplyr::mutate(metric = factor(metric,
      levels = c("Observed\n(total)", "Surrogate\ngamma", "Full-model\ngamma")))

  if (nrow(slope_data) > 0) {
    pdf(file.path(summary_dir, paste0("slope_observed_vs_gamma_", timestamp, ".pdf")),
        width = 9, height = 6)

    p <- ggplot(slope_data, aes(x = metric, y = species, group = organism)) +
      geom_ribbon(aes(ymin = lower, ymax = upper, fill = organism),
                  alpha = 0.18, color = NA) +
      geom_line(aes(color = organism), linewidth = 0.8, alpha = 0.8) +
      geom_point(aes(color = organism), size = 3) +
      geom_text(
                data = slope_data |> dplyr::filter(metric == "Observed\n(total)"),
                aes(label = paste0("n analyzed = ", scales::comma(analyzed)),
                    color = organism),
                hjust = 1, nudge_x = -0.08, size = 2.5, show.legend = FALSE) +
      ggrepel::geom_text_repel(
                data = slope_data |> dplyr::filter(metric == "Full-model\ngamma"),
                aes(label = organism, color = organism),
                hjust = 0, nudge_x = 0.08, direction = "y",
                size = 3, show.legend = FALSE, segment.size = 0.3,
                max.overlaps = Inf) +
      scale_y_log10(labels = scales::comma) +
      labs(title = "Total Amazonian Gamma Richness",
           subtitle = "Log scale; observed = total species in dataset (pre-filter). Shaded band = lower/upper bound.",
           x = NULL, y = "Species (log scale)",
           color = "Organism", fill = "Organism") +
      theme_minimal() +
      theme(legend.position = "none",
            panel.grid.minor = element_blank())
    print(p)
    dev.off()
    message("  Wrote: slope_observed_vs_gamma_", timestamp, ".pdf")
  }

  # --- Bar chart: gamma expansion ratio (gamma / observed) ---
  # Observed baseline is the pre-filter total species; each bar is annotated
  # with the observed sample size (and analyzed portion).
  ratio_plot_data <- organism_summary |>
    dplyr::filter(!is.na(gamma_best) & !is.na(pre_species) & pre_species > 0) |>
    dplyr::mutate(gamma_ratio_total = round(gamma_best / pre_species, 2)) |>
    dplyr::arrange(gamma_ratio_total)

  if (nrow(ratio_plot_data) > 0) {
    pdf(file.path(summary_dir, paste0("gamma_expansion_ratio_", timestamp, ".pdf")),
        width = 10, height = 6)

    p <- ggplot(ratio_plot_data,
                aes(y = reorder(organism, gamma_ratio_total), x = gamma_ratio_total)) +
      geom_col(width = 0.6) +
      geom_vline(xintercept = 1, linetype = "dashed", color = "red") +
      geom_text(aes(label = paste0(gamma_ratio_total, "x  (n obs = ",
                                   scales::comma(pre_species),
                                   "; analyzed = ", scales::comma(observed_spp), ")")),
                hjust = -0.05, size = 2.8) +
      scale_x_continuous(expand = expansion(mult = c(0, 0.35))) +
      labs(title = "Gamma Expansion Ratio (Estimated / Observed)",
           subtitle = "Observed = total species recorded (pre-filter). Red dashed line = 1x (no change).",
           y = NULL, x = "Gamma / Observed species") +
      theme_minimal()
    print(p)
    dev.off()
    message("  Wrote: gamma_expansion_ratio_", timestamp, ".pdf")
  }
}


# --- 6. Spatial extent overlap figure ----------------------------------------

message("\n--- Generating spatial extent overlap figure ---")
plot_spatial_extent_overlap(include = if (!is.null(include)) include else NULL,
                            group = if (!is.null(group)) group else NULL,
                            output_dir = output_dir, save_pdf = TRUE)


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
