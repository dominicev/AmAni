# AmAni.summarize.R — Cross-organism summary of extinction model results
#
# Entry points (called from AmAni.run.R after sourcing this file):
#   summarize_all()        cross-organism summary (table + figures)
#   summarize_organism()   single-organism compilation (table + figure)
#
# Reads summary CSVs, gamma richness estimates, and saved plot data from each
# organism's output directory and produces a combined summary table (CSV),
# comparison figures across organisms (PDF), and compilation figures
# regenerated from saved plot data (PDF) via the entry points. The second
# half builds the publication figures (prototype) — opt-in, see the flag at
# the bottom of this file.
#
# All function definitions live in utils/. The working directory is assumed
# to be the project root.
#
#   func_core.R           low-level: organism labels/groups, report parsing
#   func_report_plots.R   low-level: per-organism report plot builders
#   summarize_all.R       summarize_all()
#   summarize_organism.R  summarize_organism()
#   standalone_figs.R     plot_spatial_extent_overlap(), compile_workflow_diagram()
#   pub_0_settings.R      publication figures: PUB_SETTINGS
#   pub_1_collect.R       publication figures: data collectors
#   pub_2_style.R         publication figures: theme, palettes, map layers
#   pub_3_fig_sampling.R  publication figures: Figure 1 (sampling & filtering)
#   pub_4_fig_gamma*.R    publication figures: Figure 2 (gamma richness treemaps)
#   pub_5_fig_extinction.R publication figures: Figure 3 (extinction and IUCN)
#   pub_6_make_figures.R  make_publication_figures(): draws and saves all of them

library(tidyverse)
library(sf)
library(patchwork)
library(ggrepel)
library(rphylopic)

source("utils/func_core.R")
source("utils/func_report_plots.R")
source("utils/summarize_all.R")
source("utils/summarize_organism.R")
source("utils/standalone_figs.R")

source("utils/pub_0_settings.R")
source("utils/pub_1_collect.R")
source("utils/pub_2_style.R")
source("utils/pub_3_fig_sampling.R")
source("utils/pub_4_fig_gamma.R")
source("utils/pub_4_fig_gamma_uncertainty.R")
source("utils/pub_5_fig_extinction.R")
source("utils/pub_6_make_figures.R")


# =============================================================================
# PUBLICATION FIGURES (prototype)
# =============================================================================
# Rendering is opt-in: sourcing this file (e.g. from AmAni.run.R) only loads
# the functions. To draw the figures, run the block below, or set
# AMANI_PUB_FIGS=1 before sourcing so a full pipeline run regenerates them.

if (identical(Sys.getenv("AMANI_PUB_FIGS"), "1")) {

  pub <- collect_pub_data(include = PUB_SETTINGS$include,
                          group = PUB_SETTINGS$group,
                          output_dir = PUB_SETTINGS$output_dir)

  # --- Figure 1: Sampling and filtering ----------------------------------------
  # a. Records / cells / species before vs after (dumbbell, log scale)
  # b. % retained (records, cells, species)
  # c. Per-cell completeness & richness distributions before vs after
  # d. (optional) well-sampled vs super-diverse cells / spatial footprint

  retained_cells <- collect_retained_cells()
  plant_completeness_raw <- collect_plant_completeness_raw()

  fig_retained_cells(retained_cells, plant_completeness_raw)
  fig_retained_cells_combined(retained_cells, plant_completeness_raw)
  # Same map with numbered site labels: top cell per country + next 8 overall
  fig_retained_cells_combined(retained_cells, plant_completeness_raw) |>
    label_retained_cells(retained_cell_stats(retained_cells, plant_completeness_raw))
  # Pass file = NULL for an interactive rgl window instead of a snapshot
  fig_retained_cells_3d(retained_cells, plant_completeness_raw,
                        file = file.path(PUB_SETTINGS$output_dir, PUB_SETTINGS$results_dir,
                                         "fig1_retained_cells_3d.png"),
                        labels = TRUE)

  fig_sampling_retention(pub)
  fig_cell_distributions(pub)
}

# --- Figure 2: Gamma richness ------------------------------------------------------
# Treemap: area ~ richness, split observed (full colour) vs unobserved (dulled).
# value_col sets the richness used (e.g. "gamma_lower", "surrogate_best").
# The organism set is picked automatically: umbrella datasets (Hymenoptera,
# Beetles, Lepidoptera) are swapped for their alternatives when those are
# richer together (asterisked), so nested datasets don't double count.

if (identical(Sys.getenv("AMANI_PUB_FIGS"), "1")) {
  gamma_treemap_organisms(pub, value_col = "gamma_lower") |> print()
  fig_gamma_treemap(pub, value_col = "gamma_lower")
}

# --- Figure 3: Extinction and IUCN ---------------------------------------------
# a. Species lost (median, IQR, range) as count and % of gamma  fig_species_lost()
# b. Sensitivity across frequency scenarios and spatial models   fig_extinction_sensitivity()
# c. IUCN composition ribbons, split by organism                 fig_iucn_stack()
# Composite: fig_extinction(pub)

if (identical(Sys.getenv("AMANI_PUB_FIGS"), "1")) {
  fig_extinction(pub)
}

# --- Build and save all publication figures -----------------------------------------
# make_publication_figures() collects the data, draws every figure above and
# writes PDF + PNG to output/publication_figures/. Called from AmAni.run.R.
