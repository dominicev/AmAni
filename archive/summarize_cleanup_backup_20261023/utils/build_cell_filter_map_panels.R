# utils/build_cell_filter_map_panels.R
# Extracted from AmAni.summarize.R. Sourced by AmAni.summarize.R — relies on
# the packages (tidyverse, sf, patchwork, ggrepel) loaded there.

# --- Two-page cell-filtering report (per organism) --------------------------
#
# Page 1: the before/after filtering stats figure plus the cell-filter report
#         maps (animal vs plant completeness, richness, occurrences, and the
#         plant-minus-animal completeness difference).
# Page 2: plant vs animal Chao1 completeness in overlapping grid cells.
#
# Maps are rebuilt from plotdata_cell_filter_report.rds (saved by
# 02_animalAnalysis.R alongside the cell_filter_report_maps_*.pdf), so they
# match that PDF but can be arranged in a single grid.

build_cell_filter_map_panels <- function(cfr) {

  animal_cells <- cfr$animal
  plant_cells  <- cfr$plant
  delta_cells  <- cfr$delta

  outline_scale <- scale_color_manual(
    values = c(`TRUE` = "black", `FALSE` = "grey70"),
    labels = c(`TRUE` = "Retained", `FALSE` = "Filtered out"),
    name = "Animal cell"
  )

  map_theme <- theme_minimal(base_size = 9) +
    theme(plot.title = element_text(face = "bold", size = 9),
          legend.key.size = grid::unit(0.4, "cm"))

  map_pair <- function(an_var, pl_var, title, fill_name, trans = "identity") {
    list(
      ggplot(animal_cells) +
        geom_sf(aes(fill = .data[[an_var]], color = retained), linewidth = 0.4) +
        scale_fill_viridis_c(name = fill_name, trans = trans) +
        outline_scale +
        labs(title = paste("Animals:", title)) +
        map_theme,
      ggplot(plant_cells) +
        geom_sf(aes(fill = .data[[pl_var]]), color = NA) +
        scale_fill_viridis_c(name = fill_name, trans = trans) +
        labs(title = paste("Plants:", title)) +
        map_theme
    )
  }

  p_delta <- ggplot(delta_cells) +
    geom_sf(aes(fill = delta_completeness, color = retained), linewidth = 0.4) +
    scale_fill_gradient2(
      name = "Plant - animal\ncompleteness",
      low = "firebrick", mid = "white", high = "steelblue", midpoint = 0,
      limits = c(-1, 1) * max(abs(delta_cells$delta_completeness), na.rm = TRUE)) +
    outline_scale +
    labs(title = "Delta Chao1 completeness (plant - animal), animal grid cells") +
    map_theme

  c(map_pair("animal_completeness", "plant_completeness",
             "Chao1 completeness", "Completeness"),
    map_pair("animal_richness", "plant_richness",
             "species richness", "Richness"),
    map_pair("animal_occurrences", "plant_occurrences",
             "occurrences (log10)", "Occurrences", trans = "log10"),
    list(p_delta))
}
