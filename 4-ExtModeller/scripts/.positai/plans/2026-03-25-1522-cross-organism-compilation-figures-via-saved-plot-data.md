# Plan: Cross-Organism Compilation Figures via Saved Plot Data

## Approach

Saving the underlying data as RDS files alongside each PDF is the right call. This gives several advantages over stitching PDFs:
- Consistent theming/scales across organisms in compilation figures
- Ability to use shared color scales or axis limits for fair comparison
- Smaller, faster-loading data (we can drop geometry for non-map plots)
- Reproducible — regenerate any figure without re-running the pipeline

## What data to save (per organism)

Each organism's run will save a standardized `plot_data/` subdirectory inside its `animal_model/` output folder. File names are fixed (no timestamps) so the summary script can find them reliably.

| RDS file | Source script | Contents | Notes |
|----------|--------------|----------|-------|
| `plotdata_animal_richness_map.rds` | rawEndemismCalculator.R | `animal_grid_data_clipped` (sf: grid_id, animal_richness, geometry) | Clipped to Amazon boundary |
| `plotdata_animal_endemism_map.rds` | rawEndemismCalculator.R | `animal_endemism_clipped` (sf: grid_id, animal_endemism_count, geometry) | Clipped to Amazon boundary |
| `plotdata_endemism_correlation.rds` | rawEndemismCalculator.R | `combined_endemism` (df: raw_endemism, animal_endemism_count) | Drop geometry — just two columns |
| `plotdata_scaling_factors.rds` | rawEndemismCalculator.R | `scaling_df_under2` (df: scaling_factor) | Histogram data |
| `plotdata_scaled_endemism_map.rds` | rawEndemismCalculator.R | `grid_sf` (sf: grid_id, scaled_plant_endemism, geometry) | Only cells with non-NA values |
| `plotdata_simulated_endemism_dist.rds` | rawEndemismCalculator.R | list(results, med, mn, ci_low, ci_high) | Summary stats + raw vector |
| `plotdata_richness_correlation.rds` | 02_animalAnalysis.R | `richness_joined` (df: plant_richness, morphotype_richness) | Drop geometry |
| `plotdata_plant_extent_model.rds` | 02_animalAnalysis.R | list(plant_df, ci_df, y_pred, y_predUpper) | Curve + CI data |
| `plotdata_simulated_richness_map.rds` | 02_animalAnalysis.R | `sim_richness_grid_sf` (sf: grid_id, sim_richness, geometry) | Clipped if possible |

**Plant-only plots** (rawEndemismGridplot, weightedEndemismGridplot) are identical across organisms — save once per run:

| RDS file | Source script | Contents |
|----------|--------------|----------|
| `plotdata_plant_raw_endemism.rds` | rawEndemismCalculator.R | `grid_sf_clipped` (sf: grid_id, raw_endemism, geometry) |
| `plotdata_plant_weighted_endemism.rds` | rawEndemismCalculator.R | `grid_sf_clipped` (sf: grid_id, weighted_endemism, geometry) |

## Changes to source scripts

### rawEndemismCalculator.R (6 save points)

After each existing PDF write, add one `saveRDS()` call to save the plot data to the `plot_data/` subdirectory. Save only the columns needed — drop unnecessary columns and, for non-map data, drop geometry.

1. After animal_richness_map PDF → save `animal_grid_data_clipped[, c("grid_id", "animal_richness", "geometry")]`
2. After animal_endemism_map PDF → save `animal_endemism_clipped[, c("grid_id", "animal_endemism_count", "geometry")]`
3. After plant_vs_insect.corr PDF → save `sf::st_drop_geometry(combined_endemism)[, c("raw_endemism", "animal_endemism_count")]`
4. After distScalingFactors PDF → save `scaling_df_under2`
5. After animalEnd.Amaz.medianLogNorm PDF → save `grid_sf[!is.na(grid_sf$scaled_plant_endemism), c("grid_id", "scaled_plant_endemism", "geometry")]`
6. After dist.Sim.endemSpp PDF → save `list(results = scaled_results, med = med, mn = mn, ci_low = ci_low, ci_high = ci_high)`

Also save plant endemism maps (once):
7. After the grid_sf_clipped construction → save raw_endemism and weighted_endemism subsets

### 02_animalAnalysis.R (3 save points)

1. After plantPatExtentModel PDF → save `list(plant_df = plant_df, ci_df = ci_df, y_pred = y_pred, y_predUpper = y_predUpper)`
2. After plant_vs_insect.RichnessCorr PDF → save `sf::st_drop_geometry(richness_joined)[, c("plant_richness", "morphotype_richness")]`
3. After simulated_richness_map PDF → save `sim_richness_grid_sf[, c("grid_id", "sim_richness", "geometry")]`

## AmAni.summarize.R — new compilation figure section

Add a new section (section 5) that:

1. Discovers all organisms with `plot_data/` directories
2. Loads each RDS, tagging with organism name
3. Generates these compilation figures:

### Compilation figures to generate

| Figure | Type | Layout |
|--------|------|--------|
| **Animal richness maps** | sf map | Faceted by organism (facet_wrap) |
| **Animal endemism maps** | sf map | Faceted by organism |
| **Endemism correlation** (animal ~ plant) | Scatter + lm | Faceted by organism, with R² annotation |
| **Richness correlation** (animal ~ plant) | Scatter + lm | Faceted by organism, with R² annotation |
| **Scaling factor distributions** | Histogram | Faceted by organism |
| **Simulated endemic species distribution** | Histogram + vlines | Faceted by organism |
| **Plant extent model with CIs** | Line + ribbon | Faceted by organism |
| **Scaled endemism maps** | sf map | Faceted by organism |
| **Simulated richness maps** | sf map | Faceted by organism |
| **Plant endemism maps** | sf map | Single panel (shared across organisms) — include once at top |

Each figure is saved as a single multi-panel PDF in `cross_organism_summary/`.

### Note on map faceting

For sf maps, `facet_wrap` works directly with `geom_sf` when all organism data is row-bound into one sf object with an `organism` column. We'll use `scales = "fixed"` for fill to allow cross-organism comparison, and `coord_sf()` to ensure proper projection handling per facet.

## Implementation order

1. Add `plot_data/` directory creation to both source scripts (one line each)
2. Add saveRDS calls in rawEndemismCalculator.R (7 saves)
3. Add saveRDS calls in 02_animalAnalysis.R (3 saves)
4. Rewrite section 3 of AmAni.summarize.R (replace PDF copying with data loading)
5. Add section 5 to AmAni.summarize.R (compilation figures)

Total: ~13 small edits to source scripts + one major addition to AmAni.summarize.R.
