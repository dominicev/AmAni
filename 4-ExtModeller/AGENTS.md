# 4-ExtModeller — Project Memory

## Overview

This project simulates species extinction from Amazon deforestation. It uses
patch-size frequency distributions (from `3-ExtPredictor`) and simulated animal
occurrence data (from `2-ModellingAnimals`) to estimate how many species and
populations are lost under different deforestation scenarios.

Two deforestation placement models are implemented:
- **Neutral** (`01_simulate_neutral.R`): grid cells are deforested uniformly at random
- **Clustered** (`02_simulate_clustered.R`): grid cells are deforested with probability
  proportional to observed deforestation intensity (from RAISG patch centroids)

## Folder Structure

```
4-ExtModeller/
├── scripts/
│   ├── 01_simulate_neutral.R         # Spatially random deforestation
│   ├── 02_simulate_clustered.R       # Spatially clustered deforestation
│   └── utils/
│       └── simulation_functions.R    # Shared functions (supports both models)
├── output/
│   ├── neutral_Freq*.RData           # Cached neutral simulation results
│   └── clustered_Freq*.RData         # Cached clustered simulation results
├── archive/
│   └── .RData, .RDataTmp*, .Rhistory # Old session files
└── AGENTS.md
```

All scripts set their working directory to the **project root** (`4-ExtModeller/`),
so paths are relative from there.

## External Data Dependencies

| Source | Path (relative to AmAni/) | What |
|--------|--------------------------|------|
| `3-ExtPredictor/output/` | Model parameter CSVs, patch frequency tables | Patch size distributions |
| `3-ExtPredictor/output/deforestTileOutStats/` | `w*.patch_metrics.csv` | Patch centroids (for clustered weights) |
| `2-ModellingAnimals/madre de dios analysis results/` | `simulatedOccurencesForAnimal*.rds` | Simulated animal occurrences (sf) |
| `1-FilteringPlantData/` | `species_counts4_filtered.rds` | Used to construct the analysis grid |

## Key Parameters

- **Grid cell size**: 22,360 m (~22.36 km), giving ~500 km² cells in EPSG:3857
- **Total deforested area**: 454,509.4 km² (from `terra::cellSize` on `def_00_20.tif`)
- **Species**: 2,940 unique species, 27,948 population records across ~10,800 occupied grid cells
- **Single-grid endemics**: 1,051 species (35.7% of all species)
- **Simulation replicates**: 100 per scenario

## Frequency Table Scenarios

| Name | Filtering | Patches | Total area (km²) |
|------|-----------|---------|-------------------|
| `freq_over` (conservative) | Boundary-filtered + proximity dedup | 10,946 | 296,476 |
| `freq_boundary` (middle) | Boundary-filtered only | 9,376 | 256,606 |
| `freq_under` (liberal) | Centroid dedup only | 54,014 | 1,755,986 |

Note: `freq_under` total area (1.76M km²) far exceeds actual deforestation (454K km²)
due to extensive double-counting. The frequency tables have not yet been scaled to
match the actual deforested area.

## Simulation Design

### Large patches (≥ grid cell size)
- Convert each bin's total area to whole-grid-cell equivalents
- Sample that many grid cells (uniform or weighted) and remove all species in them
- Track sub-grid remainder area, fold back into small-patch model

### Small patches (< grid cell size)
- For each bin, sample `count` grid cells and compute a fractional removal
- Accumulate survival probability per grid: `p_survive = product((1 - perc)^n_hits)`
- Apply stochastic row-level removal in a single pass

### Clustered weighting (02_simulate_clustered.R)
- Spatial join of patch centroids to grid cells
- Weight = total deforested area within each grid cell + small floor
- Passed as `grid_weights` parameter through the function chain

## Results (Neutral Model, as of 2026-03-20)

| Scenario | Median spp lost | % of 2,940 | Median pops lost |
|----------|----------------|------------|-----------------|
| FreqOver | 58 | 2.0% | 1,516 |
| FreqBoundary | 50 | 1.7% | 1,305 |
| FreqUnder | 353 | 12.0% | 8,245 |

Extinctions are overwhelmingly single-grid endemics (87–100%). Multi-grid
extinctions only occur in the liberal scenario (max range: 6 grids).

## Known Issues

1. **Frequency tables not area-scaled**: Total patch area in freq tables doesn't
   match actual deforested area (454,509 km²). `freq_under` overshoots by ~4×.
2. **Deforestation is spatially independent across patches**: Each patch placement
   is independent; real deforestation expands contiguously from existing clearings.
3. **EPSG:3857 grid distortion**: Grid cells are equal in projected coordinates but
   not in true area at Amazon latitudes. `st_area()` accounts for this, but the
   22.36 km grid spacing is approximate.

## Pending Tasks

- [ ] Scale frequency tables to match 454,509 km² total deforestation
- [ ] Run clustered model and compare to neutral
- [ ] Compare which species go extinct under clustered vs neutral
- [ ] Experiment with alternative weighting schemes
- [ ] Simulate range shifts (climate change scenarios)
