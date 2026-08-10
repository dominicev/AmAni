# 3-ExtPredictor — Project Memory

## Overview

This project analyzes deforestation patch geometry across the Amazon basin using
rasterized deforestation data (2001–2020) from RAISG. The goal is to extract
patch-level metrics (area, perimeter, compactness, centroid) for large deforested
patches, which feed into downstream extinction risk modeling.

## Folder Structure

```
3-ExtPredictor/
├── data/
│   ├── def_00_20.tif                   # Binary deforestation raster
│   └── raisg-shapefiles/               # RAISG source shapefiles
├── scripts/
│   ├── 01_prep_deforestation.R         # Data prep + original grid run + resume block
│   ├── 02_run_deterministic.R          # Standalone headless runner for resume
│   ├── 03_assess_patches.R            # Downstream patch analysis
│   └── utils/
│       └── calcPatchStat.sub.R         # Core subroutine: extract_patch_metrics()
├── output/
│   ├── deforestTileOutStats/           # Deterministic run outputs (3,630 windows)
│   ├── deforestTileOutStats2/          # Single random-method test (not main analysis)
│   ├── *_model_parameters_*.csv        # Model fit outputs from 03_assess_patches.R
│   └── patch_area_frequency_*.csv      # Frequency tables from 03_assess_patches.R
├── archive/                            # Deprecated scripts, logs, one-offs
├── docs/
│   ├── TROUBLESHOOTING.md
│   └── TEMP_FILE_INFO.md
└── AGENTS.md
```

All scripts set their working directory to the **project root** (not the scripts/ folder),
so paths are relative to the root: `data/def_00_20.tif`, `output/deforestTileOutStats/`, etc.

## Data Sources

### Historical deforestation (primary analysis)
- **File**: `data/def_00_20.tif` — binary raster (1 = deforested, NA = background)
- **Origin**: Derived from `data/raisg-shapefiles/deforestacion_2001_2020_total/deforestacion2001-2020.tif`
  by thresholding (values >= 2000 → 1). Source: https://www.raisg.org/en/maps/
- **Dimensions**: 125,279 rows × 150,347 cols (~30 m resolution)
- **CRS**: Geographic (lon/lat); reprojected to EPSG:5880 (equal-area) for area calculations

### Future deforestation projections (exploratory only)
- **Directory**: `data/raisg-shapefiles/deforestacion_futura_2025/`
- **File**: `deforestacion_futura_2025.shp`
- **Key variables**: `esc_tend` (trend scenario), `esc_opti` (optimistic), `esc_pesc` (pessimistic)
- Explored in the first section of `scripts/01_prep_deforestation.R` but not used in the patch extraction pipeline

## Workflow

### Step 1: Prepare and explore (`scripts/01_prep_deforestation.R`)
- Loads and visualizes future deforestation shapefile (exploratory)
- Loads `data/def_00_20.tif` binary raster
- Contains a single random-method test call to `extract_patch_metrics()` (output in `output/deforestTileOutStats2/`)
- Contains the original deterministic grid run (windows 1–2430 processed on a 126 GB Ubuntu machine)
- Contains a resume block for remaining windows (2431–3630 + windows 303 and 533)

### Step 2: Run deterministic patch extraction (`scripts/02_run_deterministic.R`)
- Standalone headless-safe script for running the resume block via `Rscript`
- Cross-platform (Linux `/proc/meminfo` + Windows `systeminfo` for RAM detection)
- Uses absolute paths and sets worker `setwd()` to handle PSOCK cluster workers correctly

### Step 3: Assess patch results (`scripts/03_assess_patches.R`)
- Imports all `.csv` outputs from `output/deforestTileOutStats/`
- Deduplicates patches: centroid-based (exact match) then boundary-based (core-zone)
- Truncation diagnostic on largest patches
- Fits power-law and lognormal distributions to patch areas (MLE via `poweRlaw`)
- Formal model comparison using Vuong likelihood ratio test
- Produces outputs for **three filtering regimes** (matching stage 4 expectations):
  - `boundary_filtered`: core-zone dedup (principled, main analysis)
  - `underfiltered`: centroid-only dedup (retains near-duplicates from window overlap)
  - `overfiltered`: boundary-filtered + 10 km² minimum area threshold (conservative)
- Cleans stale output files before writing new ones

### Subroutine
- **`scripts/utils/calcPatchStat.sub.R`** — defines `extract_patch_metrics()`. Supports `"random"` and `"deterministic"` sampling methods. For each window: crops the raster, applies erosion→dilation morphology as a mask, labels patches, filters by minimum area (default 5 km² = 5e6 m²), computes metrics, writes CSV + GeoTIFF + GeoPackage + `.done` flag.
- **Parameters**: `output_id` (file naming prefix, defaults to `seed` for backward compatibility), `seed` (RNG seed, random method only), `method` ("random" or "deterministic").

## Deterministic Grid Configuration

- Grid tiling computed for **4500×4500** cell windows with ~50% overlap
- Actual extraction uses **4000×4000** cell windows (known mismatch; see below)
- Grid: 55 rows × 66 cols = **3,630 total windows**
- Step size: ~2,237 rows × ~2,244 cols
- Actual overlap: ~44% (intended was 50%)
- Edge strips (last 500 rows/cols) contain **zero data** — no patches missed

## Output Directory

- `output/deforestTileOutStats/` — deterministic run outputs
  - `w{id}.done` — completion flag (exists even if no patches found)
  - `w{id}.patch_metrics.csv` — patch metrics table (only if patches met threshold)
  - `w{id}.patches_filtered.tif` — filtered patch raster
  - `w{id}.patch_polygons.gpkg` — patch polygons
  - Also contains ~487 non-`w`-prefixed files from an earlier random sampling run (legacy)
- `output/deforestTileOutStats2/` — output from single random-method test call (5 files, not part of main analysis)

## Completion Status (as of 2026-03-20)

- **3,630 / 3,630** deterministic windows complete (all `.done` flags present)
- **1,005** windows produced actual output files (CSV + TIF + GPKG)
- **2,625** windows had no patches meeting the 50 km² threshold (`.done` only)

### Coverage verification
Spatial coverage was verified by checking all deforested sample points (0.25° grid) against known window extents (from TIF files + reconstructed resume-run grid). Result: **542 of 568 (95.4%)** deforested points have confirmed coverage. The 26 uncovered points are **all in original-run territory** (grid columns 6–43) where 1,720 windows found no patches ≥ 50 km² and thus produced no TIF files to verify position. The original run processed 2,428 windows ≈ 44 full columns of 55 rows, and the resume run processed the remaining 1,202. Together they cover all 3,630 grid positions. **The grid is the same — only the ID-to-position mapping was lost, not the spatial coverage.**

## Known Issues

1. **Window size mismatch**: Grid spacing assumes 4500-cell windows but extraction uses 4000. This reduces overlap from 50% to 44% but does NOT create coverage gaps (confirmed: edge strips are all NA).
2. **Legacy output files**: ~487 non-`w`-prefixed files in `output/deforestTileOutStats/` from prior random sampling run. Should be cleaned up or moved to archive.
3. **Window-ID mapping discrepancy**: The original run (windows 1–2430, on Ubuntu) used a code version whose `expand.grid()` ordering (or grid construction) differs from the current `scripts/01_prep_deforestation.R`. As a result, the `window_id` embedded in CSV/TIF filenames for windows 1–2430 does **not** correspond to the grid positions reconstructed by the current `expand.grid(start_row, start_col)` call. The resume run (windows 2431–3630) used the current code and has correct ID-to-position mapping. **Workaround in `03_assess_patches.R`**: instead of reconstructing the grid from `expand.grid()`, the script reads each window's actual geographic extent from its output `.tif` file (`w{id}.patches_filtered.tif`), which is immune to the ID ordering issue. This is critical for the boundary-based deduplication to work correctly.
4. **Morphology naming**: The `apply_morph` parameter in `extract_patch_metrics()` controls only an *additional* dilation step (patch merging). The initial erosion→dilation noise-filter mask always runs unconditionally. The parameter name is misleading but the behavior is intentional.

## Resource Settings (conservative, to avoid crashes)

- Core usage: 40% of physical cores
- RAM per worker: pre-computed in master process as 70% of (total_ram / n_workers)

## Pending Tasks

- [x] Complete the deterministic run on Ubuntu (all 3,630 windows done)
- [x] Reorganize folder structure
- [x] Refactor `03_assess_patches.R` — deduplicated fitting code, added Vuong test
- [x] Reconcile filtering regime naming with stage 4 (boundary_filtered, underfiltered, overfiltered)
- [ ] Run `scripts/03_assess_patches.R` on the complete dataset (re-run needed after refactor)
- [ ] Clean up legacy (non-`w`-prefixed) files from `output/deforestTileOutStats/`

## Changelog (2026-03-28)

- **`calcPatchStat.sub.R`**: Fixed `ifelse()` bug that returned only the first list element when `returnQ = TRUE` (now uses `if/else`). Added `output_id` parameter for explicit file naming, decoupling it from `seed` (backward compatible — defaults to `seed` if not provided).
- **`01_prep_deforestation.R`**: Fixed `setwd()` to work via both RStudio and `Rscript` (cross-platform `commandArgs()` fallback). Fixed skip-check filename typos (missing dots) and uncommented the `return(NULL)` so completed windows are properly skipped. Pre-computes worker memory limit in master process (was redundantly detected per-worker). Updated `extract_patch_metrics()` calls to use `output_id` instead of `seed` for deterministic runs.
- **`02_run_deterministic.R`**: Updated `extract_patch_metrics()` call to use `output_id`.
- **`03_assess_patches.R`**: Refactored duplicated fitting/export code into `fit_and_export_models()` function. Added formal Vuong likelihood ratio test (power-law vs lognormal). Renamed "centroid_only" → "underfiltered" and added "overfiltered" variant (boundary-filtered + 10 km² minimum) to match stage 4 expectations. Added stale output file cleanup before writing new results.

## Archived Files

- `archive/calcPatchStat_optimized.sub.R` — older, less capable version of the subroutine (random-only, in-memory morphology, no error handling). Superseded by `scripts/utils/calcPatchStat.sub.R`.
- `archive/AmazonHabitatLoss_prep.R` — v1 of the prep script
- `archive/FindingAreaOfAmazon.R`, `archive/manuallyCheckingSomeClusters.R`, `archive/reconstructingRasterFromSS.R`, `archive/test_single_window.R`, `archive/cleanup_temp.R`, `archive/restart_processing.R` — one-off exploratory/utility scripts
- `archive/AssessDefPatches.log`, `archive/nohup.out` — log files from prior runs
