# Plan: Adding Vertebrate Baseline Data Alongside Plants

## Context

The current pipeline uses plant data as the sole biodiversity baseline. Plant
richness and endemism from stages 0/1 feed into stage 2, where they serve three
purposes:

1. **Range-size distribution**: Plant species' spatial extents (# cells) are
   fitted to a lognormal, then rescaled to simulate target animal range sizes.
2. **Spatial weights**: Plant occurrence frequency per cell determines where
   simulated species originate and expand.
3. **Correlation validation**: Plant richness/endemism vs. animal
   richness/endemism validates that the proxy predicts the target taxon.

The proposal: add a vertebrate-based baseline alongside plants so that stage 2
can use both (or either) as a biodiversity proxy. The vertebrate pipeline should
be structurally parallel to stages 0 and 1 and produce outputs with an
identical schema, differing only in the taxonomic source.

## Design Decisions Needed (Before Implementation)

1. **Which vertebrate groups?** Options:
   - All Vertebrata (maximizes coverage but mixes very different ecologies)
   - Birds only (best GBIF coverage in the tropics; strongest spatial signal)
   - Birds + mammals + amphibians (comprehensive but heterogeneous)
   - Parameterized by class, so each can be run independently

   *Recommendation*: Parameterize by GBIF taxon key, so the user chooses at
   run time. Start with birds (class Aves, taxon key 212) as the initial test
   case since they have the densest Amazon coverage.

2. **Taxonomy backbone**: Stage 0 uses World Flora Online. For vertebrates,
   the GBIF backbone taxonomy is the pragmatic choice — the `species` or
   `acceptedScientificName` field from the download already carries accepted
   names. A WFO-style external matching step is unnecessary.

3. **How does stage 2 combine plant and vertebrate baselines?** Options:
   - **Side-by-side**: Run the simulation twice (once with plant weights, once
     with vertebrate weights) and report both. Simplest.
   - **Blended weights**: Combine plant + vertebrate richness per cell into a
     single weight vector (e.g., average, or R²-weighted). More ecologically
     nuanced but adds complexity.
   - **User-selected**: Let `run_animal_model()` accept a `baseline` argument
     (`"plant"`, `"vertebrate"`, `"both"`).

   *Recommendation*: Start with user-selected (default `"plant"` for backward
   compatibility). Add blended as a future option.

## Architecture

### New Directories

```
AmAni/
├── 0-CuratingPlantData/         # Existing (unchanged)
├── 0v-CuratingVertData/         # NEW — mirrors stage 0
│   ├── scripts/
│   │   └── curatingVertData.R
│   ├── output/                  # species_counts_vert.rds, etc.
│   └── AGENTS.md
├── 1-FilteringPlantData/        # Existing (unchanged)
├── 1v-FilteringVertData/        # NEW — mirrors stage 1
│   ├── scripts/
│   │   └── vertfiltering.R
│   └── output/                  # species_counts_vert_filtered.rds, etc.
├── 2-ModellingAnimals/          # Modified to accept vertebrate baseline
│   └── ...
└── AmAni.manager.R              # Modified: new functions
```

### Shared Grid

Both plant (stage 0) and vertebrate (stage 0v) pipelines must produce data on
the **identical grid** (22,360 m cells, EPSG:3857, Amazon basin boundary). The
grid is currently created inside stage 0 by:

```r
st_make_grid(amazon_3857, cellsize = 22360, what = "polygons")
```

To guarantee grid identity, the plan is to:
- Have stage 0 export the grid geometry itself (`amazon_grid.rds`)
- Have stage 0v **load** that grid rather than regenerating it
- If stage 0v runs before stage 0, it creates the grid and exports it;
  stage 0 then loads it on its next run (whichever runs first creates it)

Alternatively, extract grid creation into a shared utility script under a new
`shared/` directory.

## Implementation Steps

### Step 1: Create `0v-CuratingVertData/`

**Script: `curatingVertData.R`** — adapted from `curatingRawPlantData.R` with
these changes:

| Aspect | Stage 0 (plants) | Stage 0v (vertebrates) |
|--------|-------------------|------------------------|
| GBIF download | Plantae | Configurable (Aves, Mammalia, etc.) |
| Column selection | Same core columns (lat, lon, species, year, uncertainty) | Same |
| Quality filters | occurrenceStatus, year, coord uncertainty, duplicates | Same |
| Spatial filter | Amazon boundary (same shapefile) | Same |
| Grid | 22,360 m, EPSG:3857 | Load from stage 0 or create if absent |
| Taxonomy | WFO matching + name cleaning | GBIF backbone (`acceptedScientificName` or `species` field); simpler name cleaning (no WFO) |
| Aggregation | species_dt: grid_id × scientificName × count | Same schema |
| Endemism | weighted_endemism, raw_endemism | Same formulas |

**Outputs** (same schema as stage 0, different filenames):
- `species_counts_vert.rds` — sf: (grid_id, scientificName, count, geometry)
- `species_endemism_vert.rds` — df: (grid_id, weighted_endemism)
- `raw_species_endemism_vert.rds` — df: (grid_id, raw_endemism)

### Step 2: Create `1v-FilteringVertData/`

**Script: `vertfiltering.R`** — adapted from `plantfiltering.R`. Stage 1 is
already taxon-agnostic; the filtering logic (Q95 + 3×IQR outlier removal,
urban area exclusion, endemism recalculation) applies identically to any
taxonomic group. The main change is just input/output filenames.

**Outputs**:
- `species_counts_vert_filtered.rds` — sf: (grid_id, scientificName, count, geometry)
- `FINAL_allEndemism_vert_FullyFiltered.rds` — sf: (grid_id, raw_endemism, weighted_endemism, geometry)

### Step 3: Modify Stage 2

In `02_animalAnalysis.R` (and variants), add support for an alternative or
additional baseline. The key touchpoints:

**A. Loading baseline data** (currently lines ~139, 152):
- Add conditional loading of vertebrate files alongside plant files
- Controlled by a new `.amani_baseline` global (set by `run_animal_model()`)

**B. Computing baseline richness** (currently lines ~309-350):
- If `baseline = "vertebrate"`, use vert species_counts instead of plant
- If `baseline = "both"`, compute both and average (or R²-weighted blend)

**C. Spatial weights** (currently lines ~606-656):
- Same logic, different source data

**D. Range-size distribution fitting** (currently lines ~263-305, 520-527):
- Fit the vertebrate range-size distribution (same lognormal approach)
- Use it for the scaling factor calculation

**E. Correlation analyses** (currently lines ~1601-1707):
- Run plant-vs-animal AND vertebrate-vs-animal correlations
- Report both R² values
- If `baseline = "both"`, optionally use R² to weight the blend

**F. Export**: Same output schema (simulatedOccurencesForAnimal.rds) — the
downstream stages (3, 4) don't need to know which baseline was used.

### Step 4: Update `AmAni.manager.R`

- Add `run_vert_curation()` and `run_vert_filtering()` functions (parallel to
  the existing stage 0/1 setup, but for vertebrate data)
- Add `baseline` parameter to `run_animal_model()`:
  ```r
  run_animal_model("cockroaches", baseline = "plant")        # default
  run_animal_model("cockroaches", baseline = "vertebrate")   # vert only
  run_animal_model("cockroaches", baseline = "both")         # blended
  ```
- Add `stage_dirs$curate_vert` and `stage_dirs$filter_vert`
- Update `check_status()` to show vertebrate pipeline status
- Update `verify_upstream()` for stage 2 to check vertebrate files when
  `baseline != "plant"`

### Step 5: Update `AmAni.run.R`

Add vertebrate curation + filtering calls (run once), then use the `baseline`
argument in organism model runs as desired:

```r
# --- Vertebrate baseline (run once) ---
run_vert_curation(taxon = "Aves", data_file = "data/aves_gbif.csv")
run_vert_filtering()

# --- Then use it ---
run_animal_model("cockroaches", baseline = "both", ...)
```

## Scope and Phasing

**Phase 1** (minimum viable):
- Create 0v and 1v with scripts adapted from stages 0 and 1
- Add `baseline = "plant" | "vertebrate"` to stage 2 (no blending yet)
- Run for birds (Aves) as proof of concept

**Phase 2** (if results warrant):
- Add `baseline = "both"` with R²-weighted blending
- Parameterize 0v to run multiple vertebrate classes
- Compare plant-only vs vertebrate-only vs blended results across organisms

## Risks and Considerations

- **GBIF vertebrate coverage in the Amazon** is much sparser than plants for
  some groups (especially non-bird vertebrates). Bird data is the safest bet.
- **Ecological validity**: Whether vertebrate diversity predicts invertebrate
  diversity better than plants is an empirical question — the correlation
  analyses in stage 2 will answer this directly.
- **Grid identity**: If the Amazon boundary shapefile or grid parameters change,
  all stages must be re-run. The shared grid approach mitigates this.
- **Runtime**: Stage 0v should be fast (vertebrate data is much smaller than
  plant data). Stage 1v is trivial. The bulk of the work is in stage 2.
