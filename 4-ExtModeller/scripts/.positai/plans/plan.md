# Plan: Multi-panel composite figures + save intermediate objects

## Context

`AmAni.summarize.R` currently produces **4 individual PDF plots** (one per page)
and collects per-organism PDFs from the animal model stage. We want to:

1. Add **multi-panel composite figures** that put related information side by side
   for easier visual comparison across organisms.
2. **Save intermediate objects to disk** during the Stage 2 pipeline so that
   `AmAni.summarize.R` (and future analyses) can load structured data instead of
   parsing timestamped text files.

---

## Part A: Save intermediate objects during Stage 2 pipeline

### What to save

At the end of `ModelAnimal.R` (after both sourced scripts finish), save a single
RDS file per organism containing key model outputs. This supplements (not replaces)
the existing timestamped PDFs and text files.

**File:** `output/{organism}/animal_model/model_summary.rds`

**Contents** (a named list):

| Key | Object | Source script |
|-----|--------|--------------|
| `gamma_richness` | integer — best estimate | `02_animalAnalysis.R` |
| `gamma_upper` | integer — upper bound | `02_animalAnalysis.R` |
| `gamma_lower` | integer — lower bound | `02_animalAnalysis.R` |
| `scaling_factor` | numeric — animal_max / plant_fit_max | `02_animalAnalysis.R` |
| `n_morphotypes` | integer — unique species/morphotypes | entry script |
| `n_grid_cells` | integer — occupied grid cells | entry script |
| `endemism_r2` | numeric — R² of animal~plant endemism | `rawEndemismCalculator.R` |
| `richness_r2` | numeric — R² of animal~plant richness | `02_animalAnalysis.R` |
| `scaled_endemism_median` | numeric — median scaled endemic spp | `rawEndemismCalculator.R` |
| `scaled_endemism_ci` | numeric(2) — 5th/95th percentile | `rawEndemismCalculator.R` |
| `morphotype_distribution` | data.frame — spatial extent dist | `02_animalAnalysis.R` |

### Where to add the saving code

At the end of `ModelAnimal.R`, after `source("scripts/02_animalAnalysis.R")`,
add a block that assembles these objects from the global environment and calls
`saveRDS()`. Use `tryCatch` for each so a missing object (e.g., when endemism is
skipped) produces `NA` rather than an error.

---

## Part B: Enhance `AmAni.summarize.R`

### B1. Load structured data instead of parsing text files

Replace the current gamma richness text-file parser (section 2) with code that
reads `model_summary.rds` when available, falling back to text parsing for
organisms that ran before the RDS was added.

### B2. New composite figures using `patchwork`

Since the per-organism plots are saved as PDFs (not ggplot objects), the composite
figures will be **cross-organism comparison panels** built from the tabular data.

#### Figure 1: Extinction comparison dashboard

A single figure with 4 panels arranged via `patchwork`:
- **(a)** Gamma richness by organism (bar + error bars)
- **(b)** Species lost (FreqBoundary, dodged by model)
- **(c)** Clustered vs. neutral % reduction by organism
- **(d)** Population loss (FreqBoundary, dodged by model)

All 4 panels share the same x-axis (organism). Uses `patchwork::plot_layout`
to arrange as a 2x2 grid with a collected legend.

**Output:** `composite_extinction_dashboard_{timestamp}.pdf` (landscape, ~14x10)

#### Figure 2: Gamma richness + model validation

A 1x2 panel:
- **(a)** Gamma richness by organism (bar + CI)
- **(b)** Endemism R² and richness R² per organism (paired dot plot)

This requires `model_summary.rds` data. Only organisms with the RDS are included.

**Output:** `composite_gamma_validation_{timestamp}.pdf` (~12x5)

#### Figure 3: Scenario sensitivity heatmap

A heatmap (`geom_tile`) of median species lost with organisms on y-axis, scenarios
on x-axis, faceted by spatial model. Quick overview of how different frequency
assumptions affect each organism.

**Output:** `composite_scenario_heatmap_{timestamp}.pdf` (~10x6)

### B3. Keep existing individual plots

The existing standalone PDFs remain — the composite figures are additions, not
replacements.

---

## Implementation steps

1. **Edit `ModelAnimal.R`** — Add the `saveRDS()` block at the end of the
   `modelAnimal()` function, after the two `source()` calls.

2. **Edit `rawEndemismCalculator.R`** — Assign `endemism_r2` and
   `scaled_endemism_*` to global env at the appropriate points so `ModelAnimal.R`
   can collect them.

3. **Edit `02_animalAnalysis.R`** — Assign `richness_r2` to global env after the
   richness model is fit.

4. **Edit `AmAni.summarize.R`** — Add `library(patchwork)`, load RDS data, build
   composite figures 1-3, write to PDF.

---

## Dependencies

- `patchwork` package (for multi-panel layout). User may need to install.

---

## Not in scope

- Re-rendering per-organism maps as ggplot panels (would require saving ggplot
  objects during the pipeline, which is a larger refactor).
- Changes to Stage 4 scripts.
