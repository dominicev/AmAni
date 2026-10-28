# Option 3: Two-Phase Deforestation Pattern Separation

## Problem
`SEFD_rep_parallel` sends the full species×grid table (`sim_tbl`) to every
parallel worker. For large organisms (lepidoptera: 6.7M rows / 4.9 GB;
flies: 113 MB on disk), this causes OOM crashes even with few cores.

## Key Insight
The simulation has two logically separate phases, and only one needs species data:

- **Phase A — Generate deforestation pattern** (needs only grid IDs + weights):
  `run_large_patch_deforestation` samples which grids to fully remove;
  `run_small_patch_thinning` computes `grid_survive_prob` per remaining grid.
  These are tiny objects (~100 KB).

- **Phase B — Assess species impact** (needs the full species×grid table):
  Subset `sim_tbl` by surviving grids, apply row-level thinning,
  compute per-species range reduction and threshold counts.

## Solution
1. **Workers** (parallel, many cores, tiny memory): generate deforestation
   patterns — just `grids_to_deforest` (integer vector) and
   `grid_survive_prob` (named numeric vector). A few KB per replicate.

2. **Parent** (sequential, single-threaded): apply each pattern to the full
   `sim_tbl`. Each application takes ~1-2 seconds (subsetting + counting).

Workers need essentially zero large data. You can safely use all cores.

## Implementation
The replacement code is in `alt.simulationFunction.R` in this directory.
It is a drop-in replacement: same function signature, same return value.

To activate:
```r
source("4-ExtModeller/scripts/utils/alt.simulationFunction.R")
```

## Tradeoffs
- Phase 2 is sequential, so total wall time may be similar or slightly longer
  than the current approach on machines with enough RAM.
- On RAM-constrained machines, this is dramatically better because workers
  use negligible memory.
- Can be further optimized by parallelizing Phase 2 with mclapply (each
  apply step is independent and read-only on sim_tbl).
