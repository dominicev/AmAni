# Troubleshooting Worker Crashes

## Current Problem
Workers are crashing with "Error in unserialize(socklist[[n]]) : error reading from connection"

This means R worker processes are dying before they can return results - they're likely running out of memory or crashing during processing.

## Diagnosis Steps

### Step 1: Test if a single window can process
```r
source("test_single_window.R")
```

**If this fails:** The problem is in the processing code itself (morphological operations are too intensive)
**If this succeeds:** The problem is specific to parallelization (memory per worker, or too many workers)

### Step 2: Check your current cluster
```r
cat("Cores in use:", length(cl), "\n")
cat("Memory per worker:", round(126/length(cl), 1), "GB\n")
```

**If > 20 cores:** You have too many workers, not enough memory each

### Step 3: Properly restart with conservative settings
```r
source("restart_processing.R")
```

This will:
- Stop the old cluster
- Create new cluster with only 12-13 cores (20% of 64)
- Give each worker ~10GB RAM
- Set window size to 6000x6000

## Known Issues

### Issue 1: Too many cores (MOST LIKELY CAUSE)
**Symptoms:** Workers crash randomly, "error reading from connection"
**Cause:** Each worker gets too little memory (~3GB with 44 cores)
**Solution:** Use only 12-16 cores maximum (see restart_processing.R)

### Issue 2: Morphological operations too intensive
**Symptoms:** Even single window test fails
**Cause:** Focal operations with large kernels on 36M cells create huge temp files
**Solutions:**
- Reduce kernel sizes in calcPatchStat.sub.R (lines 141, 165)
- Skip morphological operations entirely (set apply_morph = FALSE)
- Further reduce window size to 4000x4000

### Issue 3: Temp directory fills up
**Symptoms:** Crashes after processing some windows successfully
**Cause:** Terra temp files accumulate
**Check:** `du -sh /tmp/Rtmp*`
**Solution:** Already implemented auto-cleanup in code

## Recommended Processing Strategy

### Conservative Approach (Start Here)
```r
# 1. Clean restart
source("restart_processing.R")

# 2. Test single window first
source("test_single_window.R")

# 3. If test passes, run with very few parallel workers
source("AmazonHabitatLoss_prep.v2.R")  # Now uses 12 cores
```

### If Still Crashing

**Option A: Reduce window size**
Edit line 152-153 in AmazonHabitatLoss_prep.v2.R:
```r
win_nrow <- 4000  # instead of 6000
win_ncol <- 4000
```

**Option B: Reduce morphological kernel size**
Edit lines 141, 165 in calcPatchStat.sub.R:
```r
kernel_erode <- matrix(1, nrow = 3, ncol = 3)  # instead of morph_size+4
kernel_dilate <- matrix(1, nrow = 5, ncol = 5)  # instead of morph_size+8
```

**Option C: Skip morphology entirely**
Edit line 253 in AmazonHabitatLoss_prep.v2.R:
```r
apply_morph = FALSE  # Add this parameter
```

**Option D: Use even fewer cores**
Edit restart_processing.R line 20:
```r
n_cores_use <- max(1, floor(n_cores_total * 0.15))  # Use only 15% = ~10 cores
```

## Memory Calculation

With your 126GB RAM and 64 cores:
- **Current (44 cores):** ~2.9 GB per worker ❌ TOO LOW
- **Recommended (12 cores):** ~10.5 GB per worker ✓ SAFE
- **Conservative (8 cores):** ~15.8 GB per worker ✓ VERY SAFE

Each 6000x6000 window with morphological operations needs ~5-8GB RAM per worker.

## Current Status Check
```r
# See what's been completed
completed <- list.files("deforestTileOutStats/", pattern = "\\.done$")
cat("Completed:", length(completed), "of 195 windows\n")
```
