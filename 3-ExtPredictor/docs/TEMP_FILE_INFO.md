# Terra Temp File Management

## Where Temp Files Are Stored

By default, terra stores temporary raster files in `/tmp/RtmpXXXXXX/` directories, where each R session gets its own unique directory.

**Problem with parallel processing:**
- Each worker creates its own temp directory
- With many workers and large rasters, this can use significant disk space
- Old temp directories accumulate over time (340 found on your system!)
- Currently using ~4GB in old R temp directories

## What I've Fixed

### 1. Shared Temp Directory
- All workers now share a parent temp directory
- Each worker gets its own subdirectory (by process ID) to avoid conflicts
- Located at: `/tmp/RtmpXXXXXX/terra_workers/worker_XXXXX/`

### 2. Automatic Cleanup
- Each worker cleans up its own temp files after completing a window
- Shared temp directory is deleted when parallel processing finishes
- Added `tmpFiles(remove = TRUE)` calls

### 3. Reduced Memory Pressure
- Window size reduced: 15000 → 6000 (6.25x less memory)
- Cores reduced: 50% → 25% (4x more memory per worker)
- Memory allocation increased: 0.5 → 0.7 multiplier

## Manual Cleanup

If you want to clean up old R temp directories from previous runs:

```r
source("cleanup_temp.R")
```

Or manually:
```bash
# See what's using space
du -sh /tmp/Rtmp*

# Remove old directories (be careful!)
find /tmp -maxdepth 1 -name "Rtmp*" -mtime +1 -exec rm -rf {} \;
```

## Monitoring During Processing

Watch temp file usage in real-time:
```bash
watch -n 5 'du -sh /tmp/Rtmp* 2>/dev/null | tail -5'
```

Or check just your current session:
```bash
du -sh /tmp/Rtmp*/terra_workers/
```

## If You Still Run Out of Space

1. Further reduce window size (try 4000 or even 3000)
2. Reduce number of workers even more (0.2 instead of 0.25)
3. Set a custom temp directory on a larger disk:
   ```r
   terraOptions(tempdir = "/path/to/larger/disk/terra_temp")
   ```
4. Reduce morph_size in calcPatchStat.sub.R (currently using large kernels)
