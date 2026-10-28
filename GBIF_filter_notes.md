# GBIF Occurrence Cleaning — Filter Notes
# AmAni pipeline | 00_cleanGBIF.R
# Last updated: 2026-03-27

---

## Filters currently implemented (00_cleanGBIF.R)

| # | Filter | Notes |
|---|--------|-------|
| 1 | `taxonRank == "SPECIES"` | Removes genus-, family-, and order-level records |
| 2 | Non-empty species name | Drops blank / NA species fields |
| 3 | Undetermined IDs | `sp.`, `spp.`, `cf.`, `aff.`, `nr.`, bare trailing `sp` |
| 4 | Missing coordinates | Drops rows with NA lat or lon |
| 5 | Coordinate uncertainty > 10 km | NA uncertainty records are kept (many museum records lack this field) |
| 6 | Absence records | Drops `occurrenceStatus == ABSENT` |
| 7 | High-impact GBIF issue flags | `ZERO_COORDINATE`, `COORDINATE_INVALID`, `COUNTRY_COORDINATE_MISMATCH`, `COORDINATE_OUT_OF_RANGE` |
| 8 | Whitespace standardisation | Trim + collapse internal whitespace in species name |
| 9 | Cultivated / introduced records | `establishmentMeans` in MANAGED, CULTIVATED, INTRODUCED, INVASIVE — **added 2026-03-27** |
| 10 | Implausible collection year | Year < 1800 or > current year — **added 2026-03-27** |
| 11 | Exact spatial duplicates | Same species + identical lat + lon — **added 2026-03-27** |
| 12 | Extended uncertain ID regex | `agg.`, `sensu lato`, `s.l.`, `sl`, `group`, `complex`, `?`, `×` — **added 2026-03-27** |

---

## Filters requested / under consideration

### From user wishlist

| Filter | Status | Notes |
|--------|--------|-------|
| Undetermined species (sp., cf., aff., etc.) | ✓ Implemented (filters 3 + 12) | |
| Non-specific OTUs | Partially covered | `taxonRank == SPECIES` handles most; numeric-epithet OTUs (e.g. "Aus 12") could still slip through |
| Invasive species (explicit list) | Deferred | Would require GRIIS (Global Register of Introduced and Invasive Species) list per taxon group, downloadable from gbif.org/dataset/b351a324-... Distinct from the establishmentMeans filter (filter 9), which catches records already flagged in GBIF but misses many invasives |
| Misspelled / unresolvable species names | Deferred | Best approach: `rgbif::name_backbone_checklist()` or offline GBIF backbone fuzzy match. Slow; should be cached. Flags names matching no known taxon — user can then decide whether to drop or correct |
| Occurrences with imprecise coordinates | Partially covered | Filter 5 (10 km uncertainty threshold) + filter 7 (GBIF flags). Gap: many records lack `coordinateUncertaintyInMeters` entirely; a round-coordinate proxy filter (≤ 2 decimal places ≈ 1 km) would catch more of these |
| Uncertain IDs (?, sensu lato, group, etc.) | ✓ Implemented (filter 12) | Trailing numeric morphotype suffixes (e.g. "Aus bus 1") are intentionally retained |

### Additional suggestions (from internal review)

| Filter | Priority | Notes |
|--------|----------|-------|
| Round-coordinate proxy filter | High | Flag/remove records where both lat AND lon are rounded to ≤ 2 decimal places when `coordinateUncertaintyInMeters` is NA. Approx 1 km grid at equator — still coarser than the 22 km study grid, so arguably acceptable, but worth flagging |
| Country / province centroid blacklist | Medium | The `CoordinateCleaner` R package (`cc_cen()`) maintains a database of political centroids that GBIF's own flags miss. Removes records placed exactly on capital cities, country centres, etc. |
| Study polygon spatial filter | Medium | Clip to Amazon shapefile before writing cleaned file. Records outside the polygon land on NA grid cells downstream anyway, but removing them early reduces file size and processing time |
| GBIF taxon match failure flags | Medium | `issue` column values `TAXON_MATCH_NONE` and `TAXON_MATCH_HIGHERRANK` indicate records where GBIF couldn't confidently assign taxonomy. The rank filter catches some but not all of these |
| iNaturalist-sourced records for invertebrates | Low | iNaturalist records (`datasetKey == 50c9509d-22c7-4a22-a47d-8c48425ef4a7`) are crowd-sourced and often lack expert validation for invertebrate groups. An optional flag to exclude or downweight them may be appropriate for some taxon groups |
| Duplicate removal by gbifID | Low | Deduplication by `gbifID` before any other step to handle the case of duplicate rows within the GBIF download itself (rare but possible with large multi-dataset downloads) |

---

## Notes on morphotype-style suffixes

Records like `"Aus bus 1"` or `"Aus bus 2a"` (trailing numeric/alphanumeric after a valid binomial)
are intentionally **kept** — they represent field morphotypes and are a legitimate part of some
curated datasets even if they don't correspond to formally described species.

The uncertain-ID regex (filter 12) is designed to match only recognised taxonomic uncertainty
markers (`agg.`, `group`, `complex`, `?`, `×`, etc.) and does not match trailing numbers.

---

## Parameter reference (clean_gbif)

| Parameter | Default | Description |
|-----------|---------|-------------|
| `input_file` | — | Path to raw GBIF CSV or ZIP |
| `output_file` | NULL | Derived from input if NULL |
| `max_uncertainty_m` | 10000 | Coordinate uncertainty ceiling (metres); NA records kept |
| `min_year` | 1800 | Earliest plausible collection year |
| `verbose` | TRUE | Print per-filter record counts |
