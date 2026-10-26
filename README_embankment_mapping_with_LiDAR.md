# Embankment Mapping with LiDAR

[![Language](https://img.shields.io/badge/Language-R-276DC3?logo=r)](https://www.r-project.org/)
[![Platform](https://img.shields.io/badge/Platform-DJI%20Terra%20%7C%20CloudCompare%20%7C%20ArcGIS-lightgrey)](.)
[![Data](https://img.shields.io/badge/Data-LiDAR%20%7C%20RTK%20GNSS%20%7C%20SRTM-green)](.)
[![Segments](https://img.shields.io/badge/Output-4%2C224%20segments-blue)](.)
[![License](https://img.shields.io/badge/License-Academic%20Use-yellow)](.)

---

## Overview

River embankments are critical flood-protection infrastructure. Yet systematically measuring their condition — height, width, slope — across tens of kilometres is rarely feasible with traditional survey methods.

This project builds a complete, reproducible pipeline to do exactly that. Starting from raw drone imagery and LiDAR point clouds captured in the field, it produces a spatially referenced database of **4,224 embankment segments** at 100-metre intervals, each with statistics on DEM elevation, crest width, river-side slope, and land-side slope. The pipeline runs from drone deployment to analysis-ready CSV and shapefile, and is documented end-to-end in an accompanying field and processing manual.

---

## Study Area

**Southwestern Bangladesh** — a low-lying coastal delta where embankments (locally called *polders*) are the primary defence against tidal flooding and cyclone surge. Survey grids were flown using a DJI drone with an L2 LiDAR sensor and RTK base station. The coordinate system used throughout is **WGS84 / UTM Zone 45N**.

---

## Pipeline at a Glance

```
Field (drone + RTK)
        │
        ▼
DJI Terra ──► Orthomosaic + LiDAR Point Cloud
        │
        ▼
CloudCompare ──► Ground Classification + DEM Export (50 cm)
        │
        ▼
ArcMap ──► Projection, Clipping, Mosaic → Mosaic_August.tif
        │
        ▼
R (DEM_height_correction_manual.R)
  ├─ Process 1: SRTM-based height correction
  │    └─ Resample modelled DEM to each tile (bilinear), trim top 25% of
  │       elevation values (removes embankment features), apply median dZ offset
  ├─ Process 2: Polygon footprint creation per corrected DEM tile
  ├─ Process 3: RTK ground-truth correction
  │    └─ Spatial join of RTK points to tile polygons → per-tile median dZ
  │       correction using centimetre-accurate ground control (RL values)
  └─ Process 4: Overlap-based correction for tiles without RTK coverage
       └─ Detect overlapping RTK-corrected tiles, sample 50 random overlap
          points, apply median dZ correction → Overlap_RTK_*.tif
        │
        ▼
ArcGIS (ET GeoWizards) ──► 100 m segment lines + perpendicular cross-sections
  └─ Four intersection point CSVs exported:
       top_rvr_points.csv · toe_rvr_points.csv
       top_land_points.csv · toe_land_points.csv
        │
        ▼
R (dem_elevation_extraction_with_threshold.R)
  ├─ Step 1: Inspect mosaic DEM (CRS, resolution, bounds, NoData)
  ├─ Step 2: Extract DEM elevation + confirmed UTM X/Y to all four
  │          intersection point CSVs → *_with_elev.csv
  ├─ Step 3: Compute per-cross-section geometry
  │    ├─ River side  → horizontal dist, vertical dist, slope (% and °)
  │    │                filter: horiz ≥ 0.20 m AND vert ≥ 0.20 m
  │    ├─ Land side   → same geometry and thresholds
  │    └─ Crest width → Euclidean distance top_rvr ↔ top_land; filter ≥ 0.4 m
  ├─ Step 4: Per-segment summary (mean, median, P20, P80, min, max, SD)
  │          for slope (river + land) and width
  ├─ Step 5: DEM elevation statistics per segment from embankment
  │          centreline points (Mosaic_Aug column)
  └─ Step 6: Gap-fill missing segments with 7-point centred moving average
             (original values preserved; filled rows flagged as "averaged")
        │
        ▼
Output: segment_stats.csv — 4,224 georeferenced segments
```

---

## Repository Structure

```
Embankment_mapping_with_LiDAR/
├── code/
│   ├── DEM_height_correction_manual.R              # Processes 1–4: elevation correction
│   └── dem_elevation_extraction_with_threshold.R   # Elevation extraction, geometry,
│                                                   # segment stats, gap-fill
├── data/
│   ├── segment_stats.csv                           # Final output — 4,224 segments
│   └── shapefile/
│       ├── Embankment_100m_segment_database.shp
│       ├── Embankment_100m_segment_database.dbf
│       ├── Embankment_100m_segment_database.shx
│       ├── Embankment_100m_segment_database.prj
│       ├── Embankment_100m_segment_database.sbn
│       ├── Embankment_100m_segment_database.sbx
│       └── Embankment_100m_segment_database.cpg
└── docs/
    └── Drone_survey_and_data_processing_Manual_Mahima.docx   # Full field + processing manual
```

---

## Output Dataset

The final database (`segment_stats.csv`) contains **4,224 georeferenced segments**. Each row represents one 100-metre embankment segment with the following fields:

| Field group | Columns | Description |
|---|---|---|
| Identity | `Segment`, `seg_length` | Segment ID and cumulative length |
| DEM elevation | `DEM_mean`, `DEM_median`, `DEM_p20`, `DEM_p80`, `DEM_min`, `DEM_max`, `DEM_sd` | Elevation statistics from the corrected mosaic DEM (metres) |
| Crest width | `width_mean`, `width_median`, `width_p20`, `width_p80`, `width_min`, `width_max`, `width_sd` | Euclidean distance between land-top and river-top boundary (metres) |
| River-side slope | `river_mean`, `river_median`, `river_p20`, `river_p80`, `river_min`, `river_max`, `river_sd` | Slope angle toward the river (degrees) |
| Land-side slope | `land_mean`, `land_median`, `land_p20`, `land_p80`, `land_min`, `land_max`, `land_sd` | Slope angle toward the land (degrees) |
| Data quality | `averaged`, `n` | "averaged" / "original" flag; count of valid cross-sections per segment |

**Summary statistics (DEM elevation, metres):**

| Metric | Value |
|---|---|
| Median | 3.51 m |
| 25th percentile | 3.18 m |
| 75th percentile | 3.90 m |
| Max | 7.35 m |

---

## R Scripts

### `DEM_height_correction_manual.R`

Corrects elevation offsets across multiple LiDAR-derived DEM tiles using four sequential processes. Depends on `terra` and `sf`.

**Process 1 — SRTM correction**
Lists all `.tif` files in the correction folder (excluding `modeled_dem_rf.tif`). For each tile, resamples the SRTM-modelled DEM to match the tile grid (bilinear interpolation), trims the top 25% of elevation values to exclude embankment features, and applies the median dZ difference as a flat offset. Corrected files are saved with a `corrected_` prefix.

**Process 2 — Polygon footprint creation**
Converts each SRTM-corrected raster tile to a dissolved polygon footprint (`.shp`) using `terra::as.polygons()`. These footprints are used in the next step to spatially assign RTK ground control points to their respective tiles.

**Process 3 — RTK correction**
Reads RTK ground control points from `RTK_points.shp` and spatially joins them to tile polygon footprints to identify which tile each RTK point falls within. For each tile that has RTK coverage, extracts DEM elevation at the RTK point locations, computes per-point dZ (RTK RL minus DEM value), and applies the median dZ as the correction offset. Outputs saved with an `RTK_` prefix.

**Process 4 — Overlap-based correction**
For tiles that received no RTK correction, identifies spatially overlapping RTK-corrected (`RTK_`) tiles by comparing raster extents and valid-pixel masks. Randomly samples 50 overlap points, extracts elevations from both the uncorrected target and the corrected reference, and applies the median dZ to the target. Results saved as `Overlap_RTK_*.tif` with a correction summary CSV.

---

### `dem_elevation_extraction_with_threshold.R`

A single end-to-end script that takes the four ArcGIS intersection point CSVs and the mosaicked DEM as inputs and produces the final segment statistics CSV. Depends on `terra`, `sf`, `dplyr`, `readr`, `tidyr`, `ggplot2`, `purrr`, and `zoo`.

**Step 1 — DEM inspection**
Opens `Mosaic_August.tif` with `terra::rast()` and prints CRS, dimensions, resolution, spatial extent, NoData value, and data type.

**Step 2 — Elevation extraction**
Defines `extract_dem_elevation()`, a function that converts XY coordinates to a `terra::SpatVector`, optionally reprojects from a user-specified EPSG to the DEM's CRS, performs a bounds check, and samples DEM elevation at each point using `terra::extract()`. Runs this function across all four intersection point files (`top_rvr`, `toe_rvr`, `top_land`, `toe_land`), appending `elev_*` and confirmed `x_utm` / `y_utm` columns. Outputs saved as `*_with_elev.csv`. This step replaces the earlier incorrect use of the Northing (Y) coordinate as elevation height.

**Step 3 — Per-cross-section geometry**
Joins the four enriched files on `cross_no`, `Segment`, and `seg_length`. Computes:
- **River side**: horizontal distance, vertical distance, slope in % and degrees
- **Land side**: same
- **Crest width**: Euclidean distance between river-top and land-top coordinates

Filters applied before summarising:
- River and land sides: `horizontal_dist ≥ 0.20 m` AND `vertical_dist ≥ 0.20 m`
- Width: `width ≥ 0.4 m`

A physical sanity check flags cross-sections where the toe elevation is at or above the top elevation.

**Step 4 — Segment summary (slope + width)**
Groups valid cross-sections by `Segment` and `seg_length`. Computes mean, median, P20, P80, min, max, and SD for river slope, land slope, and crest width. Uses `all_segments` (the union of all four input files) as the backbone so that segments removed by the quality filters are retained in the output with NA values rather than dropped.

**Step 5 — DEM elevation statistics**
Reads `emb_segment_1m_points.csv` (embankment centreline points with a `Mosaic_Aug` elevation column) and computes DEM mean, median, P20, P80, min, max, SD, and count per segment. Joins these to the slope/width summary.

**Step 6 — Gap-filling and export**
Applies a 7-point centred moving average (`zoo::rollapply`) to all numeric statistic columns across segments, filling only rows that were NA — original values are never overwritten. Segments that received at least one imputed value are flagged `"averaged"`; all others are marked `"original"`. Final table is sorted by `Segment` and written to `segment_stats.csv`.

---

## How to Run

### Prerequisites

```r
install.packages(c("terra", "sf", "dplyr", "readr", "tidyr", "ggplot2", "purrr", "zoo"))
```

### Steps

1. **Elevation correction** — update the folder paths in `DEM_height_correction_manual.R` to point to your local DEM and RTK directories (currently set to network paths used during the survey programme). Run Processes 1 through 4 in order.

2. **ArcGIS cross-sections** — generate 100 m segment lines, perpendicular cross-sections, and the four intersection point CSVs in ArcGIS using ET GeoWizards (see the `docs/` manual, RStudio Processing section).

3. **Mosaic the corrected DEMs** — merge the corrected tiles in ArcMap to produce `Mosaic_August.tif`.

4. **Elevation extraction and segment statistics** — update `WORK_DIR` and `DEM_PATH` in `dem_elevation_extraction_with_threshold.R`, then run the script end-to-end. Output: `summary_stats/segment_stats.csv`.

5. **Spatial join** — join `segment_stats.csv` back to the 100 m segment shapefile in ArcGIS to produce the final spatial database.

---

## Tools Used

| Tool | Purpose |
|---|---|
| DJI Terra | Point cloud reconstruction, orthomosaic generation |
| CloudCompare | Ground classification, DEM export |
| ArcMap (ArcGIS Desktop) | Projection, mosaicking, cross-section analysis |
| ET GeoWizards | Automated perpendicular cross-section generation |
| RStudio (R) | Elevation correction, DEM extraction, segment statistics |
| RTK GNSS | Centimetre-level ground control |

---

## Citation

If you use the data, scripts, or methodology from this repository, please cite:

> Tarana, M. (2025). *Embankment Mapping with LiDAR: A Drone-Based Survey and Processing Pipeline for River Embankment Condition Assessment*. GitHub repository: https://github.com/mahimatarana/Embankment_mapping_with_LiDAR

---

## License

This repository is shared for academic and research purposes. Please contact the author before reproducing or adapting the data, code, or documentation for other uses.
