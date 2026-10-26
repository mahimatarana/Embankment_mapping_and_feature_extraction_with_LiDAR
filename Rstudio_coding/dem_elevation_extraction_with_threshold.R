# =============================================================================
# DEM Elevation Extraction to Top/Toe Intersection Points  (R / RStudio version)
# =============================================================================
#
# Extracts DEM elevation (Z) and projected X/Y coordinates (UTM metres)
# from the mosaic DEM raster for each of the four top/toe intersection
# point CSVs:
#
#   top_river_points.csv  -> Top of embankment - river side
#   toe_river_points.csv  -> Toe of slope       - river side
#   top_land_points.csv   -> Top of embankment  - land side
#   toe_land_points.csv   -> Toe of slope       - land side
#
# Output: same four CSVs with two new columns added:
#   elev_<name>   - DEM elevation (metres) sampled at the exact point
#   x_utm / y_utm - confirmed UTM Easting / Northing (metres)
#
# These elevation values replace the incorrect use of _y (Northing) as
# height in slope calculations.
#
# This uses the `terra` package, which is the modern replacement for the
# old raster/rgdal stack and does not have any equivalent of the Python
# `affine` cached-property bug you were hitting in rasterio.
# =============================================================================

## ---------------------------------------------------------------------------
## 1. Packages
## ---------------------------------------------------------------------------
 install.packages(c("terra", "dplyr", "readr", "ggplot2", "tidyr", "purrr"))

library(terra)    # raster I/O + sampling (replaces rasterio)
library(dplyr)    # data wrangling (replaces pandas)
library(readr)    # fast csv read/write
library(tidyr)    # pivot for plotting
library(ggplot2)  # plotting (replaces matplotlib)
library(purrr)    # iteration helpers
library(zoo)

## ---------------------------------------------------------------------------
## EDIT THESE PATHS
## ---------------------------------------------------------------------------
WORK_DIR <- "Z:/Embankment_feature_extraction/Segmentation/Final_correction_August-2026_1"
DEM_PATH <- file.path(WORK_DIR, "Mosaic_August.tif")   # your mosaicked DEM GeoTIFF

# Input CSVs - must contain columns: top_rvr_x / top_rvr_y (or equivalent)
# csv_file -> list(x_col, y_col, elev_col_name)
POINT_FILES <- list(
  "top_rvr_points.csv" = list(x_col = "top_rvr_x",   y_col = "top_rvr_y",   elev_col = "elev_top_rvr"),
  "toe_rvr_points.csv" = list(x_col = "toe_rvr_x",   y_col = "toe_rvr_y",   elev_col = "elev_toe_rvr"),
  "top_land_points.csv"  = list(x_col = "top_land_x", y_col = "top_land_y", elev_col = "elev_top_land"),
  "toe_land_points.csv"  = list(x_col = "toe_land_x", y_col = "toe_land_y", elev_col = "elev_toe_land")
)

embankment_points<-read.csv("Z:/Embankment_feature_extraction/Segmentation/Final_correction_August-2026_1/emb_segment_1m_points.csv")
names(embankment_points)

OUTPUT_DIR <- file.path(WORK_DIR, "summary_stats")
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

cat("DEM       :", DEM_PATH, "\n")
cat("Output dir:", OUTPUT_DIR, "\n")

## ---------------------------------------------------------------------------
## 2. Inspect the DEM
## ---------------------------------------------------------------------------
dem <- rast(DEM_PATH)   # terra::rast() opens the raster (metadata only, lazy-loaded)

cat("DEM file      :", DEM_PATH, "\n")
cat("CRS           :", crs(dem, describe = TRUE)$name, "\n")
cat("Width x Height:", ncol(dem), "x", nrow(dem), "px\n")
cat("Resolution    :", paste(res(dem), collapse = " x "), "m/px\n")
cat("Bounds        :", paste(as.vector(ext(dem)), collapse = ", "), "\n")
cat("NoData value  :", NAflag(dem), "\n")
cat("Dtype         :", datatype(dem), "\n")

# Full transform-style summary in one call, if you want it:
print(dem)

## ---------------------------------------------------------------------------
## 3. Elevation Extraction Function
## ---------------------------------------------------------------------------
#' Sample DEM raster at each (x, y) point and return the data frame with
#' three new columns:
#'   x_utm            - Easting  in the DEM's CRS (metres)
#'   y_utm            - Northing in the DEM's CRS (metres)
#'   <elev_col_name>  - DEM elevation (metres) at the point
#'
#' @param df            input data frame
#' @param x_col         column name for X (Easting or Longitude)
#' @param y_col         column name for Y (Northing or Latitude)
#' @param dem           a terra SpatRaster (already opened with rast())
#' @param elev_col_name name for the new elevation column
#' @param points_epsg   EPSG of the point coordinates IF they differ from the
#'                       DEM. Set to e.g. 4326 if your CSVs are in lat/lon.
#'                       Leave NULL if points are already in the same CRS as
#'                       the DEM.
extract_dem_elevation <- function(df, x_col, y_col, dem, elev_col_name,
                                   points_epsg = NULL) {

  xs <- as.numeric(df[[x_col]])
  ys <- as.numeric(df[[y_col]])

  pts <- vect(data.frame(x = xs, y = ys), geom = c("x", "y"),
              crs = if (!is.null(points_epsg)) paste0("EPSG:", points_epsg) else crs(dem))

  # Reproject points if needed
  if (!is.null(points_epsg) && !same.crs(crs(pts), crs(dem))) {
    pts <- project(pts, crs(dem))
    cat(sprintf("  Reprojected from EPSG:%s -> %s\n", points_epsg, crs(dem, describe = TRUE)$name))
  }

  coords <- crds(pts)
  df$x_utm <- coords[, "x"]
  df$y_utm <- coords[, "y"]

  # Bounds check (points outside the DEM extent)
  e <- ext(dem)
  valid <- coords[, "x"] >= e$xmin & coords[, "x"] <= e$xmax &
           coords[, "y"] >= e$ymin & coords[, "y"] <= e$ymax
  n_invalid <- sum(!valid)
  if (n_invalid > 0) {
    cat(sprintf("  WARNING: %d point(s) fall outside the DEM extent -> NA\n", n_invalid))
  }

  # Sample elevation (terra::extract handles nodata -> NA automatically)
  elevations <- rep(NA_real_, nrow(df))
  if (any(valid)) {
    sampled <- terra::extract(dem, pts[valid, ])
    elevations[valid] <- sampled[[2]]   # column 1 = ID, column 2 = raster value
  }

  df[[elev_col_name]] <- elevations
  df
}

cat("extract_dem_elevation() defined.\n")


## ---------------------------------------------------------------------------
## 4. Extract Elevation for All Four Point Files
## ---------------------------------------------------------------------------
# Set this ONLY if your CSV coordinates are in lat/lon (EPSG:4326).
# If they are already in UTM metres (same CRS as the DEM), leave as NULL.
POINTS_EPSG <- NULL   # e.g. 4326 for lat/lon, or NULL if already UTM

results <- list()

for (csv_name in names(POINT_FILES)) {
  spec     <- POINT_FILES[[csv_name]]
  csv_path <- file.path(WORK_DIR, csv_name)

  if (!file.exists(csv_path)) {
    cat(sprintf("  SKIP - file not found: %s\n", csv_path))
    next
  }

  cat("\n", strrep("=", 60), "\n", sep = "")
  cat("Processing :", csv_name, "\n")
  cat("  X col    :", spec$x_col, "\n")
  cat("  Y col    :", spec$y_col, "\n")
  cat("  Elev col :", spec$elev_col, "\n")

  df <- read_csv(csv_path, show_col_types = FALSE)
  cat("  Rows     :", nrow(df), "\n")

  df_out <- extract_dem_elevation(
    df, spec$x_col, spec$y_col,
    dem            = dem,
    elev_col_name  = spec$elev_col,
    points_epsg    = POINTS_EPSG
  )

  n_ok  <- sum(!is.na(df_out[[spec$elev_col]]))
  n_nan <- sum(is.na(df_out[[spec$elev_col]]))
  cat(sprintf("  Extracted : %d valid  |  %d NA (outside DEM or nodata)\n", n_ok, n_nan))
  cat(sprintf("  Elev range: %.3f m  ->  %.3f m\n",
              min(df_out[[spec$elev_col]], na.rm = TRUE),
              max(df_out[[spec$elev_col]], na.rm = TRUE)))

  out_path <- file.path(OUTPUT_DIR, sub("\\.csv$", "_with_elev.csv", csv_name))
  write_csv(df_out, out_path)
  cat("  Saved     :", out_path, "\n")

  results[[csv_name]] <- df_out
}

cat(sprintf("\nDone. %d file(s) processed.\n", length(results)))

## ---------------------------------------------------------------------------
## 5. Join Top/Toe Points by cross_no & Compute Slope / Width per Cross-Section
## ---------------------------------------------------------------------------
# One point each per cross_no per file (top_rvr, toe_rvr, top_land, toe_land),
# so this is a plain 1:1 join on cross_no (kept alongside Segment/seg_length
# for readability, though cross_no alone should already be the unique key).

top_rvr   <- read_csv(file.path(OUTPUT_DIR, "top_rvr_points_with_elev.csv"),  show_col_types = FALSE)
toe_rvr   <- read_csv(file.path(OUTPUT_DIR, "toe_rvr_points_with_elev.csv"),  show_col_types = FALSE)
top_land <- read_csv(file.path(OUTPUT_DIR, "top_land_points_with_elev.csv"),   show_col_types = FALSE)
toe_land <- read_csv(file.path(OUTPUT_DIR, "toe_land_points_with_elev.csv"),   show_col_types = FALSE)

JOIN_KEYS <- c("cross_no", "Segment", "seg_length")

# ── River side: top_rvr <-> toe_rvr ───────────────────────────────────────────
river_pairs <- top_rvr %>%
  inner_join(toe_rvr, by = JOIN_KEYS, suffix = c("", "_toe_rvr")) %>%
  mutate(
    horizontal_dist_river   = sqrt((toe_rvr_x - top_rvr_x)^2 + (toe_rvr_y - top_rvr_y)^2),
    vertical_distance_river = elev_top_rvr - elev_toe_rvr,
    slope_percent_river = (vertical_distance_river / horizontal_dist_river) * 100,
    slope_degree_river  = (180 / pi) * atan(vertical_distance_river / horizontal_dist_river)
  ) %>%
  filter(!is.na(slope_degree_river)) %>%
  # filter out rows based on thresholds
  filter(horizontal_dist_river   >= 0.20,
         vertical_distance_river >= 0.20)

# ── Land side: top_land <-> toe_land ────────────────────────────────────────
land_pairs <- top_land %>%
  inner_join(toe_land, by = JOIN_KEYS, suffix = c("", "_toe_land")) %>%
  mutate(
    horizontal_dist_land   = sqrt((toe_land_x - top_land_x)^2 + (toe_land_y - top_land_y)^2),
    vertical_distance_land = elev_top_land - elev_toe_land,
    slope_percent_land = (vertical_distance_land / horizontal_dist_land) * 100,
    slope_degree_land  = (180 / pi) * atan(vertical_distance_land / horizontal_dist_land)
  ) %>%
  filter(!is.na(slope_degree_land)) %>%
  # filter out rows based on thresholds
  filter(horizontal_dist_land   >= 0.20,
         vertical_distance_land >= 0.20)

# ── Width: top_rvr <-> top_land ──────────────────────────────────────────────
width_pairs <- top_rvr %>%
  inner_join(top_land, by = JOIN_KEYS, suffix = c("", "_top_land")) %>%
  mutate(width = sqrt((top_land_x - top_rvr_x)^2 + (top_land_y - top_rvr_y)^2)) %>%
  filter(!is.na(width)) %>%
  # filter out rows based on thresholds
  filter(width >= 0.4)

cat("River pairs:", nrow(river_pairs), " | Land pairs:", nrow(land_pairs),
    " | Width pairs:", nrow(width_pairs), "\n")

# Quick physical sanity check
n_rv_invalid   <- sum(river_pairs$vertical_distance_river <= 0, na.rm = TRUE)
n_land_invalid <- sum(land_pairs$vertical_distance_land   <= 0, na.rm = TRUE)
cat(sprintf("River cross-sections with toe >= top (invalid): %d / %d\n", n_rv_invalid,   nrow(river_pairs)))
cat(sprintf("Land  cross-sections with toe >= top (invalid): %d / %d\n", n_land_invalid, nrow(land_pairs)))

write_csv(river_pairs, file.path(OUTPUT_DIR, "river_pairs_t.csv"))
write_csv(land_pairs,  file.path(OUTPUT_DIR, "land_pairs_t.csv"))
write_csv(width_pairs, file.path(OUTPUT_DIR, "width_pairs_t.csv"))

## ---------------------------------------------------------------------------
## 6. Per-Segment Summary Table
## ---------------------------------------------------------------------------
p20 <- function(x) quantile(x, 0.20, na.rm = TRUE, names = FALSE)
p80 <- function(x) quantile(x, 0.80, na.rm = TRUE, names = FALSE)

river_summary <- river_pairs %>%
  group_by(Segment, seg_length) %>%
  summarise(
    n           = n(),
    river_mean  = mean(slope_degree_river, na.rm = TRUE),
    river_median= median(slope_degree_river, na.rm = TRUE),
    river_p20   = p20(slope_degree_river),
    river_p80   = p80(slope_degree_river),
    river_min   = min(slope_degree_river, na.rm = TRUE),
    river_max   = max(slope_degree_river, na.rm = TRUE),
    river_sd    = sd(slope_degree_river, na.rm = TRUE),
    .groups = "drop"
  )

land_summary <- land_pairs %>%
  group_by(Segment, seg_length) %>%
  summarise(
    land_mean  = mean(slope_degree_land, na.rm = TRUE),
    land_median= median(slope_degree_land, na.rm = TRUE),
    land_p20   = p20(slope_degree_land),
    land_p80   = p80(slope_degree_land),
    land_min   = min(slope_degree_land, na.rm = TRUE),
    land_max   = max(slope_degree_land, na.rm = TRUE),
    land_sd    = sd(slope_degree_land, na.rm = TRUE),
    .groups = "drop"
  )

width_summary <- width_pairs %>%
  group_by(Segment, seg_length) %>%
  summarise(
    width_mean  = mean(width, na.rm = TRUE),
    width_median= median(width, na.rm = TRUE),
    width_p20   = p20(width),
    width_p80   = p80(width),
    width_min   = min(width, na.rm = TRUE),
    width_max   = max(width, na.rm = TRUE),
    width_sd    = sd(width, na.rm = TRUE),
    .groups = "drop"
  )

# ── Combine slope and width statistics ─────────────────────────────────────
segment_stats <- river_summary %>%
  full_join(
    land_summary,
    by = c("Segment", "seg_length")
  ) %>%
  full_join(
    width_summary,
    by = c("Segment", "seg_length")
  ) %>%
  arrange(Segment)

# ===========================================================================
# 6A. Make sure ALL segments are retained
# ===========================================================================
#
# Use the original point data as the source of all available segments.
# This prevents segments that were removed by the slope/width filters
# from disappearing completely from the final dataset.
#
all_segments <- bind_rows(
  top_rvr %>% select(Segment, seg_length),
  toe_rvr %>% select(Segment, seg_length),
  top_land %>% select(Segment, seg_length),
  toe_land %>% select(Segment, seg_length)
) %>%
  distinct(Segment, seg_length) %>%
  arrange(Segment)

# Join calculated statistics back to ALL segments
segment_stats_full <- all_segments %>%
  left_join(
    segment_stats,
    by = c("Segment", "seg_length")
  ) %>%
  arrange(Segment)

# ===========================================================================
# 6B. Calculate DEM statistics for each embankment segment
# ===========================================================================
#
# embankment_points must contain:
#   Segment
#   Mosaic_Aug
#
# If Mosaic_Aug is the DEM elevation extracted at embankment points,
# these statistics describe the elevation distribution within each segment.
#
embankment_stats <- embankment_points %>%
  group_by(Segment) %>%
  summarise(
    DEM_mean   = mean(Mosaic_Aug, na.rm = TRUE),
    DEM_median = median(Mosaic_Aug, na.rm = TRUE),
    DEM_p20    = quantile(Mosaic_Aug, 0.20, na.rm = TRUE, names = FALSE),
    DEM_p80    = quantile(Mosaic_Aug, 0.80, na.rm = TRUE, names = FALSE),
    DEM_min    = min(Mosaic_Aug, na.rm = TRUE),
    DEM_max    = max(Mosaic_Aug, na.rm = TRUE),
    DEM_sd     = sd(Mosaic_Aug, na.rm = TRUE),
    DEM_n      = n(),
    .groups = "drop"
  )

# Add DEM statistics
segment_stats_full <- segment_stats_full %>%
  left_join(
    embankment_stats,
    by = "Segment"
  ) %>%
  arrange(Segment)

# ===========================================================================
# 6C. Fill missing numeric values using a 7-segment moving average
# ===========================================================================
#
# Important:
# Existing values are NOT changed.
# Only NA values are replaced when a 7-segment moving average can be
# calculated.
#
# The original NA pattern is saved so we can identify which values were
# estimated rather than directly calculated.
#

numeric_cols <- names(segment_stats_full)[
  sapply(segment_stats_full, is.numeric)
]

# Do not apply moving-average filling to these identifier/length fields
numeric_fill_cols <- setdiff(
  numeric_cols,
  c("seg_length")
)

# Record which rows contained missing values BEFORE filling
missing_before <- rowSums(
  is.na(segment_stats_full[, numeric_fill_cols, drop = FALSE])
)

segment_stats_modified <- segment_stats_full %>%
  mutate(
    across(
      all_of(numeric_fill_cols),
      ~ {
        x <- .x
        
        # Calculate centred 7-segment moving average
        moving_avg <- zoo::rollapply(
          x,
          width = 7,
          align = "center",
          fill = NA_real_,
          FUN = function(z) {
            if (all(is.na(z))) {
              return(NA_real_)
            }
            
            mean(z, na.rm = TRUE)
          }
        )
        
        # Fill ONLY originally missing values
        ifelse(is.na(x), moving_avg, x)
      }
    )
  )

# Mark segments that had at least one missing value before filling
segment_stats_modified <- segment_stats_modified %>%
  mutate(
    averaged = ifelse(
      missing_before > 0,
      "averaged",
      "original"
    )
  )

# ===========================================================================
# 6D. Final ordering and export
# ===========================================================================

segment_stats_modified <- segment_stats_modified %>%
  arrange(Segment)

out_path <- file.path(
  OUTPUT_DIR,
  "segment_stats.csv"
)

write_csv(segment_stats_modified, out_path)

cat(
  "\nSaved final segment statistics:",
  out_path,
  "\n"
)

print(glimpse(segment_stats_modified))
