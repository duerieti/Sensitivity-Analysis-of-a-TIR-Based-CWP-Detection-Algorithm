library(tidyverse)
library(tmap)
library(terra)

setwd("compute_morris")

# Compute a fuzzy CWP detection probability map for one river section by
# rasterising all CWP realisations (one per Morris parameter tuple), stacking
# them into a virtual raster, and averaging pixel-wise across all bands.
# The output is a single raster where each pixel holds the fraction of parameter
# tuples for which that pixel was detected as a CWP pixel:
#   - value close to 1 → robust, parameter-insensitive detection
#   - value close to 0 → marginal, parameter-sensitive detection
compute_fuzzy_cwp <- function(
    results_dir,  # path to the folder containing one subfolder per parameter tuple
    tir_path,     # path to the TIR raster (defines extent and resolution)
    binary_dir,   # path to write the intermediate binary rasters
    vrt_path,     # path to write the virtual raster stack
    output_path   # path to write the fuzzy detection probability raster
) {

  # ── 0. LOCATE CWP REALISATIONS ──────────────────────────────────────────────
  # Each subfolder contains a shapefile of CWP polygons detected under one
  # parameter tuple. Build the full paths to all shapefiles.

  folder_names <- list.files(results_dir) %>%
    # exclude the joblist.txt file
    .[str_detect(., ".txt", negate = TRUE)]

  read_paths <- folder_names %>%
    file.path(results_dir, ., "final_polys.shp")

  # ── 1. RASTERISE EACH CWP REALISATION ────────────────────────────────────────
  # Burn each CWP shapefile into a binary raster at the extent and resolution
  # of the TIR raster: CWP pixels → 1, background → 0.

  r          <- terra::rast(tir_path)
  e          <- terra::ext(r)
  ext_string <- paste(e$xmin, e$ymin, e$xmax, e$ymax, sep = ",")
  res_string <- paste(terra::res(r), collapse = ",")

  # create the output directory for the binary rasters if it does not exist yet
  dir.create(binary_dir, showWarnings = FALSE)

  for (i in seq_along(read_paths)) {
    out_path <- file.path(binary_dir, paste0("rast_", i, ".tif"))
    system(paste0(
      'gdal vector rasterize -i ', read_paths[i], ' -o ', out_path, ' ',
      '--extent ', ext_string, ' --resolution ', res_string,
      ' --burn 1 --init 0 ',   # CWP pixels → 1, background → 0
      ' --overwrite --ot Byte --optimization RASTER'
    ))
  }

  # ── 2. STACK ALL BINARY RASTERS ──────────────────────────────────────────────
  # Combine all binary rasters into a virtual raster stack (.vrt). A VRT is a
  # lightweight pointer file — no data is copied, so this is instant regardless
  # of how many rasters are stacked. Each band corresponds to one parameter tuple.

  raster_paths        <- list.files(binary_dir, pattern = "\\.tif$", full.names = TRUE)
  raster_paths_string <- paste(raster_paths, collapse = " ")

  system(paste0(
    'gdal raster stack ', raster_paths_string, ' ', vrt_path, ' --overwrite'
  ))

  # ── 3. COMPUTE PIXEL-WISE DETECTION PROBABILITY ───────────────────────────────
  # Average across all bands using --flatten. Each pixel in the output holds
  # the fraction of parameter tuples for which that pixel was flagged as CWP.

  system(paste0(
    'gdal raster calc -i "A=', vrt_path, '" --flatten --calc "avg(A)" -o ', output_path
  ))

  # ── 4. INSPECT RESULT ─────────────────────────────────────────────────────────
  # Load and interactively visualise the result to verify before further analysis.

  fuzzy_cwp <- terra::rast(output_path)
  tmap_mode("view")
  print(tm_shape(fuzzy_cwp) + tm_raster())
}


# ── RUN FOR EACH RIVER SECTION ────────────────────────────────────────────────

compute_fuzzy_cwp(
  results_dir = "../morris_screening_obemme/results_obemme",
  tir_path    = "../data/original_data/thermal_rasters_FINAL/mean_v2obemme.tif",
  binary_dir  = "./binary_rasters_obemme",
  vrt_path    = "stacked_raster_obemme.vrt",
  output_path = "fuzzy_cwp_obemme.tif"
)

compute_fuzzy_cwp(
  results_dir = "../morris_screening_emme_v2/results_emmev2",
  tir_path    = "../data/original_data/thermal_rasters_FINAL/mean_v2emme.tif",
  binary_dir  = "./binary_rasters_upper_emme",
  vrt_path    = "stacked_raster_upper_emme.vrt",
  output_path = "fuzzy_cwp_upper_emme.tif"
)

compute_fuzzy_cwp(
  results_dir = "../morris_screening_emme_v1/results_emmev1",
  tir_path    = "../data/original_data/thermal_rasters_FINAL/mean_v01emme.tif",
  binary_dir  = "./binary_rasters_emme_v1",
  vrt_path    = "stacked_raster_emme_v1.vrt",
  output_path = "fuzzy_cwp_emmev1.tif"
)