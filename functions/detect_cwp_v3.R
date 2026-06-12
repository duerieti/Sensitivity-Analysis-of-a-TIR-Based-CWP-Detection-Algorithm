library(terra)         # handling raster data
library(sf)            # handling vector data
library(tidyverse)     # general data handling
library(exactextractr) # extracting aggregated statistics from rasters based on polygons
library(polylabelr)    # computing the pole of inaccessibility ("visual center") of polygons
library(lwgeom)        # line substring operations on sf geometries

# load the inverse distance weighting function written in C++ and make it available in R
Rcpp::sourceCpp("idw.cpp")

# set the working directory (needs to be adjusted per machine)
setwd("/home/etienne/Desktop/repos/Github_Enterprise/BSc_project/addjust_algo")

# ── HELPERS ───────────────────────────────────────────────────────────────────

# Cut the river centerline into segments of length segment_length_m and buffer
# each segment into a narrow strip (~3 pixels wide). The reference temperature
# is extracted from these strips.
produce_ref_strips <- function(rast_path, line_path, segment_length_m) {

  # load the TIR raster from the raster path
  r_tir <- terra::rast(rast_path)
  # read in the river centerline and drop any Z/M dimensions
  line  <- sf::st_read(line_path, quiet = TRUE) %>% sf::st_zm(drop = TRUE)

  # the raster must have exactly one layer
  stopifnot(terra::nlyr(r_tir) == 1)
  # set the name of the layer to "T" for temperature
  names(r_tir) <- "T"
  # if the CRS is not metric, stop — all distance calculations assume metres
  if (terra::is.lonlat(r_tir)) stop("Raster must be in a metric CRS.")

  # compute the mean pixel size in metres
  cell_m <- mean(terra::res(r_tir))

  # reproject line to raster CRS if the line has a defined CRS and it
  # differs from the raster CRS
  if (!is.na(sf::st_crs(line)) &&
      !identical(terra::crs(r_tir), sf::st_crs(line)$wkt)) {
    line <- sf::st_transform(line, terra::crs(r_tir))
  }

  # dissolve all line features into a single geometry
  g <- sf::st_union(line)
  # st_union can return a GEOMETRYCOLLECTION if the input contains mixed types;
  # extract only the LINESTRING parts in that case
  if (inherits(g, "sfc_GEOMETRYCOLLECTION"))
    g <- sf::st_collection_extract(g, "LINESTRING")
  # if the result is a MULTILINESTRING, merge it into a single LINESTRING
  # where possible
  if (inherits(g, "sfc_MULTILINESTRING")) {
    g <- sf::st_line_merge(g)
    # st_line_merge can also return a GEOMETRYCOLLECTION, so extract again
    if (inherits(g, "sfc_GEOMETRYCOLLECTION"))
      g <- sf::st_collection_extract(g, "LINESTRING")
  }

  # compute the total length of the centerline in metres
  L <- as.numeric(sf::st_length(line))
  # create normalised positions [0, 1] along the line at intervals of
  # segment_length_m. Endpoints 0 and 1 are always included.
  pos <- sort(unique(c(0, seq(0, 1, by = segment_length_m / L), 1)))
  # if there are fewer than 2 positions, segmentation is impossible
  if (length(pos) < 2) stop("Too few segments for segment_length_m = ", segment_length_m)

  # define the start and end of each segment as normalised positions
  from <- pos[seq_len(length(pos) - 1)]
  to   <- pos[2:length(pos)]
  # extract the first geometry from the line (needed for st_linesubstring)
  g1   <- sf::st_geometry(line)[1]

  # cut the centerline into line substrings, one per segment
  segments <- do.call(c, lapply(seq_along(from), function(j) {
    lwgeom::st_linesubstring(g1, from[j], to[j])
  }))

  # compute the halfwidth of the reference strip in metres (~3 pixels)
  strip_halfwidth_m <- cell_m * 3
  # buffer each segment into a narrow rectangular strip for reference
  # temperature extraction. FLAT end caps prevent overlap at segment junctions.
  # MITRE join style keeps the rectangular shape.
  strips <- sf::st_buffer(segments,
                          dist        = strip_halfwidth_m,
                          endCapStyle = "FLAT",
                          joinStyle   = "MITRE",
                          mitreLimit  = 2)

  # wrap as an sf data frame and assign each strip a unique ID
  sf::st_sf(geometry = strips) %>%
    mutate(strip_id = row_number())
}


# For each reference strip, compute the median temperature (T_ref) and the
# pole of inaccessibility (visual center). Together these form one IDW
# sample point: a location and a reference temperature.
compute_idw_points <- function(ref_strips, r_tir) {

  # repair any invalid geometries that may have arisen from buffering
  ref_strips <- sf::st_make_valid(ref_strips)

  # a slightly adjusted median function making sure that exact_extract
  # behaves like terra::extract: only cells with >= 50% coverage are included
  median_function <- function(values, coverage) {
    median(values[coverage >= 0.5], na.rm = TRUE)
  }

  # a function to compute the pole of inaccessibility ("visual center") of
  # each polygon. The pole of inaccessibility is the point inside the polygon
  # that is farthest from any edge. precision = 1 means 1 metre accuracy.
  point_function <- function(sf_obj) {
    # get coordinates of the visual center of every polygon
    centers <- polylabelr::poi(sf_obj, precision = 1)
    # for each set of coordinates, create an sf point geometry
    purrr::map(centers, ~ sf::st_point(c(.x$x, .x$y))) %>%
      # set the CRS of the new points to match the input polygons
      sf::st_sfc(crs = sf::st_crs(sf_obj))
  }

  # compute the reference temperature of each strip using the custom
  # median function
  T_ref <- exactextractr::exact_extract(r_tir, ref_strips, fun = median_function)

  # add the reference temperature as a new column to each strip
  ref_strips %>%
    tibble::add_column(T_ref = T_ref) %>%
    # compute the visual center of each strip and add as a geometry column.
    # since the strip polygon and the center point share the same row,
    # this maps the reference temperature of the strip to the point.
    mutate(center_point = point_function(ref_strips))
}


# ── MAIN FUNCTION ─────────────────────────────────────────────────────────────

# Detect cold water patches in a TIR raster of a river using IDW interpolation
# to construct a continuous reference temperature surface.
detect_cwp <- function(
    ras_path,          # path to the TIR raster
    line_path,         # path to the centerline of the river
    min_patch_area_m2 = 2,   # minimal area a CWP must have to be retained
    round_to          = 0.1, # rounding precision for the TIR raster (°C)
    delta_C           = 1.0, # temperature delta threshold for pixel flagging
    connect_diagonals = TRUE, # whether diagonally touching pixels form one polygon
    n_cores           = 4    # number of cores for OpenMP parallelism in IDW
) {

  # load the TIR raster
  r_tir  <- terra::rast(ras_path)
  # compute the mean pixel size in metres
  cell_m <- mean(terra::res(r_tir))
  # GDAL creation options: tiling improves performance of exact_extract
  co     <- "--co TILED=YES --co BLOCKXSIZE=512 --co BLOCKYSIZE=512 --co BIGTIFF=YES"

  # ── 1. PRODUCE REFERENCE STRIPS & IDW SAMPLE POINTS ───────────────────────
  # Sample 20 segment lengths from a range of 500 to 3000 metres to produce
  # overlapping strip sets. Using multiple lengths ensures the IDW surface is
  # not sensitive to any single segmentation choice.

  # draw 20 random segment lengths
  segment_lengths <- round(runif(20, min = 500, max = 3000))

  # for each segment length, cut the centerline and produce reference strips
  strip_sets <- purrr::map(
    segment_lengths,
    ~ produce_ref_strips(ras_path, line_path, segment_length_m = .x)
  )

  # for each set of strips, compute the reference temperature and the
  # center point used for IDW interpolation, then combine all into one table
  idw_samples <- purrr::map(strip_sets, ~ compute_idw_points(.x, r_tir)) %>%
    bind_rows()

  # extract the x/y coordinates of the center points as a matrix
  ref_coords <- sf::st_coordinates(idw_samples$center_point)
  # extract the reference temperatures associated with each center point
  ref_temps  <- idw_samples$T_ref

  # ── 2. BLOCK-WISE IDW INTERPOLATION ───────────────────────────────────────
  # Interpolate a continuous reference temperature surface from the sample
  # points using inverse distance weighting. The raster is processed in
  # blocks to avoid loading the entire grid into memory at once. Cells
  # outside the river (NaN in r_tir) are skipped by the C++ function.

  # create an empty output raster with the same dimensions as the TIR raster
  r_ref <- terra::rast(r_tir)
  # define blocks for block-wise writing
  bs    <- terra::blocks(r_tir)
  # start writing to disk
  terra::writeStart(r_ref, "r_ref.tif", overwrite = TRUE)

  # for every block of rows in the raster
  for (i in seq_len(bs$n)) {
    # get the cell index of the first cell in the block
    first     <- terra::cellFromRowCol(r_tir, bs$row[i], 1)
    # get the cell index of the last cell in the block
    last      <- terra::cellFromRowCol(r_tir, bs$row[i] + bs$nrows[i] - 1, terra::ncol(r_tir))
    # get the x/y coordinates of all cells in the block
    coords    <- terra::xyFromCell(r_tir, first:last)
    # get the temperature values of all cells in the block (used for NaN masking)
    cell_vals <- as.numeric(terra::values(r_tir, row = bs$row[i], nrows = bs$nrows[i]))

    # run the IDW interpolation for all cells in the block
    vals <- idw_all_cells(
      cell_x       = coords[, 1],    # x-coordinates of the cells
      cell_y       = coords[, 2],    # y-coordinates of the cells
      cell_vals    = cell_vals,       # temperature values (only used for NaN check)
      idw_points   = ref_coords,      # coordinates of the IDW reference points
      temperatures = ref_temps,       # temperatures of the IDW reference points
      power        = 2.0,            # decay power of the inverse distance weighting
      threads      = n_cores         # number of OpenMP threads
    )

    # write the interpolated values for this block to the output raster
    terra::writeValues(r_ref, vals, start = bs$row[i], nrows = bs$nrows[i])
  }

  # finalise the writing process
  terra::writeStop(r_ref)

  # ── 3. ROUND TIR RASTER ───────────────────────────────────────────────────
  # Discretise temperature to round_to steps (e.g. 0.1 °C). This is a
  # methodological decision to suppress sensor noise before pixel-level
  # comparison.

  # compute the rounding factor (e.g. round_to = 0.1 → round_factor = 10)
  round_factor <- 1 / round_to

  # using gdal as an external process: multiply by the rounding factor,
  # truncate, and divide back to achieve rounding. Write to r_tir_rounded.tif
  system(paste0(
    'gdal raster calc -i "A=', ras_path, '" ',
    '--calc "A*', round_factor, '/', round_factor, '" ',
    '-o r_tir_rounded.tif --overwrite --ot Float32 ',
    co
  ))

  # ── 4. PIXEL-LEVEL FLAGGING ───────────────────────────────────────────────
  # For every pair of corresponding pixels in A (rounded TIR) and B (reference):
  #   - if (B - A) >= delta_C, the pixel is colder than the reference → flag with 1
  #   - if the pixel in A is -9999 (nodata), exclude it
  #   - true * true = 1, any false = NaN
  # Only pixels that are both sufficiently cold AND have valid data are flagged.

  system(paste0(
    'gdal raster calc ',
    '-i "B=r_ref.tif" ',
    '-i "A=r_tir_rounded.tif" ',
    '--calc "((B - A) >= ', delta_C, ') * (A != -9999) ? 1 : NaN" ',
    '-o r_binary.tif --overwrite --ot Float32 ',
    co
  ))

  # load the binary raster of flagged pixels
  r_binary <- terra::rast("r_binary.tif")

  # ── 5. POLYGONISE ─────────────────────────────────────────────────────────
  # Convert the binary raster of flagged pixels to vector polygons.
  # dissolve = TRUE: neighbouring flagged pixels are fused into one polygon
  # (otherwise each pixel would become its own tiny square polygon).
  # eight = connect_diagonals: if TRUE, diagonally touching pixels are also
  # considered neighbours (8-connectivity instead of 4-connectivity).

  # a small tolerance value (10% of cell size) for morphological closing
  closing_tolerance <- if (connect_diagonals) cell_m * 0.1 else 0

  # convert flagged pixels to polygons, then apply morphological closing:
  # expand → union → shrink. This merges patches that are separated by
  # sub-pixel gaps from the raster-to-vector conversion.
  patches_sf <- terra::as.polygons(r_binary, dissolve = TRUE, eight = connect_diagonals) %>%
    sf::st_as_sf() %>%
    # split any MULTIPOLYGON into individual POLYGONs so each patch is a separate row
    sf::st_cast("POLYGON") %>%
    # expand slightly to close tiny gaps between adjacent polygons
    sf::st_buffer(closing_tolerance) %>%
    # merge all overlapping/touching polygons into one geometry
    sf::st_union() %>%
    # shrink back to original size, keeping merged boundaries
    sf::st_buffer(-closing_tolerance) %>%
    # split the merged result back into individual polygons
    sf::st_cast("POLYGON") %>%
    sf::st_as_sf()

  # ── 6. AREA FILTER ────────────────────────────────────────────────────────
  # All cold water patches smaller than min_patch_area_m2 are filtered out.

  patches <- patches_sf %>%
    mutate(
      # compute the area of each CWP polygon in square metres
      area_m2 = as.numeric(sf::st_area(patches_sf)),
      # assign each patch a unique ID
      ID      = row_number()
    ) %>%
    # remove patches below the minimum area threshold
    filter(area_m2 >= min_patch_area_m2)

  # ── 7. TEMPERATURE STATISTICS PER PATCH ───────────────────────────────────
  # For each patch, retrieve the mean, min, max and median temperature
  # from the original TIR raster.

  # extract temperature statistics for each patch polygon
  stats <- exactextractr::exact_extract(r_tir, patches, c("mean", "min", "max", "median")) %>%
    as_tibble() %>%
    # add the patch ID so we can join back onto the patches
    mutate(ID = patches$ID)

  # join the temperature statistics onto the patch polygons using the ID
  patches_with_stats <- patches %>%
    inner_join(stats, by = "ID") %>%
    # rename columns for cleaner output
    rename(T_min = min, T_max = max, T_med = median, T_mean = mean)

  # ── 8. POLYGON-LEVEL REFILTER ─────────────────────────────────────────────
  # The pixel flagging in step 4 operates cell-by-cell, but the polygonisation
  # and morphological closing in step 5 can merge patches whose aggregate
  # temperature no longer satisfies the threshold. This refilter validates
  # each polygon holistically against the reference surface.
  # Because the IDW surface is smooth (no discrete slab boundary jumps),
  # this filter removes fewer patches than in the zone-rasterisation approach
  # but remains a valid safeguard against false positives introduced by
  # the polygonisation step.

  # load the reference temperature raster from disk
  r_ref_disk <- terra::rast("r_ref.tif")

  # compute the reference temperature for every patch by extracting the
  # median reference temperature across each patch's extent
  T_ref_per_patch <- exactextractr::exact_extract(
    r_ref_disk, patches_with_stats,
    fun = function(values, coverage) median(values[coverage >= 0.5], na.rm = TRUE)
  )

  # add the polygon-level reference temperature and compute the temperature
  # delta, then filter out any patches where the delta is below the threshold
  patches_final <- patches_with_stats %>%
    mutate(
      T_ref  = T_ref_per_patch,
      deltaT = T_ref - T_med
    ) %>%
    filter(deltaT >= delta_C)

  # clean up all intermediate files produced by the algorithm
  file.remove(c("r_binary.tif", "r_ref.tif", "r_tir_rounded.tif"))

  return(patches_final)
}


# ── RUN ───────────────────────────────────────────────────────────────────────

ras_path  <- "../data/derived_data_products/mean_v01emme.tif"
line_path <- "../data/original_data/Centerlines_FINAL/Emme_V01.shp"

detect_cwp(ras_path, line_path, n_cores = 12)
