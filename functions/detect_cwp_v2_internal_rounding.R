library(terra)         # handling raster data
library(sf)            # handling vector data
library(tidyverse)     # general data handling
library(exactextractr) # extracting aggregated statistics from rasters based on polygons
library(lwgeom)        # line substring operations on sf geometries




# Detect cold water patches in a single TIR raster using zone rasterisation
# to construct a reference temperature surface.
detect_cwp_single <- function(
    ras_path,            # path to the TIR raster of the river
    line_path,           # path to the centerline of the river
    segment_length_m  = 500,   # step length for segmenting the centerline
    zone_halfwidth_m  = 60,    # halfwidth of the wide zone slabs used for rasterisation
    round_to          = 0.1,   # rounding precision for the TIR raster (°C)
    delta_C           = 1.0,   # temperature delta threshold for pixel flagging
    buffer_px         = 3,     # the number of pixels in the buffer strip to compute the reference temperature
    min_patch_area_m2 = 2,     # minimal surface area a CWP must have to be retained
    connect_diagonals = TRUE   # whether diagonally touching pixels form one polygon
) {

  # GDAL creation options: tiling improves the performance of exact_extract
  co <- "--co TILED=YES --co BLOCKXSIZE=512 --co BLOCKYSIZE=512 --co BIGTIFF=YES"

  # ── 0. READ INPUTS ────────────────────────────────────────────────────────

  # load the TIR raster from the raster path
  r_tir <- terra::rast(ras_path)
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

  # ── 1. ROUND TIR RASTER ───────────────────────────────────────────────────
  # Discretise temperature to round_to steps (e.g. 0.1 °C). This is a
  # methodological decision to suppress sensor noise before reference
  # temperature extraction and pixel-level comparison.

  # compute the rounding factor (e.g. round_to = 0.1 → round_factor = 10)
  round_factor <- 1 / round_to

  # using gdal as an external process: multiply by the rounding factor,
  # truncate, and divide back to achieve rounding. Write to r_tir_rounded.tif
  system(paste0(
    'gdal raster calc -i "A=', ras_path, '" ',
    '--calc "rint(A*', round_factor, ')/', round_factor, '" ',
    '-o r_tir_rounded.tif --overwrite --ot Float32 ',
    co
  ))

  # load the rounded raster
  r_tir_rounded <- terra::rast("r_tir_rounded.tif")

  # ── 2. BUILD SEGMENTS, ZONES & REFERENCE STRIPS ───────────────────────────
  # Cut the centerline into segments, then buffer each segment twice:
  #   - zone_slabs (wide, zone_halfwidth_m): used for burning zone IDs onto
  #     the raster grid so each pixel knows which reference temperature applies
  #   - ref_strips (narrow, ~3 pixels): used for extracting the reference
  #     temperature per segment

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

  # buffer each segment into a wide zone slab for spatial partitioning.
  # these wide slabs are used for rasterisation: each pixel in the raster
  # gets assigned the zone_id of the slab it falls in.
  # FLAT end caps prevent overlap at segment junctions.
  # MITRE join style keeps the rectangular shape.
  zone_slabs <- sf::st_buffer(segments,
                              dist        = zone_halfwidth_m,
                              endCapStyle = "FLAT",
                              joinStyle   = "MITRE",
                              mitreLimit  = 2)
  # wrap as an sf data frame and assign each zone slab a unique ID
  zone_slabs_sf <- sf::st_sf(geometry = zone_slabs) %>%
    mutate(zone_id = row_number())

  # write the zone slabs to disk as a shapefile so GDAL can read them
  sf::write_sf(zone_slabs_sf, "zone_slabs.shp")

  # buffer each segment into a narrow strip (~3 pixels wide) for reference
  # temperature extraction
  strip_halfwidth_m <- cell_m * buffer_px
  ref_strips <- sf::st_buffer(segments,
                              dist        = strip_halfwidth_m,
                              endCapStyle = "FLAT",
                              joinStyle   = "MITRE",
                              mitreLimit  = 2)
  # wrap as an sf data frame. zone_id matches between zone_slabs and
  # ref_strips because both are derived from the same segments.
  ref_strips_sf <- sf::st_sf(geometry = ref_strips) %>%
    mutate(zone_id = row_number())

  # ── 3. COMPUTE REFERENCE TEMPERATURES ────────────────────────────────────

  # extract the median temperature per reference strip from the rounded raster.
  # only cells with >= 50% coverage by the strip polygon are included.
  T_ref <- exactextractr::exact_extract(
    r_tir_rounded, ref_strips_sf,
    fun = function(values, coverage) median(values[coverage >= 0.5], na.rm = TRUE)
  )

  # create a lookup table mapping each zone_id to its reference temperature.
  # since zone_slabs and ref_strips share the same zone_id, this lookup can
  # be used to reclassify the zone raster.
  T_ref_lookup <- tibble(zone_id = ref_strips_sf$zone_id, T_ref = T_ref)

  # ── 4. RASTERISE ZONES ──────────────────────────────────────────────────
  # Each zone slab polygon gets rasterised. As a burn-in value the zone_id
  # is used. This way it is determinable from which zone each pixel came from.
  # Two rasterisations with opposing sort orders handle overlap conflicts:
  #   - ASC: rows burned in ascending order → higher zone_id wins (last write wins)
  #   - DESC: rows burned in descending order → lower zone_id wins

  # get the extent of the TIR raster and create a string for GDAL
  e          <- terra::ext(r_tir)
  ext_string <- paste(e$xmin, e$ymin, e$xmax, e$ymax, sep = ",")
  # get the resolution of the TIR raster and create a string for GDAL
  res_string <- paste(terra::res(r_tir), collapse = ",")

  # rasterise the zone slabs in ascending order: higher zone_id wins
  system(paste0(
    'gdal vector rasterize -i zone_slabs.shp -o zone_r_asc.tif ',
    '--dialect SQLITE --sql "SELECT * FROM zone_slabs ORDER BY zone_id ASC" ',
    '--extent ', ext_string, ' --resolution ', res_string,
    ' --overwrite --ot Int32 -a zone_id --optimization RASTER ',
    co
  ))

  # rasterise the zone slabs in descending order: lower zone_id wins
  system(paste0(
    'gdal vector rasterize -i zone_slabs.shp -o zone_r_desc.tif ',
    '--dialect SQLITE --sql "SELECT * FROM zone_slabs ORDER BY zone_id DESC" ',
    '--extent ', ext_string, ' --resolution ', res_string,
    ' --overwrite --ot Int32 -a zone_id --optimization RASTER ',
    co
  ))

  # ── 5. RECLASSIFY ZONES TO REFERENCE TEMPERATURES ─────────────────────────
  # Each cell currently holds a zone_id. Replace it with the corresponding
  # reference temperature from the lookup table.

  # write the lookup table to disk as a text file in GDAL reclassify format.
  # the format is: [from,to]=value; ... ; DEFAULT=NO_DATA
  writeLines(
    paste0(
      paste0("[", T_ref_lookup$zone_id, ",", T_ref_lookup$zone_id, "]=", T_ref_lookup$T_ref, collapse = "; "),
      "; DEFAULT=NO_DATA"
    ),
    "T_ref_lookup_gdal.txt"
  )

  # read the lookup string into memory for passing to GDAL
  lookup_str <- readLines("T_ref_lookup_gdal.txt")

  # reclassify the ascending zone raster to hold reference temperatures
  system(paste0(
    'gdal raster reclassify -i zone_r_asc.tif -o r_ref_asc.tif ',
    '--datatype Float64 --overwrite -m "', lookup_str, '" ',
    co
  ))

  # reclassify the descending zone raster to hold reference temperatures
  system(paste0(
    'gdal raster reclassify -i zone_r_desc.tif -o r_ref_desc.tif ',
    '--datatype Float64 --overwrite -m "', lookup_str, '" ',
    co
  ))

  # remove the zone rasters — they are large and no longer needed
  file.remove(c("zone_r_asc.tif", "zone_r_desc.tif"))

  # ── 6. MERGE REFERENCE SURFACES ──────────────────────────────────────────
  # Merge the two reference temperature rasters by taking the maximum value
  # per pixel. In overlap zones (where two zone slabs cover the same pixel),
  # the higher reference temperature is kept. This makes the flagging more
  # conservative: a pixel must be cold relative to the warmer of its two
  # possible references to be flagged as a CWP pixel.

  system(paste0(
    'gdal raster calc -i "A=r_ref_desc.tif" -i "B=r_ref_asc.tif" -o r_ref.tif ',
    '--calc "A > B ? A : B" ',
    '--overwrite --ot Float64 ',
    co
  ))

  # ── 7. PIXEL-LEVEL FLAGGING ───────────────────────────────────────────────
  # For every pair of corresponding pixels in A (rounded TIR) and B (reference):
  #   - if (B - A) >= delta_C, the pixel is colder than the reference → true
  #   - if the pixel in A is not -9999 (nodata) → true
  #   - true * true = 1 (flagged), any false → NaN (not flagged)
  # Only pixels that are both sufficiently cold AND have valid data are flagged.
  # B has a slightly larger extent than A, so the (A != -9999) check is needed
  # to guarantee that pixels outside A's valid extent are not falsely flagged.

  system(paste0(
    'gdal raster calc ',
    '-i "A=r_tir_rounded.tif" ',
    '-i "B=r_ref.tif" ',
    '--calc "((B - A) >= ', delta_C, ') * (A != -9999) ? 1 : NaN" ',
    '-o r_binary.tif --overwrite --ot Float32 ',
    co
  ))

  # remove the per-direction reference rasters (merged r_ref.tif is kept)
  file.remove(c("r_ref_desc.tif", "r_ref_asc.tif"))

  # load the raster of flagged pixels
  r_binary <- terra::rast("r_binary.tif")

  # ── 8. POLYGONISE ─────────────────────────────────────────────────────────
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

  # ── 9. AREA FILTER ────────────────────────────────────────────────────────
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

  # ── 10. TEMPERATURE STATISTICS PER PATCH ─────────────────────────────────
  # For each patch, retrieve the mean, min, max and median temperature
  # from the original (unrounded) TIR raster.

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

  # ── 11. POLYGON-LEVEL REFILTER ────────────────────────────────────────────
  # The pixel flagging in step 7 operates cell-by-cell, but the polygonisation
  # and morphological closing in step 8 can merge patches whose aggregate
  # temperature no longer satisfies the threshold. This refilter validates
  # each polygon holistically against the reference surface.
  # Because zone rasterisation creates discrete slab boundaries, patches
  # near those boundaries may have been flagged under a reference that
  # doesn't represent them well once aggregated at the polygon level.

  # load the merged reference temperature raster
  r_ref <- terra::rast("r_ref.tif")

  # compute the reference temperature for every patch by extracting the
  # median reference temperature across each patch's extent
  T_ref_per_patch <- exactextractr::exact_extract(
    r_ref, patches_with_stats,
    fun = function(values, coverage) median(values[coverage >= 0.5], na.rm = TRUE)
  )

  # add the polygon-level reference temperature, compute the temperature
  # delta, and filter out any patches where the delta is below the threshold
  patches_final <- patches_with_stats %>%
    mutate(
      T_ref  = T_ref_per_patch,
      deltaT = T_ref - T_med
    ) %>%
    filter(deltaT >= delta_C)

  # clean up all intermediate files produced by the algorithm
  file.remove(c("r_binary.tif", "r_ref.tif", "r_tir_rounded.tif", "T_ref_lookup_gdal.txt"))

  return(patches_final)
}


ras_path  <- "./data/original_data/thermal_rasters_FINAL/mean_v01emme.tif"
rounded_ras_path <- "./data/derived_data_products/r_rounded.tiff"
line_path <- "./data/original_data/Centerlines_FINAL/Emme_V01.shp"

cpw_test <- detect_cwp_single(
  ras_path = ras_path,
  line_path = line_path
)

library(tmap)
tmap_mode("view")
cpw_test %>% tm_shape() + tm_polygons()
