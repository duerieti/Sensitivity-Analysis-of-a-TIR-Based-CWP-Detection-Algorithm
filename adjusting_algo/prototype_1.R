library(terra) # handling rasterdata
library(sf) # handling vectordata
library(tidyverse) # generally handling data (always usefull)
library(exactextractr) # for extracting aggregated statisitcs form rasters based on polygons
library(Rcpp) # interoperability between R and C++
library(polylabelr) # to compute the "visual center of the poygons"
library(lwgeom) # 

Rcpp::sourceCpp("idw.cpp") # load the inverse distance weighting function written in C++ and make it available in R


# set the working directory (optional and needs to be addjusted)
setwd("/home/etienne/Desktop/repos/Github_Enterprise/BSc_project/addjust_algo")

# ── HELPERS ───────────────────────────────────────────────────────────────────
# a function that produces the buffer polygons to compute the reference temperature
produce_buffers <- function(rast_path, line_path, buffer_m, step_m, slab_halfwidth_m) {

  r    <- terra::rast(rast_path) # load the raster
  line <- sf::st_read(line_path, quiet = TRUE) %>% sf::st_zm(drop = TRUE) # read in the line

  stopifnot(terra::nlyr(r) == 1) # if more than one layer stop the execution
  names(r) <- "T" # set the name of the layer to "T"
  if (terra::is.lonlat(r)) stop("Raster must be in a metric CRS.") # of the crs is not metric, then stop the execution

  cell_m <- mean(terra::res(r)) # compute the mean pixel size in meters

  # Reproject line to raster CRS if line has a defined CRS and it differs from the raster CRS
  if (!is.na(sf::st_crs(line)) &&
      !identical(terra::crs(r), sf::st_crs(line)$wkt)) {
    line <- sf::st_transform(line, terra::crs(r))
  }

  g <- sf::st_union(line)
  if (inherits(g, "sfc_GEOMETRYCOLLECTION"))
    g <- sf::st_collection_extract(g, "LINESTRING")
  if (inherits(g, "sfc_MULTILINESTRING")) {
    g <- sf::st_line_merge(g)
    if (inherits(g, "sfc_GEOMETRYCOLLECTION"))
      g <- sf::st_collection_extract(g, "LINESTRING")
  }

  L   <- as.numeric(sf::st_length(line))
  pos <- sort(unique(c(0, seq(0, 1, by = step_m / L), 1)))
  npt <- length(pos)
  if (npt < 2) stop("Too few stations for step = ", step_m)

  from <- pos[seq_len(npt - 1)]
  to   <- pos[2:npt]
  g1   <- sf::st_geometry(line)[1]

  subs <- do.call(c, lapply(seq_along(from), function(j) {
    lwgeom::st_linesubstring(g1, from[j], to[j])
  }))

  slabs_buffer <- sf::st_buffer(subs,
                                dist        = buffer_m,
                                endCapStyle = "FLAT",
                                joinStyle   = "MITRE",
                                mitreLimit  = 2)

  sf::st_sf(geometry = slabs_buffer) |>
    mutate(slab_id = row_number())
}


compute_idw_points <- function(slabs, r) {

  # a slighty addjusted median function makind shure that exact_extract behaves like terra::extract
  median_function <- function(values, coverage) {
    median(values[coverage >= 0.5], na.rm = TRUE) # median for all cells, which are more than 50% covered by the polygon.
  }

  # a function to compute the "visual center" of a polygon
  point_function <- function(sf_obj) {
    centers <- polylabelr::poi(sf_obj, precision = 0.1) # get coordinates of the "visual center" of very polygon
    purrr::map(centers, ~ sf::st_point(c(.x$x, .x$y))) %>% # for each of the coordinates, create a point
      sf::st_sfc(crs = sf::st_crs(sf_obj)) # set the crs of the new points to the crs of the polygon
  }

  # compute the reference temperature of each slab using the custom median function
  slab_medians <- exactextractr::exact_extract(r, slabs, fun = median_function)

  # add the reference temperature to each slab as a new column
  slabs %>%
    tibble::add_column(slab_medians = slab_medians) %>%
    mutate(center_point = point_function(slabs)) # and then compute the "visual center" of each slab
}
# (since the slab polygon and the point share the same row in the df
# this maps the temperature of the buffer slab to the point)

# ── MAIN FUNCTION ─────────────────────────────────────────────────────────────

# detect cold water patches in a TIR raster of a river #
detect_cwp <- function(
    ras_path, # path to the TIR raster
    line_path, # path to the centerline of the river
    min_patch_area_m2 = 2, # minimal area a cwp must have
    round_to          = 0.1,
    slab_halfwidth_m  = 60,
    delta_C           = 1.0,
    connect_diagonals = TRUE,
    n_cores           = 4 # how many cores to use for parallel execution of inverse distance weighting.
) {

  r      <- terra::rast(ras_path) # load the raster
  cell_m <- mean(terra::res(r)) # compute the mean of the raster resolution
  co     <- "--co TILED=YES --co BLOCKXSIZE=512 --co BLOCKYSIZE=512 --co BIGTIFF=YES"

  # ── 1. PRODUCE SLABS & IDW POINTS ─────────────────────────────────────────

  # sample 20 step_length parameters from a range of 500 to 3000
  step_length <- round(runif(20, min = 500, max = 3000))

  # for each of the step lengths, produce buffer_slabs
  set_of_slabs <- purrr::map(
    step_length,
    ~ produce_buffers(ras_path, line_path, buffer_m = cell_m * 3, step_m = .x, slab_halfwidth_m)
  )

  # for each set of slabs compute the temperature reference points used for inverse distance weighting
  idw_points <- purrr::map(set_of_slabs, ~ compute_idw_points(.x, r)) %>%
    bind_rows()

  # get the coordinates of the center points used for inverse distance weighting
  pts   <- sf::st_coordinates(idw_points$center_point)
  # get the temperatures assoicated with the slab and therefore also with the center point
  temps <- idw_points$slab_medians

  # ── 2. BLOCK-WISE IDW ─────────────────────────────────────────────────────

  # define a new raster which will be written to
  out <- terra::rast(r)
  bs  <- terra::blocks(r) # define blocks for writing

  terra::writeStart(out, "mean_raster.tif", overwrite = TRUE) # start the writing to disk

  for (i in seq_len(bs$n)) { # for every block
    first    <- terra::cellFromRowCol(r, bs$row[i], 1) # get the first cell in the block
    last     <- terra::cellFromRowCol(r, bs$row[i] + bs$nrows[i] - 1, terra::ncol(r)) # get the last cell in the block
    coords   <- terra::xyFromCell(r, first:last) # get the coordinates of all cells between first and last cell
    cell_vals <- as.numeric(terra::values(r, row = bs$row[i], nrows = bs$nrows[i])) # get the cell values of all cells in the block

    vals <- idw_all_cells(
      cell_x       = coords[, 1], # the x-coordinates
      cell_y       = coords[, 2], # the y-coordinates
      cell_vals    = cell_vals, # temperature values of the cell
      idw_points   = pts, # coordinates of the reference points
      temperatures = temps, # the temperatures of the reference points
      power        = 2.0, # the power used in inverse distance weighting
      threads      = n_cores # how many cores should be used for parallelism
    )

    # write the computed temperatures per cell to the output raster
    terra::writeValues(out, vals, start = bs$row[i], nrows = bs$nrows[i])
  }

  # stop the writing process
  terra::writeStop(out)

  # load the produced raster
  mean_raster <- terra::rast("mean_raster.tif")



  # ── 3. ROUND RASTER ───────────────────────────────────────────────────────

  # compute the rounding factor
  rfactor <- if (!is.null(round_to) && round_to > 0) 1 / round_to else stop("round_to must be > 0")

  # using gdal as an individual process spawned through system():
  # round every pixel in the raster using the rounding factor and write to a new file
  # r_rounded.tiff
  system(paste0(
    'gdal raster calc -i "A=', ras_path, '" ',
    '--calc "A*', rfactor, '/', rfactor, '" ',
    '-o r_rounded.tiff --overwrite --ot Float32 ',
    co
  ))
  
  # load the rounded raster
  r_rounded <- terra::rast("r_rounded.tiff")

  # if the difference between reference raster and the rounded raster is larger than
  # delta_C, then flagg the cell using 1. If the difference is smaller or
  # the cell in A is a NoData cell (-9999), then set NaN
  # write the result to a raster called binary_out.tif
  system(paste0(
    'gdal raster calc ',
    '-i "B=mean_raster.tif" ',
    '-i "A=r_rounded.tiff" ',
    '--calc "((B - A) >= ', delta_C, ') * (A != -9999) ? 1 : NaN" ',
    '-o binary_out.tif --overwrite --ot Float32 ',
    co
  ))

  # load the raster
  binary <- terra::rast("binary_out.tif")

  # ── 5. POLYGONIZE ─────────────────────────────────────────────────────────

  # a small number equal to 10% of the cell size
  eps <- if (connect_diagonals) cell_m * 0.1 else 0

# Convert binary raster to vector polygons, dissolving adjacent cells into one polygon. eight = TRUE means diagonal neighbours are
# considered connected (8-connectivity). So diagonal cells are also considered as adjacent
patches_sf <- terra::as.polygons(binary, dissolve = TRUE, eight = TRUE) %>%
  sf::st_as_sf() %>%
  # Split any MULTIPOLYGON into individual POLYGONs so each patch is a separate row
  sf::st_cast("POLYGON") %>%
  # Expand slightly to close tiny gaps between adjacent polygons from raster-to-vector conversion
  sf::st_buffer(eps) %>%
  # Merge all overlapping/touching polygons into one geometry
  sf::st_union() %>%
  # Shrink back to original size, keeping merged boundaries (morphological closing)
  sf::st_buffer(-eps) %>%
  # Split merged result back into individual polygons
  sf::st_cast("POLYGON") %>%
  sf::st_as_sf()
  
  # ── 6. AREA FILTER ────────────────────────────────────────────────────────

  patches_large <- patches_sf %>%
    mutate(
      area_m2 = as.numeric(sf::st_area(.)), # get the area of each CWP polygon
      ID      = row_number() # also assign a row number
    ) %>%
    filter(area_m2 >= min_patch_area_m2) 
    # if a patch has a smaller area then beeing set, then
    # patch is removed


  # ── 7. TEMPERATURE STATISTICS PER PATCH ───────────────────────────────────
  
  # for each patch retrieve the mean, min, max and median temperature
  stats_per_poly <- exactextractr::exact_extract(r, patches_large, c("mean", "min", "max", "median")) %>%
    as_tibble() %>%
    mutate(ID = patches_large$ID) # add the patch ID from the filtered patches

  patches_large_w_stats <- patches_large %>%
    inner_join(stats_per_poly, by = "ID") %>% # join stats onto the patches using the ID
    rename(T_min = min, T_max = max, T_med = median, T_mean = mean) # some renaming for nicer output

  # ── 8. SLAB MEAN FILTER ───────────────────────────────────────────────────

  # load the reference temperature raster
  Tmean_raster <- terra::rast("mean_raster.tif")

  # what was this for again? I do not remeber
  slap_means_per_poly <- exactextractr::exact_extract(
    Tmean_raster, patches_large_w_stats,
    fun = function(values, coverage) median(values[coverage >= 0.5], na.rm = TRUE)
  )
  
  # filter them out
  patches_final <- patches_large_w_stats %>%
    mutate(
      Tmd_slb = slap_means_per_poly,
      deltaT  = Tmd_slb - T_med
    ) %>%
    filter(deltaT >= delta_C)

  # clean up the intermediate files produced by the algorithm
  file.remove(c("binary_out.tif", "mean_raster.tif", "r_rounded.tiff"))

  return(patches_final)
}


# ── RUN ───────────────────────────────────────────────────────────────────────

# define the path to the raster
ras_path  <- "../data/derived_data_products/mean_v01emme.tif"
# define the path to the line
line_path <- "../data/original_data/Centerlines_FINAL/Emme_V01.shp"

# detect the CWP in the raster
cwps <- detect_cwp(ras_path, line_path, n_cores = 12)


tmap_mode("view")
cwps %>% tm_shape() + tm_polygons()



