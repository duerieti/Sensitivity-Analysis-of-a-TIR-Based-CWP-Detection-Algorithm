
library(terra)
library(sf)
library(exactextractr)

library(tidyverse)
library(Rcpp)

setwd("/home/etienne/Desktop/repos/Github_Enterprise/BSc_project/addjust_algo")


produce_slabs <- function(
  # This function produces a set of 
  #
  rast_path, 
  line_path,
  buffer_m,
  step_m,
  slab_halfwidth_m


) {
  # ── 0. READ INPUTS ────────────────────────────────────────────────────────

r <- terra::rast(ras_path)
line <- sf::st_read(line_path, quiet = TRUE) |> sf::st_zm(TRUE, "ZM")

stopifnot(terra::nlyr(r) == 1)
names(r) <- "T"
if (terra::is.lonlat(r)) stop("Raster must be in a metric CRS.")

rres   <- terra::res(r)
cell_m <- mean(rres)

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

  # ── 1. BUILD SLABS ────────────────────────────────────────────────────────

L <- as.numeric(sf::st_length(line))
S <- step_m

pos <- sort(unique(c(0, seq(0, 1, by = S / L), 1)))
npt <- length(pos)

if (npt < 2) stop("Too few stations
  for step = ", S)

from <- pos[1:(npt - 1)]
to   <- pos[2:npt]
k    <- length(from)

g1 <- sf::st_geometry(line)[1]

subs_list <- vector("list", k)
for (j in seq_len(k)) {
  subs_list[[j]] <- lwgeom::st_linesubstring(g1, from[j], to[j])
}

subs <- do.call(c, subs_list)

slabs <- sf::st_buffer(subs,
                       dist        = slab_halfwidth_m,
                       endCapStyle = "FLAT",
                       joinStyle   = "MITRE",
                       mitreLimit  = 2)
slabs_sf <- sf::st_sf(geometry = slabs)

slabs_terra <- terra::vect(slabs_sf)

buffer_m <- cell_m * 3

slabs_buffer <- sf::st_buffer(subs,
                              dist        = buffer_m,
                              endCapStyle = "FLAT",
                              joinStyle   = "MITRE",
                              mitreLimit  = 2)

slabs_buffer_sf <- sf::st_sf(geometry = slabs_buffer) %>% mutate(slap_id = row_number())


return(

slabs_buffer_sf
  
)
}


compute_idw_points <- function(slabs, r) {

  median_function <- function(values, coverage) {
    median(values[coverage >= 0.5], na.rm = TRUE)
  }

  slab_medians <- exact_extract(r, slabs, fun = median_function)

  slabs %>%
    tibble::add_column(slab_medians = slab_medians) %>%  # preserves sf
    mutate(
      center_point = st_centroid(geometry)                # sfc_POINT column
    )
}

cppFunction('
  #include <cmath>
  #include <omp.h>

  NumericVector idw_all_cells(NumericVector cell_x,
                              NumericVector cell_y,
                              NumericMatrix idw_points,
                              NumericVector temperatures,
                              double power = 2.0,
                              int threads = 14) {
      int ncells = cell_x.size();
      int npts   = idw_points.nrow();
      NumericVector out(ncells);

      #pragma omp parallel for num_threads(threads) schedule(static)
      for (int i = 0; i < ncells; i++) {
          double x = cell_x[i];
          double y = cell_y[i];
          double wsum = 0.0, vwsum = 0.0;

          for (int j = 0; j < npts; j++) {
              double dx = idw_points(j, 0) - x;
              double dy = idw_points(j, 1) - y;
              double d  = std::sqrt(dx*dx + dy*dy);
              if (d == 0.0) { vwsum = temperatures[j]; wsum = 1.0; break; }
              double w  = 1.0 / std::pow(d, power);
              wsum  += w;
              vwsum += w * temperatures[j];
          }
          out[i] = vwsum / wsum;
      }
      return out;
  }
', plugins = "openmp")



detect_cwp <- function(
    ras_path,
    line_path,
    min_patch_area_m2 = 2,
    round_to = 0.1,
    slab_halfwidth_m = 60,
    delta_C = 1.0,
    connect_diagonals = TRUE,
    n_cores = 4
) 

{
  
  r <- terra::rast(ras_path)

  rres <- terra::res(r)
  cell_m <- mean(rres)

  co <- "--co TILED=YES --co BLOCKXSIZE=512 --co BLOCKYSIZE=512 --co BIGTIFF=YES"


  buffer_m <- cell_m * 3

  slab_lengths <- runif(10,min = 500, max = 3000) %>% round()
  set_of_slabs <- map(slab_lengths, ~ produce_slabs(ras_path, line_path, buffer_m,., slab_halfwidth_m))

  idw_points <-map(set_of_slabs, ~compute_idw_points(., r)) %>% bind_rows()



  # --- Extract points and temperatures from sf ---
  pts   <- st_coordinates(idw_points$center_point)  # n x 2 matrix
  temps <- idw_points$slab_medians



  # --- Block processing ---
  bs      <- blocks(r)
  out <- rast(r)
  out_con <- writeStart(out, "idw_out.tif", overwrite = TRUE)


  for (i in seq_len(bs$n)) {
    
    first  <- cellFromRowCol(r, bs$row[i], 1)
    last   <- cellFromRowCol(r, bs$row[i] + bs$nrows[i] - 1, ncol(r))
    coords <- xyFromCell(r, first:last)


    vals <- idw_all_cells(
      cell_x       = coords[, 1],
      cell_y       = coords[, 2],
      idw_points   = pts,
      temperatures = temps,
      power        = 2.0,
      threads      = n_cores
  )
  writeValues(out, vals, start = bs$row[i], nrows = bs$nrows[i])
  }

writeStop(out)

mean_raster <- rast("idw_out.tif")

file.remove("idw_out.tif")
  


  
mean_raster_cropped <- crop(mean_raster, r)        # crop to bounding box
r_masked  <- mask(mean_raster_cropped, r) # then mask to exact shape


writeRaster(r_masked, "mean_raster.tif")

rfactor <- if (!is.null(round_to) && round_to > 0) 1 / round_to else NA_real_

system(paste0(
  'gdal raster calc -i "A=', ras_path, '" ',
  '--calc "rint( A*', rfactor, ')/', rfactor, '" ',
  '-o r_rounded.tiff --overwrite --ot Float32 ',
  co
))
  
r_rounded <- terra::rast("r_rounded.tiff")


system(paste0(
  'gdal raster calc ',
  '-i "B=mean_raster.tif" ',
  '-i "A=../../data/r_rounded.tiff" ',
  '--calc "((B - A) >= ', delta_C, ') * (A != -9999) ? 1 : NaN" ',
  '-o binary_out.tif --overwrite --ot Float32 ',
  co
))

binary <- terra::rast("binary_out.tif")


patches_v <- terra::as.polygons(binary, dissolve = TRUE, eight = connect_diagonals)


# ── 7. POLYGONIZE ─────────────────────────────────────────────────────────

eps <- if (connect_diagonals) cell_m * 0.1 else 0

patches_sf <- patches_v %>%
  sf::st_as_sf() %>%
  sf::st_cast("POLYGON") %>%
  sf::st_buffer(eps) %>%
  sf::st_union() %>%
  sf::st_buffer(-eps) %>%
  sf::st_cast("POLYGON") %>%
  sf::st_as_sf()


  # ── 8. AREA FILTER ────────────────────────────────────────────────────────
patches_large <- patches_sf %>%
  mutate(area_m2 = st_area(.) %>% as.numeric()) %>%
  filter(area_m2 >= min_patch_area_m2) %>%
  mutate(ID = row_number())

  # ── 9. TEMPERATURE STATISTICS PER PATCH ──────────────────────────────────

print("second exact_extract")
stats_matrix <- exact_extract(r, patches_large, c("mean", "min", "max", "median"))

stats_per_poly <- stats_matrix %>%
  as_tibble() %>%
  mutate(ID = patches_large$ID)

patches_large_w_stats <- patches_large %>%
  inner_join(
    stats_per_poly,
    by = join_by(ID==ID)
  ) %>%
  rename(
    T_min=min,
    T_max=max,
    T_med=median,
    T_mean=mean,
    geometry=x

  )

Tmean_raster <- terra::rast("mean_raster.tif")


slap_means_per_poly <- exact_extract(Tmean_raster, patches_large_w_stats, fun = function(values, coverage) {
  median(values[coverage >= 0.5], na.rm = TRUE)
})

slap_means_df <- tibble(Tmd_slb = slap_means_per_poly, ID = patches_large_w_stats$ID)

patches_large_refiltered <- patches_large_w_stats %>%
  inner_join(
    slap_means_df,
    by = join_by(ID == ID)
  ) %>%
  mutate(
    deltaT = Tmd_slb - T_med
  ) %>%
  filter(deltaT >= delta_C) 


file.remove(c("binary_out.tif","mean_raster.tif"))

return(patches_large_refiltered)

}

ras_path <- file.path("../data/derived_data_products/mean_v01emme.tif")
line_path <- file.path("../data/original_data/Centerlines_FINAL/Emme_V01.shp")

detect_cwp(ras_path, line_path,n_cores = 12)

