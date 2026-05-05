library(terra)
library(sf)
library(dplyr)
library(tidyverse)
library(tictoc)
library(tmap)
library(exactextractr)
library(Rcpp)



setwd("/home/etienne/Desktop/repos/Github_Enterprise/BSc_project/addjust_algo")



ras_path <- file.path("../data/r_rounded.tiff")
line_path <- file.path("../data/Centerlines_FINAL/Emme_V01.shp")




ras_path
line_path
step_m            = 500
buffer_px         = 3
delta_C           = 1.0
min_patch_area_m2 = 2
round_to          = 0.1
slab_halfwidth_m  = 60
connect_diagonals = TRUE
union_chunk_size  = 50000L


co <- "--co TILED=YES --co BLOCKXSIZE=512 --co BLOCKYSIZE=512 --co BIGTIFF=YES"


produce_slabs <- function(
  rast_path, 
  line_path,
  step_m, 
  buffer_px, 
  slab_halfwidth_m = 60

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

if (npt < 2) stop("Too few stations for step = ", S)

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

buffer_m <- cell_m * buffer_px

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


buffer_lengths <- runif(10, 0.1,1) * 2000



set_of_slabs <- map(buffer_lengths, ~ produce_slabs(ras_path, line_path, ., buffer_px = 6))

set_of_slabs %>% length()


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


idw_points <-map(set_of_slabs, ~compute_idw_points(., r)) %>% bind_rows()


tm_shape(idw_points %>% bind_rows()) +
  tm_dots()




cppFunction('
  double compute_weighted_mean(NumericVector point,
                               NumericMatrix idw_points,
                               NumericVector temperatures,
                               double power = 2.0) {
      double x = point(0);
      double y = point(1);
      int n    = idw_points.nrow();

      double wsum = 0.0, vwsum = 0.0;

      for (int i = 0; i < n; i++) {
          double dx = idw_points(i, 0) - x;
          double dy = idw_points(i, 1) - y;
          double d  = std::sqrt(dx*dx + dy*dy);
          if (d == 0.0) return temperatures[i];
          double w  = 1.0 / std::pow(d, power);
          wsum  += w;
          vwsum += w * temperatures[i];
      }
      return vwsum / wsum;
  }
')


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

# --- Extract points and temperatures from sf ---
pts   <- st_coordinates(idw_points$center_point)  # n x 2 matrix
temps <- idw_points$slab_medians



# --- Block processing ---
terraOptions(memmax = 8)  # smaller blocks
bs      <- blocks(r)
out_con <- writeStart(r, "idw_out.tif", overwrite = TRUE)

cat("Total blocks:", bs$n, "\n")
cat("Cells per block:", bs$nrows[1] * ncol(r), "\n")

t_start <- Sys.time()

for (i in seq_len(bs$n)) {
  print(i)
  first  <- cellFromRowCol(r, bs$row[i], 1)
  last   <- cellFromRowCol(r, bs$row[i] + bs$nrows[i] - 1, ncol(r))
  coords <- xyFromCell(r, first:last)


  vals <- idw_all_cells(
    cell_x       = coords[, 1],
    cell_y       = coords[, 2],
    idw_points   = pts,
    temperatures = temps,
    power        = 2.0,
    threads      = 12
  )
  writeValues(r, vals, start = bs$row[i], nrows = bs$nrows[i])
}

writeStop(r)

for (i in seq_len(bs$n)) {
  print(i)
  first  <- cellFromRowCol(r, bs$row[i], 1)
  last   <- cellFromRowCol(r, bs$row[i] + bs$nrows[i] - 1, ncol(r))
  coords <- xyFromCell(r, first:last)

  vals <- idw_all_cells(
    cell_x       = coords[, 1],
    cell_y       = coords[, 2],
    idw_points   = pts,
    temperatures = temps,
    power        = 2.0,
    threads      = 12
  )
  writeValues(r, vals, start = bs$row[i], nrows = bs$nrows[i])
}

class(out_con)

writeStop(out_con)
vals %>% class()

# sample 100k cells from your raster
sample_cells  <- sample(ncell(r), 1000000)
sample_coords <- xyFromCell(r, sample_cells)

system.time({
  idw_all_cells(
    cell_x       = sample_coords[, 1],
    cell_y       = sample_coords[, 2],
    idw_points   = pts,
    temperatures = temps,
    power        = 2.0,
    threads      = 12
  )
})

mean_raster <- rast("idw_out.tif")




  
mean_raster_cropped <- crop(mean_raster, r)        # crop to bounding box
r_masked  <- mask(mean_raster_cropped, r) # then mask to exact shape


writeRaster(r_masked, "mean_raster.tif")



system(paste0(
  'gdal raster calc ',
  '-i "B=mean_raster.tif" ',
  '-i "A=../../data/r_rounded.tiff" ',
  '--calc "((B - A) >= ', delta_C, ') * (A != -9999) ? 1 : NaN" ',
  '-o binary_out.tif --overwrite --ot Float32 ',
  co
))


binary <- terra::rast("binary_out.tif")


patches_v <- terra::as.polygons(binary, dissolve = TRUE, eight = TRUE)


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
  filter(area_m2 >= 2) %>%
  mutate(ID = row_number())

  # ── 9. TEMPERATURE STATISTICS PER PATCH ──────────────────────────────────

tmap_mode("view")
tm_shape(patches_large ) + tm_polygons()

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
    T_mean=mean

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


sf::write_sf(patches_large_refiltered, "final_polys.shp")



