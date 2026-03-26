library(terra)
library(sf)
library(dplyr)
library(lwgeom)
library(tidyverse)
library(tictoc)
library(tmap)
library(exactextractr)


workdir <- Sys.getenv("WORKDIR")

print("dirs")
list.dirs(workdir) %>% print()

print("files")
list.files(workdir, recursive = TRUE) %>% print()



path_raster <- file.path(workdir, "data/thermal_rasters_FINAL/mean_v01emme.tif")
path_line <- file.path(workdir, "data/Centerlines_FINAL/Emme_V01.shp")


detect_cwp_single <- function(
    ras_path,
    line_path,
    step_m            = 500,
    buffer_px         = 3,
    delta_C           = 1.0,
    min_patch_area_m2 = 2,
    round_to          = 0.1,
    slab_halfwidth_m  = 60,
    connect_diagonals = TRUE,
    union_chunk_size  = 50000L
) {

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
  
rfactor <- if (!is.null(round_to) && round_to > 0) 1 / round_to else NA_real_

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
slabs_buffer_sf <- sf::st_sf(geometry = slabs_buffer)

slaps_buffer_vect <- slabs_buffer_sf %>%
  mutate(slap_id = row_number()) %>%
  terra::vect()

slabs <- sf::st_buffer(subs,
                       dist        = slab_halfwidth_m,
                       endCapStyle = "FLAT",
                       joinStyle   = "MITRE",
                       mitreLimit  = 2)
slabs_sf <- sf::st_sf(geometry = slabs)

slaps_sf_vect <- slabs_sf %>%
  mutate(slap_id = row_number()) %>%
  terra::vect()

zone_r_big <- terra::rasterize(slaps_sf_vect, r, field = "slap_id")

r_rounded <- round(r * rfactor) / rfactor

zone_r <- terra::rasterize(slaps_buffer_vect, r_rounded, field = "slap_id")

mean_raster <- terra::zonal(r_rounded, zone_r, fun = "median", na.rm = TRUE)

Tmean <- terra::classify(zone_r_big, mean_raster)

flagged_pixels <- r_rounded - Tmean
binary <- flagged_pixels <= (-1 * delta_C)

binary[binary == 0] <- NA
patches_v <- terra::as.polygons(binary, dissolve = TRUE, eight = FALSE)

eps <- if (connect_diagonals) cell_m * 0.1 else 0

patches_sf <- patches_v %>%
  sf::st_as_sf() %>%
  sf::st_cast("POLYGON") %>%
  sf::st_buffer(eps) %>%
  sf::st_union() %>%
  sf::st_buffer(-eps)

patches_large <- patches_sf %>%
  filter(as.numeric(st_area(.)) >= 2) %>%
  mutate(ID = row_number())

stats_per_poly <- terra::extract(r_rounded, patches_large) %>%
  group_by(ID) %>%
  summarise(
    mean_temp   = mean(T),
    min_temp    = min(T),
    max_temp    = max(T),
    median_temp = median(T)
  )

patches_large_with_stas <- patches_large %>%
  left_join(stats_per_poly, by = "ID") %>%
  select(!T)

slab_means <- terra::extract(Tmean, vect(patches_large_with_stas), fun = "mean") %>%
  rename(slab_mean_T = slap_id)

final_patches <- patches_large_with_stas %>%
  left_join(slab_means, by = "ID") %>%
  mutate(deltaT = slab_mean_T - mean_temp) %>%
  filter(deltaT >= delta_C)

return(final_patches)

}

detect_cwp_single(path_raster, path_line)

unlink(terra_tmp, recursive=TRUE)
