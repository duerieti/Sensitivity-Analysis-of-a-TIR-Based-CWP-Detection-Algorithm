library(terra)
library(sf)
library(dplyr)
library(tidyverse)
library(tictoc)
library(tmap)
library(exactextractr)

tic("Total Time: ")
ras_path <- file.path("../data/thermal_rasters_FINAL/mean_v01emme.tif")
line_path <- file.path("../data/Centerlines_FINAL/Emme_V01.shp")


    step_m            = 500
    buffer_px         = 3
    delta_C           = 1.0
    min_patch_area_m2 = 2
    round_to          = 0.1
    slab_halfwidth_m  = 60
    connect_diagonals = TRUE
    union_chunk_size  = 50000L


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

rfactor <- if (!is.null(round_to) && round_to > 0) 1 / round_to else NA_real_


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
slaps_buffer_sf <- sf::st_sf(geometry = slabs_buffer) %>% mutate(slap_id = row_number())

slaps_buffer_vect <- slaps_buffer_sf %>%
  mutate(slap_id = row_number()) %>%
  terra::vect()

slabs <- sf::st_buffer(subs,
                       dist        = slab_halfwidth_m,
                       endCapStyle = "FLAT",
                       joinStyle   = "MITRE",
                       mitreLimit  = 2)
slabs_sf <- sf::st_sf(geometry = slabs) %>% mutate(slap_id = row_number())
sf::write_sf(slabs_sf, "slabs_sf.shp")



  # ── 2. COMPUTE REFERENCE TEMPERATURES ────────────────────────────────────

e <- ext(r)

ext_string <- paste(e$xmin, e$ymin, e$xmax, e$ymax, sep = ",")
resolution <- res(r)
res_string <- paste(resolution, collapse = ",")


system(
  paste0(
    'gdal vector rasterize -i slabs_sf.shp -o zone_r_big.tiff --extent ', ext_string, ' --resolution ', res_string, ' --overwrite --ot Int8 -a slap_id --optimization RASTER'
  )
)


tic("tic direct gdal approach:")
system(
  paste0('gdal raster calc -i "A=../data/thermal_rasters_FINAL/mean_v01emme.tif" --calc "A*', rfactor , '/', rfactor,'"', ' -o r_rounded.tiff --overwrite --ot Float32')
)
toc()
# 118.197


r_rounded <- terra::rast("r_rounded.tiff")

tic()
slap_means <- exactextractr::exact_extract(r_rounded, slaps_buffer_sf, fun = "median")
toc()

temp_look_up_df <- tibble(slap_id = slaps_buffer_sf$slap_id, slap_means = slap_means)



  # ── 5. BURN MEDIANS INTO ZONES ────────────────────────────────────────────

writeLines(
  paste0(
    paste0("[", temp_look_up_df$slap_id, ",", temp_look_up_df$slap_id, "]=", temp_look_up_df$slap_means, collapse = "; "),
    "; DEFAULT=NO_DATA"
  ),
  "median_lookup_gdal.txt"
)

lookup_str <- readLines("median_lookup_gdal.txt")

tic()
system(paste0(
  'gdal raster reclassify -i zone_r_big.tiff  -o Tmean_raster.tiff --datatype Float64 --overwrite  -m "',
  lookup_str,
  '"'
))
toc()
# 91 seconds
# so when I also use gdal in the step for the zone_r_big creation, then this is a speedup of roughly 1/3 here



tic()
system(

paste0('gdal raster calc -i "A=r_rounded.tiff" -i "B=Tmean_raster.tiff" --calc "((A - B) >= ', delta_C, ') ? 1 : NaN" -o binary_out.tif --overwrite')
)
toc()
# 405 seconds. So roughy a doubling in performance

binary <- terra::rast("binary_out.tif")

tic()
patches_v <- terra::as.polygons(binary, dissolve = TRUE, eight = FALSE)
toc()
## Doing gdal here does not make sense i beliefe. Because I was not able to track down the command that lets you compute

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
  filter(as.numeric(st_area(.)) >= 2) %>%
  mutate(ID = row_number())


  # ── 9. TEMPERATURE STATISTICS PER PATCH ──────────────────────────────────

tic()
stats_matrix <- exact_extract(r_rounded, patches_large, c("mean", "min", "max", "median"))
toc()

stats_per_poly <- stats_matrix %>%
  as_tibble() %>%
  mutate(
      ID = patches_large$ID
  )
# 5.894 seconds -> so writing the file in chunks is a huuge speedup. So this deffinetely has to be noted
patches_large_with_stas <- patches_large %>%
  left_join(stats_per_poly, by = "ID")


slap_median_vector <- exact_extract(Tmean, patches_large_with_stas, fun = "median")

slap_median_df <- slap_median_vector %>%
  as_tibble() %>%
  mutate(
      ID = patches_large$ID
  ) %>% rename(slap_median_T = value)

final_patches <- patches_large_with_stas %>%
  left_join(slap_median_df, by = "ID") %>%
  mutate(deltaT = median- slap_median_T) %>%
  filter(deltaT >= delta_C)

toc()
