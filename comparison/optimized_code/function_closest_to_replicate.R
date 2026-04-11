library(terra)
library(sf)
library(dplyr)
library(tidyverse)
library(tictoc)
library(tmap)
library(exactextractr)

setwd("/home/etienne/Desktop/repos/Github_Enterprise/BSc_project/comparison/optimized_code")


tic("Total Time Chunked workflow: ")
ras_path <- file.path("../../data/thermal_rasters_FINAL/mean_v01emme.tif")
line_path <- file.path("../../data/Centerlines_FINAL/Emme_V01.shp")

    step_m            = 500
    buffer_px         = 3
    delta_C           = 1.0
    min_patch_area_m2 = 2
    round_to          = 0.1
    slab_halfwidth_m  = 60
    connect_diagonals = TRUE
    union_chunk_size  = 50000L

# tiling creation options — no compression, just chunked layout
co <- "--co TILED=YES --co BLOCKXSIZE=512 --co BLOCKYSIZE=512 --co BIGTIFF=YES"

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

sf::st_layers("slabs_sf.shp")

e <- ext(r)
ext_string <- paste(e$xmin, e$ymin, e$xmax, e$ymax, sep = ",")
resolution <- res(r)
res_string <- paste(resolution, collapse = ",")

system(paste0(
  'gdal vector rasterize -i slabs_sf.shp -o zone_r_big_asc.tiff ',
  '--dialect SQLITE --sql "SELECT * FROM slabs_sf ORDER BY slap_id ASC" ',
  '--extent ', ext_string, ' --resolution ', res_string,
  ' --overwrite --ot Int32 -a slap_id --optimization RASTER ',
  co
))

system(paste0(
  'gdal vector rasterize -i slabs_sf.shp -o zone_r_big_desc.tiff ',
  '--dialect SQLITE --sql "SELECT * FROM slabs_sf ORDER BY slap_id DESC" ',
  '--extent ', ext_string, ' --resolution ', res_string,
  ' --overwrite --ot Int32 -a slap_id --optimization RASTER ',
  co
))


# r_rounded: read multiple times (exact_extract x2, raster calc) → tile it
tic("tic direct gdal approach:")
system(paste0(
  'gdal raster calc -i "A=', ras_path, '"',
  ' --calc "rint(A * ', rfactor, ') / ', rfactor, '"',
  ' -o r_rounded.tiff --overwrite --ot Float64 --nodata -9999 ',
  co
))
toc()

r_rounded <- terra::rast("r_rounded.tiff")

slap_means <- exact_extract(r_rounded, slaps_buffer_sf, fun = function(values, coverage) {
  median(values[coverage >= 0.5], na.rm = TRUE)
})

#tic()
#slap_means <- exactextractr::exact_extract(r_rounded, slaps_buffer_sf, fun = "median")
#toc()

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

# Tmean_raster: read once by raster calc → tile it
tic()
system(paste0(
  'gdal raster reclassify -i zone_r_big_asc.tiff -o Tmean_raster_asc.tiff ',
  '--datatype Float64 --overwrite -m "', lookup_str, '" ',
  co
))
toc()


tic()
system(paste0(
  'gdal raster reclassify -i zone_r_big_desc.tiff -o Tmean_raster_desc.tiff ',
  '--datatype Float64 --overwrite -m "', lookup_str, '" ',
  co
))
toc()


tic("One large opperation:")
system(paste0(
  'gdal raster calc ',
  '-i "D=Tmean_raster_desc.tiff" ',
  '-i "B=Tmean_raster_asc.tiff" ',
  '-i "A=r_rounded.tiff" ',
  '--calc "((D > B ? D : B) - A >= ', delta_C, ') * (A != -9999) ? 1 : NaN" ',
  '-o binary_out.tif --overwrite --ot Float32 ',
  co
))
toc()


binary <- terra::rast("binary_out.tif")


tic()
patches_v <- terra::as.polygons(binary, dissolve = TRUE, eight = TRUE)
toc()


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

tic()
stats_matrix <- exact_extract(r, patches_large, c("mean", "min", "max", "median"))
toc()

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


Tmean_raster <- terra::rast("Tmean_raster.tiff")

tic()
slap_means_per_poly <- exact_extract(Tmean_raster, patches_large_w_stats, fun = function(values, coverage) {
  median(values[coverage >= 0.5], na.rm = TRUE)
})

toc()

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
current_polys <- sf::read_sf("../current_code/final_polys.shp")


current_polys %>% colnames()
patches_large_refiltered %>% colnames()
toc()


library(tmap)

tmap_mode("view")  # interactive; use "plot" for static
tm_shape(slabs_buffer) + tm_polygons() +
tm_shape(current_polys) +
  tm_polygons(fill = "blue") +
tm_shape(patches_large_w_stats) +
  tm_polygons(fill = "yellow") +
tm_title("Patches with stats vs. current polys") +
tm_scalebar() +
tm_compass()




