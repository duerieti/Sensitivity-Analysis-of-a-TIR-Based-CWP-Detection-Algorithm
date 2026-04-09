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

slabs_sf %>% nrow()

e <- ext(r)
ext_string <- paste(e$xmin, e$ymin, e$xmax, e$ymax, sep = ",")
resolution <- res(r)
res_string <- paste(resolution, collapse = ",")


tic()
slap_means <- exactextractr::exact_extract(r, slaps_buffer_sf, fun = "median")
toc()

Tref_vec <- round(slap_means * rfactor) / rfactor


for (j in seq_len(nrow(slabs_sf))) {
  Tref_j <- Tref_vec[j]
  if (is.na(Tref_j)) next
  
  # Rasterize slab j at FULL extent (same grid as r_rounded)
  system(paste0(
    'gdal vector rasterize -i slabs_sf.shp -o tmp_tref.tif ',
    '--dialect SQLITE --sql "SELECT * FROM slabs_sf WHERE slap_id = ', j, '" ',
    '--extent ', ext_string, ' --resolution ', res_string,
    ' --overwrite --ot Float32 --burn ', Tref_j, ' ', co
  ))
  
  # Flag: (Tref - T_rounded) >= delta, nodata where slab doesn't exist
  system(paste0(
    'gdal raster calc ',
    '-i "T=tmp_tref.tif" -i "A=r_rounded.tiff" ',
    '--calc "((T - A) >= ', delta_C, ') * (T != 0)" ',
    '-o tmp_flag.tif --overwrite --ot Byte --nodata 255 ', co
  ))
  
  # OR into accumulator
  system(paste0(
    'gdal raster calc ',
    '-i "F=flag_accum.tif" -i "N=tmp_flag.tif" ',
    '--calc "max(F, N)" ',
    '-o flag_accum_new.tif --overwrite --ot Byte --nodata 255 ', co
  ))
  
  file.rename("flag_accum_new.tif", "flag_accum.tif")
}

  # ── 2. COMPUTE REFERENCE TEMPERATURES ────────────────────────────────────

# Rasterize each slab as a separate layer in one SpatRaster stack
mask_layers <- terra::rasterize(
  terra::vect(slabs_sf),
  r,
  field = "slap_id",
  by = "slap_id"        # key: creates one band per unique slap_id
)
# Result: nlyr(mask_layers) == number of slabs
# Each band has the slap_id value inside its slab, NA outside

# Convert each band to a 1/NA binary mask
mask_binary <- !is.na(mask_layers)  # TRUE/FALSE → 1/0 per band



tic()
slap_means <- exactextractr::exact_extract(r, slaps_buffer_sf, fun = "median")
toc()

Tref_vec <- round(slap_means * rfactor) / rfactor

Tref_stack <- mask_binary * Tref_vec  

system(paste0(
  'gdal raster calc -i "A=../../data/thermal_rasters_FINAL/mean_v01emme.tif" ',
  '--calc "rint(A * ', rfactor, ') / ', rfactor, '"',
  ' -o r_rounded.tiff --overwrite --ot Float64 --nodata -9999 ',
  co
))

r_rounded <- terra::rast("r_rounded.tiff")

# ── 4. Flag: per-band (Tref_j - T_rounded >= delta) AND inside slab j ───────

delta_stack <- Tref_stack - r_rounded     # broadcasts r_rounded across all bands
flag_stack  <- (delta_stack >= delta_C) & mask_binary


# ── 5. Collapse: OR across all bands ────────────────────────────────────────

flag_final <- terra::app(flag_stack, fun = "max", na.rm = TRUE)
# 1 wherever ANY slab flagged the pixel, 0 otherwise

# Clean up: set 0 → NA for polygonization
flag_final[flag_final == 0] <- NA

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
  'gdal raster calc -i "A=../../data/thermal_rasters_FINAL/mean_v01emme.tif" ',
  '--calc "rint(A * ', rfactor, ') / ', rfactor, '"',
  ' -o r_rounded.tiff --overwrite --ot Float64 --nodata -9999 ',
  co
))
toc()


r_rounded <- terra::rast("r_rounded.tiff")



tic()
slap_means <- exactextractr::exact_extract(r, slaps_buffer_sf, fun = "median" )
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




system(paste0(
  'gdal raster calc -i "A=Tmean_raster_desc.tiff" -i "B=Tmean_raster_asc.tiff" -o Tmean_raster.tiff ',
  '--calc "A > B ? A : B" ',
  '--overwrite --ot Float32 ',
  co
))




# binary_out: read by terra::as.polygons and exact_extract → tile it
tic()
system(paste0(
  'gdal raster calc ',
  '-i "A=../../data/thermal_rasters_FINAL/mean_v01emme.tif" ',
  '-i "B=Tmean_raster.tiff" ',
  '--calc "((B - A) >= ', delta_C, ') * (A != -9999) ? 1 : NaN" ', # B has a bit of a larger extent than A, so in order to garantuee that everything works out (A != 0) is needed
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
  filter(as.numeric(st_area(.)) >= 2) %>%
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
  )



Tmean_raster <- terra::rast("Tmean_raster.tiff")

tic()
slap_means_per_poly <- exact_extract(Tmean_raster, patches_large_w_stats, fun = "median") %>% round(1)
toc()

slap_means_df <- tibble(slap_mean = slap_means_per_poly, ID = patches_large_w_stats$ID)

patches_large_w_stats <- patches_large_w_stats %>%
  inner_join(
    slap_means_df,
    by = join_by(ID == ID)
  ) %>%
  filter(slap_mean - median >= delta_C)

current_polys <- sf::read_sf("../current_code/final_polys.shp")

toc()


library(tmap)

tmap_mode("view")  # interactive; use "plot" for static

tm_shape(current_polys) +
  tm_polygons(fill = "blue") +
tm_shape(patches_large_w_stats) +
  tm_polygons(fill = "yellow") +
tm_title("Patches with stats vs. current polys") +
tm_scalebar() +
tm_compass()



