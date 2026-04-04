library(terra)
library(sf)
library(dplyr)
library(tidyverse)
library(tictoc)
library(tmap)
library(exactextractr)
install.packages("terra")


tic("Total Time Chunked workflow: ")
ras_path <- file.path("data/thermal_rasters_FINAL/mean_v01emme.tif")
line_path <- file.path("data/Centerlines_FINAL/Emme_V01.shp")

    step_m            = 500
    buffer_px         = 3
    delta_C           = 1.0
    min_patch_area_m2 = 2
    round_to          = 0.1
    slab_halfwidth_m  = 60
    connect_diagonals = TRUE
    union_chunk_size  = 50000L

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
slabs_sf <- sf::st_sf(geometry = slabs) %>% mutate(slap_id = row_number())
sf::write_sf(slabs_sf, "slabs_sf.shp")

buffer_m <- cell_m * buffer_px

slabs_buffer <- sf::st_buffer(subs,
                              dist        = buffer_m,
                              endCapStyle = "FLAT",
                              joinStyle   = "MITRE",
                              mitreLimit  = 2)
slaps_buffer_sf <- sf::st_sf(geometry = slabs_buffer) %>% mutate(slap_id = row_number())

  # ── 2. COMPUTE REFERENCE TEMPERATURES ────────────────────────────────────

tic("tic direct gdal approach:")
system(paste0(
  'gdal raster calc -i "A=data/thermal_rasters_FINAL/mean_v01emme.tif" ',
  '--calc "A*', rfactor, '/', rfactor, '" ',
  '-o r_rounded.tiff --overwrite --ot Float32  --nodata -9999 ',
  co
))
toc()

r_rounded <- terra::rast("r_rounded.tiff")

slap_means <- exactextractr::exact_extract(r_rounded, slaps_buffer_sf, fun = "median") %>% round(1)

temp_look_up_df <- tibble(slap_id = slaps_buffer_sf$slap_id, slap_means = slap_means)

  # ── 3. FLAG COLD PIXELS WITH OR LOGIC ACROSS OVERLAPPING SLABS ───────────

tic("OR-logic pixel flagging: ")


# Extract all cell-to-slab mappings — duplicates intentional for overlap zones
cell_slab_df <- exactextractr::exact_extract(
  r_rounded,
  slabs_sf,                  # use full slabs (not corridors) to match original logic
  fun = NULL,
  include_cell = TRUE,
  include_cols = "slap_id",
  progress = FALSE
) |> dplyr::bind_rows()


# Join reference temperatures, flag per (cell, slab) pair, then OR across slabs
flagged_cells <- cell_slab_df |>
  dplyr::filter(!is.na(value)) |>
  dplyr::left_join(temp_look_up_df, by = "slap_id") |>
  dplyr::mutate(flagged = (slap_means - value) >= delta_C) |>
  dplyr::group_by(cell) |>
  dplyr::summarise(flag = any(flagged, na.rm = TRUE), .groups = "drop") |>
  dplyr::filter(flag) |>
  dplyr::pull(cell)

toc()


  # ── 4. BURN FLAGGED CELLS INTO BINARY RASTER ─────────────────────────────

tic("Build binary raster: ")

xy_flagged <- terra::xyFromCell(r_rounded, flagged_cells) |>
  as.data.frame() |>
  dplyr::mutate(flag = 1L)

print("creating the vector")
pts <- terra::vect(
  xy_flagged, 
  geom  = c("x", "y"), 
  crs   = terra::crs(r_rounded)
)

terra::writeVector(pts, "points.shp", overwrite=TRUE)

print("trying the rasterisation")
terra::rasterize(
  pts,
  r_rounded,
  field    = "flag",
  filename = "binary_out.tif",
  overwrite = TRUE,
  datatype  = "INT1U",
  NAflag    = 255,
  gdal      = c("TILED=YES", "BLOCKXSIZE=512", "BLOCKYSIZE=512", "BIGTIFF=YES")
)


tic("tic direct gdal approach:")
system(paste0(
  'gdal vector rasterize -i points.shp -o binary_out --burn 1 --overwrite --datatype Int8 --nodata 255',
  co
))
toc()


toc()

getwd()

print("load")
binary <- terra::rast("binary_out.tif")



  # ── 5. POLYGONIZE ─────────────────────────────────────────────────────────

tic()
patches_v <- terra::as.polygons(binary, dissolve = TRUE, eight = TRUE)
toc()

eps <- if (connect_diagonals) cell_m * 0.1 else 0

patches_sf <- patches_v %>%
  sf::st_as_sf() %>%
  sf::st_cast("POLYGON")
  
  # ── 6. AREA FILTER ────────────────────────────────────────────────────────

patches_large <- patches_sf %>%
  mutate(area_m2 = as.numeric(st_area(.))) %>%
  filter(area_m2 >= min_patch_area_m2) %>%
  mutate(ID = row_number())

  # ── 7. TEMPERATURE STATISTICS PER PATCH ──────────────────────────────────

tic()
stats_matrix <- exact_extract(r_rounded, patches_large, c("mean", "min", "max", "median"))
toc()



stats_per_poly <- stats_matrix %>%
  as_tibble() %>%
  mutate(ID = patches_large$ID) %>%
  mutate(
    T_min = min,
    T_max = max,
    T_med = median,
    T_mean = mean
  ) %>% select(!c(min, max, median, mean))



median_of_slab <- exact_extract(, patches_large, "median", force_df = TRUE) %>%
  mutate(
    ID = row_number(),
    median = round(median,1)
) %>%
  rename(median_of_slab = median)
  


patches_large_with_stats <- patches_large %>%
  inner_join(stats_per_poly, by = join_by(ID == ID)) %>%
  inner_join(median_of_slab, by = join_by(ID == ID)) %>%
  mutate(delta_T = median_of_slab - T_med)


sf::write_sf(patches_large_with_stats, "final_polys.shp")


polys_current <- sf::read_sf("./comparison/current_code/final_polys.shp")
polys_new <- sf::read_sf("final_polys.shp")
rast_binary <- terra::rast("binary_out.tif")




tmap_mode("view")


tm_shape(polys_current, name = "Current Polygons") +
  tm_polygons(
    fill = "T_mean",
    fill.scale = tm_scale_continuous(values = "viridis"),
    fill_alpha = 0.7,
    col = "darkorange",
    lwd = 1.5,
    fill.legend = tm_legend(title = "Current Code (ID)")
  ) +
  tm_text("T_mean", size = 1.5, col = "darkorange", fontface = "bold",
          xmod = -0.002, ymod = 0.002) +   # shift UP-LEFT
tm_shape(polys_new, name = "New Polygons") +
  tm_polygons(
    fill = "T_mean",
    fill.scale = tm_scale_continuous(values = "plasma"),
    fill_alpha = 0.7,
    col = "darkblue",
    lwd = 1.5,
    fill.legend = tm_legend(title = "Optimized Code (ID)")
  ) +
  tm_text("T_mean", size = 1.5, col = "darkblue", fontface = "bold",
          xmod = 0.002, ymod = -0.002) +   # shift DOWN-RIGHT
tm_basemap(c(
    "OpenStreetMap" = "OpenStreetMap",
    "Satellite"     = "Esri.WorldImagery",
    "Topo"          = "OpenTopoMap"
  )) +
  tm_scalebar(position = c("left", "bottom")) +
  tm_compass(position = c("right", "top")) +
  tm_title("Layer Comparison: Current vs Optimized (colored by ID)")

