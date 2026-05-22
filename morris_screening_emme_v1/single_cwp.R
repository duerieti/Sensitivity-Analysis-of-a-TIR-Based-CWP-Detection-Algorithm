# detecting cold water patches in rivers
detect_cwp_single <- function(

    ras_path, # path to the TIR raster of the river
    line_path, # path to the centerline of the river
    rounded_ras_path, # path to the rounded version of the TIR raster
    step_m            = 500, # step lenght for producing the slabs following the river
    buffer_px         = 3, # number of pixels for the buffer calculation
    delta_C           = 1.0, # temperature delta for pixel flagging
    min_patch_area_m2 = 2, # minimal surface area of the CWP's
    round_to          = 0.1, # rounding precission for the raster
    slab_halfwidth_m  = 60, # halfwidth of the slabs 
    connect_diagonals = TRUE, # wheter to connect diagonally touching pixels in polygonisations to the same polygon
    union_chunk_size  = 50000L # dont know wether this is even necessairy

) {

# tiling creation options — no compression, just chunked layout
# (tiling improves the performance of exact extract)
co <- "--co TILED=YES --co BLOCKXSIZE=512 --co BLOCKYSIZE=512 --co BIGTIFF=YES"

  # ── 0. READ INPUTS ────────────────────────────────────────────────────────

# load the raster from the raster path
r <- terra::rast(ras_path)

# load the line from the line path
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

# creating the rounding factor
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

# get the extent from the TIR raster
e <- ext(r)
# create a string representation of the extent which can then be passed to GDAL
ext_string <- paste(e$xmin, e$ymin, e$xmax, e$ymax, sep = ",")

# get the resolution of the TIR raster
resolution <- res(r)
# create a string representation of the resolution which can then be passed to GDAL
res_string <- paste(resolution, collapse = ",")

## perform rasterisation of the slabs. ##
  
# each slab-polygon gets rasterized. As a burn in value the slab id is used. This way it is 
# determinable from which slab each pixel came from. if two slabs overlap the one with the
# higher slab id will win.

system(paste0(
  'gdal vector rasterize -i slabs_sf.shp -o zone_r_big_asc.tiff ',
  '--dialect SQLITE --sql "SELECT * FROM slabs_sf ORDER BY slap_id ASC" ',
  '--extent ', ext_string, ' --resolution ', res_string,
  ' --overwrite --ot Int32 -a slap_id --optimization RASTER ',
  co
))


## perform rasterisation of the slabs. ##
  
# each slab-polygon gets rasterized. As a burn in value the slab id is used. This way it is 
# determinable from which slab each pixel came from. if two slabs overlap the one with the
# lower slab id will win.
  
system(paste0(
  'gdal vector rasterize -i slabs_sf.shp -o zone_r_big_desc.tiff ',
  '--dialect SQLITE --sql "SELECT * FROM slabs_sf ORDER BY slap_id DESC" ',
  '--extent ', ext_string, ' --resolution ', res_string,
  ' --overwrite --ot Int32 -a slap_id --optimization RASTER ',
  co
))

# note: These two operations will produce the exact same raster with the only difference that the issue of overlapping rasters
# is handled in exactly opposing ways.
# for one version of the reasterisation, the slab production is handled in a way where slab with the larger id will win in conflicing cases
# for the other version, the slab witht the lower id will win.

# load the rounded raster
r_rounded <- terra::rast(rounded_ras_path)

# extract the reference temperatures from the rounded raster using the buffer slabs
slap_means <- exact_extract(r_rounded, slaps_buffer_sf, fun = function(values, coverage) {
  median(values[coverage >= 0.5], na.rm = TRUE)
})

# create a look up structure as a dataframe. One row holds the reference temperature and the other the slab id,
# of the slab_buffer, the temperature was derived with.
temp_look_up_df <- tibble(slap_id = slaps_buffer_sf$slap_id, slap_means = slap_means)
  

# write the lookup table to disc as a text file. The structure has to be very specific, so that it can be used
# in gdal for rasterisation.

writeLines(
  paste0(
    paste0("[", temp_look_up_df$slap_id, ",", temp_look_up_df$slap_id, "]=", temp_look_up_df$slap_means, collapse = "; "),
    "; DEFAULT=NO_DATA"
  ),
  "median_lookup_gdal.txt"
)


  # ── 5. BURN MEDIANS INTO ZONES ────────────────────────────────────────────

# the slabs and the buffer_slabs shared the exact same IDs. Now the look up table can be used to reclassify the raster
# each cell holds the ID, from the slab it was produced. By this makes it possible to assign the reference temperatures
# computed with the buffer slab that corresponds to the slab to the pixels


# read the string into memory
lookup_str <- readLines("median_lookup_gdal.txt")

# Reclassify the ascending version of the raster to hold the reference temperatures
system(paste0(
  'gdal raster reclassify -i zone_r_big_asc.tiff -o Tmean_raster_asc.tiff ',
  '--datatype Float64 --overwrite -m "', lookup_str, '" ',
  co
))

# reclassify the descending version of the raster to hold the reference temperatures
system(paste0(
  'gdal raster reclassify -i zone_r_big_desc.tiff -o Tmean_raster_desc.tiff ',
  '--datatype Float64 --overwrite -m "', lookup_str, '" ',
  co
))


# Do some cleaning up of raster files. This is necessairy as these raster have quite a large size
# especially when then running the function in parallel, this can lead to tremendous spikes in disk space (in the orders of hunderst of GB)
file.remove(c("zone_r_big_asc.tiff", "zone_r_big_desc.tiff"))

# this merges the two currently existing reference temperature rasters.
# this is just for sorting out the overlapping regions.
# for an overlapping region the the larger temperature value for a specific pixel will always be preffered. This is to make shure that
# a pixel that belongs to two slabs will always be flagged as a cold water pixel if it has a sufficiently slower temperature to be accounted
# to the cwp-pixels in at least one of the slabs.

system(paste0(
  'gdal raster calc -i "A=Tmean_raster_desc.tiff" -i "B=Tmean_raster_asc.tiff" -o Tmean_raster.tiff ',
  '--calc "A > B ? A : B" ',
  '--overwrite --ot Float64 ',
  co
))


# perform the pixel flagging

# for every pair of corresponding pixels in A and B:
  # if the difference between the temperatures in the rounded raster and the temperature in the reference raster is larger than
  # then (B - A) >= delta_C becomes true. 
  # if the pixel in A is not -9999 (Na value), then (A != -9999) becomes true
  # true * true is true, ture * false is false, false * true is false, false * false is false
  # therefore only if the temperature difference is larger than delta_C and the value in A is not a Na value the pixel is flaged with 1
  # else the pixel is set to NaN.

system(paste0(
  'gdal raster calc ',
  '-i "A=', rounded_ras_path, '" ',
  '-i "B=Tmean_raster.tiff" ',
  '--calc "((B - A) >= ', delta_C, ') * (A != -9999) ? 1 : NaN" ', # B has a bit of a larger extent than A, so in order to garantuee that everything works out (A != -9999) is needed
  '-o binary_out.tif --overwrite --ot Float32 ',
  co
))
  
# again remove some of the intermediate files which are no longer needed
file.remove(c("Tmean_raster_desc.tiff", "Tmean_raster_asc.tiff"))

# load the raster of flagged pixels
binary <- terra::rast("binary_out.tif")

# convert binary raster into polygons. 
# dissolve = TRUE : neighboring pixels are fused to the same polygon (else each pixel would become its separate polygon)

patches_v <- terra::as.polygons(binary, dissolve = TRUE, eight = TRUE)


# ── 7. POLYGONIZE ─────────────────────────────────────────────────────────

  
# a small value 10% of the cell length/width
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

# all cold water patches smaller than min_patch_area_m2 are filtered out
patches_large <- patches_sf %>%
  mutate(area_m2 = st_area(.) %>% as.numeric()) %>%
  filter(area_m2 >= min_patch_area_m2) %>%
  mutate(ID = row_number())

  # ── 9. TEMPERATURE STATISTICS PER PATCH ──────────────────────────────────

# computing temperature statistics for the patches
stats_matrix <- exact_extract(r, patches_large, c("mean", "min", "max", "median"))

# adding the patch ID to the stats matrix
stats_per_poly <- stats_matrix %>%
  as_tibble() %>%
  mutate(ID = patches_large$ID)


# joining the patch statistics onto the patches using the id
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

# loading the reference temperature raster
Tmean_raster <- terra::rast("Tmean_raster.tiff")

# computing the refference temperature for every patch again 
slap_means_per_poly <- exact_extract(Tmean_raster, patches_large_w_stats, fun = function(values, coverage) {
  median(values[coverage >= 0.5], na.rm = TRUE)
})
  

# creating a tibble of the reference temperature values together with the slab ID's
slap_means_df <- tibble(Tmd_slb = slap_means_per_poly, ID = patches_large_w_stats$ID)


# joining the reference temperature values onto the patches using the id and filtering all
# polygons that do not fullfill the filter condition
patches_large_refiltered <- patches_large_w_stats %>%
  inner_join(
    slap_means_df,
    by = join_by(ID == ID)
  ) %>%
  mutate(
    deltaT = Tmd_slb - T_med
  ) %>%
  filter(deltaT >= delta_C) 

# remove the intermediate files
file.remove(c("binary_out.tif","Tmean_raster.tiff"))

return(patches_large_refiltered)
  
}


