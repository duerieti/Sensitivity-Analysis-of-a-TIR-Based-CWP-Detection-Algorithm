library(tidyverse)
library(tmap)
library(terra)

getwd()
setwd("compute_morris")

file_and_folder_names <- list.files("./results_upper_emme")
folder_names <- file_and_folder_names %>%
  .[str_detect(., ".txt", negate = TRUE)]

read_paths <- folder_names %>% file.path("./results_upper_emme", ., "final_polys.shp")

r <- terra::rast("../data/thermal_rasters_FINAL/mean_v2emme.tif")
e <- ext(r)
ext_string <- paste(e$xmin, e$ymin, e$xmax, e$ymax, sep = ",")
resolution <- res(r)
res_string <- paste(resolution, collapse = ",")

dir.create("./binary_rasters_upper_emme", showWarnings = FALSE)



for (i in seq_along(read_paths)) {
  out_path <- file.path("./binary_rasters_upper_emme", paste0("rast_", i, ".tif"))
  system(paste0(
    'gdal vector rasterize -i ', read_paths[i], ' -o ', out_path, ' ',
    '--extent ', ext_string, ' --resolution ', res_string,
    ' --burn 1 --init 0 ',
    ' --overwrite --ot Byte --optimization RASTER'
  ))
}



test_r <- rast("./binary_rasters_upper_emme/rast_2.tif")

polys_test <- sf::read_sf("./results_upper_emme/job_2_7_880/final_polys.shp")



# Step 1: Build a virtual raster stack (instant, no data copied)
raster_paths <- list.files("./binary_rasters_upper_emme", pattern = "\\.tif$", full.names = TRUE)


raster_paths_string <- paste0(raster_paths, collapse = " ")




system(
  paste0(
    'gdal raster stack ',raster_paths_string, ' stacked_raster_upper_emme.vrt --overwrite'
  )
  
)


system(
  paste0(
    'gdal raster calc -i "A=stacked_raster_upper_emme.vrt" --flatten --calc "avg(A)" -o fuzzy_cwp_upper_emme.tif'
  )
  
)


sets <- terra::rast("fuzzy_cwp_upper_emme.tif")

tmap_mode("view")
tm_shape(sets) +
  tm_raster()
