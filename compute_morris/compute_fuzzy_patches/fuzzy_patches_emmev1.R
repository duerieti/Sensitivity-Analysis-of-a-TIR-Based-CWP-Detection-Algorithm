library(tidyverse)
library(tmap)
library(terra)

getwd()
setwd("compute_morris")

file_and_folder_names <- list.files("../../morris_screening_emme_v1/results")
folder_names <- file_and_folder_names %>%
  .[str_detect(., ".txt", negate = TRUE)]

read_paths <- folder_names %>% file.path("../../morris_screening_emme_v1/results", ., "final_polys.shp")

r <- terra::rast("../data/original_data/thermal_rasters_FINAL/mean_v01emme.tif")
e <- ext(r)
ext_string <- paste(e$xmin, e$ymin, e$xmax, e$ymax, sep = ",")
resolution <- res(r)
res_string <- paste(resolution, collapse = ",")

dir.create("./binary_rasters_emme_v1", showWarnings = FALSE)

for (i in seq_along(read_paths)) {
  out_path <- file.path("./binary_rasters_emme_v1", paste0("rast_", i, ".tif"))
  system(paste0(
    'gdal vector rasterize -i ', read_paths[i], ' -o ', out_path, ' ',
    '--extent ', ext_string, ' --resolution ', res_string,
    ' --burn 1 --init 0 ',
    ' --overwrite --ot Byte --optimization RASTER'
  ))
}



# Step 1: Build a virtual raster stack (instant, no data copied)
raster_paths <- list.files("./binary_rasters_emme_v1", pattern = "\\.tif$", full.names = TRUE)


raster_paths_string <- paste0(raster_paths, collapse = " ")


system(
  paste0(
    'gdal raster stack ',raster_paths_string, ' stacked_raster.vrt --overwrite'
  )
  
)


system(
  paste0(
    'gdal raster calc -i "A=stacked_raster.vrt" --flatten --calc "avg(A)" -o fuzzy_cwp_emmev1.tif'
  )
  
)
