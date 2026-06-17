library(terra) # for handling rasterdata
library(sf) # for handling vectordata
library(tidyverse) # generally always usefull pipe opperator, mutate etc.
library(tictoc) # for timing
library(tmap) # ploting data
library(exactextractr) # extracting statistics from rasters based on polygon footprints

# load the CWP detetion algorithm
source("../../../functions/detect_cwp_v3.R")

terraOptions(memmax=19)
tic() # start the timer
patches_new <- detect_cwp(
    ras_path = file.path("../../../data/original_data/thermal_rasters_FINAL/mean_v01emme.tif"), 
    line_path = file.path("../../../data/original_data/Centerlines_FINAL/Emme_V01.shp"),
    delta_C           = 0.5,
    min_patch_area_m2 = 2
)

# write the result to the folder
sf::st_write(patches_new, "final_polys.shp", append = FALSE)

toc() # stop the timer

