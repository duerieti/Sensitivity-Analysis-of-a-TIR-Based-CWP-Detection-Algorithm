library(terra) # for handling rasterdata
library(sf) # for handling vectordata
library(tidyverse) # generally always usefull pipe opperator, mutate etc.
library(tictoc) # for timing
library(tmap) # ploting data
library(exactextractr) # extracting statistics from rasters based on polygon footprints

# load the CWP detetion algorithm
source("../../../functions/detect_cwp_v2_internal_rounding.R")
terraOptions(memmax=19)
tic() # start the timer
patches_new <- detect_cwp_single(
    ras_path = file.path("../../../data/original_data/thermal_rasters_FINAL/mean_v01emme.tif"), 
    line_path = file.path("../../../data/original_data/Centerlines_FINAL/Emme_V01.shp"),
    segment_length_m  = 1200,   # step length for segmenting the centerline
    zone_halfwidth_m  = 60,    # halfwidth of the wide zone slabs used for rasterisation
    buffer_px         = 3,
    delta_C           = 1.0,
    min_patch_area_m2 = 2,
    round_to          = 0.1,   # rounding precision for the TIR raster (°C)
    connect_diagonals = TRUE   # whether diagonally touching pixels form one polygon
)

# write the result to the folder
sf::st_write(patches_new, "final_polys.shp", append = FALSE)

toc() # stop the timer

