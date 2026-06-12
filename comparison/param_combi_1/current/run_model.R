library(terra) # for handling rasterdata
library(sf) # for handling vectordata
library(tidyverse) # generally always usefull pipe opperator, mutate etc.
library(tictoc)    # for timing
library(tmap) # ploting data
library(exactextractr) # extracting statistics from rasters based on polygon footprints

# load the CWP detetion algorithm
source("../../../functions/function_current.R")

terraOptions(memmax=19)

tic() # start timer
patches_new <- detect_cwp_single(
    ras_path = file.path("../../../data/original_data/thermal_rasters_FINAL/mean_v01emme.tif"), 
    line_path = file.path("../../../data/original_data/Centerlines_FINAL/Emme_V01.shp"),
    step_m            = 200,       # slab length along centerline (m)
    slab_halfwidth_m  = 60,
    buffer_px         = 4,         # corridor half-width for reference T (pixels)
    delta_C           = 0.5,       # temperature anomaly threshold (deg C)
    min_patch_area_m2 = 2,         # drop patches smaller than this
    round_to          = 0.1,       # round temperatures to this precision
    connect_diagonals = TRUE      # 8-connectivity when dissolving patch
    )

# write the result to the folder
sf::st_write(patches_new, "final_polys.shp", append = FALSE)

toc() # stop timer

