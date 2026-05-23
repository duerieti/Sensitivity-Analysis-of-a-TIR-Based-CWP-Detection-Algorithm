library(terra) # for handling rasterdata
library(sf) # for handling vectordata
library(tidyverse) # generally always usefull for pipe operator, mutate, etc.
library(tictoc) # for timing
library(tmap) # ploting data
library(exactextractr) # extracting statistics from rasters based on polygon footprints

# load the CWP detection algorithm
source("../functions/detect_cwp_v2.R")

# parse the arguments passed to the R script and store in a vector
args <- commandArgs(trailingOnly = TRUE)

# first argument is a simple index
index 	  <- as.integer(args[1])
# second argument is the number of pixels in the buffer
buffer_px <- as.numeric(args[2])
# third argument is the step length in the slab production process
step_m    <- as.numeric(args[3])

# set the memory terra is allowed to consume maximally
terraOptions(memmax=8.9)

# job_id is simply the list of arguments passed seperated by _
job_id <- paste(args, collapse = "_")

# create a file path pointing to a folder where all intermediate results will be written
# and the final result will be writen to
job_dir <- file.path(getwd(), "results", paste0("job_", job_id))

# create this directory
dir.create(job_dir, recursive = TRUE, showWarnings = FALSE)

# set the path of the interpreter to this directory
setwd(job_dir)

# devine the path to the raster
ras_path <- file.path("../../../data/thermal_rasters_FINAL/mean_v2emme.tif")
# define the path to the rounded version of the raster
# (this is a bit optimised for the scanning runs. The raster only needs to be rounded
# once. Its not necessairy that every worker rounds the raster internally, when this can be
# be done one a one-time-basis)
rounded_ras_path <- file.path("../../../data/r_rounded_aligned.tiff") # needs to be changed too
# define the path to the line
line_path <- file.path("../../../data/Centerlines_FINAL/Emme_V02.shp")

# start the function with the arguments passed from GNU parallel
patches_new <- detect_cwp_single(
    ras_path, line_path, rounded_ras_path,
    segment_length_m = step_m,
    buffer_px        = buffer_px
)

sf::st_write(patches_new, "final_polys.shp", append = FALSE)
