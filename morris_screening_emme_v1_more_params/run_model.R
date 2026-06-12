library(terra) # for handling rasterdata
library(sf) # for handling vectordata
library(tidyverse) # generally always usefull. pipe operator, mutate etc.
library(tictoc) # for timing
library(tmap) # ploting data
library(exactextractr) # extracting statistics from rasters based on polygon footprints

# load the CWP detection algorithm 
source("../functions/detect_cwp_v2.R")

# parse the arguments passed to the R script and store it in a vector
args <- commandArgs(trailingOnly = TRUE)

# first argument is a simple index
index     <- as.integer(args[1])
# second argument is the number of pixels in the buffer
buffer_px <- as.numeric(args[2])
# third argument is the step length in the slab production process
step_m    <- as.numeric(args[3])
# forth one is the temperature delta used for pixel flagging
delta_t   <- as.numeric(args[4])

# set the memory terra is allowed to consume maxiammly
terraOptions(memmax = 14)

# job_id is simply the list of arguments passed separated by _
job_id  <- paste(args, collapse = "_")

# create a file path pointing to a folder where all intermediate results will be written
# and the final results will be written to.
job_dir <- file.path(getwd(), "results", paste0("job_", job_id))

# if the directory path allready exists, then the job allready has been completed. => terminate the process
if (dir.exists(job_dir)) {
  message("Job directory already exists, skipping: ", job_dir)
  quit(save = "no", status = 0)
}

# create the job directory
dir.create(job_dir, recursive = TRUE, showWarnings = FALSE)

# set the path of the interpreter to this directory
setwd(job_dir)

# define the path to the raster
ras_path         <- file.path("../../../data/original_data/thermal_rasters_FINAL/mean_v01emme.tif")
# define the path to the rounded version of the raster
# (this is a bit optimised for the scanning runs. The raster only needs to be rounded
# once. Its not necessairy that every worker rounds the raster internally, when this can be
# be done one a one-time-basis)
rounded_ras_path <- file.path("../../../data/derived_data_products/mean_v01emme_rounded.tif")
# define the path to the line
line_path        <- file.path("../../../data/original_data/Centerlines_FINAL/Emme_V01.shp")


patches_new <- detect_cwp_single(
  ras_path, line_path, rounded_ras_path,
  segment_length_m  = step_m,
  buffer_px         = buffer_px,
  delta_C           = delta_t
)

# write the results to the job directory
sf::st_write(patches_new, "final_polys.shp", append = FALSE)
