library(terra) # for handling rasterdata
library(sf) # for handling vectordata
library(tidyverse) # generally always usefull. pipe operator, mutate etc.
library(tictoc) # for timing
library(tmap) # for ploting data
library(exactextractr) # extracting statistics from rasters based on polygon footprints

# load the CWP detection algorithm
source("../functions/detect_cwp_v2.R")

# fetch the possitional that are passed to this R script
args <- commandArgs(trailingOnly = TRUE)
# convert the arguments to nummeric types and store in variable
index     <- as.integer(args[1])
buffer_px <- as.numeric(args[2])
step_m    <- as.numeric(args[3])
delta_t   <- as.numeric(args[4])

# set a quota to how much memory terra is allowed to use
terraOptions(memmax = 9)

# create a id for the job. This is just a string containing the arguments that where passed to the script eg "1_7_760_1.5"
job_id  <- paste(args, collapse = "_")
# create a filepath describing the directory, where the final results will be stored and intermediate results will be written to. The folder name contains the job id e.g. "job_1_7_760_1.5"
job_dir <- file.path(getwd(), "results", paste0("job_", job_id))

# if the directory allready exists, then skipt the computation entirely
if (dir.exists(job_dir)) {
  message("Job directory already exists, skipping: ", job_dir)
  quit(save = "no", status = 0)
}
# else the directory will be created
dir.create(job_dir, recursive = TRUE, showWarnings = FALSE)

# move the interpreter path to the directory (all following paths will now be relative to the job directory
setwd(job_dir)

# read the mean temperature raster fo the "emme v2" segment
ras_path         <- file.path("../../../data/original_data/thermal_rasters_FINAL/mean_v2emme.tif")

# read the roundev version of the temperature raster
rounded_ras_path <- file.path("../../../data/derived_data_products/mean_v2emme_rounded.tif")

# read the centerline for the "emme v2" segment
line_path        <- file.path("../../../data/original_data/Centerlines_FINAL/Emme_V02.shp")

# Run the CWP detetion algorithm with the parammeters that where passed to the R script by GNU parallel

patches_new <- detect_cwp_single(
  ras_path, line_path, rounded_ras_path,
  segment_length_m     = step_m,
  buffer_px            = buffer_px,
  delta_C              = delta_t
)

# write the result to the job directory
sf::st_write(patches_new, "final_polys.shp", append = FALSE)
