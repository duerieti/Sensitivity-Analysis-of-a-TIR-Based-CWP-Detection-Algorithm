library(terra)
library(sf)
library(dplyr)
library(tidyverse)
library(tictoc)
library(tmap)
library(exactextractr)

source("single_cwp.R")

args <- commandArgs(trailingOnly = TRUE)


index 	  <- as.integer(args[1])
buffer_px <- as.numeric(args[2])
step_m    <- as.numeric(args[3])

terraOptions(memmax=19)


job_id <- paste(args, collapse = "_")

job_dir <- file.path(getwd(), "results", paste0("job_", job_id))
dir.create(job_dir, recursive = TRUE, showWarnings = FALSE)
setwd(job_dir)


ras_path <- file.path("../../../data/thermal_rasters_FINAL/mean_v2obemme.tif")
rounded_ras_path <- file.path("../../../data/r_v2obemme_aligned.tif") # needs to be changed too
line_path <- file.path("../../../data/Centerlines_FINAL/Obere_Emme_V02.shp")

patches_new <- detect_cwp_single(
    ras_path, line_path, rounded_ras_path,
    step_m           = step_m,
    buffer_px        = buffer_px,
)

sf::st_write(patches_new, "final_polys.shp", append = FALSE)
