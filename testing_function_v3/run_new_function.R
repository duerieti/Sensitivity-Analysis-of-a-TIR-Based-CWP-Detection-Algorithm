setwd("testing_function_v3")

ras_path  <- "../data/original_data/thermal_rasters_FINAL/mean_v01emme.tif"
line_path <- "../data/original_data/Centerlines_FINAL/Emme_V01.shp"

source("../functions/detect_cwp_v3.R")

cwp_new <- detect_cwp(ras_path, line_path, n_cores = 12)


sf::write_sf(cwp_new, "cwp_of_function_v3.shp")

tmap_mode("view")
tm_shape(cwp_new) + tm_polygons()

source("../functions/detect_cwp_v2_internal_rounding.R")

cwp_odl <-  detect_cwp_single(ras_path, line_path)


tmap_mode("view")
tm_shape(cwp_new) + tm_polygons(col = "red", alpha = 0.5) +
  tm_shape(cwp_odl) + tm_polygons(col = "blue", alpha = 0.5)


source("../functions/jaccard.R")
jaccard_similarity(cwp_new, cwp_odl)
