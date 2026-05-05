library(tidyverse)
library(tmap)
library(sf)

getwd()
setwd("./compute_morris")

bbox_ilfis <- st_bbox(c(
  xmin =2623723,
  ymin =1199977,
  xmax =2623793,
  ymax =1200424

))

bbox_zollbrueck <- st_bbox(c(
  xmin =2622726,
  ymin =1201647,
  xmax =2623656,
  ymax =1203076

))

bbox_schuepbachkanal <- st_bbox(c(
  xmin =2622675,
  ymin =1197452,
  xmax =2622838,
  ymax =1197726
))


bbox_schuepbachkanal <- st_bbox(c(
  xmin =2622675,
  ymin =1197452,
  xmax =2622838,
  ymax =1197726
))


bbox_schuepbach <- st_bbox(c(
  xmin =2622828,
  ymin =1196960,
  xmax =2623227,
  ymax =1197490
))


crop_to_regions <- function(polys) {

  segment_ilfis                <- st_crop(polys, bbox_ilfis) %>% mutate(location = "ilfis")
  segment_zollbrueck           <- st_crop(polys, bbox_zollbrueck) %>% mutate(location = "zollbrueck")
  segment_schuepbachkanal      <- st_crop(polys, bbox_schuepbachkanal) %>% mutate(location = "schuepbachkanal")
  segment_schuepbach           <- st_crop(polys, bbox_schuepbach) %>% mutate(location = "schuepbach")

  list(segment_ilfis, segment_zollbrueck,segment_schuepbachkanal,segment_schuepbach)
}

compute_lumped_statistic <- function(poly_region) {
  lumped_stats <- poly_region %>%
    sf::st_drop_geometry() %>%
    mutate(
      weighted_T_mean = T_mean * area_m2
    ) %>%
    summarise(
      total_area        = sum(area_m2),
      total_averaged_temperature = sum(weighted_T_mean) / total_area,
      morris_row_index  = unique(morris_row_index),
      location     = unique(location)
    )
  
  return(lumped_stats)
}

process_polygon <- function(folder_name, base_path) {
  index     <- as.integer(str_split(folder_name, "_")[[1]][2])
  read_path <- file.path(base_path, folder_name, "final_polys.shp")


  region_polys <- read_path %>%
    sf::read_sf() %>%
    mutate(morris_row_index = index) %>%
    crop_to_regions()

  lumped_stats_per_region <- map(region_polys, compute_lumped_statistic) %>%
    bind_rows()

  return(lumped_stats_per_region)
}



file_and_folder_names <- list.files("./results_upper_emme")

folder_names <- file_and_folder_names %>%
  .[str_detect(., ".txt", negate = TRUE)]



lumped_stats_per_location <- map(
  folder_names,
  ~process_polygon(., "results_upper_emme")
) %>% bind_rows()


write_csv(lumped_stats_per_location, "lumped_stats_upper_emme.csv")



