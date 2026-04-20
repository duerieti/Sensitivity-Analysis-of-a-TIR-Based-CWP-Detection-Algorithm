library(tidyverse)
library(tmap)
library(sf)

crop_to_regions <- function(polys) {
  segment_limpbach  <- st_crop(polys, bbox_limpbach)  %>% mutate(location = "limpbach")
  segment_grundbach    <- st_crop(polys, bbox_grundbach)    %>% mutate(location = "grundbach")
  segment_hasle     <- st_crop(polys, bbox_hasle)     %>% mutate(location = "hasle")
  segment_emmepark  <- st_crop(polys, bbox_emmepark)  %>% mutate(location = "emmepark")
  list(segment_limpbach, segment_grundbach, segment_hasle, segment_emmepark)
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

bbox_limpbach <- st_bbox(c(
  xmin = 2607841,
  ymin = 1223113,
  xmax = 2608876,
  ymax = 1223590
), crs = st_crs(2056))

bbox_grundbach <- st_bbox(c(
  xmin = 2607970,
  ymin = 1220090,
  xmax = 2608178,
  ymax = 1220820
), crs = st_crs(2056))

bbox_hasle <- st_bbox(c(
  xmin = 2616141,
  ymin = 1207441,
  xmax = 2616417,
  ymax = 1207772
), crs = st_crs(2056))

bbox_emmepark <- st_bbox(c(
  xmin = 2607871,
  ymin = 1222261,
  xmax = 2608128,
  ymax = 1222754
), crs = st_crs(2056))

file_and_folder_names <- list.files("./results")

folder_names <- file_and_folder_names %>%
  .[str_detect(., ".txt", negate = TRUE)]

lumped_stats_per_location <- map(
  folder_names,
  ~process_polygon(., "results")
) %>% bind_rows()


write_csv(lumped_stats_per_location, "lumped_stats_per_location.csv")


install.packages("sensitivity")
