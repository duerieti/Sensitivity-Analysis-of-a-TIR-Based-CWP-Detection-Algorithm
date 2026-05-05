library(tidyverse)
library(tmap)
library(sf)

getwd()
setwd("./compute_morris")




classified_cwp_data <- read_delim("cwp_annotated.csv") %>%
  separate_wider_delim(cols = Bbox, delim = ",", names = c("xmin", "ymin", "xmax", "ymax")) %>%
  mutate(
    xmin = as.numeric(xmin),
    ymin = as.numeric(ymin),
    xmax = as.numeric(xmax),
    ymax = as.numeric(ymax)
  )


crop_to_regions <- function(polys, df) {
n_row <- df %>% nrow()


cropped_list <- vector("list", n_row)

  for (i in seq(1,n_row,1)) {
    row <- df[i,  ]


    bbox <- st_bbox(
      c(
      xmin = row$xmin,
      ymin = row$ymin,
      xmax = row$xmax,
      ymax = row$ymax),
    crs = st_crs(2056)
    )

    croped_region <- st_crop(polys, bbox)
    if (nrow(croped_region) != 0) {
      cropped_list[[i]] <- croped_region %>% 
        mutate(
        Class = row$Class,
        identifier = i
      )
    }
  }

  return(cropped_list)
  
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
      identifier = unique(identifier),
      Class = unique(Class)
    )
  
  return(lumped_stats)
}

process_polygon <- function(folder_name, base_path, df) {
  index     <- as.integer(str_split(folder_name, "_")[[1]][2])
  read_path <- file.path(base_path, folder_name, "final_polys.shp")


  region_polys <- read_path %>%
    sf::read_sf() %>%
    mutate(morris_row_index = index) %>%
    crop_to_regions(df) %>% 
    compact()

  lumped_stats_per_region <- map(region_polys, compute_lumped_statistic) %>%
    bind_rows()

  return(lumped_stats_per_region)
}



file_and_folder_names <- list.files("./results_lower_emme")

folder_names <- file_and_folder_names %>%
  .[str_detect(., ".txt", negate = TRUE)]

lumped_stats_per_location_lower_emme <- map(
  folder_names,
  ~process_polygon(., "results_lower_emme", classified_cwp_data)
) %>% bind_rows()


file_and_folder_names <- list.files("./results_upper_emme")

folder_names <- file_and_folder_names %>%
  .[str_detect(., ".txt", negate = TRUE)]

lumped_stats_per_location_upper_emme <- map(
  folder_names,
  ~process_polygon(., "results_upper_emme", classified_cwp_data)
) %>% bind_rows()


file_and_folder_names <- list.files("./results_obemme")

folder_names <- file_and_folder_names %>%
  .[str_detect(., ".txt", negate = TRUE)]

lumped_stats_per_location_obemme <- map(
  folder_names,
  ~process_polygon(., "results_obemme", classified_cwp_data)
) %>% bind_rows()






lumped_stats_per_location <- bind_rows(lumped_stats_per_location_lower_emme, lumped_stats_per_location_upper_emme,lumped_stats_per_location_obemme)



write_csv(lumped_stats_per_location, "lumped_stats_emme.csv")



  
