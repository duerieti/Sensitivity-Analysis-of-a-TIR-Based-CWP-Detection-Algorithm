library(tidyverse)
library(tmap)
library(sf)


setwd("./compute_morris")



data <- read_delim("cwp_annotated.csv")

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
        identifier = row$Index
      )
    }

    else {
      cropped_list[[i]] <- tibble(
        total_area = 0,
        morris_row_index = unique(polys$morris_row_index),
        identifier = row$Index,
        Class = row$Class
      )
    
    }
  }

  return(cropped_list)
  
}




compute_lumped_statistic <- function(poly_region) {

  if ("sf" %in% class(poly_region)) {
    lumped_stats <- poly_region %>%
      sf::st_drop_geometry() %>%
      summarise(
        total_area        = sum(area_m2),
        morris_row_index  = unique(morris_row_index),
        identifier = unique(identifier),
        Class = unique(Class)
      )
  }

  else {
    lumped_stats <- poly_region
  }

  
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



file_and_folder_names <- list.files("./results_emmev1_big")

folder_names <- file_and_folder_names %>%
  .[str_detect(., ".txt", negate = TRUE)]


emme_v1_data <- classified_cwp_data %>% filter(Dataset == "emme_v1")

lumped_stats_per_location_emme_v1_big <- map(
  folder_names,
  ~process_polygon(., "results_emmev1_big", emme_v1_data)
) %>% bind_rows()


file_and_folder_names <- list.files("./results_emme_v2_big")

folder_names <- file_and_folder_names %>%
  .[str_detect(., ".txt", negate = TRUE)]


emme_v2_data <- classified_cwp_data %>% filter(Dataset == "emme_v2")

lumped_stats_per_location_emme_v2_big <- map(
  folder_names,
  ~process_polygon(., "results_emme_v2_big", emme_v2_data)
) %>% bind_rows()



lumped_stats_per_location <- bind_rows(lumped_stats_per_location_emme_v2_big, lumped_stats_per_location_emme_v1_big)



write_csv(lumped_stats_per_location, "lumped_stats_emme_big.csv")



  
