library(tidyverse)
library(tmap)
library(sf)

# set the interpreter path to this folder
setwd("./compute_morris")



# read the anotation data of the cold water patches 
# seperate the coordinate string into four different columns
classified_cwp_data <- read_delim("cwp_annotated.csv") %>%
  separate_wider_delim(cols = Bbox, delim = ",", names = c("xmin", "ymin", "xmax", "ymax")) %>%
  mutate(
    xmin = as.numeric(xmin),
    ymin = as.numeric(ymin),
    xmax = as.numeric(xmax),
    ymax = as.numeric(ymax)
  )



crop_to_regions <- function(polys, df) {

  # get the number of rows in the dataframe
  n_row <- df %>% nrow()

  # create a list of length n_row
  cropped_list <- vector("list", n_row)

  # or every row
  for (i in seq(1,n_row,1)) {
    row <- df[i,  ] # get the row

    # prepare a bounding box if the xmin, xmax, ymin, ymax columns in the row
    bbox <- st_bbox(
      c(
      xmin = row$xmin,
      ymin = row$ymin,
      xmax = row$xmax,
      ymax = row$ymax),
    crs = st_crs(2056) # swiss coordinate system
    )

    # crop the cwp simple feature collection to the bounding box
    # this has the effect that only the cwp lying within that bounding box is retained
    # and all other cwps are discarded

    croped_region <- st_crop(polys, bbox)
    
    # if the filtered sfc of the cwps is not empty
    if (nrow(croped_region) != 0) {
      # add the class of the cwp in the bounding box (Tributary / Non Tributary)
      # and the index of the cwp in the bounding box
      cropped_list[[i]] <- croped_region %>% 
        mutate(
        Class = row$Class,
        identifier = row$Index
      )
    }
    # else if the sfc is empty since the cwp is not even present in this realisation of the cwp detection algorithm
    # which can happen depending on the parameters beeing set, especially due to the delta_C for flagging beeing used
    else {
      # pass a tibble with a total area of zero (the cwp is not there at all, so the area is zero)
      cropped_list[[i]] <- tibble(
        total_area = 0,
        morris_row_index = unique(polys$morris_row_index), # get the index of the parameter tupple which the cwp realisation belongs to
        identifier = row$Index, # also add the index of the bounding box
        Class = row$Class # also add the class of the bounding box
      )
    
    }
  }

  return(cropped_list)
  
}




compute_lumped_statistic <- function(poly_region) {

  # if poly_region is a simple feature object
  # compute the total area of the cwp
  # add the morris row index of the parameter tupple corresponding to the
  # realisation of the cwps
  # ad the iddentifier of the bounding box
  # add the class of the bounding box (classification of the cwp)
  if ("sf" %in% class(poly_region)) {
    lumped_stats <- poly_region %>%
      sf::st_drop_geometry() %>% # drop geometry to go to a data.frame
      summarise(
        total_area        = sum(area_m2),
        morris_row_index  = unique(morris_row_index),
        identifier = unique(identifier),
        Class = unique(Class)
      )
  }

  # else if it is something else (in this case a tibble)
  else {
    lumped_stats <- poly_region # just assign that directly
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


file_and_folder_names <- list.files("./results_obemme_big")

folder_names <- file_and_folder_names %>%
  .[str_detect(., ".txt", negate = TRUE)]

obemme_data <- classified_cwp_data %>% filter(Dataset == "obemme")

lumped_stats_per_location_obemme_big <- map(
  folder_names,
  ~process_polygon(., "results_obemme_big", obemme_data )
) %>% bind_rows()




lumped_stats_per_location <- bind_rows(lumped_stats_per_location_emme_v2_big, lumped_stats_per_location_emme_v1_big, lumped_stats_per_location_obemme_big)



write_csv(lumped_stats_per_location, "lumped_stats_emme_big.csv")



  
