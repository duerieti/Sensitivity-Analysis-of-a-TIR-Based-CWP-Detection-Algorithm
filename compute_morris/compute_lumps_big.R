library(tidyverse)
library(tmap)
library(sf)

setwd("./compute_morris")

# ── 0. READ ANNOTATION DATA ───────────────────────────────────────────────────
# The annotation data contains bounding boxes of manually identified and
# classified cold water patches (Tributary / Non-Tributary). Each bounding box
# defines a region of interest within which the detected CWP polygons will be
# evaluated. The bounding box coordinates are stored as a single comma-separated
# string and need to be split into four separate columns.

classified_cwp_data <- read_delim("cwp_annotated.csv") %>%
  # split the bounding box coordinate string into four separate columns
  separate_wider_delim(cols = Bbox, delim = ",", names = c("xmin", "ymin", "xmax", "ymax")) %>%
  mutate(
    xmin = as.numeric(xmin),
    ymin = as.numeric(ymin),
    xmax = as.numeric(xmax),
    ymax = as.numeric(ymax)
  )


# ── HELPERS ───────────────────────────────────────────────────────────────────

# For a set of detected CWP polygons and a dataframe of annotated bounding boxes,
# crop the polygons to each bounding box. This isolates the detected CWPs at each
# manually annotated location. If no CWP was detected within a bounding box
# (which can happen depending on the parameter tuple, especially delta_C), a
# zero-area placeholder is returned. This because if a cwp is not present in a result, that
# actually coresponds to an area of for this cwp.

crop_to_regions <- function(polys, df) {

  n_row <- df %>% nrow()
  cropped_list <- vector("list", n_row)

  for (i in seq(1, n_row, 1)) {
    # get the i-th annotated bounding box
    row <- df[i, ]

    # construct an sf bounding box from the annotation coordinates.
    # CRS 2056 is the Swiss metric coordinate system (LV95).
    bbox <- st_bbox(
      c(xmin = row$xmin, ymin = row$ymin, xmax = row$xmax, ymax = row$ymax),
      crs = st_crs(2056)
    )

    # retain only the detected CWP polygons that fall within this bounding box
    cropped_region <- st_crop(polys, bbox)

    # if at least one CWP polygon was detected within the bounding box
    if (nrow(cropped_region) != 0) {
      # tag the detected polygons with the class (Tributary / Non-Tributary)
      # and the identifier of the annotated bounding box
      cropped_list[[i]] <- cropped_region %>%
        mutate(
          Class      = row$Class,
          identifier = row$Index
        )
    } else {
      # no CWP was detected at this location for this parameter tuple.
      # return a zero-area placeholder. This because if a cwp is not present in a result, that
      # actually coresponds to an area of for this cwp.

      cropped_list[[i]] <- tibble(
        total_area       = 0,
        morris_row_index = unique(polys$morris_row_index), # parameter tuple index
        identifier       = row$Index,
        Class            = row$Class
      )
    }
  }

  return(cropped_list)
}


# Reduce the detected CWP polygons within one annotated bounding box to a single
# scalar: the total detected area. This is the model output metric used to compute
# Morris elementary effects. If the input is already a zero-area placeholder tibble
# (no detection), it is passed through unchanged.
compute_lumped_statistic <- function(poly_region) {

  if ("sf" %in% class(poly_region)) {
    # sum the area of all detected CWP polygons within the bounding box and
    # carry the parameter tuple index, bounding box identifier and class
    lumped_stats <- poly_region %>%
      sf::st_drop_geometry() %>%
      summarise(
        total_area       = sum(area_m2),
        morris_row_index = unique(morris_row_index),
        identifier       = unique(identifier),
        Class            = unique(Class)
      )
  } else {
    # input is already a zero-area placeholder tibble — pass through unchanged
    lumped_stats <- poly_region
  }

  return(lumped_stats)
}


# For one results folder (one Morris parameter tuple), read the detected CWP
# polygons, crop them to the annotated bounding boxes, and compute the total
# detected area per bounding box. The folder name encodes the Morris row index the realisations
# of the parameters in this specific run (a.k.a the parameter tuple)
# which is used to link results back to the parameter tuple they were produced with.
process_polygon <- function(folder_name, base_path, df) {

  # extract the Morris row index from the folder name (format: results_X)
  index <- as.integer(str_split(folder_name, "_")[[1]][2])

  # build the path to the detected CWP shapefile for this parameter tuple
  read_path <- file.path(base_path, folder_name, "final_polys.shp")

  # read the detected CWP polygons, tag them with the Morris row index,
  # crop to each annotated bounding box, and drop empty list elements
  region_polys <- read_path %>%
    sf::read_sf() %>%
    mutate(morris_row_index = index) %>%
    crop_to_regions(df) %>%
    compact()

  # compute the total detected area per annotated bounding box
  lumped_stats_per_region <- map(region_polys, compute_lumped_statistic) %>%
    bind_rows()

  return(lumped_stats_per_region)
}


# ── 1. PROCESS RESULTS PER RIVER SECTION ─────────────────────────────────────
# For each river section, list all result folders (one per Morris parameter
# tuple), filter the annotation data to that section, and compute the total
# detected CWP area per annotated location per parameter tuple.

# ── emme v1 ───────────────────────────────────────────────────────────────────

# list all result folders for emme v1, excluding the joblist.txt file
folder_names <- list.files("../morris_screening_emme_v1_more_params/results") %>%
  .[str_detect(., ".txt", negate = TRUE)]

# filter annotation data to emme v1 bounding boxes
emme_v1_data <- classified_cwp_data %>% filter(Dataset == "emme_v1")

# process all parameter tuples for emme v1
lumped_stats_emme_v1 <- map(
  folder_names,
  ~ process_polygon(., "../morris_screening_emme_v1_more_params/results", emme_v1_data)
) %>% bind_rows()

# ── emme v2 ───────────────────────────────────────────────────────────────────

folder_names <- list.files("../morris_screening_emme_v2_more_params/results") %>%
  .[str_detect(., ".txt", negate = TRUE)]

emme_v2_data <- classified_cwp_data %>% filter(Dataset == "emme_v2")

lumped_stats_emme_v2 <- map(
  folder_names,
  ~ process_polygon(., "../morris_screening_emme_v2_more_params/results", emme_v2_data)
) %>% bind_rows()

# ── obere emme ────────────────────────────────────────────────────────────────

folder_names <- list.files("../morris_screening_obemme_more_params/results") %>%
  .[str_detect(., ".txt", negate = TRUE)]

obemme_data <- classified_cwp_data %>% filter(Dataset == "obemme")

lumped_stats_obemme <- map(
  folder_names,
  ~ process_polygon(., "../morris_screening_obemme_more_params/results", obemme_data)
) %>% bind_rows()


# ── 2. COMBINE AND WRITE ──────────────────────────────────────────────────────
# Combine results from all three river sections into one table. Each row
# represents one annotated CWP location under one Morris parameter tuple,
# with the total detected area as the scalar model output.

lumped_stats_all <- bind_rows(lumped_stats_emme_v1, lumped_stats_emme_v2, lumped_stats_obemme)
write_csv(lumped_stats_all, "lumped_stats_emme_big.csv")