library(tidyverse)
library(tmap)
library(sf)



# ── 0. READ ANNOTATION DATA ───────────────────────────────────────────────────
# (unchanged)

classified_cwp_data <- read_delim("compute_morris/cwp_annotated.csv") %>%
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
# zero-area / zero-deltaT placeholder is returned. This is because "no CWP
# detected" corresponds to both area = 0 and deltaT = 0 (no thermal anomaly).

crop_to_regions <- function(polys, df) {

  n_row <- df %>% nrow()
  cropped_list <- vector("list", n_row)

  for (i in seq(1, n_row, 1)) {
    row <- df[i, ]

    bbox <- st_bbox(
      c(xmin = row$xmin, ymin = row$ymin, xmax = row$xmax, ymax = row$ymax),
      crs = st_crs(2056)
    )

    cropped_region <- st_crop(polys, bbox)

    if (nrow(cropped_region) != 0) {
      cropped_list[[i]] <- cropped_region %>%
        mutate(
          Class      = row$Class,
          identifier = row$Index
        )
    } else {
      # no CWP was detected at this location for this parameter tuple.
      # area = 0 (no patch), deltaT = 0 (no thermal anomaly) — both are
      # the natural "neutral" values corresponding to "nothing detected".
      cropped_list[[i]] <- tibble(
        total_area       = 0,
        mean_deltaT      = 0,
        morris_row_index = unique(polys$morris_row_index), # parameter tuple index
        identifier       = row$Index,
        Class            = row$Class
      )
    }
  }

  return(cropped_list)
}


# Reduce the detected CWP polygons within one annotated bounding box to two
# scalars: the total detected area, and the area-weighted mean temperature
# delta (deltaT) across all detected polygons. These are the model output
# metrics used to compute Morris elementary effects. If the input is already
# a zero-placeholder tibble (no detection), it is passed through unchanged.
compute_lumped_statistic <- function(poly_region) {

  if ("sf" %in% class(poly_region)) {
    lumped_stats <- poly_region %>%
      sf::st_drop_geometry() %>%
      summarise(
        total_area       = sum(area_m2),
        # area-weighted mean deltaT: larger polygons contribute proportionally
        # more to the location's overall thermal-deficit characterization
        mean_deltaT      = sum(deltaT * area_m2) / sum(area_m2),
        morris_row_index = unique(morris_row_index),
        identifier       = unique(identifier),
        Class            = unique(Class)
      )
  } else {
    # input is already a zero-placeholder tibble — pass through unchanged
    lumped_stats <- poly_region
  }

  return(lumped_stats)
}


# (process_polygon unchanged — it just calls the two helpers above)
process_polygon <- function(folder_name, base_path, df) {

  index <- as.integer(str_split(folder_name, "_")[[1]][2])

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

# ── 1. PROCESS RESULTS PER RIVER SECTION ─────────────────────────────────────
# (unchanged — total_area AND mean_deltaT now flow through automatically)

folder_names <- list.files("./morris_screening_emme_v1_more_params/results_emmev1_big") %>%
  .[str_detect(., ".txt", negate = TRUE)]

emme_v1_data <- classified_cwp_data %>% filter(Dataset == "emme_v1")

lumped_stats_emme_v1 <- map(
  folder_names,
  ~ process_polygon(., "./morris_screening_emme_v1_more_params/results_emmev1_big", emme_v1_data)
) %>% bind_rows()

folder_names <- list.files("./morris_screening_emme_v2_more_params/results_emmev2_big") %>%
  .[str_detect(., ".txt", negate = TRUE)]

emme_v2_data <- classified_cwp_data %>% filter(Dataset == "emme_v2")

lumped_stats_emme_v2 <- map(
  folder_names,
  ~ process_polygon(., "./morris_screening_emme_v2_more_params/results_emmev2_big", emme_v2_data)
) %>% bind_rows()

folder_names <- list.files("./morris_screening_obemme_more_params/results_obemme_big") %>%
  .[str_detect(., ".txt", negate = TRUE)]

obemme_data <- classified_cwp_data %>% filter(Dataset == "obemme")

lumped_stats_obemme <- map(
  folder_names,
  ~ process_polygon(., "./morris_screening_obemme_more_params/results_obemme_big", obemme_data)
) %>% bind_rows()


# ── 2. COMBINE AND WRITE ──────────────────────────────────────────────────────
# Each row now carries both total_area and mean_deltaT per
# (annotated location × Morris parameter tuple).


lumped_stats_all <- bind_rows(lumped_stats_emme_v1, lumped_stats_emme_v2, lumped_stats_obemme)

lumped_stats_all

write_csv(lumped_stats_all, "lumped_stats_emme_big.csv")
