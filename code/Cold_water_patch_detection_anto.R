# --- Cold-water patch detection & sensitivity analysis -----------------------
# Dependencies
rm(list = ls())

library(terra)
library(sf)
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(stringr)
library(readr)
library(lwgeom)
library(qs)
library(grDevices)




options(sf_use_s2 = FALSE)

# ---- USER INPUTS ------------------------------------------------------------
#EMME_V1 (Emme UNterlauf)
setwd("D:/Daten_anto/Thermo_bern/automatisation/v1_emme")
ras_path  <- r"(s:\pools\n\N-IUNR-Allgemein\Zentren\Integrative Ökologie\FG Ökohydrologie\1. Projekte\3. Aktuell\Thermo_BE_2024_2025\8_thermal_analysis\thermal_rasters_FINAL\mean_v01emme.tif)"      # <-  raster with water pixels only
line_path <-  r"(s:\pools\n\N-IUNR-Allgemein\Zentren\Integrative Ökologie\FG Ökohydrologie\1. Projekte\3. Aktuell\Thermo_BE_2024_2025\8_thermal_analysis\Centerlines_FINAL\Emme_V01.shp)"      # <- centerline
out_dir      <- "output/cold_patches_out_linear_buffer_emmev1_01"


#urtene_V1
setwd("D:/Daten_anto/Thermo_bern/automatisation/v1_urtene/")
ras_path  <- r"(s:\pools\n\N-IUNR-Allgemein\Zentren\Integrative Ökologie\FG Ökohydrologie\1. Projekte\3. Aktuell\Thermo_BE_2024_2025\8_thermal_analysis\thermal_rasters_FINAL\mean_v1urtene.tif)"      # <- raster with water pixels only
line_path <-  r"(s:\pools\n\N-IUNR-Allgemein\Zentren\Integrative Ökologie\FG Ökohydrologie\1. Projekte\3. Aktuell\Thermo_BE_2024_2025\8_thermal_analysis\Centerlines_FINAL\Urtene_V01.shp)"      # <- centerline
out_dir      <- "output/cold_patches_out_linear_buffer_urtenev1_01"


#grüene_V2
setwd("D:/Daten_anto/Thermo_bern/automatisation/v2_emmgrue/")
ras_path  <- r"(s:\pools\n\N-IUNR-Allgemein\Zentren\Integrative Ökologie\FG Ökohydrologie\1. Projekte\3. Aktuell\Thermo_BE_2024_2025\8_thermal_analysis\thermal_rasters_FINAL\mean_v2gruen.tif)"      # <-  raster with water pixels only
line_path <-  r"(s:\pools\n\N-IUNR-Allgemein\Zentren\Integrative Ökologie\FG Ökohydrologie\1. Projekte\3. Aktuell\Thermo_BE_2024_2025\8_thermal_analysis\Centerlines_FINAL\Gruene_V02.shp)"      # <-  centerline
out_dir      <- "output/cold_patches_out_linear_buffer_gruene_01"


#emme_V2  (Emme Oberlauf 2)
setwd("D:/Daten_anto/Thermo_bern/automatisation/v2_emmgrue/")
ras_path  <- r"(s:\pools\n\N-IUNR-Allgemein\Zentren\Integrative Ökologie\FG Ökohydrologie\1. Projekte\3. Aktuell\Thermo_BE_2024_2025\8_thermal_analysis\thermal_rasters_FINAL\mean_v2emme.tif)"      # <- your raster with water pixels only
line_path <-  r"(s:\pools\n\N-IUNR-Allgemein\Zentren\Integrative Ökologie\FG Ökohydrologie\1. Projekte\3. Aktuell\Thermo_BE_2024_2025\8_thermal_analysis\Centerlines_FINAL\Emme_V02.shp)"      # <- your centerline
out_dir      <- "output/cold_patches_out_linear_buffer_emmev2_01"

#ilf_V2
setwd("D:/Daten_anto/Thermo_bern/automatisation/v2_ilf/")
ras_path  <- r"(s:\pools\n\N-IUNR-Allgemein\Zentren\Integrative Ökologie\FG Ökohydrologie\1. Projekte\3. Aktuell\Thermo_BE_2024_2025\8_thermal_analysis\thermal_rasters_FINAL\mean_v02ilfis.tif)"      # <- your raster with water pixels only
line_path <-  r"(s:\pools\n\N-IUNR-Allgemein\Zentren\Integrative Ökologie\FG Ökohydrologie\1. Projekte\3. Aktuell\Thermo_BE_2024_2025\8_thermal_analysis\Centerlines_FINAL\Ilfis_V02.shp)"      # <- your centerline
out_dir      <- "output/cold_patches_out_linear_buffer_ilfv2_01"


#obereemme_V2 (Emme Oberlauf 1)
setwd("D:/Daten_anto/Thermo_bern/automatisation/v2_obereemme/")
ras_path  <- r"(s:\pools\n\N-IUNR-Allgemein\Zentren\Integrative Ökologie\FG Ökohydrologie\1. Projekte\3. Aktuell\Thermo_BE_2024_2025\8_thermal_analysis\thermal_rasters_FINAL\mean_v2obemme.tif)"      # <- your raster with water pixels only
line_path <-  r"(s:\pools\n\N-IUNR-Allgemein\Zentren\Integrative Ökologie\FG Ökohydrologie\1. Projekte\3. Aktuell\Thermo_BE_2024_2025\8_thermal_analysis\Centerlines_FINAL\Obere_Emme_V02.shp)"      # <- your centerline
out_dir      <- "output/cold_patches_out_linear_buffer_v2_obereemme_01"




steps_m      <- c(100, 200, 500, 1000, 2000, 5000)                 # sampling step along centerline (m)
buffers_px   <- c(3)                          #  buffer radius (pixels)
deltas_C     <- c(0.5, 1.0, 2.0, 3.0)                   # ΔT thresholds [°C]
min_area_m2  <- 2                                  # min patch area
handle_small <- "drop"                             # "drop" OR "merge"


# ---- CORE FUNCTIONS ----------------------------------------------------------


# --- helper for pretty, timestamped logs ------------------------------------
# --- tiny logger --------------------------------------------------------------
log_msg <- function(fmt, ...) message(sprintf("[%s] %s",
                                              format(Sys.time(), "%H:%M:%S"), sprintf(fmt, ...)))


detect_cold_patches <- function(
    ras_path, line_path,
    steps = c(100, 200, 500, 1000),
    buffers_px = c(1,3,5,10),
    deltas = c(0.5,1,2,3),
    min_patch_area_m2 = 2,
    handle_small = c("drop","merge"),
    round_to = 0.1,                 # 0 or NULL disables rounding
    slab_halfwidth_m = 60,          # <- default lateral half-width (meters)
    out_dir = "output/cold_patches_out_linear_buffer",
    connect_diagonals = TRUE,       # 8-neighbor connectivity
    union_chunk_size = 50000L       # rectangles per chunk before partial dissolve
) {
  # logger
  log_msg <- function(fmt, ...) message(sprintf("[%s] %s",
                                                format(Sys.time(), "%H:%M:%S"), sprintf(fmt, ...)))
  
  options(sf_use_s2 = FALSE)
  handle_small <- match.arg(handle_small)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  shp_dir <- file.path(out_dir, "shapefiles"); dir.create(shp_dir, showWarnings = FALSE)
  
  # --- read inputs
  log_msg("Reading raster: %s", ras_path)
  r <- terra::rast(ras_path)
  stopifnot(terra::nlyr(r) == 1)
  names(r) <- "T"
  if (terra::is.lonlat(r)) stop("Raster must be in a metric CRS.")
  rres  <- terra::res(r); cell_m <- mean(rres)
  
  log_msg("Reading centerline: %s", line_path)
  line <- sf::st_read(line_path, quiet = TRUE) |> sf::st_zm(TRUE, "ZM")
  if (!is.na(sf::st_crs(line)) && !identical(terra::crs(r), sf::st_crs(line)$wkt)) {
    log_msg("Reprojecting centerline to raster CRS"); line <- sf::st_transform(line, terra::crs(r))
  }
  g <- sf::st_union(line)
  if (inherits(g, "sfc_GEOMETRYCOLLECTION")) g <- sf::st_collection_extract(g, "LINESTRING")
  if (inherits(g, "sfc_MULTILINESTRING")) {
    g <- sf::st_line_merge(g)
    if (inherits(g, "sfc_GEOMETRYCOLLECTION")) g <- sf::st_collection_extract(g, "LINESTRING")
  }
  if (!inherits(g, "sfc_LINESTRING")) stop("Centerline must resolve to a LINESTRING.")
  line <- sf::st_sf(geometry = g)
  
  rfactor <- if (!is.null(round_to) && round_to > 0) 1 / round_to else NA_real_
  
  total_combos <- length(steps) * length(buffers_px) * length(deltas)
  combo_i <- 0L
  log_msg("Grid: steps=%s | buffers_px=%s | deltas=%s | combos=%d",
          paste(steps, collapse=","), paste(buffers_px, collapse=","), paste(deltas, collapse=","), total_combos)
  
  all_summaries <- list(); all_areas <- list(); written_shps <- character(0)
  
  # --- helper: vector-only polygonization → single-part polygons
  polygonize_flagged_vector <- function(flagged_cells, r, chunk_size = 50000L,
                                        connect_diag = TRUE) {
    if (!length(flagged_cells)) {
      return(sf::st_sf(geometry = sf::st_sfc(), crs = sf::st_crs(terra::crs(r))))
    }
    res    <- terra::res(r)
    crs_wk <- terra::crs(r)
    
    make_rects_sfc <- function(xc, yc, rx, ry, crs_wk) {
      polys <- vector("list", length(xc))
      for (i in seq_along(xc)) {
        x <- xc[i]; y <- yc[i]
        m <- matrix(c(
          x - rx/2, y - ry/2,
          x + rx/2, y - ry/2,
          x + rx/2, y + ry/2,
          x - rx/2, y + ry/2,
          x - rx/2, y - ry/2
        ), ncol = 2, byrow = TRUE)
        polys[[i]] <- sf::st_polygon(list(m))
      }
      sf::st_sfc(polys, crs = sf::st_crs(crs_wk))
    }
    
    n       <- length(flagged_cells)
    chunks  <- split(seq_len(n), ceiling(seq_len(n) / chunk_size))
    parts_sv <- vector("list", length(chunks))  # dissolved per chunk
    
    # 1) dissolve per chunk with terra::aggregate()
    for (ci in seq_along(chunks)) {
      idx <- chunks[[ci]]
      xy  <- terra::xyFromCell(r, flagged_cells[idx])
      sfc <- make_rects_sfc(xy[,1], xy[,2], res[1], res[2], crs_wk)
      
      v <- terra::vect(sf::st_sf(geometry = sfc))
      v$grp <- 1L
      vd <- terra::aggregate(v, by = "grp")
      parts_sv[[ci]] <- vd
      message(sprintf("     chunk %d/%d: %d cells → %d multipart piece(s)",
                      ci, length(chunks), length(idx), nrow(vd)))
    }
    
    # 2) global dissolve (sf::st_union), with optional diagonal connection
    parts_sf <- lapply(parts_sv, sf::st_as_sf)
    stitched <- do.call(rbind, parts_sf)
    geom     <- sf::st_geometry(stitched)
    geom     <- try(sf::st_make_valid(geom), silent = TRUE)
    eps      <- if (connect_diag) min(res) * 0.1 else 0
    
    if (eps > 0) {
      merged <- sf::st_union(sf::st_buffer(geom, eps))
      merged <- sf::st_buffer(merged, -eps)
    } else {
      merged <- sf::st_union(geom)
    }
    
    # 3) explode to single-parts: ensure **POLYGON** geometry (no multiparts)
    polys <- sf::st_collection_extract(merged, "POLYGON", warn = FALSE)
    polys <- suppressWarnings(sf::st_cast(polys, "POLYGON"))   # guarantees sfc_POLYGON
    sf::st_sf(geometry = polys)
  }
  
  # --- main loops
  for (S in steps) {
    L <- as.numeric(sf::st_length(line))
    if (L < 2*S) { log_msg("Skip step=%d (line too short)", S); next }
    
    pos <- sort(unique(c(0, seq(0, 1, by = S/L), 1)))
    pts <- sf::st_line_sample(line, sample = pos) |> sf::st_cast("POINT")
    pts_sf <- sf::st_as_sf(pts)
    tvec <- pos; npt <- length(tvec)
    if (npt < 2) { log_msg("Skip step=%d (stations=%d)", S, npt); next }
    log_msg("Step %d m: stations=%d", S, npt)
    
    # slabs: i..i+1 (length ≈ step), lateral half-width = slab_halfwidth_m
    from <- tvec[1:(npt-1)]; to <- tvec[2:npt]; k <- length(from)
    g1 <- sf::st_geometry(line)[1]
    log_msg("Building %d slabs (length≈%d m, halfwidth=%.2f m)", k, S, slab_halfwidth_m)
    
    subs_list <- vector("list", k)
    for (j in seq_len(k)) subs_list[[j]] <- lwgeom::st_linesubstring(g1, from[j], to[j])
    subs  <- do.call(c, subs_list)
    slabs <- sf::st_buffer(subs, dist = slab_halfwidth_m, endCapStyle="FLAT", joinStyle="MITRE", mitreLimit=2)
    slabs_sf <- sf::st_sf(geometry = slabs)
    
    # pixels per slab 
    log_msg("Extracting slab pixels (full raster) …")
    slab_df <- terra::extract(r, terra::vect(slabs_sf), cells = TRUE)
    slab_df <- slab_df[!is.na(slab_df$T), , drop = FALSE]
    if (!nrow(slab_df)) { log_msg("No slab pixels at step=%d → skipping this step", S); next }
    if (!is.na(rfactor)) slab_df$T <- round(slab_df$T * rfactor) / rfactor
    by_slab <- split(slab_df[, c("ID","cell","T")], slab_df$ID)
    log_msg("Slab rows=%d (unique cells=%d)", nrow(slab_df), length(unique(slab_df$cell)))
    
    # medians from buffers at station i (left endpoint of slab)
    # medians from a LINEAR CORRIDOR buffer around each slab (not a point buffer)
    for (Rpx in buffers_px) {
      Rm <- Rpx * cell_m
      log_msg("  Linear buffer %d px (~%.2f m): medians per %d slabs …", Rpx, Rm, k)
      
      median_round <- function(x) {
        x <- x[!is.na(x)]
        if (!length(x)) return(NA_real_)
        if (!is.na(rfactor)) x <- round(x * rfactor) / rfactor
        stats::median(x)
      }
      
      # compute Tmed per slab using a corridor buffer around the line substring
      Tmed <- numeric(k)
      for (j in seq_len(k)) {
        # corridor around slab j (half-width = Rm)
        # end caps flat so corridor length matches the slab; joinStyle can be MITRE/ROUND
        corridor_j <- sf::st_buffer(
          sf::st_sfc(subs[j], crs = sf::st_crs(line)),
          dist = Rm, endCapStyle = "FLAT", joinStyle = "MITRE", mitreLimit = 2
        )
        # extract raster values inside corridor
        vals_df <- terra::extract(r, terra::vect(sf::st_sf(geometry = corridor_j)))
        v <- vals_df[[ names(r)[1] ]]
        # median with optional rounding
        Tmed[j] <- median_round(v)
      }
      
      log_msg("  Medians ready (NA=%d)", sum(is.na(Tmed)))
      
      for (dC in deltas) {
        combo_i  <- combo_i + 1L
        combo_tag <- sprintf("S%sm_B%spx_D%0.1fC", S, Rpx, dC)
        log_msg("   [%d/%d] ΔT=%.1f°C → %s", combo_i, total_combos, dC, combo_tag)
        
        # flag slab pixels using corridor median at slab j
        flagged <- integer(0L)
        for (j in seq_along(by_slab)) {
          tm <- Tmed[j]; if (is.na(tm)) next
          dfj <- by_slab[[j]]
          keep <- dfj$cell[(tm - dfj$T) >= dC]
          if (length(keep)) flagged <- c(flagged, keep)
        }
        flagged <- unique(flagged)
        log_msg("   Flagged cells: %d", length(flagged))
        
        
        # polygonize → global dissolve → explode to single parts
        if (!length(flagged)) {
          pol_sf <- sf::st_sf(geometry = sf::st_sfc(), crs = sf::st_crs(line))
          log_msg("   No polygons for %s", combo_tag)
        } else {
          pol_sf <- polygonize_flagged_vector(flagged, r,
                                              chunk_size = union_chunk_size,
                                              connect_diag = connect_diagonals)
          log_msg("   Patches (single-part) created: %d", nrow(pol_sf))
        }
        
        # compute areas AFTER single-part explode
        if (inherits(pol_sf, "sfc")) pol_sf <- sf::st_sf(geometry = pol_sf, crs = sf::st_crs(line))
        if (nrow(pol_sf) > 0) {
          pol_sf$area_m2 <- tryCatch(
            as.numeric(sf::st_area(pol_sf)),
            error = function(e) as.numeric(terra::expanse(terra::vect(pol_sf), unit = "m", transform = FALSE))
          )
          
          # drop/merge small AFTER dissolve+explode
          if (handle_small == "drop") {
            pol_sf <- dplyr::filter(pol_sf, area_m2 >= min_patch_area_m2)
            log_msg("   Dropped patches < %.2f m²; kept: %d", min_patch_area_m2, nrow(pol_sf))
          } else {
            big   <- dplyr::filter(pol_sf, area_m2 >= min_patch_area_m2)
            small <- dplyr::filter(pol_sf, area_m2 <  min_patch_area_m2)
            if (nrow(small) > 0 && nrow(big) == 0) {
              pol_sf <- dplyr::slice_max(small, area_m2, n = 1)
              log_msg("   No big patches → kept only largest small patch")
            } else if (nrow(small) > 0 && nrow(big) > 0) {
              nearest_idx <- sf::st_nearest_feature(small, big)
              small$merge_id <- nearest_idx
              big$merge_id   <- seq_len(nrow(big))
              pol_sf <- dplyr::bind_rows(big, small) |>
                dplyr::group_by(merge_id) |>
                dplyr::summarise(area_m2 = sum(area_m2), .groups = "drop")
              log_msg("   Merged %d small patches → now %d patches", nrow(small), nrow(pol_sf))
            } else {
              pol_sf <- big
            }
          }
        }
        
        
        # --- Add temperature statistics for each polygon ---
        if (nrow(pol_sf) > 0) {
          log_msg("   Calculating temperature statistics for %d poligons …", nrow(pol_sf))
          
          # convert raster in 'SpatRaster' and poligons in 'SpatVector'
          pol_v <- terra::vect(pol_sf)
          
          # Extract the temperature values within each polygon
          vals_list <- terra::extract(r, pol_v)
          vals_list <- split(vals_list[[names(r)[1]]], vals_list$ID)
          
          # helper funtion for calculating statistics
          get_stats <- function(x) {
            x <- x[!is.na(x)]
            if (!length(x)) return(rep(NA_real_, 4))
            c(min(x), max(x), stats::median(x), mean(x))
          }
          
          stats_mat <- t(vapply(vals_list, get_stats, numeric(4)))
          colnames(stats_mat) <- c("T_min", "T_max", "T_med", "T_mean")
          
          pol_sf$T_min  <- stats_mat[, "T_min"]
          pol_sf$T_max  <- stats_mat[, "T_max"]
          pol_sf$T_med  <- stats_mat[, "T_med"]
          pol_sf$T_mean <- stats_mat[, "T_mean"]
          
          # Assign to each polygon the median of the corresponding slab
          # (computed from Tmed[j] of the selected pixels)
          
          # Identify the pixels contained within each polygon
          cell_extract <- terra::extract(r, pol_v, cells = TRUE)
          cell_extract <- cell_extract[!is.na(cell_extract$cell), ]
          
          
          # Map each cell to the slab it originates from (using by_slab)
          cell_to_slab_list <- lapply(seq_along(by_slab), function(j) {
            df <- by_slab[[j]]
            if (nrow(df) == 0) return(NULL)        # salta slab vuoti
            data.frame(cell = df$cell, slab_id = j)
          })
          cell_to_slab <- do.call(rbind, cell_to_slab_list)
          
          
          merged <- merge(cell_extract, cell_to_slab, by = "cell", all.x = TRUE)
          
          # For each polygon, compute the median of the corresponding Tmed[slab_id]
          slab_meds <- tapply(merged$slab_id, merged$ID, function(ids) {
            if (all(is.na(ids))) return(NA_real_)
            stats::median(Tmed[ids], na.rm = TRUE)
          })
          
          pol_sf$Tmed_slab <- as.numeric(slab_meds[as.character(seq_len(nrow(pol_sf)))])
          
          # deltaT = difference between the polygon T_med and Tmed_slab
          pol_sf$deltaT <- pol_sf$Tmed_slab - pol_sf$T_med 
          
          pol_sf$delta_thresh <- dC
          pol_sf <- dplyr::filter(pol_sf, deltaT >= dC)
          
          
          log_msg("   Temperature statistics added to the features")
        }
        
        
        
        # write one shapefile per combo (skip empty)
        if (nrow(pol_sf) > 0) {
          shp_path <- file.path(shp_dir, paste0(combo_tag, ".shp"))
          if (file.exists(shp_path)) unlink(sub("\\.shp$", ".*", shp_path)) # clear old set
          ok <- try(sf::st_write(pol_sf, shp_path, driver = "ESRI Shapefile", quiet = TRUE), silent = TRUE)
          if (inherits(ok, "try-error")) {
            log_msg("   WARN: failed to write shapefile %s", shp_path)
          } else {
            log_msg("   Wrote shapefile: %s", shp_path)
            written_shps <- c(written_shps, shp_path)
          }
        } else {
          log_msg("   Skipped shapefile write (no polygons)")
        }
        
        # summaries
        if (nrow(pol_sf) > 0) {
          all_summaries[[length(all_summaries)+1]] <- tibble::tibble(
            step_m = S, buffer_px = Rpx, delta = dC,
            n_patches = nrow(pol_sf),
            total_area_m2 = sum(pol_sf$area_m2),
            median_area_m2 = stats::median(pol_sf$area_m2),
            mean_area_m2 = mean(pol_sf$area_m2)
          )
          all_areas[[length(all_areas)+1]] <-
            dplyr::transmute(pol_sf, step_m = S, buffer_px = Rpx, delta = dC, area_m2 = area_m2)
        } else {
          all_summaries[[length(all_summaries)+1]] <- tibble::tibble(
            step_m = S, buffer_px = Rpx, delta = dC,
            n_patches = 0, total_area_m2 = 0,
            median_area_m2 = NA_real_, mean_area_m2 = NA_real_
          )
        }
        
      } # deltas
    }   # buffers
    
    log_msg("Finished step %d m", S)
  }     # steps
  
  # --- outputs: CSVs + plots (TIFF + PDF)
  summary_tbl <- dplyr::bind_rows(all_summaries)
  areas_tbl   <- dplyr::bind_rows(all_areas)
  
  log_msg("Writing summaries to %s", out_dir)
  readr::write_csv(summary_tbl, file.path(out_dir, "summary_patch_counts_and_area.csv"))
  if (nrow(areas_tbl)) readr::write_csv(areas_tbl, file.path(out_dir, "summary_patch_areas_long.csv"))
  

  
  log_msg("Done. Shapefiles in: %s", shp_dir)
  invisible(list(summary = summary_tbl, areas = areas_tbl, shp_dir = shp_dir,
                 shapefiles = written_shps, out_dir = out_dir))
}





# ---- RUN --------------------------------------------------------------------
res <- detect_cold_patches(
  ras_path       = ras_path,
  line_path      = line_path,
  steps          = steps_m,
  buffers_px     = buffers_px,
  deltas         = deltas_C,
  min_patch_area_m2 = min_area_m2,
  handle_small   = handle_small,
  out_dir        = out_dir
)



if (F) {
  #test with selected parameter combinations ----
  steps_m      <- c(100)                 # sampling step along centerline (m)
  buffers_px   <- c(10)                    # circular buffer radius (pixels)
  deltas_C     <- c( 2, 3)                   # ΔT thresholds [°C]
  min_area_m2  <- 2                                  # min patch area
  handle_small <- "drop"                             # "drop" OR "merge"
  out_dir      <- "output/cold_patches_out"
}



#Save and export the results -----
qsave(res, file=paste(out_dir,"/Results_CWP_detection.qs", sep=""))









if (F){
  pr=982
  
  # 2. Import results-------
  #out_dir      <- "output/cold_patches_out"
  res <- qread(file=paste(out_dir,"/Results_CWP_detection.qs", sep=""))
}
# 3. Plot the results ---------



# ==== Sensitivity plots driven by `res` =======================================
# Uses: res$summary (counts/areas by combo), res$areas (per-patch areas)
# Output: <res$out_dir>/plots_sensitivity/*.tiff|*.pdf


plot_dir <- file.path(out_dir, "plots_sensitivity")
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

sum0 <- res$summary %>% 
  dplyr::filter(!is.na(n_patches)) %>%
  dplyr::mutate(
    step_m    = as.integer(step_m),
    buffer_px = as.integer(buffer_px),
    delta     = as.numeric(delta),
    step_f    = factor(step_m),
    buf_f     = factor(buffer_px),
    dlt_f     = factor(delta)
  )

# Helper saver (TIFF + PDF)

save_tiff_pdf <- function(plot, path_no_ext, width_in, height_in, dpi = 300) {
  # TIFF
  grDevices::tiff(
    filename = paste0(path_no_ext, ".tiff"),
    width = width_in, height = height_in, units = "in",
    res = dpi, compression = "lzw"
  )
  print(plot)
  grDevices::dev.off()
  
  # PDF
  grDevices::pdf(
    file = paste0(path_no_ext, ".pdf"),
    width = width_in, height = height_in,
    useDingbats = FALSE
  )
  print(plot)
  grDevices::dev.off()
}


#grDevices::tiff(p_ecdf) and dev.off().

# 
# 1) Variance-explained bars (ANOVA main effects): which parameter explains more?
#    (balanced grid → Type I SS ~ main-effect contribution)
#



sum0 <- res$summary

# --- Varianzanteile für drei Zielgrössen --------------------------------------
metrics <- c("n_patches", "total_area_m2", "mean_area_m2","median_area_m2" )



# Detect study area automatically
# (expects path like ".../automatisation/<study_area>/...")

# Studiengebiet automatisch erkennen (aus Pfad: ".../automatisation/<gebiet>/...")

wd <- getwd()
parts <- strsplit(wd, .Platform$file.sep)[[1]]
auto_idx <- which(parts == "automatisation")

if (length(auto_idx) == 1 && length(parts) > auto_idx) {
  study_area <- parts[auto_idx + 1]
} else {
  study_area <- "unbekannt"
}

study_area_file  <- tolower(study_area)
study_area_title <- paste0(toupper(substring(study_area, 1, 1)), substring(study_area, 2))
message("Erkanntes Studiengebiet: ", study_area_title)


# Sichere Speicherfunktion (für ggplot / patchwork / gridExtra)

save_plot_pdf_tiff <- function(plot_obj, file_base, width = 10, height = 4, dpi = 300) {
  pdf_file  <- paste0(file_base, ".pdf")
  tiff_file <- paste0(file_base, ".tiff")
  dir.create(dirname(pdf_file), recursive = TRUE, showWarnings = FALSE)
  
  .draw_plot <- function(p) {
    if (inherits(p, "ggplot") || inherits(p, "patchwork")) {
      print(p)
    } else if (inherits(p, "grob") || inherits(p, "gTree")) {
      grid::grid.draw(p)
    } else if (is.list(p) && all(c("p1", "p2") %in% names(p))) {
      print(p$p1); print(p$p2)
    } else {
      try(print(p), silent = TRUE)
    }
  }
  
  grDevices::pdf(pdf_file, width = width, height = height)
  .draw_plot(plot_obj)
  grDevices::dev.off()
  
  grDevices::tiff(tiff_file, width = width, height = height,
                  units = "in", res = dpi, compression = "lzw")
  .draw_plot(plot_obj)
  grDevices::dev.off()
  
  message("✅ Gespeichert: ", basename(pdf_file), " und ", basename(tiff_file))
}

# Variance-explained bars (ANOVA main effects)
# Automatically skip factors with only one unique value



# Parameter-Check: nur Faktoren mit mehr als einem Wert

params_all <- c("step_m", "buffer_px", "delta")
params_vary <- params_all[sapply(sum0[params_all], function(x) length(unique(x)) > 1)]

message("Verwendete Parameter: ", paste(params_vary, collapse = ", "))


# ANOVA pro Kennzahl (inkl. mittlere Fläche)

fx_list <- lapply(metrics, function(m) {
  f <- as.formula(paste(m, "~", paste(sprintf("factor(%s)", params_vary), collapse = " + ")))
  fit <- lm(f, data = sum0)
  ss  <- anova(fit)
  den <- sum(ss[seq_along(params_vary), "Sum Sq"])
  
  data.frame(
    param    = params_vary,
    var_expl = if (den > 0) ss[seq_along(params_vary), "Sum Sq"] / den else NA_real_,
    metric   = m,
    row.names = NULL
  )
})

# Bisherige Metriken um mittlere Fläche erweitern
fx <- dplyr::bind_rows(fx_list) |>
  dplyr::mutate(
    metric = dplyr::recode(metric,
                           n_patches       = "Anzahl der CWP",
                           total_area_m2   = "Gesamtfläche",
                           mean_area_m2    = "Mittlere Fläche",
                           median_area_m2  = "Medianfläche"
    )
  )


# Plot

p_fx <- ggplot2::ggplot(
  fx, ggplot2::aes(x = reorder(param, -var_expl), y = var_expl)
) +
  ggplot2::geom_col(width = 0.7) +
  ggplot2::scale_y_continuous(labels = scales::percent) +
  ggplot2::facet_wrap(~ metric, nrow = 1, scales = "free_y") +
  ggplot2::labs(
    x = "Parameter",
    y = "% erklärte Varianz (Haupteffekte)",
    title = ""
  ) +
  ggplot2::theme_minimal(base_size = 12)


# Speichern

out_dir_plots <- file.path(out_dir, "plots_sensitivity")
dir.create(out_dir_plots, showWarnings = FALSE, recursive = TRUE)

out_file <- file.path(
  out_dir_plots,
  paste0("sensitivity_variance_explained_all_metrics_",pr, study_area_file)
)

save_tiff_pdf(p_fx, out_file, 15, 4)
pr=pr+1










#4. definitives plot ----




# Optionale Kombination nebeneinander

.combine_side_by_side <- function(..., title = NULL) {
  plots <- list(...)
  if (requireNamespace("patchwork", quietly = TRUE)) {
    comb <- Reduce(`+`, plots) + patchwork::plot_layout(ncol = length(plots), guides = "collect")
    comb <- comb & ggplot2::theme(legend.position = "bottom")
    if (!is.null(title)) comb <- comb + patchwork::plot_annotation(title = title)
    return(comb)
  } else if (requireNamespace("gridExtra", quietly = TRUE)) {
    grob <- do.call(gridExtra::arrangeGrob, c(plots, list(ncol = length(plots))))
    grid::grid.newpage(); grid::grid.draw(grob)
    return(grob)
  } else {
    return(plots)
  }
}


# Datenvorbereitung

steps_levels  <- sort(unique(sum0$step_m))
deltas_levels <- sort(unique(sum0$delta))
has_buffer <- "buffer_px" %in% names(sum0) && length(unique(sum0$buffer_px)) > 1

if (has_buffer) {
  buffers_levels <- sort(unique(sum0$buffer_px))
  sum0 <- sum0 %>%
    dplyr::mutate(
      step_f   = factor(step_m,    levels = steps_levels),
      buffer_f = factor(buffer_px, levels = buffers_levels)
    )
  message("Buffer erkannt mit ", length(buffers_levels), " Stufen.")
} else {
  sum0 <- sum0 %>%
    dplyr::mutate(step_f = factor(step_m, levels = steps_levels))
  message("⚠️ Kein variierender Buffer – Buffer aus der Darstellung entfernt.")
}


# Plot 1: Anzahl der CWP vs ΔT

if (has_buffer) {
  p_count <- ggplot2::ggplot(sum0, ggplot2::aes(
    x = delta, y = n_patches,
    color = step_f, linetype = buffer_f,
    group = interaction(step_f, buffer_f)
  ))
} else {
  p_count <- ggplot2::ggplot(sum0, ggplot2::aes(
    x = delta, y = n_patches,
    color = step_f, group = step_f
  ))
}

p_count <- p_count +
  ggplot2::geom_line(size = 0.9, na.rm = TRUE) +
  ggplot2::geom_point(size = 1.2, na.rm = TRUE) +
  ggplot2::scale_x_continuous(breaks = deltas_levels) +
  ggplot2::labs(
    x = expression(Delta*T~"[°C]"),
    y = "Anzahl der CWP",
    color = "Step [m]",
    linetype = if (has_buffer) "Puffer (Pixel)" else NULL,
    title = " "
  ) +
  ggplot2::theme_minimal(base_size = 12)


# Plot 2: Gesamtfläche vs ΔT

sum_tot <- sum0 %>% dplyr::filter(!is.na(total_area_m2))

if (has_buffer) {
  p_total_area <- ggplot2::ggplot(sum_tot, ggplot2::aes(
    x = delta, y = total_area_m2,
    color = step_f, linetype = buffer_f,
    group = interaction(step_f, buffer_f)
  ))
} else {
  p_total_area <- ggplot2::ggplot(sum_tot, ggplot2::aes(
    x = delta, y = total_area_m2,
    color = step_f, group = step_f
  ))
}

p_total_area <- p_total_area +
  ggplot2::geom_line(size = 0.9, na.rm = TRUE) +
  ggplot2::geom_point(size = 1.2, na.rm = TRUE) +
  ggplot2::scale_x_continuous(breaks = deltas_levels) +
  ggplot2::labs(
    x = expression(Delta*T~"[°C]"),
    y = "Gesamtfläche der CWP [m²]",
    color = "Step [m]",
    linetype = if (has_buffer) "Puffer (Pixel)" else NULL,
    title = " "
  ) +
  ggplot2::theme_minimal(base_size = 12)


# Plot 3: Mittlere Fläche vs ΔT

sum_mean <- sum0 %>% dplyr::filter(!is.na(mean_area_m2))

if (has_buffer) {
  p_mean_area <- ggplot2::ggplot(sum_mean, ggplot2::aes(
    x = delta, y = mean_area_m2,
    color = step_f, linetype = buffer_f,
    group = interaction(step_f, buffer_f)
  ))
} else {
  p_mean_area <- ggplot2::ggplot(sum_mean, ggplot2::aes(
    x = delta, y = mean_area_m2,
    color = step_f, group = step_f
  ))
}

p_mean_area <- p_mean_area +
  ggplot2::geom_line(size = 0.9, na.rm = TRUE) +
  ggplot2::geom_point(size = 1.2, na.rm = TRUE) +
  ggplot2::scale_x_continuous(breaks = deltas_levels) +
  ggplot2::labs(
    x = expression(Delta*T~"[°C]"),
    y = "Mittlere CWP-Flächengrösse [m²]",
    color = "Step [m]",
    linetype = if (has_buffer) "Puffer (Pixel)" else NULL,
    title = " "
  ) +
  ggplot2::theme_minimal(base_size = 12)


# Plot 4: Medianfläche vs ΔT

sum_med <- sum0 %>% dplyr::filter(!is.na(median_area_m2))

if (has_buffer) {
  p_med_area <- ggplot2::ggplot(sum_med, ggplot2::aes(
    x = delta, y = median_area_m2,
    color = step_f, linetype = buffer_f,
    group = interaction(step_f, buffer_f)
  ))
} else {
  p_med_area <- ggplot2::ggplot(sum_med, ggplot2::aes(
    x = delta, y = median_area_m2,
    color = step_f, group = step_f
  ))
}

p_med_area <- p_med_area +
  ggplot2::geom_line(size = 0.9, na.rm = TRUE) +
  ggplot2::geom_point(size = 1.2, na.rm = TRUE) +
  ggplot2::scale_x_continuous(breaks = deltas_levels) +
  ggplot2::labs(
    x = expression(Delta*T~"[°C]"),
    y = "Median der CWP-Flächengrösse [m²]",
    color = "Step [m]",
    linetype = if (has_buffer) "Puffer (Pixel)" else NULL,
    title = " "
  ) +
  ggplot2::theme_minimal(base_size = 12)


# Kombination & Export (alle vier nebeneinander)

if (!requireNamespace("cowplot", quietly = TRUE)) install.packages("cowplot")

# Extract legend from one representative plot
legend_plot <- cowplot::get_legend(
  p_count +
    ggplot2::theme(
      legend.position = "right",
      legend.title = ggplot2::element_text(size = 10),
      legend.text  = ggplot2::element_text(size = 9)
    )
)

# Remove legends from all plots
p_count      <- p_count      + ggplot2::theme(legend.position = "none")
p_total_area <- p_total_area + ggplot2::theme(legend.position = "none")
p_mean_area  <- p_mean_area  + ggplot2::theme(legend.position = "none")
p_med_area   <- p_med_area   + ggplot2::theme(legend.position = "none")

# Combine 4 plots horizontally
plots_row <- cowplot::plot_grid(
  p_count, p_total_area, p_med_area, p_mean_area,
  nrow = 1, align = "v"
)

# Add legend to the right of the combined plots
combined_4 <- cowplot::plot_grid(
  plots_row, legend_plot,
  ncol = 2, rel_widths = c(1, 0.12)
)

# # Add overall title
# title_grob <- grid::textGrob(
#   paste0("Sensitivität – Anzahl und Flächenkennwerte vs ΔT (", study_area_title, ")"),
#   gp = grid::gpar(fontsize = 14, fontface = "bold")
# )
#combined_4 <- cowplot::plot_grid(title_grob, combined_4, ncol = 1, rel_heights = c(0.1, 1))

out_base_4 <- file.path(plot_dir, paste0("linien_anzahl_gesamt_mittel_median_nach_delta_", pr, study_area_file))
save_plot_pdf_tiff(combined_4, out_base_4, width = 20, height = 5)
pr=pr+1










