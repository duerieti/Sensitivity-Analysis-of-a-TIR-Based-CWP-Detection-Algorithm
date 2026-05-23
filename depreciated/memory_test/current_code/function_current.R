library(terra)
library(sf)
library(dplyr)
library(lwgeom)
library(tidyverse)
library(tictoc)
library(tmap)


path_raster <- "../../../data/thermal_rasters_FINAL/mean_v01emme.tif"
path_line <- "../../../data/Centerlines_FINAL/Emme_V01.shp"

terraOptions(memmax = 9.8)


polygonize_flagged_vector <- function(flagged_cells, r,
                                      chunk_size   = 50000L,
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

  n      <- length(flagged_cells)
  chunks <- split(seq_len(n), ceiling(seq_len(n) / chunk_size))
  parts_sv <- vector("list", length(chunks))

  for (ci in seq_along(chunks)) {
    idx <- chunks[[ci]]
    xy  <- terra::xyFromCell(r, flagged_cells[idx])
    sfc <- make_rects_sfc(xy[, 1], xy[, 2], res[1], res[2], crs_wk)

    v       <- terra::vect(sf::st_sf(geometry = sfc))
    v$grp   <- 1L
    vd      <- terra::aggregate(v, by = "grp")
    parts_sv[[ci]] <- vd
  }

  parts_sf <- lapply(parts_sv, sf::st_as_sf)
  stitched <- do.call(rbind, parts_sf)
  geom     <- sf::st_geometry(stitched)
  geom     <- try(sf::st_make_valid(geom), silent = TRUE)

  eps <- if (connect_diag) min(res) * 0.1 else 0

  if (eps > 0) {
    merged <- sf::st_union(sf::st_buffer(geom, eps))
    merged <- sf::st_buffer(merged, -eps)
  } else {
    merged <- sf::st_union(geom)
  }

  polys <- sf::st_collection_extract(merged, "POLYGON", warn = FALSE)
  polys <- suppressWarnings(sf::st_cast(polys, "POLYGON"))
  sf::st_sf(geometry = polys)
}


# =============================================================================
detect_cwp_single <- function(
    ras_path,
    line_path,
    step_m            = 500,       # slab length along centerline (m)
    buffer_px         = 3,         # corridor half-width for reference T (pixels)
    delta_C           = 1.0,       # temperature anomaly threshold (deg C)
    min_patch_area_m2 = 2,         # drop patches smaller than this
    round_to          = 0.1,       # round temperatures to this precision
    slab_halfwidth_m  = 60,        # lateral half-width of slab polygons (m)
    connect_diagonals = TRUE,      # 8-connectivity when dissolving patches
    union_chunk_size  = 50000L     # max cells per dissolve chunk
) {

  timings <- numeric()
  tic <- function() proc.time()[["elapsed"]]

  # ── 0. READ INPUTS ────────────────────────────────────────────────────────
  t0 <- tic()

  r <- terra::rast(ras_path)
  stopifnot(terra::nlyr(r) == 1)
  names(r) <- "T"
  if (terra::is.lonlat(r)) stop("Raster must be in a metric CRS.")
  rres   <- terra::res(r)
  cell_m <- mean(rres)

  line <- sf::st_read(line_path, quiet = TRUE) |> sf::st_zm(TRUE, "ZM")
  if (!is.na(sf::st_crs(line)) &&
      !identical(terra::crs(r), sf::st_crs(line)$wkt)) {
    line <- sf::st_transform(line, terra::crs(r))
  }
  g <- sf::st_union(line)
  if (inherits(g, "sfc_GEOMETRYCOLLECTION"))
    g <- sf::st_collection_extract(g, "LINESTRING")
  if (inherits(g, "sfc_MULTILINESTRING")) {
    g <- sf::st_line_merge(g)
    if (inherits(g, "sfc_GEOMETRYCOLLECTION"))
      g <- sf::st_collection_extract(g, "LINESTRING")
  }
  if (!inherits(g, "sfc_LINESTRING"))
    stop("Centerline must resolve to a LINESTRING.")
  line <- sf::st_sf(geometry = g)

  rfactor <- if (!is.null(round_to) && round_to > 0) 1 / round_to else NA_real_

  timings[["read_inputs"]] <- tic() - t0

  # ── 1. BUILD SLABS ────────────────────────────────────────────────────────
  t0 <- tic()

  L   <- as.numeric(sf::st_length(line))
  S   <- step_m
  pos <- sort(unique(c(0, seq(0, 1, by = S / L), 1)))
  npt <- length(pos)
  if (npt < 2) stop("Too few stations for step = ", S)

  from <- pos[1:(npt - 1)]
  to   <- pos[2:npt]
  k    <- length(from)
  g1   <- sf::st_geometry(line)[1]

  subs_list <- vector("list", k)
  for (j in seq_len(k)) {
    subs_list[[j]] <- lwgeom::st_linesubstring(g1, from[j], to[j])
  }
  subs  <- do.call(c, subs_list)
  slabs <- sf::st_buffer(subs,
                          dist        = slab_halfwidth_m,
                          endCapStyle = "FLAT",
                          joinStyle   = "MITRE",
                          mitreLimit  = 2)
  slabs_sf <- sf::st_sf(geometry = slabs)

  timings[["build_slabs"]] <- tic() - t0

  # ── 2. EXTRACT SLAB PIXELS ────────────────────────────────────────────────
  t0 <- tic()

  slab_df <- terra::extract(r, terra::vect(slabs_sf), cells = TRUE)
  slab_df <- slab_df[!is.na(slab_df$T), , drop = FALSE]
  if (!nrow(slab_df)) stop("No pixels found inside slabs.")
  if (!is.na(rfactor)) slab_df$T <- round(slab_df$T * rfactor) / rfactor
  by_slab <- split(slab_df[, c("ID", "cell", "T")], slab_df$ID)

  timings[["extract_slab_pixels"]] <- tic() - t0

  # ── 3. COMPUTE REFERENCE TEMPERATURES (median per slab corridor) ───────
  t0 <- tic()

  Rm <- buffer_px * cell_m

  median_round <- function(x) {
    x <- x[!is.na(x)]
    if (!length(x)) return(NA_real_)
    if (!is.na(rfactor)) x <- round(x * rfactor) / rfactor
    stats::median(x)
  }

  Tmed <- numeric(k)
  for (j in seq_len(k)) {
    corridor_j <- sf::st_buffer(
      sf::st_sfc(subs[j], crs = sf::st_crs(line)),
      dist        = Rm,
      endCapStyle = "FLAT",
      joinStyle   = "MITRE",
      mitreLimit  = 2
    )
    vals_df  <- terra::extract(r, terra::vect(sf::st_sf(geometry = corridor_j)))
    v        <- vals_df[[names(r)[1]]]
    Tmed[j]  <- median_round(v)
  }

  timings[["compute_ref_temps"]] <- tic() - t0

  # ── 4. FLAG COLD PIXELS ───────────────────────────────────────────────────
  t0 <- tic()

  flagged <- integer(0L)
  for (j in seq_along(by_slab)) {
    tm <- Tmed[j]
    if (is.na(tm)) next
    dfj  <- by_slab[[j]]
    keep <- dfj$cell[(tm - dfj$T) >= delta_C]
    if (length(keep)) flagged <- c(flagged, keep)
  }
  flagged <- unique(flagged)





  timings[["flag_pixels"]] <- tic() - t0

  # ── 5. POLYGONIZE FLAGGED CELLS ────────────────────────────────────────────
  t0 <- tic()

  if (!length(flagged)) {
    pol_sf <- sf::st_sf(geometry = sf::st_sfc(), crs = sf::st_crs(line))
  } else {
    pol_sf <- polygonize_flagged_vector(flagged, r,
                                         chunk_size   = union_chunk_size,
                                         connect_diag = connect_diagonals)
  }

  timings[["polygonize"]] <- tic() - t0

  # ── 6. AREA FILTER (drop small patches) ───────────────────────────────────
  t0 <- tic()

  if (inherits(pol_sf, "sfc"))
    pol_sf <- sf::st_sf(geometry = pol_sf, crs = sf::st_crs(line))

  if (nrow(pol_sf) > 0) {
    pol_sf$area_m2 <- tryCatch(
      as.numeric(sf::st_area(pol_sf)),
      error = function(e)
        as.numeric(terra::expanse(terra::vect(pol_sf), unit = "m", transform = FALSE))
    )
    pol_sf <- dplyr::filter(pol_sf, area_m2 >= min_patch_area_m2)
  }

  timings[["area_filter"]] <- tic() - t0

  # ── 7. TEMPERATURE STATISTICS PER PATCH ────────────────────────────────────
  t0 <- tic()

  if (nrow(pol_sf) > 0) {
    pol_v     <- terra::vect(pol_sf)
    vals_list <- terra::extract(r, pol_v)
    vals_list <- split(vals_list[[names(r)[1]]], vals_list$ID)

    get_stats <- function(x) {
      x <- x[!is.na(x)]
      if (!length(x)) return(rep(NA_real_, 4))
      c(min(x), max(x), stats::median(x), mean(x))
    }
    stats_mat           <- t(vapply(vals_list, get_stats, numeric(4)))
    colnames(stats_mat) <- c("T_min", "T_max", "T_med", "T_mean")
    pol_sf$T_min  <- stats_mat[, "T_min"]
    pol_sf$T_max  <- stats_mat[, "T_max"]
    pol_sf$T_med  <- stats_mat[, "T_med"]
    pol_sf$T_mean <- stats_mat[, "T_mean"]

    # Assign slab reference temperature to each patch
    cell_extract <- terra::extract(r, pol_v, cells = TRUE)
    cell_extract <- cell_extract[!is.na(cell_extract$cell), ]

    cell_to_slab_list <- lapply(seq_along(by_slab), function(j) {
      df <- by_slab[[j]]
      if (nrow(df) == 0) return(NULL)
      data.frame(cell = df$cell, slab_id = j)
    })
    cell_to_slab <- do.call(rbind, cell_to_slab_list)

    merged_df <- merge(cell_extract, cell_to_slab, by = "cell", all.x = TRUE)

    slab_meds <- tapply(merged_df$slab_id, merged_df$ID, function(ids) {
      if (all(is.na(ids))) return(NA_real_)
      stats::median(Tmed[ids], na.rm = TRUE)
    })
    pol_sf$Tmed_slab <- as.numeric(slab_meds[as.character(seq_len(nrow(pol_sf)))])

    pol_sf$deltaT       <- pol_sf$Tmed_slab - pol_sf$T_med
    pol_sf$delta_thresh <- delta_C

    # Second filter: patch-level delta must meet threshold
    pol_sf <- dplyr::filter(pol_sf, deltaT >= delta_C)
  }

  timings[["patch_stats"]] <- tic() - t0

  # ── TOTAL ──────────────────────────────────────────────────────────────────
  timings[["TOTAL"]] <- sum(timings)

  list(patches = pol_sf, timings = timings)


  sf::write_sf(pol_sf, "final_polys.shp")
}



tic("Current Function:")
detect_cwp_single(path_raster,path_line)
toc()

