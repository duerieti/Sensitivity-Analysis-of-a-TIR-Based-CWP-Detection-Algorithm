
# =============================================================================
# CWP DETECTION — OPTIMISATION LAB
# =============================================================================
# This script extracts the core cold-water-patch detection logic from the
# original sensitivity-analysis code and exposes it as a single-call function
# `detect_cwp_single()`.  No parameter-grid loop — you call it once with one
# set of parameters and get back the result + timing breakdown.
#
# Purpose: benchmark, profile, and optimise each stage of the detection
# pipeline without the overhead / complexity of the outer parameter sweep.
# =============================================================================

rm(list = ls())

# ── LIBRARIES ────────────────────────────────────────────────────────────────
library(terra)
library(sf)
library(dplyr)
library(lwgeom)

options(sf_use_s2 = FALSE)

# ── DATA PATHS (Emme reaches — relative to project root) ────────────────────



# =============================================================================
# CORE FUNCTION: detect_cwp_single()
# =============================================================================
# Runs the full CWP detection pipeline for ONE set of parameters.
# Returns a list with the result sf + a named vector of timings (seconds).
#
# This is the function you will profile and optimise.
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

  from <- pos[1:(npt - 1)]fc_GEOMETRY of length 2; first list element: List of 1
##   ..$ : num [1:11845, 1:2] -738336 -738357 -738363 -738368 -738378 ...
##   ..- attr(*, "class")= chr  "XY" "POLYGON" "sfg"
##  - attr(*, "sf_column")= chr "."
##  - attr(*, "agr")= Factor w/ 3 levels "constant","aggregate",..: NA
##   ..- attr(*, "names")= chr "centrum"

Lo and behold - a sf data frame with just two rows, one for center and one for the outskirts.

And while it was not my original intention the outskirts area is a single spatial object consisting of three separate polygons. Which I found rather neat.

Tagged: data-viz, spatial
JLA Data

© 2018 / Powered by Hugo

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
}


# =============================================================================
# CORE FUNCTION: detect_cwp_single()  —  ZONAL STATS REFACTOR
# =============================================================================
# Key changes vs original:
#
#   Step 2 (extract slab pixels):
#     OLD → terra::extract() materialises a giant flat table (ID, T, cell)
#           one row per pixel × slab — memory spike, bad_alloc risk
#     NEW → terra::rasterize() burns slab IDs into a raster (one integer per
#           pixel), then terra::zonal() computes medians entirely in C++.
#           Only a tiny k-row table ever enters R memory.
#
#   Step 3 (per-slab corridor median loop):
#     OLD → k separate terra::extract() calls inside a for-loop
#     NEW → GONE. Medians come directly from the zonal stats table above.
#           The corridor buffer distinction is preserved by using the slab
#           polygons themselves (same geometry, same intent).
#           NOTE: if you need a NARROWER corridor than slab_halfwidth_m for
#           the reference temperature, see the optional block at the end of
#           step 2 — you can rasterize a separate narrow-corridor layer.
#
#   Step 4 (flag cold pixels):
#     OLD → R-level for-loop comparing T_pixel vs Tmed[slab] row by row
#     NEW → pure raster algebra in C++:
#             Tref_r  = raster where each pixel holds its slab's median T
#             delta_r = Tref_r - r   (difference everywhere)
#             cold_r  = delta_r >= delta_C  (logical raster)
#           Cell indices extracted with terra::cells() — only flagged pixels
#           enter R memory, not the full pixel table.
#
#   Steps 5–7 (polygonize, area filter, patch stats):
#     UNCHANGED — these operate on small patch polygons, not the full raster.
# =============================================================================

detect_cwp_single_zonal <- function(
    ras_path,
    line_path,
    step_m            = 500,
    buffer_px         = 3,
    delta_C           = 1.0,
    min_patch_area_m2 = 2,
    round_to          = 0.1,
    slab_halfwidth_m  = 60,
    connect_diagonals = TRUE,
    union_chunk_size  = 50000L
) {

  timings <- numeric()
  tic <- function() proc.time()[["elapsed"]]


  # ── 0. READ INPUTS ─────────────────────────────────────────────────────────
  # Unchanged from original.
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

  # Crop raster to river corridor immediately — reduces all subsequent
  # raster operations to only the pixels that could ever be inside a slab.
  # This is the single most impactful memory/speed improvement.
  corridor_full <- sf::st_buffer(line, dist = slab_halfwidth_m + 10)
  r <- terra::crop(r, terra::vect(corridor_full))
  r <- terra::mask(r, terra::vect(corridor_full))

  timings[["read_inputs"]] <- tic() - t0


  # ── 1. BUILD SLABS ─────────────────────────────────────────────────────────
  # Unchanged from original — produces slabs_sf (sf of POLYGON geometries).
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
  for (j in seq_len(k))
    subs_list[[j]] <- lwgeom::st_linesubstring(g1, from[j], to[j])
  subs  <- do.call(c, subs_list)

  slabs <- sf::st_buffer(subs,
                         dist        = slab_halfwidth_m,
                         endCapStyle = "FLAT",
                         joinStyle   = "MITRE",
                         mitreLimit  = 2)
  slabs_sf      <- sf::st_sf(geometry = slabs)
  slabs_sf$slabID <- seq_len(k)   # explicit integer ID column — needed for rasterize

  timings[["build_slabs"]] <- tic() - t0


  # ── 2. RASTERIZE SLABS + ZONAL MEDIANS ────────────────────────────────────
  # NEW — replaces both the big terra::extract() AND the per-slab corridor loop.
  #
  # terra::rasterize(): burns the slabID value into each pixel that falls
  #   inside that slab. Pixels in NO slab get NA. Pixels in overlapping slabs
  #   get the LAST slab's ID (overlap is minimal with flat end caps).
  #   Output: a single-band integer raster, same extent/resolution as r.
  #
  # terra::zonal(): for each unique zone value (= slabID), computes the
  #   requested statistic over the corresponding pixels in r.
  #   Output: a small 2-column data frame — one row per slab.
  #   Everything stays in terra's C++ layer; no giant table in R memory.
  t0 <- tic()

  # Optional: if you want the reference temperature from a NARROWER corridor
  # (matching the original buffer_px logic), build narrow corridor polygons
  # here and rasterize those instead of slabs_sf.
  Rm <- buffer_px * cell_m
  if (Rm < slab_halfwidth_m) {
    # Narrow corridor for reference temperature only
    corridors <- sf::st_buffer(subs,
                               dist        = Rm,
                               endCapStyle = "FLAT",
                               joinStyle   = "MITRE",
                               mitreLimit  = 2)
    corridors_sf         <- sf::st_sf(geometry = corridors)
    corridors_sf$slabID  <- seq_len(k)
    zone_r <- terra::rasterize(terra::vect(corridors_sf), r, field = "slabID")
  } else {
    # Corridor same width as slab — rasterize slabs directly
    zone_r <- terra::rasterize(terra::vect(slabs_sf), r, field = "slabID")
  }

  # Apply optional rounding to r before computing medians
  r_med <- r
  if (!is.na(rfactor)) r_med <- round(r_med * rfactor) / rfactor

  # Compute median temperature per slab zone — result is a tiny data frame:
  #   col 1: slabID (zone value)
  #   col 2: T      (median temperature for that slab)
  zonal_df <- terra::zonal(r_med, zone_r, fun = "median", na.rm = TRUE)
  colnames(zonal_df) <- c("slabID", "Tmed")

  # Build Tmed vector indexed by slabID (1..k) for fast lookup below
  # Slabs with no pixels get NA (e.g. slab entirely outside water mask)
  Tmed <- rep(NA_real_, k)
  Tmed[zonal_df$slabID] <- zonal_df$Tmed

  timings[["zonal_ref_temps"]] <- tic() - t0


  # ── 3. BUILD T_ref RASTER + FLAG COLD PIXELS ──────────────────────────────
  # NEW — replaces the R-level for-loop over by_slab.
  #
  # terra::classify(): remaps zone_r pixel values using a lookup table.
  #   Each pixel that had value slabID now gets that slab's Tmed value.
  #   Result: Tref_r — a raster where every pixel holds its local reference T.
  #
  # Raster algebra then computes the anomaly and flags cold pixels,
  # entirely in C++ with no R-level loops.
  t0 <- tic()

  # Lookup table: from slabID → Tmed value
  # classify() expects a matrix with columns [from, to]
  lookup <- cbind(from = zonal_df$slabID,
                  to   = zonal_df$Tmed)
  Tref_r <- terra::classify(zone_r, lookup)

  # Difference raster: positive where pixel is COLDER than its slab median
  delta_r <- Tref_r - r

  # Logical raster: TRUE (1) where anomaly meets or exceeds threshold
  cold_r  <- delta_r >= delta_C

  # Extract cell indices of flagged pixels — only TRUE cells enter R memory
  flagged <- terra::cells(cold_r, 1)   # returns matrix with column "cell"
  print(flagged)
  print(class(flagged))
  flagged <- as.integer(flagged[, "cell"])

  
  timings[["flag_pixels"]] <- tic() - t0


  # ── 4. POLYGONIZE FLAGGED CELLS ────────────────────────────────────────────
  # Unchanged — polygonize_flagged_vector() still works on cell indices.
  t0 <- tic()

  if (!length(flagged)) {
    pol_sf <- sf::st_sf(geometry = sf::st_sfc(), crs = sf::st_crs(line))
  } else {
    pol_sf <- polygonize_flagged_vector(flagged, r,
                                        chunk_size   = union_chunk_size,
                                        connect_diag = connect_diagonals)
  }

  timings[["polygonize"]] <- tic() - t0


  # ── 5. AREA FILTER ─────────────────────────────────────────────────────────
  # Unchanged.
  t0 <- tic()

  if (inherits(pol_sf, "sfc"))
    pol_sf <- sf::st_sf(geometry = pol_sf, crs = sf::st_crs(line))

  if (nrow(pol_sf) > 0) {
    pol_sf$area_m2 <- tryCatch(
      as.numeric(sf::st_area(pol_sf)),
      error = function(e)
        as.numeric(terra::expanse(terra::vect(pol_sf),
                                  unit = "m", transform = FALSE))
    )
    pol_sf <- dplyr::filter(pol_sf, area_m2 >= min_patch_area_m2)
  }

  timings[["area_filter"]] <- tic() - t0


  # ── 6. TEMPERATURE STATISTICS PER PATCH ────────────────────────────────────
  # Largely unchanged — operates on small patch polygons so memory is fine.
  # One change: Tmed_slab is now looked up from the Tmed vector via zone_r
  # instead of reconstructing cell_to_slab from the old by_slab list.
  t0 <- tic()

  if (nrow(pol_sf) > 0) {
    pol_v     <- terra::vect(pol_sf)

    # Per-patch temperature stats (min, max, median, mean)
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

    # Assign each patch its slab reference temperature.
    # NEW approach: extract zone_r values (= slabIDs) inside each patch,
    # then look up Tmed[slabID]. No need to rebuild cell_to_slab from scratch.
    zone_vals     <- terra::extract(zone_r, pol_v)   # slabID per pixel per patch
    zone_vals_sp  <- split(zone_vals[["slabID"]], zone_vals$ID)

    pol_sf$Tmed_slab <- vapply(zone_vals_sp, function(ids) {
      ids <- ids[!is.na(ids)]
      if (!length(ids)) return(NA_real_)
      # median of the Tmed values of all slabs represented in this patch
      stats::median(Tmed[ids], na.rm = TRUE)
    }, numeric(1))

    # Confirmed delta: how much colder is the patch median vs its reference?
    pol_sf$deltaT       <- pol_sf$Tmed_slab - pol_sf$T_med
    pol_sf$delta_thresh <- delta_C

    # Second filter: patch-level delta must meet threshold
    pol_sf <- dplyr::filter(pol_sf, deltaT >= delta_C)
  }

  timings[["patch_stats"]] <- tic() - t0


  # ── TOTAL ───────────────────────────────────────────────────────────────────
  timings[["TOTAL"]] <- sum(timings[names(timings) != "TOTAL"])

  list(patches = pol_sf, timings = timings)
}


# =============================================================================
# HELPER: polygonize_flagged_vector()
# =============================================================================
# Converts flagged raster cell indices into dissolved vector polygons.
# Identical to the original — extracted here so detect_cwp_single() can call it.
# =============================================================================
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
# TIMING UTILITY: benchmark a single run and print a breakdown
# =============================================================================
print_timings <- function(result) {
  t <- result$timings
  total <- t[["TOTAL"]]
  cat("\n=== CWP Detection Timing Breakdown ===\n")
  for (nm in names(t)) {
    pct <- sprintf("(%5.1f%%)", t[[nm]] / total * 100)
    cat(sprintf("  %-25s %7.3f s  %s\n", nm, t[[nm]], if (nm != "TOTAL") pct else ""))
  }
  cat(sprintf("\n  Patches found: %d\n", nrow(result$patches)))
}


# =============================================================================
# RUN: Single detection on Emme V01 with default parameters
# =============================================================================
message("\n====== Running CWP detection on Emme V01 ======")
res_v01 <- detect_cwp_single(
  ras_path  = datasets$emme_v01$ras_path,
  line_path = datasets$emme_v01$line_path,
  step_m    = 500,
  buffer_px = 3,
  delta_C   = 1.0
)

print_timings(res_ob)
2 = FALSE)

# ── DATA PATHS (Emme reaches — relative to project root) ────────────────────

datasets <- list(
  emme_v01 = list(
    ras_path  = file.path("data/thermal_rasters_FINAL/mean_v01emme.tif"),
    line_path = file.path("data/Centerlines_FINAL/Emme_V01.shp")
  )
)



# Quick sanity check — do the files exist?
for (nm in names(datasets)) {
  ok_r <- file.exists(datasets[[nm]]$ras_path)
  ok_l <- file.exists(datasets[[nm]]$line_path)
  message(sprintf("  %-18s  raster: %s  |  centerline: %s", nm,
                  ifelse(ok_r, "OK", "MISSING"), ifelse(ok_l, "OK", "MISSING")))
}




message("\n====== Running CWP detection on Emme V01 ======")
res_v01 <- detect_cwp_single_zonal(
  ras_path  = "./data/emme_v0_quarter.tif",
  line_path = "./data/quarter_line.shp",
  step_m    = 500,
  buffer_px = 3,
  delta_C   = 1.0
)


print_timings(res_ob)


# =============================================================================
# BENCHMARK: Run N repetitions to get stable timing
# =============================================================================
benchmark_cwp <- function(ras_path, line_path, n_reps = 5, ...) {
  cat(sprintf("\nBenchmarking %d repetitions...\n", n_reps))
  all_timings <- list()
  for (i in seq_len(n_reps)) {
    res <- detect_cwp_single(ras_path, line_path, ...)
    all_timings[[i]] <- res$timings
    cat(sprintf("  rep %d: %.3f s  (%d patches)\n", i, res$timings[["TOTAL"]], nrow(res$patches)))
  }
  # Aggregate: median timing per stage
  stages <- names(all_timings[[1]])
  med <- sapply(stages, function(s) median(sapply(all_timings, `[[`, s)))
  cat("\n  Median timings:\n")
  for (s in stages) {
    cat(sprintf("    %-25s %7.3f s\n", s, med[[s]]))
  }
  invisible(all_timings)
}

# Uncomment to run a proper benchmark:
# benchmark_cwp(datasets$emme_v01$ras_path, datasets$emme_v01$line_path, n_reps = 5)
