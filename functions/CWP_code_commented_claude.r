# =============================================================================
# COLD WATER PATCH (CWP) DETECTION — THERMAL RIVER ANALYSIS
# =============================================================================
# PURPOSE:
#   Detect anomalously cold zones ("cold water patches") along river reaches
#   using thermal infrared raster data. Cold patches are areas where the local
#   water temperature is significantly lower than the surrounding river reach —
#   often caused by cold groundwater inflows, tributary confluences, or upwelling.
#
# KEY CONCEPT:
#   For every combination of parameters (step × buffer × ΔT), the script:
#     1. Cuts the river centerline into segments ("slabs")
#     2. Computes a local reference temperature (median) per slab
#     3. Flags pixels that are ≥ ΔT colder than their slab's reference temp
#     4. Polygonizes the flagged pixels into discrete patches
#     5. Saves one shapefile per parameter combination + a sensitivity summary
#
# SCIENTIFIC CONTEXT:
#   Part of project Thermo_BE_2024_2025 — thermal analysis of rivers in the
#   canton of Bern, Switzerland. Cold refugia are ecologically important for
#   cold-water fish (e.g. trout) under heat stress.
# =============================================================================

# Clear the workspace to avoid conflicts with previous sessions
rm(list = ls())

# ── LIBRARIES ─────────────────────────────────────────────────────────────────
library(terra)    # raster reading, extraction, reprojection (replaces raster pkg)
library(sf)       # vector geometry (points, lines, polygons), CRS handling
library(dplyr)    # data manipulation (filter, mutate, summarise, etc.)
library(tidyr)    # data reshaping (not heavily used here but loaded)
library(purrr)    # functional programming helpers (map, lapply alternatives)
library(ggplot2)  # plotting sensitivity results
library(stringr)  # string manipulation utilities
library(readr)    # fast CSV writing (write_csv)
library(lwgeom)   # st_linesubstring: cut a line at fractional positions
library(qs)       # fast R object serialisation (qsave/qread) — much faster than saveRDS
library(grDevices)# TIFF and PDF device drivers for saving plots

# Disable spherical geometry engine — use flat (planar) geometry instead.
# Important because the raster is in a projected (metric) CRS, not geographic.
options(sf_use_s2 = FALSE)


# =============================================================================
# USER INPUTS — set working directory and file paths per study area
# =============================================================================
# Each block below configures one river reach. Only the LAST executed block
# takes effect — the others are overwritten. In practice, you run one block
# at a time (comment out the rest) before calling detect_cold_patches().

# --- Emme lower reach (V1) ---
setwd("D:/Daten_anto/Thermo_bern/automatisation/v1_emme")
ras_path  <- r"(s:\pools\n\...\mean_v01emme.tif)"   # thermal raster (water pixels only)
line_path <- r"(s:\pools\n\...\Emme_V01.shp)"        # river centerline shapefile
out_dir   <- "output/cold_patches_out_linear_buffer_emmev1_01"

# --- Urtene V1 ---
setwd("D:/Daten_anto/Thermo_bern/automatisation/v1_urtene/")
ras_path  <- r"(s:\pools\n\...\mean_v1urtene.tif)"
line_path <- r"(s:\pools\n\...\Urtene_V01.shp)"
out_dir   <- "output/cold_patches_out_linear_buffer_urtenev1_01"

# --- Grüene V2 ---
setwd("D:/Daten_anto/Thermo_bern/automatisation/v2_emmgrue/")
ras_path  <- r"(s:\pools\n\...\mean_v2gruen.tif)"
line_path <- r"(s:\pools\n\...\Gruene_V02.shp)"
out_dir   <- "output/cold_patches_out_linear_buffer_gruene_01"

# --- Emme upper reach 2 (V2) ---
setwd("D:/Daten_anto/Thermo_bern/automatisation/v2_emmgrue/")
ras_path  <- r"(s:\pools\n\...\mean_v2emme.tif)"
line_path <- r"(s:\pools\n\...\Emme_V02.shp)"
out_dir   <- "output/cold_patches_out_linear_buffer_emmev2_01"

# --- Ilfis V2 ---
setwd("D:/Daten_anto/Thermo_bern/automatisation/v2_ilf/")
ras_path  <- r"(s:\pools\n\...\mean_v02ilfis.tif)"
line_path <- r"(s:\pools\n\...\Ilfis_V02.shp)"
out_dir   <- "output/cold_patches_out_linear_buffer_ilfv2_01"

# --- Obere Emme (upper reach 1, V2) ---
setwd("D:/Daten_anto/Thermo_bern/automatisation/v2_obereemme/")
ras_path  <- r"(s:\pools\n\...\mean_v2obemme.tif)"
line_path <- r"(s:\pools\n\...\Obere_Emme_V02.shp)"
out_dir   <- "output/cold_patches_out_linear_buffer_v2_obereemme_01"


# ── PARAMETER GRID ────────────────────────────────────────────────────────────
# The sensitivity analysis loops over ALL combinations of these three parameters.
# Total runs = length(steps_m) × length(buffers_px) × length(deltas_C)
# With these defaults: 6 × 1 × 4 = 24 combinations per river.

steps_m    <- c(100, 200, 500, 1000, 2000, 5000)
# ↑ Slab length along the centerline (meters).
#   Controls the spatial resolution of the LOCAL reference temperature.
#   Small step (100m) → many short slabs → T_ref adapts tightly to local conditions.
#   Large step (5000m) → few long slabs → T_ref averages over a long reach.

buffers_px <- c(3)
# ↑ Half-width of the corridor used to compute the reference temperature median,
#   expressed in PIXELS (converted to meters internally: Rm = Rpx × cell_size).
#   Buffer of 3px at 0.5m resolution = 1.5m corridor on each side of the line.
#   A narrow buffer → reference temp comes from pixels very close to the centerline.
#   A wide buffer → reference temp includes more of the river cross-section.

deltas_C   <- c(0.5, 1.0, 2.0, 3.0)
# ↑ Temperature anomaly threshold ΔT [°C].
#   A pixel is flagged as "cold" if: T_ref − T_pixel ≥ ΔT
#   Low ΔT (0.5°C) → catches subtle cold anomalies, more patches, risk of false positives.
#   High ΔT (3.0°C) → only strong cold inflows are detected, fewer but more certain patches.

min_area_m2  <- 2        # Minimum patch area in m² — tiny pixel artifacts are discarded.
handle_small <- "drop"   # What to do with sub-threshold patches: "drop" removes them,
                          # "merge" attaches them to the nearest larger patch.


# =============================================================================
# HELPER: LOGGER
# =============================================================================
# Prints a timestamped message to the console, useful for tracking progress
# during long multi-parameter runs.
log_msg <- function(fmt, ...) message(sprintf("[%s] %s",
                                              format(Sys.time(), "%H:%M:%S"),
                                              sprintf(fmt, ...)))


# =============================================================================
# MAIN FUNCTION: detect_cold_patches()
# =============================================================================
detect_cold_patches <- function(
    ras_path,              # path to the thermal raster (.tif)
    line_path,             # path to the river centerline (.shp)
    steps          = c(100, 200, 500, 1000),   # slab lengths to test (m)
    buffers_px     = c(1, 3, 5, 10),           # corridor half-widths to test (px)
    deltas         = c(0.5, 1, 2, 3),          # ΔT thresholds to test (°C)
    min_patch_area_m2 = 2,                     # discard patches smaller than this
    handle_small   = c("drop", "merge"),       # how to handle small patches
    round_to       = 0.1,   # round temperatures to nearest 0.1°C before computing medians
                             # reduces sensitivity to raster floating-point noise
    slab_halfwidth_m = 60,  # lateral half-width of each slab polygon (m)
                             # determines how far from the centerline pixels are included
    out_dir        = "output/cold_patches_out_linear_buffer",
    connect_diagonals = TRUE,    # use 8-neighbour connectivity when dissolving patches
                                  # TRUE: diagonally touching pixels merge into one patch
                                  # FALSE: only cardinal neighbours merge (4-connectivity)
    union_chunk_size = 50000L    # max rectangles per dissolve chunk (memory management)
                                  # large rasters with many flagged cells are processed in
                                  # batches to avoid RAM exhaustion during st_union()
) {

  # Local copy of logger (so the function is self-contained)
  log_msg <- function(fmt, ...) message(sprintf("[%s] %s",
                                                format(Sys.time(), "%H:%M:%S"),
                                                sprintf(fmt, ...)))

  options(sf_use_s2 = FALSE)
  handle_small <- match.arg(handle_small)  # validate argument, error if neither "drop" nor "merge"

  # Create output directories (won't fail if they already exist)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  shp_dir <- file.path(out_dir, "shapefiles")
  dir.create(shp_dir, showWarnings = FALSE)


  # ── READ AND VALIDATE INPUTS ───────────────────────────────────────────────

  log_msg("Reading raster: %s", ras_path)
  r <- terra::rast(ras_path)
  stopifnot(terra::nlyr(r) == 1)   # must be single-band (one temperature layer)
  names(r) <- "T"                   # rename band to "T" for consistent access later
  if (terra::is.lonlat(r)) stop("Raster must be in a metric CRS.")
  # ↑ Distances (slab lengths, buffer widths) are in meters — a geographic CRS
  #   (lat/lon degrees) would silently produce wrong results.
  rres   <- terra::res(r)           # [xres, yres] in meters
  cell_m <- mean(rres)              # average cell size — used to convert px → meters

  log_msg("Reading centerline: %s", line_path)
  line <- sf::st_read(line_path, quiet = TRUE) |>
    sf::st_zm(TRUE, "ZM")           # drop Z/M coordinates if present (2D only needed)

  # Reproject centerline to match the raster CRS if they differ
  if (!is.na(sf::st_crs(line)) &&
      !identical(terra::crs(r), sf::st_crs(line)$wkt)) {
    log_msg("Reprojecting centerline to raster CRS")
    line <- sf::st_transform(line, terra::crs(r))
  }

  # The shapefile may contain multiple features (e.g. separate line segments).
  # We dissolve them into a single continuous LINESTRING for consistent sampling.
  g <- sf::st_union(line)
  if (inherits(g, "sfc_GEOMETRYCOLLECTION"))
    g <- sf::st_collection_extract(g, "LINESTRING")
  if (inherits(g, "sfc_MULTILINESTRING")) {
    g <- sf::st_line_merge(g)       # stitch segments into one line where endpoints touch
    if (inherits(g, "sfc_GEOMETRYCOLLECTION"))
      g <- sf::st_collection_extract(g, "LINESTRING")
  }
  if (!inherits(g, "sfc_LINESTRING"))
    stop("Centerline must resolve to a LINESTRING.")
  line <- sf::st_sf(geometry = g)

  # Pre-compute rounding factor once (avoids repeated if-checks inside loops)
  # rfactor = 10 when round_to = 0.1  →  round(x * 10) / 10
  rfactor <- if (!is.null(round_to) && round_to > 0) 1 / round_to else NA_real_

  # Report total number of parameter combinations to process
  total_combos <- length(steps) * length(buffers_px) * length(deltas)
  combo_i      <- 0L
  log_msg("Grid: steps=%s | buffers_px=%s | deltas=%s | combos=%d",
          paste(steps, collapse = ","), paste(buffers_px, collapse = ","),
          paste(deltas, collapse = ","), total_combos)

  # Accumulators for results across all combos
  all_summaries <- list()   # one tibble row per combo
  all_areas     <- list()   # per-patch area rows (long format)
  written_shps  <- character(0)  # paths to successfully written shapefiles


  # ============================================================================
  # INNER HELPER: polygonize_flagged_vector()
  # ============================================================================
  # PURPOSE:
  #   Convert a set of flagged raster cell indices into dissolved vector polygons.
  #   Works entirely in vector space (no rasterize → polygonize raster operation).
  #
  # WHY CHUNKED?
  #   If thousands of cells are flagged, calling st_union() on all of them at once
  #   can exhaust RAM. We dissolve in batches (chunks), then do a final global union.
  #
  # WHY DIAGONAL CONNECTION?
  #   By default, two pixels sharing only a corner (diagonal neighbour) are treated
  #   as disconnected patches. With connect_diagonals=TRUE, we expand each pixel by
  #   a tiny epsilon before dissolving so that diagonal neighbours touch and merge.
  #   After dissolving, we shrink back by the same epsilon to restore original shape.
  # ============================================================================
  polygonize_flagged_vector <- function(flagged_cells, r,
                                        chunk_size  = 50000L,
                                        connect_diag = TRUE) {

    # Return an empty sf object if no cells are flagged
    if (!length(flagged_cells)) {
      return(sf::st_sf(geometry = sf::st_sfc(),
                       crs      = sf::st_crs(terra::crs(r))))
    }

    res    <- terra::res(r)     # pixel dimensions [xres, yres]
    crs_wk <- terra::crs(r)

    # Sub-helper: build a list of square polygon geometries from cell centres.
    # Each cell centre (xc, yc) becomes a closed rectangle of size rx × ry.
    make_rects_sfc <- function(xc, yc, rx, ry, crs_wk) {
      polys <- vector("list", length(xc))
      for (i in seq_along(xc)) {
        x <- xc[i]; y <- yc[i]
        # 5-point closed ring: SW → SE → NE → NW → SW
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
    # Split cell indices into batches of chunk_size
    chunks <- split(seq_len(n), ceiling(seq_len(n) / chunk_size))
    parts_sv <- vector("list", length(chunks))

    # ── STEP 1: Dissolve within each chunk ────────────────────────────────────
    # For each batch: build pixel rectangles → convert to terra SpatVector →
    # dissolve all into one multipolygon using terra::aggregate().
    # terra::aggregate() is faster than sf::st_union() for large datasets.
    for (ci in seq_along(chunks)) {
      idx <- chunks[[ci]]
      xy  <- terra::xyFromCell(r, flagged_cells[idx])   # cell index → (x, y) coordinates
      sfc <- make_rects_sfc(xy[, 1], xy[, 2], res[1], res[2], crs_wk)

      v       <- terra::vect(sf::st_sf(geometry = sfc))
      v$grp   <- 1L                              # all cells in same group → dissolve together
      vd      <- terra::aggregate(v, by = "grp") # dissolve: adjacent rects merge
      parts_sv[[ci]] <- vd
      message(sprintf("     chunk %d/%d: %d cells → %d multipart piece(s)",
                      ci, length(chunks), length(idx), nrow(vd)))
    }

    # ── STEP 2: Global dissolve across all chunks ─────────────────────────────
    # Each chunk produced a partially-dissolved multipolygon. Now union them all
    # together to merge patches that span chunk boundaries.
    parts_sf <- lapply(parts_sv, sf::st_as_sf)
    stitched <- do.call(rbind, parts_sf)
    geom     <- sf::st_geometry(stitched)
    geom     <- try(sf::st_make_valid(geom), silent = TRUE)  # fix any geometry errors

    # Diagonal connection trick:
    # eps = 10% of minimum pixel dimension (tiny, sub-pixel expansion)
    # Buffer OUT by eps → diagonal neighbours now touch → st_union merges them
    # Buffer IN by eps → restore original pixel footprint
    eps <- if (connect_diag) min(res) * 0.1 else 0

    if (eps > 0) {
      merged <- sf::st_union(sf::st_buffer(geom, eps))
      merged <- sf::st_buffer(merged, -eps)
    } else {
      merged <- sf::st_union(geom)
    }

    # ── STEP 3: Explode multipart → single-part polygons ─────────────────────
    # st_union returns one big MULTIPOLYGON. We split it so each spatially
    # disconnected patch becomes its own row — needed for per-patch area/stats.
    polys <- sf::st_collection_extract(merged, "POLYGON", warn = FALSE)
    polys <- suppressWarnings(sf::st_cast(polys, "POLYGON"))  # guarantee sfc_POLYGON
    sf::st_sf(geometry = polys)
  }


  # ============================================================================
  # OUTER LOOP: iterate over slab step lengths
  # ============================================================================
  for (S in steps) {

    L <- as.numeric(sf::st_length(line))  # total centerline length in meters

    # Skip if the line is shorter than 2 slab lengths
    # (would produce only 0 or 1 slabs — not meaningful)
    if (L < 2 * S) { log_msg("Skip step=%d (line too short)", S); next }

    # Sample fractional positions along the line [0, 1].
    # pos = 0 is the start, pos = 1 is the end. Steps are evenly spaced at S/L intervals.
    # sort(unique(c(0, ..., 1))) ensures start and end are always included.
    pos <- sort(unique(c(0, seq(0, 1, by = S / L), 1)))
    pts <- sf::st_line_sample(line, sample = pos) |> sf::st_cast("POINT")
    pts_sf <- sf::st_as_sf(pts)
    tvec <- pos
    npt  <- length(tvec)
    if (npt < 2) { log_msg("Skip step=%d (stations=%d)", S, npt); next }
    log_msg("Step %d m: stations=%d", S, npt)

    # Build slab polygons: one per consecutive pair of positions.
    # Each slab = a short line segment buffered laterally to slab_halfwidth_m.
    # endCapStyle="FLAT" → flat ends at the station points (no rounded caps),
    #   so adjacent slabs tile perfectly without gaps or overlaps.
    from <- tvec[1:(npt - 1)]   # start fractions for each slab
    to   <- tvec[2:npt]         # end fractions for each slab
    k    <- length(from)        # number of slabs
    g1   <- sf::st_geometry(line)[1]
    log_msg("Building %d slabs (length≈%d m, halfwidth=%.2f m)", k, S, slab_halfwidth_m)

    subs_list <- vector("list", k)
    for (j in seq_len(k)) {
      # Cut the centerline between positions from[j] and to[j]
      # Returns a LINESTRING covering exactly that fraction of the river
      subs_list[[j]] <- lwgeom::st_linesubstring(g1, from[j], to[j])
    }
    subs  <- do.call(c, subs_list)   # combine into an sfc (list of geometries)
    slabs <- sf::st_buffer(subs,     # buffer each substring laterally
                           dist         = slab_halfwidth_m,
                           endCapStyle  = "FLAT",
                           joinStyle    = "MITRE",   # sharp corners at bends
                           mitreLimit   = 2)
    slabs_sf <- sf::st_sf(geometry = slabs)

    # Extract ALL raster pixels that fall within any slab.
    # slab_df columns: ID (slab index 1..k), cell (raster cell number), T (temperature)
    # This is done ONCE per step (not per buffer/delta) for efficiency.
    log_msg("Extracting slab pixels (full raster) …")
    slab_df <- terra::extract(r, terra::vect(slabs_sf), cells = TRUE)
    slab_df <- slab_df[!is.na(slab_df$T), , drop = FALSE]  # remove NA (non-water) pixels
    if (!nrow(slab_df)) {
      log_msg("No slab pixels at step=%d → skipping this step", S); next
    }

    # Round temperatures to reduce floating-point noise
    if (!is.na(rfactor)) slab_df$T <- round(slab_df$T * rfactor) / rfactor

    # Split the big extraction table by slab ID → list of data frames, one per slab.
    # by_slab[[j]] contains the cell indices and temperatures of all pixels in slab j.
    by_slab <- split(slab_df[, c("ID", "cell", "T")], slab_df$ID)
    log_msg("Slab rows=%d (unique cells=%d)", nrow(slab_df), length(unique(slab_df$cell)))


    # ── INNER LOOP 1: buffer width (controls reference temperature corridor) ──
    for (Rpx in buffers_px) {

      Rm <- Rpx * cell_m   # convert pixels → meters
      log_msg("  Linear buffer %d px (~%.2f m): medians per %d slabs …", Rpx, Rm, k)

      # Helper: median of a vector, with optional rounding and NA removal
      median_round <- function(x) {
        x <- x[!is.na(x)]
        if (!length(x)) return(NA_real_)
        if (!is.na(rfactor)) x <- round(x * rfactor) / rfactor
        stats::median(x)
      }

      # Compute the reference temperature Tmed[j] for each slab j.
      # The reference is the MEDIAN temperature inside a narrow corridor
      # (half-width = Rm) around the slab's centerline substring.
      #
      # WHY MEDIAN instead of mean?
      #   The median is robust to outliers (including the cold patches we're looking
      #   for). Using the mean would pull T_ref down near cold inflows, which would
      #   reduce their detected ΔT and potentially hide them. The median stays stable.
      #
      # NOTE: the corridor (Rm) is typically much narrower than the slab (slab_halfwidth_m).
      #   slab_halfwidth_m = 60m: defines which pixels belong to this river reach.
      #   Rm = buffer_px × cell_size: defines which pixels inform T_ref.
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
      log_msg("  Medians ready (NA=%d)", sum(is.na(Tmed)))


      # ── INNER LOOP 2: temperature threshold ─────────────────────────────────
      for (dC in deltas) {

        combo_i   <- combo_i + 1L
        combo_tag <- sprintf("S%sm_B%spx_D%0.1fC", S, Rpx, dC)  # e.g. "S100m_B3px_D1.0C"
        log_msg("   [%d/%d] ΔT=%.1f°C → %s", combo_i, total_combos, dC, combo_tag)

        # ── FLAG COLD PIXELS ────────────────────────────────────────────────
        # For each slab j: compare every pixel temperature against Tmed[j].
        # Flag the pixel if  T_ref − T_pixel ≥ dC  (i.e. pixel is at least dC colder).
        # We collect the CELL INDICES of flagged pixels (not their coordinates yet).
        flagged <- integer(0L)
        for (j in seq_along(by_slab)) {
          tm  <- Tmed[j]
          if (is.na(tm)) next           # skip slabs with no reference temperature
          dfj  <- by_slab[[j]]
          keep <- dfj$cell[(tm - dfj$T) >= dC]   # flagging condition
          if (length(keep)) flagged <- c(flagged, keep)
        }
        # Remove duplicates: a cell overlapping two slabs might be flagged twice
        flagged <- unique(flagged)
        log_msg("   Flagged cells: %d", length(flagged))


        # ── POLYGONIZE ──────────────────────────────────────────────────────
        # Convert the set of flagged cell indices into dissolved vector polygons.
        # Output: sf data frame with one row per spatially disconnected cold patch.
        if (!length(flagged)) {
          pol_sf <- sf::st_sf(geometry = sf::st_sfc(), crs = sf::st_crs(line))
          log_msg("   No polygons for %s", combo_tag)
        } else {
          pol_sf <- polygonize_flagged_vector(flagged, r,
                                              chunk_size   = union_chunk_size,
                                              connect_diag = connect_diagonals)
          log_msg("   Patches (single-part) created: %d", nrow(pol_sf))
        }


        # ── COMPUTE PATCH AREAS ─────────────────────────────────────────────
        # Areas are computed AFTER dissolving and splitting, so each row is one
        # spatially contiguous patch. Fallback to terra::expanse() if st_area fails.
        if (inherits(pol_sf, "sfc"))
          pol_sf <- sf::st_sf(geometry = pol_sf, crs = sf::st_crs(line))

        if (nrow(pol_sf) > 0) {
          pol_sf$area_m2 <- tryCatch(
            as.numeric(sf::st_area(pol_sf)),
            error = function(e)
              as.numeric(terra::expanse(terra::vect(pol_sf),
                                        unit = "m", transform = FALSE))
          )

          # ── HANDLE SMALL PATCHES ──────────────────────────────────────────
          # Option "drop": simply remove patches below min_patch_area_m2.
          #   Good default — avoids single-pixel noise in results.
          # Option "merge": attach small patches to the spatially nearest big patch.
          #   Useful if you want to preserve the total area even for patchy regions.
          if (handle_small == "drop") {
            pol_sf <- dplyr::filter(pol_sf, area_m2 >= min_patch_area_m2)
            log_msg("   Dropped patches < %.2f m²; kept: %d", min_patch_area_m2, nrow(pol_sf))
          } else {
            big   <- dplyr::filter(pol_sf, area_m2 >= min_patch_area_m2)
            small <- dplyr::filter(pol_sf, area_m2 <  min_patch_area_m2)
            if (nrow(small) > 0 && nrow(big) == 0) {
              # Edge case: ALL patches are small — keep only the largest one
              pol_sf <- dplyr::slice_max(small, area_m2, n = 1)
              log_msg("   No big patches → kept only largest small patch")
            } else if (nrow(small) > 0 && nrow(big) > 0) {
              # Find the nearest big patch for each small patch and merge them
              nearest_idx  <- sf::st_nearest_feature(small, big)
              small$merge_id <- nearest_idx
              big$merge_id   <- seq_len(nrow(big))
              pol_sf <- dplyr::bind_rows(big, small) |>
                dplyr::group_by(merge_id) |>
                dplyr::summarise(area_m2 = sum(area_m2), .groups = "drop")
              log_msg("   Merged %d small patches → now %d patches", nrow(small), nrow(pol_sf))
            } else {
              pol_sf <- big  # no small patches — nothing to do
            }
          }
        }


        # ── TEMPERATURE STATISTICS PER PATCH ───────────────────────────────
        # For each surviving patch polygon, extract the raster values inside it
        # and compute descriptive statistics. Also compute deltaT (confirmed ΔT)
        # and apply a SECOND filter to remove patches whose actual ΔT < threshold.
        #
        # WHY A SECOND FILTER?
        #   The first flagging used per-pixel comparison (T_ref − T_pixel ≥ dC).
        #   After dissolving, a polygon may include some pixels that just barely
        #   missed the threshold or were merged in from adjacent slabs. The second
        #   filter ensures the PATCH MEDIAN is truly dC colder than its slab reference.
        if (nrow(pol_sf) > 0) {
          log_msg("   Calculating temperature statistics for %d polygons …", nrow(pol_sf))

          pol_v <- terra::vect(pol_sf)

          # Extract temperatures for all pixels within each polygon
          vals_list <- terra::extract(r, pol_v)
          vals_list <- split(vals_list[[names(r)[1]]], vals_list$ID)

          # Compute min, max, median, mean for each polygon
          get_stats <- function(x) {
            x <- x[!is.na(x)]
            if (!length(x)) return(rep(NA_real_, 4))
            c(min(x), max(x), stats::median(x), mean(x))
          }
          stats_mat         <- t(vapply(vals_list, get_stats, numeric(4)))
          colnames(stats_mat) <- c("T_min", "T_max", "T_med", "T_mean")
          pol_sf$T_min  <- stats_mat[, "T_min"]
          pol_sf$T_max  <- stats_mat[, "T_max"]
          pol_sf$T_med  <- stats_mat[, "T_med"]
          pol_sf$T_mean <- stats_mat[, "T_mean"]

          # Assign each polygon its slab's reference temperature (Tmed_slab).
          # A polygon may span multiple slabs (especially with coarse steps),
          # so we take the MEDIAN of all the Tmed[j] values of its constituent cells.
          cell_extract <- terra::extract(r, pol_v, cells = TRUE)
          cell_extract <- cell_extract[!is.na(cell_extract$cell), ]

          # Build a lookup table: cell index → slab index
          cell_to_slab_list <- lapply(seq_along(by_slab), function(j) {
            df <- by_slab[[j]]
            if (nrow(df) == 0) return(NULL)
            data.frame(cell = df$cell, slab_id = j)
          })
          cell_to_slab <- do.call(rbind, cell_to_slab_list)

          # Join: for each cell inside a polygon, find its slab
          merged_df <- merge(cell_extract, cell_to_slab, by = "cell", all.x = TRUE)

          # For each polygon: take the median of Tmed values across its cells' slabs
          slab_meds <- tapply(merged_df$slab_id, merged_df$ID, function(ids) {
            if (all(is.na(ids))) return(NA_real_)
            stats::median(Tmed[ids], na.rm = TRUE)
          })
          pol_sf$Tmed_slab <- as.numeric(slab_meds[as.character(seq_len(nrow(pol_sf)))])

          # deltaT: how much colder is this patch compared to its local reference?
          # Positive value = patch is colder than reference (expected for CWP).
          pol_sf$deltaT       <- pol_sf$Tmed_slab - pol_sf$T_med
          pol_sf$delta_thresh <- dC   # record which threshold was used

          # SECOND FILTER: keep only polygons where the actual patch-level ΔT ≥ threshold
          pol_sf <- dplyr::filter(pol_sf, deltaT >= dC)
          log_msg("   Temperature statistics added to the features")
        }


        # ── WRITE SHAPEFILE ─────────────────────────────────────────────────
        # One .shp per parameter combination, named by the combo tag.
        # Overwrites any previous file for the same combo.
        if (nrow(pol_sf) > 0) {
          shp_path <- file.path(shp_dir, paste0(combo_tag, ".shp"))
          if (file.exists(shp_path))
            unlink(sub("\\.shp$", ".*", shp_path))  # delete old sidecar files (.dbf, .prj, etc.)
          ok <- try(sf::st_write(pol_sf, shp_path,
                                 driver = "ESRI Shapefile", quiet = TRUE),
                    silent = TRUE)
          if (inherits(ok, "try-error")) {
            log_msg("   WARN: failed to write shapefile %s", shp_path)
          } else {
            log_msg("   Wrote shapefile: %s", shp_path)
            written_shps <- c(written_shps, shp_path)
          }
        } else {
          log_msg("   Skipped shapefile write (no polygons)")
        }


        # ── ACCUMULATE SUMMARY STATISTICS ───────────────────────────────────
        # One row per parameter combination, recording patch count and area metrics.
        # These feed the sensitivity plots later.
        if (nrow(pol_sf) > 0) {
          all_summaries[[length(all_summaries) + 1]] <- tibble::tibble(
            step_m         = S,
            buffer_px      = Rpx,
            delta          = dC,
            n_patches      = nrow(pol_sf),
            total_area_m2  = sum(pol_sf$area_m2),
            median_area_m2 = stats::median(pol_sf$area_m2),
            mean_area_m2   = mean(pol_sf$area_m2)
          )
          # Also save individual patch areas in long format (for distribution plots)
          all_areas[[length(all_areas) + 1]] <-
            dplyr::transmute(pol_sf,
                             step_m    = S,
                             buffer_px = Rpx,
                             delta     = dC,
                             area_m2   = area_m2)
        } else {
          # Record zeros for combos that produced no patches
          all_summaries[[length(all_summaries) + 1]] <- tibble::tibble(
            step_m         = S, buffer_px = Rpx, delta = dC,
            n_patches      = 0, total_area_m2 = 0,
            median_area_m2 = NA_real_, mean_area_m2 = NA_real_
          )
        }

      } # ── end delta loop
    }   # ── end buffer loop

    log_msg("Finished step %d m", S)
  }     # ── end step loop


  # ── WRITE SUMMARY CSVs ────────────────────────────────────────────────────
  # summary_patch_counts_and_area.csv  — one row per combo (wide format)
  # summary_patch_areas_long.csv       — one row per patch (long format)
  summary_tbl <- dplyr::bind_rows(all_summaries)
  areas_tbl   <- dplyr::bind_rows(all_areas)

  log_msg("Writing summaries to %s", out_dir)
  readr::write_csv(summary_tbl,
                   file.path(out_dir, "summary_patch_counts_and_area.csv"))
  if (nrow(areas_tbl))
    readr::write_csv(areas_tbl,
                     file.path(out_dir, "summary_patch_areas_long.csv"))

  log_msg("Done. Shapefiles in: %s", shp_dir)

  # Return all results invisibly (caller can assign to a variable or ignore)
  invisible(list(
    summary    = summary_tbl,
    areas      = areas_tbl,
    shp_dir    = shp_dir,
    shapefiles = written_shps,
    out_dir    = out_dir
  ))
}


# =============================================================================
# RUN
# =============================================================================
# Trigger the full detection for the currently configured river reach.
# Results are saved to disk AND returned as a list (assigned to `res`).
res <- detect_cold_patches(
  ras_path          = ras_path,
  line_path         = line_path,
  steps             = steps_m,
  buffers_px        = buffers_px,
  deltas            = deltas_C,
  min_patch_area_m2 = min_area_m2,
  handle_small      = handle_small,
  out_dir           = out_dir
)

# Quick test block (wrapped in if(F) so it never runs accidentally).
# Uncomment and set if(T) to test a single parameter combination quickly.
if (F) {
  steps_m      <- c(100)
  buffers_px   <- c(10)
  deltas_C     <- c(2, 3)
  min_area_m2  <- 2
  handle_small <- "drop"
  out_dir      <- "output/cold_patches_out"
}

# Save the full results object to disk using qs (faster than RDS, smaller file)
qsave(res, file = paste0(out_dir, "/Results_CWP_detection.qs"))


# =============================================================================
# RELOAD RESULTS (optional)
# =============================================================================
# The if(F) block below is a manual reload path — useful when re-running only
# the plotting section without re-detecting patches.
if (F) {
  pr  <- 982   # plot counter (used for naming output files sequentially)
  res <- qread(file = paste0(out_dir, "/Results_CWP_detection.qs"))
}


# =============================================================================
# SENSITIVITY PLOTS
# =============================================================================
# These plots visualise how the CWP detection results change across the
# parameter grid. The goal is to identify which parameter drives the most
# variability — and to find a "stable zone" where results are robust.

plot_dir <- file.path(out_dir, "plots_sensitivity")
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

# Prepare the summary table: cast columns to correct types and add factor versions
sum0 <- res$summary |>
  dplyr::filter(!is.na(n_patches)) |>
  dplyr::mutate(
    step_m    = as.integer(step_m),
    buffer_px = as.integer(buffer_px),
    delta     = as.numeric(delta),
    step_f    = factor(step_m),    # factor version for ggplot colour/group mapping
    buf_f     = factor(buffer_px),
    dlt_f     = factor(delta)
  )

# ── HELPER: save a plot as both TIFF and PDF ───────────────────────────────
save_tiff_pdf <- function(plot, path_no_ext, width_in, height_in, dpi = 300) {
  grDevices::tiff(paste0(path_no_ext, ".tiff"),
                  width = width_in, height = height_in, units = "in",
                  res = dpi, compression = "lzw")
  print(plot)
  grDevices::dev.off()

  grDevices::pdf(paste0(path_no_ext, ".pdf"),
                 width = width_in, height = height_in, useDingbats = FALSE)
  print(plot)
  grDevices::dev.off()
}


# ── PLOT 1: Variance explained by each parameter (ANOVA approach) ──────────
# Run a simple linear model for each output metric (n_patches, total_area, etc.)
# using the three parameters as factors. The ANOVA sum-of-squares tells us
# what proportion of the total variance each parameter accounts for.
#
# This answers: "Does step, buffer, or ΔT matter most for the results?"
# A parameter with ~90% variance explained is the dominant driver.
# A parameter with <5% can be fixed at any reasonable value without much effect.

sum0 <- res$summary

metrics <- c("n_patches", "total_area_m2", "mean_area_m2", "median_area_m2")

# Auto-detect study area name from the working directory path
# Expects: ".../automatisation/<study_area>/..."
wd       <- getwd()
parts    <- strsplit(wd, .Platform$file.sep)[[1]]
auto_idx <- which(parts == "automatisation")
study_area <- if (length(auto_idx) == 1 && length(parts) > auto_idx)
                parts[auto_idx + 1]
              else
                "unbekannt"

study_area_file  <- tolower(study_area)
study_area_title <- paste0(toupper(substring(study_area, 1, 1)),
                            substring(study_area, 2))
message("Erkanntes Studiengebiet: ", study_area_title)

# Extended save helper that handles ggplot, patchwork, and grid grobs
save_plot_pdf_tiff <- function(plot_obj, file_base, width = 10, height = 4,
                               dpi = 300) {
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

  message("Gespeichert: ", basename(pdf_file), " und ", basename(tiff_file))
}

# Only include parameters that actually vary in this run
# (if only one buffer value was tested, it can't explain any variance)
params_all  <- c("step_m", "buffer_px", "delta")
params_vary <- params_all[sapply(sum0[params_all],
                                 function(x) length(unique(x)) > 1)]
message("Verwendete Parameter: ", paste(params_vary, collapse = ", "))

# Fit one LM per metric, extract ANOVA SS, compute fraction of total SS
fx_list <- lapply(metrics, function(m) {
  f   <- as.formula(paste(m, "~",
                          paste(sprintf("factor(%s)", params_vary),
                                collapse = " + ")))
  fit <- lm(f, data = sum0)
  ss  <- anova(fit)
  den <- sum(ss[seq_along(params_vary), "Sum Sq"])   # total SS across all parameters
  data.frame(
    param    = params_vary,
    var_expl = if (den > 0)
                 ss[seq_along(params_vary), "Sum Sq"] / den
               else NA_real_,
    metric   = m,
    row.names = NULL
  )
})

# Combine and relabel metrics into German (project convention)
fx <- dplyr::bind_rows(fx_list) |>
  dplyr::mutate(
    metric = dplyr::recode(metric,
      n_patches      = "Anzahl der CWP",
      total_area_m2  = "Gesamtfläche",
      mean_area_m2   = "Mittlere Fläche",
      median_area_m2 = "Medianfläche"
    )
  )

# Bar chart: one facet per metric, bars sorted by variance explained
p_fx <- ggplot2::ggplot(fx, ggplot2::aes(x = reorder(param, -var_expl),
                                          y = var_expl)) +
  ggplot2::geom_col(width = 0.7) +
  ggplot2::scale_y_continuous(labels = scales::percent) +
  ggplot2::facet_wrap(~ metric, nrow = 1, scales = "free_y") +
  ggplot2::labs(x = "Parameter",
                y = "% erklärte Varianz (Haupteffekte)",
                title = "") +
  ggplot2::theme_minimal(base_size = 12)

out_dir_plots <- file.path(out_dir, "plots_sensitivity")
dir.create(out_dir_plots, showWarnings = FALSE, recursive = TRUE)

save_tiff_pdf(p_fx,
              file.path(out_dir_plots,
                        paste0("sensitivity_variance_explained_all_metrics_",
                               pr, study_area_file)),
              15, 4)
pr <- pr + 1


# =============================================================================
# PLOT 2–5: Line plots of patch metrics vs ΔT, coloured by step length
# =============================================================================
# These four plots show directly how each output metric changes as ΔT increases,
# with one line per step value (and linetype per buffer if multiple buffers tested).
#
# Expected pattern:
#   n_patches and total_area should DECREASE as ΔT increases
#   (stricter threshold → fewer qualifying patches).
#   If lines for different step values converge, the step parameter doesn't
#   matter much. If they diverge wildly, step is a sensitive choice.

# Helper for optional side-by-side combination using patchwork or gridExtra
.combine_side_by_side <- function(..., title = NULL) {
  plots <- list(...)
  if (requireNamespace("patchwork", quietly = TRUE)) {
    comb <- Reduce(`+`, plots) +
      patchwork::plot_layout(ncol = length(plots), guides = "collect")
    comb <- comb & ggplot2::theme(legend.position = "bottom")
    if (!is.null(title))
      comb <- comb + patchwork::plot_annotation(title = title)
    return(comb)
  } else if (requireNamespace("gridExtra", quietly = TRUE)) {
    grob <- do.call(gridExtra::arrangeGrob, c(plots, list(ncol = length(plots))))
    grid::grid.newpage(); grid::grid.draw(grob)
    return(grob)
  } else {
    return(plots)
  }
}

# Prepare factor levels for consistent ordering in legends
steps_levels  <- sort(unique(sum0$step_m))
deltas_levels <- sort(unique(sum0$delta))

# Check whether buffer varies (if only one value, it adds nothing to the plot)
has_buffer <- "buffer_px" %in% names(sum0) &&
              length(unique(sum0$buffer_px)) > 1

if (has_buffer) {
  buffers_levels <- sort(unique(sum0$buffer_px))
  sum0 <- sum0 |> dplyr::mutate(
    step_f   = factor(step_m,    levels = steps_levels),
    buffer_f = factor(buffer_px, levels = buffers_levels)
  )
} else {
  sum0 <- sum0 |> dplyr::mutate(step_f = factor(step_m, levels = steps_levels))
}

# Shared aesthetic mapping helper: returns the correct aes() depending on
# whether buffer varies (uses linetype) or not (single linetype)
make_aes <- function(y_var) {
  if (has_buffer)
    ggplot2::aes(x = delta, y = !!ggplot2::sym(y_var),
                 color = step_f, linetype = buffer_f,
                 group = interaction(step_f, buffer_f))
  else
    ggplot2::aes(x = delta, y = !!ggplot2::sym(y_var),
                 color = step_f, group = step_f)
}

# Plot 1: number of detected patches vs ΔT
p_count <- ggplot2::ggplot(sum0, make_aes("n_patches")) +
  ggplot2::geom_line(size = 0.9, na.rm = TRUE) +
  ggplot2::geom_point(size = 1.2, na.rm = TRUE) +
  ggplot2::scale_x_continuous(breaks = deltas_levels) +
  ggplot2::labs(x      = expression(Delta * T ~ "[°C]"),
                y      = "Anzahl der CWP",
                color  = "Step [m]",
                linetype = if (has_buffer) "Puffer (Pixel)" else NULL,
                title  = " ") +
  ggplot2::theme_minimal(base_size = 12)

# Plot 2: total patch area vs ΔT
p_total_area <- ggplot2::ggplot(
    dplyr::filter(sum0, !is.na(total_area_m2)),
    make_aes("total_area_m2")
  ) +
  ggplot2::geom_line(size = 0.9, na.rm = TRUE) +
  ggplot2::geom_point(size = 1.2, na.rm = TRUE) +
  ggplot2::scale_x_continuous(breaks = deltas_levels) +
  ggplot2::labs(x = expression(Delta * T ~ "[°C]"),
                y = "Gesamtfläche der CWP [m²]",
                color = "Step [m]",
                linetype = if (has_buffer) "Puffer (Pixel)" else NULL,
                title = " ") +
  ggplot2::theme_minimal(base_size = 12)

# Plot 3: mean patch area vs ΔT
p_mean_area <- ggplot2::ggplot(
    dplyr::filter(sum0, !is.na(mean_area_m2)),
    make_aes("mean_area_m2")
  ) +
  ggplot2::geom_line(size = 0.9, na.rm = TRUE) +
  ggplot2::geom_point(size = 1.2, na.rm = TRUE) +
  ggplot2::scale_x_continuous(breaks = deltas_levels) +
  ggplot2::labs(x = expression(Delta * T ~ "[°C]"),
                y = "Mittlere CWP-Flächengrösse [m²]",
                color = "Step [m]",
                linetype = if (has_buffer) "Puffer (Pixel)" else NULL,
                title = " ") +
  ggplot2::theme_minimal(base_size = 12)

# Plot 4: median patch area vs ΔT
p_med_area <- ggplot2::ggplot(
    dplyr::filter(sum0, !is.na(median_area_m2)),
    make_aes("median_area_m2")
  ) +
  ggplot2::geom_line(size = 0.9, na.rm = TRUE) +
  ggplot2::geom_point(size = 1.2, na.rm = TRUE) +
  ggplot2::scale_x_continuous(breaks = deltas_levels) +
  ggplot2::labs(x = expression(Delta * T ~ "[°C]"),
                y = "Median der CWP-Flächengrösse [m²]",
                color = "Step [m]",
                linetype = if (has_buffer) "Puffer (Pixel)" else NULL,
                title = " ") +
  ggplot2::theme_minimal(base_size = 12)


# ── COMBINE AND EXPORT ────────────────────────────────────────────────────
# Extract a shared legend from p_count, then remove individual legends
# so the combined 4-panel figure has one clean legend on the right.

if (!requireNamespace("cowplot", quietly = TRUE)) install.packages("cowplot")

legend_plot <- cowplot::get_legend(
  p_count + ggplot2::theme(
    legend.position = "right",
    legend.title    = ggplot2::element_text(size = 10),
    legend.text     = ggplot2::element_text(size = 9)
  )
)

# Strip legends from individual panels
p_count      <- p_count      + ggplot2::theme(legend.position = "none")
p_total_area <- p_total_area + ggplot2::theme(legend.position = "none")
p_mean_area  <- p_mean_area  + ggplot2::theme(legend.position = "none")
p_med_area   <- p_med_area   + ggplot2::theme(legend.position = "none")

# Arrange the four panels side by side
plots_row <- cowplot::plot_grid(
  p_count, p_total_area, p_med_area, p_mean_area,
  nrow = 1, align = "v"
)

# Attach the shared legend on the right (12% of total width)
combined_4 <- cowplot::plot_grid(
  plots_row, legend_plot,
  ncol = 2, rel_widths = c(1, 0.12)
)

# Save the combined figure as both PDF and TIFF
out_base_4 <- file.path(
  plot_dir,
  paste0("linien_anzahl_gesamt_mittel_median_nach_delta_", pr, study_area_file)
)
save_plot_pdf_tiff(combined_4, out_base_4, width = 20, height = 5)
pr <- pr + 1
