library(terra)
library(sf)
library(dplyr)
library(lwgeom)
library(tictoc)
library(parallel)


# ── Per-slab worker (defined at top level so mclapply can find it) ────────────
# Each call receives one slab + its buffer and does the full detection pipeline.
# The rounded raster is NOT passed as a SpatRaster (C++ pointer, unsafe to fork);
# instead every worker re-reads from the same temp file on disk.
.process_one_slab <- function(args) {

  suppressPackageStartupMessages({
    library(terra)
    library(sf)
    library(dplyr)
  })

  # Reconstruct sf objects from the serialised sfc geometries
  slab_sf <- sf::st_sf(geometry = args$slab_geom)
  buf_sf  <- sf::st_sf(geometry = args$buf_geom)

  r <- terra::rast(args$raster_path)
  names(r) <- "T"

  # Crop to extents (masking slab; buffer kept unmasked for median extraction)
  r_slab <- terra::crop(r, terra::vect(slab_sf), mask = TRUE)
  r_buf  <- terra::crop(r, terra::vect(buf_sf),  mask = FALSE)

  # Reference temperature: median of all pixels inside the buffer strip
  ref_vals <- terra::extract(r_buf, terra::vect(buf_sf))
  ref_temp <- median(ref_vals[[2]], na.rm = TRUE)
  if (is.na(ref_temp)) return(NULL)

  # Flag pixels colder than ref_temp - delta_C
  binary <- (r_slab - ref_temp) <= (-args$delta_C)
  binary[binary == 0] <- NA

  if (all(is.na(terra::values(binary)))) return(NULL)

  patches_v <- terra::as.polygons(binary, dissolve = TRUE, eight = args$connect_d)
  if (length(patches_v) == 0) return(NULL)

  # Polygonise — optional diagonal-pixel merging via small expand/shrink
  eps <- if (args$connect_d) args$cell_m * 0.1 else 0
  patches_sf <- tryCatch({
    patches_v |>
      sf::st_as_sf()           |>
      sf::st_cast("POLYGON")   |>
      sf::st_buffer(eps)       |>
      sf::st_union()           |>
      sf::st_buffer(-eps)      |>
      sf::st_cast("POLYGON")   |>
      sf::st_as_sf()
  }, error = function(e) NULL)

  if (is.null(patches_sf) || nrow(patches_sf) == 0) return(NULL)

  # Area filter; assign local row ID for the stat join below
  patches_sf <- patches_sf |>
    dplyr::select(-dplyr::any_of("T")) |>
    dplyr::filter(as.numeric(sf::st_area(.)) >= args$min_area) |>
    dplyr::mutate(ID = dplyr::row_number())

  if (nrow(patches_sf) == 0) return(NULL)

  # Temperature statistics per patch
  stats <- terra::extract(r_slab, terra::vect(patches_sf)) |>
    dplyr::group_by(ID) |>
    dplyr::summarise(
      mean_temp   = mean(T,   na.rm = TRUE),
      min_temp    = min(T,    na.rm = TRUE),
      max_temp    = max(T,    na.rm = TRUE),
      median_temp = median(T, na.rm = TRUE),
      .groups = "drop"
    )

  patches_sf |>
    dplyr::left_join(stats, by = "ID") |>
    dplyr::mutate(
      slab_id  = args$j,
      ref_temp = ref_temp,
      deltaT   = ref_temp - mean_temp
    ) |>
    dplyr::filter(deltaT >= args$delta_C) |>
    dplyr::select(-ID)
}


detect_cwp_parallel <- function(
    ras_path,
    line_path,
    step_m            = 500,
    buffer_px         = 3,
    delta_C           = 1.0,
    min_patch_area_m2 = 2,
    round_to          = 0.1,
    slab_halfwidth_m  = 60,
    connect_diagonals = TRUE,
    n_cores           = max(1L, parallel::detectCores() - 1L)
) {

  timings <- numeric()
  tic_fn  <- function() proc.time()[["elapsed"]]

  # ── 0. READ INPUTS ──────────────────────────────────────────────────────────
  t0 <- tic_fn()

  r    <- terra::rast(ras_path)
  line <- sf::st_read(line_path, quiet = TRUE) |> sf::st_zm(TRUE, "ZM")

  stopifnot(terra::nlyr(r) == 1)
  names(r) <- "T"
  if (terra::is.lonlat(r)) stop("Raster must be in a metric CRS.")

  cell_m <- mean(terra::res(r))

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

  rfactor <- if (!is.null(round_to) && round_to > 0) 1 / round_to else NA_real_

  # Round raster once, write to temp file — workers re-read from disk.
  # SpatRasters are backed by C++ pointers and cannot be forked safely.
  r_rounded  <- round(r * rfactor) / rfactor
  tmp_raster <- tempfile(fileext = ".tif")
  terra::writeRaster(r_rounded, tmp_raster, overwrite = TRUE)
  on.exit(unlink(tmp_raster), add = TRUE)
  rm(r, r_rounded)

  timings[["read_inputs"]] <- tic_fn() - t0

  # ── 1. BUILD SLABS ──────────────────────────────────────────────────────────
  t0 <- tic_fn()

  L   <- as.numeric(sf::st_length(line))
  pos <- sort(unique(c(0, seq(0, 1, by = step_m / L), 1)))
  npt <- length(pos)

  if (npt < 2) stop("Too few stations for step = ", step_m)

  from <- pos[seq_len(npt - 1)]
  to   <- pos[seq_len(npt - 1) + 1]
  k    <- length(from)

  g1        <- sf::st_geometry(line)[1]
  subs_list <- vector("list", k)
  for (j in seq_len(k))
    subs_list[[j]] <- lwgeom::st_linesubstring(g1, from[j], to[j])
  subs <- do.call(c, subs_list)

  slabs     <- sf::st_buffer(subs, dist = slab_halfwidth_m,
                             endCapStyle = "FLAT", joinStyle = "MITRE", mitreLimit = 2)
  slabs_buf <- sf::st_buffer(subs, dist = cell_m * buffer_px,
                             endCapStyle = "FLAT", joinStyle = "MITRE", mitreLimit = 2)

  timings[["build_slabs"]] <- tic_fn() - t0

  # ── 2. DISPATCH ONE SLAB PER CORE ───────────────────────────────────────────
  t0 <- tic_fn()

  # Build plain-list argument bundles — fully serialisable, no terra handles.
  # sfc objects are ordinary R lists with a CRS attribute; they fork cleanly.
  slab_args <- lapply(seq_len(k), function(j) {
    list(
      j           = j,
      slab_geom   = slabs[j],
      buf_geom    = slabs_buf[j],
      raster_path = tmp_raster,
      delta_C     = delta_C,
      min_area    = min_patch_area_m2,
      cell_m      = cell_m,
      connect_d   = connect_diagonals
    )
  })

  n_cores <- min(n_cores, k)

  if (.Platform$OS.type == "windows") {
    # PSOCK cluster: workers are fresh R sessions — need to export the helper
    cl <- parallel::makeCluster(n_cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterExport(cl, varlist = ".process_one_slab", envir = environment())
    results <- parallel::parLapply(cl, slab_args, .process_one_slab)
  } else {
    # Fork-based: workers inherit the process image, helper already available
    results <- parallel::mclapply(slab_args, .process_one_slab,
                                  mc.cores      = n_cores,
                                  mc.preschedule = TRUE)
  }

  timings[["process_slabs_parallel"]] <- tic_fn() - t0

  # ── 3. COMBINE RESULTS ──────────────────────────────────────────────────────
  t0 <- tic_fn()

  non_null <- Filter(Negate(is.null), results)

  if (length(non_null) == 0) {
    final_patches <- sf::st_sf(geometry = sf::st_sfc(crs = sf::st_crs(line)))
  } else {
    final_patches <- do.call(rbind, non_null) |>
      dplyr::mutate(ID = dplyr::row_number())
  }

  timings[["combine_results"]] <- tic_fn() - t0
  timings[["TOTAL"]]           <- sum(timings)

  list(patches = final_patches, timings = timings)
}


path_raster <- file.path("../data/thermal_rasters_FINAL/mean_v01emme.tif")
path_line   <- file.path("../data/Centerlines_FINAL/Emme_V01.shp")

detect_cwp_parallel(path_raster, path_line)
