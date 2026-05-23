import numpy as np
import geopandas as gpd
import rasterio
import rasterio.features
import rasterio.transform
from rasterio.features import shapes
from shapely.geometry import shape, mapping
from shapely.ops import unary_union, linemerge
import shapely
import time
import warnings
from scipy.ndimage import label, binary_dilation
import dask.array as da
import dask_geopandas  # optional, for parallel vector ops

# ─────────────────────────────────────────────────────────────────────────────
def detect_cwp_single(
    ras_path,
    line_path,
    step_m            = 500,
    buffer_px         = 3,
    delta_C           = 1.0,
    min_patch_area_m2 = 2,
    round_to          = 0.1,
    slab_halfwidth_m  = 60,
    connect_diagonals = True,
    use_dask          = True,       # toggle parallel raster ops
    dask_chunks       = 512,        # chunk size for dask array
):
    timings = {}
    rfactor = (1 / round_to) if round_to and round_to > 0 else None

    # ── 0. READ INPUTS ────────────────────────────────────────────────────────
    t0 = time.perf_counter()

    src       = rasterio.open(ras_path)
    transform = src.transform
    crs       = src.crs
    nrows, ncols = src.height, src.width
    nodata    = src.nodata
    cell_m    = (abs(transform.a) + abs(transform.e)) / 2

    # read as dask array for parallel ops, or plain numpy
    if use_dask:
        import dask.array as da
        rdata = da.from_array(
            src.read(1).astype("float32"),
            chunks=(dask_chunks, dask_chunks)
        )
    else:
        rdata = src.read(1).astype("float32")

    if nodata is not None:
        rdata = da.where(rdata == nodata, np.nan, rdata) \
                if use_dask else \
                np.where(rdata == nodata, np.nan, rdata)

    gdf_line = gpd.read_file(line_path).to_crs(crs)
    merged   = linemerge(unary_union(gdf_line.geometry))
    if merged.geom_type == "MultiLineString":
        merged = linemerge(merged)
    assert merged.geom_type == "LineString"

    timings["read_inputs"] = time.perf_counter() - t0

    # ── 1. BUILD SLABS ────────────────────────────────────────────────────────
    t0 = time.perf_counter()

    L        = merged.length
    positions = np.unique(np.concatenate([[0], np.arange(0, L, step_m), [L]]))
    k        = len(positions) - 1
    buffer_m = buffer_px * cell_m

    slabs, slabs_buffer = [], []
    for j in range(k):
        seg  = shapely.ops.substring(merged, positions[j], positions[j + 1])
        slab = seg.buffer(slab_halfwidth_m, cap_style=2,
                          join_style=3, mitre_limit=2)
        buf  = seg.buffer(buffer_m, cap_style=2,
                          join_style=3, mitre_limit=2)
        slabs.append(slab)
        slabs_buffer.append(buf)

    # GeoDataFrames with slab_id (1-indexed to match rasterize background=0)
    gdf_slabs  = gpd.GeoDataFrame(
        {"slab_id": np.arange(1, k + 1)},
        geometry=slabs, crs=crs
    )
    gdf_buffer = gpd.GeoDataFrame(
        {"slab_id": np.arange(1, k + 1)},
        geometry=slabs_buffer, crs=crs
    )

    timings["build_slabs"] = time.perf_counter() - t0

    # ── 2. RASTERIZE SLABS → zone raster ─────────────────────────────────────
    # equivalent of terra::rasterize(slaps_sf_vect, r, field="slap_id")
    t0 = time.perf_counter()

    zone_arr = rasterio.features.rasterize(
        shapes   = zip(gdf_slabs.geometry.apply(mapping),
                       gdf_slabs["slab_id"]),
        out_shape = (nrows, ncols),
        transform = transform,
        fill      = 0,          # background = 0
        dtype     = "int32"
    )

    timings["rasterize_slabs"] = time.perf_counter() - t0

    # ── 3. ROUND RASTER ───────────────────────────────────────────────────────
    # equivalent of r_rounded <- round(r * rfactor) / rfactor
    # this is where dask parallelises transparently
    t0 = time.perf_counter()

    if rfactor is not None:
        r_rounded = da.round(rdata * rfactor) / rfactor \
                    if use_dask else \
                    np.round(rdata * rfactor) / rfactor
    else:
        r_rounded = rdata

    # materialise if dask — needed for geometry_mask ops below
    r_np = r_rounded.compute() if use_dask else r_rounded

    timings["round_raster"] = time.perf_counter() - t0

    # ── 4. COMPUTE MEDIAN PER SLAB CORRIDOR ───────────────────────────────────
    # equivalent of terra::extract(r_rounded, slabs_buffer_vect, fun="median")
    t0 = time.perf_counter()

    # vectorised: build one mask per buffer polygon, compute median
    # parallelisable via joblib if k is large
    from joblib import Parallel, delayed

    def _slab_median(geom, r_np, transform, nrows, ncols, rfactor):
        mask = rasterio.features.geometry_mask(
            [mapping(geom)],
            out_shape = (nrows, ncols),
            transform = transform,
            invert    = True
        )
        vals = r_np[mask & ~np.isnan(r_np)]
        if not len(vals):
            return np.nan
        m = float(np.median(vals))
        return np.round(m * rfactor) / rfactor if rfactor else m

    # parallel across slabs — this is embarrassingly parallel
    Tmed = Parallel(n_jobs=-1)(          # n_jobs=-1 → all cores
        delayed(_slab_median)(geom, r_np, transform, nrows, ncols, rfactor)
        for geom in gdf_buffer.geometry
    )
    Tmed = np.array(Tmed)

    timings["compute_ref_temps"] = time.perf_counter() - t0

    # ── 5. BURN MEDIANS INTO ZONE RASTER ─────────────────────────────────────
    # equivalent of terra::classify(zone_r_big, mean_raster)
    # lookup table: slab_id → median temperature
    t0 = time.perf_counter()

    # build lookup array: index = slab_id, value = Tmed
    lut       = np.full(k + 1, np.nan, dtype="float32")
    lut[1:]   = Tmed.astype("float32")          # slab_id is 1-indexed

    # vectorised reclassify — pure numpy, very fast
    valid     = zone_arr > 0
    Tmean_arr = np.full((nrows, ncols), np.nan, dtype="float32")
    Tmean_arr[valid] = lut[zone_arr[valid]]

    timings["burn_medians"] = time.perf_counter() - t0

    # ── 6. FLAG COLD PIXELS ───────────────────────────────────────────────────
    # equivalent of flagged_pixels <- r_rounded - Tmean; binary <- flagged <= -delta_C
    # dask parallelises this if use_dask=True
    t0 = time.perf_counter()

    if use_dask:
        r_da     = da.from_array(r_np,      chunks=(dask_chunks, dask_chunks))
        Tmean_da = da.from_array(Tmean_arr, chunks=(dask_chunks, dask_chunks))
        diff     = r_da - Tmean_da
        cold     = (diff <= -delta_C).compute()
    else:
        diff = r_np - Tmean_arr
        cold = diff <= -delta_C

    cold &= ~np.isnan(r_np) & ~np.isnan(Tmean_arr)

    timings["flag_pixels"] = time.perf_counter() - t0

    # ── 7. POLYGONIZE ─────────────────────────────────────────────────────────
    # equivalent of terra::as.polygons + st_buffer/union/cast dance
    t0 = time.perf_counter()

    struct = np.ones((3, 3)) if connect_diagonals else \
             np.array([[0,1,0],[1,1,1],[0,1,0]])

    labeled, _ = label(cold.astype("int32"), structure=struct)

    patch_geoms = [
        shape(geom)
        for geom, val in shapes(labeled.astype("int32"), transform=transform)
        if val > 0
    ]

    if connect_diagonals and patch_geoms:
        eps = cell_m * 0.1
        patch_geoms = [g.buffer(eps) for g in patch_geoms]
        merged_geom = unary_union(patch_geoms)
        patch_geoms = [g.buffer(-eps) for g in
                       (merged_geom.geoms
                        if merged_geom.geom_type == "MultiPolygon"
                        else [merged_geom])]

    gdf_patches = gpd.GeoDataFrame(geometry=patch_geoms, crs=crs) \
                  if patch_geoms else \
                  gpd.GeoDataFrame(geometry=[], crs=crs)

    timings["polygonize"] = time.perf_counter() - t0

    # ── 8. AREA FILTER ────────────────────────────────────────────────────────
    t0 = time.perf_counter()

    if len(gdf_patches):
        gdf_patches["area_m2"] = gdf_patches.geometry.area
        gdf_patches = gdf_patches[
            gdf_patches["area_m2"] >= min_patch_area_m2
        ].reset_index(drop=True)
        gdf_patches["ID"] = gdf_patches.index + 1

    timings["area_filter"] = time.perf_counter() - t0

    # ── 9. TEMPERATURE STATISTICS PER PATCH ──────────────────────────────────
    # equivalent of terra::extract + group_by + summarise
    t0 = time.perf_counter()

    if len(gdf_patches):

        def _patch_stats(geom, r_np, Tmean_arr, transform, nrows, ncols):
            mask = rasterio.features.geometry_mask(
                [mapping(geom)],
                out_shape = (nrows, ncols),
                transform = transform,
                invert    = True
            )
            vals      = r_np[mask & ~np.isnan(r_np)]
            tmean_vals = Tmean_arr[mask & ~np.isnan(Tmean_arr)]

            if not len(vals):
                return dict(mean_temp=np.nan, min_temp=np.nan,
                            max_temp=np.nan, median_temp=np.nan,
                            slab_mean_T=np.nan)
            return dict(
                mean_temp   = float(np.mean(vals)),
                min_temp    = float(np.min(vals)),
                max_temp    = float(np.max(vals)),
                median_temp = float(np.median(vals)),
                slab_mean_T = float(np.mean(tmean_vals))
                              if len(tmean_vals) else np.nan
            )

        # parallel across patches
        stats_list = Parallel(n_jobs=-1)(
            delayed(_patch_stats)(row.geometry, r_np, Tmean_arr,
                                  transform, nrows, ncols)
            for _, row in gdf_patches.iterrows()
        )

        stats_df = gpd.pd.DataFrame(stats_list)
        for col in stats_df.columns:
            gdf_patches[col] = stats_df[col].values

        gdf_patches["deltaT"] = gdf_patches["slab_mean_T"] - \
                                 gdf_patches["mean_temp"]
        gdf_patches["delta_thresh"] = delta_C

        gdf_patches = gdf_patches[
            gdf_patches["deltaT"] >= delta_C
        ].reset_index(drop=True)

    timings["patch_stats"] = time.perf_counter() - t0
    timings["TOTAL"]       = sum(timings.values())

    return {"patches": gdf_patches, "timings": timings}


# ── run ───────────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    result = detect_cwp_single(
        "../data/thermal_rasters_FINAL/mean_v01emme.tif",
        "../data/Centerlines_FINAL/Emme_V01.shp",
        use_dask    = True,
        dask_chunks = 512
    )
    print(result["timings"])
    print(result["patches"])
