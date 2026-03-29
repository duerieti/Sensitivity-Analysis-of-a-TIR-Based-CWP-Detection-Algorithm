# Step 1: rasterize + reclassify → write tmean to disk
gdal raster pipeline \
  ! read ../data/test_dir/zone_r_big.tif \
  ! reclassify -m "$(cat ../data/test_dir/median_lookup_gdal.txt)" \
  ! write ../data/test_dir/tmean.tif --overwrite

# Step 2: calc using both r.tif and tmean.tif
gdal raster calc \
  -i "A=../data/thermal_rasters_FINAL/mean_v01emme.tif" \
  -i "B=../data/test_dir/tmean.tif" \
  --calc "((A - B) <= -0.07) ? (A - B) : NaN" \
  --no-check-extent \
  -o ../data/test_dir/binary_out.tif --overwrite
