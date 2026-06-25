#!/bin/bash
#SBATCH --job-name=run_new_function
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=04:00:00
#SBATCH --partition=earth-3
#SBATCH --mem=10G
module load USS/2022
module load gcc/9.4.0-pe5.34
module load lsfm-init-miniconda/1.0.0
conda activate R_env
     gdal raster calc -i "A=../original_data/thermal_rasters_FINAL/mean_v2emme.tif" \
     --calc "rint(A*10)/10" \
     -o mean_v2emme_rounded.tif --overwrite --ot Float32 --nodata -9999 --propagate-nodata \
     --co TILED=YES --co BLOCKXSIZE=512 --co BLOCKYSIZE=512 --co BIGTIFF=YES --co COMPRESS=DEFLATE --co PREDICTOR=2
