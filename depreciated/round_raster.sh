#!/bin/bash
#
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


gdal raster reproject \
    --bbox 2623352.294469989836216,1183505.412900009884685,2638338.164469989836216,1197093.121100009884685 \
    --resolution 0.206889999565001,0.206889999798000 \
    -i ../data/r_v2obemme_rounded.tif \
    -o ../data/r_v2obemme_aligned.tif \
    --overwrite \
    --src-nodata -3.3999999521443642e+38 \
    --dst-nodata -9999
