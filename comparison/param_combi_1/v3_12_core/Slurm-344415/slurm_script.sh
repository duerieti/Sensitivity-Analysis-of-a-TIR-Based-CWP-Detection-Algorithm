#!/bin/bash
#BATCH --job-name=morris_screening
#SBATCH --time=0-08:00:00
#SBATCH --partition=earth-3
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=100G

module load USS/2022
module load gcc/9.4.0-pe5.34 parallel/20221022-pe5.34
module load lsfm-init-miniconda/1.0.0

# Use hardcoded conda path instead of dynamic conda info --base
source /cfs/software/uss/2022/spack/linux-rocky8-x86_64/gcc-9.4.0/miniconda3-4.12.0-7riiqd3uthjsmxqkabjrfhfxhyh5epcl/etc/profile.d/conda.sh
conda activate R_env

# Force R to use module-loaded GCC 9.4
mkdir -p ~/.R
echo "CC=$(which gcc)" > ~/.R/Makevars
echo "CXX=$(which g++)" >> ~/.R/Makevars
echo "CXX17=$(which g++)" >> ~/.R/Makevars
echo "FC=$(which gfortran)" >> ~/.R/Makevars

Rscript run_model.R
