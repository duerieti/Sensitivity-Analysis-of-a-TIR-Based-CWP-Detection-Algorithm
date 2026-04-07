#!/bin/bash
#
#SBATCH --job-name=run_new_function
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=04:00:00
#SBATCH --partition=earth-3
#SBATCH --mem=50G

module load USS/2022
module load gcc/9.4.0-pe5.34
module load lsfm-init-miniconda/1.0.0

conda activate R_env

Rscript function_terra.R


