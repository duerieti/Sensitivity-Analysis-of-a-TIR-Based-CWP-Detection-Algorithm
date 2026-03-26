#!/bin/bash
#
#SBATCH --job-name=run_new_function
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=04:00:00
#SBATCH --partition=earth-1
#SBATCH --mem=150G

module load USS/2022
module load gcc/9.4.0-pe5.34
module load lsfm-init-miniconda/1.0.0

WORKDIR=/data/scratch/duerieti/$SLURM_JOB_ID
mkdir -p $WORKDIR

trap "rm -rf $WORKDIR" EXIT

cp -r /net/beegfs/earth/scratch/duerieti/BSc_project/data $WORKDIR/


export WORKDIR

conda activate R_env

Rscript function_extract.R


