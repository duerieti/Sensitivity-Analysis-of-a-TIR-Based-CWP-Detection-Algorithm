#!/bin/bash
#SBATCH --job-name=morris_screening
#SBATCH --time=0-02:00:00
#SBATCH --partition=earth-3
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=15
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=11G


module load USS/2022
module load gcc/9.4.0-pe5.34 parallel/20221022-pe5.34
module load lsfm-init-miniconda/1.0.0

source $(conda info --base)/etc/profile.d/conda.sh
conda activate R_env



mkdir -p results

parallel --colsep ' ' \
         --joblog results/joblog.txt \
         --eta \
         -j $SLURM_NTASKS \
         Rscript run_model.R {1} {2} {3} \
         :::: params_testing.txt
