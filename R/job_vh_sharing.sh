#! /bin/bash
#SBATCH --time=02:00:00
#SBATCH --account=ctb-tpoisot
#SBATCH --array=1-127
#SBATCH --nodes=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=100G

module load r/4.1.2

export R_LIBS=~/local/R_libs/
Rscript R/viral_host_sharing_corrected.R
