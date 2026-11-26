#!/bin/bash
#SBATCH --job-name=traj_summary
#SBATCH --output=traj_summary_%j.out
#SBATCH --time=02:00:00
#SBATCH --mem=48G
#SBATCH --account=YOUR_ACCOUNT
#SBATCH --mail-type=END,FAIL
source "${SLURM_SUBMIT_DIR:-$PWD}/config.sh"   # site paths (see config.sh)
module load gcc/12.3.0 R/4.4.0
cd ${WHI_DIR}/10_Rebuild_v2
Rscript 09_trajgwas_summary.R
