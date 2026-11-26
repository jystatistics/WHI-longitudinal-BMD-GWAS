#!/bin/bash
#SBATCH --job-name=trajgwas_v2
#SBATCH --output=trajgwas_v2_%A_%a.out
#SBATCH --time=48:00:00
#SBATCH --mem=48G
#SBATCH --cpus-per-task=4
#SBATCH --account=YOUR_ACCOUNT
#SBATCH --mail-type=FAIL
#SBATCH --array=1-22
# 2026-11-24. One array task per chromosome.
# Submit: sbatch 07_trajgwas_run.sh <group> <outcome> <phenotype csv>
#   e.g.  sbatch 07_trajgwas_run.sh black hipbmd ${WHI_DIR}/01_Black/pheno_black_hipbmd_data.csv
set -euo pipefail
source "${SLURM_SUBMIT_DIR:-$PWD}/config.sh"   # site paths (see config.sh)
export JULIA_NUM_THREADS=$SLURM_CPUS_PER_TASK
export COLUMNS=400
# same julia as the earlier TrajGWAS runs (whatever is on PATH); printed so the log shows which one ran
echo "julia: $(command -v julia || echo NOT FOUND)  $(julia --version 2>/dev/null)"
cd ${WHI_DIR}/10_Rebuild_v2
julia --threads $SLURM_CPUS_PER_TASK 07_trajgwas.jl "$1" "$2" "$SLURM_ARRAY_TASK_ID" "$3"
