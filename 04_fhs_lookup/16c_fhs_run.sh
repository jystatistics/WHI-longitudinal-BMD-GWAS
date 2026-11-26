#!/bin/bash
#SBATCH --job-name=fhs_traj
#SBATCH --output=fhs_traj_%A_%a.out
#SBATCH --time=04:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=2
#SBATCH --account=YOUR_ACCOUNT
#SBATCH --mail-type=FAIL
#SBATCH --array=1-3
set -euo pipefail
source "${SLURM_SUBMIT_DIR:-$PWD}/config.sh"   # site paths (see config.sh)
export COLUMNS=400
F=${FHS_DIR}; O=$F/10_Lookup_v2
case $SLURM_ARRAY_TASK_ID in
  1) o=fnbmd; d=01_FNBMD; y=femoral_neck_bmd;;
  2) o=trbmd; d=02_TRBMD; y=trochanteric_bmd;;
  3) o=lsbmd; d=03_LSBMD; y=lumbar_spine_l24_bmd;;
esac
echo "julia: $(command -v julia || echo NOT FOUND)  $(julia --version 2>/dev/null)"
cd ${WHI_DIR}/10_Rebuild_v2
julia 16_fhs_trajgwas.jl fhs_$o $y $F/$d/fhs_pheno_female.csv $O/lookup_$o $O/pca_$o.eigenvec $O
