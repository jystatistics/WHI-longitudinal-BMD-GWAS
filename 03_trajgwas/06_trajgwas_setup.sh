#!/bin/bash
#SBATCH --job-name=trajgwas_setup
#SBATCH --output=trajgwas_setup_%j.out
#SBATCH --time=06:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --account=YOUR_ACCOUNT
#SBATCH --mail-type=END,FAIL
# 2026-11-24. One-time setup for the unified TrajGWAS (run after the final merge is built).
# 2026-11-25: uses WHI_v4_r08 = R2>=0.8 variants; 912 PAGE/SHARE duplicates (04b) and 134 GARNET-WHISP/PAGE
#             cross-ID duplicates (04e) and 152 Z2-based duplicates (04f) removed -> WHI_v5_r08
#  1. cohort_by_iid.txt : IID -> genotype sub-study (GARNET, WHIMS, WHISP, PAGE, SHARE), used as a covariate
#  2. geno/chrN         : WHI_v5_r08 (R2>=0.8; three duplicate passes removed; 04b, 04e, 04f) split by chromosome, so each array task scans one chromosome
set -euo pipefail
source "${SLURM_SUBMIT_DIR:-$PWD}/config.sh"   # site paths (see config.sh)
B=${WHI_DIR}/10_Rebuild_v2
cd $B
# cohort_by_fidiid.txt is written by 04e_dedup_cross_id.sh (FID IID cohort for every per-cohort file)
awk 'NR==FNR{c[$1" "$2]=$3; next} {print $2, c[$1" "$2]}' cohort_by_fidiid.txt WHI_v5_r08.fam > cohort_by_iid.txt
[ -z "$(awk 'NF<2' cohort_by_iid.txt)" ] || { echo "women without a cohort label"; exit 1; }
[ -z "$(cut -d' ' -f1 cohort_by_iid.txt | sort | uniq -d)" ] || { echo "duplicated IIDs in cohort_by_iid.txt"; exit 1; }
echo "cohort_by_iid.txt: $(wc -l < cohort_by_iid.txt) women"; cut -d' ' -f2 cohort_by_iid.txt | sort | uniq -c
mkdir -p geno trajgwas
for i in $(seq 1 22); do
  plink --memory 56000 --bfile WHI_v5_r08 --chr $i --make-bed --out geno/chr$i > /dev/null
  awk '{$1=$2; print}' geno/chr$i.fam > geno/chr$i.fam.tmp && mv geno/chr$i.fam.tmp geno/chr$i.fam   # FID = IID for TrajGWAS
  echo "chr$i: $(wc -l < geno/chr$i.bim) variants"
done
