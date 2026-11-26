#!/bin/bash
#SBATCH --job-name=share_merge
#SBATCH --output=share_merge_%x_%j.out
#SBATCH --error=share_merge_%x_%j.err
#SBATCH --time=06:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --account=YOUR_ACCOUNT
#SBATCH --mail-type=END,FAIL
# 2026-11-23. Merge the 22 per-chromosome PLINK sets from 03_SHARE_download_convert.sh
# into SHARE_${CG}_R2 and write the combined r2pass.txt that 04_merge_all_cohorts.sh reads.
# Submit right after the array job, so it starts only if all 22 tasks succeed:
#   J=$(sbatch --parsable --export=ALL,CG=c1 03_SHARE_download_convert.sh)
#   sbatch --dependency=afterok:${J} --export=ALL,CG=c1 03b_SHARE_merge_chr.sh
set -euo pipefail
source "${SLURM_SUBMIT_DIR:-$PWD}/config.sh"   # site paths (see config.sh)
CG=${CG:?set CG=c1 or c2}
OUT=${WHI_DBGAP_DIR}/03_NHLBI/phg000061.v2.p1.WHI_NHLBI.genotype-calls-matrixfmt.${CG}/subject_level_filtered_PLINK_sets/05_converted_R2
cd ${OUT}
for i in $(seq 1 22); do [ -s chr${i}.bed ] || { echo "missing chr${i}.bed"; exit 1; }; done
cat r2pass_chr{1..22}.txt > r2pass.txt
: > mergelist.txt; for i in $(seq 2 22); do echo "chr${i}" >> mergelist.txt; done
plink --memory 56000 --bfile chr1 --merge-list mergelist.txt --make-bed --out SHARE_${CG}_R2
echo "SHARE_${CG}_R2: $(wc -l < SHARE_${CG}_R2.bim) variants, $(wc -l < SHARE_${CG}_R2.fam) people"
# SPTBN1 rs4233949 (hg19 2:54659707) check, by position so the ID format does not matter
plink --memory 8000 --bfile SHARE_${CG}_R2 --chr 2 --from-bp 54659707 --to-bp 54659707 --freq --out sptbn1_check
cat sptbn1_check.frq
# after checking the counts, the per-chromosome sets can be removed:  rm chr{1..22}.{bed,bim,fam,log,nosex}
