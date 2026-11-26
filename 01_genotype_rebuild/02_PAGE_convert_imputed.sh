#!/bin/bash
#SBATCH --job-name=page_convert
#SBATCH --output=page_convert_%A_%a.out
#SBATCH --error=page_convert_%A_%a.err
#SBATCH --time=24:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --account=YOUR_ACCOUNT
#SBATCH --mail-type=END,FAIL
#SBATCH --array=0-2
# DRAFT 2026-11-21. Converts PAGE c0/c1/c2 Michigan (HRC r1.1) dose VCFs to PLINK, keeping only
# variants with imputation R2 >= 0.3 in that batch. Replaces merged_c1/merged_c2, which were built
# from ~680k array SNPs instead of the imputed data. Writes only to 02_PAGE/05_Rebuild/.
# Submit:  sbatch 02_PAGE_convert_imputed.sh      (runs c0, c1, c2 as array tasks 0,1,2)
# NOTE: drop c0 (use --array=1-2) if c0 turns out to be a no-consent group.
set -euo pipefail
source "${SLURM_SUBMIT_DIR:-$PWD}/config.sh"   # site paths (see config.sh)
GROUPS_=(c0 c1 c2)
CG=${GROUPS_[$SLURM_ARRAY_TASK_ID]}
P=${WHI_DBGAP_DIR}/02_PAGE
D=${P}/phg000910.v1.PAGE_WHI.genotype-calls-matrixfmt.${CG}/filtered_subject_level_set/02_after_imputation
OUT=${P}/05_Rebuild/${CG}
mkdir -p ${OUT}; cd ${OUT}

# 1. Variants with R2 >= 0.3 (Minimac info: col1 SNP id, col7 Rsq)
: > r2pass.txt
for i in $(seq 1 22); do
  zcat ${D}/chr$i.info.gz | awk 'NR>1 && $7>=0.3 {print $1}' >> r2pass.txt
done
echo "${CG}: variants with R2>=0.3: $(wc -l < r2pass.txt)"

# 2. Per-chromosome conversion (hard calls from GT, same as the other cohorts)
: > mergelist.txt
for i in $(seq 1 22); do
  plink --memory 56000 --vcf ${D}/chr$i.dose.vcf.gz --extract r2pass.txt --make-bed --out chr$i
  [ $i -gt 1 ] && echo "chr$i" >> mergelist.txt
done
plink --memory 56000 --bfile chr1 --merge-list mergelist.txt --make-bed --out PAGE_${CG}_R2
rm -f chr[0-9]*.bed chr[0-9]*.bim chr[0-9]*.fam chr[0-9]*.log chr[0-9]*.nosex

# 3. Report
echo "${CG}: final variants $(wc -l < PAGE_${CG}_R2.bim), people $(wc -l < PAGE_${CG}_R2.fam), chr21 $(awk '$1==21' PAGE_${CG}_R2.bim | wc -l)"
plink --memory 56000 --bfile PAGE_${CG}_R2 --snp 2:54659707:C:G --freq --out sptbn1_check || echo "SPTBN1 not present (R2<0.3?)"
cat sptbn1_check.frq 2>/dev/null || true
zcat ${D}/chr2.info.gz | grep -m1 '^2:54659707' || true
