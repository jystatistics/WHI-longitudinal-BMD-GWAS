#!/bin/bash
#SBATCH --job-name=share_fixmerge
#SBATCH --output=share_fixmerge_%x_%j.out
#SBATCH --error=share_fixmerge_%x_%j.err
#SBATCH --time=08:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --account=YOUR_ACCOUNT
#SBATCH --mail-type=END,FAIL
# 2026-11-23. Replaces 03b. Michigan IDs are chr:pos only, so multi-allelic sites share an ID and the merge
# fails ("3+ alleles"). Rename each variant to chr:pos:REF:ALT (same format as PAGE/WHIMS/WHISP/GARNET)
# using the info file, keep only records with their own R2 >= 0.3, then merge chr1-22.
set -euo pipefail
source "${SLURM_SUBMIT_DIR:-$PWD}/config.sh"   # site paths (see config.sh)
CG=${CG:?set CG=c1 or c2}
B=${WHI_DBGAP_DIR}/03_NHLBI/phg000061.v2.p1.WHI_NHLBI.genotype-calls-matrixfmt.${CG}/subject_level_filtered_PLINK_sets
D=${B}/04_after_imputation_hg19; OUT=${B}/05_converted_R2
cd ${OUT}
rm -f SHARE_${CG}_R2-merge.missnp SHARE_${CG}_R2-merge.fam
: > mergelist.txt
for i in $(seq 1 22); do
  [ -s chr${i}.bed ] && [ -s ${D}/chr${i}.info.gz ] || { echo "missing chr${i}.bed or chr${i}.info.gz"; exit 1; }
  [ -s chr${i}.bim.orig ] || cp chr${i}.bim chr${i}.bim.orig
  zcat ${D}/chr${i}.info.gz | awk -v c=${i} -v OFS='\t' '
    /^#/ {next}  $1=="SNP" {next}
    { if ($8 ~ /R2=/) { pos=$2; ref=$4; alt=$5; r=-1; if (match($8,/R2=[0-9.eE-]+/)) r=substr($8,RSTART+3,RLENGTH-3)+0 }
      else { split($1,s,":"); pos=s[2]; ref=$2; alt=$3; r=$7+0 }
      k=(ref<alt) ? pos":"ref":"alt : pos":"alt":"ref
      print k, c":"pos":"ref":"alt, r }' > idmap_chr${i}.txt
  : > r2pass_ids_chr${i}.txt
  awk -v OFS='\t' -v rf=r2pass_ids_chr${i}.txt 'NR==FNR { id[$1]=$2; r[$1]=$3; split($1,p,":"); m[p[1]":"p[2]]++; m[p[1]":"p[3]]++; mk[p[1]":"p[2]]=$1; mk[p[1]":"p[3]]=$1; next }
    { k=($5<$6) ? $4":"$5":"$6 : $4":"$6":"$5
      if ($5=="0" && m[$4":"$6]==1) k=mk[$4":"$6]
      if ((k in id) && !(k in seen)) { $2=id[k]; seen[k]=1; if (r[k]>=0.3) print $2 > rf }
      else { $2=$1":"$4":"$6":"$5":nomatch"FNR }
      print }' idmap_chr${i}.txt chr${i}.bim.orig > chr${i}.bim
  plink --memory 56000 --bfile chr${i} --extract r2pass_ids_chr${i}.txt --make-bed --out chr${i}_id > /dev/null
  echo "chr${i}: $(wc -l < chr${i}.bim.orig) in, $(grep -c nomatch chr${i}.bim || true) unmatched, $(wc -l < chr${i}_id.bim) kept (own R2>=0.3)"
  rm -f idmap_chr${i}.txt
  [ ${i} -gt 1 ] && echo "chr${i}_id" >> mergelist.txt
done
cat r2pass_ids_chr{1..22}.txt > r2pass.txt
plink --memory 56000 --bfile chr1_id --merge-list mergelist.txt --make-bed --out SHARE_${CG}_R2
echo "SHARE_${CG}_R2: $(wc -l < SHARE_${CG}_R2.bim) variants, $(wc -l < SHARE_${CG}_R2.fam) people"
plink --memory 8000 --bfile SHARE_${CG}_R2 --chr 2 --from-bp 54659707 --to-bp 54659707 --freq --out sptbn1_check
cat sptbn1_check.frq
