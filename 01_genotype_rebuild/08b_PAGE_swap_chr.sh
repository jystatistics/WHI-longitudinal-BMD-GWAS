#!/bin/bash
#SBATCH --job-name=page_swap
#SBATCH --output=page_swap_%j.out
#SBATCH --time=06:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --account=YOUR_ACCOUNT
#SBATCH --mail-type=END,FAIL
# 2026-11-24. Replace the re-imputed chromosomes in PAGE_cX_R2 -> PAGE_cX_R2v2 (old files untouched),
# and write r2pass_v2.txt. Run after every task of 08_PAGE_reimputed_download.sh has finished.
# 2026-11-25: the re-imputed .fam files hold the joined Michigan sample name (FID_IID, read with --double-id),
# while PAGE_cX_R2 split it into FID and IID. Map them back to the PAGE_cX_R2 IDs before merging.
set -euo pipefail
source "${SLURM_SUBMIT_DIR:-$PWD}/config.sh"   # site paths (see config.sh)
P=${WHI_DBGAP_DIR}/02_PAGE/05_Rebuild
for CG in c1 c2; do
  [ $CG = c1 ] && CHRS="6 9 17 19" || CHRS="6 14 20 21 22"
  R=$P/reimpute/$CG/results; cd $P/$CG
  for i in $CHRS; do [ -s $R/chr$i.bim ] || { echo "missing $R/chr$i.bim"; exit 1; }; done
  for i in $CHRS; do   # FID_IID -> FID IID as in PAGE_cX_R2.fam (safe to run twice)
    awk 'NR==FNR{j[$1"_"$2]=$1" "$2; ok[$1" "$2]; next}
         { k=$1" "$2; if (!(k in ok)) { if ($2 in j) { split(j[$2],a," "); $1=a[1]; $2=a[2] } else bad++ } print }
         END{ if (bad) { print bad" samples not found in PAGE_cX_R2.fam" > "/dev/stderr"; exit 1 } }' \
      PAGE_${CG}_R2.fam $R/chr$i.fam > $R/chr$i.fam.tmp
    mv $R/chr$i.fam.tmp $R/chr$i.fam
    echo "$CG chr$i: $(awk 'NR==FNR{k[$1" "$2];next} ($1" "$2 in k)' PAGE_${CG}_R2.fam $R/chr$i.fam | wc -l) of $(wc -l < $R/chr$i.fam) re-imputed samples match PAGE_${CG}_R2.fam"
  done
  plink --memory 56000 --bfile PAGE_${CG}_R2 --not-chr $(echo $CHRS | tr ' ' ',') --make-bed --out tmp_rest > /dev/null
  for i in $CHRS; do echo $R/chr$i; done > swap_list.txt
  plink --memory 56000 --bfile tmp_rest --merge-list swap_list.txt --make-bed --out PAGE_${CG}_R2v2 > /dev/null
  awk -F: -v chrs=" $CHRS " 'index(chrs, " "$1" ")==0' r2pass.txt > r2pass_v2.txt
  for i in $CHRS; do cat $R/r2pass_chr$i.txt >> r2pass_v2.txt; done
  echo "== PAGE $CG: $(wc -l < PAGE_${CG}_R2v2.bim) variants, $(wc -l < PAGE_${CG}_R2v2.fam) people; r2pass_v2 $(wc -l < r2pass_v2.txt)"
  [ $(wc -l < PAGE_${CG}_R2v2.fam) -eq $(wc -l < PAGE_${CG}_R2.fam) ] || { echo "People count changed ($(wc -l < PAGE_${CG}_R2.fam) -> $(wc -l < PAGE_${CG}_R2v2.fam)): sample IDs do not match"; exit 1; }
  awk -F: '{n[$1]++} END{for(c=1;c<=22;c++) printf "%d:%dk ", c, n[c]/1000; print ""}' r2pass_v2.txt
  plink --memory 56000 --bfile PAGE_${CG}_R2v2 --snp 2:54659707:C:G --freq --out sptbn1_v2 > /dev/null; cat sptbn1_v2.frq
done
