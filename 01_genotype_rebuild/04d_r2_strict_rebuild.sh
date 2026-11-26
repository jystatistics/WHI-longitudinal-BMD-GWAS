#!/bin/bash
#SBATCH --job-name=r2strict
#SBATCH --output=r2strict_%j.out
#SBATCH --time=24:00:00
#SBATCH --mem=160G
#SBATCH --cpus-per-task=8
#SBATCH --account=YOUR_ACCOUNT
#SBATCH --mail-type=END,FAIL
# 2026-11-25. Stricter imputation-quality rule for hard-called genotypes (duplicate PAGE/SHARE check showed
# ~21% hard-call discordance for common SNPs at R2>=0.3). Builds WHI_v3_r08 from the per-cohort *_v2 files:
#  - R2 >= RT (default 0.8) in every non-WHISP batch imputed for that chromosome
#  - WHISP women kept where R2 >= RT in both WHISP batches, otherwise missing
#  - 2,649 unique women (PAGE copies of the 912 PAGE/SHARE duplicates removed)
#  (2026-11-25 fix: the awk threshold variable was named RT, which is gawk's built-in record terminator,
#   so every variant passed; renamed THR.)
#  - PCA on well-imputed, LD-pruned common SNPs outside long-range LD regions; reports PC outliers and relatedness
# Submit: sbatch 04d_r2_strict_rebuild.sh        (or  sbatch --export=ALL,RT=0.9 04d_r2_strict_rebuild.sh)
set -euo pipefail
source "${SLURM_SUBMIT_DIR:-$PWD}/config.sh"   # site paths (see config.sh)
RT=${RT:-0.8}; TAG=r$(echo $RT | tr -d '.')        # 0.8 -> r08
W=${WHI_DBGAP_DIR}
cd ${WHI_DIR}/10_Rebuild_v2
mkdir -p r2_$TAG
r2list () {  # $1 = folder with chr*.info.gz, $2 = output name, $3 = optional list of chromosomes
  local out=r2_$TAG/$2.txt; : > $out
  for i in ${3:-$(seq 1 22)}; do
    if [ ! -s $1/chr$i.info.gz ]; then echo "  $2: no chr$i.info.gz"; continue; fi
    zcat $1/chr$i.info.gz | awk -v OFS=: -v THR=$RT '/^#/{next} $1=="SNP"{next}
      $8 ~ /R2=/ { c=$1; sub(/^chr/,"",c); if (match($8,/R2=[0-9.eE-]+/)) { r=substr($8,RSTART+3,RLENGTH-3); if (r+0>=THR) print c,$2,$4,$5 }; next }
      { split($1,s,":"); c=s[1]; sub(/^chr/,"",c); if ($7+0>=THR) print c,s[2],$2,$3 }' >> $out
  done
  echo "$2: $(wc -l < $out) variants R2>=$RT; SPTBN1 listed: $(grep -c '^2:54659707:C:G$' $out || true)"
}
r2list ${W}/04_GARNET/phg000139.v1.GARNET_WHI.genotype-calls-matrixfmt.c1/subject_level_filtered_PLINK_sets/02_after_imputation/01_no_GRS GARNET_c1
r2list ${W}/04_GARNET/phg000139.v1.GARNET_WHI.genotype-calls-matrixfmt.c2/subject_level_filtered_PLINK_sets/02_after_imputation GARNET_c2
r2list ${W}/09_WHIMS_GWAS/phg000396.v1.NHLBI_WHI_v9_WHIMS.genotype-calls-matrixfmt.c1.HMB-IRB/02_after_imputation WHIMS_c1
r2list ${W}/09_WHIMS_GWAS/phg000396.v1.NHLBI_WHI_v9_WHIMS.genotype-calls-matrixfmt.c2.HMB-IRB-NPU/02_after_imputation WHIMS_c2
r2list ${W}/05_WHISP/phg000685.v1.ESP_WHISP_6800.genotype-calls-vcf.c1/01_after_imputation WHISP_c1
r2list ${W}/05_WHISP/phg000685.v1.ESP_WHISP_6800.genotype-calls-vcf.c2/01_after_imputation WHISP_c2
for c in c1 c2; do
  r2list ${W}/03_NHLBI/phg000061.v2.p1.WHI_NHLBI.genotype-calls-matrixfmt.$c/subject_level_filtered_PLINK_sets/04_after_imputation_hg19 SHARE_$c
done
# PAGE: original imputation, except the re-imputed chromosomes
r2list ${W}/02_PAGE/phg000910.v1.PAGE_WHI.genotype-calls-matrixfmt.c1/filtered_subject_level_set/02_after_imputation PAGE_c1 "1 2 3 4 5 7 8 10 11 12 13 14 15 16 18 20 21 22"
mv r2_$TAG/PAGE_c1.txt r2_$TAG/PAGE_c1_orig.txt; r2list ${W}/02_PAGE/05_Rebuild/reimpute/c1/results PAGE_c1 "6 9 17 19"; cat r2_$TAG/PAGE_c1_orig.txt >> r2_$TAG/PAGE_c1.txt; rm r2_$TAG/PAGE_c1_orig.txt
r2list ${W}/02_PAGE/phg000910.v1.PAGE_WHI.genotype-calls-matrixfmt.c2/filtered_subject_level_set/02_after_imputation PAGE_c2 "1 2 3 4 5 7 8 9 10 11 12 13 15 16 17 18 19"
mv r2_$TAG/PAGE_c2.txt r2_$TAG/PAGE_c2_orig.txt; r2list ${W}/02_PAGE/05_Rebuild/reimpute/c2/results PAGE_c2 "6 14 20 21 22"; cat r2_$TAG/PAGE_c2_orig.txt >> r2_$TAG/PAGE_c2.txt; rm r2_$TAG/PAGE_c2_orig.txt
for b in PAGE_c1 PAGE_c2; do echo "$b total: $(wc -l < r2_$TAG/$b.txt)"; done

# same rule as 04_merge_all_cohorts.sh
CORE=$(ls r2_$TAG/*.txt | grep -v '/WHISP_')
for f in $CORE; do cut -d: -f1 $f | sort -u | sed "s/^/$(basename $f .txt) /"; done > r2_${TAG}_chroms.txt
awk '{n[$2]++} END{for(c in n) print c, n[c]}' r2_${TAG}_chroms.txt | sort -n > r2_${TAG}_nbatch.txt
cat $CORE | sort | uniq -c | awk 'NR==FNR{nb[$1]=$2; next} {split($2,a,":"); if ($1==nb[a[1]]) print $2}' r2_${TAG}_nbatch.txt - > keep_$TAG.txt
cat r2_$TAG/WHISP_c1.txt r2_$TAG/WHISP_c2.txt | sort | uniq -d | awk 'NR==FNR{k[$1];next} ($1 in k)' keep_$TAG.txt - > keep_${TAG}_whisp.txt
echo "R2>=$RT in every core batch: $(wc -l < keep_$TAG.txt); also in both WHISP batches: $(wc -l < keep_${TAG}_whisp.txt)"
echo "per chromosome, all / WHISP:"; for k in keep_$TAG.txt keep_${TAG}_whisp.txt; do awk -F: '{n[$1]++} END{for(c=1;c<=22;c++) printf "%d:%d ", c, n[c]; print ""}' $k; done

# subset the per-cohort files and merge (PAGE duplicates removed)
for f in GARNET_v2 WHIMS_v2 PAGE_c1_v2 PAGE_c2_v2 SHARE_c1_v2 SHARE_c2_v2; do
  plink --memory 150000 --bfile $f --extract keep_$TAG.txt --remove remove_page_dups.txt --make-bed --out ${f%_v2}_$TAG > /dev/null
done
plink --memory 150000 --bfile WHISP_v2 --extract keep_${TAG}_whisp.txt --make-bed --out WHISP_$TAG > /dev/null
for f in GARNET WHIMS WHISP PAGE_c1 PAGE_c2 SHARE_c1 SHARE_c2; do echo "${f}_$TAG: $(wc -l < ${f}_$TAG.bim) variants, $(wc -l < ${f}_$TAG.fam) women"; done
printf "%s\n" WHIMS_$TAG WHISP_$TAG PAGE_c1_$TAG PAGE_c2_$TAG SHARE_c1_$TAG SHARE_c2_$TAG > mergelist_$TAG.txt
plink --memory 150000 --bfile GARNET_$TAG --merge-list mergelist_$TAG.txt --make-bed --out WHI_v3_$TAG > /dev/null
[ ! -f WHI_v3_$TAG-merge.missnp ] || { echo "merge conflicts: $(wc -l < WHI_v3_$TAG-merge.missnp)"; exit 1; }
echo "WHI_v3_$TAG: $(wc -l < WHI_v3_$TAG.bim) variants, $(wc -l < WHI_v3_$TAG.fam) women, unique IIDs $(awk '{print $2}' WHI_v3_$TAG.fam | sort -u | wc -l)"
plink --memory 150000 --bfile WHI_v3_$TAG --snp 2:54659707:C:G --freq --out sptbn1_$TAG > /dev/null 2>&1 && cat sptbn1_$TAG.frq || echo "SPTBN1 (2:54659707:C:G) not in the R2>=$RT set"

# PCA on well-imputed common SNPs, excluding long-range LD regions (hg19)
printf "6 25000000 33500000 MHC\n8 8000000 12000000 inv8p23\n17 40900000 45000000 inv17q21\n2 129000000 137000000 LCT\n11 87500000 90500000 r11\n" > longLD_hg19.txt
plink --memory 150000 --threads 8 --bfile WHI_v3_$TAG --maf 0.05 --geno 0.02 --hwe 1e-10 --exclude range longLD_hg19.txt --indep-pairwise 200 50 0.1 --out prune_$TAG > /dev/null
echo "PCA SNPs: $(wc -l < prune_$TAG.prune.in); per chr: $(awk -F: '{n[$1]++} END{for(c=1;c<=22;c++) printf "%d:%d ", c, n[c]}' prune_$TAG.prune.in)"
plink --memory 150000 --threads 8 --bfile WHI_v3_$TAG --extract prune_$TAG.prune.in --pca 20 --out WHI_PCA_v3_$TAG > /dev/null
awk '{for(j=3;j<=12;j++){s[j]+=$j; q[j]+=$j*$j}; n++; v[n]=$0} END{for(j=3;j<=12;j++){m[j]=s[j]/n; sd[j]=sqrt(q[j]/n-m[j]^2)}
  for(i=1;i<=n;i++){split(v[i],x," "); for(j=3;j<=12;j++) if((x[j]-m[j])/sd[j]>6||(x[j]-m[j])/sd[j]<-6) out[j-2]++}
  printf "PC outliers (|z|>6) by PC1-10:"; for(p=1;p<=10;p++) printf " %d", out[p]; print ""}' WHI_PCA_v3_$TAG.eigenvec
plink --memory 150000 --threads 8 --bfile WHI_v3_$TAG --extract prune_$TAG.prune.in --genome --min 0.2 --out rel_$TAG > /dev/null
awk 'NR>1 {if($10>=0.9) d++; else if($10>=0.354) f++; else if($10>=0.177) s++} END{print "pairs (PLINK PI_HAT): >=0.9 " d+0 ", 0.354-0.9 " f+0 ", 0.2-0.354 " s+0}' rel_$TAG.genome
command -v plink2 >/dev/null && { plink2 --bfile WHI_v3_$TAG --extract prune_$TAG.prune.in --make-king-table --king-table-filter 0.0884 --out king_$TAG > /dev/null;
  awk 'NR>1 {k=$NF; if(k>=0.354) d++; else if(k>=0.177) f++; else s++} END{print "KING kinship: >=0.354 (dup/MZ) " d+0 ", 1st degree " f+0 ", 2nd degree " s+0}' king_$TAG.kin0; } || echo "plink2 not on PATH: KING skipped"
echo "Done."
