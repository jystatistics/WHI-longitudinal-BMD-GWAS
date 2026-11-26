#!/bin/bash
#SBATCH --job-name=merge_whi_v2
#SBATCH --output=merge_whi_v2_%j.out
#SBATCH --error=merge_whi_v2_%j.err
#SBATCH --time=48:00:00
#SBATCH --mem=200G
#SBATCH --cpus-per-task=6
#SBATCH --account=YOUR_ACCOUNT
#SBATCH --mail-type=END,FAIL
# 2026-11-24: WHISP handled separately (missing for WHISP women where it fails). 2026-11-23 (draft 2026-11-21). GARNET c1 confirmed: built from 02_after_imputation/01_no_GRS (chr*.c1.log). Rebuilds the 5-cohort WHI genotype set (replaces WHI_merged_all).
#  - SHARE and PAGE from the rebuilt, R2-filtered files; PAGE c0 (NRUP, no consent) excluded
#  - rule: keep variants with imputation R2 >= 0.3 in every non-WHISP batch imputed for that chromosome; WHISP genotypes kept only where R2 >= 0.3 in both WHISP batches
#  - same participants as before (matched on IID to the old filtered .fam files)
#  - PCA on LD-pruned common SNPs (the old PCA used all 39M variants, unpruned)
# Writes to ${WHI_DIR}/10_Rebuild_v2/ ; old files untouched.
set -euo pipefail
source "${SLURM_SUBMIT_DIR:-$PWD}/config.sh"   # site paths (see config.sh)
W=${WHI_DBGAP_DIR}
OUT=${WHI_DIR}/10_Rebuild_v2
mkdir -p ${OUT}/r2; cd ${OUT}

# ---------- 1. R2>=0.3 lists for each imputation batch ----------
r2list () {  # $1 = folder with chr*.info.gz, $2 = output name
  # IDs are rebuilt as chr:pos:REF:ALT from the columns (not the ID column), so every batch uses the same format
  if [ -f r2/$2.done ]; then echo "$2: already built, $(wc -l < r2/$2.txt) variants"; return; fi
  : > r2/$2.txt
  for i in $(seq 1 22); do
    if [ ! -s $1/chr$i.info.gz ]; then echo "WARNING: $2 has no chr$i.info.gz; chr$i will be missing for this batch"; continue; fi
    zcat $1/chr$i.info.gz | awk -v OFS=: '/^#/{next} $1=="SNP"{next}
      $8 ~ /R2=/ { c=$1; sub(/^chr/,"",c); if (match($8,/R2=[0-9.eE-]+/)) { r=substr($8,RSTART+3,RLENGTH-3); if (r+0>=0.3) print c,$2,$4,$5 }; next }
      { split($1,s,":"); c=s[1]; sub(/^chr/,"",c); if ($7+0>=0.3) print c,s[2],$2,$3 }' >> r2/$2.txt
  done
  echo "$2: $(wc -l < r2/$2.txt) variants R2>=0.3; first ID $(head -1 r2/$2.txt); SPTBN1 listed: $(grep -c '^2:54659707:C:G$' r2/$2.txt || true)"
  touch r2/$2.done
}
r2list ${W}/04_GARNET/phg000139.v1.GARNET_WHI.genotype-calls-matrixfmt.c1/subject_level_filtered_PLINK_sets/02_after_imputation/01_no_GRS GARNET_c1
r2list ${W}/04_GARNET/phg000139.v1.GARNET_WHI.genotype-calls-matrixfmt.c2/subject_level_filtered_PLINK_sets/02_after_imputation GARNET_c2
r2list ${W}/09_WHIMS_GWAS/phg000396.v1.NHLBI_WHI_v9_WHIMS.genotype-calls-matrixfmt.c1.HMB-IRB/02_after_imputation WHIMS_c1
r2list ${W}/09_WHIMS_GWAS/phg000396.v1.NHLBI_WHI_v9_WHIMS.genotype-calls-matrixfmt.c2.HMB-IRB-NPU/02_after_imputation WHIMS_c2
r2list ${W}/05_WHISP/phg000685.v1.ESP_WHISP_6800.genotype-calls-vcf.c1/01_after_imputation WHISP_c1
r2list ${W}/05_WHISP/phg000685.v1.ESP_WHISP_6800.genotype-calls-vcf.c2/01_after_imputation WHISP_c2
cp ${W}/02_PAGE/05_Rebuild/c1/r2pass_v2.txt r2/PAGE_c1.txt
cp ${W}/02_PAGE/05_Rebuild/c2/r2pass_v2.txt r2/PAGE_c2.txt
for c in c1 c2; do cp ${W}/03_NHLBI/phg000061.v2.p1.WHI_NHLBI.genotype-calls-matrixfmt.$c/subject_level_filtered_PLINK_sets/05_converted_R2/r2pass.txt r2/SHARE_$c.txt; done
# Core batches (all except WHISP), per chromosome: keep variants with R2>=0.3 in every core batch imputed for that chromosome
# WHISP (exome-based imputation, many Michigan-excluded chunks, no chr11/chr16) is handled separately below:
# a variant that fails in WHISP is kept for everyone else and is missing for the WHISP women.
CORE=$(ls r2/*.txt | grep -v '/WHISP_')
for f in $CORE; do cut -d: -f1 $f | sort -u | sed "s/^/$(basename $f .txt) /"; done > r2_chroms.txt
awk '{n[$2]++} END{for(c in n) print c, n[c]}' r2_chroms.txt | sort -n > r2_nbatch_per_chr.txt
echo "Core batches imputed per chromosome (chr nbatches):"; cat r2_nbatch_per_chr.txt | tr '\n' ' '; echo
cat $CORE | sort | uniq -c | awk 'NR==FNR{nb[$1]=$2; next} {split($2,a,":"); if ($1==nb[a[1]]) print $2}' r2_nbatch_per_chr.txt - > keep_variants.txt
echo "Variants with R2>=0.3 in every core batch imputed for that chromosome: $(wc -l < keep_variants.txt)"
cat r2/WHISP_c1.txt r2/WHISP_c2.txt | sort | uniq -d > whisp_pass.txt
awk 'NR==FNR{k[$1];next} ($1 in k)' whisp_pass.txt keep_variants.txt > keep_whisp.txt
echo "Of these, also R2>=0.3 in both WHISP batches (kept for WHISP women): $(wc -l < keep_whisp.txt)"
echo "Kept variants per chromosome, all / WHISP:"; for k in keep_variants.txt keep_whisp.txt; do awk -F: '{n[$1]++} END{for(c=1;c<=22;c++) printf "%d:%d ", c, n[c]; print ""}' $k; done
for f in r2/PAGE_c1.txt r2/PAGE_c2.txt r2/SHARE_c1.txt r2/SHARE_c2.txt; do echo "$f: $(wc -l < $f); first ID $(head -1 $f); SPTBN1 listed: $(grep -c '^2:54659707:C:G$' $f || true)"; done
[ $(wc -l < keep_variants.txt) -ge 1000000 ] || { echo "Too few shared variants: ID formats probably differ between batches (see first IDs above)"; exit 1; }

# ---------- 2. Same participants as the published analysis (match on IID) ----------
keepids () {  # $1 old fam, $2 new bfile, $3 name
  awk 'NR==FNR{k[$2];next} ($2 in k){print $1,$2}' $1 $2.fam > keep_$3.txt
  echo "$3: old $(wc -l < $1)  matched $(wc -l < keep_$3.txt)"
}
G=${W}/04_GARNET/01_Combined_C1_C2/GARNET_merged_filtered_rearranged
M=${W}/09_WHIMS_GWAS/01_Combined_C1_C2/WHIMS_merged_filtered_rearranged
S=${W}/05_WHISP/01_Combined_C1_C2/WHISP_merged_filtered_rearranged
P1=${W}/02_PAGE/05_Rebuild/c1/PAGE_c1_R2v2
P2=${W}/02_PAGE/05_Rebuild/c2/PAGE_c2_R2v2
H1=${W}/03_NHLBI/phg000061.v2.p1.WHI_NHLBI.genotype-calls-matrixfmt.c1/subject_level_filtered_PLINK_sets/05_converted_R2/SHARE_c1_R2
H2=${W}/03_NHLBI/phg000061.v2.p1.WHI_NHLBI.genotype-calls-matrixfmt.c2/subject_level_filtered_PLINK_sets/05_converted_R2/SHARE_c2_R2
keepids ${W}/02_PAGE/04_Merged/merged_page_cleaned.fam ${P1} PAGE_c1
keepids ${W}/02_PAGE/04_Merged/merged_page_cleaned.fam ${P2} PAGE_c2
keepids ${W}/03_NHLBI/02_C1_C2_Merged/SHARE_merged_filtered_rearranged.fam ${H1} SHARE_c1
keepids ${W}/03_NHLBI/02_C1_C2_Merged/SHARE_merged_filtered_rearranged.fam ${H2} SHARE_c2

# ---------- 3. Subset every cohort to the common variant set, then merge ----------
plink --memory 180000 --bfile ${G} --extract keep_variants.txt --make-bed --out GARNET_v2
plink --memory 180000 --bfile ${M} --extract keep_variants.txt --make-bed --out WHIMS_v2
plink --memory 180000 --bfile ${S} --extract keep_whisp.txt --make-bed --out WHISP_v2
plink --memory 180000 --bfile ${P1} --keep keep_PAGE_c1.txt --extract keep_variants.txt --make-bed --out PAGE_c1_v2
plink --memory 180000 --bfile ${P2} --keep keep_PAGE_c2.txt --extract keep_variants.txt --make-bed --out PAGE_c2_v2
plink --memory 180000 --bfile ${H1} --keep keep_SHARE_c1.txt --extract keep_variants.txt --make-bed --out SHARE_c1_v2
plink --memory 180000 --bfile ${H2} --keep keep_SHARE_c2.txt --extract keep_variants.txt --make-bed --out SHARE_c2_v2
for f in GARNET_v2 WHIMS_v2 WHISP_v2 PAGE_c1_v2 PAGE_c2_v2 SHARE_c1_v2 SHARE_c2_v2; do
  ref=keep_variants.txt; [ $f = WHISP_v2 ] && ref=keep_whisp.txt
  echo "$f: $(wc -l < $f.bim) of $(wc -l < $ref) variants in $ref, $(wc -l < $f.fam) people; first bim ID $(head -1 $f.bim | cut -f2)"
  [ $(wc -l < $f.bim) -ge $(( $(wc -l < $ref) * 9 / 10 )) ] || { echo "$f is missing >10% of shared variants: its .bim IDs are probably not chr:pos:REF:ALT"; exit 1; }
done
printf "%s\n" WHIMS_v2 WHISP_v2 PAGE_c1_v2 PAGE_c2_v2 SHARE_c1_v2 SHARE_c2_v2 > mergelist.txt
plink --memory 180000 --bfile GARNET_v2 --merge-list mergelist.txt --make-bed --out WHI_merged_v2
if [ -f WHI_merged_v2-merge.missnp ]; then echo "Merge conflicts: $(wc -l < WHI_merged_v2-merge.missnp) variants"; exit 1; fi
echo "Merged: $(wc -l < WHI_merged_v2.bim) variants, $(wc -l < WHI_merged_v2.fam) people"

# ---------- 4. Checks ----------
plink --memory 180000 --bfile WHI_merged_v2 --missing --out WHI_merged_v2_miss
plink --memory 180000 --bfile WHI_merged_v2 --snp 2:54659707:C:G --freq --out sptbn1_v2; cat sptbn1_v2.frq

# ---------- 5. PCA on LD-pruned common autosomal SNPs ----------
plink --memory 180000 --bfile WHI_merged_v2 --maf 0.05 --geno 0.02 --hwe 1e-10 --indep-pairwise 200 50 0.2 --out prune
plink --memory 180000 --bfile WHI_merged_v2 --extract prune.prune.in --pca 20 --out WHI_PCA_v2
echo "Done."
