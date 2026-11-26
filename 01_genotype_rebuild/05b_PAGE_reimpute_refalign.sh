#!/bin/bash
#SBATCH --job-name=page_refalign
#SBATCH --output=page_refalign_%j.out
#SBATCH --time=03:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --account=YOUR_ACCOUNT
#SBATCH --mail-type=END,FAIL
# 2026-11-24. Michigan 2.0.11 rejected the PAGE re-imputation uploads: "More than 100 allele switches"
# (c1: 6,390 SNPs had REF/ALT the other way round from HRC). The older server accepted this, the new one does not.
# Fix (same as SHARE, 01b): force the HRC reference allele to be the VCF REF, then sort and index.
# Output: 02_PAGE/05_Rebuild/reimpute/cX/vcf_refaligned/  -> upload these files instead.
set -euo pipefail
source "${SLURM_SUBMIT_DIR:-$PWD}/config.sh"   # site paths (see config.sh)
module load bcftools/1.17
HRC=${TOOLS_DIR}/HRC.r1-1.GRCh37.wgs.mac5.sites.tab
for CG in c1 c2; do
  [ $CG = c1 ] && CHRS="6 9 17 19" || CHRS="6 14 20 21 22"
  cd ${WHI_DBGAP_DIR}/02_PAGE/05_Rebuild/reimpute/$CG
  mkdir -p vcf_refaligned
  echo "===== PAGE $CG"
  for i in $CHRS; do
    plink --memory 28000 --vcf PAGE_${CG}-chr$i.vcf.gz --double-id --keep-allele-order --make-bed --out tmp_chr$i > /dev/null
  done
  cat $(for i in $CHRS; do echo tmp_chr$i.bim; done) > all.bim
  # SNP ID <tab> HRC REF allele, matched on chr:pos where the SNP's two alleles are HRC's REF/ALT
  awk 'NR==FNR{k=$1":"$4; id[k]=$2; a[k]=$5; b[k]=$6; next}
       FNR>1 { k=$1":"$2; if ((k in id) && !(k in done) && (($4==a[k] && $5==b[k]) || ($4==b[k] && $5==a[k]))) { print id[k]"\t"$4; done[k]=1 } }' all.bim $HRC > hrc_ref_alleles.txt
  echo "SNPs: $(wc -l < all.bim)   with HRC REF found: $(wc -l < hrc_ref_alleles.txt)"
  for i in $CHRS; do
    plink --memory 28000 --bfile tmp_chr$i --a2-allele hrc_ref_alleles.txt 2 1 --recode vcf-iid bgz --out tmp_ref_chr$i > /dev/null
    bcftools sort tmp_ref_chr$i.vcf.gz -Oz -o vcf_refaligned/PAGE_${CG}-chr$i.vcf.gz 2>/dev/null
    bcftools index -t vcf_refaligned/PAGE_${CG}-chr$i.vcf.gz
    bcftools query -f '%ID\t%REF\n' vcf_refaligned/PAGE_${CG}-chr$i.vcf.gz | \
      awk -v c=$i 'NR==FNR{r[$1]=$2;next} ($1 in r){n++; if(r[$1]==$2)m++} END{printf "chr%s: REF matches HRC for %d of %d HRC SNPs; samples: ", c, m, n}' hrc_ref_alleles.txt -
    bcftools query -l vcf_refaligned/PAGE_${CG}-chr$i.vcf.gz | wc -l
  done
done
echo "Upload reimpute/c1/vcf_refaligned/*.vcf.gz and reimpute/c2/vcf_refaligned/*.vcf.gz (not the .tbi files)."
