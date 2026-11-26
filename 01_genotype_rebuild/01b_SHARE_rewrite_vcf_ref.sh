#!/bin/bash
#SBATCH --job-name=share_vcf_ref
#SBATCH --output=share_vcf_ref_%A_%a.out
#SBATCH --error=share_vcf_ref_%A_%a.err
#SBATCH --time=4:00:00
#SBATCH --mem=32G
#SBATCH --account=YOUR_ACCOUNT
#SBATCH --mail-type=END,FAIL
#SBATCH --array=1-2
# 2026-11-22. Michigan QC reported "Allele switch: 654,227 / Match: 0": every SNP had REF/ALT swapped.
# Cause: the HRC check script sets the HRC reference allele as PLINK A1, but PLINK 1.9 writes A2 as
# the VCF REF. Fix: look up the HRC reference allele by chr:pos and force it to be A2 (= VCF REF).
set -euo pipefail
source "${SLURM_SUBMIT_DIR:-$PWD}/config.sh"   # site paths (see config.sh)
module load bcftools/1.17
CG=c${SLURM_ARRAY_TASK_ID}
OUT=${WHI_DBGAP_DIR}/03_NHLBI/phg000061.v2.p1.WHI_NHLBI.genotype-calls-matrixfmt.${CG}/subject_level_filtered_PLINK_sets/03_rebuild_TOPMed
cd ${OUT}
HRC=${TOOLS_DIR}/HRC.r1-1.GRCh37.wgs.mac5.sites.tab
cat SHARE_${CG}_b37_qc-updated-chr*.bim > all_updated.bim
# ID <tab> HRC REF allele, matched on chr:pos (HRC tab: CHROM POS ID REF ALT ...)
# keep the HRC row whose REF/ALT are the SNP's two alleles; one line per SNP
awk 'NR==FNR{k=$1":"$4; id[k]=$2; a[k]=$5; b[k]=$6; next}
     FNR>1 { k=$1":"$2; if ((k in id) && !(k in done) && (($4==a[k] && $5==b[k]) || ($4==b[k] && $5==a[k]))) { print id[k]"\t"$4; done[k]=1 } }' all_updated.bim ${HRC} > hrc_ref_alleles.txt
REFA=hrc_ref_alleles.txt
echo "SNPs: $(wc -l < all_updated.bim)   with HRC REF found: $(wc -l < ${REFA})"
mkdir -p vcf_refaligned
for i in $(seq 1 22); do
  plink --memory 28000 --bfile SHARE_${CG}_b37_qc-updated-chr$i --a2-allele ${REFA} 2 1 \
        --recode vcf-iid bgz --out tmp_chr$i
  bcftools sort tmp_chr$i.vcf.gz -Oz -o vcf_refaligned/SHARE_${CG}_chr$i.vcf.gz
  bcftools index -t vcf_refaligned/SHARE_${CG}_chr$i.vcf.gz
  rm -f tmp_chr$i.vcf.gz tmp_chr$i.log tmp_chr$i.nosex
done
# Check: VCF REF should equal the HRC reference allele for (almost) every SNP
bcftools query -f '%ID\t%REF\n' vcf_refaligned/SHARE_${CG}_chr22.vcf.gz | \
  awk 'NR==FNR{r[$1]=$2;next} ($1 in r){n++; if(r[$1]==$2)m++} END{print "chr22 REF matches HRC:",m,"of",n}' ${REFA} -
echo "Upload vcf_refaligned/*.vcf.gz (not vcf/) to Michigan."
