#!/bin/bash
#SBATCH --job-name=share_prep
#SBATCH --output=share_prep_%j.out
#SBATCH --error=share_prep_%j.err
#SBATCH --time=12:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --account=YOUR_ACCOUNT
#SBATCH --mail-type=END,FAIL
# DRAFT v2 2026-11-21 (liftOver version; no Oxford strand file needed).
# SHARE (Affy 6.0) .bim has build-36 (hg18) positions; the earlier HRC (GRCh37) check matched only
# 13,341 SNPs, so only ~1,950 SNPs were imputed. This job lifts positions hg18 -> hg19, drops A/T and
# C/G SNPs, aligns strand/alleles to HRC GRCh37 (no frequency filter), and writes per-chr VCFs
# for the TOPMed Imputation Server.
# Before running, on the LOGIN node:
#   cd ${TOOLS_DIR}
#   wget https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/liftOver && chmod +x liftOver
#   wget https://hgdownload.soe.ucsc.edu/goldenPath/hg18/liftOver/hg18ToHg19.over.chain.gz
# Submit:  sbatch --export=CG=c1 01_SHARE_prep_for_TOPMed.sh   (then CG=c2)
set -euo pipefail
source "${SLURM_SUBMIT_DIR:-$PWD}/config.sh"   # site paths (see config.sh)
module load bcftools/1.17

CG=${CG:-c1}
BASE=${WHI_DBGAP_DIR}/03_NHLBI/phg000061.v2.p1.WHI_NHLBI.genotype-calls-matrixfmt.${CG}/subject_level_filtered_PLINK_sets
IN=${BASE}/NHLBI_SHARE_WHI_subject_level_${CG}
OUT=${BASE}/03_rebuild_TOPMed
REF=${TOOLS_DIR}
for f in ${REF}/liftOver ${REF}/hg18ToHg19.over.chain.gz ${REF}/HRC-1000G-check-bim.pl ${REF}/HRC.r1-1.GRCh37.wgs.mac5.sites.tab; do
  [ -s "$f" ] || { echo "Missing $f"; exit 1; }
done
mkdir -p ${OUT}/vcf; cd ${OUT}

# 1. liftOver hg18 -> hg19 (autosomes + X)
awk '$1>=1 && $1<=23 && $4>0 {c=($1==23)?"X":$1; print "chr"c"\t"$4-1"\t"$4"\t"$2}' ${IN}.bim > hg18.bed
${REF}/liftOver hg18.bed ${REF}/hg18ToHg19.over.chain.gz hg19.bed unmapped.bed
awk '$1 ~ /^chr([0-9]+|X)$/ {c=$1; sub("chr","",c); if(c=="X")c=23; print $4, c, $3}' hg19.bed > lifted.txt
echo "hg18 SNPs: $(wc -l < hg18.bed)   lifted: $(wc -l < lifted.txt)"
awk '{print $1}' lifted.txt > keep.snps
awk '{print $1, $2}' lifted.txt > chr.txt
awk '{print $1, $3}' lifted.txt > pos.txt
plink --bfile ${IN} --extract keep.snps --make-bed --out tmp1
plink --bfile tmp1 --update-chr chr.txt --make-bed --out tmp2
plink --bfile tmp2 --update-map pos.txt --make-bed --out tmp3

# 2. Drop strand-ambiguous A/T, C/G SNPs; basic QC
awk '($5=="A"&&$6=="T")||($5=="T"&&$6=="A")||($5=="C"&&$6=="G")||($5=="G"&&$6=="C"){print $2}' tmp3.bim > palin.snps
plink --bfile tmp3 --exclude palin.snps --chr 1-22 \
      --geno 0.05 --mind 0.05 --maf 0.01 --hwe 1e-20 \
      --make-bed --out SHARE_${CG}_b37_qc
plink --bfile SHARE_${CG}_b37_qc --freq --out SHARE_${CG}_b37_qc
rm -f tmp1.* tmp2.* tmp3.*

# 3. Align strand / ref-alt to HRC GRCh37. -t 1 = no allele-frequency exclusion
#    (HRC is mostly European; frequency filtering removes good SNPs in Black/Hispanic women)
perl ${REF}/HRC-1000G-check-bim.pl -b SHARE_${CG}_b37_qc.bim -f SHARE_${CG}_b37_qc.frq \
     -r ${REF}/HRC.r1-1.GRCh37.wgs.mac5.sites.tab -h -t 1
grep -E 'Position Matches|No Match|Total in bim|Strand to change|Non Matching' LOG-SHARE_${CG}_b37_qc-HRC.txt || true
sh Run-plink.sh

# 4. Report: expect ~600k+ SNPs total; chr21 ~8-10k
cat SHARE_${CG}_b37_qc-updated-chr*.bim | wc -l
wc -l < SHARE_${CG}_b37_qc-updated-chr21.bim

# 5. Per-chromosome sorted, bgzipped VCFs for upload
for i in $(seq 1 22); do
  plink --bfile SHARE_${CG}_b37_qc-updated-chr$i --recode vcf-iid bgz --keep-allele-order --out tmp_chr$i
  bcftools sort tmp_chr$i.vcf.gz -Oz -o vcf/SHARE_${CG}_chr$i.vcf.gz
  bcftools index -t vcf/SHARE_${CG}_chr$i.vcf.gz
  rm tmp_chr$i.vcf.gz
done
echo "Done. Upload vcf/*.vcf.gz to the Michigan server with the same panel used for the other cohorts."
