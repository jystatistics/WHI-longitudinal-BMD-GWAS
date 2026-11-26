#!/bin/bash
#SBATCH --job-name=share_dl
#SBATCH --output=share_dl_%x_%A_%a.out
#SBATCH --error=share_dl_%x_%A_%a.err
#SBATCH --time=12:00:00
#SBATCH --mem=48G
#SBATCH --cpus-per-task=4
#SBATCH --account=YOUR_ACCOUNT
#SBATCH --mail-type=FAIL
#SBATCH --array=1-22%2
# 2026-11-22. Per chromosome: download Michigan zip -> extract dose + info only -> delete zip ->
# keep the analysis participants and R2>=0.3 variants -> PLINK -> delete dose VCF.
# Never holds more than ~2 chromosomes' zips on disk (project quota is 500 GB).
# Before submitting (once per consent group, see instructions):
#   ~/mis_dl_${CG}.txt  : download commands from Michigan (curl -sL <link> > file), all jobs for that group
#   ~/.mis_pw_${CG}     : the result password(s), one per line (chmod 600)
# Submit:  sbatch --export=ALL,CG=c1 03_SHARE_download_convert.sh     (then CG=c2)
set -euo pipefail
source "${SLURM_SUBMIT_DIR:-$PWD}/config.sh"   # site paths (see config.sh)
CG=${CG:?set CG=c1 or c2}
i=${SLURM_ARRAY_TASK_ID}
B=${WHI_DBGAP_DIR}/03_NHLBI/phg000061.v2.p1.WHI_NHLBI.genotype-calls-matrixfmt.${CG}/subject_level_filtered_PLINK_sets
D=${B}/04_after_imputation_hg19; OUT=${B}/05_converted_R2
mkdir -p ${D} ${OUT}
# participants used in the published analysis (IID), FID=IID after --double-id
[ -s ${OUT}/keep_ids.txt ] || awk '{print $2, $2}' ${WHI_DBGAP_DIR}/03_NHLBI/02_C1_C2_Merged/SHARE_merged_filtered_rearranged.fam > ${OUT}/keep_ids.txt

cd ${D}
# 1. download this chromosome's zip (the curl/wget line in the Michigan list that names chr_${i}.zip;
#    the list also has "echo Downloading file chr_N.zip" lines, which are skipped)
echo "Running: $(grep -E "chr_${i}\.zip" ~/mis_dl_${CG}.txt | grep -vE "^[[:space:]]*echo" | head -1 | sed -E 's#(https?://[^/]+/)[^ ]*#\1...#g')"
CMD=$(grep -E "chr_${i}\.zip" ~/mis_dl_${CG}.txt | grep -vE "^[[:space:]]*echo" | head -1)
[ -n "$CMD" ] || { echo "No download line for chr_${i}.zip in ~/mis_dl_${CG}.txt"; exit 1; }
eval "$CMD"
[ -s chr_${i}.zip ] || { echo "Download failed for chr ${i}"; exit 1; }

# 2. extract only dose + info (try each password; chr10 of c1 has its own)
ok=0
while IFS= read -r PW; do
  if 7z x -y -p"$PW" chr_${i}.zip "chr${i}.dose.vcf.gz" "chr${i}.info.gz" >/dev/null 2>&1; then ok=1; break; fi
done < ~/.mis_pw_${CG}
[ $ok -eq 1 ] || { echo "Could not unzip chr ${i} with any password"; exit 1; }
rm -f chr_${i}.zip

# 3. R2>=0.3 list (MIS2 info is VCF-style with R2= in INFO; older format has Rsq in col 7)
zcat chr${i}.info.gz | awk '/^#/{next} $1=="SNP"{next}
  $8 ~ /R2=/ { if (match($8,/R2=[0-9.eE-]+/)) { r=substr($8,RSTART+3,RLENGTH-3); if (r+0>=0.3) print $3 }; next }
  $7+0>=0.3 {print $1}' > ${OUT}/r2pass_chr${i}.txt

# 4. convert, keeping analysis participants and R2>=0.3 variants
plink --memory 40000 --vcf chr${i}.dose.vcf.gz --double-id --keep ${OUT}/keep_ids.txt \
      --extract ${OUT}/r2pass_chr${i}.txt --make-bed --out ${OUT}/chr${i}
rm -f chr${i}.dose.vcf.gz
echo "chr${i}: R2>=0.3 $(wc -l < ${OUT}/r2pass_chr${i}.txt), kept $(wc -l < ${OUT}/chr${i}.bim) variants, $(wc -l < ${OUT}/chr${i}.fam) people"
