# Site-specific paths. Edit these once, keep this file in the directory you submit jobs from,
# and run `source config.sh` before running any R or Julia script by hand.
# Every shell/SLURM script sources it with:  source "${SLURM_SUBMIT_DIR:-$PWD}/config.sh"

# Root of the analysis working area (rebuilt genotype files, TrajGWAS output, phenotype files).
# The pipeline writes to $WHI_DIR/10_Rebuild_v2.
export WHI_DIR=/path/to/WHI_analysis

# Folder for the phenotype files written by 02_phenotype/build_pheno_files.R
# (01_Black/pheno_black_hipbmd_data.csv, ...) and read by the TrajGWAS step.
export PHENO_DIR=$WHI_DIR/02_phenotype

# Root of the WHI dbGaP downloads (one sub-folder per genotype study, plus 08_Phenotype).
export WHI_DBGAP_DIR=/path/to/dbGaP/WHI

# Root of the Framingham Heart Study look-up area (FHS genotype + phenotype files per BMD site).
export FHS_DIR=/path/to/FHS_analysis

# Folder holding reference files and tools: HRC.r1-1.GRCh37.wgs.mac5.sites.tab, liftOver,
# hg18ToHg19.over.chain.gz, and a Julia installation with TrajGWAS.
export TOOLS_DIR=/path/to/tools
