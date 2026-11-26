# Build the long-format WHI phenotype files used by TrajGWAS, directly from the dbGaP tables.
# One row per woman per DXA visit:  FID IID dbGaP_Subject_ID time <outcome> AGE HEIGHTX WEIGHTX RACE
#   time     = DXA visit year (HIPVY / SPNVY: 0 = baseline, then 1, 3, 6, 9 ...)
#   hipbmd   = HIPBMDCOR  (total hip BMD, scanner-corrected, g/cm2)   from bmdhip_rel2
#   spinebmd = SPNBMDCOR  (total spine BMD, scanner-corrected, g/cm2) from bmdspine_rel2
#   AGE, RACE         = baseline age and self-reported race/ethnicity       from f2_rel1
#   HEIGHTX, WEIGHTX  = baseline height (cm) and weight (kg), visit year 0  from f80_rel1
# Genotype IDs (fam IID) are linked to dbGaP_Subject_ID through the WHI Sample and Subject files (SAMPLE_ID, SUBJID, or
# SOURCE_SUBJECT_ID for GARNET).
# A woman with more than one genotype record keeps the first one in the fam file (duplicates are resolved later, 04b-04f).
# Usage (after `source config.sh`):  Rscript build_pheno_files.R [fam file] [output folder]
suppressPackageStartupMessages(library(data.table))
P   <- Sys.getenv("WHI_DBGAP_DIR"); stopifnot(nzchar(P))
args <- commandArgs(trailingOnly = TRUE)
FAM <- if (length(args) >= 1) args[1] else file.path(Sys.getenv("WHI_DIR"), "WHI_merged_all.fam")
OUT <- if (length(args) >= 2) args[2] else Sys.getenv("PHENO_DIR", file.path(Sys.getenv("WHI_DIR"), "02_phenotype"))
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

# read every dbGaP file whose name matches a pattern (all consent groups), header line starts with dbGaP_Subject_ID
rd <- function(pat, cols = NULL) {
  fs <- list.files(P, pattern = paste0(pat, ".*\\.txt(\\.gz)?$"), full.names = TRUE, recursive = TRUE)
  fs <- fs[!grepl("data_dict|var_report", fs)]
  fs <- fs[!(grepl("\\.gz$", fs) & sub("\\.gz$", "", fs) %in% fs)]   # plain .txt wins over a .gz copy
  if (!length(fs)) stop("no dbGaP file matching ", pat, " under ", P)
  d <- rbindlist(lapply(fs, function(f) fread(f, skip = "dbGaP_Subject_ID", sep = "\t", quote = "", colClasses = "character",
                                              select = cols)), fill = TRUE)
  cat(sprintf("%-22s %d file(s), %d rows\n", pat, length(fs), nrow(d))); unique(d)
}
num <- function(x) suppressWarnings(as.numeric(x))

# ---- baseline covariates
f2  <- rd("f2_rel1", c("dbGaP_Subject_ID", "AGE", "RACE"))[, .(dbGaP_Subject_ID, AGE = num(AGE), RACE = as.integer(num(RACE)))]
f2  <- f2[!duplicated(dbGaP_Subject_ID)]
f80 <- rd("f80_rel1", c("dbGaP_Subject_ID", "F80VY", "F80DAYS", "HEIGHTX", "WEIGHTX"))
f80 <- f80[num(F80VY) == 0][order(dbGaP_Subject_ID, num(F80DAYS))][!duplicated(dbGaP_Subject_ID),
          .(dbGaP_Subject_ID, HEIGHTX = num(HEIGHTX), WEIGHTX = num(WEIGHTX))]

# ---- repeated BMD (one value per woman and visit year; the earliest scan if a year has two)
hip <- rd("bmdhip_rel2", c("dbGaP_Subject_ID", "HIPVY", "HIPDAYS", "HIPBMDCOR"))[
         , .(dbGaP_Subject_ID, time = num(HIPVY), days = num(HIPDAYS), hipbmd = num(HIPBMDCOR))]
spn <- rd("bmdspine_rel2", c("dbGaP_Subject_ID", "SPNVY", "SPNDAYS", "SPNBMDCOR"))[
         , .(dbGaP_Subject_ID, time = num(SPNVY), days = num(SPNDAYS), spinebmd = num(SPNBMDCOR))]
first_scan <- function(d, y) {
  d <- d[!is.na(get(y)) & get(y) > 0 & !is.na(time)][order(dbGaP_Subject_ID, time, days)]
  n0 <- nrow(d); d <- d[!duplicated(d, by = c("dbGaP_Subject_ID", "time"))]
  cat(sprintf("%s: %d scans, %d repeat scans in the same visit year dropped\n", y, n0, n0 - nrow(d))); d[, days := NULL][]
}
hip <- first_scan(hip, "hipbmd"); spn <- first_scan(spn, "spinebmd")

# ---- genotype IID -> dbGaP_Subject_ID
idcol <- function(d, rx) grep(rx, names(d), ignore.case = TRUE, value = TRUE)[1]
smp <- rd("Sample\\.MULTI"); sc <- idcol(smp, "^SAMPLE_?ID$")
sbj <- rd("Subject\\.MULTI"); uc <- idcol(sbj, "^SUBJID$")
# GARNET genotype IDs are the SOURCE_SUBJECT_ID of the Subject file; the other studies use SAMPLE_ID or SUBJID
map <- unique(rbind(smp[, .(id = get(sc), dbGaP_Subject_ID)], sbj[, .(id = get(uc), dbGaP_Subject_ID)],
                    sbj[, .(id = SOURCE_SUBJECT_ID, dbGaP_Subject_ID)]))
map <- map[!is.na(id) & nzchar(id)][!duplicated(id)]
fam <- fread(FAM, header = FALSE, colClasses = "character")[, .(FID = V1, IID = V2)]
fam <- merge(fam[, ord := .I], map, by.x = "IID", by.y = "id", all.x = TRUE)[order(ord)][, ord := NULL]
cat(sprintf("genotype records: %d; linked to a dbGaP subject: %d; unique women: %d\n",
            nrow(fam), sum(!is.na(fam$dbGaP_Subject_ID)), uniqueN(na.omit(fam$dbGaP_Subject_ID))))
fam <- fam[!is.na(dbGaP_Subject_ID) & !duplicated(dbGaP_Subject_ID)]

base <- merge(merge(fam, f2, by = "dbGaP_Subject_ID", all.x = TRUE), f80, by = "dbGaP_Subject_ID", all.x = TRUE)
groups <- list(black = c(3L, "01_Black/pheno_black"), hispanic = c(4L, "03_Hispanic/pheno_hispanic"),
               ai_an = c(1L, "04_AI_AN/pheno_ai"))
for (g in names(groups)) for (s in c("hip", "spine")) {
  y <- paste0(s, "bmd"); bmd <- if (s == "hip") hip else spn
  d <- merge(base[RACE == as.integer(groups[[g]][1])], bmd, by = "dbGaP_Subject_ID")
  setcolorder(d, c("FID", "IID", "dbGaP_Subject_ID", "time", y, "AGE", "HEIGHTX", "WEIGHTX", "RACE"))
  setorder(d, IID, time)
  f <- file.path(OUT, paste0(groups[[g]][2], "_", y, "_data.csv")); dir.create(dirname(f), FALSE, TRUE)
  fwrite(d, f)
  cat(sprintf("%-9s %-8s: %5d women, %5d visits -> %s\n", g, y, uniqueN(d$IID), nrow(d), f))
}
