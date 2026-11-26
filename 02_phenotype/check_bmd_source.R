# 2026-11-24. Which WHI BMD variable did the Aim 2 phenotype files use? Prints summary numbers only.
# Compares hipbmd / spinebmd in the old pheno files with every BMD column of the dbGaP tables
# bmdhip_rel2 / bmdspine_rel2, matched on SUBJID and exam year (time == HIPVY / SPNVY).
suppressPackageStartupMessages(library(data.table))
P <- file.path(Sys.getenv("WHI_DBGAP_DIR"), "08_Phenotype")
W <- Sys.getenv("WHI_DIR")
rdg <- function(pat) rbindlist(lapply(list.files(P, pattern = pat, full.names = TRUE),
                                      fread, skip = "dbGaP_Subject_ID"), fill = TRUE)
ref <- list(hip = rdg("bmdhip_rel2.*txt$"), spine = rdg("bmdspine_rel2.*txt$"))
vy  <- c(hip = "HIPVY", spine = "SPNVY")
for (s in names(ref)) cat(s, "table:", nrow(ref[[s]]), "rows; BMD columns:",
                          grep("BMD", names(ref[[s]]), value = TRUE), "\n")
grps <- c(black = "01_Black/pheno_black", hispanic = "03_Hispanic/pheno_hispanic", ai_an = "04_AI_AN/pheno_ai")
for (g in names(grps)) for (s in c("hip", "spine")) {
  f <- file.path(W, paste0(grps[g], "_", s, "bmd_data.csv"))
  if (!file.exists(f)) { cat("\n", g, s, ": file not found\n"); next }
  old <- fread(f)[, .(SUBJID = as.character(SUBJID), time = as.numeric(time), y = as.numeric(get(paste0(s, "bmd"))))][!is.na(y)]
  r <- copy(ref[[s]]); r[, SUBJID := as.character(SUBJID)]; setnames(r, vy[s], "VY")
  m <- merge(old, r, by.x = c("SUBJID", "time"), by.y = c("SUBJID", "VY"))
  cat(sprintf("\n== %s %s: %d old visits, %d matched to dbGaP on SUBJID + exam year\n", g, s, nrow(old), nrow(m)))
  for (c in grep("BMD", names(r), value = TRUE)) {
    x <- suppressWarnings(as.numeric(m[[c]])); ok <- !is.na(x) & x > 0
    if (sum(ok) < 10) next
    cat(sprintf("  %-12s n=%5d  cor=%.4f  mean(old-new)=%+.4f  identical(|diff|<0.0005)=%5.1f%%\n",
                c, sum(ok), cor(m$y[ok], x[ok]), mean(m$y[ok] - x[ok]), 100 * mean(abs(m$y[ok] - x[ok]) < 5e-4)))
  }
}
