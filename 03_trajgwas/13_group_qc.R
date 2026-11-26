# 2026-11-25. Group-specific QC for the TrajGWAS results.
# The maf/hwepval columns written by TrajGWAS are computed over the whole genotype file (all 2,363 women), not the group
# analysed, so QC must use each set's own women: MAF >= 0.01, minor-allele count >= 20, HWE P >= 1e-6 (within group).
# Writes per set: summary/grp_<set>.frq/.hwe (plink), lambda table, lead SNPs (+/-500 kb) for the mean and variability
# tests, and the Wald list summary/wald_leads_<set>.txt. Summary statistics only.
suppressMessages(library(data.table))
R   <- Sys.getenv("WHI_DIR")
B   <- file.path(R, "10_Rebuild_v2"); TG <- file.path(B, "trajgwas"); S <- file.path(TG, "summary")
sets <- data.table(set = c("black_hipbmd", "black_spinebmd", "hispanic_hipbmd", "hispanic_spinebmd", "ai_an_hipbmd"),
                   pheno = file.path(R, c("01_Black/pheno_black_hipbmd_data.csv", "01_Black/pheno_black_spinebmd_data.csv",
                                          "03_Hispanic/pheno_hispanic_hipbmd_data.csv", "03_Hispanic/pheno_hispanic_spinebmd_data.csv",
                                          "04_AI_AN/pheno_ai_hipbmd_data.csv")),
                   outcome = c("hipbmd", "spinebmd", "hipbmd", "spinebmd", "hipbmd"))
fl  <- strsplit(trimws(readLines(file.path(B, "WHI_v5_r08.fam"))), "[[:space:]]+")   # read as text: no NA conversion of IDs
cat("fam lines:", length(fl), "| lines with < 6 fields:", sum(lengths(fl) < 6), "\n")
fam <- data.table(FID = vapply(fl, `[`, "", 1), IID = vapply(fl, `[`, "", 2))
bim <- fread(file.path(B, "WHI_v5_r08.bim"), header = FALSE, select = c(1, 2, 4, 5, 6), col.names = c("chr", "SNP", "pos", "a1", "a2"),
             colClasses = list(character = c(2, 5, 6)))
bim[, key := paste(as.integer(chr), pos, pmin(a1, a2), pmax(a1, a2), sep = ":")][, c("chr", "pos", "a1", "a2") := NULL]
al  <- fread(file.path(B, "alias_all.txt"), header = FALSE, colClasses = "character"); amap <- setNames(al$V2, al$V1)
resolve <- function(id) { for (k in 1:5) { j <- !is.na(amap[id]); if (!any(j)) break; id[j] <- amap[id[j]] }; id }
lam <- function(p, df = 1) { p <- p[is.finite(p) & p > 0 & p <= 1]; if (!length(p)) return(NA_real_); median(qchisq(p, df, lower.tail = FALSE)) / qchisq(0.5, df) }
lead <- function(d, col, k = 10, win = 5e5) {
  x <- d[is.finite(get(col))][order(get(col))][1:min(.N, 5000)]; keep <- integer(0)
  for (i in seq_len(nrow(x))) { if (length(keep) >= k) break
    if (!any(x$chr[keep] == x$chr[i] & abs(x$pos[keep] - x$pos[i]) < win)) keep <- c(keep, i) }
  x[keep]
}
L <- list(); N <- list()
for (r in seq_len(nrow(sets))) {
  s <- sets$set[r]; ph <- fread(sets$pheno[r]); ph[, IID := resolve(as.character(IID))]
  ph <- ph[is.finite(suppressWarnings(as.numeric(get(sets$outcome[r]))))]
  ids <- intersect(unique(ph$IID), fam$IID)
  kp <- fam[IID %in% ids]; stopifnot(nrow(kp) > 0, all(nzchar(kp$FID)), all(nzchar(kp$IID)))
  writeLines(paste(kp$FID, kp$IID), file.path(S, paste0("keep_", s, ".txt")))
  cat(s, ": pheno women", uniqueN(ph$IID), "| in genotype file", nrow(kp), "\n")
  rc <- system(sprintf("plink --bfile %s --keep %s --freq --hardy --out %s > /dev/null", file.path(B, "WHI_v5_r08"),
                 file.path(S, paste0("keep_", s, ".txt")), file.path(S, paste0("grp_", s))))
  if (rc != 0) stop("plink failed for ", s, "; see ", file.path(S, paste0("grp_", s, ".log")))
  fq <- fread(file.path(S, paste0("grp_", s, ".frq")), select = c("SNP", "MAF", "NCHROBS"))
  hw <- fread(file.path(S, paste0("grp_", s, ".hwe")), select = c("SNP", "TEST", "P"))[TEST %in% c("ALL", "ALL(NP)"), .(SNP, hwe_grp = P)]
  d <- rbindlist(lapply(file.path(TG, s, sprintf("chr%d.pval.txt", 1:22)), fread, select = c("chr", "pos", "snpid", "allele1", "allele2", "betapval", "betadir", "taupval", "taudir"), colClasses = list(character = c("snpid", "allele1", "allele2"))))
  for (v in c("betapval", "taupval")) set(d, j = v, value = suppressWarnings(as.numeric(d[[v]])))
  d[, snpid := as.character(snpid)]
  cat("  ID check: pval snpid e.g.", head(d$snpid, 2), "| frq SNP e.g.", head(fq$SNP, 2), "| exact ID matches", sum(d$snpid %in% fq$SNP), "\n")
  g <- hw[fq, on = "SNP"][bim, on = "SNP", nomatch = NULL][, .(key, maf_grp = MAF, NCHROBS, hwe_grp)]   # frq/hwe -> chr:pos:alleles
  d[, key := paste(as.integer(chr), pos, pmin(allele1, allele2), pmax(allele1, allele2), sep = ":")]
  g <- g[!duplicated(key) & !key %in% key[duplicated(key)]]
  d <- g[d, on = "key"]; cat("  matched by chr:pos:alleles", sum(!is.na(d$maf_grp)), "of", nrow(d), "\n")
  d[, mac := maf_grp * NCHROBS]
  cat("  non-missing: MAF", sum(!is.na(d$maf_grp)), "| HWE", sum(!is.na(d$hwe_grp)), "| MAF>=.01", sum(d$maf_grp >= 0.01, na.rm = TRUE), "| MAC>=20", sum(d$mac >= 20, na.rm = TRUE), "| HWE>=1e-6", sum(d$hwe_grp >= 1e-6, na.rm = TRUE), "\n")
  q <- d[maf_grp >= 0.01 & mac >= 20 & hwe_grp >= 1e-6]
  if (nrow(q) == 0) stop("no SNPs passed group QC for ", s, " - join problem, nothing written")
  q[, bin := cut(maf_grp, c(0, 0.05, 0.2, 0.51), right = FALSE, labels = c("MAF <.05", "MAF .05-.20", "MAF >=.20"))]
  N[[s]] <- data.table(set = s, women = length(ids), tested = nrow(d), pass_QC = nrow(q),
                       lam_beta = lam(q$betapval), lam_tau = lam(q$taupval),
                       n_beta_5e8 = sum(q$betapval < 5e-8, na.rm = TRUE), n_tau_5e8 = sum(q$taupval < 5e-8, na.rm = TRUE),
                       n_beta_1e5 = sum(q$betapval < 1e-5, na.rm = TRUE), n_tau_1e5 = sum(q$taupval < 1e-5, na.rm = TRUE),
                       min_beta = min(q$betapval, na.rm = TRUE), min_tau = min(q$taupval, na.rm = TRUE))
  L[[s]] <- q[, .(n = .N, lam_beta = lam(betapval), lam_tau = lam(taupval)), by = bin][order(bin)][, set := s]
  lb <- lead(q, "betapval")[, test := "mean (beta)"]; lt <- lead(q, "taupval")[, test := "variability (tau)"]
  ld <- rbind(lb, lt)[, set := s]
  fwrite(ld[, .(set, test, chr, pos, snpid, allele1, allele2, maf_grp, mac, hwe_grp, betapval, betadir, taupval, taudir)], file.path(S, paste0("leads_", s, ".csv")))
  writeLines(unique(ld$snpid), file.path(S, paste0("wald_leads_", s, ".txt")))
  fwrite(q[, .(key)], file.path(S, paste0("qcpass_", s, ".txt")), col.names = FALSE)
  r2 <- q[chr == 2 & pos >= 54400000 & pos <= 55200000]
  cat(sprintf("%s: %d women; SPTBN1 region after group QC: %d SNPs, min beta p %.3g, min tau p %.3g\n", s, length(ids), nrow(r2), min(r2$betapval, na.rm = TRUE), min(r2$taupval, na.rm = TRUE)))
  rm(d, q, fq, hw); gc()
}
options(width = 250)
NN <- rbindlist(N); fwrite(NN, file.path(S, "groupqc_summary.csv")); fwrite(rbindlist(L), file.path(S, "groupqc_lambda_by_maf.csv"))
cat("\nGroup-specific QC (MAF >= .01, MAC >= 20, HWE P >= 1e-6 within group):\n"); print(NN, digits = 3)
cat("\nLambda by group MAF:\n"); print(rbindlist(L)[, .(set, bin, n, lam_beta, lam_tau)], digits = 3)
cat("\nLead SNPs (+/-500 kb), top 3 per test:\n")
for (s in sets$set) print(fread(file.path(S, paste0("leads_", s, ".csv")))[, head(.SD, 3), by = test][, .(set, test, snpid, maf_grp, betapval, taupval)], digits = 3)
writeLines(format(Sys.time()), file.path(S, "groupqc_done.txt"))
