# 2026-11-25. Calibration check: genomic inflation by MAF bin for each test, and top hits for the separate
# (SPA-calibrated) mean (beta) and variability (tau) tests. Summary statistics only.
suppressMessages(library(data.table))
B <- file.path(Sys.getenv("WHI_DIR"), "10_Rebuild_v2/trajgwas")
sets <- c("black_hipbmd", "black_spinebmd", "hispanic_hipbmd", "hispanic_spinebmd", "ai_an_hipbmd")
lam <- function(p, df) { p <- p[is.finite(p) & p > 0 & p <= 1]; if (!length(p)) return(NA_real_); median(qchisq(p, df, lower.tail = FALSE)) / qchisq(0.5, df) }
L <- list(); H <- list()
for (s in sets) {
  d <- rbindlist(lapply(file.path(B, s, sprintf("chr%d.pval.txt", 1:22)), fread))
  for (v in c("maf", "hwepval", "betapval", "taupval", "jointpval")) set(d, j = v, value = suppressWarnings(as.numeric(d[[v]])))
  d <- d[maf >= 0.01 & hwepval >= 1e-6 & is.finite(jointpval)]
  d[, bin := cut(maf, c(0.01, 0.05, 0.2, 0.51), right = FALSE, labels = c("MAF .01-.05", "MAF .05-.20", "MAF >=.20"))]
  L[[s]] <- d[, .(n = .N, lam_beta = lam(betapval, 1), lam_tau = lam(taupval, 1), lam_joint = lam(jointpval, 2),
                  joint_vs_parts = mean(jointpval < pchisq(qchisq(betapval, 1, lower.tail = FALSE) + qchisq(taupval, 1, lower.tail = FALSE), 2, lower.tail = FALSE) / 100, na.rm = TRUE)),
              by = bin][order(bin)][, set := s]
  H[[s]] <- rbind(head(d[order(betapval)], 5)[, test := "beta (mean)"], head(d[order(taupval)], 5)[, test := "tau (variability)"])[, set := s]
  rm(d); gc()
}
options(width = 250)
cat("Genomic inflation by MAF bin (QC: MAF >= .01, HWE P >= 1e-6).\n")
cat("joint_vs_parts = share of SNPs whose joint P is >100x smaller than a 2-df test built from the beta and tau P (calibration check)\n")
print(rbindlist(L)[, .(set, bin, n, lam_beta, lam_tau, lam_joint, joint_vs_parts)], digits = 3)
cat("\nTop 5 SNPs by the separate SPA-calibrated tests:\n")
print(rbindlist(H)[, .(set, test, snpid, maf, betapval, taupval, jointpval)], digits = 3)
