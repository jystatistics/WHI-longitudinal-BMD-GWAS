# 2026-11-25. Genomic inflation (overall and by MAF) for the TrajGWAS output folders named on the command line.
suppressMessages(library(data.table))
B <- file.path(Sys.getenv("WHI_DIR"), "10_Rebuild_v2/trajgwas")
lam <- function(p, df) { p <- p[is.finite(p) & p > 0 & p <= 1]; if (!length(p)) return(NA_real_); median(qchisq(p, df, lower.tail = FALSE)) / qchisq(0.5, df) }
out <- list()
for (s in commandArgs(TRUE)) {
  fs <- file.path(B, s, sprintf("chr%d.pval.txt", 1:22)); fs <- fs[file.exists(fs)]
  d <- rbindlist(lapply(fs, fread))
  for (v in c("maf", "hwepval", "betapval", "taupval", "jointpval")) set(d, j = v, value = suppressWarnings(as.numeric(d[[v]])))
  d <- d[maf >= 0.01 & hwepval >= 1e-6 & is.finite(jointpval)]
  d[, bin := cut(maf, c(0.01, 0.05, 0.2, 0.51), right = FALSE, labels = c("MAF .01-.05", "MAF .05-.20", "MAF >=.20"))]
  out[[s]] <- rbind(d[, .(bin = "all", n = .N, lam_beta = lam(betapval, 1), lam_tau = lam(taupval, 1), lam_joint = lam(jointpval, 2),
                          min_beta = min(betapval, na.rm = TRUE), min_tau = min(taupval, na.rm = TRUE))],
                    d[, .(n = .N, lam_beta = lam(betapval, 1), lam_tau = lam(taupval, 1), lam_joint = lam(jointpval, 2),
                          min_beta = min(betapval, na.rm = TRUE), min_tau = min(taupval, na.rm = TRUE)), by = bin][order(bin)][, bin := as.character(bin)])[, `:=`(set = s, chr_files = length(fs))]
  rm(d); gc()
}
options(width = 250); print(rbindlist(out, use.names = TRUE)[, .(set, chr_files, bin, n, lam_beta, lam_tau, lam_joint, min_beta, min_tau)], digits = 3)
