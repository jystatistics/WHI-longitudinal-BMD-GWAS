# 2026-11-25. Summarise genome-wide TrajGWAS output (summary statistics only; no individual-level data).
# Per set: completeness, QC counts, genomic inflation, top hits, SPTBN1 region, QQ + Manhattan plots.
suppressMessages({ library(data.table); library(ragg) })
B    <- file.path(Sys.getenv("WHI_DIR"), "10_Rebuild_v2/trajgwas")
OUT  <- file.path(B, "summary"); dir.create(OUT, showWarnings = FALSE)
sets <- c("black_hipbmd", "black_spinebmd", "hispanic_hipbmd", "hispanic_spinebmd", "ai_an_hipbmd")
MAFMIN <- 0.01; HWEMIN <- 1e-6
lam <- function(p, df) { p <- p[is.finite(p) & p > 0 & p <= 1]; median(qchisq(p, df, lower.tail = FALSE)) / qchisq(0.5, df) }
summ <- list(); tops <- list(); sptb <- list()
for (s in sets) {
  fs <- file.path(B, s, sprintf("chr%d.pval.txt", 1:22)); have <- file.exists(fs)
  cat(sprintf("== %s: %d/22 chromosome files%s\n", s, sum(have), if (all(have)) "" else paste0(" (missing chr ", paste(which(!have), collapse = ","), ")")))
  d <- rbindlist(lapply(fs[have], fread))
  for (v in c("maf", "hwepval", "betapval", "taupval", "jointpval")) set(d, j = v, value = suppressWarnings(as.numeric(d[[v]])))
  q <- d[maf >= MAFMIN & hwepval >= HWEMIN & is.finite(jointpval)]
  q5 <- q[maf >= 0.05]
  summ[[s]] <- data.table(set = s, variants_tested = nrow(d), after_QC = nrow(q),
    lambda_joint = lam(q$jointpval, 2), lambda_beta = lam(q$betapval, 1), lambda_tau = lam(q$taupval, 1),
    lambda_joint_maf05 = lam(q5$jointpval, 2), tau_p_zero = sum(q$taupval == 0, na.rm = TRUE),
    n_joint_5e8 = sum(q$jointpval < 5e-8), n_joint_1e6 = sum(q$jointpval < 1e-6),
    n_beta_5e8 = sum(q$betapval < 5e-8, na.rm = TRUE), n_tau_5e8 = sum(q$taupval < 5e-8, na.rm = TRUE))
  tops[[s]] <- cbind(set = s, head(q[order(jointpval)], 30))
  r <- q[chr == 2 & pos >= 54400000 & pos <= 55200000]   # SPTBN1 (hg19 chr2:54.68-54.90 Mb) +/- ~300 kb
  sptb[[s]] <- data.table(set = s, n_variants = nrow(r),
    min_joint = if (nrow(r)) min(r$jointpval) else NA, snp_min_joint = if (nrow(r)) r$snpid[which.min(r$jointpval)] else NA,
    min_beta = if (nrow(r)) min(r$betapval, na.rm = TRUE) else NA, min_tau = if (nrow(r)) min(r$taupval, na.rm = TRUE) else NA,
    rs4233949_present = any(d$chr == 2 & d$pos == 54659707))
  # plots: joint p; thin p > 0.01 to keep files small
  pl <- q[, .(chr, pos, p = jointpval)][p < 0.01 | runif(.N) < 0.05]
  setorder(pl, chr, pos); off <- pl[, .(m = max(pos)), by = chr][, cum := cumsum(as.numeric(m)) - m]
  pl <- off[pl, on = "chr"][, x := pos + cum]
  agg_png(file.path(OUT, paste0(s, "_manhattan.png")), width = 2400, height = 900, res = 200)
  plot(pl$x, -log10(pl$p), pch = 20, cex = 0.4, col = ifelse(pl$chr %% 2 == 1, "grey25", "steelblue"), xaxt = "n",
       xlab = "Chromosome", ylab = "-log10(joint P)", main = s)
  axis(1, at = off$cum + off$m / 2, labels = off$chr, cex.axis = 0.7); abline(h = -log10(c(5e-8, 1e-5)), col = c("red", "orange"), lty = 2)
  dev.off()
  agg_png(file.path(OUT, paste0(s, "_qq.png")), width = 1200, height = 1200, res = 200)
  o <- -log10(sort(q$jointpval)); e <- -log10(ppoints(length(o)))
  k <- c(which(o > 2), sample(which(o <= 2), min(20000, sum(o <= 2))))
  plot(e[k], o[k], pch = 20, cex = 0.4, xlab = "Expected -log10(P)", ylab = "Observed -log10(joint P)",
       main = sprintf("%s  lambda(joint) = %.3f", s, summ[[s]]$lambda_joint)); abline(0, 1, col = "red")
  dev.off()
  rm(d, q, q5, pl); gc()
}
S <- rbindlist(summ); TT <- rbindlist(tops, fill = TRUE); P <- rbindlist(sptb)
fwrite(S, file.path(OUT, "trajgwas_summary.csv")); fwrite(TT, file.path(OUT, "trajgwas_top30_joint.csv")); fwrite(P, file.path(OUT, "trajgwas_SPTBN1_region.csv"))
options(width = 250); cat(sprintf("\nQC: MAF >= %g, HWE P >= %g\n", MAFMIN, HWEMIN)); print(S, digits = 3); cat("\nSPTBN1 region:\n"); print(P, digits = 3)
cat("\nTop 5 joint hits per set:\n"); print(TT[, head(.SD, 5), by = set][, .(set, snpid, maf, betapval, taupval, jointpval)], digits = 3)
cat("\nWritten to", OUT, "\n")
