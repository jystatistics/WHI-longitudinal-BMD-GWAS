# 2026-11-25. Revised GWAS figures from the rebuilt TrajGWAS (summary statistics only).
#   Fig 1  hip   Manhattan: rows = mean (beta) / within-subject variability (tau); columns = Black, Hispanic, AI/AN
#   Fig 2  spine Manhattan: Black, Hispanic (AI/AN spine not analysed: null model could not be estimated)
#   Fig S3 / S4: matching QQ plots with genomic inflation (lambda)
# QC (group-specific, from 13_group_qc.R): MAF >= 0.01, MAC >= 20, HWE P >= 1e-6. P-values = TrajGWAS score test with saddlepoint approximation.
# Output: TIFF (LZW, 300 dpi, 7.5 in wide; PLOS) + PNG preview in trajgwas/figures/.
suppressMessages({ library(data.table); library(ragg) })
B   <- file.path(Sys.getenv("WHI_DIR"), "10_Rebuild_v2/trajgwas")
OUT <- file.path(B, "figures"); dir.create(OUT, showWarnings = FALSE)
lab <- c(black = "Black", hispanic = "Hispanic", ai_an = "AI/AN")
lam <- function(p) { p <- p[is.finite(p) & p > 0 & p <= 1]; median(qchisq(p, 1, lower.tail = FALSE)) / qchisq(0.5, 1) }
load_set <- function(s) {   # group-specific QC from 13_group_qc.R: MAF >= 0.01, MAC >= 20, HWE P >= 1e-6 within the group
  S  <- file.path(B, "summary")
  fq <- fread(file.path(S, paste0("grp_", s, ".frq")), select = c("SNP", "MAF", "NCHROBS"))
  hw <- fread(file.path(S, paste0("grp_", s, ".hwe")), select = c("SNP", "TEST", "P"))[TEST %in% c("ALL", "ALL(NP)"), .(SNP, P)]
  lf <- fread(file.path(S, paste0("leads_", s, ".csv")), nrows = 0)   # stops here if group QC did not finish
  ok <- fread(file.path(S, paste0("qcpass_", s, ".txt")), header = FALSE)$V1
  d <- rbindlist(lapply(file.path(B, s, sprintf("chr%d.pval.txt", 1:22)), fread, select = c("chr", "pos", "allele1", "allele2", "betapval", "taupval"), colClasses = list(character = c("allele1", "allele2"))))
  for (v in c("betapval", "taupval")) set(d, j = v, value = suppressWarnings(as.numeric(d[[v]])))
  d[paste(as.integer(chr), pos, pmin(allele1, allele2), pmax(allele1, allele2), sep = ":") %chin% ok]
}
manh <- function(d, col, main) {
  x <- d[is.finite(get(col)) & get(col) > 0, .(chr, pos, p = get(col))]
  set.seed(1); x <- x[p < 0.05 | runif(.N) < 0.1]                       # thin the uninformative bulk for file size
  off <- x[, .(m = max(pos)), by = chr][order(chr)][, cum := cumsum(as.numeric(m)) - m]
  x <- off[x, on = "chr"][, xx := pos + cum]
  ymax <- max(8, ceiling(-log10(min(x$p))) + 0.5)
  plot(x$xx, -log10(x$p), pch = 20, cex = 0.25, col = ifelse(x$chr %% 2 == 1, "#3B4A6B", "#8FA3C7"), xaxt = "n", xaxs = "i",
       ylim = c(0, ymax), xlab = "", ylab = expression(-log[10](italic(P))), main = main, cex.main = 0.85, cex.lab = 0.8, cex.axis = 0.7, las = 1)
  axis(1, at = off$cum + off$m / 2, labels = ifelse(off$chr %in% c(1:11, 13, 15, 18, 22), off$chr, ""), cex.axis = 0.55, tick = FALSE, line = -0.7, gap.axis = -1)
  mtext("Chromosome", side = 1, line = 1.1, cex = 0.6)
  abline(h = -log10(5e-8), col = "#C0392B", lty = 2, lwd = 0.8); abline(h = -log10(1e-5), col = "grey50", lty = 3, lwd = 0.8)
}
qq <- function(p, main) {
  p <- sort(p[is.finite(p) & p > 0 & p <= 1]); n <- length(p); e <- -log10(ppoints(n)); o <- -log10(p)
  k <- c(which(o > 2), round(seq(1, n, length.out = 20000))); k <- unique(k[k >= 1 & k <= n])
  i <- round(10^seq(0, log10(n), length.out = 300)); ci <- -log10(cbind(qbeta(0.025, i, n - i + 1), qbeta(0.975, i, n - i + 1)))
  plot(e[k], o[k], type = "n", xlab = expression(Expected ~ -log[10](italic(P))), ylab = expression(Observed ~ -log[10](italic(P))),
       main = main, cex.main = 0.85, cex.lab = 0.8, cex.axis = 0.7, las = 1, xlim = c(0, max(e)), ylim = c(0, max(max(o), max(e), 8)))
  polygon(c(-log10(i / (n + 1)), rev(-log10(i / (n + 1)))), c(ci[, 1], rev(ci[, 2])), col = "grey88", border = NA)
  points(e[k], o[k], pch = 20, cex = 0.3, col = "#3B4A6B"); abline(0, 1, col = "#C0392B", lwd = 0.8)
  legend("topleft", bty = "n", cex = 0.75, legend = as.expression(bquote(lambda[GC] == .(sprintf("%.2f", lam(p))))))
}
draw <- function(file, groups, outcome, fun) {
  for (dev in c("tiff", "png")) {
    f <- file.path(OUT, paste0(file, ".", dev)); w <- 7.5; h <- 5.2
    if (dev == "tiff") agg_tiff(f, width = w, height = h, units = "in", res = 300, compression = "lzw") else agg_png(f, width = w, height = h, units = "in", res = 150)
    par(mfrow = c(2, length(groups)), mar = c(3.0, 3.2, 1.8, 0.6), mgp = c(1.9, 0.5, 0), tcl = -0.25)
    for (row in c("betapval", "taupval")) for (g in groups) {
      d <- DATA[[paste0(g, "_", outcome)]]
      ttl <- sprintf("%s%s: %s", lab[[g]], if (g == "black") "" else " (exploratory)", if (row == "betapval") "mean BMD" else "BMD variability")
      fun(d, row, ttl)
    }
    dev.off()
  }
  cat("written", file, "\n")
}
sets <- c("black_hipbmd", "hispanic_hipbmd", "ai_an_hipbmd", "black_spinebmd", "hispanic_spinebmd")
DATA <- setNames(lapply(sets, load_set), sets)
draw("Fig1_Manhattan_hip",  c("black", "hispanic", "ai_an"), "hipbmd",   manh)
draw("Fig2_Manhattan_spine", c("black", "hispanic"),         "spinebmd", manh)
draw("FigS3_QQ_hip",        c("black", "hispanic", "ai_an"), "hipbmd",   function(d, col, main) qq(d[[col]], main))
draw("FigS4_QQ_spine",      c("black", "hispanic"),          "spinebmd", function(d, col, main) qq(d[[col]], main))
cat("\nFigures in", OUT, "\n")
