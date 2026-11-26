# 2026-11-25. FHS look-up, step 2: drop relatives (FHS is family based; TrajGWAS assumes independent women).
# PLINK --genome on LD-pruned SNPs; pairs with PI_HAT >= 0.1875 (2nd degree or closer); greedily drop the woman in most pairs.
# Args: keep list, .genome file, output keep list.
a <- commandArgs(TRUE); suppressMessages(library(data.table))
k <- fread(a[1], header = FALSE, colClasses = "character"); setnames(k, c("FID", "IID"))
g <- if (file.exists(a[2]) && file.size(a[2]) > 0) fread(a[2], colClasses = list(character = c("FID1", "IID1", "FID2", "IID2"))) else data.table()
p <- if (nrow(g)) g[PI_HAT >= 0.1875, .(a = IID1, b = IID2)] else data.table(a = character(), b = character())
cat("women", nrow(k), "| related pairs (PI_HAT >= 0.1875):", nrow(p), "\n")
drop <- character()
while (nrow(p)) {
  deg <- table(c(p$a, p$b)); w <- names(deg)[which.max(deg)]
  drop <- c(drop, w); p <- p[a != w & b != w]
}
u <- k[!IID %in% drop]
writeLines(paste(u$FID, u$IID), a[3]); cat("removed", length(drop), "-> unrelated women:", nrow(u), "\n")
