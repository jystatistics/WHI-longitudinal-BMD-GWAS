# 2026-11-25. Nearest gene (hg19 RefSeq, UCSC refGene) for the lead SNPs and rs4233949. Summary statistics only.
suppressMessages(library(data.table))
S <- file.path(Sys.getenv("WHI_DIR"), "10_Rebuild_v2/trajgwas/summary")
rf <- file.path(dirname(S), "refGene_hg19.txt.gz")
if (!file.exists(rf)) download.file("https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/refGene.txt.gz", rf, mode = "wb")
g <- fread(rf, header = FALSE, select = c(2, 3, 5, 6, 13), col.names = c("acc", "chrom", "start", "end", "gene"))
g <- g[grepl("^chr([0-9]+)$", chrom) & grepl("^NM_", acc)][, chr := as.integer(sub("chr", "", chrom))]
g <- g[, .(start = min(start) + 1L, end = max(end)), by = .(chr, gene)]           # one span per protein-coding gene
x <- fread(file.path(S, "lead_snps_wald.csv"), colClasses = list(character = c("snpid", "allele1", "allele2")))
x <- rbind(x, data.table(set = "SPTBN1", test = "rs4233949", chr = 2L, pos = 54659707L, snpid = "2:54659707"), fill = TRUE)
near <- function(c, p) {
  h <- g[chr == c]; inside <- h[start <= p & end >= p]
  if (nrow(inside)) return(paste(sort(unique(inside$gene)), collapse = "/"))
  d <- pmin(abs(h$start - p), abs(h$end - p)); i <- which.min(d)
  sprintf("%s (%s kb %s)", h$gene[i], round(d[i] / 1000), ifelse(h$start[i] > p, "upstream", "downstream"))
}
x[, nearest_gene := mapply(near, chr, pos)]
fwrite(x[set != "SPTBN1"], file.path(S, "lead_snps_wald_genes.csv"))
print(x[, .(set, test, snpid, nearest_gene)][order(set, test)], nrows = 200)
