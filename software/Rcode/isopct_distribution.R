## ----- isopct_distribution
## <<isopct_distribution.R>>

## Plot distribution of RSEM-estimated isoform percentages

args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

print(referencefile)
print(outdir)

bnm <- basename(referencefile)

estimates <- read.delim(referencefile, header = TRUE, as.is = TRUE)

m <- max(table(estimates$gene_id))

IsoPct <- t(sapply(unique(estimates$gene_id), function(w) {
  tmp <- estimates[estimates$gene_id == w, ]
  tmp2 <- rep(0, m)
  tmp2[1:nrow(tmp)] <- sort(tmp$IsoPct, decreasing = TRUE)
  tmp2
}))

pdf(paste0(outdir, "/", bnm, "_isoform_percentages.pdf"))
boxplot(IsoPct[rowSums(IsoPct) > 0, 1:5], 
        xlab = "Isoform number", ylab = "Percentage",
        main = "Estimated isoform percentages, expressed genes",
        cex.main = 0.9)
dev.off()

f <- function(n) {
  a <- exp(-n + 0.3^(n - 1))
  100*a/sum(a)
}

pdf(paste0(outdir, "/", bnm, "_isoform_percentages_per_ntranscript.pdf"), 
    width = 7, height = 10.5)
par(mfrow = c(3, 2))
for (i in 2:7) {
  boxplot(IsoPct[rowSums(IsoPct > 0) == i, 1:i], 
          xlab = "Isoform number", ylab = "Percentage",
          main = paste0("Estimated isoform percentages for genes with ", i, 
                        " expressed isoforms\n", 
                        "Red line = Previous estimate"),
          cex.main = 0.75)
  lines(1:i, f(1:i), col = "red", lwd = 2)
}
dev.off()

pdf(paste0(outdir, "/", bnm, "_isoform_percentages_per_ntranscript_noline.pdf"), 
    width = 7, height = 10.5)
par(mfrow = c(3, 2))
for (i in 2:7) {
  boxplot(IsoPct[rowSums(IsoPct > 0) == i, 1:i], 
          xlab = "Isoform number", ylab = "Percentage",
          main = paste0("Estimated isoform percentages for genes with ", i, 
                        " expressed isoforms\n"),
          cex.main = 0.75)
}
dev.off()

## Calculate differences between relative abundance of consecutive isoforms
pdf(paste0(outdir, "/", bnm, "_isoform_percentage_differences_per_ntranscript.pdf"), 
    width = 7, height = 10.5)
par(mfrow = c(3, 2))
for (i in 2:7) {
  tmp <- IsoPct[rowSums(IsoPct > 0) == i, 1:i]
  tmp <- -t(diff(t(tmp)))
  colnames(tmp) <- paste0((1:(i - 1)), "-", (2:i))
  boxplot(tmp[, 1:(i - 1)], xlab = "Isoform pair",
          ylab = "Relative abundance difference", 
          main = paste0("Relative abundance difference between consecutive\n",
                        "isoforms for genes with ", i, " expressed isoforms."),
          cex.main = 0.75)
}
dev.off()

pdf(paste0(outdir, "/", bnm, "_isoform_percentage_differences_per_ntranscript_min2w10.pdf"), 
    width = 7, height = 10.5)
par(mfrow = c(3, 2))
for (i in 2:7) {
  tmp <- IsoPct[rowSums(IsoPct > 0) == i & rowSums(IsoPct > 10) > 1, 1:i]
  tmp <- -t(diff(t(tmp)))
  colnames(tmp) <- paste0((1:(i - 1)), "-", (2:i))
  boxplot(tmp[, 1:(i - 1)], xlab = "Isoform pair",
          ylim = c(0, 100), 
          ylab = "Relative abundance difference", 
          main = paste0("Relative abundance difference between consecutive\n",
                        "isoforms for genes with ", i, " expressed isoforms\n", 
                        "(at least 2 with abundance fraction > 10%)"),
          cex.main = 0.75)
}
dev.off()