## ----- plot_dexseq_results
## <<plot_dexseq_results.R>>

## Generate informative plots relating gene characteristics to DEXSeq results
## Input: Rdata object with the results from DEXSeq (dxd, res and pgq objects)
## Input: output basename
## Input: truth table
args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

print(path_to_Rdata_object)
print(output_basename)
print(path_to_truth_table)
print(method_name)

library(DEXSeq)
library(dplyr)
library(Hmisc)

load(path_to_Rdata_object)

## If the object is an ExonCountSet (DEXSeq 1.8.0), change some of the 
## column headers so it works with the code below
if (class(dxd) == "ExonCountSet") {
  idx <- match(c("padjust", "geneID", "meanBase", "exonID"), colnames(res))
  colnames(res)[idx] <- c("padj", "groupID", "exonBaseMean", "featureID")
}

truth <- read.delim(path_to_truth_table, header = TRUE, as.is = TRUE)
ds <- truth$gene[which(truth$ds_status == 1)]
nonds <- truth$gene[which(truth$ds_status == 0)]

print(table(res$padj < 0.05, useNA = "ifany"))
print(table(pgq < 0.05, useNA = "ifany"))

fp <- intersect(nonds, names(pgq)[which(pgq <= 0.05)])
tp <- intersect(ds, names(pgq)[which(pgq <= 0.05)])
tn <- intersect(nonds, names(pgq)[which(pgq > 0.05)])
fn <- intersect(ds, names(pgq)[which(pgq > 0.05)])

## Mean-variance plot
if (class(dxd) == "DEXSeqDataSet") {
  conds <- data.frame(colData(dxd))[colData(dxd)$exon == "this", "condition"]
  pdf(paste0(output_basename, "_meanvar_cond1.pdf"))
  plot(apply(featureCounts(dxd)[, conds == levels(factor(conds))[1]], 1, mean),
       apply(featureCounts(dxd)[, conds == levels(factor(conds))[1]], 1, var),
       log = "xy", xlab = "Mean", ylab = "Variance", main = method_name)
  abline(0, 1, col = "red", lwd = 3)
  dev.off()

  ## Mean-dispersion estimate
  pdf(paste0(output_basename, "_meandispest.pdf"))
  DEXSeq::plotDispEsts(dxd)
  dev.off()
}

## Density of log2(dispersion)
dfp.disp <- density(log2(res$dispersion[res$groupID %in% fp]), na.rm = TRUE)
dtn.disp <- density(log2(res$dispersion[res$groupID %in% tn]), na.rm = TRUE)
dtp.disp <- density(log2(res$dispersion[res$groupID %in% tp]), na.rm = TRUE)
dfn.disp <- density(log2(res$dispersion[res$groupID %in% fn]), na.rm = TRUE)
pdf(paste0(output_basename, "_density_exonbindispersion.pdf"))
plot(dfp.disp, xlim = range(c(dfp.disp$x, dtn.disp$x, dtp.disp$x, dfn.disp$x)),
     ylim = range(c(dfp.disp$y, dtn.disp$y, dtp.disp$y, dfn.disp$y)), 
     col = "black", lwd = 2, xlab = "log2(exon bin dispersion)", 
     sub = "", main = method_name)
lines(dtn.disp, col = "blue", lwd = 2)
lines(dtp.disp, col = "green", lwd = 2)
lines(dfn.disp, col = "red", lwd = 2)
legend("topleft", legend = c("FP", "TN", "TP", "FN"),
       col = c("black", "blue", "green", "red"), 
       lty = 1, lwd = 2)
dev.off()

## Plot density of log2(exonBaseMean)
dfp.expr <- density(log2(res$exonBaseMean[res$groupID %in% fp] + 1), na.rm = TRUE, from = 0)
dtn.expr <- density(log2(res$exonBaseMean[res$groupID %in% tn] + 1), na.rm = TRUE, from = 0)
dtp.expr <- density(log2(res$exonBaseMean[res$groupID %in% tp] + 1), na.rm = TRUE, from = 0)
dfn.expr <- density(log2(res$exonBaseMean[res$groupID %in% fn] + 1), na.rm = TRUE, from = 0)
pdf(paste0(output_basename, "_density_exonbinexpression.pdf"))
plot(dfp.expr, xlim = range(c(dfp.expr$x, dtn.expr$x, dtp.expr$x, dfn.expr$x)),
     ylim = range(c(dfp.expr$y, dtn.expr$y, dtp.expr$y, dfn.expr$y)), 
     col = "black", lwd = 2, xlab = "log2(exon bin count + 1)", 
     sub = "", main = method_name)
lines(dtn.expr, col = "blue", lwd = 2)
lines(dtp.expr, col = "green", lwd = 2)
lines(dfn.expr, col = "red", lwd = 2)
legend("topright", legend = c("FP", "TN", "TP", "FN"),
       col = c("black", "blue", "green", "red"), 
       lty = 1, lwd = 2)
dev.off()

## Calculate the variance of log2(dispersion) values within each gene
res$genomicData <- NULL
res$transcripts <- NULL
G <- as.data.frame(res) %>% group_by(groupID) %>%
  summarise(vld = var(log2(dispersion), na.rm = TRUE),
            vex = var(log2(exonBaseMean + 1), na.rm = TRUE),
            nbrbin = length(dispersion),
            maxex = max(log2(exonBaseMean + 1), na.rm = TRUE),
            minex = min(log2(exonBaseMean + 1), na.rm = TRUE),
            rangex = maxex - minex)
G$nbrbinbin <- cut2(G$nbrbin, g = 6)
G$status <- NA
G$status[match(fp, G$groupID)] <- "FP"
G$status[match(tn, G$groupID)] <- "TN"
G$status[match(tp, G$groupID)] <- "TP"
G$status[match(fn, G$groupID)] <- "FN"
G$diff_IsoPct <- truth$diff_IsoPct[match(G$groupID, truth$gene)]

## Cut the variance(log2(dispersion))
G$vldbin <- cut2(G$vld, g = 6)

## Subset G
Gfp <- subset(G, status == "FP")
Gtn <- subset(G, status == "TN")
Gtp <- subset(G, status == "TP")
Gfn <- subset(G, status == "FN")

## Plot density of highest exon expression within gene
dfp.maxex <- density(Gfp$maxex, na.rm = TRUE)
dtn.maxex <- density(Gtn$maxex, na.rm = TRUE)
dtp.maxex <- density(Gtp$maxex, na.rm = TRUE)
dfn.maxex <- density(Gfn$maxex, na.rm = TRUE)
pdf(paste0(output_basename, "_density_max_exon_expression.pdf"))
plot(dfp.maxex,
     xlim = range(c(dfp.maxex$x, dtn.maxex$x, dtp.maxex$x, dfn.maxex$x)),
     ylim = range(c(dfp.maxex$y, dtn.maxex$y, dtp.maxex$y, dfn.maxex$y)),
     col = "black", lwd = 2,
     xlab = "Expression of most highly expressed bin",
     sub = "", main = method_name)
lines(dtn.maxex, col = "blue", lwd = 2)
lines(dtp.maxex, col = "green", lwd = 2)
lines(dfn.maxex, col = "red", lwd = 2)
legend("topright", legend = c("FP", "TN", "TP", "FN"),
       col = c("black", "blue", "green", "red"),
       lty = 1, lwd = 2)
dev.off()

## Plot density of lowest exon expression within gene
dfp.minex <- density(Gfp$minex, na.rm = TRUE)
dtn.minex <- density(Gtn$minex, na.rm = TRUE)
dtp.minex <- density(Gtp$minex, na.rm = TRUE)
dfn.minex <- density(Gfn$minex, na.rm = TRUE)
pdf(paste0(output_basename, "_density_min_exon_expression.pdf"))
plot(dfp.minex,
     xlim = range(c(dfp.minex$x, dtn.minex$x, dtp.minex$x, dfn.minex$x)),
     ylim = range(c(dfp.minex$y, dtn.minex$y, dtp.minex$y, dfn.minex$y)),
     col = "black", lwd = 2,
     xlab = "Expression of most lowly expressed bin",
     sub = "", main = method_name)
lines(dtn.minex, col = "blue", lwd = 2)
lines(dtp.minex, col = "green", lwd = 2)
lines(dfn.minex, col = "red", lwd = 2)
legend("topright", legend = c("FP", "TN", "TP", "FN"),
       col = c("black", "blue", "green", "red"),
       lty = 1, lwd = 2)
dev.off()

## Plot density of exon expression range within gene
dfp.rangex <- density(Gfp$rangex, na.rm = TRUE)
dtn.rangex <- density(Gtn$rangex, na.rm = TRUE)
dtp.rangex <- density(Gtp$rangex, na.rm = TRUE)
dfn.rangex <- density(Gfn$rangex, na.rm = TRUE)
pdf(paste0(output_basename, "_density_range_exon_expression.pdf"))
plot(dfp.rangex,
     xlim = range(c(dfp.rangex$x, dtn.rangex$x, dtp.rangex$x, dfn.rangex$x)),
     ylim = range(c(dfp.rangex$y, dtn.rangex$y, dtp.rangex$y, dfn.rangex$y)),
     col = "black", lwd = 2,
     xlab = "Range of bin expression values",
     sub = "", main = method_name)
lines(dtn.rangex, col = "blue", lwd = 2)
lines(dtp.rangex, col = "green", lwd = 2)
lines(dfn.rangex, col = "red", lwd = 2)
legend("topright", legend = c("FP", "TN", "TP", "FN"),
       col = c("black", "blue", "green", "red"),
       lty = 1, lwd = 2)
dev.off()

## Plot density of diff_isopct (the dominance of the top isoform)
dfp.isopct <- density(Gfp$diff_IsoPct, na.rm = TRUE, from = 0, to = 1)
dtn.isopct <- density(Gtn$diff_IsoPct, na.rm = TRUE, from = 0, to = 1)
dtp.isopct <- density(Gtp$diff_IsoPct, na.rm = TRUE, from = 0, to = 1)
dfn.isopct <- density(Gfn$diff_IsoPct, na.rm = TRUE, from = 0, to = 1)
pdf(paste0(output_basename, "_density_diff_isopct.pdf"))
plot(dfp.isopct, 
     xlim = range(c(dfp.isopct$x, dtn.isopct$x, dtp.isopct$x, dfn.isopct$x)),
     ylim = range(c(dfp.isopct$y, dtn.isopct$y, dtp.isopct$y, dfn.isopct$y)), 
     col = "black", lwd = 2, 
     xlab = "relative abundance difference between two most abundant isoforms", 
     sub = "", main = method_name)
lines(dtn.isopct, col = "blue", lwd = 2)
lines(dtp.isopct, col = "green", lwd = 2)
lines(dfn.isopct, col = "red", lwd = 2)
legend("topright", legend = c("FP", "TN", "TP", "FN"),
       col = c("black", "blue", "green", "red"), 
       lty = 1, lwd = 2)
dev.off()

## Plot density of variance(log2(dispersion))
dfp.dispvar <- density(Gfp$vld, na.rm = TRUE, from = 0)
dtn.dispvar <- density(Gtn$vld, na.rm = TRUE, from = 0)
dtp.dispvar <- density(Gtp$vld, na.rm = TRUE, from = 0)
dfn.dispvar <- density(Gfn$vld, na.rm = TRUE, from = 0)
pdf(paste0(output_basename, "_density_variance_exonbindispersion.pdf"))
plot(dfp.dispvar, 
     xlim = range(c(dfp.dispvar$x, dtn.dispvar$x, dtp.dispvar$x, dfn.dispvar$x)),
     ylim = range(c(dfp.dispvar$y, dtn.dispvar$y, dtp.dispvar$y, dfn.dispvar$y)), 
     col = "black", lwd = 2, 
     xlab = "within-gene variance(log2(exon bin dispersion))", 
     sub = "", main = method_name)
lines(dtn.dispvar, col = "blue", lwd = 2)
lines(dtp.dispvar, col = "green", lwd = 2)
lines(dfn.dispvar, col = "red", lwd = 2)
legend("topright", legend = c("FP", "TN", "TP", "FN"),
       col = c("black", "blue", "green", "red"), 
       lty = 1, lwd = 2)
dev.off()

## Plot density of variance(log expression)
dfp.dispvar <- density(Gfp$vex, na.rm = TRUE, from = 0)
dtn.dispvar <- density(Gtn$vex, na.rm = TRUE, from = 0)
dtp.dispvar <- density(Gtp$vex, na.rm = TRUE, from = 0)
dfn.dispvar <- density(Gfn$vex, na.rm = TRUE, from = 0)
pdf(paste0(output_basename, "_density_variance_exonbinexpression.pdf"))
plot(dfp.dispvar, 
     xlim = range(c(dfp.dispvar$x, dtn.dispvar$x, dtp.dispvar$x, dfn.dispvar$x)),
     ylim = range(c(dfp.dispvar$y, dtn.dispvar$y, dtp.dispvar$y, dfn.dispvar$y)), 
     col = "black", lwd = 2, 
     xlab = "within-gene variance(log2(exon bin count + 1))", 
     sub = "", main = method_name)
lines(dtn.dispvar, col = "blue", lwd = 2)
lines(dtp.dispvar, col = "green", lwd = 2)
lines(dfn.dispvar, col = "red", lwd = 2)
legend("topright", legend = c("FP", "TN", "TP", "FN"),
       col = c("black", "blue", "green", "red"), 
       lty = 1, lwd = 2)
dev.off()

## Plot density of the number of exon bins
dfp.nbrbin <- density(Gfp$nbrbin, na.rm = TRUE, from = 0)
dtn.nbrbin <- density(Gtn$nbrbin, na.rm = TRUE, from = 0)
dtp.nbrbin <- density(Gtp$nbrbin, na.rm = TRUE, from = 0)
dfn.nbrbin <- density(Gfn$nbrbin, na.rm = TRUE, from = 0)
pdf(paste0(output_basename, "_density_nbrexonbins.pdf"))
plot(dfp.nbrbin, 
     xlim = range(c(dfp.nbrbin$x, dtn.nbrbin$x, dtp.nbrbin$x, dfn.nbrbin$x)),
     ylim = range(c(dfp.nbrbin$y, dtn.nbrbin$y, dtp.nbrbin$y, dfn.nbrbin$y)), 
     col = "black", lwd = 2, 
     xlab = "number of exon bins", 
     sub = "", main = method_name)
lines(dtn.nbrbin, col = "blue", lwd = 2)
lines(dtp.nbrbin, col = "green", lwd = 2)
lines(dfn.nbrbin, col = "red", lwd = 2)
legend("topright", legend = c("FP", "TN", "TP", "FN"),
       col = c("black", "blue", "green", "red"), 
       lty = 1, lwd = 2)
dev.off()

## Plot number of exon bins vs variance in dispersion
pdf(paste0(output_basename, "_nbrexonbins_vs_variance_exonbindispersion.pdf"))
boxplot(Gfp$vld ~ Gfp$nbrbinbin, col = "grey", 
        xlab = "Number of exon bins", las = 2, 
        ylab = "within-gene variance(log2(dispersion))",
        main = paste0(method_name, ", FP (n = ", nrow(Gfp), ")"))
boxplot(Gtn$vld ~ Gtn$nbrbinbin, col = "lightblue", 
        xlab = "Number of exon bins", las = 2,  
        ylab = "within-gene variance(log2(dispersion))",
        main = paste0(method_name, ", TN (n = ", nrow(Gtn), ")"))
boxplot(Gtp$vld ~ Gtp$nbrbinbin, col = "lightgreen", 
        xlab = "Number of exon bins", las = 2,  
        ylab = "within-gene variance(log2(dispersion))",
        main = paste0(method_name, ", TP (n = ", nrow(Gtp), ")"))
boxplot(Gfn$vld ~ Gfn$nbrbinbin, col = "pink", 
        xlab = "Number of exon bins", las = 2,  
        ylab = "within-gene variance(log2(dispersion))",
        main = paste0(method_name, ", FN (n = ", nrow(Gfn), ")"))
dev.off()

## Plot status vs variance in log dispersion
counts <- table(G$status, G$vldbin)
pdf(paste0(output_basename, "_status_variance_exonbindispersion.pdf"))
par(mar = c(8, 4, 4, 4), mgp = c(6, 1, 0), xpd = TRUE)
barplot(counts, col = c("grey", "pink", "lightblue", "lightgreen"), 
        legend.text = rownames(counts), las = 2, 
        xlab = "within-gene variance(log2(dispersion))",
        args.legend = list(x = 8.3), main = method_name)
par(mar = c(5, 4, 4, 2), mgp = c(3, 1, 0), xpd = FALSE)
dev.off()

## sessionInfo
print(sessionInfo())


