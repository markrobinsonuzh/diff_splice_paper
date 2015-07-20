## ----- plot_simulation_details
## <<plot_simulation_details.R>>

## Plot the TPMs of genes/isoforms across all pairs of samples, and color by ds/de status
args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

print(path_to_sim_details)
print(output_pdf)

library(dplyr)
library(scales)

details <- read.delim(path_to_sim_details, header = TRUE, as.is = TRUE)

pdf(output_pdf, width = 7, height = 14)
par(mfrow = c(4, 2))
idx <- grep("_isoformTPM$", colnames(details), value = TRUE)
for (i1 in 1:(length(idx) - 1)) {
	for (i2 in (i1 + 1):(length(idx))) {
		df <- details[, c(idx[i1], idx[i2], "gene_ds_status", "gene_de_status", 
		"transcript_ds_status", "gene_id")]
		df$gene_ds_de_status = paste0(df$gene_ds_status, ":", df$gene_de_status)
		df$transcript_ds_de_status = paste0(df$transcript_ds_status, ":", df$gene_de_status)
		
		plot(df[, idx[i1]], df[, idx[i2]], log = "xy", pch = ".", cex = 1, 
		col = factor(df$gene_de_status), xlab = idx[i1], ylab = idx[i2], 
		main = "TPM (transcript), colored by gene_de_status")
		legend("topleft", legend = levels(factor(df$gene_de_status)), 
		col = 1:nlevels(factor(df$gene_de_status)), pch = 19)
		
		plot(df[, idx[i1]], df[, idx[i2]], log = "xy", pch = ".", cex = 1, 
		col = factor(df$gene_ds_status), xlab = idx[i1], ylab = idx[i2],
		main = "TPM (transcript), colored by gene_ds_status")
		legend("topleft", legend = levels(factor(df$gene_ds_status)), 
		col = 1:nlevels(factor(df$gene_ds_status)), pch = 19)
		
		plot(df[, idx[i1]], df[, idx[i2]], log = "xy", pch = ".", cex = 1, 
		col = factor(df$transcript_ds_status), xlab = idx[i1], ylab = idx[i2],
		main = "TPM (transcript), colored by transcript_ds_status")
		legend("topleft", legend = levels(factor(df$transcript_ds_status)), 
		col = 1:nlevels(factor(df$transcript_ds_status)), pch = 19)
		
		plot(df[, idx[i1]], df[, idx[i2]], log = "xy", pch = ".", cex = 1, 
		col = factor(df$gene_ds_de_status), xlab = idx[i1], ylab = idx[i2],
		main = "TPM (transcript), colored by gene_ds_de_status")
		legend("topleft", legend = levels(factor(df$gene_ds_de_status)), 
		col = 1:nlevels(factor(df$gene_ds_de_status)), pch = 19)
		
		plot(df[, idx[i1]], df[, idx[i2]], log = "xy", pch = ".", cex = 1, 
		col = factor(df$transcript_ds_de_status), xlab = idx[i1], ylab = idx[i2],
		main = "TPM (transcript), colored by transcript_ds_de_status")
		legend("topleft", legend = levels(factor(df$transcript_ds_de_status)), 
		col = 1:nlevels(factor(df$transcript_ds_de_status)), pch = 19)
		
		df2 <- df %>% group_by(gene_id) %>% summarise_(tpm1 = paste0("sum(", idx[i1], ")"), 
		tpm2 = paste0("sum(", idx[i2], ")"), ds_status = "max(as.numeric(gene_ds_status))", 
		de_status = "max(as.numeric(gene_de_status))", 
		ds_de_status = "unique(gene_ds_de_status)")
		
		plot(df2$tpm1, df2$tpm2, log = "xy", pch = ".", cex = 1, 
		col = factor(df2$ds_status), xlab = idx[i1], ylab = idx[i2],
		main = "TPM (gene), colored by gene_ds_status")
		legend("topleft", legend = levels(factor(df2$ds_status)), 
		col = 1:nlevels(factor(df2$ds_status)), pch = 19)
		
		plot(df2$tpm1, df2$tpm2, log = "xy", pch = ".", cex = 1, 
		col = factor(df2$de_status), xlab = idx[i1], ylab = idx[i2],
		main = "TPM (gene), colored by gene_de_status")
		legend("topleft", legend = levels(factor(df2$de_status)), 
		col = 1:nlevels(factor(df2$de_status)), pch = 19)
		
		plot(df2$tpm1, df2$tpm2, log = "xy", pch = ".", cex = 1, 
		col = factor(df2$ds_de_status), xlab = idx[i1], ylab = idx[i2],
		main = "TPM (gene), colored by gene_ds_de_status")
		legend("topleft", legend = levels(factor(df2$ds_de_status)), 
		col = 1:nlevels(factor(df2$ds_de_status)), pch = 19)
	}
}
xdenods <- subset(details, gene_ds_status == 0 & gene_de_status == 1)
xdenods <- xdenods %>% mutate(geneCount_gr1 = 1/3 * (s1_geneCount + s2_geneCount
                                                     + s3_geneCount), 
                              geneCount_gr2 = 1/3 * (s4_geneCount + s5_geneCount 
                                                     + s6_geneCount),
                              isoformCount_gr1 = 1/3 * (s1_isoformCount + 
                                                          s2_isoformCount + 
                                                          s3_isoformCount), 
                              isoformCount_gr2 = 1/3 * (s4_isoformCount + 
                                                          s5_isoformCount + 
                                                          s6_isoformCount),
                              isoformTPM_gr1 = 1/3 * (s1_isoformTPM + 
                                                        s2_isoformTPM + 
                                                        s3_isoformTPM), 
                              isoformTPM_gr2 = 1/3 * (s4_isoformTPM + 
                                                        s5_isoformTPM + 
                                                        s6_isoformTPM))

xnodenods <- subset(details, gene_ds_status == 0 & gene_de_status == 0)
xnodenods <- xnodenods %>% mutate(geneCount_gr1 = 1/3 * (s1_geneCount + s2_geneCount
                                                         + s3_geneCount), 
                                  geneCount_gr2 = 1/3 * (s4_geneCount + s5_geneCount 
                                                         + s6_geneCount),
                                  isoformCount_gr1 = 1/3 * (s1_isoformCount + 
                                                              s2_isoformCount + 
                                                              s3_isoformCount), 
                                  isoformCount_gr2 = 1/3 * (s4_isoformCount + 
                                                              s5_isoformCount + 
                                                              s6_isoformCount),
                                  isoformTPM_gr1 = 1/3 * (s1_isoformTPM + 
                                                            s2_isoformTPM + 
                                                            s3_isoformTPM), 
                                  isoformTPM_gr2 = 1/3 * (s4_isoformTPM + 
                                                            s5_isoformTPM + 
                                                            s6_isoformTPM))

d1 <- density(log2(xdenods$geneCount_gr2/xdenods$geneCount_gr1), na.rm = TRUE)
d2 <- density(log2(xnodenods$geneCount_gr2/xnodenods$geneCount_gr1), na.rm = TRUE)
plot(d1, xlim = range(c(d1$x, d2$x)), ylim = range(c(d1$y, d2$y)), lwd = 2,
     main = "Average gene count per group", sub = "", 
     xlab = "log(group2average/group1average)")
lines(d2, col = "red", lwd = 2)
legend("topright", legend = c("DE, no DS", "no DE, no DS"), 
       col = c("black", "red"), lty = 1)

d1 <- density(log2(xdenods$isoformCount_gr2/xdenods$isoformCount_gr1), na.rm = TRUE)
d2 <- density(log2(xnodenods$isoformCount_gr2/xnodenods$isoformCount_gr1), na.rm = TRUE)
plot(d1, xlim = range(c(d1$x, d2$x)), ylim = range(c(d1$y, d2$y)), lwd = 2,
     main = "Average isoform count per group", sub = "", 
     xlab = "log(group2average/group1average)")
lines(d2, col = "red", lwd = 2)
legend("topright", legend = c("DE, no DS", "no DE, no DS"), 
       col = c("black", "red"), lty = 1)

d1 <- density(log2(xdenods$isoformTPM_gr2/xdenods$isoformTPM_gr1), na.rm = TRUE)
d2 <- density(log2(xnodenods$isoformTPM_gr2/xnodenods$isoformTPM_gr1), na.rm = TRUE)
plot(d1, xlim = range(c(d1$x, d2$x)), ylim = range(c(d1$y, d2$y)), lwd = 2,
     main = "Average isoform TPM per group", sub = "", 
     xlab = "log(group2average/group1average)")
lines(d2, col = "red", lwd = 2)
legend("topright", legend = c("DE, no DS", "no DE, no DS"), 
       col = c("black", "red"), lty = 1)

dev.off()