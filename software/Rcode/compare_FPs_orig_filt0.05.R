## Generate lists of 'new' and 'eliminated' FPs when filtering gtf, for plotting
library(DEXSeq)

basedir <- "/home/Shared/data/alt_splicing_simulations/Simulation5_Charlotte/hsapiens"
truth <- read.delim(paste0(basedir, "/no_diffexpression/non_null_simulation", 
                           "/3_truth/truth_human_non_null.txt"), 
                    header = TRUE, as.is = TRUE)
load(paste0(basedir, "/no_diffexpression/non_null_simulation/4_results/", 
            "dexseq_htseq_nomerge.Rdata"))
res_orig <- res
pgq_orig <- pgq
load(paste0(basedir, "/no_diffexpression/non_null_simulation/4_results/", 
            "INCOMPLETE_ATLEAST5/dexseq_htseq_nomerge_atleast5.Rdata"))
res_filt <- res
pgq_filt <- pgq

fp_filt <- intersect(names(pgq_filt)[pgq_filt < 0.05], 
                     truth$gene[truth$ds_status == 0])
fp_orig <- intersect(names(pgq_orig)[pgq_orig < 0.05], 
                     truth$gene[truth$ds_status == 0])

fp_filt_not_orig <- setdiff(fp_filt, fp_orig)
fp_orig_not_filt <- setdiff(fp_orig, fp_filt)

FP_filt_not_orig <- data.frame(gene = fp_filt_not_orig, 
                               green_orig = NA,
                               red_orig = NA,
                               green_filt = NA,
                               red_filt = NA,
                               stringsAsFactors = FALSE)
for (i in 1:length(FP_filt_not_orig$gene)) {
  g <- FP_filt_not_orig$gene[i]
  spv <- (subset(res_filt, groupID == g & pvalue < 0.05))$featureID
  sadj <- (subset(res_filt, groupID == g & padj < 0.05))$featureID
  spv <- setdiff(spv, sadj)
  FP_filt_not_orig[i, "green_filt"] <- paste(spv, collapse = ",")
  FP_filt_not_orig[i, "red_filt"] <- paste(sadj, collapse = ",")
  
  spv <- (subset(res_orig, groupID == g & pvalue < 0.05))$featureID
  sadj <- (subset(res_orig, groupID == g & padj < 0.05))$featureID
  spv <- setdiff(spv, sadj)
  FP_filt_not_orig[i, "green_orig"] <- paste(spv, collapse = ",")
  FP_filt_not_orig[i, "red_orig"] <- paste(sadj, collapse = ",")
}
FP_filt_not_orig[FP_filt_not_orig == ""] <- NA

write.table(FP_filt_not_orig, 
            file = paste0(basedir, "/figures/coverage_plots/FP_filt_0.05_not_orig.txt"), 
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")

FP_orig_not_filt <- data.frame(gene = fp_orig_not_filt, 
                               green_orig = NA,
                               red_orig = NA,
                               green_filt = NA,
                               red_filt = NA,
                               stringsAsFactors = FALSE)
for (i in 1:length(FP_orig_not_filt$gene)) {
  g <- FP_orig_not_filt$gene[i]
  spv <- (subset(res_filt, groupID == g & pvalue < 0.05))$featureID
  sadj <- (subset(res_filt, groupID == g & padj < 0.05))$featureID
  spv <- setdiff(spv, sadj)
  FP_orig_not_filt[i, "green_filt"] <- paste(spv, collapse = ",")
  FP_orig_not_filt[i, "red_filt"] <- paste(sadj, collapse = ",")
  
  spv <- (subset(res_orig, groupID == g & pvalue < 0.05))$featureID
  sadj <- (subset(res_orig, groupID == g & padj < 0.05))$featureID
  spv <- setdiff(spv, sadj)
  FP_orig_not_filt[i, "green_orig"] <- paste(spv, collapse = ",")
  FP_orig_not_filt[i, "red_orig"] <- paste(sadj, collapse = ",")
}
FP_orig_not_filt[FP_orig_not_filt == ""] <- NA

write.table(FP_orig_not_filt, 
            file = paste0(basedir, "/figures/coverage_plots/FP_orig_not_filt_0.05.txt"), 
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")






