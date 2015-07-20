## ----- dexseq_1.8_function
## <<dexseq_1.8_function.R>>

run_dexseq_1.8 <- function(count_files, conditions, flattened_file = NULL,
                           out_basename, method_name) {
  
  library(DEXSeq)
  
  sample_table <- data.frame(condition = conditions)
  
  dxd <- read.HTSeqCounts(count_files, 
                          design = sample_table)
  
  dxd <- estimateSizeFactors(dxd)
  dxd <- estimateDispersions(dxd, nCores = 4)
  dxd <- fitDispersionFunction(dxd)
  dxd <- testForDEU(dxd, nCores = 4)
  res <- DEUresultTable(dxd)
  pgq <- perGeneQValue(dxd, p = "pvalue")
  nm <- names(pgq)
  pgq <- pmin(1, pgq)
  names(pgq) <- nm
  
  save(dxd, res, pgq, file = paste0(out_basename, ".Rdata"))
  tmp <- cbind(gene = names(pgq), "adjP" = pgq)
  colnames(tmp)[colnames(tmp) == "adjP"] <- paste0(method_name, ":adjP")
  write.table(tmp, 
              file = paste0(out_basename, ".txt"), 
              col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
  
  print(sessionInfo())
}
