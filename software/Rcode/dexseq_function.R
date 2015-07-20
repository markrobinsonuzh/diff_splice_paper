## ----- dexseq_function
## <<dexseq_function.R>>

## Run DEXSeq on a given set of count files
run_dexseq <- function(count_files, conditions, flattened_file = NULL,
                       out_basename, method_name) {
  
  library(DEXSeq)
  BPPARAM <- MulticoreParam(workers = 6)
  
  ## Generate sample table
  sample_table <- data.frame(sample = 1:length(conditions),
                             condition = conditions)
  
  ## Generate DEXSeqDataSet
  message("Generating DEXSeqDataSet...")
  dxd <- DEXSeqDataSetFromHTSeq(count_files, 
                                sampleData = sample_table,
                                design = ~ sample + exon + condition:exon,
                                flattenedfile = flattened_file)
  
  ## Estimate size factors and dispersions
  dxd <- estimateSizeFactors(dxd)
  message("Estimating dispersions...")
  dxd <- estimateDispersions(dxd, BPPARAM = BPPARAM)
  
  ## Run the test and get results for exon bins and genes
  message("Running test...")
  dxd <- testForDEU(dxd, BPPARAM = BPPARAM)
  
  message("Summarizing results on gene level...")
  res <- DEXSeqResults(dxd)
  pgq <- perGeneQValue(res, p = "pvalue")
  
  ## Save results
  save(dxd, res, pgq, file = paste0(out_basename, ".Rdata"))
  tmp <- cbind(gene = names(pgq), "adjP" = pgq)
  colnames(tmp)[colnames(tmp) == "adjP"] <- paste0(method_name, ":adjP")
  write.table(tmp, 
              file = paste0(out_basename, ".txt"), 
              col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
  
  print(sessionInfo())
}
