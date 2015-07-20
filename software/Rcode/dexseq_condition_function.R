## ----- dexseq_condition_function
## <<dexseq_condition_function.R>>

## Run DEXSeq on a given set of count files
run_dexseq_condition <- function(count_files, conditions, flattened_file = NULL,
                                 out_basename, method_name) {
  
  library(DEXSeq)
  BPPARAM <- MulticoreParam(workers = 4)
  
  ## Generate sample table
  sample_table <- data.frame(sample = 1:length(conditions),
                             condition = conditions)
  
  ## Genreate DEXSeqDataSet
  dxd <- DEXSeqDataSetFromHTSeq(count_files, 
                                sampleData = sample_table,
                                design = ~ condition + exon + condition:exon,
                                flattenedfile = flattened_file)
  
  ## Estimate size factors and dispersions
  dxd <- estimateSizeFactors(dxd)
  dxd <- estimateDispersions(dxd, formula = ~ condition + exon + condition:exon,
                             BPPARAM = BPPARAM)
  
  ## Run the test and get results for exon bins and genes
  dxd <- testForDEU(dxd, fullModel = ~ condition + exon + condition:exon,
                    reducedModel = ~ condition + exon, BPPARAM = BPPARAM)
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
