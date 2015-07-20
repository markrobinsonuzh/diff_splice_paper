## ----- splicinggraph_function
## <<splicinggraph_function.R>>

## Count reads with SplicingGraph
splicinggraph <- function(gtf.file, output.dir, bamfiles, samplenames) {
  library(SplicingGraphs)
  library(Rsamtools)
  library(GenomicAlignments)
  
  ## Generate splicing graphs
  txdb  <- makeTxDbFromGFF(gtf.file, format = "gtf")
  
  sg <- SplicingGraphs(txdb)
  
  ## Define which reads are going to be imported from bam files
  flag0 <- scanBamFlag(isSecondaryAlignment = FALSE,
                       isNotPassingQualityControls = FALSE,
                       isDuplicate = FALSE)
  param0 <- ScanBamParam(flag = flag0)
  
  ## Read bam files and count reads
  for (i in seq_along(bamfiles)) {
    bam_file <- bamfiles[i]
    cat("Processing run ", names(bam_file), " ... ", sep = "")
    galp <- readGAlignmentPairs(bam_file, use.names = TRUE, param = param0)
    sg1 <- assignReads(sg, reads = galp, sample.name = paste0("sample", i))
    cat("OK\n")
    sg_counts <- countReads(sg1)
    sg_counts <- sg_counts[, -which(colnames(sg_counts) == "ex_or_in")]
    sg_counts <- as.data.frame(sg_counts)
    write.table(sg_counts, paste0(output.dir, "/splicinggraph", samplenames[i], ".txt"),
                row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
  }
  
  print(sessionInfo())
}