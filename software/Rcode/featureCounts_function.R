## ----- featureCounts_function
## <<featureCounts_function.R>>

## Count reads mapping to each exon using featureCounts (Rsubread package)
fcounts <- function(gtf, bamfiles, output.dir, sample.names) {
  library(Rsubread)
  
  fc <- featureCounts(bamfiles, annot.ext = gtf, 
                      isGTFAnnotationFile = TRUE, GTF.featureType = "exon", 
                      GTF.attrType = "exon_id", useMetaFeatures = FALSE, 
                      isPairedEnd = TRUE, allowMultiOverlap = TRUE, 
                      strandSpecific = 0, requireBothEndsMapped = TRUE)
  fc$counts <- cbind(rownames(fc$counts), fc$counts)
  print(head(fc$counts))
  
  for (i in 2:ncol(fc$counts)) {
    write.table(fc$counts[, c(1, i)], paste0(output.dir, "/featureCounts",
                                             sample.names[i - 1], ".txt"),
                row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)
  }
}
