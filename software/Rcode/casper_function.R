## ----- casper_function
## <<casper_function.R>>

## Run casper to count reads aligning to exon paths
casper <- function(gtf, bamfiles, genome, read_length, output.dir, sample_names) {
  library(Rsamtools)
  library(GenomicFeatures)
  library(rtracklayer)
  library(edgeR)
  library(casper)
  library(dplyr)
  
  genDB <- makeTxDbFromGFF(gtf, format = "gtf")
  genomeDb <- procGenome(genDB = genDB, genome = genome, mc.cores = 10)
  
  ## Define which reads that will be imported from the bam files and 
  ## considered for the counting
  what <- scanBamWhat()
  what <- what[!(what %in% c("seq", "qual"))]
  flag <- scanBamFlag(isPaired = TRUE, hasUnmappedMate = FALSE)
  param <- ScanBamParam(flag = flag, what = what, tag = "XS")
  
  ## Create list of path counts for the bam files
  CTS <- list()
  for (fl in 1:length(bamfiles)) {
##  CTS <- lapply(1:length(bamfiles), function(fl) {
    ## Read bam file
    bam <- scanBam(file = bamfiles[fl], param = param)[[1]]
    bam <- rmShortInserts(bam, isizeMin = 50)
    pbam <- procBam(bam)
    
    ## Compute fragment start and length distributions
    distrs <- getDistrs(genomeDb, bam = bam, readLength = read_length)
    
    ## Compute counts for exon paths
    pc <- pathCounts(pbam, DB = genomeDb, mc.cores = 16)
    
    b <- as.data.frame(unlist(pc@counts[[1]]))
    d <- data.frame(rownames(b), b)
    colnames(d) <- c("path_number", "counts")
    d$path_number <- as.character(d$path_number)
    d$gene_number <- strsplit2(d$path_number, "[.]")[, 1]
    colnames(d)[3] <- "island_id"
    
    ## Calculate expression
    eset <- calcExp(distrs = distrs, genomeDB = genomeDb, pc = pc, 
                    readLength = read_length, rpkm = FALSE, mc.cores = 16)
    counts <- fData(eset)
    ## Merge the gene names of the same island
    counts <- data.frame(counts %>% group_by(island_id) %>% 
      summarise(gene_id = paste(unique(gene_id), collapse = "+")) %>% ungroup(), 
      stringsAsFactors = FALSE)
    counts <- merge(counts, d)
    ##counts  <- counts[!duplicated(counts$path_number), ]
    counts$path_id <- paste0(counts$gene_id, ":", counts$path_number)
    counts <- counts[, c("path_id", "counts")]
    colnames(counts)[colnames(counts) == "counts"] <- paste0("counts_", 
                                                             sample_names[fl])
    CTS[[fl]] <- counts
  }
  #)
  
  tmp <- CTS[[1]]
  if (length(CTS) > 1) {
    for (i in 2:length(CTS)) {
      tmp <- merge(tmp, CTS[[i]], all = TRUE)
    }
  }
  tmp[is.na(tmp)] <- 0
  
  for (i in 2:ncol(tmp)) {
    write.table(tmp[, c(1, i)], file = paste0(output.dir, "/casper", 
                                              gsub("counts_", "", 
                                                   colnames(tmp)[i]), ".txt"),
                col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
  }
  
  print(sessionInfo())
}
