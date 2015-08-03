## ----- tophatjunction_run
## <<tophatjunction_run.R>>

## Read input arguments
args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

print(gtf.file)
print(tophat.output.folder)
print(output.dir)
print(pattern)

## Help function to list subdirectories of a given directory
## From http://stackoverflow.com/questions/4749783/how-to-obtain-a-list-of-directories-within-a-directory-like-list-files-but-i
list.dirs <- function(path = ".", pattern = NULL, all.dirs = FALSE,
                      full.names = FALSE, ignore.case = FALSE) {
  all <- list.files(path, pattern, all.dirs,
                    full.names = TRUE, recursive = FALSE, ignore.case)
  dirs <- all[file.info(all)$isdir]
  if(isTRUE(full.names))
    return(dirs)
  else
    return(basename(dirs))
}

## List subdirectories in tophat directory and define input bed files
dir.names <- list.dirs(tophat.output.folder, pattern = pattern, full.names = FALSE)
sample.names <- gsub("sample", "", dir.names)
junction.files <- paste0(tophat.output.folder, "/", dir.names, "/junctions.bed")

keep <- which(file.exists(junction.files))
dir.names <- dir.names[keep]
sample.names <- sample.names[keep]
junction.files <- junction.files[keep]

print(junction.files)
print(sample.names)

## Generate junction read count files
if (length(junction.files) > 0) {
  library(GenomicRanges)
  library(GenomicFeatures)
  library(rtracklayer)
  library(dplyr)
  
  gtf <- makeTxDbFromGFF(gtf.file, format = "gtf")
  gtf2 <- import(gtf.file)
  
  ## Generate gene/transcript correspondence
  gene_tx_mapping <- mcols(gtf2)[, c("gene_id", "transcript_id")] %>% unique()
  
  ## Define junctions from gtf file
  tx <- exonsBy(gtf, "tx", use.names = TRUE)
  junctions <- lapply(tx, function(w) GenomicRanges::setdiff(range(w), w) + 1)
  junctions2 <- unlist(GRangesList(junctions))
  mcols(junctions2) <- data.frame(transcript = names(junctions2))
  names(junctions2) <- NULL
  junctions2 <- as.data.frame(junctions2)
  junctions2$gene <- gene_tx_mapping$gene_id[match(junctions2$transcript, 
                                                   gene_tx_mapping$transcript_id)]
  junctions2$position <- paste0(junctions2$seqnames, ":", junctions2$start,
                                "-", junctions2$end, junctions2$strand)
  junctions2 <- as.data.frame(junctions2  %>% distinct(position) %>% 
                                group_by(gene) %>% 
                                mutate(juncnr = 1:length(transcript)) %>% ungroup())
  junctions2$juncname <- paste(junctions2$gene, junctions2$juncnr, sep = ":")
  
  ## Extract number of mapped reads for each junction from TopHat files
  for (jf in 1:length(junction.files)) {
    tophat <- read.delim(junction.files[jf], skip = 1, header = FALSE, as.is = TRUE)
    tophat$add <- as.numeric(sapply(strsplit(tophat[, 11], ","), .subset, 1))
    tophat$subtr <- as.numeric(sapply(strsplit(tophat[, 11], ","), .subset, 2))
    tophat[, 2] <- tophat[, 2] + tophat$add
    tophat[, 3] <- tophat[, 3] - tophat$subtr + 1
    tophat$position <- paste0(tophat[, 1], ":", tophat[, 2], "-", tophat[, 3], tophat[, 6])
    
    final_table <- junctions2[, c("juncname", "position")]
    final_table$count <- tophat[match(final_table$position, tophat$position), 5]
    final_table[is.na(final_table)] <- 0
    
    write.table(final_table[, c("juncname", "count")], 
                file = paste0(output.dir, "/tophatjunc", sample.names[jf], ".txt"),
                row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)
  }
}