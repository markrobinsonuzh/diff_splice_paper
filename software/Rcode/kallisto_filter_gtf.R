## ----- kallisto_filter_gtf
## <<kallisto_filter_gtf.R>>

## Read input arguments
args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

library(dplyr)

print(path_to_kallisto_output)
print(path_to_conversion_table)
print(path_to_isoform_results)
print(pattern)
print(path_to_gtf)
print(min_abundance_perc)
print(output_gtf)
print(output_missingtx)

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

## List subdirectories in kallisto directory
dir_names <- list.dirs(path_to_kallisto_output, pattern = pattern, full.names = FALSE)
sample_names <- gsub("sample", "", dir_names)
abundance_files <- paste0(path_to_kallisto_output, "/", dir_names, "/abundance.txt")

keep <- which(file.exists(abundance_files))
dir_names <- dir_names[keep]
sample_names <- sample_names[keep]
abundance_files <- abundance_files[keep]

print(abundance_files)
print(sample_names)

conv <- read.delim(path_to_conversion_table, header = FALSE, as.is = TRUE, 
                   sep = " ")
colnames(conv) <- c("kallisto_id", "transcript_id")

genetranscript <- read.delim(path_to_isoform_results, header = TRUE, as.is = TRUE)

if (length(abundance_files) > 0) {
  ABUNDANCE <- lapply(abundance_files, function(af) {
    fl <- read.delim(af, header = TRUE, as.is = TRUE)
    fl$transcript <- conv$transcript_id[match(fl$target_id, conv$kallisto_id)]
    fl$gene <- genetranscript$gene_id[match(fl$transcript, genetranscript$transcript_id)]
    fl$genetr <- paste(fl$gene, fl$transcript, sep = ":")
    fl
  })
  X <- ABUNDANCE[[1]][, c("genetr", "tpm")]
  if (length(ABUNDANCE) > 1) {
    for (i in 2:length(ABUNDANCE)) {
      X <- merge(X, ABUNDANCE[[i]][, c("genetr", "tpm")], by = "genetr", 
                 all = TRUE)
    }
  }
  colnames(X) <- c("transcript", sample_names)
  X$gene <- sapply(strsplit(X$transcript, ":"), .subset, 1)
  X$transcript <- sapply(strsplit(X$transcript, ":"), .subset, 2)
  tmp <- X %>% group_by(gene) %>% mutate_each(funs(perc = ./sum(.)), -transcript)
  tmp[is.na(tmp)] <- 0
  tmp <- as.data.frame(tmp)
  tmp$gene <- NULL
  tmp$maxfrac <- apply(tmp[, sample_names], 1, max)
  
  keep_tr <- tmp$transcript[tmp$maxfrac >= min_abundance_perc/100]
  remove_tr <- tmp$transcript[tmp$maxfrac < min_abundance_perc/100]
  
  library(rtracklayer)
  gtf_file <- import(path_to_gtf)

  new_gtf <- gtf_file[-which(gtf_file$transcript_id %in% remove_tr), ]
  new_gtf$exon_number <- as.character(new_gtf$exon_number)
  
  ## Write to files
  write.table(data.frame(transcript = remove_tr, 
                         status = genetranscript$transcript_ds_status[match(remove_tr, 
                                                                            genetranscript$transcript_id)]),
              file = output_missingtx, row.names = FALSE, col.names = TRUE,
              quote = FALSE, sep = "\t")
  
  export(new_gtf, output_gtf)
  
}


