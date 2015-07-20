## ----- kallisto_summarise_run
## <<kallisto_summarise_run.R>>

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
print(output_dir)

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

## List subdirectories in tophat directory and define input bam files
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
    fl$est_counts <- round(fl$est_counts)
    fl
    })
  X <- ABUNDANCE[[1]][, c("genetr", "est_counts")]
  if (length(ABUNDANCE) > 1) {
    for (i in 2:length(ABUNDANCE)) {
      X <- merge(X, ABUNDANCE[[i]][, c("genetr", "est_counts")], by = "genetr", 
                 all = TRUE)
    }
  }
  colnames(X) <- c("gene", sample_names)
  for (i in 2:ncol(X)) {
    write.table(X[, c(1, i)], 
                file = paste0(output_dir, "/kallisto", sample_names[i - 1], ".txt"),
                row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)
  }
}


