## ----- rsem_summarise_run
## <<rsem_summarise_run.R>>

args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

print(path_to_rsem)
print(output_directory)
print(pattern)

library(dplyr)

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

## List subdirectories in RSEM directory and define input files
dir.names <- list.dirs(path_to_rsem, pattern = pattern, full.names = FALSE)
sample.names <- gsub("sample", "", dir.names)
isores.files <- paste0(path_to_rsem, "/", dir.names, "/", dir.names, 
                       ".isoforms.results")

keep <- which(file.exists(isores.files))
dir.names <- dir.names[keep]
sample.names <- sample.names[keep]
isores.files <- isores.files[keep]

print(isores.files)
print(sample.names)

if (length(isores.files) > 0) {
  L <- lapply(isores.files, function(w) {
    tmp <- read.delim(w, header = TRUE, as.is = TRUE)
    data.frame(tr = paste(tmp$gene_id, tmp$transcript_id, sep = ":"), 
               count = round(tmp$expected_count))
  })
}
names(L) <- sample.names
if (length(L) > 0) {
  cm <- L[[sample.names[1]]]
}
if (length(L) > 1) {
  for (i in 2:length(L)) {
    cm <- merge(cm, L[[sample.names[i]]], by = "tr", all = TRUE)
  }
}
cm[is.na(cm)] <- 0
colnames(cm) <- c("transcript", sample.names)

for (i in 2:ncol(cm)) {
  write.table(cm[, c(1, i)], file = paste0(output_directory, "/rsem", 
                                           sample.names[i - 1], ".txt"),
              row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)
}
