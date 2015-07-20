## ----- featureCounts_run
## <<featureCounts_run.R>>

## Read input arguments
args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

print(path_to_featurecounts_fcn)
print(gtf)
print(output.dir)
print(path_to_tophat)
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

## List subdirectories in tophat directory and define input bam files
dir.names <- list.dirs(path_to_tophat, pattern = pattern, full.names = FALSE)
sample.names <- gsub("sample", "", dir.names)
bam.files <- paste0(path_to_tophat, "/", dir.names, "/accepted_hits.bam")

keep <- which(file.exists(bam.files))
dir.names <- dir.names[keep]
sample.names <- sample.names[keep]
bam.files <- bam.files[keep]

print(bam.files)
print(sample.names)

## Run featureCounts
if (length(bam.files) > 0) {
  source(path_to_featurecounts_fcn)
  fcounts(gtf, bam.files, output.dir, sample.names)
}