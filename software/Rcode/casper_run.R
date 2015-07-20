## ----- casper_run
## <<casper_run.R>>

## Read input arguments
args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

print(path_to_casper_fcn)
print(path_to_gtf)
print(path_to_tophat)
print(pattern)
print(genome)
print(read_length)
print(output.dir)

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
dir_names <- list.dirs(path_to_tophat, pattern = pattern, full.names = FALSE)
sample_names <- gsub("sample", "", dir_names)
bam_files <- paste0(path_to_tophat, "/", dir_names, "/accepted_hits.bam")

keep <- which(file.exists(bam_files))
dir_names <- dir_names[keep]
sample_names <- sample_names[keep]
bam_files <- bam_files[keep]

print(bam_files)
print(sample_names)

## Run casper to count reads for paths
if (length(bam_files) > 0) {
  source(path_to_casper_fcn)
  casper(gtf = path_to_gtf, bamfiles = bam_files, genome = genome,
         read_length = read_length, output.dir = output.dir,
         sample_names = sample_names)
}