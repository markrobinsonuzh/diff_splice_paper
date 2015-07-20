## ----- compare_complete_incomplete
## <<compare_complete_incomplete.R>>

## Read input arguments
args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

print(path_to_counts)
print(output_txt_file)

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

## List subdirectories in count directory and define input bam files
dir_names <- list.dirs(path_to_counts, full.names = FALSE)
dir_names <- setdiff(dir_names, c("INCOMPLETE_MISSING20", 
                                  "INCOMPLETE_ATLEAST5", 
                                  "INCOMPLETE_ATLEAST10", 
                                  "INCOMPLETE_ATLEAST15", 
                                  "INCOMPLETE_ATLEAST25", 
                                  "rMATS", "cuffdiff",
                                  "INCOMPLETE_KALLISTOEST"))

T <- lapply(c(dir_names, paste0("INCOMPLETE_MISSING20/", dir_names, "_missing20"),
              paste0("INCOMPLETE_ATLEAST5/", dir_names, "_atleast5"),
              paste0("INCOMPLETE_ATLEAST10/", dir_names, "_atleast10"),
              paste0("INCOMPLETE_ATLEAST15/", dir_names, "_atleast15"),
              paste0("INCOMPLETE_ATLEAST25/", dir_names, "_atleast25"),
              paste0("INCOMPLETE_KALLISTOEST/", dir_names, "_kallistoest_atleast5"),
              paste0("INCOMPLETE_KALLISTOEST/", dir_names, "_kallistoest_atleast10"),
              paste0("INCOMPLETE_KALLISTOEST/", dir_names, "_kallistoest_atleast15"),
              paste0("INCOMPLETE_KALLISTOEST/", dir_names, "_kallistoest_atleast25")), 
            function(i) {
              print(i)
              f <- list.files(paste(path_to_counts, i, sep = "/"), pattern = "\\.txt$", 
                              full.names = TRUE, recursive = FALSE)
              if (length(f) > 0) {
                b <- t(sapply(f, function(fl) {
                  tmp <- read.delim(fl, header = FALSE, as.is = TRUE)
                  w <- grep("^_", tmp[, 1])
                  print(w)
                  if (length(w) > 0)
                    tmp <- tmp[-w, ]
                  c(nbrbin = nrow(tmp), totcount = sum(tmp[, 2]))
                }))
                print(dim(b))
                rownames(b) <- NULL
                b
              } else {
                data.frame(nbrbin = NA, totcount = NA)
              }
            })

names(T) <- c(paste0(dir_names, "_complete"), paste0(dir_names, "_incomplete"),
              paste0(dir_names, "_atleast5"), paste0(dir_names, "_atleast10"),
              paste0(dir_names, "_atleast15"), paste0(dir_names, "atleast25"),
              paste0(dir_names, "_kallistoest_atleast5"),
              paste0(dir_names, "_kallistoest_atleast10"),
              paste0(dir_names, "_kallistoest_atleast15"),
              paste0(dir_names, "_kallistoest_atleast25"))
S <- lapply(T, function(w) {
  c(nbrbin = mean(w[, "nbrbin"]), totcount = mean(w[, "totcount"]))
})
R <- do.call(rbind, S)
R <- round(R)

write.table(R, file = output_txt_file, row.names = TRUE, col.names = TRUE, 
            sep = "\t", quote = FALSE)