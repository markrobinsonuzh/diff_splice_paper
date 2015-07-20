## ----- subset_miso_counts_assignable
## <<subset_miso_counts_assignable.R>>

## Read input arguments
args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

print(path_to_input_files)
print(path_to_output_files)

## Exclude all rows corresponding to reads not compatible with any isoforms

input_files <- list.files(path_to_input_files, pattern = "\\.txt$", 
                          full.names = FALSE)
for (fl in input_files) {
  x <- read.delim(paste0(path_to_input_files, "/", fl), header = FALSE, 
                  as.is = TRUE)
  comb <- sapply(strsplit(x[, 1], ":"), .subset, 2)
  keep <- which(gsub(",", "", gsub("0", "", comb)) != "")
  x <- x[keep, ]
  write.table(x, file = paste0(path_to_output_files, "/", fl), quote = FALSE, 
              sep = "\t", row.names = FALSE, col.names = FALSE)
}