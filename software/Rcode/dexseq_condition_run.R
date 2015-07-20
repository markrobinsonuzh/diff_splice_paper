## ----- dexseq_condition_run
## <<dexseq_condition_run.R>>

## Read input arguments
args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

print(path_to_dexseq_fcn)
print(path_to_count_files)
print(conditions)
print(flattened_file)
print(out_basename)
print(method_name)

## Retrieve count files from the provided directory
count_files <- list.files(path_to_count_files, pattern = ".txt$", 
                          recursive = FALSE, full.names = TRUE)

print(count_files)

## Run DEXSeq on the count files
source(path_to_dexseq_fcn)
run_dexseq_condition(count_files, conditions, flattened_file,
                     out_basename, method_name)