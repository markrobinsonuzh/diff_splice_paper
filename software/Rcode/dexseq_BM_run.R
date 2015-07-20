## ----- dexseq_BM_run
## <<dexseq_BM_run.R>>

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

count_files <- list.files(path_to_count_files, pattern = ".txt$", 
                          recursive = FALSE, full.names = TRUE)

source(path_to_dexseq_fcn)
run_dexseq_BM(count_files, conditions, flattened_file,
              out_basename, method_name)