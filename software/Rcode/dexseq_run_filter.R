## ----- dexseq_run
## <<dexseq_run.R>>

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
if (!is.null(filter_bin_count))
  filter_bin_count <- as.numeric(filter_bin_count)
print(filter_bin_count)
if (!is.null(filter_bin_variance))
  filter_bin_variance <- as.numeric(filter_bin_variance)
print(filter_bin_variance)
if (!is.null(filter_perc))
  filter_perc <- as.numeric(filter_perc)
print(filter_perc)
if (!is.null(filter_perc_perkb))
  filter_perc_perkb <- as.numeric(filter_perc_perkb)
print(filter_perc_perkb)
if (!is.null(filter_perc_pertx))
  filter_perc_pertx <- as.numeric(filter_perc_pertx)
print(filter_perc_pertx)
if (!is.null(filter_perc_perkb_pertx))
  filter_perc_perkb_pertx <- as.numeric(filter_perc_perkb_pertx)
print(filter_perc_perkb_pertx)

## Retrieve count files from the provided directory
count_files <- list.files(path_to_count_files, pattern = ".txt$", 
                          recursive = FALSE, full.names = TRUE)

print(count_files)

## Run DEXSeq on the count files
source(path_to_dexseq_fcn)
run_dexseq(count_files = count_files, conditions = conditions, 
           flattened_file = flattened_file,
           out_basename = out_basename, method_name = method_name, 
           filter_bin_count = filter_bin_count,
           filter_bin_variance = filter_bin_variance, 
           filter_perc = filter_perc, filter_perc_perkb = filter_perc_perkb, 
           filter_perc_pertx = filter_perc_pertx,
           filter_perc_perkb_pertx = filter_perc_perkb_pertx)


