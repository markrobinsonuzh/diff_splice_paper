## ----- generate_truth_table_run
## <<generate_truth_table_run.R>>

## Read input arguments
args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

print(path_to_generate_truth_file)
print(path_to_final_summary)
print(out.file)
print(astalavista.file)
print(gtf.file)
print(flattened.gtf.file)
print(missing.annot.file)

## Read simulation summary file
final_summary <- read.delim(path_to_final_summary, header = TRUE, as.is = TRUE)

## Generate truth file
source(path_to_generate_truth_file)
generate_truth_table(final_summary = final_summary, out.file = out.file, 
                     astalavista.file = astalavista.file,
                     gtf.file = gtf.file, flattened.gtf.file = flattened.gtf.file,
                     missing.annot.file = missing.annot.file) 