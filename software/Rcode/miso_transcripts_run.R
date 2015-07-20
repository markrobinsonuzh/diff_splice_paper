## ----- miso_transcripts_run
## <<miso_transcripts_run.R>>

## Read input arguments
args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

print(path_to_miso_fcn)
print(miso_output_dir)
print(counts_output_dir)

## Summarise MISO data
source(path_to_miso_fcn)
miso_transcript(miso_output_dir, counts_output_dir)
