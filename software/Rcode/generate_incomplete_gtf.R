## ----- generate_incomplete_gtf
## <<generate_incomplete_gtf.R>>

## Read input arguments
args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

print(path_to_gtf)
print(path_to_truth)
print(remove_fraction)
print(output_gtf)
print(seed)
print(output_missingtx)

set.seed(seed)

library(rtracklayer)

## Read gtf file and truth
gtf_file <- import(path_to_gtf)
truth <- read.delim(path_to_truth, header = TRUE, as.is = TRUE)

## Determine which transcripts to exclude
diffspliced <- truth$transcript_id[truth$transcript_ds_status == 1]
notdiffspliced <- truth$transcript_id[truth$transcript_ds_status == 0]
tx_remove <- c(sample(diffspliced, length(diffspliced) * remove_fraction),
               sample(notdiffspliced, length(notdiffspliced) * remove_fraction))

## Subset gtf file
new_gtf <- gtf_file[-which(gtf_file$transcript_id %in% tx_remove), ]

## Write to files
write.table(data.frame(transcript = tx_remove, 
                       status = truth$transcript_ds_status[match(tx_remove, 
                                                                 truth$transcript_id)]),
            file = output_missingtx, row.names = FALSE, col.names = TRUE,
            quote = FALSE, sep = "\t")

export(new_gtf, output_gtf)
