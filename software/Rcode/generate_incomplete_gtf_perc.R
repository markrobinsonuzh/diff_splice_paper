## ----- generate_incomplete_gtf
## <<generate_incomplete_gtf.R>>

## Read input arguments
args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

print(path_to_gtf)
print(path_to_truth)
print(min_abundance_perc)
print(output_gtf)
print(output_missingtx)

min_abundance_perc <- as.numeric(min_abundance_perc)

library(rtracklayer)

## Read gtf file and truth
gtf_file <- import(path_to_gtf)
truth <- read.delim(path_to_truth, header = TRUE, as.is = TRUE)

## Determine which transcripts to exclude
tx_remove <- truth$transcript_id[intersect(which(truth$IsoPct < min_abundance_perc),
                                           which(truth$IsoPct2 < min_abundance_perc))]

## Subset gtf file
new_gtf <- gtf_file[-which(gtf_file$transcript_id %in% tx_remove), ]

new_gtf$exon_number <- as.character(new_gtf$exon_number)

## Write to files
write.table(data.frame(transcript = tx_remove, 
                       status = truth$transcript_ds_status[match(tx_remove, 
                                                                 truth$transcript_id)]),
            file = output_missingtx, row.names = FALSE, col.names = TRUE,
            quote = FALSE, sep = "\t")

export(new_gtf, output_gtf)
