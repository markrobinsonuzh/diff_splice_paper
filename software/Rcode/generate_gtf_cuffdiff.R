## ----- generate_gtf_cuffdiff
## <<generate_gtf_cuffdiff.R>>

## Generate gtf file where the gene_name is replaced by gene_id (Ensembl ID)
args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

print(input_gtf)
print(output_gtf)

library(GenomicRanges)
library(rtracklayer)

gtf <- import(input_gtf)
gtf$gene_name <- gtf$gene_id

gtf$exon_number <- as.character(gtf$exon_number)

export(object = gtf, con = output_gtf, format = "gtf")