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

## Import gtf file
gtf <- import(input_gtf)
## Change gene symbol to Ensembl ID
gtf$gene_name <- gtf$gene_id

## Keep exon number as character in output
gtf$exon_number <- as.character(gtf$exon_number)

export(object = gtf, con = output_gtf, format = "gtf")