## ----- generate_gtf_exontocds
## <<generate_gtf_exontocds.R>>

## Generate gtf file where all CDSs are taken out, and all exons are renamed as CDSs
args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

print(input_gtf)
print(output_gtf)

library(GenomicRanges)
library(rtracklayer)

## Import gtf
gtf <- import(input_gtf)

## Exclude CDS entries
gtf <- subset(gtf, type != "CDS")

## Rename exons as CDSs, and replace protein IDs with transcript IDs
gtf$type[gtf$type == "exon"] <- "CDS"
gtf$protein_id <- gtf$transcript_id

## Replace gene symbol with Ensembl ID
gtf$gene_name <- gtf$gene_id

## Keep exon number as character in output
gtf$exon_number <- as.character(gtf$exon_number)

export(object = gtf, con = output_gtf, format = "gtf")