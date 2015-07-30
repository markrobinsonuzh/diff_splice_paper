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

gtf <- import(input_gtf)
gtf <- subset(gtf, type != "CDS")
gtf$type[gtf$type == "exon"] <- "CDS"
gtf$protein_id <- gtf$transcript_id
gtf$gene_name <- gtf$gene_id

gtf$exon_number <- as.character(gtf$exon_number)

export(object = gtf, con = output_gtf, format = "gtf")