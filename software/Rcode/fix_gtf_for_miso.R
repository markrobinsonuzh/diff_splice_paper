## ----- fix_gtf_for_miso
## <<fix_gtf_for_miso.R>>

## Read input arguments
args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

print(input_gtf)
print(output_gtf)

library(rtracklayer)

## Convert exon number column to characters, otherwise they are left out 
## in the conversion to gff3 by the perl script
gtf <- import(input_gtf)
gtf$exon_number <- as.character(gtf$exon_number)
export(gtf, output_gtf, format = "gtf")