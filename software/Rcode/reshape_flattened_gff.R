## ----- reshape_flattened_gff
## <<reshape_flattened_gff.R>>

## Reshape flattened gff file from DEXSeq_prepare_annotation so that it works 
## with featureCounts

args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

print(input_gff_file)
print(output_gff_file)

library(rtracklayer)

x <- as.data.frame(import(input_gff_file))
x$type <- as.character(x$type)
x$type[x$type == "exonic_part"] <- "exon"
x$group <- as.character(x$group)
x$group <- sapply(x$group, function(w) {
  a <- strsplit(w, "; ")[[1]]
  gid <- gsub("gene_id ", "", a[grep("gene_id", a)])
  eid <- gsub("exonic_part_number ", "", a[grep("exonic_part_number", a)])
  gexid <- paste0("exon_id ", gid, ":", eid)
  gexid <- gsub("\":\"", ":", gexid)
  paste(c(a, gexid), collapse = "; ")
})

export(x, con = output_gff_file, format = "GFF")
