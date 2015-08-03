## ----- generate_flattened_gtf
## <<generate_flattened_gtf.R>>

## Generate manually flattened gtf files to use with featureCounts
args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

print(input_gtf)
print(output_gtf)
print(ignore_strand)

library(GenomicRanges)
library(rtracklayer)

## Import gtf file
gtf0 <- import(input_gtf)

## Keep only exons
idx <- mcols(gtf0)$type == "exon" 
gtf0 <- gtf0[idx]
  
# Split by gene
gtf.g <- split(gtf0, mcols(gtf0)$gene_id)

# Flatten per gene
gtf.g.flat <- disjoin(gtf.g)

gtf.g.flat <- unlist(gtf.g.flat)
mcols(gtf.g.flat)$gene_id <- names(gtf.g.flat)

# Flatten overall
gtf.o.flat <- disjoin(gtf.g.flat, ignore.strand = ignore_strand)

# Remove overlapping parts
co <- countOverlaps(gtf.o.flat, gtf.g.flat)
gtf.o.flat <- gtf.o.flat[co == 1]
  
# Find the gene names
fo <- findOverlaps(gtf.o.flat, gtf.g.flat)
mcols(gtf.o.flat)$gene_id <- mcols(gtf.g.flat)$gene_id[subjectHits(fo)]
mcols(gtf.o.flat)$transcript_id <- mcols(gtf.o.flat)$gene_id
mcols(gtf.o.flat)$type <- "exon"

gtf.o.flat.sort <- gtf.o.flat[order(mcols(gtf.o.flat)$gene_id, decreasing = FALSE)]

exon.id <- split(mcols(gtf.o.flat.sort)$gene_id, mcols(gtf.o.flat.sort)$gene_id)

exon.id.new <- unlist(lapply(exon.id, function(g){ 
  seq(1, length(g)) 
}))

mcols(gtf.o.flat.sort)$exon_number <- exon.id.new
mcols(gtf.o.flat.sort)$exon_id <- paste0(mcols(gtf.o.flat.sort)$gene_id,
                                         ":", sprintf( "%03d", exon.id.new))
export(gtf.o.flat.sort, output_gtf, format = "gtf")

## Fix gff file for DEXSeq. Note! Still need to change the strand to . instead of *
gtf.o.flat.sort$type[gtf.o.flat.sort$type == "exon"] <- "exonic_part"
mcols(gtf.o.flat.sort)$exonic_part_number <- sprintf("%03d", gtf.o.flat.sort$exon_number)
tmp <- gtf.o.flat.sort
mcols(tmp)$group <- paste0("transcripts \"", tmp$transcript_id, 
                           "\"; exonic_part_number \"", tmp$exonic_part_number, 
                           "\"; gene_id \"", tmp$gene_id, "\"")
  
export(tmp, gsub("\\.gtf$", ".gff", output_gtf), format = "gff")


