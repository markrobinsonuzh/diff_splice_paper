## ----- generate_renamed_gtf_for_featurecounts
## <<generate_renamed_gtf_for_featurecounts.R>>

## Read input arguments
args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

print(input_gtf)
print(output_gtf)

library(GenomicRanges)
library(rtracklayer)
library(dplyr)

gtf0 <- import(input_gtf)

gtf0 <- as.data.frame(gtf0[gtf0$type == "exon"])

gtf1 <- gtf0 %>% group_by(gene_id) %>% 
  mutate(exon_id = paste0(gene_id, ":", sprintf("%03d", 1:length(gene_id)))) %>%
  ungroup()
gtf1 <- as.data.frame(gtf1)
gtf2 <- makeGRangesFromDataFrame(gtf1, keep.extra.columns = TRUE)

export(gtf2, con = output_gtf, format = "gtf")