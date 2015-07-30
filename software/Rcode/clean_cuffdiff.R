## ----- clean_cuffdiff
## <<clean_cuffdiff.R>>

## Clean the result files returned by cuffdiff

## Read input arguments
args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

print(path_to_cuffdiff_result)
print(path_to_gtf_file)
print(output_file)
print(method_name)

library(rtracklayer)

gtf <- import(path_to_gtf_file)
gtf <- as.data.frame(gtf)

res <- read.delim(path_to_cuffdiff_result, header = TRUE, as.is = TRUE)

gene_id <- sapply(res$gene, function(g) {
  g <- strsplit(g, split = ",")[[1]]
  g2 <- unique(gtf$gene_id[gtf$gene_name %in% g])
  paste(g2, collapse = "+")
})

res$gene_ensembl <- gene_id
df <- res[, c("gene_ensembl", "q_value")]

## For the few genes where there are many representations (which happens if there are isoforms that don't overlap any other isoform), keep the smallest FDR.
library(dplyr)
df2 <- as.data.frame(df %>% group_by(gene_ensembl) %>% summarise(q_value = min(q_value)) %>% ungroup())

colnames(df2) <- c("gene", paste0(method_name, ":adjP"))

write.table(df2, file = output_file, sep = "\t", row.names = FALSE,
            col.names = TRUE, quote = FALSE)