## ----- clean_cuffdiff_ensembl
## <<clean_cuffdiff_ensembl.R>>

## Clean the result files returned by cuffdiff when used on a gtf file 
## where the gene names are the Ensembl IDs already

## Read input arguments
args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

print(path_to_cuffdiff_result)
print(output_file)
print(method_name)

res <- read.delim(path_to_cuffdiff_result, header = TRUE, as.is = TRUE)

df <- res[, c("gene", "q_value")]
## For the few genes where there are many representations (which happens if there are isoforms that don't overlap any other isoform), keep the smallest FDR.
library(dplyr)
df2 <- as.data.frame(df %>% group_by(gene) %>% summarise(q_value = min(q_value)) %>% ungroup())

colnames(df2) <- c("gene", paste0(method_name, ":adjP"))

write.table(df2, file = output_file, sep = "\t", row.names = FALSE,
            col.names = TRUE, quote = FALSE)