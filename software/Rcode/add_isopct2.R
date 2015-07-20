## ----- add_isopct2
## <<add_isopct2.R>>

## Add the column IsoPct2 (isoform percentages in condition 2) 
## to the simulation details file. 

## Read input arguments
args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

print(input_file)
print(output_file)

tmp <- read.delim(input_file, header = TRUE, as.is = TRUE)

mutate_IsoPct <- function(w, ref) {
  o <- order(ref, decreasing = TRUE)[1:2]
  w[rev(o)] <- w[o]
  w
}

library(dplyr)

tmp_ds <- subset(tmp, gene_ds_status == 1)
tmp_nonds <- subset(tmp, gene_ds_status == 0)

tmp_ds <- tmp_ds %>% group_by(gene_id) %>% 
  mutate(IsoPct2 = mutate_IsoPct(IsoPct, IsoPct))
tmp_nonds$IsoPct2 <- tmp_nonds$IsoPct

tmp <- rbind(tmp_ds, tmp_nonds)
tmp <- tmp[order(tmp$gene_id), ]

write.table(tmp, file = output_file, row.names = FALSE, 
            col.names = TRUE, quote = FALSE, sep = "\t")
