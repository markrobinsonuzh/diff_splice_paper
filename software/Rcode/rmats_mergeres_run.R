## ----- rmats_mergeres_run
## <<rmats_mergeres_run.R>>

## Put together rMATS results
args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

print(mats_dir)
print(output_file)
print(method_name)

library(DEXSeq)

## Function 'perGeneQvalueExact' from DEXSeq (1.14.0)
perGeneQValueExact = function(pGene, theta, geneSplit) {
  stopifnot(length(pGene) == length(geneSplit))
  
  ## Compute the numerator \sum_{i=1}^M 1-(1-theta)^{n_i}
  ## Below we first identify the summands which are the same
  ## (because they have the same n_i), then do the sum via the
  ## mapply
  numExons     = listLen(geneSplit)
  tab          = tabulate(numExons)
  notZero      = (tab>0)
  numerator    = mapply(function(m, n) m * (1 - (1-theta)^n),
                        m = tab[notZero],
                        n = which(notZero))
  numerator    = rowSums(numerator)
  
  ## Compute the denominator: for each value of theta, the number
  ## of genes with pGene <= theta[i].
  ## Note that in cut(..., right=TRUE), the intervals are
  ## right-closed (left open) intervals.
  bins   = cut(pGene, breaks = c(-Inf, as.vector(theta)), right = TRUE, include.lowest = TRUE)
  counts = tabulate(bins, nbins = nlevels(bins))
  denom  = cumsum(counts)
  stopifnot(denom[length(denom)]==length(pGene))
  
  return(numerator/denom)
}

## On target and junction reads
a <- read.table(paste0(mats_dir, "/A3SS.MATS.ReadsOnTargetAndJunctionCounts.txt"), 
                header = TRUE)[, c("GeneID", "PValue", "FDR")]
b <- read.table(paste0(mats_dir, "/A5SS.MATS.ReadsOnTargetAndJunctionCounts.txt"), 
                header = TRUE)[, c("GeneID", "PValue", "FDR")]
d <- read.table(paste0(mats_dir, "/MXE.MATS.ReadsOnTargetAndJunctionCounts.txt"), 
                header = TRUE)[, c("GeneID", "PValue", "FDR")]
e <- read.table(paste0(mats_dir, "/RI.MATS.ReadsOnTargetAndJunctionCounts.txt"), 
                header = TRUE)[, c("GeneID", "PValue", "FDR")]
f <- read.table(paste0(mats_dir, "/SE.MATS.ReadsOnTargetAndJunctionCounts.txt"), 
                header = TRUE)[, c("GeneID", "PValue", "FDR")]

all <- rbind(a, b, d, e, f)
colnames(all)[1] <- "gene_id"
all <- all[order(all$gene_id), ]

## Summarise on gene level
geneSplit = split(seq(along = all$gene_id), all$gene_id)
pGene = sapply(geneSplit, function(i) min(all$PValue[i]))
stopifnot(all(is.finite(pGene)))
theta = unique(sort(pGene))
q = perGeneQValueExact(pGene, theta, geneSplit)
qres_P = rep(NA_real_, length(pGene))
qres_P = q[match(pGene, theta)]
qres_P = pmin(1, qres_P)
names(qres_P) = names(geneSplit)
stopifnot(!any(is.na(qres_P)))

## Write to file
tmp <- data.frame(gene = names(qres_P), "adjP" = qres_P)
colnames(tmp)[colnames(tmp) == "adjP"] <- paste0(method_name, "_target_junction:adjP")
write.table(tmp, file = gsub("\\.txt", "_ontarget_junction.txt", output_file), 
            col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

## ----------------------------------------------------------------- ##
## Only junction reads
a <- read.table(paste0(mats_dir, "/A3SS.MATS.JunctionCountOnly.txt"), 
                header = TRUE)[, c("GeneID", "PValue", "FDR")]
b <- read.table(paste0(mats_dir, "/A5SS.MATS.JunctionCountOnly.txt"), 
                header = TRUE)[, c("GeneID", "PValue", "FDR")]
d <- read.table(paste0(mats_dir, "/MXE.MATS.JunctionCountOnly.txt"), 
                header = TRUE)[, c("GeneID", "PValue", "FDR")]
e <- read.table(paste0(mats_dir, "/RI.MATS.JunctionCountOnly.txt"), 
                header = TRUE)[, c("GeneID", "PValue", "FDR")]
f <- read.table(paste0(mats_dir, "/SE.MATS.JunctionCountOnly.txt"), 
                header = TRUE)[, c("GeneID", "PValue", "FDR")]

all <- rbind(a, b, d, e, f)
colnames(all)[1] <- "gene_id"
all <- all[order(all$gene_id), ]

## Summarise on gene level
geneSplit = split(seq(along = all$gene_id), all$gene_id)
pGene = sapply(geneSplit, function(i) min(all$PValue[i]))
stopifnot(all(is.finite(pGene)))
theta = unique(sort(pGene))
q = perGeneQValueExact(pGene, theta, geneSplit)
qres_P = rep(NA_real_, length(pGene))
qres_P = q[match(pGene, theta)]
qres_P = pmin(1, qres_P)
names(qres_P) = names(geneSplit)
stopifnot(!any(is.na(qres_P)))

## Write to file
tmp <- data.frame(gene = names(qres_P), "adjP" = qres_P)
colnames(tmp)[colnames(tmp) == "adjP"] <- paste0(method_name, "_junction:adjP")
write.table(tmp, file = gsub("\\.txt", "_junction.txt", output_file), 
            col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")



sessionInfo()
