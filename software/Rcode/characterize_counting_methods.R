## ----- characterize_counting_methods
## <<characterize_counting_methods.R>>

## Compare the various counting methods in terms of the number of counting 
## bins per gene and the count distributions (using a single sample).

## Read input arguments
args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

print(path_to_count_tables)
print(count_methods)
print(output_file)

L <- lapply(count_methods, function(cm) {
  fl <- list.files(paste0(path_to_count_tables, "/", cm), 
                   pattern = "1.txt$", full.names = TRUE)
  print(fl)
  tmp <- read.delim(fl, header = FALSE, as.is = TRUE)
  colnames(tmp) <- c("bin", "count")
  tmp$gene <- sapply(strsplit(tmp$bin, ":"), .subset, 1)
  tmp
})

methods <- c("dexseq", "dexseq_nomerge", 
             "featurecounts_flat",  "featurecounts_noflat", 
             "casper", "miso_assignable", "splicinggraph",
             "tophat_junc", "kallisto")
method_names <- c("DEXSeq-default", "DEXSeq-noaggreg",
                  "featureCounts-flat", "featureCounts-exon",
                  "casper", "MISO", "SplicingGraph", 
                  "TopHat-junctions", "kallisto")
names(method_names) <- methods
method_cols <- c(rgb(240, 228, 66, maxColorValue = 255),
                 rgb(0, 158, 115, maxColorValue = 255),
                 rgb(230, 159, 0, maxColorValue = 255),
                 rgb(86, 180, 233, maxColorValue = 255),
                 rgb(0, 0, 0, maxColorValue = 255),
                 rgb(0, 114, 178, maxColorValue = 255),
                 rgb(213, 94, 0, maxColorValue = 255),
                 rgb(204, 121, 167, maxColorValue = 255),
                 rgb(0, 255, 0, maxColorValue = 255))
names(method_cols) <- method_names
colors <- method_cols[method_names[match(count_methods, names(method_names))]]
names(L) <- method_names[match(count_methods, methods)]
##colors <- method_cols[count_methods]
##colors <- c("black", "red", "green", "blue", "cyan", "orange", "grey", 
##            "pink", "magenta")

library(dplyr)
## Calculate density of number of counting bins per gene
dens_log10nbins_pergene <- lapply(L, function(w) {
  tmp <- w %>% group_by(gene) %>% summarise(nbin = length(bin))
  density(log10(tmp$nbin), from = 0)
})
rx <- range(sapply(dens_log10nbins_pergene, function(w) range(w$x)))
ry <- range(sapply(dens_log10nbins_pergene, function(w) range(w$y)))
pdf(output_file)
plot(dens_log10nbins_pergene[[1]], xlim = rx, ylim = ry, col = colors[1], 
     xlab = "log10(number of bins per gene)", sub = "", lwd = 3, main = "", 
     cex.lab = 1.5)
for (i in 2:length(dens_log10nbins_pergene)) {
  lines(dens_log10nbins_pergene[[i]], col = colors[i], lwd = 3)
}
legend("topright", names(dens_log10nbins_pergene), lwd = 3, 
       col = colors[1:length(dens_log10nbins_pergene)])

## Calculate and plot density of count per bin + 1
dens_log2countperbin <- lapply(L, function(w) {
  density(log2(w$count + 1), from = 0)
})
rx <- range(sapply(dens_log2countperbin, function(w) range(w$x)))
ry <- range(sapply(dens_log2countperbin, function(w) range(w$y)))
plot(dens_log2countperbin[[1]], xlim = rx, ylim = ry, col = colors[1], 
     xlab = "log2(count per bin + 1)", sub = "", lwd = 3, main = "", 
     cex.lab = 1.5)
for (i in 2:length(dens_log2countperbin)) {
  lines(dens_log2countperbin[[i]], col = colors[i], lwd = 3)
}
legend("topright", names(dens_log2countperbin), lwd = 3, 
       col = colors[1:length(dens_log2countperbin)])

## Calculate and plot density of count per gene + 1
dens_log2countpergene <- lapply(L, function(w) {
  tmp <- w %>% group_by(gene) %>% summarise(countgene = sum(count))
  density(log2(tmp$countgene + 1), from = 0)
})
rx <- range(sapply(dens_log2countpergene, function(w) range(w$x)))
ry <- range(sapply(dens_log2countpergene, function(w) range(w$y)))
plot(dens_log2countpergene[[1]], xlim = rx, ylim = ry, col = colors[1], 
     xlab = "log2(count per gene + 1)", sub = "", lwd = 3, main = "", 
     cex.lab = 1.5)
for (i in 2:length(dens_log2countpergene)) {
  lines(dens_log2countpergene[[i]], col = colors[i], lwd = 3)
}
legend("topright", names(dens_log2countpergene), lwd = 3, 
       col = colors[1:length(dens_log2countpergene)])
dev.off()
