## Illustrate the different "flattening" options
## Use Gviz to illustrate the gene models for two overlapping genes, and 
## then show the bins created by each flattening approach

## Plot the coverage for each sample and a given gene
## Read input arguments
args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

print(path_to_gene_model_gtf_file)
print(paths_to_gff_files)
print(paths_to_gtf_files)
print(all_gene_to_display)
print(output_filename)
print(plot_height)
print(ext_before)
print(ext_after)

library(combinat)
library(Gviz)
options(ucscChromosomeNames = FALSE)
## Generate genome axis track (genome coordinates)
gtr <- GenomeAxisTrack()

## Load gene models
library(rtracklayer)
genemodels <- import(path_to_gene_model_gtf_file)
idx <- match(c("transcript_id", "gene_id", "exon_id"), colnames(mcols(genemodels)))
colnames(mcols(genemodels))[idx] <- c("transcript", "gene", "exon")
mcols(genemodels)$symbol <- mcols(genemodels)$transcript
genemodels <- subset(genemodels, type == "exon")

## Load individual gff/gtf files
fullbins <- list()
for (fl in names(paths_to_gff_files)) {
  tmp <- import(paths_to_gff_files[fl])
  tmp <- subset(tmp, type == "exonic_part")
  mcols(tmp)$exon <- sapply(as.character(mcols(tmp)$group), 
                            function(s) {
                              m <- strsplit(s, ";")[[1]]
                              m <- m[grep("exonic_part_number", m)]
                              gsub("\"", "", 
                                   gsub("^\"", "\"E", 
                                        gsub(" exonic_part_number ", "", m)))
                            })
  mcols(tmp)$transcript <- sapply(as.character(mcols(tmp)$group), 
                                  function(s) {
                                    m <- strsplit(s, ";")[[1]]
                                    m <- m[grep("gene_id", m)]
                                    gsub("\"", "", 
                                         gsub(" gene_id ", "", m))
                                  })
  mcols(tmp)$gene <- mcols(tmp)$transcript
  mcols(tmp)$symbol <- mcols(tmp)$transcript
  mcols(tmp)$feature <- mcols(tmp)$transcript
  fullbins[[fl]] <- tmp
}

for (fl in names(paths_to_gtf_files)) {
  tmp <- import(paths_to_gtf_files[fl])
  idx <- match(c("transcript_id", "gene_id", "exon_id"), colnames(mcols(tmp)))
  colnames(mcols(tmp))[idx] <- c("transcript", "gene", "exon")
  mcols(tmp)$symbol <- mcols(tmp)$transcript
  tmp <- subset(tmp, type == "exon")
  mcols(tmp)$feature <- mcols(tmp)$gene
  mcols(tmp)$group <- mcols(tmp)$gene
  fullbins[[fl]] <- tmp
}

pdf(output_filename, height = plot_height)
for (gene_to_display in all_gene_to_display) {
  ## Split the input gene complex name and create all possible combinations
  tmpgenes <- strsplit(gene_to_display, "\\+")[[1]]
  genes_to_display <- c()
  for (i in 1:length(tmpgenes)) {  ## number of elements to select
    cbn <- cbind(combn(length(tmpgenes), i))
    for (j in 1:ncol(cbn)) {
      pms <- lapply(permn(length(cbn[, j])), function(w) cbn[w, j])
      for (k in 1:length(pms)) {
        genes_to_display <- c(genes_to_display, paste(tmpgenes[pms[[k]]], collapse = "+"))
      }
    }
  }
  
  ## Subset to genes/complexes of interest
  gm <- subset(genemodels, gene %in% genes_to_display)
  show_chr <- unique(seqnames(gm))[1]
  min_coord <- min(start(gm)) - ext_before*(max(end(gm)) - min(start(gm)))
  max_coord <- max(end(gm)) + ext_after*(max(end(gm)) - min(start(gm)))
  grtr <- GeneRegionTrack(gm, showId = FALSE, name = "Gene model",
                          col = NULL, fill = "orange", cex.title = 1)
  group(grtr) <- mcols(gm)$transcript
  feature(grtr) <- mcols(gm)$gene
  
  bins <- lapply(fullbins, function(w) subset(w, gene %in% genes_to_display))
  
  ## Keep only gene combinations that appear in at least one file
  genes_to_display <- union(gm$gene, unlist(lapply(bins, function(w) w$gene)))
  
  ## Assign colors to the different genes/complexes
  geneColors <- c("lightblue", "green", "pink", "grey", "cyan",
                  "magenta", "yellow")[1:length(genes_to_display)]
  names(geneColors) <- genes_to_display
  
  backgroundColors <- c("#FFFEDB", "white", "#FFFEDB", "white")
  
  ## Assemble coloring arguments
  myoptions <- ""
  for (g in genes_to_display) {
    myoptions <- paste0(myoptions, "'", g, "'= '", geneColors[g], "', ")
  }
  myoptions <- gsub(", $", "", myoptions)
  
  ## Generate tracks for the different bins
  binTracks <- lapply(1:length(bins), function(i) {
    assign(paste0("bintr", i), 
           GeneRegionTrack(bins[[i]], showId = FALSE, 
                           name = names(bins)[i], 
                           col = "black", cex.title = 0.5,
                           background.panel = backgroundColors[i]))
  })
  
  ## Plot
  tracks <- c(grtr, gtr, binTracks)
  eval(parse(text = paste("plotTracks(tracks, 
                                   main = gene_to_display, cex.main = 1,
                                   from = min_coord, to = max_coord, 
                                   min.width = 0, min.distance = 0, collapse = FALSE,", 
                          myoptions, ")")))
  
}  
dev.off()
