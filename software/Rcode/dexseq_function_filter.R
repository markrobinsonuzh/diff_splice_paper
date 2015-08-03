## ----- dexseq_function
## <<dexseq_function.R>>

## Run DEXSeq on a given set of count files
run_dexseq <- function(count_files, conditions, flattened_file = NULL,
                       out_basename, method_name, filter_bin_count = NULL,
                       filter_bin_variance = NULL, 
                       filter_perc = NULL, filter_perc_perkb = NULL, 
                       filter_perc_pertx = NULL, 
                       filter_perc_perkb_pertx = NULL) {
  
  if (length(sum(!is.null(c(filter_bin_count, filter_perc, filter_perc_perkb, 
                            filter_perc_pertx, filter_perc_perkb_pertx)))) > 1)
    stop("Only one filter criterion is allowed")
  
  library(DEXSeq)
  BPPARAM <- MulticoreParam(workers = 6)
  
  ## Generate sample table
  sample_table <- data.frame(sample = 1:length(conditions),
                             condition = conditions)
  
  ## Generate DEXSeqDataSet
  message("Generating DEXSeqDataSet...")
  dxd <- DEXSeqDataSetFromHTSeq(count_files, 
                                sampleData = sample_table,
                                design = ~ sample + exon + condition:exon,
                                flattenedfile = flattened_file)
  
  dxd <- estimateSizeFactors(dxd)
  
  ## Filter by normalized bin counts
  if (!is.null(filter_bin_count)) {
    norm_counts <- featureCounts(dxd, normalized = TRUE)
    keep <- which(rowSums(norm_counts) >= filter_bin_count)
    dxd <- dxd[keep, ]
  }
  
  ## Filter by variance of log2(normalized count + 1)
  if (!is.null(filter_bin_variance)) {
    norm_counts <- featureCounts(dxd, normalized = TRUE)
    vars <- apply(log2(norm_counts + 1), 1, var)
    keep <- which(vars >= filter_bin_variance)
    dxd <- dxd[keep, ]
  }
  
  ## Filter by relative abundance of bins (without accounting for bin size)
  if (!is.null(filter_perc)) {
    norm_count <- featureCounts(dxd, normalized = TRUE)
    norm_count <- as.data.frame(norm_count, stringsAsFactors = FALSE)
    norm_count$gene <- sapply(strsplit(rownames(norm_count), ":"),
                              .subset, 1)
    library(dplyr)
    count_perc <- as.data.frame(norm_count %>% group_by(gene) %>% 
                                  mutate_each(funs(perc = . / sum(.))))
    count_perc[is.na(count_perc)] <- 0
    rownames(count_perc) <- rownames(norm_count)
    count_perc$gene <- NULL
    max_count_perc <- apply(count_perc, 1, max)
    
    keep <- names(max_count_perc)[which(max_count_perc >= filter_perc)]
    dxd <- dxd[rownames(dxd) %in% keep, ]
  }
  
  ## Filter by relative abundance of bins after accounting for bin size
  if (!is.null(filter_perc_perkb)) {
    norm_count <- featureCounts(dxd, normalized = TRUE)
    library(rtracklayer)
    x <- import(flattened_file)
    x <- subset(x, type == "exonic_part")
    L <- width(x)
    names(L) <- sapply(as.character(x$group), function(i) {
      tmp <- strsplit(i, ";")[[1]]
      tmp2 <- tmp[grep("gene_id", tmp)]
      tmp3 <- tmp[grep("exonic_part_number", tmp)]
      tmp2 <- gsub("\"", "", gsub(" gene_id ", "", tmp2))
      tmp3 <- gsub("\"", "", gsub(" exonic_part_number ", "", tmp3))
      paste0(tmp2, ":E", tmp3)
    })
    L <- L[rownames(norm_count)]
    norm_count_perkb <- sweep(norm_count, 1, L, "/") * 1e3
    norm_count_perkb <- as.data.frame(norm_count_perkb, stringsAsFactors = FALSE)
    norm_count_perkb$gene <- sapply(strsplit(rownames(norm_count_perkb), ":"),
                                    .subset, 1)
    library(dplyr)
    count_perc <- as.data.frame(norm_count_perkb %>% group_by(gene) %>% 
                                  mutate_each(funs(perc = . / sum(.))))
    count_perc[is.na(count_perc)] <- 0
    rownames(count_perc) <- rownames(norm_count_perkb)
    count_perc$gene <- NULL
    max_count_perc <- apply(count_perc, 1, max)
    
    keep <- names(max_count_perc)[which(max_count_perc >= filter_perc_perkb)]
    dxd <- dxd[rownames(dxd) %in% keep, ]
  }
  
  ## Filter by relative abundance per isoform of bins without accounting for bin size
  if (!is.null(filter_perc_pertx)) {
    norm_count <- featureCounts(dxd, normalized = TRUE)
    library(rtracklayer)
    x <- import(flattened_file)
    x <- subset(x, type == "exonic_part")
    nbr_isoforms <- t(sapply(as.character(x$group), function(i) {
      tmp <- strsplit(i, ";")[[1]]
      tmp1 <- tmp[grep("transcripts", tmp)]
      tmp2 <- tmp[grep("gene_id", tmp)]
      tmp3 <- tmp[grep("exonic_part_number", tmp)]
      tmp1 <- gsub("\"", "", gsub("transcripts ", "", tmp1))
      tmp1 <- length(strsplit(tmp1, "\\+")[[1]])
      tmp2 <- gsub("\"", "", gsub(" gene_id ", "", tmp2))
      tmp3 <- gsub("\"", "", gsub(" exonic_part_number ", "", tmp3))
      c(paste0(tmp2, ":E", tmp3), tmp1)
    }))
    rn <- nbr_isoforms[, 1]
    nbr_isoforms <- as.numeric(nbr_isoforms[, 2])
    names(nbr_isoforms) <- rn
    nbr_isoforms <- nbr_isoforms[rownames(norm_count)]
    
    norm_count_pertx <- sweep(norm_count, 1, nbr_isoforms, "/")
    norm_count_pertx <- as.data.frame(norm_count_pertx, stringsAsFactors = FALSE)
    norm_count_pertx$gene <- sapply(strsplit(rownames(norm_count_pertx), ":"),
                                    .subset, 1)
    
    library(dplyr)
    count_perc <- as.data.frame(norm_count_pertx %>% group_by(gene) %>% 
                                  mutate_each(funs(perc = . / sum(.))))
    count_perc[is.na(count_perc)] <- 0
    rownames(count_perc) <- rownames(norm_count_pertx)
    count_perc$gene <- NULL
    max_count_perc <- apply(count_perc, 1, max)
    
    keep <- names(max_count_perc)[which(max_count_perc >= filter_perc_pertx)]
    dxd <- dxd[rownames(dxd) %in% keep, ]
  }
  
  ## Filter by relative abundance per isoform of bins after accounting for bin size
  if (!is.null(filter_perc_perkb_pertx)) {
    norm_count <- featureCounts(dxd, normalized = TRUE)
    library(rtracklayer)
    x <- import(flattened_file)
    x <- subset(x, type == "exonic_part")
    L <- width(x)
    names(L) <- sapply(as.character(x$group), function(i) {
      tmp <- strsplit(i, ";")[[1]]
      tmp2 <- tmp[grep("gene_id", tmp)]
      tmp3 <- tmp[grep("exonic_part_number", tmp)]
      tmp2 <- gsub("\"", "", gsub(" gene_id ", "", tmp2))
      tmp3 <- gsub("\"", "", gsub(" exonic_part_number ", "", tmp3))
      paste0(tmp2, ":E", tmp3)
    })
    L <- L[rownames(norm_count)]
    norm_count_perkb <- sweep(norm_count, 1, L, "/") * 1e3
    
    nbr_isoforms <- t(sapply(as.character(x$group), function(i) {
      tmp <- strsplit(i, ";")[[1]]
      tmp1 <- tmp[grep("transcripts", tmp)]
      tmp2 <- tmp[grep("gene_id", tmp)]
      tmp3 <- tmp[grep("exonic_part_number", tmp)]
      tmp1 <- gsub("\"", "", gsub("transcripts ", "", tmp1))
      tmp1 <- length(strsplit(tmp1, "\\+")[[1]])
      tmp2 <- gsub("\"", "", gsub(" gene_id ", "", tmp2))
      tmp3 <- gsub("\"", "", gsub(" exonic_part_number ", "", tmp3))
      c(paste0(tmp2, ":E", tmp3), tmp1)
    }))
    rn <- nbr_isoforms[, 1]
    nbr_isoforms <- as.numeric(nbr_isoforms[, 2])
    names(nbr_isoforms) <- rn
    nbr_isoforms <- nbr_isoforms[rownames(norm_count_perkb)]
    
    norm_count_perkb_pertx <- sweep(norm_count_perkb, 1, nbr_isoforms, "/")
    norm_count_perkb_pertx <- as.data.frame(norm_count_perkb_pertx, stringsAsFactors = FALSE)
    norm_count_perkb_pertx$gene <- 
      sapply(strsplit(rownames(norm_count_perkb_pertx), ":"), .subset, 1)
    
    library(dplyr)
    count_perc <- as.data.frame(norm_count_perkb_pertx %>% group_by(gene) %>% 
                                  mutate_each(funs(perc = . / sum(.))))
    count_perc[is.na(count_perc)] <- 0
    rownames(count_perc) <- rownames(norm_count_perkb_pertx)
    count_perc$gene <- NULL
    max_count_perc <- apply(count_perc, 1, max)
    
    keep <- names(max_count_perc)[which(max_count_perc >= filter_perc_perkb_pertx)]
    dxd <- dxd[rownames(dxd) %in% keep, ]
  }
  
  print(dim(featureCounts(dxd)))
  print(mean(colSums(featureCounts(dxd))))
  
  ## Estimate size factors and dispersions
  message("Estimating dispersions...")
  dxd <- estimateDispersions(dxd, BPPARAM = BPPARAM)
  
  ## Run the test and get results for exon bins and genes
  message("Running test...")
  dxd <- testForDEU(dxd, BPPARAM = BPPARAM)
  
  message("Summarizing results on gene level...")
  res <- DEXSeqResults(dxd)
  pgq <- perGeneQValue(res, p = "pvalue")
  
  ## Save results
  save(dxd, res, pgq, file = paste0(out_basename, ".Rdata"))
  tmp <- cbind(gene = names(pgq), "adjP" = pgq)
  colnames(tmp)[colnames(tmp) == "adjP"] <- paste0(method_name, ":adjP")
  write.table(tmp, 
              file = paste0(out_basename, ".txt"), 
              col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
  
  print(sessionInfo())
}
