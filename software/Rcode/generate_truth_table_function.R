## ----- generate_truth_table_function
## <<generate_truth_table_function.R>>

## Generate truth table for genes based on simulation summary
generate_truth_table <- function(final_summary, out.file, 
                                 astalavista.file = NULL,
                                 gtf.file = NULL,
                                 flattened.gtf.file = NULL,
                                 missing.annot.file = NULL) {
  library(GenomicRanges)
  library(GenomicFeatures)
  library(Hmisc)
  library(dplyr)
  library(rtracklayer)
  
  if (!is.null(missing.annot.file)) {
    mfile <- read.delim(missing.annot.file, header = TRUE)
    final_summary$missing_tr <- 0
    final_summary$missing_tr[match(mfile$transcript, 
                                   final_summary$transcript_id)] <- 1
  
    gene_table <- final_summary %>% group_by(gene_id) %>%
      summarise(ds_status = unique(gene_ds_status),
                de_status = unique(gene_de_status), 
                ds_de_status = paste0(ds_status, ":", de_status), 
                TPM = sum(TPM),
                nbr_isoforms = unique(nbr_isoforms),
                nbr_isoforms_above5 = length(which(IsoPct >= 5)),
                nbr_isoforms_above10 = length(which(IsoPct >= 10)),
                nbr_isoforms_above15 = length(which(IsoPct >= 15)),
                nbr_isoforms_above25 = length(which(IsoPct >= 25)),
                diff_IsoPct = unique(diff_IsoPct),
                nbr_missing_ds_tr = sum(missing_tr[transcript_ds_status == 1]),
                nbr_missing_nonds_tr = sum(missing_tr[transcript_ds_status == 0]),
                nbr_missing_tr = sum(missing_tr),
                nbr_remaining_tr = nbr_isoforms - nbr_missing_tr)
    gene_table$nbr_missing_tr_cat <- cut2(gene_table$nbr_missing_tr, 
                                          cuts = c(0, 1, 2, 3, 4))
    gene_table$nbr_remaining_tr_cat <- cut2(gene_table$nbr_remaining_tr,
                                            cuts = c(0, 1, 2, 3, 4))
  } else {
    gene_table <- final_summary %>% group_by(gene_id) %>%
      summarise(ds_status = unique(gene_ds_status),
                de_status = unique(gene_de_status), 
                ds_de_status = paste0(ds_status, ":", de_status), 
                TPM = sum(TPM),
                nbr_isoforms = unique(nbr_isoforms),
                nbr_isoforms_above5 = length(which(IsoPct >= 5)),
                nbr_isoforms_above10 = length(which(IsoPct >= 10)),
                nbr_isoforms_above15 = length(which(IsoPct >= 15)),
                nbr_isoforms_above25 = length(which(IsoPct >= 25)),
                diff_IsoPct = unique(diff_IsoPct))
  }
  gene_table$TPM_2 = cut2(gene_table$TPM, g = 2)
  gene_table$TPM_5 = cut2(gene_table$TPM, g = 5)
  gene_table$TPM_10 = cut2(gene_table$TPM, g = 10)
  gene_table$nbr_isoforms_2 = cut2(gene_table$nbr_isoforms, g = 2)
  gene_table$nbr_isoforms_5 = cut2(gene_table$nbr_isoforms, g = 5)
  gene_table$nbr_isoforms_10 = cut2(gene_table$nbr_isoforms, g = 10)
  gene_table$diff_IsoPct_2 = cut2(gene_table$diff_IsoPct, g = 2)
  gene_table$diff_IsoPct_5 = cut2(gene_table$diff_IsoPct, g = 5)
  gene_table$diff_IsoPct_10 = cut2(gene_table$diff_IsoPct, g = 10)
  gene_table$diff_IsoPct_3eq <- cut2(gene_table$diff_IsoPct, 
                                     cuts = c(0, 1/3, 2/3, 1))
  
  if (!is.null(gtf.file)) {
    ## Calculate the number of basepairs by which the two most abundant
    ## isoforms differ. Only for the differentially spliced genes. To
    ## do it for all genes, just remove the subset when creating tmp.
    tmp <- subset(final_summary, transcript_ds_status == 1) %>% 
      group_by(gene_id) %>%
      summarise(topisoform = transcript_id[which.max(IsoPct)],
                secisoform = ifelse(length(transcript_id) > 1, 
                                    transcript_id[order(IsoPct, 
                                                        decreasing = TRUE)[2]], 
                                    "NA"))
    txdb <- makeTxDbFromGFF(gtf.file, format = "gtf")
    ebt <- exonsBy(txdb, "tx", use.names = TRUE)
    ebt2 <- ebt[setdiff(c(tmp$topisoform, tmp$secisoform), "NA")]
    tmp$diffbp <- sapply(1:nrow(tmp), function(i) {
      if (any(tmp[i, ] == "NA")) 0
      else {
        a <- ebt2[[unlist(tmp[i, "topisoform"])]]
        b <- ebt2[[unlist(tmp[i, "secisoform"])]]
        sum(width(GenomicRanges::union(a, b))) - 
          sum(width(GenomicRanges::intersect(a, b)))
      }
    })
    tmp$diffbp_2 <- cut2(tmp$diffbp, g = 2)
    tmp$diffbp_5 <- cut2(tmp$diffbp, g = 5)
    tmp$diffbp_10 <- cut2(tmp$diffbp, g = 10)
    gene_table <- merge(gene_table, tmp[, c("gene_id", "diffbp", 
                                            "diffbp_2", "diffbp_5",
                                            "diffbp_10")], 
                        by = "gene_id", all = TRUE)
  }
  
  if (!is.null(astalavista.file)) {
    ## Add the type of splicing event for the differentially spliced genes
    astalavista <- read.delim(astalavista.file, header = FALSE, as.is = TRUE)
    av <- t(sapply(astalavista[, 9], function(w) {
      tmp <- unlist(strsplit(w, ";"))
      a <- gsub(" ", "", gsub("structure ", "", tmp[grep("structure", tmp)]))
      b <- gsub("transcript_id ", "", tmp[grep("transcript_id", tmp)])
      b1 <- paste(sort(unlist(strsplit(b, ","))), collapse = ";")
      c(a, b1)
    }))
    rownames(av) <- NULL
    colnames(av) <- c("as_type", "transcripts")
    av <- as.data.frame(av, stringsAsFactors = FALSE)
    
    av1 <- av[-grep("/", av[, 2]), ]
    av2 <- av[grep("/", av[, 2]), ]
    
    ## The "transcripts" column may combine multiple transcript combinations
    myfun = function(mystr) {
      mm = lapply(strsplit(mystr, ";"), strsplit, split = "/")
      mm2 = expand.grid(mm[[1]][[1]], mm[[1]][[2]], stringsAsFactors = FALSE)
      sapply(1:nrow(mm2), function(i) paste(sort(mm2[i, ]), collapse = ";"))
    }
    tmp_trs <- lapply(av2[, 2], myfun)
    names(tmp_trs) <- av2[, 1]
    tmp_trs <- stack(tmp_trs)
    colnames(tmp_trs) <- c("transcripts", "as_type1")
    tmp_trs <- tmp_trs %>% group_by(transcripts) %>% 
      summarise(as_type = paste(as_type1, collapse = ":"))
    tmp_trs <- tmp_trs[, c("as_type", "transcripts")]
    av <- rbind(av1, tmp_trs)
    
    ds_transcripts <- subset(final_summary, transcript_ds_status == 1) %>% 
      group_by(gene_id) %>%
      summarise(transcripts = paste(sort(transcript_id), collapse = ";"))
    ds_transcripts$as_type <- av$as_type[match(ds_transcripts$transcripts, av$transcripts)]
    gene_table <- merge(gene_table, ds_transcripts[, c("gene_id", "as_type")], 
                        by.x = "gene_id", by.y = "gene_id", all = TRUE)
    
    gene_table$as_type_common <- gene_table$as_type
    gene_table$as_type_common[!gene_table$as_type_common %in% 
                                c(NA, "0,1-2^", "1^,2^", "1-,2-", 
                                  "0,1^2-")] <- "other"
    gene_table$as_type_common[gene_table$as_type_common == "0,1-2^"] <- "skipped exon"
    gene_table$as_type_common[gene_table$as_type_common == "1^,2^"] <- "alternative donors"
    gene_table$as_type_common[gene_table$as_type_common == "1-,2-"] <- "alternative acceptors"
    gene_table$as_type_common[gene_table$as_type_common == "0,1^2-"] <- "retained intron"
    
  }
  
  if (!is.null(flattened.gtf.file)) {
    ## Extract the number of exon bins for each gene
    flatfile <- import(flattened.gtf.file)
    flatfile <- subset(flatfile, type == "exonic_part")
    genes <- sapply(strsplit(as.character(flatfile$group), ";"), .subset, 3)
    genes <- gsub(" gene_id \"", "", genes)
    genes <- gsub("\"", "", genes)
    tbl <- as.data.frame(table(genes), stringsAsFactors = FALSE)
    colnames(tbl)[colnames(tbl) == "Freq"] <- "nbrexonbins"
    tbl$nbrexonbins_2 <- cut2(tbl$nbrexonbins, g = 2)
    tbl$nbrexonbins_5 <- cut2(tbl$nbrexonbins, g = 5)
    tbl$nbrexonsbins_10 <- cut2(tbl$nbrexonbins, g = 10)
    gene_table <- merge(gene_table, tbl, by.x = "gene_id",
                        by.y = "genes", all = TRUE)
  }
  
  idx <- which(colnames(gene_table) == "gene_id")
  colnames(gene_table)[idx] <- "gene"
  
  write.table(gene_table, file = out.file, 
              row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
}