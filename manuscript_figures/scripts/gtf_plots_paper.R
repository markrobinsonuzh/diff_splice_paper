## ----- gtf_plots_paper
## <<gtf_plots_paper.R>>

path_to_human_gtf <- "/home/Shared/data/alt_splicing_simulations/Simulation5_Charlotte/hsapiens/reference_files/Homo_sapiens.GRCh37.71.primary_assembly.protein_coding.gtf"
path_to_drosophila_gtf <- "/home/Shared/data/alt_splicing_simulations/Simulation5_Charlotte/drosophila/reference_files/Drosophila_melanogaster.BDGP5.70.protein_coding.gtf"

path_to_human_gff <- paste0("/home/Shared/data/alt_splicing_simulations/", 
                            "Simulation5_Charlotte/hsapiens/reference_files/", 
                            "Homo_sapiens.GRCh37.71.primary_assembly.protein_coding.flattened.gff")
path_to_drosophila_gff <- paste0("/home/Shared/data/alt_splicing_simulations/", 
                                 "Simulation5_Charlotte/drosophila/", 
                                 "reference_files/Drosophila_melanogaster.BDGP5.70.protein_coding.flattened.gff")
rsem_isoforms_drosophila <- 
  read.delim(paste0("/home/Shared/data/alt_splicing_simulations/", 
                    "Simulation5_Charlotte/drosophila/reference_files/", 
                    "rsem_model/SRR1501444.isoforms.results"), 
             header = TRUE, as.is = TRUE)
rsem_isoforms_human <- 
  read.delim(paste0("/home/Shared/data/alt_splicing_simulations/", 
                    "Simulation5_Charlotte/hsapiens/reference_files/", 
                    "rsem_model/SRR493366.isoforms.results"), 
             header = TRUE, as.is = TRUE)

library(rtracklayer)
library(dplyr)
library(ggplot2)
library(reshape2)
library(GenomicFeatures)

## Read gtf/gff files
gtf_human <- import(path_to_human_gtf)
gtf_drosophila <- import(path_to_drosophila_gtf)
gff_human <- import(path_to_human_gff)
gff_drosophila <- import(path_to_drosophila_gff)

## Keep only exons
gtf_human <- subset(gtf_human, type == "exon")
gtf_drosophila <- subset(gtf_drosophila, type == "exon")
gff_human <- subset(gff_human, type == "exonic_part")
gff_drosophila <- subset(gff_drosophila, type == "exonic_part")

## Extract gene ID from gff annotation
mcols(gff_human)$gene_id <- sapply(as.character(mcols(gff_human)$group), function(w) {
  a <- strsplit(w, ";")[[1]]
  a <- a[grep("gene_id", a)]
  gsub(" ", "", gsub("gene_id", "", a))
})
mcols(gff_drosophila)$gene_id <- sapply(as.character(mcols(gff_drosophila)$group), function(w) {
  a <- strsplit(w, ";")[[1]]
  a <- a[grep("gene_id", a)]
  gsub(" ", "", gsub("gene_id", "", a))
})

## Number of transcripts per gene
ntxbygene_human <- as.data.frame(mcols(gtf_human)) %>% 
  group_by(gene_id) %>% summarise(ntx = length(unique(transcript_id)),
                                  organism = "Human",
                                  measure = "Number of transcripts per gene")
ntxbygene_drosophila <- as.data.frame(mcols(gtf_drosophila)) %>% 
  group_by(gene_id) %>% summarise(ntx = length(unique(transcript_id)),
                                  organism = "Drosophila",
                                  measure = "Number of transcripts per gene")
ntxbygene <- rbind(ntxbygene_human, ntxbygene_drosophila)
colnames(ntxbygene)[match(c("gene_id", "ntx"), colnames(ntxbygene))] <- 
  c("group_id", "value")

## Number of exons by transcript
nexbytx_human <- as.data.frame(mcols(gtf_human)) %>% 
  group_by(transcript_id) %>% summarise(nex = length(unique(exon_number)),
                                        organism = "Human",
                                        measure = "Number of exons per transcript")
nexbytx_drosophila <- as.data.frame(mcols(gtf_drosophila)) %>% 
  group_by(transcript_id) %>% summarise(nex = length(unique(exon_number)),
                                        organism = "Drosophila",
                                        measure = "Number of exons per transcript")
nexbytx <- rbind(nexbytx_human, nexbytx_drosophila)
colnames(nexbytx)[match(c("transcript_id", "nex"), colnames(nexbytx))] <- 
  c("group_id", "value")

## Transcript length
tmp <- split(gtf_human, gtf_human$transcript_id)
tx_length_human <- sum(width(tmp))
tmp <- split(gtf_drosophila, gtf_drosophila$transcript_id)
tx_length_drosophila <- sum(width(tmp))
tx_length <- rbind(data.frame(transcript_id = names(tx_length_human),
                              tx_length = tx_length_human, 
                              organism = "Human",
                              measure = "Transcript length (bp)", 
                              stringsAsFactors = FALSE),
                   data.frame(transcript_id = names(tx_length_drosophila),
                              tx_length = tx_length_drosophila, 
                              organism = "Drosophila",
                              measure = "Transcript length (bp)",
                              stringsAsFactors = FALSE))
colnames(tx_length)[match(c("transcript_id", "tx_length"), colnames(tx_length))] <- 
  c("group_id", "value")

## Exon length
tmp <- gtf_human[-which(duplicated(gtf_human$exon_id))]
exon_length_human <- width(tmp)
tmp <- gtf_drosophila[-which(duplicated(gtf_drosophila$exon_id))]
exon_length_drosophila <- width(tmp)
exon_length <- rbind(data.frame(exon_id = 1:length(exon_length_human),
                                exon_length = exon_length_human, 
                                organism = "Human",
                                measure = "Exon length (bp)", 
                                stringsAsFactors = FALSE),
                     data.frame(exon_id = 1:length(exon_length_drosophila),
                                exon_length = exon_length_drosophila, 
                                organism = "Drosophila",
                                measure = "Exon length (bp)",
                                stringsAsFactors = FALSE))
colnames(exon_length)[match(c("exon_id", "exon_length"), colnames(exon_length))] <- 
  c("group_id", "value")

## Bins/genes in complexes
n_bins_human <- length(gff_human)
n_bins_in_complex_human <- length(gff_human[grep("\\+", mcols(gff_human)$gene_id), ])
n_bins_drosophila <- length(gff_drosophila)
n_bins_in_complex_drosophila <- length(gff_drosophila[grep("\\+", mcols(gff_drosophila)$gene_id), ])

n_genes_human <- length(unique(mcols(gff_human)$gene_id))
n_complexes_human <- length(unique(mcols(gff_human)$gene_id[grep("\\+", mcols(gff_human)$gene_id)]))
n_genes_drosophila <- length(unique(mcols(gff_drosophila)$gene_id))
n_complexes_drosophila <- length(unique(mcols(gff_drosophila)$gene_id[grep("\\+", mcols(gff_drosophila)$gene_id)]))

## Number of transcripts with 0/non0 reads
n_tx_0_drosophila <- length(which(rsem_isoforms_drosophila$TPM == 0))
n_tx_non0_drosophila <- length(which(rsem_isoforms_drosophila$TPM != 0))
n_tx_0_human <- length(which(rsem_isoforms_human$TPM == 0))
n_tx_non0_human <- length(which(rsem_isoforms_human$TPM != 0))

## Summary tables
summary_drosophila <- c(nbr_genes = length(ntxbygene_drosophila$gene_id),
                        nbr_transcripts = length(nexbytx_drosophila$transcript_id),
                        nbr_exons = length(gtf_drosophila),
                        nbr_unique_exons = length(unique(gtf_drosophila$exon_id)),
                        min_tx_by_gene = min(ntxbygene_drosophila$ntx),
                        median_tx_by_gene = median(ntxbygene_drosophila$ntx),
                        mean_tx_by_gene = mean(ntxbygene_drosophila$ntx),
                        max_tx_by_gene = max(ntxbygene_drosophila$ntx),
                        min_ex_by_tx = min(nexbytx_drosophila$nex),
                        median_ex_by_tx = median(nexbytx_drosophila$nex),
                        mean_ex_by_tx = mean(nexbytx_drosophila$nex),
                        max_ex_by_tx = max(nexbytx_drosophila$nex),
                        nbr_bins = n_bins_drosophila,
                        nbr_bins_in_complex = n_bins_in_complex_drosophila,
                        n_genes_complexes = n_genes_drosophila,
                        n_complexes = n_complexes_drosophila,
                        n_tx_0 = n_tx_0_drosophila,
                        n_tx_non0 = n_tx_non0_drosophila,
                        frac_tx_0 = n_tx_0_drosophila/(n_tx_0_drosophila + 
                                                         n_tx_non0_drosophila))
summary_human <- c(nbr_genes = length(ntxbygene_human$gene_id),
                   nbr_transcripts = length(nexbytx_human$transcript_id),
                   nbr_exons = length(gtf_human),
                   nbr_unique_exons = length(unique(gtf_human$exon_id)),
                   min_tx_by_gene = min(ntxbygene_human$ntx),
                   median_tx_by_gene = median(ntxbygene_human$ntx),
                   mean_tx_by_gene = mean(ntxbygene_human$ntx),
                   max_tx_by_gene = max(ntxbygene_human$ntx),
                   min_ex_by_tx = min(nexbytx_human$nex),
                   median_ex_by_tx = median(nexbytx_human$nex),
                   mean_ex_by_tx = mean(nexbytx_human$nex),
                   max_ex_by_tx = max(nexbytx_human$nex),
                   nbr_bins = n_bins_human,
                   nbr_bins_in_complex = n_bins_in_complex_human,
                   n_genes_complexes = n_genes_human,
                   n_complexes = n_complexes_human,
                   n_tx_0 = n_tx_0_human,
                   n_tx_non0 = n_tx_non0_human,
                   frac_tx_0 = n_tx_0_human/(n_tx_0_human + 
                                               n_tx_non0_human))

write.table(data.frame(summary_drosophila), 
            file = paste0("/home/Shared/data/alt_splicing_simulations/",
                          "Simulation5_Charlotte/manuscript_figures/", 
                          "figures/summary_drosophila.txt"), 
            row.names = TRUE, col.names = FALSE, quote = FALSE, sep = "\t")
write.table(data.frame(summary_human), 
            file = paste0("/home/Shared/data/alt_splicing_simulations/",
                          "Simulation5_Charlotte/manuscript_figures/", 
                          "figures/summary_human.txt"), 
            row.names = TRUE, col.names = FALSE, quote = FALSE, sep = "\t")

## Plots
pdf(paste0("/home/Shared/data/alt_splicing_simulations/Simulation5_Charlotte/", 
           "manuscript_figures/figures/gtf_characterisation.pdf"))
df <- rbind(ntxbygene, nexbytx, tx_length, exon_length)
print(ggplot(df, aes(x = organism, y = value, fill = organism)) + geom_violin() + 
        scale_y_log10() + facet_wrap(~measure, scales = "free") + 
        ylab("") + xlab("") + #coord_flip() + 
        scale_fill_manual(values = c("blue", "red"), guide = FALSE) + 
        guides(fill = guide_legend(nrow = 1, size = 15)) + 
        theme(legend.position = "bottom", 
              panel.background = element_rect(fill = NA, colour = "black"),
              panel.border = element_rect(fill = NA, colour = "black", size = 1),
              panel.grid.minor.x = element_blank(),
              panel.grid.minor.y = element_blank(),
              axis.text.x = element_text(size = 15),
              axis.text.y = element_text(size = 15),
              axis.title.x = element_text(size = 15),
              axis.title.y = element_text(size = 15),
              strip.text = element_text(size = 12),
              strip.background = element_rect(fill = NA, colour = "black")))

print(ggplot(df, aes(x = value, col = organism)) + 
        geom_path(stat = "density", size = 2, alpha = 0.6) + 
        scale_x_log10() + facet_wrap(~measure, scales = "free") + 
        ylab("Density") + xlab("") + #coord_flip() + 
        scale_color_manual(values = c("blue", "red"), guide = FALSE) + 
        guides(col = guide_legend(nrow = 1, size = 15)) + 
        theme(legend.position = "bottom", 
              panel.background = element_rect(fill = NA, colour = "black"),
              panel.border = element_rect(fill = NA, colour = "black", size = 1),
              panel.grid.minor.x = element_blank(),
              panel.grid.minor.y = element_blank(),
              axis.text.x = element_text(size = 15),
              axis.text.y = element_text(size = 15),
              axis.title.x = element_text(size = 15),
              axis.title.y = element_text(size = 15),
              strip.text = element_text(size = 12),
              strip.background = element_rect(fill = NA, colour = "black")))
dev.off()

