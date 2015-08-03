## ----- characterize_rmats_fp
## <<characterize_rmats_fp.R>>

args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

print(mats_dir_human)
print(mats_dir_drosophila)
print(truth_file_human)
print(truth_file_drosophila)
print(output_pdf)

## ------------------------ ON TARGET AND JUNCTION ---------------------- ##
## HUMAN
a <- read.table(paste0(mats_dir_human, "/A3SS.MATS.ReadsOnTargetAndJunctionCounts.txt"), 
                header = TRUE)[, c("GeneID", "PValue", "FDR", "IncLevel1",
                                   "IncLevel2", "IncLevelDifference")]
b <- read.table(paste0(mats_dir_human, "/A5SS.MATS.ReadsOnTargetAndJunctionCounts.txt"), 
                header = TRUE)[, c("GeneID", "PValue", "FDR", "IncLevel1",
                                   "IncLevel2", "IncLevelDifference")]
d <- read.table(paste0(mats_dir_human, "/MXE.MATS.ReadsOnTargetAndJunctionCounts.txt"), 
                header = TRUE)[, c("GeneID", "PValue", "FDR", "IncLevel1",
                                   "IncLevel2", "IncLevelDifference")]
e <- read.table(paste0(mats_dir_human, "/RI.MATS.ReadsOnTargetAndJunctionCounts.txt"), 
                header = TRUE)[, c("GeneID", "PValue", "FDR", "IncLevel1",
                                   "IncLevel2", "IncLevelDifference")]
f <- read.table(paste0(mats_dir_human, "/SE.MATS.ReadsOnTargetAndJunctionCounts.txt"), 
                header = TRUE)[, c("GeneID", "PValue", "FDR", "IncLevel1",
                                   "IncLevel2", "IncLevelDifference")]

x_human <- rbind(a, b, d, e, f)
x_human$IncLevel1_mean <- sapply(as.character(x_human$IncLevel1), function(w) {
  mean(as.numeric(strsplit(w, ",")[[1]]), na.rm = TRUE)
})
x_human$IncLevel2_mean <- sapply(as.character(x_human$IncLevel2), function(w) {
  mean(as.numeric(strsplit(w, ",")[[1]]), na.rm = TRUE)
})

truth_human <- read.delim(truth_file_human, header = TRUE, as.is = TRUE)

m_human <- merge(x_human, truth_human[, c("gene", "ds_status")], 
                 by.x = "GeneID", by.y = "gene", all.x = TRUE)
m_human$signif <- as.numeric(m_human$FDR < 0.05)

m_human$status <- NA
m_human$status[intersect(which(m_human$signif == 1), 
                         which(m_human$ds_status == 1))] <- "TP"
m_human$status[intersect(which(m_human$signif == 1), 
                         which(m_human$ds_status == 0))] <- "FP"
m_human$status[intersect(which(m_human$signif == 0), 
                         which(m_human$ds_status == 1))] <- "FN"
m_human$status[intersect(which(m_human$signif == 0), 
                         which(m_human$ds_status == 0))] <- "TN"

m_human$IncLevelmean <- (m_human$IncLevel1_mean + m_human$IncLevel2_mean)/2
m_human$organism <- "Human"

## DROSOPHILA
a <- read.table(paste0(mats_dir_drosophila, 
                       "/A3SS.MATS.ReadsOnTargetAndJunctionCounts.txt"), 
                header = TRUE)[, c("GeneID", "PValue", "FDR", "IncLevel1",
                                   "IncLevel2", "IncLevelDifference")]
b <- read.table(paste0(mats_dir_drosophila, 
                       "/A5SS.MATS.ReadsOnTargetAndJunctionCounts.txt"), 
                header = TRUE)[, c("GeneID", "PValue", "FDR", "IncLevel1",
                                   "IncLevel2", "IncLevelDifference")]
d <- read.table(paste0(mats_dir_drosophila, 
                       "/MXE.MATS.ReadsOnTargetAndJunctionCounts.txt"), 
                header = TRUE)[, c("GeneID", "PValue", "FDR", "IncLevel1",
                                   "IncLevel2", "IncLevelDifference")]
e <- read.table(paste0(mats_dir_drosophila, 
                       "/RI.MATS.ReadsOnTargetAndJunctionCounts.txt"), 
                header = TRUE)[, c("GeneID", "PValue", "FDR", "IncLevel1",
                                   "IncLevel2", "IncLevelDifference")]
f <- read.table(paste0(mats_dir_drosophila, 
                       "/SE.MATS.ReadsOnTargetAndJunctionCounts.txt"), 
                header = TRUE)[, c("GeneID", "PValue", "FDR", "IncLevel1",
                                   "IncLevel2", "IncLevelDifference")]

x_drosophila <- rbind(a, b, d, e, f)
x_drosophila$IncLevel1_mean <- sapply(as.character(x_drosophila$IncLevel1), function(w) {
  mean(as.numeric(strsplit(w, ",")[[1]]), na.rm = TRUE)
})
x_drosophila$IncLevel2_mean <- sapply(as.character(x_drosophila$IncLevel2), function(w) {
  mean(as.numeric(strsplit(w, ",")[[1]]), na.rm = TRUE)
})

truth_drosophila <- read.delim(truth_file_drosophila, header = TRUE, as.is = TRUE)

m_drosophila <- merge(x_drosophila, truth_drosophila[, c("gene", "ds_status")], 
                 by.x = "GeneID", by.y = "gene", all.x = TRUE)
m_drosophila$signif <- as.numeric(m_drosophila$FDR < 0.05)

m_drosophila$status <- NA
m_drosophila$status[intersect(which(m_drosophila$signif == 1), 
                         which(m_drosophila$ds_status == 1))] <- "TP"
m_drosophila$status[intersect(which(m_drosophila$signif == 1), 
                         which(m_drosophila$ds_status == 0))] <- "FP"
m_drosophila$status[intersect(which(m_drosophila$signif == 0), 
                         which(m_drosophila$ds_status == 1))] <- "FN"
m_drosophila$status[intersect(which(m_drosophila$signif == 0), 
                         which(m_drosophila$ds_status == 0))] <- "TN"

m_drosophila$IncLevelmean <- (m_drosophila$IncLevel1_mean + m_drosophila$IncLevel2_mean)/2
m_drosophila$organism <- "Drosophila"

## PLOT
m <- rbind(m_human, m_drosophila)
library(ggplot2)
pdf(gsub("\\.pdf", "_target_junction.pdf", output_pdf), width = 7, height = 7)
print(ggplot(subset(m, status %in% c("FP", "TN")), 
             aes_string(x = "IncLevelmean", group = "status", col = "status")) + 
        geom_line(stat = "density", size = 2, adjust = 0.1) + 
        scale_color_discrete(name = "") +
        facet_wrap(~organism, nrow = 2) + 
        ##facet_wrap(as.formula(paste("organism ~ ", split_variable)), scales = "free") + 
        xlab("Average inclusion level") + 
        expand_limits(y = 0) + 
        theme(legend.position = "bottom", 
              panel.background = element_rect(fill = NA, colour = "black"),
              panel.border = element_rect(fill = NA, colour = "black", size = 1),
              panel.grid.minor.x = element_blank(),
              panel.grid.minor.y = element_blank(),
              axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, 
                                         size = 15),
              axis.text.y = element_text(size = 15),
              axis.title.x = element_text(size = 15),
              axis.title.y = element_text(size = 15),
              strip.text = element_text(size = 15),
              strip.background = element_rect(fill = NA, colour = "black")))

library(MASS)
library(RColorBrewer)
rf <- colorRampPalette(rev(brewer.pal(11,'Spectral')))
r <- rf(32)
k <- kde2d(subset(m, status == "TN")$IncLevel1_mean, 
           subset(m, status == "TN")$IncLevel2_mean, h = 0.1)
image(k, col = r, main = "True negatives")
k <- kde2d(subset(m, status == "FP")$IncLevel1_mean, 
           subset(m, status == "FP")$IncLevel2_mean, h = 0.1)
image(k, col = r, main = "False positives")
k <- kde2d(subset(m, status == "TP")$IncLevel1_mean, 
           subset(m, status == "TP")$IncLevel2_mean, h = 0.1)
image(k, col = r, main = "True positives")
dev.off()


## ----------------------------- ONLY JUNCTION --------------------------- ##
## HUMAN
a <- read.table(paste0(mats_dir_human, "/A3SS.MATS.JunctionCountOnly.txt"), 
                header = TRUE)[, c("GeneID", "PValue", "FDR", "IncLevel1",
                                   "IncLevel2", "IncLevelDifference")]
b <- read.table(paste0(mats_dir_human, "/A5SS.MATS.JunctionCountOnly.txt"), 
                header = TRUE)[, c("GeneID", "PValue", "FDR", "IncLevel1",
                                   "IncLevel2", "IncLevelDifference")]
d <- read.table(paste0(mats_dir_human, "/MXE.MATS.JunctionCountOnly.txt"), 
                header = TRUE)[, c("GeneID", "PValue", "FDR", "IncLevel1",
                                   "IncLevel2", "IncLevelDifference")]
e <- read.table(paste0(mats_dir_human, "/RI.MATS.JunctionCountOnly.txt"), 
                header = TRUE)[, c("GeneID", "PValue", "FDR", "IncLevel1",
                                   "IncLevel2", "IncLevelDifference")]
f <- read.table(paste0(mats_dir_human, "/SE.MATS.JunctionCountOnly.txt"), 
                header = TRUE)[, c("GeneID", "PValue", "FDR", "IncLevel1",
                                   "IncLevel2", "IncLevelDifference")]

x_human <- rbind(a, b, d, e, f)
x_human$IncLevel1_mean <- sapply(as.character(x_human$IncLevel1), function(w) {
  mean(as.numeric(strsplit(w, ",")[[1]]), na.rm = TRUE)
})
x_human$IncLevel2_mean <- sapply(as.character(x_human$IncLevel2), function(w) {
  mean(as.numeric(strsplit(w, ",")[[1]]), na.rm = TRUE)
})

truth_human <- read.delim(truth_file_human, header = TRUE, as.is = TRUE)

m_human <- merge(x_human, truth_human[, c("gene", "ds_status")], 
                 by.x = "GeneID", by.y = "gene", all.x = TRUE)
m_human$signif <- as.numeric(m_human$FDR < 0.05)

m_human$status <- NA
m_human$status[intersect(which(m_human$signif == 1), 
                         which(m_human$ds_status == 1))] <- "TP"
m_human$status[intersect(which(m_human$signif == 1), 
                         which(m_human$ds_status == 0))] <- "FP"
m_human$status[intersect(which(m_human$signif == 0), 
                         which(m_human$ds_status == 1))] <- "FN"
m_human$status[intersect(which(m_human$signif == 0), 
                         which(m_human$ds_status == 0))] <- "TN"

m_human$IncLevelmean <- (m_human$IncLevel1_mean + m_human$IncLevel2_mean)/2
m_human$organism <- "Human"

## DROSOPHILA
a <- read.table(paste0(mats_dir_drosophila, 
                       "/A3SS.MATS.JunctionCountOnly.txt"), 
                header = TRUE)[, c("GeneID", "PValue", "FDR", "IncLevel1",
                                   "IncLevel2", "IncLevelDifference")]
b <- read.table(paste0(mats_dir_drosophila, 
                       "/A5SS.MATS.JunctionCountOnly.txt"), 
                header = TRUE)[, c("GeneID", "PValue", "FDR", "IncLevel1",
                                   "IncLevel2", "IncLevelDifference")]
d <- read.table(paste0(mats_dir_drosophila, 
                       "/MXE.MATS.JunctionCountOnly.txt"), 
                header = TRUE)[, c("GeneID", "PValue", "FDR", "IncLevel1",
                                   "IncLevel2", "IncLevelDifference")]
e <- read.table(paste0(mats_dir_drosophila, 
                       "/RI.MATS.JunctionCountOnly.txt"), 
                header = TRUE)[, c("GeneID", "PValue", "FDR", "IncLevel1",
                                   "IncLevel2", "IncLevelDifference")]
f <- read.table(paste0(mats_dir_drosophila, 
                       "/SE.MATS.JunctionCountOnly.txt"), 
                header = TRUE)[, c("GeneID", "PValue", "FDR", "IncLevel1",
                                   "IncLevel2", "IncLevelDifference")]

x_drosophila <- rbind(a, b, d, e, f)
x_drosophila$IncLevel1_mean <- sapply(as.character(x_drosophila$IncLevel1), function(w) {
  mean(as.numeric(strsplit(w, ",")[[1]]), na.rm = TRUE)
})
x_drosophila$IncLevel2_mean <- sapply(as.character(x_drosophila$IncLevel2), function(w) {
  mean(as.numeric(strsplit(w, ",")[[1]]), na.rm = TRUE)
})

truth_drosophila <- read.delim(truth_file_drosophila, header = TRUE, as.is = TRUE)

m_drosophila <- merge(x_drosophila, truth_drosophila[, c("gene", "ds_status")], 
                      by.x = "GeneID", by.y = "gene", all.x = TRUE)
m_drosophila$signif <- as.numeric(m_drosophila$FDR < 0.05)

m_drosophila$status <- NA
m_drosophila$status[intersect(which(m_drosophila$signif == 1), 
                              which(m_drosophila$ds_status == 1))] <- "TP"
m_drosophila$status[intersect(which(m_drosophila$signif == 1), 
                              which(m_drosophila$ds_status == 0))] <- "FP"
m_drosophila$status[intersect(which(m_drosophila$signif == 0), 
                              which(m_drosophila$ds_status == 1))] <- "FN"
m_drosophila$status[intersect(which(m_drosophila$signif == 0), 
                              which(m_drosophila$ds_status == 0))] <- "TN"

m_drosophila$IncLevelmean <- (m_drosophila$IncLevel1_mean + m_drosophila$IncLevel2_mean)/2
m_drosophila$organism <- "Drosophila"

## PLOT
m <- rbind(m_human, m_drosophila)
library(ggplot2)
pdf(gsub("\\.pdf", "_junction.pdf", output_pdf), width = 7, height = 7)
print(ggplot(subset(m, status %in% c("FP", "TN")), 
             aes_string(x = "IncLevelmean", group = "status", col = "status")) + 
        geom_line(stat = "density", size = 2, adjust = 0.1) + 
        scale_color_discrete(name = "") +
        facet_wrap(~organism, nrow = 2) + 
        ##facet_wrap(as.formula(paste("organism ~ ", split_variable)), scales = "free") + 
        xlab("Average inclusion level") + 
        expand_limits(y = 0) + 
        theme(legend.position = "bottom", 
              panel.background = element_rect(fill = NA, colour = "black"),
              panel.border = element_rect(fill = NA, colour = "black", size = 1),
              panel.grid.minor.x = element_blank(),
              panel.grid.minor.y = element_blank(),
              axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, 
                                         size = 15),
              axis.text.y = element_text(size = 15),
              axis.title.x = element_text(size = 15),
              axis.title.y = element_text(size = 15),
              strip.text = element_text(size = 15),
              strip.background = element_rect(fill = NA, colour = "black")))

library(MASS)
library(RColorBrewer)
rf <- colorRampPalette(rev(brewer.pal(11,'Spectral')))
r <- rf(32)
k <- kde2d(subset(m, status == "TN")$IncLevel1_mean, 
           subset(m, status == "TN")$IncLevel2_mean, h = 0.1)
image(k, col = r, main = "True negatives")
k <- kde2d(subset(m, status == "FP")$IncLevel1_mean, 
           subset(m, status == "FP")$IncLevel2_mean, h = 0.1)
image(k, col = r, main = "False positives")
k <- kde2d(subset(m, status == "TP")$IncLevel1_mean, 
           subset(m, status == "TP")$IncLevel2_mean, h = 0.1)
image(k, col = r, main = "True positives")
dev.off()
