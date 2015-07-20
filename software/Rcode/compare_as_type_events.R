## ----- compare_as_type_events
## <<compare_as_type_events.R>>

## Read input arguments
args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

print(path_to_astalavista)
print(path_to_simulation_design)
print(output_pdf)

## Read the AStalavista gtf file
astalavista <- read.delim(path_to_astalavista, header = FALSE, as.is = TRUE)
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
colnames(tmp_trs) <- c("transcripts", "as_type")
tmp_trs <- tmp_trs[, c("as_type", "transcripts")]
av <- rbind(av1, tmp_trs)
## Now, av is a many-to-many mapping between transcript pairs and AS events

## Read the simulation design file (contains isoform percentages and 
## which isoforms are simulated to be differentially used)
simdes <- read.delim(path_to_simulation_design, header = TRUE, as.is = TRUE)

library(dplyr)
## Extract the two isoforms with highest abundance for each gene (with >1 isoform)
simdes_top2 <- simdes %>% group_by(gene_id) %>% 
  filter(transcript_id %in% transcript_id[order(IsoPct, decreasing = TRUE)[1:2]]) %>%
  summarise(transcripts = paste(sort(transcript_id), collapse = ";"))
simdes_top2 <- simdes_top2[grep(";", simdes_top2$transcripts), ]
simdes_top2 <- merge(simdes_top2, av, all = TRUE)
simdes_top2 <- subset(simdes_top2, !is.na(gene_id))

## Extract the differentially used isoform pairs
simdes_diffsp <- subset(simdes, transcript_ds_status == 1) %>%
  group_by(gene_id) %>%
  summarise(transcripts = paste(sort(transcript_id), collapse = ";"))
simdes_diffsp <- merge(simdes_diffsp, av, all = TRUE)
simdes_diffsp <- subset(simdes_diffsp, !is.na(gene_id))

combfun <- function(inp_vec) {
  ## Create all combinations of the elements in the input vector
  a <- c()
  for (i in 1:(length(inp_vec) - 1)) {
    for (j in (i + 1):length(inp_vec)) {
      a <- c(a, paste(sort(c(inp_vec[i], inp_vec[j])), collapse = ";"))
    }
  }
  a
}

## Generate all pairs of isoforms expressed above 10%
simdes_tmp <- subset(simdes, IsoPct > 10)
tb <- names(table(simdes_tmp$gene_id)[table(simdes_tmp$gene_id) != 1])
simdes_tmp <- simdes_tmp[simdes_tmp$gene_id %in% tb, ]
simdes_above10 <- lapply(unique(simdes_tmp$gene_id), function(w) {
  combfun(simdes_tmp[simdes_tmp$gene_id == w, "transcript_id"])
})
names(simdes_above10) <- unique(simdes_tmp$gene_id)
simdes_above10 <- stack(simdes_above10)
colnames(simdes_above10) <- c("transcripts", "gene_id")
simdes_above10 <- merge(simdes_above10, av, all = TRUE)
simdes_above10 <- subset(simdes_above10, !is.na(gene_id))

DiffUsed <- sort(table(simdes_diffsp$as_type)/nrow(simdes_diffsp), 
                 decreasing = TRUE)
Top2 <- sort(table(simdes_top2$as_type)/nrow(simdes_top2), 
                  decreasing = TRUE)
All <- sort(table(av$as_type)/nrow(av), decreasing = TRUE)
Above10 <- sort(table(simdes_above10$as_type)/nrow(simdes_above10),
                decreasing = TRUE)

M <- merge(as.data.frame(All), as.data.frame(Top2), by = 0, all = TRUE)
M <- merge(M, as.data.frame(DiffUsed), by.x = "Row.names", by.y = 0, all = TRUE)
M <- merge(M, as.data.frame(Above10), by.x = "Row.names", by.y = 0, all = TRUE)
keep <- union(union(names(DiffUsed[1:20]), names(Top2[1:20])), 
              union(names(All[1:20]), names(Above10[1:20])))
M <- subset(M, Row.names %in% keep)
M <- M[order(rowSums(M[, -1]), decreasing = TRUE), ]
colnames(M)[colnames(M) == "Row.names"] <- "as_type"
M$as_type <- factor(M$as_type, levels = rev(M$as_type))

library(reshape2)
library(ggplot2)
library(data.table)
pdf(output_pdf)
B <- melt(data.table(M), id.vars = "as_type")
B$value <- signif(B$value, digits = 3)
ggplot(B, aes(variable, as_type, fill = value)) + geom_tile() + 
  geom_text(aes(label = value)) + 
  scale_fill_gradient2(limits = c(min(B$value), 1.5*max(B$value))) + 
  scale_x_discrete(expand = c(0, 0)) + 
  scale_y_discrete(expand = c(0, 0)) + 
  xlab("") + ylab("AS event type") + 
  theme(panel.background = element_rect(fill = NA, colour = NA),
        legend.position = "none",
        axis.ticks = element_blank(),
        plot.title = element_text(size = 20),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 10))
dev.off()