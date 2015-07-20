## ----- plot_simulation_details_2
## <<plot_simulation_details_2.R>>

## Plot the TPMs of genes/isoforms across all pairs of samples, and color by ds/de status
args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

print(path_to_sim_details)
print(output_pdf)

library(dplyr)
library(ggplot2)
library(scales)

details <- read.delim(path_to_sim_details, header = TRUE, as.is = TRUE)

details$gr1_isoformTPM <- 1/3*(details$s1_isoformTPM + 
                                 details$s2_isoformTPM + 
                                 details$s3_isoformTPM)
details$gr2_isoformTPM <- 1/3*(details$s4_isoformTPM + 
                                 details$s5_isoformTPM + 
                                 details$s6_isoformTPM)
details$transcript_ds_cat <- ifelse(details$transcript_ds_status == 0,
                                    "non-differentially used", 
                                    "differentially used")
x <- details %>% group_by(gene_id) %>% 
  summarise(s1_geneTPM = sum(s1_isoformTPM),
            s2_geneTPM = sum(s2_isoformTPM),
            s3_geneTPM = sum(s3_isoformTPM),
            s4_geneTPM = sum(s4_isoformTPM),
            s5_geneTPM = sum(s5_isoformTPM),
            s6_geneTPM = sum(s6_isoformTPM),
            gr1_geneTPM = 1/3*(s1_geneTPM + 
                                 s2_geneTPM + 
                                 s3_geneTPM),
            gr2_geneTPM = 1/3*(s4_geneTPM + 
                                 s5_geneTPM + 
                                 s6_geneTPM),
            s1_geneCount = sum(s1_isoformCount),
            s2_geneCount = sum(s2_isoformCount),
            s3_geneCount = sum(s3_isoformCount),
            s4_geneCount = sum(s4_isoformCount),
            s5_geneCount = sum(s5_isoformCount),
            s6_geneCount = sum(s6_isoformCount),
            gr1_geneCount = 1/3*(s1_geneCount + 
                                   s2_geneCount + 
                                   s3_geneCount),
            gr2_geneCount = 1/3*(s4_geneCount + 
                                   s5_geneCount + 
                                   s6_geneCount),
            status_ds = gene_ds_status[1]) 

x$status_ds <- ifelse(x$status_ds == 0, "non-differentially spliced", 
                      "differentially spliced")
pdf(output_pdf)
print(ggplot(x, aes(x = gr1_geneCount, y = gr2_geneCount, col = status_ds)) + 
        geom_point(alpha = 1/2) + scale_x_log10() + scale_y_log10() + 
        scale_color_manual(values = c("red", "black"), name = "") +
        guides(colour = guide_legend(override.aes = list(size = 5))) + 
        xlab("Average gene count, condition 1") + 
        ylab("Average gene count, condition 2") +
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
              legend.text=element_text(size = 12),
              strip.background = element_rect(fill = NA, colour = "black")) + 
        ggtitle(""))

print(ggplot(x, aes(x = gr1_geneTPM, y = gr2_geneTPM, col = status_ds)) + 
        geom_point(alpha = 1/2) + scale_x_log10() + scale_y_log10() + 
        scale_color_manual(values = c("red", "black"), name = "") +
        guides(colour = guide_legend(override.aes = list(size = 5))) + 
        xlab("Average gene TPM, condition 1") + 
        ylab("Average gene TPM, condition 2") +
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
              legend.text=element_text(size = 12),
              strip.background = element_rect(fill = NA, colour = "black")) + 
        ggtitle(""))

print(ggplot(details, aes(x = gr1_isoformTPM, y = gr2_isoformTPM, 
                          col = transcript_ds_cat)) + 
        geom_point(alpha = 1/2) + scale_x_log10() + scale_y_log10() + 
        scale_color_manual(values = c("red", "black"), name = "") +
        guides(colour = guide_legend(override.aes = list(size = 5))) + 
        xlab("Average isoform TPM, condition 1") + 
        ylab("Average isoform TPM, condition 2") +
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
              legend.text=element_text(size = 12),
              strip.background = element_rect(fill = NA, colour = "black")) + 
        ggtitle(""))
dev.off()