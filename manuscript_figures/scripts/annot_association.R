## ----- annot_association
## <<annot_association.R>>

library(ggplot2)

truth_human <- read.delim("/home/Shared/data/alt_splicing_simulations/Simulation5_Charlotte/hsapiens/no_diffexpression/non_null_simulation/3_truth/truth_human_non_null_missing20.txt", 
                          header = TRUE, as.is = TRUE)
truth_drosophila <- read.delim("/home/Shared/data/alt_splicing_simulations/Simulation5_Charlotte/drosophila/no_diffexpression/non_null_simulation/3_truth/truth_drosophila_non_null_missing20.txt", 
                          header = TRUE, as.is = TRUE)
truth_human$organism <- "Human"
truth_drosophila$organism <- "Drosophila"
truth <- rbind(truth_human, truth_drosophila)

pdf("/home/Shared/data/alt_splicing_simulations/Simulation5_Charlotte/manuscript_figures/figures/annot_association.pdf", width = 12)
print(ggplot(subset(truth, TPM != 0), aes(x = nbr_isoforms, y = TPM)) + 
        geom_point() + scale_x_log10() + scale_y_log10() + 
        stat_smooth() + facet_wrap(~organism) + 
        xlab("Number of isoforms") + ylab("Gene TPM") + 
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

print(ggplot(subset(truth, TPM != 0 & diff_IsoPct != 0), aes(x = diff_IsoPct, y = TPM)) + 
        geom_point() + scale_y_log10() + scale_x_log10() + 
        stat_smooth() + facet_wrap(~organism) + 
        xlab("Difference in relative abundance between two most dominant isoforms") + 
        ylab("Gene TPM") + 
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



