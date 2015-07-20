## ----- generate_figures
## <<generate_figures.R>>

library(reshape2)
library(Hmisc)
library(ggplot2)
source("plot_function.R")

basedir <- "/home/Shared/data/alt_splicing_simulations/Simulation5_Charlotte"

## Generate figures for paper
output_directory <- paste0(basedir, "/manuscript_figures/figures")
output_directory_truth <- paste0(basedir, "/manuscript_figures/result_summaries")

## Define path to truth and result files
path_to_results_drosophila <- paste0(basedir, "/drosophila/no_diffexpression/", 
                                     "non_null_simulation/4_results")
path_to_truth_drosophila <- paste0(basedir, "/drosophila/no_diffexpression/", 
                                   "non_null_simulation/3_truth/", 
                                   "truth_drosophila_non_null_missing20.txt")

path_to_results_human <- paste0(basedir, "/hsapiens/no_diffexpression/", 
                                "non_null_simulation/4_results")
path_to_truth_human <- paste0(basedir, "/hsapiens/no_diffexpression/", 
                              "non_null_simulation/3_truth/", 
                              "truth_human_non_null_missing20.txt")

## Define colors and names for method display
methods <- c("dexseq_htseq", "dexseq_htseq_nomerge", 
             "dexseq_featurecounts_flat",  "dexseq_featurecounts_noflat", 
             "dexseq_casper", "dexseq_miso_assignable", "dexseq_splicinggraph",
             "dexseq_tophat_junc", "dexseq_kallisto", "cuffdiff", "rMATS_junction")
method_names <- c("DEXSeq-default", "DEXSeq-noaggreg", "featureCounts-flat", 
                  "featureCounts-exon", "casper", "MISO", "SplicingGraph", 
                  "TopHat-junctions", "kallisto", "cuffdiff", "rMATS")
method_cols <- c(rgb(240, 228, 66, maxColorValue = 255),
                 rgb(0, 158, 115, maxColorValue = 255),
                 rgb(230, 159, 0, maxColorValue = 255),
                 rgb(86, 180, 233, maxColorValue = 255),
                 rgb(0, 0, 0, maxColorValue = 255),
                 rgb(0, 114, 178, maxColorValue = 255),
                 rgb(213, 94, 0, maxColorValue = 255),
                 rgb(204, 121, 167, maxColorValue = 255),
                 rgb(0, 255, 0, maxColorValue = 255),
                 rgb(114, 178, 12, maxColorValue = 255),
                 rgb(123, 0, 119, maxColorValue = 255))

## Define colors and names for different types of filtering
## 5%
methods_5_drosophila <- c("dexseq_htseq_nomerge",
                         "INCOMPLETE_ATLEAST5/dexseq_htseq_nomerge_atleast5",
                         "INCOMPLETE_KALLISTOEST/dexseq_htseq_nomerge_kallistoest_atleast5",
                         "FILTERING/dexseq_htseq_nomerge_binfilt_count_4.018573",
                         "FILTERING/dexseq_htseq_nomerge_binfilt_perc_0.002398478",
                         "FILTERING/dexseq_htseq_nomerge_binfilt_perc_perkb_0.004674837",
                         "FILTERING/dexseq_htseq_nomerge_binfilt_perc_pertx_0.005971978",
                         "FILTERING/dexseq_htseq_nomerge_binfilt_perc_perkb_pertx_0.009185838",
                         "FILTERING/dexseq_htseq_nomerge_binfilt_variance_0.02936164")
methods_5_human <- c("dexseq_htseq_nomerge",
                    "INCOMPLETE_ATLEAST5/dexseq_htseq_nomerge_atleast5",
                    "INCOMPLETE_KALLISTOEST/dexseq_htseq_nomerge_kallistoest_atleast5",
                    "FILTERING/dexseq_htseq_nomerge_binfilt_count_39.10969",
                    "FILTERING/dexseq_htseq_nomerge_binfilt_perc_0.01506810",
                    "FILTERING/dexseq_htseq_nomerge_binfilt_perc_perkb_0.009028732",
                    "FILTERING/dexseq_htseq_nomerge_binfilt_perc_pertx_0.01463351",
                    "FILTERING/dexseq_htseq_nomerge_binfilt_perc_perkb_pertx_0.01032703",
                    "FILTERING/dexseq_htseq_nomerge_binfilt_variance_0.1500974")
method_names_5 <- c("DEXSeq-noaggreg", "DEXSeq-noaggreg (RSEM >5%)",
                    "DEXSeq-noaggreg (kallisto >5%)", "DEXSeq-noaggreg (f1)",
                    "DEXSeq-noaggreg (f2)", "DEXSeq-noaggreg (f3)", 
                    "DEXSeq-noaggreg (f4)", "DEXSeq-noaggreg (f5)",
                    "DEXSeq-noaggreg (f6)")
method_cols_5 <- c(rgb(0, 158, 115, maxColorValue = 255),
                   rgb(0, 0, 0, maxColorValue = 255),
                   rgb(230, 159, 0, maxColorValue = 255),
                   rgb(86, 180, 233, maxColorValue = 255),
                   rgb(0, 114, 178, maxColorValue = 255),
                   rgb(213, 94, 0, maxColorValue = 255),
                   rgb(204, 121, 167, maxColorValue = 255),
                   rgb(0, 255, 0, maxColorValue = 255),
                   rgb(114, 178, 12, maxColorValue = 255))

## 10%
methods_10_drosophila <- c("dexseq_htseq_nomerge",
                          "INCOMPLETE_ATLEAST10/dexseq_htseq_nomerge_atleast10",
                          "INCOMPLETE_KALLISTOEST/dexseq_htseq_nomerge_kallistoest_atleast10",
                          "FILTERING/dexseq_htseq_nomerge_binfilt_count_11.67996",
                          "FILTERING/dexseq_htseq_nomerge_binfilt_perc_0.006635892",
                          "FILTERING/dexseq_htseq_nomerge_binfilt_perc_perkb_0.011305494",
                          "FILTERING/dexseq_htseq_nomerge_binfilt_perc_pertx_0.013075295",
                          "FILTERING/dexseq_htseq_nomerge_binfilt_perc_perkb_pertx_0.015232710",
                          "FILTERING/dexseq_htseq_nomerge_binfilt_variance_0.04505319")
methods_10_human <- c("dexseq_htseq_nomerge",
                     "INCOMPLETE_ATLEAST10/dexseq_htseq_nomerge_atleast10",
                     "INCOMPLETE_KALLISTOEST/dexseq_htseq_nomerge_kallistoest_atleast10",
                     "FILTERING/dexseq_htseq_nomerge_binfilt_count_74.73187",
                     "FILTERING/dexseq_htseq_nomerge_binfilt_perc_0.02045075",
                     "FILTERING/dexseq_htseq_nomerge_binfilt_perc_perkb_0.012398071",
                     "FILTERING/dexseq_htseq_nomerge_binfilt_perc_pertx_0.01828132",
                     "FILTERING/dexseq_htseq_nomerge_binfilt_perc_perkb_pertx_0.01357155",
                     "FILTERING/dexseq_htseq_nomerge_binfilt_variance_0.1773393")
method_names_10 <- c("DEXSeq-noaggreg", "DEXSeq-noaggreg (RSEM >10%)",
                     "DEXSeq-noaggreg (kallisto >10%)", "DEXSeq-noaggreg (f1)",
                     "DEXSeq-noaggreg (f2)", "DEXSeq-noaggreg (f3)", 
                     "DEXSeq-noaggreg (f4)", "DEXSeq-noaggreg (f5)",
                     "DEXSeq-noaggreg (f6)")
method_cols_10 <- c(rgb(0, 158, 115, maxColorValue = 255),
                    rgb(0, 0, 0, maxColorValue = 255),
                    rgb(230, 159, 0, maxColorValue = 255),
                    rgb(86, 180, 233, maxColorValue = 255),
                    rgb(0, 114, 178, maxColorValue = 255),
                    rgb(213, 94, 0, maxColorValue = 255),
                    rgb(204, 121, 167, maxColorValue = 255),
                    rgb(0, 255, 0, maxColorValue = 255),
                    rgb(114, 178, 12, maxColorValue = 255))

## 15%
methods_15_drosophila <- c("dexseq_htseq_nomerge",
                           "INCOMPLETE_ATLEAST15/dexseq_htseq_nomerge_atleast15",
                           "INCOMPLETE_KALLISTOEST/dexseq_htseq_nomerge_kallistoest_atleast15",
                           "FILTERING/dexseq_htseq_nomerge_binfilt_count_20.56182",
                           "FILTERING/dexseq_htseq_nomerge_binfilt_perc_0.011371713",
                           "FILTERING/dexseq_htseq_nomerge_binfilt_perc_perkb_0.017371318",
                           "FILTERING/dexseq_htseq_nomerge_binfilt_perc_pertx_0.018757327",
                           "FILTERING/dexseq_htseq_nomerge_binfilt_perc_perkb_pertx_0.020156113",
                           "FILTERING/dexseq_htseq_nomerge_binfilt_variance_0.05583263")
methods_15_human <- c("dexseq_htseq_nomerge",
                      "INCOMPLETE_ATLEAST15/dexseq_htseq_nomerge_atleast15",
                      "INCOMPLETE_KALLISTOEST/dexseq_htseq_nomerge_kallistoest_atleast15",
                      "FILTERING/dexseq_htseq_nomerge_binfilt_count_108.93658",
                      "FILTERING/dexseq_htseq_nomerge_binfilt_perc_0.02355890",
                      "FILTERING/dexseq_htseq_nomerge_binfilt_perc_perkb_0.014812590",
                      "FILTERING/dexseq_htseq_nomerge_binfilt_perc_pertx_0.02063856",
                      "FILTERING/dexseq_htseq_nomerge_binfilt_perc_perkb_pertx_0.01594104",
                      "FILTERING/dexseq_htseq_nomerge_binfilt_variance_0.2069248")
method_names_15 <- c("DEXSeq-noaggreg", "DEXSeq-noaggreg (RSEM >15%)",
                     "DEXSeq-noaggreg (kallisto >15%)", "DEXSeq-noaggreg (f1)",
                     "DEXSeq-noaggreg (f2)", "DEXSeq-noaggreg (f3)", 
                     "DEXSeq-noaggreg (f4)", "DEXSeq-noaggreg (f5)",
                     "DEXSeq-noaggreg (f6)")
method_cols_15 <- c(rgb(0, 158, 115, maxColorValue = 255),
                    rgb(0, 0, 0, maxColorValue = 255),
                    rgb(230, 159, 0, maxColorValue = 255),
                    rgb(86, 180, 233, maxColorValue = 255),
                    rgb(0, 114, 178, maxColorValue = 255),
                    rgb(213, 94, 0, maxColorValue = 255),
                    rgb(204, 121, 167, maxColorValue = 255),
                    rgb(0, 255, 0, maxColorValue = 255),
                    rgb(114, 178, 12, maxColorValue = 255))

## 25%
methods_25_drosophila <- c("dexseq_htseq_nomerge",
                           "INCOMPLETE_ATLEAST25/dexseq_htseq_nomerge_atleast25",
                           "INCOMPLETE_KALLISTOEST/dexseq_htseq_nomerge_kallistoest_atleast25",
                           "FILTERING/dexseq_htseq_nomerge_binfilt_count_38.93041",
                           "FILTERING/dexseq_htseq_nomerge_binfilt_perc_0.019841270",
                           "FILTERING/dexseq_htseq_nomerge_binfilt_perc_perkb_0.026573817",
                           "FILTERING/dexseq_htseq_nomerge_binfilt_perc_pertx_0.026543975",
                           "FILTERING/dexseq_htseq_nomerge_binfilt_perc_perkb_pertx_0.028110773",
                           "FILTERING/dexseq_htseq_nomerge_binfilt_variance_0.07173037")
methods_25_human <- c("dexseq_htseq_nomerge",
                      "INCOMPLETE_ATLEAST25/dexseq_htseq_nomerge_atleast25",
                      "INCOMPLETE_KALLISTOEST/dexseq_htseq_nomerge_kallistoest_atleast25",
                      "FILTERING/dexseq_htseq_nomerge_binfilt_count_174.93327",
                      "FILTERING/dexseq_htseq_nomerge_binfilt_perc_0.02811174",
                      "FILTERING/dexseq_htseq_nomerge_binfilt_perc_perkb_0.018823718",
                      "FILTERING/dexseq_htseq_nomerge_binfilt_perc_pertx_0.02430609",
                      "FILTERING/dexseq_htseq_nomerge_binfilt_perc_perkb_pertx_0.01988352",
                      "FILTERING/dexseq_htseq_nomerge_binfilt_variance_0.2545928")
method_names_25 <- c("DEXSeq-noaggreg", "DEXSeq-noaggreg (RSEM >25%)",
                     "DEXSeq-noaggreg (kallisto >25%)", "DEXSeq-noaggreg (f1)",
                     "DEXSeq-noaggreg (f2)", "DEXSeq-noaggreg (f3)", 
                     "DEXSeq-noaggreg (f4)", "DEXSeq-noaggreg (f5)",
                     "DEXSeq-noaggreg (f6)")
method_cols_25 <- c(rgb(0, 158, 115, maxColorValue = 255),
                    rgb(0, 0, 0, maxColorValue = 255),
                    rgb(230, 159, 0, maxColorValue = 255),
                    rgb(86, 180, 233, maxColorValue = 255),
                    rgb(0, 114, 178, maxColorValue = 255),
                    rgb(213, 94, 0, maxColorValue = 255),
                    rgb(204, 121, 167, maxColorValue = 255),
                    rgb(0, 255, 0, maxColorValue = 255),
                    rgb(114, 178, 12, maxColorValue = 255))

## Define colors and names for DEXSeq-noaggreg with isoform filtering
methods_filtering <- c("dexseq_htseq_nomerge", 
                       "INCOMPLETE_ATLEAST5/dexseq_htseq_nomerge_atleast5",
                       "INCOMPLETE_ATLEAST10/dexseq_htseq_nomerge_atleast10",
                       "INCOMPLETE_ATLEAST15/dexseq_htseq_nomerge_atleast15",
                       "INCOMPLETE_ATLEAST25/dexseq_htseq_nomerge_atleast25")
method_names_filtering <- c("DEXSeq-noaggreg", "DEXSeq-noaggreg (>5%)",
                            "DEXSeq-noaggreg (>10%)", "DEXSeq-noaggreg (>15%)",
                            "DEXSeq-noaggreg (>25%)")
method_cols_filtering <- c(rgb(0, 158, 115, maxColorValue = 255),
                           rgb(0, 0, 0, maxColorValue = 255),
                           rgb(230, 159, 0, maxColorValue = 255),
                           rgb(86, 180, 233, maxColorValue = 255),
                           rgb(0, 114, 178, maxColorValue = 255))

## Define colors and names for missing annotation display
methods_missing20 <- c("dexseq_htseq", "INCOMPLETE_MISSING20/dexseq_htseq_missing20", 
                       "dexseq_htseq_nomerge", "INCOMPLETE_MISSING20/dexseq_htseq_nomerge_missing20", 
                       "dexseq_featurecounts_flat", "INCOMPLETE_MISSING20/dexseq_featurecounts_flat_missing20", 
                       "dexseq_featurecounts_noflat", "INCOMPLETE_MISSING20/dexseq_featurecounts_noflat_missing20", 
                       "dexseq_casper", "INCOMPLETE_MISSING20/dexseq_casper_missing20", 
                       "dexseq_miso_assignable", "INCOMPLETE_MISSING20/dexseq_miso_assignable_missing20", 
                       "dexseq_splicinggraph", "INCOMPLETE_MISSING20/dexseq_splicinggraph_missing20",
                       "dexseq_tophat_junc", "INCOMPLETE_MISSING20/dexseq_tophat_junc_missing20", 
                       "dexseq_kallisto", "INCOMPLETE_MISSING20/dexseq_kallisto_missing20")
method_names_missing20 <- c("DEXSeq-default", "DEXSeq-default, incomplete",
                            "DEXSeq-noaggreg", "DEXSeq-noaggreg, incomplete", 
                            "featureCounts-flat", "featureCounts-flat, incomplete", 
                            "featureCounts-exon", "featureCounts-exon, incomplete", 
                            "casper", "casper, incomplete", 
                            "MISO", "MISO, incomplete", 
                            "SplicingGraph", "SplicingGraph, incomplete", 
                            "TopHat-junctions", "TopHat-junctions, incomplete", 
                            "kallisto", "kallisto, incomplete")
method_colnames_missing20 <- c("DEXSeq-default", "DEXSeq-default",
                               "DEXSeq-noaggreg", "DEXSeq-noaggreg", 
                               "featureCounts-flat", "featureCounts-flat", 
                               "featureCounts-exon", "featureCounts-exon", 
                               "casper", "casper", 
                               "MISO", "MISO", "SplicingGraph", "SplicingGraph", 
                               "TopHat-junctions", "TopHat-junctions", 
                               "kallisto", "kallisto")
a <- 0
method_cols_missing20 <- c(rgb(240, 228, 66, maxColorValue = 255),
                           rgb(240 + a*(255 - 240), 228 + a*(255 - 228), 
                               66 + 1*(255 - 66), maxColorValue = 255),
                           rgb(0, 158, 115, maxColorValue = 255),
                           rgb(0 + a*(255 - 0), 158 + a*(255 - 158), 
                               115 + a*(255 - 115), maxColorValue = 255),
                           rgb(230, 159, 0, maxColorValue = 255),
                           rgb(230 + a*(255 - 230), 159 + a*(255 - 159), 
                               0 + a*(255 - 0), maxColorValue = 255),
                           rgb(86, 180, 233, maxColorValue = 255),
                           rgb(86 + a*(255 - 86), 180 + a*(255 - 180), 
                               233 + a*(255 - 233), maxColorValue = 255),
                           rgb(0, 0, 0, maxColorValue = 255),
                           rgb(0 + a*(255 - 0), 0 + a*(255 - 0), 
                               0 + a*(255 - 0), maxColorValue = 255),
                           rgb(0, 114, 178, maxColorValue = 255),
                           rgb(0 + a*(255 - 0), 114 + a*(255 - 114), 
                               178 + a*(255 - 178), maxColorValue = 255),
                           rgb(213, 94, 0, maxColorValue = 255),
                           rgb(213 + a*(255 - 213), 94 + a*(255 - 94), 
                               0 + a*(255 - 0), maxColorValue = 255),
                           rgb(204, 121, 167, maxColorValue = 255),
                           rgb(204 + a*(255 - 204), 121 + a*(255 - 121), 
                               167 + a*(255 - 167), maxColorValue = 255),
                           rgb(0, 255, 0, maxColorValue = 255),
                           rgb(0 + a*(255 - 0), 255 + a*(255 - 255), 
                               0 + a*(255 - 0), maxColorValue = 255))
method_pchnames_missing20 <- rep(c("complete", "incomplete"), 9)
method_pch_missing20 <- rep(c(21, 24), 9)
method_pch2_missing20 <- rep(c(19, 17), 9)

## Determine thresholds and load truth files
thresholds <- c(0.01, 0.05, 0.1)
truth_drosophila <- read.delim(path_to_truth_drosophila, header = TRUE, 
                               as.is = TRUE)
truth_human <- read.delim(path_to_truth_human, header = TRUE, as.is = TRUE)

## ========================================================================= ##
## OVERALL PERFORMANCE
truth_drosophila$all <- paste0("(n = ", length(which(!is.na(truth_drosophila$ds_status))), ")")
truth_human$all <- paste0("(n = ", length(which(!is.na(truth_human$ds_status))), ")")

## Original data
plot_fdr_tpr_paper(path_to_results_drosophila = path_to_results_drosophila, 
                   path_to_results_human = path_to_results_human,
                   truth_drosophila = truth_drosophila, 
                   truth_human = truth_human, 
                   methods = methods, method_names = method_names, 
                   method_cols = method_cols, thresholds = thresholds, 
                   split_variable = "all",
                   output_filename = paste0(output_directory, "/fdr_tpr_overall.pdf"),
                   pointsize = 3, stripsize = 10,
                   axistitlesize = 15, axistextsize = 10, fig_height = 5) 

## (consider only genes for which we have both truth and result, e.g. no complexes)
plot_fdr_tpr_paper(path_to_results_drosophila = path_to_results_drosophila, 
                   path_to_results_human = path_to_results_human,
                   truth_drosophila = truth_drosophila, 
                   truth_human = truth_human, 
                   methods = methods, method_names = method_names, 
                   method_cols = method_cols, thresholds = thresholds, 
                   split_variable = "all",
                   output_filename = paste0(output_directory, "/fdr_tpr_overall_onlyshared.pdf"),
                   pointsize = 3, stripsize = 10, fig_height = 5, 
                   axistitlesize = 15, axistextsize = 10, only_shared = TRUE) 

plot_tpr_paper(path_to_results_drosophila = path_to_results_drosophila, 
               path_to_results_human = path_to_results_human,
               truth_drosophila = truth_drosophila, 
               truth_human = truth_human,  
               methods = methods, method_names = method_names, 
               method_cols = method_cols, thresholds = thresholds, 
               split_variable = "all",
               output_filename = paste0(output_directory, "/tpr_overall.pdf"),
               pointsize = 3, stripsize = 10,
               axistitlesize = 15, axistextsize = 10) 

plot_perftable(path_to_results = path_to_results_drosophila, 
               truth = truth_drosophila, 
               methods = methods, method_names = method_names, 
               method_cols = method_cols, threshold = 0.05, 
               split_variable = "all",
               output_filename = paste0(output_directory, "/perf_overall_drosophila.pdf"))

plot_perftable(path_to_results = path_to_results_human, 
               truth = truth_human, 
               methods = methods, method_names = method_names, 
               method_cols = method_cols, threshold = 0.05, 
               split_variable = "all",
               output_filename = paste0(output_directory, "/perf_overall_human.pdf"))

plot_perftable(path_to_results = path_to_results_human, 
               truth = truth_human, 
               methods = methods, method_names = method_names, 
               method_cols = method_cols, threshold = 0.05, 
               split_variable = "all", only_shared = TRUE, 
               output_filename = paste0(output_directory, "/perf_overall_human_onlyshared.pdf"))

## Filtering
plot_fdr_tpr_paper(path_to_results_drosophila = path_to_results_drosophila, 
                   path_to_results_human = path_to_results_human,
                   truth_drosophila = truth_drosophila, 
                   truth_human = truth_human, 
                   methods = methods_filtering, 
                   method_names = method_names_filtering, 
                   method_cols = method_cols_filtering, 
                   thresholds = thresholds, 
                   split_variable = "all",
                   output_filename = paste0(output_directory, 
                                            "/fdr_tpr_overall_filtering.pdf"),
                   pointsize = 3, stripsize = 10,
                   axistitlesize = 15, axistextsize = 10,
                   legend_nrow = 2, fig_height = 5) 

## (consider only genes for which we have both truth and result, so only the genes that actually have a q-value)
plot_fdr_tpr_paper(path_to_results_drosophila = path_to_results_drosophila, 
                   path_to_results_human = path_to_results_human,
                   truth_drosophila = truth_drosophila, 
                   truth_human = truth_human, 
                   methods = methods_filtering, 
                   method_names = method_names_filtering, 
                   method_cols = method_cols_filtering, 
                   thresholds = thresholds, 
                   split_variable = "all",
                   output_filename = paste0(output_directory, 
                                            "/fdr_tpr_overall_filtering_onlyshared.pdf"),
                   pointsize = 3, stripsize = 10,
                   axistitlesize = 15, axistextsize = 10,
                   legend_nrow = 2, only_shared = TRUE, fig_height = 5) 

plot_tpr_paper(path_to_results_drosophila = path_to_results_drosophila, 
               path_to_results_human = path_to_results_human,
               truth_drosophila = truth_drosophila, 
               truth_human = truth_human,  
               methods = methods_filtering, 
               method_names = method_names_filtering, 
               method_cols = method_cols_filtering, thresholds = thresholds, 
               split_variable = "all",
               output_filename = paste0(output_directory, "/tpr_overall_filtering.pdf"),
               pointsize = 3, stripsize = 10,
               axistitlesize = 15, axistextsize = 10) 

plot_perftable(path_to_results = path_to_results_drosophila, 
               truth = truth_drosophila, 
               methods = methods_filtering, 
               method_names = method_names_filtering, 
               method_cols = method_cols_filtering, threshold = 0.05, 
               split_variable = "all",
               output_filename = paste0(output_directory,
                                        "/perf_overall_filtering_drosophila.pdf"))

plot_perftable(path_to_results = path_to_results_human, 
               truth = truth_human, 
               methods = methods_filtering, 
               method_names = method_names_filtering, 
               method_cols = method_cols_filtering, threshold = 0.05, 
               split_variable = "all",
               output_filename = paste0(output_directory, 
                                        "/perf_overall_filtering_human.pdf"))

## Different types of filtering
plot_fdr_tpr_paper(path_to_results_drosophila = path_to_results_drosophila, 
                   path_to_results_human = path_to_results_human,
                   truth_drosophila = truth_drosophila, 
                   truth_human = truth_human, 
                   methods_human = methods_5_human,
                   methods_drosophila = methods_5_drosophila,
                   method_names = method_names_5, 
                   method_cols = method_cols_5, 
                   thresholds = thresholds, 
                   split_variable = "all",
                   output_filename = paste0(output_directory, 
                                            "/fdr_tpr_overall_filtering_types_5.pdf"),
                   pointsize = 3, stripsize = 10,
                   axistitlesize = 15, axistextsize = 10,
                   legend_nrow = 3, fig_height = 5) 

plot_fdr_tpr_paper(path_to_results_drosophila = path_to_results_drosophila, 
                   path_to_results_human = path_to_results_human,
                   truth_drosophila = truth_drosophila, 
                   truth_human = truth_human, 
                   methods_human = methods_10_human,
                   methods_drosophila = methods_10_drosophila,
                   method_names = method_names_10, 
                   method_cols = method_cols_10, 
                   thresholds = thresholds, 
                   split_variable = "all",
                   output_filename = paste0(output_directory, 
                                            "/fdr_tpr_overall_filtering_types_10.pdf"),
                   pointsize = 3, stripsize = 10,
                   axistitlesize = 15, axistextsize = 10,
                   legend_nrow = 3, fig_height = 5) 

plot_fdr_tpr_paper(path_to_results_drosophila = path_to_results_drosophila, 
                   path_to_results_human = path_to_results_human,
                   truth_drosophila = truth_drosophila, 
                   truth_human = truth_human, 
                   methods_human = methods_15_human,
                   methods_drosophila = methods_15_drosophila,
                   method_names = method_names_15, 
                   method_cols = method_cols_15, 
                   thresholds = thresholds, 
                   split_variable = "all",
                   output_filename = paste0(output_directory, 
                                            "/fdr_tpr_overall_filtering_types_15.pdf"),
                   pointsize = 3, stripsize = 10,
                   axistitlesize = 15, axistextsize = 10,
                   legend_nrow = 3, fig_height = 5) 

plot_fdr_tpr_paper(path_to_results_drosophila = path_to_results_drosophila, 
                   path_to_results_human = path_to_results_human,
                   truth_drosophila = truth_drosophila, 
                   truth_human = truth_human, 
                   methods_human = methods_25_human,
                   methods_drosophila = methods_25_drosophila,
                   method_names = method_names_25, 
                   method_cols = method_cols_25, 
                   thresholds = thresholds, 
                   split_variable = "all",
                   output_filename = paste0(output_directory, 
                                            "/fdr_tpr_overall_filtering_types_25.pdf"),
                   pointsize = 3, stripsize = 10,
                   axistitlesize = 15, axistextsize = 10,
                   legend_nrow = 3, fig_height = 5) 


## Incomplete annotation
plot_fdr_tpr_paper_pch(path_to_results_drosophila = path_to_results_drosophila, 
                       path_to_results_human = path_to_results_human,
                       truth_drosophila = truth_drosophila, 
                       truth_human = truth_human, 
                       methods = methods_missing20, 
                       method_names = method_names_missing20, 
                       method_cols = method_cols_missing20, 
                       method_pch = method_pch_missing20, 
                       method_pchnames = method_pchnames_missing20,
                       method_colnames = method_colnames_missing20,
                       thresholds = thresholds, 
                       split_variable = "all",
                       output_filename = paste0(output_directory, 
                                                "/fdr_tpr_overall_incomplete_pch.pdf"),
                       pointsize = 3, stripsize = 10,
                       axistitlesize = 15, axistextsize = 10,
                       legend_nrow = 2, fig_height = 8) 

plot_tpr_paper_pch(path_to_results_drosophila = path_to_results_drosophila, 
                   path_to_results_human = path_to_results_human,
                   truth_drosophila = truth_drosophila, 
                   truth_human = truth_human, 
                   methods = methods_missing20, 
                   method_names = method_names_missing20, 
                   method_cols = method_cols_missing20, 
                   method_pch = method_pch2_missing20, 
                   method_pchnames = method_pchnames_missing20,
                   method_colnames = method_colnames_missing20,
                   thresholds = thresholds, 
                   split_variable = "all",
                   output_filename = paste0(output_directory, 
                                            "/tpr_overall_incomplete_pch.pdf"),
                   pointsize = 3, stripsize = 10,
                   axistitlesize = 15, axistextsize = 10,
                   legend_nrow = 2) 

plot_perftable(path_to_results = path_to_results_drosophila, 
               truth = truth_drosophila, 
               methods = methods_missing20, 
               method_names = method_names_missing20, 
               method_cols = method_cols_missing20, threshold = 0.05, 
               split_variable = "all",
               output_filename = paste0(output_directory, 
                                        "/perf_overall_incomplete_drosophila.pdf"))

plot_perftable(path_to_results = path_to_results_human, 
               truth = truth_human, 
               methods = methods_missing20, 
               method_names = method_names_missing20, 
               method_cols = method_cols_missing20, threshold = 0.05, 
               split_variable = "all",
               output_filename = paste0(output_directory, 
                                        "/perf_overall_incomplete_human.pdf"))

## SPLIT BY NUMBER OF MISSING TRANSCRIPTS
## Number of missing ds transcripts (note that only true ds genes can miss ds transcripts)
truth_drosophila$nbr_missing_2 <- sapply(truth_drosophila$nbr_missing_ds_tr, function(i) {
  paste0(i, " (n = ", length(which(truth_drosophila$nbr_missing_ds_tr == i)), ")")
})
truth_human$nbr_missing_2 <- sapply(truth_human$nbr_missing_ds_tr, function(i) {
  paste0(i, " (n = ", length(which(truth_human$nbr_missing_ds_tr == i)), ")")
})
truth_drosophila$nbr_missing_ds_tr2 <- paste(truth_drosophila$nbr_missing_ds_tr, "du tx missing")
truth_human$nbr_missing_ds_tr2 <- paste(truth_human$nbr_missing_ds_tr, "du tx missing")

plot_tpr_paper_pch(path_to_results_drosophila = path_to_results_drosophila, 
                   path_to_results_human = path_to_results_human,
                   truth_drosophila = truth_drosophila, 
                   truth_human = truth_human, 
                   methods = methods_missing20, 
                   method_names = method_names_missing20, 
                   method_cols = method_cols_missing20, 
                   method_pch = method_pch2_missing20, 
                   method_pchnames = method_pchnames_missing20,
                   method_colnames = method_colnames_missing20,
                   thresholds = thresholds, 
                   split_variable = "nbr_missing_2",
                   output_filename = paste0(output_directory, 
                                            "/tpr_nbr_missing_ds_incomplete_pch.pdf"),
                   pointsize = 1.5, stripsize = 10,
                   axistitlesize = 15, axistextsize = 7,
                   legend_nrow = 2, fig_height = 7) 

plot_perftable(path_to_results = path_to_results_drosophila, 
               truth = truth_drosophila, 
               methods = methods_missing20, 
               method_names = method_names_missing20, 
               method_cols = method_cols_missing20, threshold = 0.05, 
               split_variable = "nbr_missing_ds_tr2",
               output_filename = paste0(output_directory, 
                                        "/perf_nbr_missing_ds_incomplete_drosophila.pdf"))

plot_perftable(path_to_results = path_to_results_human, 
               truth = truth_human, 
               methods = methods_missing20, 
               method_names = method_names_missing20, 
               method_cols = method_cols_missing20, threshold = 0.05, 
               split_variable = "nbr_missing_ds_tr2",
               output_filename = paste0(output_directory, 
                                        "/perf_nbr_missing_ds_incomplete_human.pdf"))

## Number of missing ds and nonds transcripts
truth_drosophila$nbr_missing_nonds_2 <- cut2(truth_drosophila$nbr_missing_nonds_tr,
                                             cuts = c(0, 1, 2))
truth_human$nbr_missing_nonds_2 <- cut2(truth_human$nbr_missing_nonds_tr,
                                        cuts = c(0, 1, 2))
truth_drosophila$nbr_missing_combined <- 
  droplevels(interaction(paste(truth_drosophila$nbr_missing_ds_tr, "ds"), 
                         paste(truth_drosophila$nbr_missing_nonds_2, "nonds")))
truth_human$nbr_missing_combined <- 
  droplevels(interaction(paste(truth_human$nbr_missing_ds_tr, "ds"),
                         paste(truth_human$nbr_missing_nonds_2, "nonds")))
truth_drosophila$nbr_missing_4 <- sapply(truth_drosophila$nbr_missing_combined, function(i) {
  paste0(i, " (n = ", length(which(truth_drosophila$nbr_missing_combined == i)), ")")
})
truth_human$nbr_missing_4 <- sapply(truth_human$nbr_missing_combined, function(i) {
  paste0(i, " (n = ", length(which(truth_human$nbr_missing_combined == i)), ")")
})
plot_fdr_tpr_paper_pch(path_to_results_drosophila = path_to_results_drosophila, 
                       path_to_results_human = NULL,
                       truth_drosophila = truth_drosophila, 
                       truth_human = NULL, 
                       methods = methods_missing20, 
                       method_names = method_names_missing20, 
                       method_cols = method_cols_missing20, 
                       method_pch = method_pch_missing20, 
                       method_pchnames = method_pchnames_missing20,
                       method_colnames = method_colnames_missing20,
                       thresholds = thresholds, 
                       split_variable = "nbr_missing_4",
                       output_filename = paste0(output_directory, 
                                                "/fdr_tpr_nbr_missing_ds_nonds", 
                                                "_incomplete_pch_drosophila.pdf"),
                       pointsize = 2, stripsize = 10,
                       axistitlesize = 10, axistextsize = 10,
                       legend_nrow = 2) 

plot_tpr_paper_pch(path_to_results_drosophila = path_to_results_drosophila, 
                   path_to_results_human = NULL,
                   truth_drosophila = truth_drosophila, 
                   truth_human = NULL, 
                   methods = methods_missing20, 
                   method_names = method_names_missing20, 
                   method_cols = method_cols_missing20, 
                   method_pch = method_pch2_missing20, 
                   method_pchnames = method_pchnames_missing20,
                   method_colnames = method_colnames_missing20,
                   thresholds = thresholds, 
                   split_variable = "nbr_missing_4",
                   output_filename = paste0(output_directory, 
                                            "/tpr_nbr_missing_ds_nonds", 
                                            "_incomplete_pch_drosophila.pdf"),
                   pointsize = 1.5, stripsize = 10,
                   axistitlesize = 15, axistextsize = 7,
                   legend_nrow = 2) 

plot_perftable(path_to_results = path_to_results_drosophila, 
               truth = truth_drosophila, 
               methods = methods_missing20, 
               method_names = method_names_missing20, 
               method_cols = method_cols_missing20, threshold = 0.05, 
               split_variable = "nbr_missing_combined",
               output_filename = paste0(output_directory, 
                                        "/perf_nbr_missing_ds_nonds", 
                                        "_incomplete_drosophila.pdf"))

plot_fdr_tpr_paper_pch(path_to_results_drosophila = NULL, 
                       path_to_results_human = path_to_results_human,
                       truth_drosophila = NULL, 
                       truth_human = truth_human, 
                       methods = methods_missing20, 
                       method_names = method_names_missing20, 
                       method_cols = method_cols_missing20, 
                       method_pch = method_pch_missing20, 
                       method_pchnames = method_pchnames_missing20,
                       method_colnames = method_colnames_missing20,
                       thresholds = thresholds, 
                       split_variable = "nbr_missing_4",
                       output_filename = paste0(output_directory, 
                                                "/fdr_tpr_nbr_missing_ds_nonds", 
                                                "_incomplete_pch_human.pdf"),
                       pointsize = 2, stripsize = 10,
                       axistitlesize = 10, axistextsize = 10,
                       legend_nrow = 2) 

plot_tpr_paper_pch(path_to_results_drosophila = NULL, 
                   path_to_results_human = path_to_results_human,
                   truth_drosophila = NULL, 
                   truth_human = truth_human, 
                   methods = methods_missing20, 
                   method_names = method_names_missing20, 
                   method_cols = method_cols_missing20, 
                   method_pch = method_pch2_missing20, 
                   method_pchnames = method_pchnames_missing20,
                   method_colnames = method_colnames_missing20,
                   thresholds = thresholds, 
                   split_variable = "nbr_missing_4",
                   output_filename = paste0(output_directory, 
                                            "/tpr_nbr_missing_ds_nonds", 
                                            "_incomplete_pch_human.pdf"),
                   pointsize = 1.5, stripsize = 10,
                   axistitlesize = 15, axistextsize = 7,
                   legend_nrow = 2) 

plot_perftable(path_to_results = path_to_results_human, 
               truth = truth_human, 
               methods = methods_missing20, 
               method_names = method_names_missing20, 
               method_cols = method_cols_missing20, threshold = 0.05, 
               split_variable = "nbr_missing_combined",
               output_filename = paste0(output_directory, 
                                        "/perf_nbr_missing_ds_nonds", 
                                        "_incomplete_human.pdf"))

## Number of remaining transcripts
truth_drosophila$nbr_remaining_2 <- cut2(truth_drosophila$nbr_remaining_tr,
                                         cuts = c(0, 1, 2))
truth_human$nbr_remaining_2 <- cut2(truth_human$nbr_remaining_tr,
                                    cuts = c(0, 1, 2))
truth_drosophila$nbr_remaining_3 <- sapply(truth_drosophila$nbr_remaining_2, function(i) {
  paste0(i, " (n = ", length(which(truth_drosophila$nbr_remaining_2 == i)), ")")
})
truth_human$nbr_remaining_3 <- sapply(truth_human$nbr_remaining_2, function(i) {
  paste0(i, " (n = ", length(which(truth_human$nbr_remaining_2 == i)), ")")
})
plot_fdr_tpr_paper_pch(path_to_results_drosophila = path_to_results_drosophila, 
                       path_to_results_human = path_to_results_human,
                       truth_drosophila = truth_drosophila, 
                       truth_human = truth_human, 
                       methods = methods_missing20, 
                       method_names = method_names_missing20, 
                       method_cols = method_cols_missing20, 
                       method_pch = method_pch_missing20, 
                       method_pchnames = method_pchnames_missing20,
                       method_colnames = method_colnames_missing20,
                       thresholds = thresholds, 
                       split_variable = "nbr_remaining_3",
                       output_filename = paste0(output_directory, 
                                                "/fdr_tpr_nbr_rem", 
                                                "_incomplete_pch.pdf"),
                       pointsize = 2, stripsize = 10,
                       axistitlesize = 15, axistextsize = 10,
                       legend_nrow = 2) 

plot_tpr_paper_pch(path_to_results_drosophila = path_to_results_drosophila, 
                   path_to_results_human = path_to_results_human,
                   truth_drosophila = truth_drosophila, 
                   truth_human = truth_human, 
                   methods = methods_missing20, 
                   method_names = method_names_missing20, 
                   method_cols = method_cols_missing20, 
                   method_pch = method_pch2_missing20, 
                   method_pchnames = method_pchnames_missing20,
                   method_colnames = method_colnames_missing20,
                   thresholds = thresholds, 
                   split_variable = "nbr_remaining_3",
                   output_filename = paste0(output_directory, 
                                            "/tpr_nbr_rem", 
                                            "_incomplete_pch.pdf"),
                   pointsize = 2, stripsize = 10,
                   axistitlesize = 15, axistextsize = 7,
                   legend_nrow = 2) 

plot_perftable(path_to_results = path_to_results_drosophila, 
               truth = truth_drosophila, 
               methods = methods_missing20, 
               method_names = method_names_missing20, 
               method_cols = method_cols_missing20, threshold = 0.05, 
               split_variable = "nbr_remaining_2",
               output_filename = paste0(output_directory, 
                                        "/perf_nbr_rem", 
                                        "_incomplete_drosophila.pdf"))

plot_perftable(path_to_results = path_to_results_human, 
               truth = truth_human, 
               methods = methods_missing20, 
               method_names = method_names_missing20, 
               method_cols = method_cols_missing20, threshold = 0.05, 
               split_variable = "nbr_remaining_2",
               output_filename = paste0(output_directory, 
                                        "/perf_nbr_rem", 
                                        "_incomplete_human.pdf"))

## Number of missing ds and remaining transcripts
truth_drosophila$nbr_missing_combined_2 <- 
  droplevels(interaction(paste(truth_drosophila$nbr_missing_ds_tr, "ds"), 
                         paste(truth_drosophila$nbr_remaining_2, "rem")))
truth_human$nbr_missing_combined_2 <- 
  droplevels(interaction(paste(truth_human$nbr_missing_ds_tr, "ds"),
                         paste(truth_human$nbr_remaining_2, "rem")))
truth_drosophila$nbr_missing_5 <- sapply(truth_drosophila$nbr_missing_combined_2, function(i) {
  paste0(i, " (n = ", length(which(truth_drosophila$nbr_missing_combined_2 == i)), ")")
})
truth_human$nbr_missing_5 <- sapply(truth_human$nbr_missing_combined_2, function(i) {
  paste0(i, " (n = ", length(which(truth_human$nbr_missing_combined_2 == i)), ")")
})
plot_fdr_tpr_paper_pch(path_to_results_drosophila = path_to_results_drosophila, 
                       path_to_results_human = NULL,
                       truth_drosophila = subset(truth_drosophila, 
                                                 !(nbr_missing_combined_2 %in% 
                                                     c("0 ds. 0 rem", "0 ds. 1 rem"))), 
                       truth_human = NULL, 
                       methods = methods_missing20, 
                       method_names = method_names_missing20, 
                       method_cols = method_cols_missing20, 
                       method_pch = method_pch_missing20, 
                       method_pchnames = method_pchnames_missing20,
                       method_colnames = method_colnames_missing20,
                       thresholds = thresholds, 
                       split_variable = "nbr_missing_5",
                       output_filename = paste0(output_directory, 
                                                "/fdr_tpr_nbr_missing_ds_nbr_rem", 
                                                "_incomplete_pch_drosophila.pdf"),
                       pointsize = 2, stripsize = 10,
                       axistitlesize = 15, axistextsize = 10,
                       legend_nrow = 2) 

plot_tpr_paper_pch(path_to_results_drosophila = path_to_results_drosophila, 
                   path_to_results_human = NULL,
                   truth_drosophila = subset(truth_drosophila, 
                                             !(nbr_missing_combined_2 %in% 
                                                 c("0 ds. 0 rem", "0 ds. 1 rem"))), 
                   truth_human = NULL, 
                   methods = methods_missing20, 
                   method_names = method_names_missing20, 
                   method_cols = method_cols_missing20, 
                   method_pch = method_pch2_missing20, 
                   method_pchnames = method_pchnames_missing20,
                   method_colnames = method_colnames_missing20,
                   thresholds = thresholds, 
                   split_variable = "nbr_missing_5",
                   output_filename = paste0(output_directory, 
                                            "/tpr_nbr_missing_ds_nbr_rem", 
                                            "_incomplete_pch_drosophila.pdf"),
                   pointsize = 1.5, stripsize = 10,
                   axistitlesize = 15, axistextsize = 7,
                   legend_nrow = 2) 

plot_perftable(path_to_results = path_to_results_drosophila, 
               truth = truth_drosophila, 
               methods = methods_missing20, 
               method_names = method_names_missing20, 
               method_cols = method_cols_missing20, threshold = 0.05, 
               split_variable = "nbr_missing_combined_2",
               output_filename = paste0(output_directory, 
                                        "/perf_nbr_missing_ds_nbr_rem", 
                                        "_incomplete_drosophila.pdf"))

plot_fdr_tpr_paper_pch(path_to_results_drosophila = NULL, 
                       path_to_results_human = path_to_results_human,
                       truth_drosophila = NULL, 
                       truth_human = subset(truth_human, 
                                            !(nbr_missing_combined_2 %in% 
                                                c("0 ds. 0 rem", "0 ds. 1 rem"))), 
                       methods = methods_missing20, 
                       method_names = method_names_missing20, 
                       method_cols = method_cols_missing20, 
                       method_pch = method_pch_missing20, 
                       method_pchnames = method_pchnames_missing20,
                       method_colnames = method_colnames_missing20,
                       thresholds = thresholds, 
                       split_variable = "nbr_missing_5",
                       output_filename = paste0(output_directory, 
                                                "/fdr_tpr_nbr_missing_ds_nbr_rem", 
                                                "_incomplete_pch_human.pdf"),
                       pointsize = 2, stripsize = 10,
                       axistitlesize = 15, axistextsize = 10,
                       legend_nrow = 2) 

plot_tpr_paper_pch(path_to_results_drosophila = NULL, 
                   path_to_results_human = path_to_results_human,
                   truth_drosophila = NULL, 
                   truth_human = subset(truth_human, 
                                        !(nbr_missing_combined_2 %in% 
                                            c("0 ds. 0 rem", "0 ds. 1 rem"))), 
                   methods = methods_missing20, 
                   method_names = method_names_missing20, 
                   method_cols = method_cols_missing20, 
                   method_pch = method_pch2_missing20, 
                   method_pchnames = method_pchnames_missing20,
                   method_colnames = method_colnames_missing20,
                   thresholds = thresholds, 
                   split_variable = "nbr_missing_5",
                   output_filename = paste0(output_directory, 
                                            "/tpr_nbr_missing_ds_nbr_rem", 
                                            "_incomplete_pch_human.pdf"),
                   pointsize = 1.5, stripsize = 10,
                   axistitlesize = 15, axistextsize = 7,
                   legend_nrow = 2) 

plot_perftable(path_to_results = path_to_results_human, 
               truth = truth_human, 
               methods = methods_missing20, 
               method_names = method_names_missing20, 
               method_cols = method_cols_missing20, threshold = 0.05, 
               split_variable = "nbr_missing_combined_2",
               output_filename = paste0(output_directory, 
                                        "/perf_nbr_missing_ds_nbr_rem", 
                                        "_incomplete_human.pdf"))

## SPLIT BY ISOFORM PERCENTAGE DIFFERENCE
## Original data
truth_drosophila$diff_IsoPct_3 <- cut2(truth_drosophila$diff_IsoPct, 
                                       cuts = c(0, 1/3, 2/3, 1))
truth_human$diff_IsoPct_3 <- cut2(truth_human$diff_IsoPct, 
                                  cuts = c(0, 1/3, 2/3, 1))
truth_drosophila$diff_IsoPct3 <- sapply(truth_drosophila$diff_IsoPct_3, function(i) {
  paste0(i, " (n = ", length(which(truth_drosophila$diff_IsoPct_3 == i)), ")")
})
truth_human$diff_IsoPct3 <- sapply(truth_human$diff_IsoPct_3, function(i) {
  paste0(i, " (n = ", length(which(truth_human$diff_IsoPct_3 == i)), ")")
})
truth_drosophila$diff_IsoPct3_2 <- sapply(truth_drosophila$diff_IsoPct_3, function(i) {
  paste0(i, " (n = ", length(which(truth_drosophila$diff_IsoPct_3 == i)), ", n.ds = ",
         length(intersect(which(truth_drosophila$diff_IsoPct_3 == i), 
                          which(truth_drosophila$ds_status == 1))), ")")
})
truth_human$diff_IsoPct3_2 <- sapply(truth_human$diff_IsoPct_3, function(i) {
  paste0(i, " (n = ", length(which(truth_human$diff_IsoPct_3 == i)), ", n.ds = ",
         length(intersect(which(truth_human$diff_IsoPct_3 == i), 
                          which(truth_human$ds_status == 1))), ")")
})
plot_fdr_tpr_paper(path_to_results_drosophila = path_to_results_drosophila, 
                   path_to_results_human = path_to_results_human,
                   truth_drosophila = truth_drosophila, 
                   truth_human = truth_human, 
                   methods = methods, method_names = method_names, 
                   method_cols = method_cols, thresholds = thresholds, 
                   split_variable = "diff_IsoPct3",
                   output_filename = paste0(output_directory, "/fdr_tpr_diffisopct.pdf"),
                   pointsize = 3, stripsize = 10,
                   axistitlesize = 15, axistextsize = 10) 

plot_fdr_tpr_paper(path_to_results_drosophila = path_to_results_drosophila, 
                   path_to_results_human = path_to_results_human,
                   truth_drosophila = truth_drosophila, 
                   truth_human = truth_human, 
                   methods = methods, method_names = method_names, 
                   method_cols = method_cols, thresholds = thresholds, 
                   split_variable = "diff_IsoPct3_2",
                   output_filename = paste0(output_directory, "/fdr_tpr_diffisopct_withnds.pdf"),
                   pointsize = 3, stripsize = 9,
                   axistitlesize = 15, axistextsize = 10) 

plot_tpr_paper(path_to_results_drosophila = path_to_results_drosophila, 
               path_to_results_human = path_to_results_human,
               truth_drosophila = truth_drosophila, 
               truth_human = truth_human, 
               methods = methods, method_names = method_names, 
               method_cols = method_cols, thresholds = thresholds, 
               split_variable = "diff_IsoPct3",
               output_filename = paste0(output_directory, "/tpr_diffisopct.pdf"),
               pointsize = 3, stripsize = 10,
               axistitlesize = 15, axistextsize = 10) 

plot_perftable(path_to_results = path_to_results_drosophila, 
               truth = truth_drosophila, 
               methods = methods, method_names = method_names, 
               method_cols = method_cols, threshold = 0.05, 
               split_variable = "diff_IsoPct_3",
               output_filename = paste0(output_directory, "/perf_diffisopct_drosophila.pdf"))

plot_perftable(path_to_results = path_to_results_human, 
               truth = truth_human, 
               methods = methods, method_names = method_names, 
               method_cols = method_cols, threshold = 0.05, 
               split_variable = "diff_IsoPct_3",
               output_filename = paste0(output_directory, "/perf_diffisopct_human.pdf"))

plot_counting_characteristics(path_to_results_drosophila = path_to_results_drosophila,
                              path_to_results_human = path_to_results_human, 
                              truth_drosophila = truth_drosophila,
                              truth_human = truth_human, 
                              showmeasure = "dispvar", 
                              labelmeasure = "within-gene variance(log2(dispersion))", 
                              method = "dexseq_htseq_nomerge", 
                              method_name = "DEXSeq-noaggreg", 
                              threshold = 0.05, split_variable = "diff_IsoPct3",
                              output_filename = paste0(output_directory, "/counting_diffisopct_dexseq_noaggreg.pdf"), 
                              stripsize = 10,
                              axistitlesize = 15, axistextsize = 10,
                              legend_nrow = 2, linewidth = 1.5)

## Filtering
plot_fdr_tpr_paper(path_to_results_drosophila = path_to_results_drosophila, 
                   path_to_results_human = path_to_results_human,
                   truth_drosophila = truth_drosophila, 
                   truth_human = truth_human, 
                   methods = methods_filtering, 
                   method_names = method_names_filtering, 
                   method_cols = method_cols_filtering, thresholds = thresholds, 
                   split_variable = "diff_IsoPct3",
                   output_filename = paste0(output_directory, 
                                            "/fdr_tpr_diffisopct_filtering.pdf"),
                   pointsize = 3, stripsize = 10,
                   axistitlesize = 15, axistextsize = 10) 

plot_tpr_paper(path_to_results_drosophila = path_to_results_drosophila, 
               path_to_results_human = path_to_results_human,
               truth_drosophila = truth_drosophila, 
               truth_human = truth_human, 
               methods = methods_filtering, 
               method_names = method_names_filtering, 
               method_cols = method_cols_filtering, thresholds = thresholds, 
               split_variable = "diff_IsoPct3",
               output_filename = paste0(output_directory, 
                                        "/tpr_diffisopct_filtering.pdf"),
               pointsize = 3, stripsize = 10,
               axistitlesize = 15, axistextsize = 10) 

plot_perftable(path_to_results = path_to_results_drosophila, 
               truth = truth_drosophila, 
               methods = methods_filtering, 
               method_names = method_names_filtering, 
               method_cols = method_cols, threshold = 0.05, 
               split_variable = "diff_IsoPct_3",
               output_filename = paste0(output_directory, 
                                        "/perf_diffisopct_filtering_drosophila.pdf"))

plot_perftable(path_to_results = path_to_results_human, 
               truth = truth_human, 
               methods = methods_filtering, 
               method_names = method_names_filtering, 
               method_cols = method_cols_filtering, threshold = 0.05, 
               split_variable = "diff_IsoPct_3",
               output_filename = paste0(output_directory, 
                                        "/perf_diffisopct_filtering_human.pdf"))

plot_counting_characteristics(path_to_results_drosophila = path_to_results_drosophila,
                              path_to_results_human = path_to_results_human, 
                              truth_drosophila = truth_drosophila,
                              truth_human = truth_human, 
                              showmeasure = "dispvar", 
                              labelmeasure = "within-group variance(log2(dispersion))", 
                              method = "INCOMPLETE_ATLEAST10/dexseq_htseq_nomerge_atleast10", 
                              method_name = "DEXSeq-noaggreg (>10%)", 
                              threshold = 0.05, split_variable = "diff_IsoPct3",
                              output_filename = paste0(output_directory, "/counting_diffisopct_dexseq_noaggreg_atleast10.pdf"), 
                              stripsize = 10,
                              axistitlesize = 15, axistextsize = 10,
                              legend_nrow = 2, linewidth = 1.5)

## Different types of filtering
plot_fdr_tpr_paper(path_to_results_drosophila = path_to_results_drosophila, 
                   path_to_results_human = path_to_results_human,
                   truth_drosophila = truth_drosophila, 
                   truth_human = truth_human, 
                   methods_human = methods_5_human,
                   methods_drosophila = methods_5_drosophila,
                   method_names = method_names_5, 
                   method_cols = method_cols_5, 
                   thresholds = thresholds, 
                   split_variable = "diff_IsoPct3",
                   output_filename = paste0(output_directory, 
                                            "/fdr_tpr_diffisopct_filtering_types_5.pdf"),
                   pointsize = 3, stripsize = 10,
                   axistitlesize = 15, axistextsize = 10,
                   legend_nrow = 3) 

plot_fdr_tpr_paper(path_to_results_drosophila = path_to_results_drosophila, 
                   path_to_results_human = path_to_results_human,
                   truth_drosophila = truth_drosophila, 
                   truth_human = truth_human, 
                   methods_human = methods_10_human,
                   methods_drosophila = methods_10_drosophila,
                   method_names = method_names_10, 
                   method_cols = method_cols_10, 
                   thresholds = thresholds, 
                   split_variable = "diff_IsoPct3",
                   output_filename = paste0(output_directory, 
                                            "/fdr_tpr_diffisopct_filtering_types_10.pdf"),
                   pointsize = 3, stripsize = 10,
                   axistitlesize = 15, axistextsize = 10,
                   legend_nrow = 3) 

plot_fdr_tpr_paper(path_to_results_drosophila = path_to_results_drosophila, 
                   path_to_results_human = path_to_results_human,
                   truth_drosophila = truth_drosophila, 
                   truth_human = truth_human, 
                   methods_human = methods_15_human,
                   methods_drosophila = methods_15_drosophila,
                   method_names = method_names_15, 
                   method_cols = method_cols_15, 
                   thresholds = thresholds, 
                   split_variable = "diff_IsoPct3",
                   output_filename = paste0(output_directory, 
                                            "/fdr_tpr_diffisopct_filtering_types_15.pdf"),
                   pointsize = 3, stripsize = 10,
                   axistitlesize = 15, axistextsize = 10,
                   legend_nrow = 3) 

plot_fdr_tpr_paper(path_to_results_drosophila = path_to_results_drosophila, 
                   path_to_results_human = path_to_results_human,
                   truth_drosophila = truth_drosophila, 
                   truth_human = truth_human, 
                   methods_human = methods_25_human,
                   methods_drosophila = methods_25_drosophila,
                   method_names = method_names_25, 
                   method_cols = method_cols_25, 
                   thresholds = thresholds, 
                   split_variable = "diff_IsoPct3",
                   output_filename = paste0(output_directory, 
                                            "/fdr_tpr_diffisopct_filtering_types_25.pdf"),
                   pointsize = 3, stripsize = 10,
                   axistitlesize = 15, axistextsize = 10,
                   legend_nrow = 3) 



## SPLIT BY (GENE) EXPRESSION LEVEL
truth_drosophila$exprlevel <- "low"
truth_drosophila$exprlevel[truth_drosophila$TPM > 
                             median(truth_drosophila$TPM[truth_drosophila$ds_status == 1])] <- "high"
truth_drosophila$exprlevel <- factor(truth_drosophila$exprlevel, levels = c("low", "high"))
truth_human$exprlevel <- "low"
truth_human$exprlevel[truth_human$TPM > median(truth_human$TPM[truth_human$ds_status == 1])] <- "high"
truth_human$exprlevel <- factor(truth_human$exprlevel, levels = c("low", "high"))
truth_drosophila$exprlevel2 <- sapply(truth_drosophila$exprlevel, function(i) {
  paste0(i, " (n = ", length(which(truth_drosophila$exprlevel == i)), ")")
})
truth_human$exprlevel2 <- sapply(truth_human$exprlevel, function(i) {
  paste0(i, " (n = ", length(which(truth_human$exprlevel == i)), ")")
})
## Original data
plot_fdr_tpr_paper(path_to_results_drosophila = path_to_results_drosophila, 
                   path_to_results_human = path_to_results_human,
                   truth_drosophila = truth_drosophila, 
                   truth_human = truth_human, 
                   methods = methods, method_names = method_names, 
                   method_cols = method_cols, thresholds = thresholds, 
                   split_variable = "exprlevel2",
                   output_filename = paste0(output_directory, "/fdr_tpr_TPM.pdf"),
                   pointsize = 3, stripsize = 10,
                   axistitlesize = 15, axistextsize = 10) 

plot_tpr_paper(path_to_results_drosophila = path_to_results_drosophila, 
               path_to_results_human = path_to_results_human,
               truth_drosophila = truth_drosophila, 
               truth_human = truth_human, 
               methods = methods, method_names = method_names, 
               method_cols = method_cols, thresholds = thresholds, 
               split_variable = "exprlevel2",
               output_filename = paste0(output_directory, "/tpr_TPM.pdf"),
               pointsize = 3, stripsize = 10,
               axistitlesize = 15, axistextsize = 10) 

plot_perftable(path_to_results = path_to_results_drosophila, 
               truth = truth_drosophila, 
               methods = methods, method_names = method_names, 
               method_cols = method_cols, threshold = 0.05, 
               split_variable = "exprlevel",
               output_filename = paste0(output_directory, "/perf_TPM_drosophila.pdf"))

plot_perftable(path_to_results = path_to_results_human, 
               truth = truth_human, 
               methods = methods, method_names = method_names, 
               method_cols = method_cols, threshold = 0.05, 
               split_variable = "exprlevel",
               output_filename = paste0(output_directory, "/perf_TPM_human.pdf"))

## Filtering
plot_fdr_tpr_paper(path_to_results_drosophila = path_to_results_drosophila, 
                   path_to_results_human = path_to_results_human,
                   truth_drosophila = truth_drosophila, 
                   truth_human = truth_human, 
                   methods = methods_filtering, 
                   method_names = method_names_filtering, 
                   method_cols = method_cols_filtering, thresholds = thresholds, 
                   split_variable = "exprlevel2",
                   output_filename = paste0(output_directory, 
                                            "/fdr_tpr_TPM_filtering.pdf"),
                   pointsize = 3, stripsize = 10,
                   axistitlesize = 15, axistextsize = 10) 

plot_tpr_paper(path_to_results_drosophila = path_to_results_drosophila, 
               path_to_results_human = path_to_results_human,
               truth_drosophila = truth_drosophila, 
               truth_human = truth_human, 
               methods = methods_filtering, 
               method_names = method_names_filtering, 
               method_cols = method_cols_filtering, thresholds = thresholds, 
               split_variable = "exprlevel2",
               output_filename = paste0(output_directory, 
                                        "/tpr_TPM_filtering.pdf"),
               pointsize = 3, stripsize = 10,
               axistitlesize = 15, axistextsize = 10) 

plot_perftable(path_to_results = path_to_results_drosophila, 
               truth = truth_drosophila, 
               methods = methods_filtering, 
               method_names = method_names_filtering, 
               method_cols = method_cols, threshold = 0.05, 
               split_variable = "exprlevel",
               output_filename = paste0(output_directory, 
                                        "/perf_TPM_filtering_drosophila.pdf"))

plot_perftable(path_to_results = path_to_results_human, 
               truth = truth_human, 
               methods = methods_filtering, 
               method_names = method_names_filtering, 
               method_cols = method_cols_filtering, threshold = 0.05, 
               split_variable = "exprlevel",
               output_filename = paste0(output_directory, 
                                        "/perf_TPM_filtering_human.pdf"))

# ## SPLIT BY NUMBER OF ISOFORMS WITH ABUNDANCE > 5%
# truth_drosophila$nbr_isoforms_a5 <- cut2(truth_drosophila$nbr_isoforms_above5, 
#                                          cuts = c(2, 3, max(truth_drosophila$nbr_isoforms_above5)),
#                                          oneval = FALSE)
# truth_human$nbr_isoforms_a5 <- cut2(truth_human$nbr_isoforms_above5, 
#                                     cuts = c(2, 3, max(truth_human$nbr_isoforms_above5)),
#                                     oneval = FALSE)
# truth_drosophila$nbr_isoforms_a5_2 <- sapply(truth_drosophila$nbr_isoforms_a5, function(i) {
#   paste0(i, " (n = ", length(which(truth_drosophila$nbr_isoforms_a5 == i)), ")")
# })
# truth_human$nbr_isoforms_a5_2 <- sapply(truth_human$nbr_isoforms_a5, function(i) {
#   paste0(i, " (n = ", length(which(truth_human$nbr_isoforms_a5 == i)), ")")
# })
# ## Original data
# plot_fdr_tpr_paper(path_to_results_drosophila = path_to_results_drosophila, 
#                    path_to_results_human = path_to_results_human,
#                    truth_drosophila = truth_drosophila, 
#                    truth_human = truth_human,  
#                    methods = methods, method_names = method_names, 
#                    method_cols = method_cols, thresholds = thresholds, 
#                    split_variable = "nbr_isoforms_a5_2",
#                    output_filename = paste0(output_directory, "/fdr_tpr_nbr_isoforms_above5.pdf"),
#                    pointsize = 3, stripsize = 10,
#                    axistitlesize = 15, axistextsize = 10) 
# 
# plot_tpr_paper(path_to_results_drosophila = path_to_results_drosophila, 
#                path_to_results_human = path_to_results_human,
#                truth_drosophila = truth_drosophila, 
#                truth_human = truth_human,  
#                methods = methods, method_names = method_names, 
#                method_cols = method_cols, thresholds = thresholds, 
#                split_variable = "nbr_isoforms_a5_2",
#                output_filename = paste0(output_directory, "/tpr_nbr_isoforms_above5.pdf"),
#                pointsize = 3, stripsize = 10,
#                axistitlesize = 15, axistextsize = 10) 
# 
# plot_perftable(path_to_results = path_to_results_drosophila, 
#                truth = truth_drosophila, 
#                methods = methods, method_names = method_names, 
#                method_cols = method_cols, threshold = 0.05, 
#                split_variable = "nbr_isoforms_a5",
#                output_filename = paste0(output_directory, "/perf_nbr_isoforms_above5_drosophila.pdf"))
# 
# plot_perftable(path_to_results = path_to_results_human, 
#                truth = truth_human, 
#                methods = methods, method_names = method_names, 
#                method_cols = method_cols, threshold = 0.05, 
#                split_variable = "nbr_isoforms_a5",
#                output_filename = paste0(output_directory, "/perf_nbr_isoforms_above5_human.pdf"))
# 
# ## Filtering
# plot_fdr_tpr_paper(path_to_results_drosophila = path_to_results_drosophila, 
#                    path_to_results_human = path_to_results_human,
#                    truth_drosophila = truth_drosophila, 
#                    truth_human = truth_human, 
#                    methods = methods_filtering, 
#                    method_names = method_names_filtering, 
#                    method_cols = method_cols_filtering, thresholds = thresholds, 
#                    split_variable = "nbr_isoforms_a5_2",
#                    output_filename = paste0(output_directory, 
#                                             "/fdr_tpr_nbr_isoforms_above5_filtering.pdf"),
#                    pointsize = 3, stripsize = 10,
#                    axistitlesize = 15, axistextsize = 10) 
# 
# plot_tpr_paper(path_to_results_drosophila = path_to_results_drosophila, 
#                path_to_results_human = path_to_results_human,
#                truth_drosophila = truth_drosophila, 
#                truth_human = truth_human, 
#                methods = methods_filtering, 
#                method_names = method_names_filtering, 
#                method_cols = method_cols_filtering, thresholds = thresholds, 
#                split_variable = "nbr_isoforms_a5_2",
#                output_filename = paste0(output_directory, 
#                                         "/tpr_nbr_isoforms_above5_filtering.pdf"),
#                pointsize = 3, stripsize = 10,
#                axistitlesize = 15, axistextsize = 10) 
# 
# plot_perftable(path_to_results = path_to_results_drosophila, 
#                truth = truth_drosophila, 
#                methods = methods_filtering, 
#                method_names = method_names_filtering, 
#                method_cols = method_cols, threshold = 0.05, 
#                split_variable = "nbr_isoforms_a5",
#                output_filename = paste0(output_directory, 
#                                         "/perf_nbr_isoforms_above5_filtering_drosophila.pdf"))
# 
# plot_perftable(path_to_results = path_to_results_human, 
#                truth = truth_human, 
#                methods = methods_filtering, 
#                method_names = method_names_filtering, 
#                method_cols = method_cols_filtering, threshold = 0.05, 
#                split_variable = "nbr_isoforms_a5",
#                output_filename = paste0(output_directory, 
#                                         "/perf_nbr_isoforms_above5_filtering_human.pdf"))
# 
# ## SPLIT BY NUMBER OF ISOFORMS WITH ABUNDANCE > 10%
# truth_drosophila$nbr_isoforms_a10 <- cut2(truth_drosophila$nbr_isoforms_above10, 
#                                          cuts = c(2, 3, max(truth_drosophila$nbr_isoforms_above10)),
#                                          oneval = FALSE)
# truth_human$nbr_isoforms_a10 <- cut2(truth_human$nbr_isoforms_above10, 
#                                     cuts = c(2, 3, max(truth_human$nbr_isoforms_above10)),
#                                     oneval = FALSE)
# truth_drosophila$nbr_isoforms_a10_2 <- sapply(truth_drosophila$nbr_isoforms_a10, function(i) {
#   paste0(i, " (n = ", length(which(truth_drosophila$nbr_isoforms_a10 == i)), ")")
# })
# truth_human$nbr_isoforms_a10_2 <- sapply(truth_human$nbr_isoforms_a10, function(i) {
#   paste0(i, " (n = ", length(which(truth_human$nbr_isoforms_a10 == i)), ")")
# })
# ## Original data
# plot_fdr_tpr_paper(path_to_results_drosophila = path_to_results_drosophila, 
#                    path_to_results_human = path_to_results_human,
#                    truth_drosophila = truth_drosophila, 
#                    truth_human = truth_human,  
#                    methods = methods, method_names = method_names, 
#                    method_cols = method_cols, thresholds = thresholds, 
#                    split_variable = "nbr_isoforms_a10_2",
#                    output_filename = paste0(output_directory, "/fdr_tpr_nbr_isoforms_above10.pdf"),
#                    pointsize = 3, stripsize = 10,
#                    axistitlesize = 15, axistextsize = 10) 
# 
# plot_tpr_paper(path_to_results_drosophila = path_to_results_drosophila, 
#                path_to_results_human = path_to_results_human,
#                truth_drosophila = truth_drosophila, 
#                truth_human = truth_human,  
#                methods = methods, method_names = method_names, 
#                method_cols = method_cols, thresholds = thresholds, 
#                split_variable = "nbr_isoforms_a10_2",
#                output_filename = paste0(output_directory, "/tpr_nbr_isoforms_above10.pdf"),
#                pointsize = 3, stripsize = 10,
#                axistitlesize = 15, axistextsize = 10) 
# 
# plot_perftable(path_to_results = path_to_results_drosophila, 
#                truth = truth_drosophila, 
#                methods = methods, method_names = method_names, 
#                method_cols = method_cols, threshold = 0.05, 
#                split_variable = "nbr_isoforms_a10",
#                output_filename = paste0(output_directory, "/perf_nbr_isoforms_above10_drosophila.pdf"))
# 
# plot_perftable(path_to_results = path_to_results_human, 
#                truth = truth_human, 
#                methods = methods, method_names = method_names, 
#                method_cols = method_cols, threshold = 0.05, 
#                split_variable = "nbr_isoforms_a10",
#                output_filename = paste0(output_directory, "/perf_nbr_isoforms_above10_human.pdf"))
# 
# ## Filtering
# plot_fdr_tpr_paper(path_to_results_drosophila = path_to_results_drosophila, 
#                    path_to_results_human = path_to_results_human,
#                    truth_drosophila = truth_drosophila, 
#                    truth_human = truth_human, 
#                    methods = methods_filtering, 
#                    method_names = method_names_filtering, 
#                    method_cols = method_cols_filtering, thresholds = thresholds, 
#                    split_variable = "nbr_isoforms_a10_2",
#                    output_filename = paste0(output_directory, 
#                                             "/fdr_tpr_nbr_isoforms_above10_filtering.pdf"),
#                    pointsize = 3, stripsize = 10,
#                    axistitlesize = 15, axistextsize = 10) 
# 
# plot_tpr_paper(path_to_results_drosophila = path_to_results_drosophila, 
#                path_to_results_human = path_to_results_human,
#                truth_drosophila = truth_drosophila, 
#                truth_human = truth_human, 
#                methods = methods_filtering, 
#                method_names = method_names_filtering, 
#                method_cols = method_cols_filtering, thresholds = thresholds, 
#                split_variable = "nbr_isoforms_a10_2",
#                output_filename = paste0(output_directory, 
#                                         "/tpr_nbr_isoforms_above10_filtering.pdf"),
#                pointsize = 3, stripsize = 10,
#                axistitlesize = 15, axistextsize = 10) 
# 
# plot_perftable(path_to_results = path_to_results_drosophila, 
#                truth = truth_drosophila, 
#                methods = methods_filtering, 
#                method_names = method_names_filtering, 
#                method_cols = method_cols, threshold = 0.05, 
#                split_variable = "nbr_isoforms_a10",
#                output_filename = paste0(output_directory, 
#                                         "/perf_nbr_isoforms_above10_filtering_drosophila.pdf"))
# 
# plot_perftable(path_to_results = path_to_results_human, 
#                truth = truth_human, 
#                methods = methods_filtering, 
#                method_names = method_names_filtering, 
#                method_cols = method_cols_filtering, threshold = 0.05, 
#                split_variable = "nbr_isoforms_a10",
#                output_filename = paste0(output_directory, 
#                                         "/perf_nbr_isoforms_above10_filtering_human.pdf"))
# 
# ## SPLIT BY NUMBER OF ISOFORMS WITH ABUNDANCE > 15%
# truth_drosophila$nbr_isoforms_a15 <- cut2(truth_drosophila$nbr_isoforms_above15, 
#                                           cuts = c(2, 3, max(truth_drosophila$nbr_isoforms_above15)),
#                                           oneval = FALSE)
# truth_human$nbr_isoforms_a15 <- cut2(truth_human$nbr_isoforms_above15, 
#                                      cuts = c(2, 3, max(truth_human$nbr_isoforms_above15)),
#                                      oneval = FALSE)
# truth_drosophila$nbr_isoforms_a15_2 <- sapply(truth_drosophila$nbr_isoforms_a15, function(i) {
#   paste0(i, " (n = ", length(which(truth_drosophila$nbr_isoforms_a15 == i)), ")")
# })
# truth_human$nbr_isoforms_a15_2 <- sapply(truth_human$nbr_isoforms_a15, function(i) {
#   paste0(i, " (n = ", length(which(truth_human$nbr_isoforms_a15 == i)), ")")
# })
# ## Original data
# plot_fdr_tpr_paper(path_to_results_drosophila = path_to_results_drosophila, 
#                    path_to_results_human = path_to_results_human,
#                    truth_drosophila = truth_drosophila, 
#                    truth_human = truth_human,  
#                    methods = methods, method_names = method_names, 
#                    method_cols = method_cols, thresholds = thresholds, 
#                    split_variable = "nbr_isoforms_a15_2",
#                    output_filename = paste0(output_directory, "/fdr_tpr_nbr_isoforms_above15.pdf"),
#                    pointsize = 3, stripsize = 10,
#                    axistitlesize = 15, axistextsize = 10) 
# 
# plot_tpr_paper(path_to_results_drosophila = path_to_results_drosophila, 
#                path_to_results_human = path_to_results_human,
#                truth_drosophila = truth_drosophila, 
#                truth_human = truth_human,  
#                methods = methods, method_names = method_names, 
#                method_cols = method_cols, thresholds = thresholds, 
#                split_variable = "nbr_isoforms_a15_2",
#                output_filename = paste0(output_directory, "/tpr_nbr_isoforms_above15.pdf"),
#                pointsize = 3, stripsize = 10,
#                axistitlesize = 15, axistextsize = 10) 
# 
# plot_perftable(path_to_results = path_to_results_drosophila, 
#                truth = truth_drosophila, 
#                methods = methods, method_names = method_names, 
#                method_cols = method_cols, threshold = 0.05, 
#                split_variable = "nbr_isoforms_a15",
#                output_filename = paste0(output_directory, "/perf_nbr_isoforms_above15_drosophila.pdf"))
# 
# plot_perftable(path_to_results = path_to_results_human, 
#                truth = truth_human, 
#                methods = methods, method_names = method_names, 
#                method_cols = method_cols, threshold = 0.05, 
#                split_variable = "nbr_isoforms_a15",
#                output_filename = paste0(output_directory, "/perf_nbr_isoforms_above15_human.pdf"))
# 
# ## Filtering
# plot_fdr_tpr_paper(path_to_results_drosophila = path_to_results_drosophila, 
#                    path_to_results_human = path_to_results_human,
#                    truth_drosophila = truth_drosophila, 
#                    truth_human = truth_human, 
#                    methods = methods_filtering, 
#                    method_names = method_names_filtering, 
#                    method_cols = method_cols_filtering, thresholds = thresholds, 
#                    split_variable = "nbr_isoforms_a15_2",
#                    output_filename = paste0(output_directory, 
#                                             "/fdr_tpr_nbr_isoforms_above15_filtering.pdf"),
#                    pointsize = 3, stripsize = 10,
#                    axistitlesize = 15, axistextsize = 10) 
# 
# plot_tpr_paper(path_to_results_drosophila = path_to_results_drosophila, 
#                path_to_results_human = path_to_results_human,
#                truth_drosophila = truth_drosophila, 
#                truth_human = truth_human, 
#                methods = methods_filtering, 
#                method_names = method_names_filtering, 
#                method_cols = method_cols_filtering, thresholds = thresholds, 
#                split_variable = "nbr_isoforms_a15_2",
#                output_filename = paste0(output_directory, 
#                                         "/tpr_nbr_isoforms_above15_filtering.pdf"),
#                pointsize = 3, stripsize = 10,
#                axistitlesize = 15, axistextsize = 10) 
# 
# plot_perftable(path_to_results = path_to_results_drosophila, 
#                truth = truth_drosophila, 
#                methods = methods_filtering, 
#                method_names = method_names_filtering, 
#                method_cols = method_cols, threshold = 0.05, 
#                split_variable = "nbr_isoforms_a15",
#                output_filename = paste0(output_directory, 
#                                         "/perf_nbr_isoforms_above15_filtering_drosophila.pdf"))
# 
# plot_perftable(path_to_results = path_to_results_human, 
#                truth = truth_human, 
#                methods = methods_filtering, 
#                method_names = method_names_filtering, 
#                method_cols = method_cols_filtering, threshold = 0.05, 
#                split_variable = "nbr_isoforms_a15",
#                output_filename = paste0(output_directory, 
#                                         "/perf_nbr_isoforms_above15_filtering_human.pdf"))
# 
# ## SPLIT BY NUMBER OF ISOFORMS WITH ABUNDANCE > 25%
# truth_drosophila$nbr_isoforms_a25 <- cut2(truth_drosophila$nbr_isoforms_above25, 
#                                           cuts = c(2, 3, max(truth_drosophila$nbr_isoforms_above25)),
#                                           oneval = FALSE)
# truth_human$nbr_isoforms_a25 <- cut2(truth_human$nbr_isoforms_above25, 
#                                      cuts = c(2, 3, max(truth_human$nbr_isoforms_above25)),
#                                      oneval = FALSE)
# truth_drosophila$nbr_isoforms_a25_2 <- sapply(truth_drosophila$nbr_isoforms_a25, function(i) {
#   paste0(i, " (n = ", length(which(truth_drosophila$nbr_isoforms_a25 == i)), ")")
# })
# truth_human$nbr_isoforms_a25_2 <- sapply(truth_human$nbr_isoforms_a25, function(i) {
#   paste0(i, " (n = ", length(which(truth_human$nbr_isoforms_a25 == i)), ")")
# })
# ## Original data
# plot_fdr_tpr_paper(path_to_results_drosophila = path_to_results_drosophila, 
#                    path_to_results_human = path_to_results_human,
#                    truth_drosophila = truth_drosophila, 
#                    truth_human = truth_human,  
#                    methods = methods, method_names = method_names, 
#                    method_cols = method_cols, thresholds = thresholds, 
#                    split_variable = "nbr_isoforms_a25_2",
#                    output_filename = paste0(output_directory, "/fdr_tpr_nbr_isoforms_above25.pdf"),
#                    pointsize = 3, stripsize = 10,
#                    axistitlesize = 15, axistextsize = 10) 
# 
# plot_tpr_paper(path_to_results_drosophila = path_to_results_drosophila, 
#                path_to_results_human = path_to_results_human,
#                truth_drosophila = truth_drosophila, 
#                truth_human = truth_human,  
#                methods = methods, method_names = method_names, 
#                method_cols = method_cols, thresholds = thresholds, 
#                split_variable = "nbr_isoforms_a25_2",
#                output_filename = paste0(output_directory, "/tpr_nbr_isoforms_above25.pdf"),
#                pointsize = 3, stripsize = 10,
#                axistitlesize = 15, axistextsize = 10) 
# 
# plot_perftable(path_to_results = path_to_results_drosophila, 
#                truth = truth_drosophila, 
#                methods = methods, method_names = method_names, 
#                method_cols = method_cols, threshold = 0.05, 
#                split_variable = "nbr_isoforms_a25",
#                output_filename = paste0(output_directory, "/perf_nbr_isoforms_above25_drosophila.pdf"))
# 
# plot_perftable(path_to_results = path_to_results_human, 
#                truth = truth_human, 
#                methods = methods, method_names = method_names, 
#                method_cols = method_cols, threshold = 0.05, 
#                split_variable = "nbr_isoforms_a25",
#                output_filename = paste0(output_directory, "/perf_nbr_isoforms_above25_human.pdf"))
# 
# ## Filtering
# plot_fdr_tpr_paper(path_to_results_drosophila = path_to_results_drosophila, 
#                    path_to_results_human = path_to_results_human,
#                    truth_drosophila = truth_drosophila, 
#                    truth_human = truth_human, 
#                    methods = methods_filtering, 
#                    method_names = method_names_filtering, 
#                    method_cols = method_cols_filtering, thresholds = thresholds, 
#                    split_variable = "nbr_isoforms_a25_2",
#                    output_filename = paste0(output_directory, 
#                                             "/fdr_tpr_nbr_isoforms_above25_filtering.pdf"),
#                    pointsize = 3, stripsize = 10,
#                    axistitlesize = 15, axistextsize = 10) 
# 
# plot_tpr_paper(path_to_results_drosophila = path_to_results_drosophila, 
#                path_to_results_human = path_to_results_human,
#                truth_drosophila = truth_drosophila, 
#                truth_human = truth_human, 
#                methods = methods_filtering, 
#                method_names = method_names_filtering, 
#                method_cols = method_cols_filtering, thresholds = thresholds, 
#                split_variable = "nbr_isoforms_a25_2",
#                output_filename = paste0(output_directory, 
#                                         "/tpr_nbr_isoforms_above25_filtering.pdf"),
#                pointsize = 3, stripsize = 10,
#                axistitlesize = 15, axistextsize = 10) 
# 
# plot_perftable(path_to_results = path_to_results_drosophila, 
#                truth = truth_drosophila, 
#                methods = methods_filtering, 
#                method_names = method_names_filtering, 
#                method_cols = method_cols, threshold = 0.05, 
#                split_variable = "nbr_isoforms_a25",
#                output_filename = paste0(output_directory, 
#                                         "/perf_nbr_isoforms_above25_filtering_drosophila.pdf"))
# 
# plot_perftable(path_to_results = path_to_results_human, 
#                truth = truth_human, 
#                methods = methods_filtering, 
#                method_names = method_names_filtering, 
#                method_cols = method_cols_filtering, threshold = 0.05, 
#                split_variable = "nbr_isoforms_a25",
#                output_filename = paste0(output_directory, 
#                                         "/perf_nbr_isoforms_above25_filtering_human.pdf"))

## SPLIT BY NUMBER OF BINS (DEXSeq-noaggreg)
plot_counting_characteristics(path_to_results_drosophila = path_to_results_drosophila,
                              path_to_results_human = path_to_results_human, 
                              truth_drosophila = truth_drosophila,
                              truth_human = truth_human, 
                              showmeasure = "dispvar", 
                              labelmeasure = "within-group variance(log2(dispersion))", 
                              method = "dexseq_htseq_nomerge", 
                              method_name = "DEXSeq-noaggreg", 
                              threshold = 0.05, split_variable = "nbrbinbin_2",
                              output_filename = paste0(output_directory, "/counting_nbrbins_dexseq_noaggreg.pdf"), 
                              stripsize = 10,
                              axistitlesize = 15, axistextsize = 10,
                              legend_nrow = 2, linewidth = 1.5)

## SPLIT BY NUMBER OF ISOFORMS
truth_drosophila$nbr_isoforms_cat <- cut2(truth_drosophila$nbr_isoforms, 
                                          cuts = c(2, 4, 10, max(truth_drosophila$nbr_isoforms)))
truth_human$nbr_isoforms_cat <- cut2(truth_human$nbr_isoforms, 
                                     cuts = c(2, 4, 10, max(truth_human$nbr_isoforms)))
truth_drosophila$nbr_isoforms_cat2 <- sapply(truth_drosophila$nbr_isoforms_cat, function(i) {
  paste0(i, " (n = ", length(which(truth_drosophila$nbr_isoforms_cat == i)), ")")
})
truth_human$nbr_isoforms_cat2 <- sapply(truth_human$nbr_isoforms_cat, function(i) {
  paste0(i, " (n = ", length(which(truth_human$nbr_isoforms_cat == i)), ")")
})

truth_drosophila$nbr_isoforms_cat3 <- sapply(truth_drosophila$nbr_isoforms_cat, function(i) {
  paste0(i, " (n = ", length(which(truth_drosophila$nbr_isoforms_cat == i)), ", n.ds = ",
         length(intersect(which(truth_drosophila$nbr_isoforms_cat == i), 
                          which(truth_drosophila$ds_status == 1))), ")")
})
truth_human$nbr_isoforms_cat3 <- sapply(truth_human$nbr_isoforms_cat, function(i) {
  paste0(i, " (n = ", length(which(truth_human$nbr_isoforms_cat == i)), ", n.ds = ",
         length(intersect(which(truth_human$nbr_isoforms_cat == i), 
                          which(truth_human$ds_status == 1))), ")")
})

## Original data
plot_fdr_tpr_paper(path_to_results_drosophila = path_to_results_drosophila, 
                   path_to_results_human = path_to_results_human,
                   truth_drosophila = droplevels(subset(truth_drosophila, nbr_isoforms_cat != " 1")), 
                   truth_human = droplevels(subset(truth_human, nbr_isoforms_cat != " 1")),  
                   methods = methods, method_names = method_names, 
                   method_cols = method_cols, thresholds = thresholds, 
                   split_variable = "nbr_isoforms_cat2",
                   output_filename = paste0(output_directory, "/fdr_tpr_nbr_isoforms.pdf"),
                   pointsize = 3, stripsize = 10,
                   axistitlesize = 15, axistextsize = 10) 

plot_fdr_tpr_paper(path_to_results_drosophila = path_to_results_drosophila, 
                   path_to_results_human = path_to_results_human,
                   truth_drosophila = droplevels(subset(truth_drosophila, nbr_isoforms_cat != " 1")), 
                   truth_human = droplevels(subset(truth_human, nbr_isoforms_cat != " 1")),  
                   methods = methods, method_names = method_names, 
                   method_cols = method_cols, thresholds = thresholds, 
                   split_variable = "nbr_isoforms_cat3",
                   output_filename = paste0(output_directory, "/fdr_tpr_nbr_isoforms_withnds.pdf"),
                   pointsize = 3, stripsize = 10,
                   axistitlesize = 15, axistextsize = 10) 

plot_tpr_paper(path_to_results_drosophila = path_to_results_drosophila, 
               path_to_results_human = path_to_results_human,
               truth_drosophila = droplevels(subset(truth_drosophila, nbr_isoforms_cat != " 1")), 
               truth_human = droplevels(subset(truth_human, nbr_isoforms_cat != " 1")),  
               methods = methods, method_names = method_names, 
               method_cols = method_cols, thresholds = thresholds, 
               split_variable = "nbr_isoforms_cat2",
               output_filename = paste0(output_directory, "/tpr_nbr_isoforms.pdf"),
               pointsize = 3, stripsize = 10,
               axistitlesize = 15, axistextsize = 10) 

plot_perftable(path_to_results = path_to_results_drosophila, 
               truth = truth_drosophila, 
               methods = methods, method_names = method_names, 
               method_cols = method_cols, threshold = 0.05, 
               split_variable = "nbr_isoforms_cat",
               output_filename = paste0(output_directory, "/perf_nbr_isoforms_drosophila.pdf"))

plot_perftable(path_to_results = path_to_results_human, 
               truth = truth_human, 
               methods = methods, method_names = method_names, 
               method_cols = method_cols, threshold = 0.05, 
               split_variable = "nbr_isoforms_cat",
               output_filename = paste0(output_directory, "/perf_nbr_isoforms_human.pdf"))

## Filtering
plot_fdr_tpr_paper(path_to_results_drosophila = path_to_results_drosophila, 
                   path_to_results_human = path_to_results_human,
                   truth_drosophila = droplevels(subset(truth_drosophila, nbr_isoforms_cat != " 1")), 
                   truth_human = droplevels(subset(truth_human, nbr_isoforms_cat != " 1")), 
                   methods = methods_filtering, 
                   method_names = method_names_filtering, 
                   method_cols = method_cols_filtering, thresholds = thresholds, 
                   split_variable = "nbr_isoforms_cat2",
                   output_filename = paste0(output_directory, 
                                            "/fdr_tpr_nbr_isoforms_filtering.pdf"),
                   pointsize = 3, stripsize = 10,
                   axistitlesize = 15, axistextsize = 10) 

plot_tpr_paper(path_to_results_drosophila = path_to_results_drosophila, 
               path_to_results_human = path_to_results_human,
               truth_drosophila = droplevels(subset(truth_drosophila, nbr_isoforms_cat != " 1")), 
               truth_human = droplevels(subset(truth_human, nbr_isoforms_cat != " 1")), 
               methods = methods_filtering, 
               method_names = method_names_filtering, 
               method_cols = method_cols_filtering, thresholds = thresholds, 
               split_variable = "nbr_isoforms_cat2",
               output_filename = paste0(output_directory, 
                                        "/tpr_nbr_isoforms_filtering.pdf"),
               pointsize = 3, stripsize = 10,
               axistitlesize = 15, axistextsize = 10) 

plot_perftable(path_to_results = path_to_results_drosophila, 
               truth = truth_drosophila, 
               methods = methods_filtering, 
               method_names = method_names_filtering, 
               method_cols = method_cols, threshold = 0.05, 
               split_variable = "nbr_isoforms_cat",
               output_filename = paste0(output_directory, 
                                        "/perf_nbr_isoforms_filtering_drosophila.pdf"))

plot_perftable(path_to_results = path_to_results_human, 
               truth = truth_human, 
               methods = methods_filtering, 
               method_names = method_names_filtering, 
               method_cols = method_cols_filtering, threshold = 0.05, 
               split_variable = "nbr_isoforms_cat",
               output_filename = paste0(output_directory, 
                                        "/perf_nbr_isoforms_filtering_human.pdf"))

## SPLIT BY AS TYPE
truth_drosophila$as_type_2 <- truth_drosophila$as_type
truth_drosophila$as_type_2[!truth_drosophila$as_type_2 %in% 
                             c(NA, "0,1-2^", "1^,2^", "1-,2-", "0,1^2-", "1-2^,3-4^")] <- "other"
truth_drosophila$as_type_2[truth_drosophila$as_type_2 == "0,1-2^"] <- "skipped exon"
truth_drosophila$as_type_2[truth_drosophila$as_type_2 == "1^,2^"] <- "alternative donors"
truth_drosophila$as_type_2[truth_drosophila$as_type_2 == "1-,2-"] <- "alternative acceptors"
truth_drosophila$as_type_2[truth_drosophila$as_type_2 == "0,1^2-"] <- "retained intron"
truth_drosophila$as_type_2[truth_drosophila$as_type_2 == "1-2^,3-4^"] <- "mutually exclusive exons"

truth_human$as_type_2 <- truth_human$as_type
truth_human$as_type_2[!truth_human$as_type_2 %in% 
                             c(NA, "0,1-2^", "1^,2^", "1-,2-", "0,1^2-", "1-2^,3-4^")] <- "other"
truth_human$as_type_2[truth_human$as_type_2 == "0,1-2^"] <- "skipped exon"
truth_human$as_type_2[truth_human$as_type_2 == "1^,2^"] <- "alternative donors"
truth_human$as_type_2[truth_human$as_type_2 == "1-,2-"] <- "alternative acceptors"
truth_human$as_type_2[truth_human$as_type_2 == "0,1^2-"] <- "retained intron"
truth_human$as_type_2[truth_human$as_type_2 == "1-2^,3-4^"] <- "mutually exclusive exons"

truth_drosophila$as_type3 <- sapply(truth_drosophila$as_type_2, function(i) {
  paste0(i, " (n = ", length(which(truth_drosophila$as_type_2 == i)), ")")
})
truth_human$as_type3 <- sapply(truth_human$as_type_2, function(i) {
  paste0(i, " (n = ", length(which(truth_human$as_type_2 == i)), ")")
})

plot_tpr_paper(path_to_results_drosophila = path_to_results_drosophila, 
               path_to_results_human = NULL,
               truth_drosophila = droplevels(subset(truth_drosophila, !is.na(as_type_2))), 
               truth_human = NULL,  
               methods = methods, method_names = method_names, 
               method_cols = method_cols, thresholds = thresholds, 
               split_variable = "as_type3",
               output_filename = paste0(output_directory, "/tpr_as_type.pdf"),
               pointsize = 3, stripsize = 8,
               axistitlesize = 15, axistextsize = 10) 

plot_perftable(path_to_results = path_to_results_drosophila, 
               truth = truth_drosophila, 
               methods = methods, method_names = method_names, 
               method_cols = method_cols, threshold = 0.05, 
               split_variable = "as_type_2",
               output_filename = paste0(output_directory, "/perf_as_type_drosophila.pdf"))

plot_perftable(path_to_results = path_to_results_human, 
               truth = truth_human, 
               methods = methods, method_names = method_names, 
               method_cols = method_cols, threshold = 0.05, 
               split_variable = "as_type_2",
               output_filename = paste0(output_directory, "/perf_as_type_human.pdf"))

## SPLIT BY NUMBER OF DIFFERING BP BETWEEN DIFF USED ISOFORMS
truth_drosophila$diffbp_3 <- cut2(truth_drosophila$diffbp, 
                                  cuts = c(100, 1000, max(truth_drosophila$diffbp, 
                                                          na.rm = TRUE)))
truth_human$diffbp_3 <- cut2(truth_human$diffbp, 
                             cuts = c(100, 1000, max(truth_human$diffbp, 
                                                     na.rm = TRUE)))
truth_drosophila$diffbp4 <- sapply(truth_drosophila$diffbp_3, function(i) {
  paste0(i, " (n = ", length(which(truth_drosophila$diffbp_3 == i)), ")")
})
truth_human$diffbp4 <- sapply(truth_human$diffbp_3, function(i) {
  paste0(i, " (n = ", length(which(truth_human$diffbp_3 == i)), ")")
})

## Original data
plot_tpr_paper(path_to_results_drosophila = path_to_results_drosophila, 
               path_to_results_human = path_to_results_human,
               truth_drosophila = droplevels(subset(truth_drosophila, !is.na(diffbp_3))), 
               truth_human = droplevels(subset(truth_human, !is.na(diffbp_3))),  
               methods = methods, method_names = method_names, 
               method_cols = method_cols, thresholds = thresholds, 
               split_variable = "diffbp4",
               output_filename = paste0(output_directory, "/tpr_diffbp.pdf"),
               pointsize = 3, stripsize = 10,
               axistitlesize = 15, axistextsize = 10) 

plot_perftable(path_to_results = path_to_results_drosophila, 
               truth = truth_drosophila, 
               methods = methods, method_names = method_names, 
               method_cols = method_cols, threshold = 0.05, 
               split_variable = "diffbp_3",
               output_filename = paste0(output_directory, "/perf_diffbp_drosophila.pdf"))

plot_perftable(path_to_results = path_to_results_human, 
               truth = truth_human, 
               methods = methods, method_names = method_names, 
               method_cols = method_cols, threshold = 0.05, 
               split_variable = "diffbp_3",
               output_filename = paste0(output_directory, "/perf_diffbp_human.pdf"))

## Filtering
plot_tpr_paper(path_to_results_drosophila = path_to_results_drosophila, 
               path_to_results_human = path_to_results_human,
               truth_drosophila = droplevels(subset(truth_drosophila, !is.na(diffbp_3))), 
               truth_human = droplevels(subset(truth_human, !is.na(diffbp_3))), 
               methods = methods_filtering, 
               method_names = method_names_filtering, 
               method_cols = method_cols_filtering, thresholds = thresholds, 
               split_variable = "diffbp4",
               output_filename = paste0(output_directory, 
                                        "/tpr_diffbp_filtering.pdf"),
               pointsize = 3, stripsize = 10,
               axistitlesize = 15, axistextsize = 10) 

plot_perftable(path_to_results = path_to_results_drosophila, 
               truth = truth_drosophila, 
               methods = methods_filtering, 
               method_names = method_names_filtering, 
               method_cols = method_cols, threshold = 0.05, 
               split_variable = "diffbp_3",
               output_filename = paste0(output_directory, 
                                        "/perf_diffbp_filtering_drosophila.pdf"))

plot_perftable(path_to_results = path_to_results_human, 
               truth = truth_human, 
               methods = methods_filtering, 
               method_names = method_names_filtering, 
               method_cols = method_cols_filtering, threshold = 0.05, 
               split_variable = "diffbp_3",
               output_filename = paste0(output_directory, 
                                        "/perf_diffbp_filtering_human.pdf"))

## SPLIT BY DIFFERENTIAL EXPRESSION STATUS
# truth_human_de <- read.delim(paste0(basedir, "/hsapiens/with_diffexpression/", 
#                                     "non_null_simulation/3_truth/", 
#                                     "truth_human_non_null.txt"),
#                              header = TRUE, as.is = TRUE)
# path_to_results_human_de <- paste0(basedir, "/hsapiens/with_diffexpression/", 
#                                    "non_null_simulation/4_results")
# 
# truth_human_de$de_status2 <- truth_human_de$de_status
# truth_human_de$de_status2[truth_human_de$de_status2 == 0] <- "not DE"
# truth_human_de$de_status2[truth_human_de$de_status2 == 1] <- "DE"
# truth_human_de$de_status3 <- sapply(truth_human_de$de_status2, function(i) {
#   paste0(i, " (n = ", length(which(truth_human_de$de_status2 == i)), ")")
# })
# 
# plot_fdr_tpr_paper(path_to_results_drosophila = NULL, 
#                    path_to_results_human = path_to_results_human_de,
#                    truth_drosophila = NULL, 
#                    truth_human = truth_human_de,  
#                    methods = methods, method_names = method_names, 
#                    method_cols = method_cols, thresholds = thresholds, 
#                    split_variable = "de_status3",
#                    output_filename = paste0(output_directory, "/fdr_tpr_diffexpstatus.pdf"),
#                    pointsize = 3, stripsize = 15,
#                    axistitlesize = 15, axistextsize = 10) 
# 
# plot_tpr_paper(path_to_results_drosophila = NULL, 
#                path_to_results_human = path_to_results_human_de,
#                truth_drosophila = NULL, 
#                truth_human = truth_human_de,  
#                methods = methods, method_names = method_names, 
#                method_cols = method_cols, thresholds = thresholds, 
#                split_variable = "de_status3",
#                output_filename = paste0(output_directory, "/tpr_diffexpstatus.pdf"),
#                pointsize = 3, stripsize = 15,
#                axistitlesize = 15, axistextsize = 10) 
# 
# plot_perftable(path_to_results = path_to_results_human, 
#                truth = truth_human_de, 
#                methods = methods, method_names = method_names, 
#                method_cols = method_cols, threshold = 0.05, 
#                split_variable = "de_status2",
#                output_filename = paste0(output_directory, "/perf_diffexpstatus_human.pdf"))

## GENERATE TRUTH FILES WITH ANNOTATIONS USED ABOVE
truth_human_save <- truth_human[, c("gene", "ds_status", "nbr_missing_ds_tr", 
                                    "nbr_missing_combined", "nbr_remaining_2", 
                                    "nbr_missing_combined_2", "diff_IsoPct_3", 
                                    "exprlevel", "nbr_isoforms_cat", 
                                    "as_type_2", "diffbp_3")]
colnames(truth_human_save) <- c("gene", "ds_status", "nbr_missing_ds_transcripts", 
                                "nbr_missing_ds_nonds_transcripts",
                                "nbr_remaining_transcripts", 
                                "nbr_missing_ds_nbr_remaining_transcripts",
                                "diff_rel_abundance", "gene_expression",
                                "nbr_isoforms", "AS_event_type", 
                                "differentiating_basepairs")

truth_drosophila_save <- truth_drosophila[, c("gene", "ds_status", "nbr_missing_ds_tr", 
                                              "nbr_missing_combined", "nbr_remaining_2", 
                                              "nbr_missing_combined_2", "diff_IsoPct_3", 
                                              "exprlevel", "nbr_isoforms_cat", 
                                              "as_type_2", "diffbp_3")]
colnames(truth_drosophila_save) <- c("gene", "ds_status", "nbr_missing_ds_transcripts", 
                                     "nbr_missing_ds_nonds_transcripts",
                                     "nbr_remaining_transcripts", 
                                     "nbr_missing_ds_nbr_remaining_transcripts",
                                     "diff_rel_abundance", "gene_expression",
                                     "nbr_isoforms", "AS_event_type", 
                                     "differentiating_basepairs")

write.table(truth_human_save, file = paste0(output_directory_truth, "/", 
                                            gsub("\\.txt", "_ms.txt", 
                                                 basename(path_to_truth_human))),
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
write.table(truth_drosophila_save, 
            file = paste0(output_directory_truth, "/", 
                          gsub("\\.txt", "_ms.txt",
                               basename(path_to_truth_drosophila))),
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
