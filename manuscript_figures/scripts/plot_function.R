calc_FDR <- function(methods, method_names, method_cols, organism, 
                     thresholds, truth, path_to_results, method_pch = NULL,
                     only_shared = FALSE) {
  out <- as.data.frame(do.call(rbind, lapply(methods, function(w) {
    tmp <- read.delim(paste0(path_to_results, "/", w, ".txt"), 
                      header = TRUE, as.is = TRUE, check.names = FALSE)
    if (isTRUE(only_shared)) {
      ist <- intersect(truth$gene, tmp$gene)
      tmp <- subset(tmp, gene %in% ist)
      truth <- subset(truth, gene %in% ist)
    } else {
      tmp <- subset(tmp, gene %in% truth$gene)
    }
    sapply(thresholds, function(i) {
      length(intersect(truth$gene[which(truth$ds_status == 0)],
                       tmp[which(tmp[, 2] <= i), 1]))/
        length(which(tmp[, 2] <= i))
    })
  })))
  colnames(out) <- thresholds
  rownames(out) <- methods
  out$method_name <- factor(method_names, levels = method_names)
  out$method_col <- method_cols
##  out$method_col <- factor(method_cols, levels = method_cols[!duplicated(method_cols)])
  if (!is.null(method_pch)) out$method_pch <- method_pch
  out$organism <- organism
  out
}

calc_TPR <- function(methods, method_names, method_cols, organism, 
                     thresholds, truth, path_to_results, method_pch = NULL,
                     only_shared = FALSE) {
  out <- as.data.frame(do.call(rbind, lapply(methods, function(w) {
    tmp <- read.delim(paste0(path_to_results, "/", w, ".txt"), 
                      header = TRUE, as.is = TRUE, check.names = FALSE)
    if (isTRUE(only_shared)) {
      ist <- intersect(truth$gene, tmp$gene)
      tmp <- subset(tmp, gene %in% ist)
      truth <- subset(truth, gene %in% ist)
    } else {
      tmp <- subset(tmp, gene %in% truth$gene)
    }
    sapply(thresholds, function(i) {
      length(intersect(truth$gene[which(truth$ds_status == 1)],
                       tmp[which(tmp[, 2] <= i), 1]))/
        length(which(truth$ds_status == 1))
    })
  })))
  colnames(out) <- thresholds
  rownames(out) <- methods
  out$method_name <- factor(method_names, levels = method_names)
  out$method_col <- method_cols
##  out$method_col <- factor(method_cols, levels = method_cols[!duplicated(method_cols)])
  if (!is.null(method_pch)) out$method_pch <- method_pch
  out$organism <- organism
  out
}

plot_counting_characteristics <- function(path_to_results_drosophila,
                                          path_to_results_human, 
                                          truth_drosophila,
                                          truth_human, 
                                          showmeasure = "dispvar", 
                                          labelmeasure = "within-gene variance(log2(dispersion))", 
                                          method, method_name, 
                                          threshold, split_variable,
                                          output_filename, stripsize = 20,
                                          axistitlesize = 20, axistextsize = 15,
                                          legend_nrow = 2, linewidth = 2) {
  ## Read result file
  library(DEXSeq)
  library(dplyr)
  library(Hmisc)
  ## Drosophila
  if (!is.null(path_to_results_drosophila) && !is.null(truth_drosophila)) {
    load(paste0(path_to_results_drosophila, "/", method, ".Rdata"))
    res$genomicData <- NULL
    res$transcripts <- NULL
    G_dr <- as.data.frame(res) %>% group_by(groupID) %>%
      summarise(dispvar = var(log2(dispersion), na.rm = TRUE),
                exprvar = var(log2(exonBaseMean + 1), na.rm = TRUE),
                nbrbin = length(dispersion),
                maxexpr = max(log2(exonBaseMean + 1), na.rm = TRUE),
                minexpr = min(log2(exonBaseMean + 1), na.rm = TRUE),
                rangexpr = maxexpr - minexpr)
    pgq <- data.frame(pgq)
    G_dr <- merge(G_dr, pgq, by.x = "groupID", by.y = 0, all = TRUE)
    G_dr$pgq[is.na(G_dr$pgq)] <- 1
    G_dr$nbrbinbin <- cut2(G_dr$nbrbin, cuts = c(10, 40))
    G_dr$nbrbinbin_2 <- sapply(G_dr$nbrbinbin, function(i) {
      paste0(i, " (n = ", length(which(G_dr$nbrbinbin == i)), ")")
    })
    G_dr$nbrbinbin_2[G_dr$nbrbinbin_2 == "NA (n = 0)"] <- NA
    G_dr <- merge(G_dr, truth_drosophila, by.x = "groupID", by.y = "gene", all = TRUE)
    G_dr$status <- NA
    G_dr$status[intersect(which(G_dr$ds_status == 1), which(G_dr$pgq <= threshold))] <- "TP"
    G_dr$status[intersect(which(G_dr$ds_status == 0), which(G_dr$pgq <= threshold))] <- "FP"
    G_dr$status[intersect(which(G_dr$ds_status == 1), which(G_dr$pgq > threshold))] <- "FN"
    G_dr$status[intersect(which(G_dr$ds_status == 0), which(G_dr$pgq > threshold))] <- "TN"
    G_dr$organism <- "Drosophila"
  } else {
    G_dr <- NULL
  }
  ## Human
  if (!is.null(path_to_results_human) && !is.null(truth_human)) {
    load(paste0(path_to_results_human, "/", method, ".Rdata"))
    res$genomicData <- NULL
    res$transcripts <- NULL
    G_hs <- as.data.frame(res) %>% group_by(groupID) %>%
      summarise(dispvar = var(log2(dispersion), na.rm = TRUE),
                exprvar = var(log2(exonBaseMean + 1), na.rm = TRUE),
                nbrbin = length(dispersion),
                maxexpr = max(log2(exonBaseMean + 1), na.rm = TRUE),
                minexpr = min(log2(exonBaseMean + 1), na.rm = TRUE),
                rangexpr = maxexpr - minexpr)
    pgq <- data.frame(pgq)
    G_hs <- merge(G_hs, pgq, by.x = "groupID", by.y = 0, all = TRUE)
    G_hs$pgq[is.na(G_hs$pgq)] <- 1
    G_hs$nbrbinbin <- cut2(G_hs$nbrbin, cuts = c(10, 40))
    G_hs$nbrbinbin_2 <- sapply(G_hs$nbrbinbin, function(i) {
      paste0(i, " (n = ", length(which(G_hs$nbrbinbin == i)), ")")
    })
    G_hs$nbrbinbin_2[G_hs$nbrbinbin_2 == "NA (n = 0)"] <- NA
    G_hs <- merge(G_hs, truth_human, by.x = "groupID", by.y = "gene", all = TRUE)
    G_hs$status <- NA
    G_hs$status[intersect(which(G_hs$ds_status == 1), which(G_hs$pgq <= threshold))] <- "TP"
    G_hs$status[intersect(which(G_hs$ds_status == 0), which(G_hs$pgq <= threshold))] <- "FP"
    G_hs$status[intersect(which(G_hs$ds_status == 1), which(G_hs$pgq > threshold))] <- "FN"
    G_hs$status[intersect(which(G_hs$ds_status == 0), which(G_hs$pgq > threshold))] <- "TN"
    G_hs$organism <- "Human"
  }  else{
    G_hs <- NULL
  }
  G <- rbind(G_dr, G_hs)
  if (split_variable == "nbrbinbin_2") {
    G <- subset(G, !is.na(G$nbrbinbin_2))
  }
  pdf(output_filename, width = 10)
  print(ggplot(G, aes_string(x = showmeasure, group = "status", col = "status")) + 
          geom_line(stat = "density", size = linewidth) + 
          scale_color_discrete(name = "") +
          facet_wrap(as.formula(paste("organism ~ ", split_variable)), scales = "free") + 
          xlab(labelmeasure) + 
          theme(legend.position = "bottom", 
                panel.background = element_rect(fill = NA, colour = "black"),
                panel.border = element_rect(fill = NA, colour = "black", size = 1),
                panel.grid.minor.x = element_blank(),
                panel.grid.minor.y = element_blank(),
                axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, 
                                           size = axistextsize),
                axis.text.y = element_text(size = axistextsize),
                axis.title.x = element_text(size = axistitlesize),
                axis.title.y = element_text(size = axistitlesize),
                strip.text = element_text(size = stripsize),
                strip.background = element_rect(fill = NA, colour = "black")) + 
          ggtitle(method_name))
  dev.off()
}
  

plot_fdr_tpr_paper <- function(path_to_results_drosophila, 
                               path_to_results_human,
                               truth_drosophila, truth_human, 
                               methods = NULL, method_names, 
                               methods_drosophila = NULL,
                               methods_human = NULL, 
                               method_cols, thresholds, split_variable,
                               output_filename,
                               pointsize = 3, stripsize = 20,
                               axistitlesize = 20, axistextsize = 15,
                               legend_nrow = 2, fig_height = 7, 
                               fig_width = 10, only_shared = FALSE) {
  
  if (is.null(methods_drosophila)) {
    methods_drosophila <- methods
  }
  if (is.null(methods_human)) {
    methods_human <- methods
  }
  ## Calculate FDR 
  if (is.null(split_variable)) {
    if (!is.null(path_to_results_drosophila) && (!is.null(truth_drosophila))) {
      FDR_drosophila <- calc_FDR(methods = methods_drosophila, method_names = method_names, 
                                 method_cols = method_cols, organism = "Drosophila", 
                                 thresholds = thresholds, truth = truth_drosophila,
                                 path_to_results = path_to_results_drosophila, 
                                 only_shared = only_shared)
      TPR_drosophila <- calc_TPR(methods = methods_drosophila, method_names = method_names, 
                                 method_cols = method_cols, organism = "Drosophila", 
                                 thresholds = thresholds, truth = truth_drosophila,
                                 path_to_results = path_to_results_drosophila,
                                 only_shared = only_shared)
    } else {
      FDR_drosophila <- NULL
      TPR_drosophila <- NULL
    }
    if (!is.null(path_to_results_human) && (!is.null(truth_human))) {
      FDR_human <- calc_FDR(methods = methods_human, method_names = method_names, 
                            method_cols = method_cols, organism = "Human", 
                            thresholds = thresholds, truth = truth_human,
                            path_to_results = path_to_results_human,
                            only_shared = only_shared)
      TPR_human <- calc_TPR(methods = methods_human, method_names = method_names, 
                            method_cols = method_cols, organism = "Human", 
                            thresholds = thresholds, truth = truth_human,
                            path_to_results = path_to_results_human,
                            only_shared = only_shared)
    } else {
      FDR_human <- NULL
      TPR_human <- NULL
    }
  } else {
    if (!is.null(path_to_results_drosophila) && (!is.null(truth_drosophila))) {
      FDR_drosophila <- lapply(split(truth_drosophila, 
                                     truth_drosophila[, split_variable]), 
                               function(w) calc_FDR(methods = methods_drosophila, 
                                                    method_names = method_names, 
                                                    method_cols = method_cols, 
                                                    organism = "Drosophila", 
                                                    thresholds = thresholds, 
                                                    truth = w,
                                                    path_to_results = path_to_results_drosophila,
                                                    only_shared = only_shared))
      TPR_drosophila <- lapply(split(truth_drosophila, 
                                     truth_drosophila[, split_variable]), 
                               function(w) calc_TPR(methods = methods_drosophila, 
                                                    method_names = method_names, 
                                                    method_cols = method_cols, 
                                                    organism = "Drosophila", 
                                                    thresholds = thresholds, 
                                                    truth = w,
                                                    path_to_results = path_to_results_drosophila,
                                                    only_shared = only_shared))
    } else {
      FDR_drosophila <- NULL
      TPR_drosophila <- NULL
    }
    if (!is.null(path_to_results_human) && (!is.null(truth_human))) {
      FDR_human <- lapply(split(truth_human, 
                                truth_human[, split_variable]), 
                          function(w) calc_FDR(methods = methods_human, 
                                               method_names = method_names, 
                                               method_cols = method_cols, 
                                               organism = "Human", 
                                               thresholds = thresholds, 
                                               truth = w,
                                               path_to_results = path_to_results_human,
                                               only_shared = only_shared))
      TPR_human <- lapply(split(truth_human, 
                                truth_human[, split_variable]), 
                          function(w) calc_TPR(methods = methods_human, 
                                               method_names = method_names, 
                                               method_cols = method_cols, 
                                               organism = "Human", 
                                               thresholds = thresholds, 
                                               truth = w,
                                               path_to_results = path_to_results_human,
                                               only_shared = only_shared))
    } else {
      FDR_human <- NULL
      TPR_human <- NULL
    }
  }
  if (!is.null(FDR_drosophila))
    FDR_drosophila <- melt(FDR_drosophila, variable.name = "threshold", 
                           value.name = "FDR")
  if (!is.null(FDR_human))
    FDR_human <- melt(FDR_human, variable.name = "threshold", 
                      value.name = "FDR")
  if (!is.null(TPR_drosophila))
    TPR_drosophila <- melt(TPR_drosophila, variable.name = "threshold", 
                           value.name = "TPR")
  if (!is.null(TPR_human))
    TPR_human <- melt(TPR_human, variable.name = "threshold", 
                      value.name = "TPR")
  
  ## Put results together
  if (!is.null(FDR_drosophila))
    FDR_TPR_drosophila <- merge(FDR_drosophila, TPR_drosophila)
  else
    FDR_TPR_drosophila <- NULL
  if (!is.null(FDR_human))
    FDR_TPR_human <- merge(FDR_human, TPR_human)
  else
    FDR_TPR_human <- NULL
  FDR_TPR <- rbind(FDR_TPR_drosophila, FDR_TPR_human)
  FDR_TPR$threshold <- as.numeric(as.character(FDR_TPR$threshold))
  FDR_TPR$fill_color <- FDR_TPR$method_col
  FDR_TPR$fill_color[FDR_TPR$FDR > FDR_TPR$threshold] <- "white"
  fill_color <- unique(FDR_TPR$fill_color)
  names(fill_color) <- fill_color
  uniq <- !duplicated(FDR_TPR$method_name)
  col_color <- FDR_TPR$method_col[uniq]
  names(col_color) <- FDR_TPR$method_name[uniq]
  
  write.table(FDR_TPR, file = gsub("\\.pdf$", ".txt", output_filename), 
              row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
  if (is.null(split_variable)) {
    pdf(output_filename, width = fig_width, height = fig_height)
    print(ggplot(FDR_TPR, aes(x = FDR, y = TPR, group = method_name)) + 
            geom_vline(xintercept = seq(0, 1, 0.1), 
                       colour = "lightgrey", linetype = "dashed") + 
            geom_vline(xintercept = thresholds, linetype = "dashed") + 
            geom_path(size = 1, aes(colour = method_name)) + 
            facet_wrap(~ organism) + 
            geom_point(size = pointsize + 1, 
                       aes(colour = method_name), shape = 19) + 
            geom_point(size = pointsize, 
                       aes(fill = fill_color, colour = method_name), shape = 21) + 
            scale_fill_manual(values = fill_color, guide = FALSE) + 
            scale_color_manual(values = col_color, name = "") +
            guides(col = guide_legend(nrow = legend_nrow)) + 
            ylim(0, 1) +
            scale_x_continuous(breaks = c(0, thresholds, seq(0.1, 1, 0.1)),
                               labels = c("", thresholds, seq(0.1, 1, 0.1)),
                               limits = c(0, 1)) + 
            theme(legend.position = "bottom", 
                  panel.background = element_rect(fill = NA, colour = "black"),
                  panel.border = element_rect(fill = NA, colour = "black", size = 1),
                  panel.grid.minor.x = element_blank(),
                  panel.grid.minor.y = element_blank(),
                  axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, 
                                             size = axistextsize),
                  axis.text.y = element_text(size = axistextsize),
                  axis.title.x = element_text(size = axistitlesize),
                  axis.title.y = element_text(size = axistitlesize),
                  strip.text = element_text(size = stripsize),
                  strip.background = element_rect(fill = NA, colour = "black")) + 
            ggtitle(""))
    dev.off()
  } else {
    pdf(output_filename, width = fig_width, height = fig_height)
    print(ggplot(FDR_TPR, aes(x = FDR, y = TPR, group = method_name)) + 
            geom_vline(xintercept = seq(0, 1, 0.1), 
                       colour = "lightgrey", linetype = "dashed") + 
            geom_vline(xintercept = thresholds, linetype = "dashed") + 
            geom_path(size = 1, aes(colour = method_name)) + 
            facet_wrap(organism ~ L1) + 
            geom_point(size = pointsize + 1, 
                       aes(colour = method_name), shape = 19) + 
            geom_point(size = pointsize, 
                       aes(fill = fill_color, colour = method_name), shape = 21) + 
            scale_fill_manual(values = fill_color, guide = FALSE) + 
            scale_color_manual(values = col_color, name = "") +
            guides(col = guide_legend(nrow = legend_nrow)) + 
            ylim(0, 1) +
            scale_x_continuous(breaks = c(0, thresholds, seq(0.1, 1, 0.1)),
                               labels = c("", thresholds, seq(0.1, 1, 0.1)),
                               limits = c(0, 1)) + 
            theme(legend.position = "bottom", 
                  panel.background = element_rect(fill = NA, colour = "black"),
                  panel.border = element_rect(fill = NA, colour = "black", size = 1),
                  panel.grid.minor.x = element_blank(),
                  panel.grid.minor.y = element_blank(),
                  axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, 
                                             size = axistextsize),
                  axis.text.y = element_text(size = axistextsize),
                  axis.title.x = element_text(size = axistitlesize),
                  axis.title.y = element_text(size = axistitlesize),
                  strip.text = element_text(size = stripsize),
                  strip.background = element_rect(fill = NA, colour = "black")) + 
            ggtitle(""))
    dev.off()
  }
}

## With different pch
plot_fdr_tpr_paper_pch <- function(path_to_results_drosophila, 
                                   path_to_results_human,
                                   truth_drosophila, truth_human, 
                                   methods = NULL, method_names, 
                                   methods_drosophila = NULL,
                                   methods_human = NULL, 
                                   method_pchnames, method_colnames, 
                                   method_pch, 
                                   method_cols, thresholds, split_variable,
                                   output_filename,
                                   pointsize = 3, stripsize = 20,
                                   axistitlesize = 20, axistextsize = 15,
                                   legend_nrow = 2, fig_height = 7,
                                   fig_width = 10, only_shared = FALSE) {
  
  if (is.null(methods_drosophila)) {
    methods_drosophila <- methods
  }
  if (is.null(methods_human)) {
    methods_human <- methods
  }
  
  ## Calculate FDR 
  if (is.null(split_variable)) {
    if (!is.null(path_to_results_drosophila) && (!is.null(truth_drosophila))) {
      FDR_drosophila <- calc_FDR(methods = methods_drosophila, method_names = method_names, 
                                 method_cols = method_colnames, 
                                 organism = "Drosophila", 
                                 thresholds = thresholds, 
                                 truth = truth_drosophila,
                                 path_to_results = path_to_results_drosophila, 
                                 method_pch = method_pchnames,
                                 only_shared = only_shared)
      TPR_drosophila <- calc_TPR(methods = methods_drosophila, method_names = method_names, 
                                 method_cols = method_colnames, 
                                 organism = "Drosophila", 
                                 thresholds = thresholds, 
                                 truth = truth_drosophila,
                                 path_to_results = path_to_results_drosophila, 
                                 method_pch = method_pchnames,
                                 only_shared = only_shared)
    } else {
      FDR_drosophila <- NULL
      TPR_drosophila <- NULL
    }
    if (!is.null(path_to_results_human) && (!is.null(truth_human))) {
      FDR_human <- calc_FDR(methods = methods_human, method_names = method_names, 
                            method_cols = method_colnames, 
                            organism = "Human", 
                            thresholds = thresholds, 
                            truth = truth_human,
                            path_to_results = path_to_results_human, 
                            method_pch = method_pchnames,
                            only_shared = only_shared)
      TPR_human <- calc_TPR(methods = methods_human, method_names = method_names, 
                            method_cols = method_colnames, 
                            organism = "Human", 
                            thresholds = thresholds, 
                            truth = truth_human,
                            path_to_results = path_to_results_human, 
                            method_pch = method_pchnames,
                            only_shared = only_shared)
    } else {
      FDR_human <- NULL
      TPR_human <- NULL
    }
  } else {
    if (!is.null(path_to_results_drosophila) && (!is.null(truth_drosophila))) {
      FDR_drosophila <- lapply(split(truth_drosophila, 
                                     truth_drosophila[, split_variable]), 
                               function(w) calc_FDR(methods = methods_drosophila, 
                                                    method_names = method_names, 
                                                    method_cols = method_colnames, 
                                                    organism = "Drosophila", 
                                                    thresholds = thresholds, 
                                                    truth = w,
                                                    path_to_results = path_to_results_drosophila,
                                                    method_pch = method_pchnames,
                                                    only_shared = only_shared))
      TPR_drosophila <- lapply(split(truth_drosophila, 
                                     truth_drosophila[, split_variable]), 
                               function(w) calc_TPR(methods = methods_drosophila, 
                                                    method_names = method_names, 
                                                    method_cols = method_colnames, 
                                                    organism = "Drosophila", 
                                                    thresholds = thresholds, 
                                                    truth = w,
                                                    path_to_results = path_to_results_drosophila,
                                                    method_pch = method_pchnames,
                                                    only_shared = only_shared))
    } else {
      FDR_drosophila <- NULL
      TPR_drosophila <- NULL
    }
    if (!is.null(path_to_results_human) && (!is.null(truth_human))) {
      FDR_human <- lapply(split(truth_human, 
                                truth_human[, split_variable]), 
                          function(w) calc_FDR(methods = methods_human, 
                                               method_names = method_names, 
                                               method_cols = method_colnames, 
                                               organism = "Human", 
                                               thresholds = thresholds, 
                                               truth = w,
                                               path_to_results = path_to_results_human,
                                               method_pch = method_pchnames,
                                               only_shared = only_shared))
      TPR_human <- lapply(split(truth_human, 
                                truth_human[, split_variable]), 
                          function(w) calc_TPR(methods = methods_human, 
                                               method_names = method_names, 
                                               method_cols = method_colnames, 
                                               organism = "Human", 
                                               thresholds = thresholds, 
                                               truth = w,
                                               path_to_results = path_to_results_human,
                                               method_pch = method_pchnames,
                                               only_shared = only_shared))
    } else {
      FDR_human <- NULL
      TPR_human <- NULL
    }
  }
  if (!is.null(FDR_drosophila))
    FDR_drosophila <- melt(FDR_drosophila, variable.name = "threshold", 
                           value.name = "FDR", measure.vars = as.character(thresholds))
  if (!is.null(FDR_human))
    FDR_human <- melt(FDR_human, variable.name = "threshold", 
                      value.name = "FDR", measure.vars = as.character(thresholds))
  if (!is.null(TPR_drosophila))
    TPR_drosophila <- melt(TPR_drosophila, variable.name = "threshold", 
                           value.name = "TPR", measure.vars = as.character(thresholds))
  if (!is.null(TPR_human))
    TPR_human <- melt(TPR_human, variable.name = "threshold", 
                      value.name = "TPR", measure.vars = as.character(thresholds))

  ## Put results together
  if (!is.null(FDR_drosophila))
    FDR_TPR_drosophila <- merge(FDR_drosophila, TPR_drosophila)
  else
    FDR_TPR_drosophila <- NULL
  if (!is.null(FDR_human))
    FDR_TPR_human <- merge(FDR_human, TPR_human)
  else
    FDR_TPR_human <- NULL
  FDR_TPR <- rbind(FDR_TPR_drosophila, FDR_TPR_human)
  FDR_TPR$threshold <- as.numeric(as.character(FDR_TPR$threshold))
  FDR_TPR$fill_color <- method_cols[match(FDR_TPR$method_col, method_colnames)]
  FDR_TPR$fill_color[FDR_TPR$FDR > FDR_TPR$threshold] <- "white"
  fill_color <- unique(FDR_TPR$fill_color)
  names(fill_color) <- fill_color
  uniq <- !duplicated(method_colnames)
  col_color <- method_cols[uniq]
  names(col_color) <- method_colnames[uniq]
  uniq2 <- !duplicated(method_pchnames)
  pchs <- method_pch[uniq2]
  names(pchs) <- method_pchnames[uniq2]

  write.table(FDR_TPR, file = gsub("\\.pdf$", ".txt", output_filename), 
              row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
  if (is.null(split_variable)) {
    pdf(output_filename, width = fig_width, height = fig_height)
    print(ggplot(FDR_TPR, aes(x = FDR, y = TPR, group = method_name)) + 
            geom_vline(xintercept = seq(0, 1, 0.1), 
                       colour = "lightgrey", linetype = "dashed") + 
            geom_vline(xintercept = thresholds, linetype = "dashed") + 
            geom_path(size = 1, aes(colour = method_col)) + 
            facet_wrap(~ organism) + 
            geom_point(size = pointsize + 1, 
                       aes(colour = method_col, shape = method_pch)) + 
            geom_point(size = pointsize + 0.5, 
                       aes(colour = method_col, shape = method_pch)) + 
            geom_point(size = pointsize, 
                       aes(fill = fill_color, colour = method_col, shape = method_pch)) + 
            scale_fill_manual(values = fill_color, guide = FALSE) + 
            scale_color_manual(values = col_color, name = "") +
            scale_shape_manual(values = pchs, name = "") + 
            guides(col = guide_legend(nrow = legend_nrow)) + 
            ylim(0, 1) +
            scale_x_continuous(breaks = c(0, thresholds, seq(0.1, 1, 0.1)),
                               labels = c("", thresholds, seq(0.1, 1, 0.1)),
                               limits = c(0, 1)) + 
            theme(legend.position = "bottom", 
                  panel.background = element_rect(fill = NA, colour = "black"),
                  panel.border = element_rect(fill = NA, colour = "black", size = 1),
                  panel.grid.minor.x = element_blank(),
                  panel.grid.minor.y = element_blank(),
                  axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, 
                                             size = axistextsize),
                  axis.text.y = element_text(size = axistextsize),
                  axis.title.x = element_text(size = axistitlesize),
                  axis.title.y = element_text(size = axistitlesize),
                  strip.text = element_text(size = stripsize),
                  strip.background = element_rect(fill = NA, colour = "black")) + 
            ggtitle(""))
    dev.off()
  } else {
    pdf(output_filename, width = fig_width, height = fig_height)
    print(ggplot(FDR_TPR, aes(x = FDR, y = TPR, group = method_name)) + 
            geom_vline(xintercept = seq(0, 1, 0.1), 
                       colour = "lightgrey", linetype = "dashed") + 
            geom_vline(xintercept = thresholds, linetype = "dashed") + 
            geom_path(size = 1, aes(colour = method_col)) + 
            facet_wrap(organism ~ L1) + 
            geom_point(size = pointsize + 1, 
                       aes(colour = method_col, shape = method_pch)) + 
            geom_point(size = pointsize + 0.5, 
                       aes(colour = method_col, shape = method_pch)) + 
            geom_point(size = pointsize, 
                       aes(fill = fill_color, colour = method_col, shape = method_pch)) + 
            scale_fill_manual(values = fill_color, guide = FALSE) + 
            scale_color_manual(values = col_color, name = "") +
            scale_shape_manual(values = pchs, name = "") + 
            guides(col = guide_legend(nrow = legend_nrow)) + 
            ylim(0, 1) +
            scale_x_continuous(breaks = c(0, thresholds, seq(0.1, 1, 0.1)),
                               labels = c("", thresholds, seq(0.1, 1, 0.1)),
                               limits = c(0, 1)) + 
            theme(legend.position = "bottom", 
                  panel.background = element_rect(fill = NA, colour = "black"),
                  panel.border = element_rect(fill = NA, colour = "black", size = 1),
                  panel.grid.minor.x = element_blank(),
                  panel.grid.minor.y = element_blank(),
                  axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, 
                                             size = axistextsize),
                  axis.text.y = element_text(size = axistextsize),
                  axis.title.x = element_text(size = axistitlesize),
                  axis.title.y = element_text(size = axistitlesize),
                  strip.text = element_text(size = stripsize),
                  strip.background = element_rect(fill = NA, colour = "black")) + 
            ggtitle(""))
    dev.off()
  }
}

plot_tpr_paper <- function(path_to_results_drosophila, 
                           path_to_results_human,
                           truth_drosophila, truth_human, 
                           methods = NULL, method_names, 
                           methods_drosophila = NULL,
                           methods_human = NULL,
                           method_cols, thresholds, split_variable,
                           output_filename,
                           pointsize = 3, stripsize = 20,
                           axistitlesize = 20, axistextsize = 15,
                           legend_nrow = 2, only_shared = FALSE,
                           fig_width = 11, fig_height = 7) {
  
  if (is.null(methods_drosophila)) {
    methods_drosophila <- methods
  }
  if (is.null(methods_human)) {
    methods_human <- methods
  }
  
  ## Calculate TPR
  if (is.null(split_variable)) {
    if (!is.null(truth_drosophila) && !is.null(path_to_results_drosophila))
      TPR_drosophila <- calc_TPR(methods = methods_drosophila, method_names = method_names, 
                                 method_cols = method_cols, 
                                 organism = "Drosophila", 
                                 thresholds = thresholds, truth = truth_drosophila,
                                 path_to_results = path_to_results_drosophila,
                                 only_shared = only_shared)
    else
      TPR_drosophila <- NULL
    if (!is.null(truth_human) && !is.null(path_to_results_human))
      TPR_human <- calc_TPR(methods = methods_human, method_names = method_names, 
                            method_cols = method_cols, organism = "Human", 
                            thresholds = thresholds, truth = truth_human,
                            path_to_results = path_to_results_human, 
                            only_shared = only_shared)
    else
      TPR_human <- NULL
  } else {
    if (!is.null(truth_drosophila) && !is.null(path_to_results_drosophila))
      TPR_drosophila <- lapply(split(truth_drosophila, 
                                     truth_drosophila[, split_variable]), 
                               function(w) calc_TPR(methods = methods_drosophila, 
                                                    method_names = method_names, 
                                                    method_cols = method_cols, 
                                                    organism = "Drosophila", 
                                                    thresholds = thresholds, 
                                                    truth = w,
                                                    path_to_results = path_to_results_drosophila,
                                                    only_shared = only_shared))
    else 
      TPR_drosophila <- NULL
    if (!is.null(truth_human) && !is.null(path_to_results_human))
      TPR_human <- lapply(split(truth_human, 
                                truth_human[, split_variable]), 
                          function(w) calc_TPR(methods = methods_human, 
                                               method_names = method_names, 
                                               method_cols = method_cols, 
                                               organism = "Human", 
                                               thresholds = thresholds, 
                                               truth = w,
                                               path_to_results = path_to_results_human,
                                               only_shared = only_shared))
    else
      TPR_human <- NULL
  }
  if (!is.null(TPR_drosophila))
    TPR_drosophila <- melt(TPR_drosophila, variable.name = "threshold", 
                           value.name = "TPR")
  if (!is.null(TPR_human))
    TPR_human <- melt(TPR_human, variable.name = "threshold", 
                      value.name = "TPR")
  
  ## Put results together
  TPR <- rbind(TPR_drosophila, TPR_human)
  TPR$threshold <- as.numeric(as.character(TPR$threshold))
  uniq <- !duplicated(TPR$method_name)
  col_color <- TPR$method_col[uniq]
  names(col_color) <- TPR$method_name[uniq]
  
  write.table(TPR, file = gsub("\\.pdf$", ".txt", output_filename), 
              row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
  if (is.null(split_variable)) {
    pdf(output_filename, width = fig_width, height = fig_height)
    print(ggplot(TPR, aes(x = TPR, y = method_name, group = method_name)) + 
            geom_path(size = 1, aes(colour = method_name)) + 
            facet_wrap(~ organism) + 
            geom_point(size = pointsize + 1, 
                       aes(colour = method_name), shape = 19) + 
            scale_color_manual(values = col_color, name = "") +
            guides(col = guide_legend(nrow = legend_nrow)) + 
            xlim(0, 1) + ylab("") +
            theme(legend.position = "bottom", 
                  panel.grid.major.x = element_line(colour = "lightgrey", 
                                                    linetype = "dotted"),
                  panel.grid.minor.x = element_blank(),
                  panel.grid.minor.y = element_blank(),
                  panel.background = element_rect(fill = NA, colour = "black"),
                  panel.border = element_rect(fill = NA, colour = "black", size = 1),
                  axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, 
                                             size = axistextsize),
                  axis.text.y = element_text(size = axistextsize),
                  axis.title.x = element_text(size = axistitlesize),
                  axis.title.y = element_text(size = axistitlesize),
                  strip.text = element_text(size = stripsize),
                  strip.background = element_rect(fill = NA, colour = "black")) + 
            ggtitle(""))
    dev.off()
  } else {
    pdf(output_filename, width = fig_width, height = fig_height)
    print(ggplot(TPR, aes(x = TPR, y = method_name, group = method_name)) + 
            geom_path(size = 1, aes(colour = method_name)) + 
            facet_wrap(organism ~ L1) + 
            geom_point(size = pointsize + 1, 
                       aes(colour = method_name), shape = 19) + 
            scale_color_manual(values = col_color, name = "") +
            guides(col = guide_legend(nrow = legend_nrow)) + 
            xlim(0, 1) + ylab("") + 
            theme(legend.position = "bottom", 
                  panel.grid.major.x = element_line(colour = "lightgrey", 
                                                    linetype = "dotted"),
                  panel.grid.minor.x = element_blank(),
                  panel.grid.minor.y = element_blank(),
                  panel.background = element_rect(fill = NA, colour = "black"),
                  panel.border = element_rect(fill = NA, colour = "black", size = 1),
                  axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, 
                                             size = axistextsize),
                  axis.text.y = element_text(size = axistextsize),
                  axis.title.x = element_text(size = axistitlesize),
                  axis.title.y = element_text(size = axistitlesize),
                  strip.text = element_text(size = stripsize),
                  strip.background = element_rect(fill = NA, colour = "black")) + 
            ggtitle(""))
    dev.off()
  }
}

plot_tpr_paper_pch <- function(path_to_results_drosophila, 
                               path_to_results_human,
                               truth_drosophila, truth_human, 
                               methods = NULL, method_names, method_pchnames,
                               methods_drosophila = NULL,
                               methods_human = NULL,
                               method_colnames, method_pch, 
                               method_cols, thresholds, split_variable,
                               output_filename,
                               pointsize = 3, stripsize = 20,
                               axistitlesize = 20, axistextsize = 15,
                               legend_nrow = 2, only_shared = FALSE,
                               fig_height = 7, fig_width = 11) {
  
  if (is.null(methods_drosophila)) {
    methods_drosophila <- methods
  }
  if (is.null(methods_human)) {
    methods_human <- methods
  }
  
  ## Calculate TPR
  if (is.null(split_variable)) {
    if (!is.null(truth_drosophila) && !is.null(path_to_results_drosophila))
      TPR_drosophila <- calc_TPR(methods = methods_drosophila, method_names = method_names, 
                                 method_cols = method_colnames, 
                                 organism = "Drosophila", 
                                 thresholds = thresholds, 
                                 truth = truth_drosophila,
                                 path_to_results = path_to_results_drosophila, 
                                 method_pch = method_pchnames,
                                 only_shared = only_shared)
    else
      TPR_drosophila <- NULL
    if (!is.null(truth_human) && !is.null(path_to_results_human))
      TPR_human <- calc_TPR(methods = methods_human, method_names = method_names, 
                            method_cols = method_colnames, 
                            organism = "Human", 
                            thresholds = thresholds, truth = truth_human,
                            path_to_results = path_to_results_human, 
                            method_pch = method_pchnames,
                            only_shared = only_shared)
    else
      TPR_human <- NULL
  } else {
    if (!is.null(truth_drosophila) && !is.null(path_to_results_drosophila))
      TPR_drosophila <- lapply(split(truth_drosophila, 
                                     truth_drosophila[, split_variable]), 
                               function(w) calc_TPR(methods = methods, 
                                                    method_names = method_names, 
                                                    method_cols = method_colnames, 
                                                    organism = "Drosophila", 
                                                    thresholds = thresholds, 
                                                    truth = w,
                                                    path_to_results = path_to_results_drosophila,
                                                    method_pch = method_pchnames,
                                                    only_shared = only_shared))
    else 
      TPR_drosophila <- NULL
    if (!is.null(truth_human) && !is.null(path_to_results_human))
      TPR_human <- lapply(split(truth_human, 
                                truth_human[, split_variable]), 
                          function(w) calc_TPR(methods = methods, 
                                               method_names = method_names, 
                                               method_cols = method_colnames,
                                               organism = "Human", 
                                               thresholds = thresholds, 
                                               truth = w,
                                               path_to_results = path_to_results_human,
                                               method_pch = method_pchnames,
                                               only_shared = only_shared))
    else
      TPR_human <- NULL
  }
  if (!is.null(TPR_drosophila))
    TPR_drosophila <- melt(TPR_drosophila, variable.name = "threshold", 
                           value.name = "TPR", measure.vars = as.character(thresholds))
  if (!is.null(TPR_human))
    TPR_human <- melt(TPR_human, variable.name = "threshold", 
                      value.name = "TPR", measure.vars = as.character(thresholds))
  
  ## Put results together
  TPR <- rbind(TPR_drosophila, TPR_human)
  TPR$threshold <- as.numeric(as.character(TPR$threshold))
  uniq <- !duplicated(method_colnames)
  col_color <- method_cols[uniq]
  names(col_color) <- method_colnames[uniq]
  uniq2 <- !duplicated(method_pchnames)
  pchs <- method_pch[uniq2]
  names(pchs) <- method_pchnames[uniq2]
#   
#   uniq <- !duplicated(TPR$method_name)
#   col_color <- TPR$method_col[uniq]
#   names(col_color) <- TPR$method_name[uniq]
  
  write.table(TPR, file = gsub("\\.pdf$", ".txt", output_filename), 
              row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
  if (is.null(split_variable)) {
    pdf(output_filename, width = fig_width, height = fig_height)
    print(ggplot(TPR, aes(x = TPR, y = method_name, group = method_name)) + 
            geom_path(size = 1, aes(colour = method_col)) + 
            facet_wrap(~ organism) + 
            geom_point(size = pointsize + 1, 
                       aes(colour = method_col, shape = method_pch)) + 
            scale_color_manual(values = col_color, name = "") +
            scale_shape_manual(values = pchs, name = "") + 
            guides(col = guide_legend(nrow = legend_nrow)) + 
            xlim(0, 1) + ylab("") +
            theme(legend.position = "bottom", 
                  panel.grid.major.x = element_line(colour = "lightgrey", 
                                                    linetype = "dotted"),
                  panel.grid.minor.x = element_blank(),
                  panel.grid.minor.y = element_blank(),
                  panel.background = element_rect(fill = NA, colour = "black"),
                  panel.border = element_rect(fill = NA, colour = "black", size = 1),
                  axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, 
                                             size = axistextsize),
                  axis.text.y = element_text(size = axistextsize),
                  axis.title.x = element_text(size = axistitlesize),
                  axis.title.y = element_text(size = axistitlesize),
                  strip.text = element_text(size = stripsize),
                  strip.background = element_rect(fill = NA, colour = "black")) + 
            ggtitle(""))
    dev.off()
  } else {
    pdf(output_filename, width = fig_width, height = fig_height)
    print(ggplot(TPR, aes(x = TPR, y = method_name, group = method_name)) + 
            geom_path(size = 1, aes(colour = method_col)) + 
            facet_wrap(organism ~ L1) + 
            geom_point(size = pointsize + 1, 
                       aes(colour = method_col, shape = method_pch)) + 
            scale_color_manual(values = col_color, name = "") +
            scale_shape_manual(values = pchs, name = "") + 
            guides(col = guide_legend(nrow = legend_nrow)) + 
            xlim(0, 1) + ylab("") + 
            theme(legend.position = "bottom", 
                  panel.grid.major.x = element_line(colour = "lightgrey", 
                                                    linetype = "dotted"),
                  panel.grid.minor.x = element_blank(),
                  panel.grid.minor.y = element_blank(),
                  panel.background = element_rect(fill = NA, colour = "black"),
                  panel.border = element_rect(fill = NA, colour = "black", size = 1),
                  axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, 
                                             size = axistextsize),
                  axis.text.y = element_text(size = axistextsize),
                  axis.title.x = element_text(size = axistitlesize),
                  axis.title.y = element_text(size = axistitlesize),
                  strip.text = element_text(size = stripsize),
                  strip.background = element_rect(fill = NA, colour = "black")) + 
            ggtitle(""))
    dev.off()
  }
}

plot_perftable <- function(path_to_results, truth, 
                           methods, method_names, 
                           method_cols, threshold, split_variable,
                           output_filename, 
                           only_shared = FALSE) {
  
  pdf(output_filename, width = 12, height = 5)
  trt <- split(truth, truth[, split_variable])
  for (m in seq_along(methods)) {
    A <- do.call(rbind, lapply(trt, function(w1) {
      tmp1 <- read.delim(paste0(path_to_results, "/", methods[m], ".txt"), 
                         header = TRUE, as.is = TRUE, check.names = FALSE)
      if (isTRUE(only_shared)) {
        ist <- intersect(tmp1$gene, w1$gene)
        tmp <- subset(tmp1, gene %in% ist)
        w <- subset(w1, gene %in% ist)
      } else {
        w <- w1
        tmp <- tmp1
      }
      ds <- length(which(w$ds_status == 1))
      nonds <- length(which(w$ds_status == 0))
      
      TP <- length(intersect(w$gene[w$ds_status == 1], 
                             tmp$gene[tmp[, 2] <= threshold]))
      FP <- length(intersect(w$gene[w$ds_status == 0], 
                             tmp$gene[tmp[, 2] <= threshold]))
      FN <- length(intersect(w$gene[w$ds_status == 1], 
                             tmp$gene[tmp[, 2] > threshold]))
      TN <- length(intersect(w$gene[w$ds_status == 0], 
                             tmp$gene[tmp[, 2] > threshold]))
      uncl.P <- length(intersect(w$gene[!(w$ds_status %in% c(0, 1))],
                                 tmp$gene[tmp[, 2] > threshold]))
      uncl.N <- length(intersect(w$gene[!(w$ds_status %in% c(0, 1))],
                                 tmp$gene[tmp[, 2] > threshold]))
      tot.called <- length(intersect(w$gene,
                                     tmp$gene[which(!is.na(tmp[, 2]))]))
      ##tot.called <- TP + FP + FN + TN
      tot.status <- ds + nonds
      c(ds = ds, nonds = nonds, tot.status = tot.status, 
        TP = TP, FP = FP, FN = FN, TN = TN, uncl.P = uncl.P,
        uncl.N = uncl.N, tot.called = tot.called)
    }))
    tmp1 <- read.delim(paste0(path_to_results, "/", methods[m], ".txt"), 
                       header = TRUE, as.is = TRUE, check.names = FALSE)
    A <- rbind(A, c(ds = "-", nonds = "-", tot.status = "-",
                    TP = "-", FP = "-", FN = "-", TN = "-",
                    uncl.P = length(intersect(tmp1$gene[tmp1[, 2] <= threshold],
                                              setdiff(tmp1$gene, truth$gene))),
                    uncl.N = length(intersect(tmp1$gene[tmp1[, 2] > threshold],
                                              setdiff(tmp1$gene, truth$gene))),
                    tot.called = length(intersect(tmp1$gene[tmp1[, 2] <= threshold],
                                                  setdiff(tmp1$gene, truth$gene))) + 
                      length(intersect(tmp1$gene[tmp1[, 2] > threshold],
                                       setdiff(tmp1$gene, truth$gene)))))
    B <- melt(A)
    B$fillcolor <- "white"
    B$fillcolor[B$Var2 == "ds"] <- "lightgreen"
    B$fillcolor[B$Var2 == "nonds"] <- "pink"
    fc <- unique(B$fillcolor)
    names(fc) <- fc
    print(ggplot(B, aes(Var2, Var1, fill = fillcolor)) + 
            geom_tile() + 
            geom_text(aes(label = value), color = "black", size = 7) + 
            scale_fill_manual(values = fc, name = "") + 
            scale_x_discrete(expand = c(0, 0)) + 
            scale_y_discrete(expand = c(0, 0)) + 
            xlab("") + ylab(split_variable) + 
            geom_vline(xintercept = 3.5, linetype = "dashed") + 
            theme(panel.background = element_rect(fill = NA, colour = NA),
                  legend.position = "none",
                  axis.ticks = element_blank(),
                  plot.title = element_text(size = 20),
                  axis.text.x = element_text(size = 15),
                  axis.text.y = element_text(size = 15),
                  axis.title.y = element_text(size = 15)) + 
            ggtitle(method_names[m]))
  }
  dev.off()
}
  
  